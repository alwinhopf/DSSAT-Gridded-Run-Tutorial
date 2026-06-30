# =============================================================================
# validate_against_observed.R — observed vs. original-DSSAT vs. pipeline yields.
#
# For each DSSAT example experiment that ships with OBSERVED data (a FileA, e.g.
# UFGA8201.MZA), this compares three pathways of end-of-season grain yield (HWAM):
#
#   1. OBSERVED      — measured HWAM per treatment, read from the FileA.
#   2. ORIGINAL      — DSSAT run of the experiment with its OWN local weather
#                      (.WTH) and soil, i.e. the calibrated "textbook" run.
#   3. PIPELINE      — the same experiment run at its own coordinates but with
#                      this pipeline's GRIDDED weather x soil inputs (the full
#                      sensitivity sweep), so you can also see which gridded
#                      input best reproduces the observed / locally-simulated yield.
#
# It writes a tidy long table, a metrics table (RMSE, nRMSE, bias, Willmott d,
# modelling efficiency, R2), and comparison figures.
#
# USAGE
#   Rscript validate_against_observed.R                 # all configured experiments
#   Rscript validate_against_observed.R UFGA8201        # one (or several) by name
#   Rscript validate_against_observed.R --dry-run       # resolve metadata only, run nothing
#
# Config is the small block below (paths, crop, which experiments, parallelism).
# The pipeline sweep is delegated to run_experiment.R, and the original run uses
# the DSSAT R package — the same engine the pipeline itself uses.
# =============================================================================

suppressWarnings(suppressMessages({
  for (p in c("yaml")) if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
}))

# ---- Configuration --------------------------------------------------------
PROJECT_ROOT <- tryCatch({
  if (requireNamespace("this.path", quietly = TRUE)) dirname(this.path::this.path()) else getwd()
}, error = function(e) getwd())
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, mustWork = FALSE)
setwd(PROJECT_ROOT)

DSSAT_DIR   <- normalizePath(file.path(PROJECT_ROOT, "..", "DSSAT48"), mustWork = FALSE)
DSSAT_EXE   <- Sys.getenv("DSSAT_EXE", unset = file.path(DSSAT_DIR, "dscsm048"))
CROP_DIR    <- file.path(DSSAT_DIR, "Maize")
CROP_EXT    <- "MZ"                       # MZX experiment / MZA observed
OUT_DIR     <- file.path(PROJECT_ROOT, "validation")
MAX_PARALLEL <- 2                         # passed through to the pipeline sweep
MODEL_CODE  <- "MZCER048"                 # CERES-Maize, for the original run

# Experiments to validate (must each have a matching .MZA observed file).
EXPERIMENTS <- c("UFGA8201","BRPI0202","FLSC8101","GAGR0201","GHWA0401",
                 "IBWA8301","IUAF9901","SIAZ9501","SIAZ9601","EBPL8501")

args      <- commandArgs(trailingOnly = TRUE)
DRY_RUN   <- any(args %in% c("--dry-run","-n"))
CHECK_ONLY<- any(args %in% c("--check"))        # report source runnability, then exit
sel       <- args[!grepl("^-", args)]
if (length(sel)) EXPERIMENTS <- sel

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Helpers: parsing -----------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a
yy_to_year <- function(yy) { yy <- as.integer(yy); ifelse(yy > 30, 1900L + yy, 2000L + yy) }

# Observed FileA -> data.frame(treatment, HWAM_obs)
parse_filea <- function(path) {
  ln <- readLines(path, warn = FALSE)
  h  <- grep("^@TRNO", ln)[1]
  if (is.na(h)) return(NULL)
  hdr  <- strsplit(trimws(sub("^@", "", ln[h])), "\\s+")[[1]]
  rows <- list(); i <- h + 1
  while (i <= length(ln) && grepl("^\\s*[0-9]", ln[i])) {
    rows[[length(rows) + 1]] <- strsplit(trimws(ln[i]), "\\s+")[[1]]; i <- i + 1
  }
  if (!length(rows)) return(NULL)
  m  <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(m) <- hdr[seq_len(ncol(m))]
  num <- function(x) suppressWarnings(as.numeric(sub("\\.$", "", x)))
  data.frame(treatment = as.integer(m$TRNO), HWAM_obs = num(m$HWAM))
}

# FileX -> metadata: treatments, planting year, WSTA, soil id, and coordinates
# (preferred from the matching weather-file header, which is unambiguous).
parse_filex <- function(exp) {
  fx <- file.path(CROP_DIR, paste0(exp, ".", CROP_EXT, "X"))
  ln <- readLines(fx, warn = FALSE)
  sect <- function(name) grep(paste0("^\\*", name), ln)[1]
  # treatments
  t0 <- sect("TREATMENTS"); trts <- integer()
  if (!is.na(t0)) for (i in (t0 + 2):length(ln)) {
    if (grepl("^\\*", ln[i])) break
    if (grepl("^\\s*[0-9]", ln[i])) trts <- c(trts, as.integer(strsplit(trimws(ln[i]), "\\s+")[[1]][1]))
  }
  # planting year: PDATE is the 2nd token (after the P level) under *PLANTING DETAILS
  p0 <- sect("PLANTING"); pyear <- NA_integer_
  if (!is.na(p0)) for (i in (p0 + 1):length(ln)) {
    if (grepl("^\\s*[0-9]", ln[i])) { pyear <- yy_to_year(substr(strsplit(trimws(ln[i]), "\\s+")[[1]][2], 1, 2)); break }
  }
  # WSTA + soil id from the *FIELDS data line
  f0 <- sect("FIELDS"); wsta <- NA; soil <- NA
  if (!is.na(f0)) for (i in (f0 + 1):length(ln)) {
    if (grepl("^\\s*[0-9]", ln[i])) {
      tok <- strsplit(trimws(ln[i]), "\\s+")[[1]]
      wsta <- tok[3]
      sid  <- grep("^[A-Z]{2}[A-Z0-9]{8}$", tok, value = TRUE)
      soil <- if (length(sid)) tail(sid, 1) else NA
      break
    }
  }
  # coordinates: prefer the weather-file header (labelled LAT/LONG)
  wprefix <- sub("\\s+$", "", wsta)
  wf <- list.files(file.path(DSSAT_DIR, "Weather"),
                   pattern = paste0("^", substr(wprefix, 1, 4), ".*\\.WTH$"), full.names = TRUE)
  wf <- c(file.path(DSSAT_DIR, "Weather", paste0(exp, ".WTH")), wf)
  wf <- wf[file.exists(wf)][1]
  lat <- NA; lon <- NA
  if (!is.na(wf)) {
    wln <- readLines(wf, warn = FALSE); hh <- grep("@ *INSI", wln)[1]
    if (!is.na(hh) && hh + 1 <= length(wln)) {
      v <- suppressWarnings(as.numeric(strsplit(trimws(wln[hh + 1]), "\\s+")[[1]]))
      v <- v[is.finite(v)]; if (length(v) >= 2) { lat <- v[1]; lon <- v[2] }
    }
  }
  list(exp = exp, filex = fx, treatments = trts, n_treatments = length(trts),
       planting_year = pyear, wsta = wsta, wprefix = wprefix, soil_id = soil,
       lat = lat, lon = lon)
}

# ---- Helpers: region & feasible gridded sources ---------------------------
# ALL_SOURCES: TRUE  = every weather/soil source that is geographically &
#                      temporally feasible for the site (the full comparison).
#              FALSE = a minimal recommended pair (faster, fewer runs).
ALL_SOURCES <- TRUE
PRUNE_UNRUNNABLE <- TRUE        # drop sources whose deps/keys/files are missing (preflight)
has_cds_credentials <- function() {
  file.exists(path.expand("~/.cdsapirc")) || nzchar(Sys.getenv("CDSAPI_RC")) ||
    nzchar(Sys.getenv("CDSAPI_KEY"))
}
HAS_CDS <- has_cds_credentials()
if (!HAS_CDS && ALL_SOURCES && requireNamespace("dssatutils", quietly = TRUE)) {
  HAS_CDS <- tryCatch({
    dssatutils::setup_cds_credentials(prompt = interactive(), quiet = TRUE)
    has_cds_credentials()
  }, error = function(e) {
    if (interactive()) message("Copernicus CDS setup skipped: ", conditionMessage(e))
    FALSE
  })
}
HAS_HWSD <- file.exists(file.path(PROJECT_ROOT, "HWSD", "HWSD2.bil")) || nzchar(Sys.getenv("HWSD_RASTER_FILE"))

region_of <- function(lat, lon) {
  if (is.na(lat) || is.na(lon)) return(NA_character_)
  if (lon >= -125 && lon <= -66 && lat >= 24 && lat <= 50) return("US")     # contiguous US
  if (lon >= -160 && lon <= -154 && lat >= 18 && lat <= 23) return("US")     # Hawaii
  if (lon >= -145 && lon <= -52 && lat >= 14 && lat <= 60) return("NAM")     # rest of N. America
  "GLOBAL"
}
is_conus <- function(lat, lon) lon >= -125 && lon <= -66 && lat >= 24 && lat <= 50

# Every gridded source feasible for a site, honouring each source's coverage:
#   DAYMET (N.Am, 1980+) | GRIDMET (CONUS, 1979+) | OPEN_METEO (global, 1940+)
#   NASA_POWER (global, 1984+) | NASA_POWER_CHIRPS (50S-50N, 1984+)
#   NASA_POWER_CHIRPS_V3 (60S-60N, 1984+) | AGERA5 (global, key)
#   SSURGO/SOILGRIDS_10K (US) | SOILGRIDS_ONLINE (global) | HWSD (global, files)
select_sources <- function(region, lat, lon, year) {
  if (is.na(region)) return(NULL)
  if (!ALL_SOURCES) {                                  # minimal recommended pair
    wx <- "OPEN_METEO"
    if (region %in% c("US","NAM") && year >= 1980) wx <- c("DAYMET", wx)
    if (year >= 1984) wx <- c(wx, "NASA_POWER")
    soil <- if (region == "US") c("SSURGO","SOILGRIDS_10K") else "SOILGRIDS_ONLINE"
    return(list(weather = unique(wx), soil = soil))
  }
  conus <- is_conus(lat, lon); wx <- character(); soil <- character()
  if (region %in% c("US","NAM") && year >= 1980)  wx <- c(wx, "DAYMET")
  if (conus && year >= 1979)                       wx <- c(wx, "GRIDMET")
  if (year >= 1940)                                wx <- c(wx, "OPEN_METEO")
  if (year >= 1984)                                wx <- c(wx, "NASA_POWER")
  if (year >= 1984 && abs(lat) <= 50)              wx <- c(wx, "NASA_POWER_CHIRPS")
  if (year >= 1984 && abs(lat) <= 60)              wx <- c(wx, "NASA_POWER_CHIRPS_V3")
  if (HAS_CDS && year >= 1979)                     wx <- c(wx, "AGERA5")
  if (region == "US")                              soil <- c(soil, "SSURGO", "SOILGRIDS_10K")
  soil <- c(soil, "SOILGRIDS_ONLINE")
  if (HAS_HWSD)                                    soil <- c(soil, "HWSD")
  if (!length(wx) || !length(soil)) return(NULL)
  list(weather = unique(wx), soil = unique(soil))
}

# ---- Metrics (simulated vs observed) --------------------------------------
metrics <- function(obs, sim) {
  ok <- is.finite(obs) & is.finite(sim) & obs > 0
  o <- obs[ok]; s <- sim[ok]; n <- length(o)
  if (n < 1) return(data.frame(n = 0, RMSE = NA, nRMSE_pct = NA, MBE = NA, d = NA, EF = NA, R2 = NA))
  rmse <- sqrt(mean((s - o)^2)); ob <- mean(o)
  d  <- 1 - sum((s - o)^2) / sum((abs(s - ob) + abs(o - ob))^2)
  ef <- 1 - sum((o - s)^2) / sum((o - ob)^2)
  r2 <- if (n > 1 && sd(o) > 0 && sd(s) > 0) cor(o, s)^2 else NA
  data.frame(n = n, RMSE = round(rmse, 1), nRMSE_pct = round(100 * rmse / ob, 1),
             MBE = round(mean(s - o), 1), d = round(d, 3), EF = round(ef, 3),
             R2 = round(r2, 3))
}

# ---- Pathway 2: original DSSAT run (native local weather + soil) -----------
run_original <- function(meta) {
  if (!requireNamespace("DSSAT", quietly = TRUE)) { message("  [original] DSSAT R package missing — skipped."); return(NULL) }
  suppressMessages(library(DSSAT))
  rd <- file.path(OUT_DIR, "original", meta$exp); unlink(rd, recursive = TRUE); dir.create(rd, recursive = TRUE)
  old <- getwd(); on.exit(setwd(old), add = TRUE)
  fx_base <- basename(meta$filex)
  file.copy(meta$filex, file.path(rd, fx_base), overwrite = TRUE)
  # native weather (all files sharing the station prefix) + soil profile
  wfs <- list.files(file.path(DSSAT_DIR, "Weather"),
                    pattern = paste0("^", substr(meta$wprefix, 1, 4)), full.names = TRUE)
  wfs <- c(wfs, file.path(DSSAT_DIR, "Weather", paste0(meta$exp, ".WTH")))
  file.copy(unique(wfs[file.exists(wfs)]), rd, overwrite = TRUE)
  sols <- list.files(file.path(DSSAT_DIR, "Soil"), pattern = "\\.SOL$", full.names = TRUE)
  hit  <- if (!is.na(meta$soil_id)) sols[vapply(sols, function(f) any(grepl(meta$soil_id, readLines(f, warn = FALSE), fixed = TRUE)), logical(1))] else character()
  file.copy(if (length(hit)) hit else file.path(DSSAT_DIR, "Soil", "SOIL.SOL"), rd, overwrite = TRUE)
  setwd(rd)
  options(DSSAT.CSM = DSSAT_EXE)
  ok <- tryCatch({
    write_dssbatch(x = fx_base, trtno = meta$treatments, file_name = file.path(getwd(), "DSSBatch.V48"))
    run_dssat(); TRUE
  }, error = function(e) { message("  [original] run failed: ", conditionMessage(e)); FALSE })
  if (!ok || !file.exists("Summary.OUT")) return(NULL)
  s <- tryCatch(DSSAT::read_dssat("Summary.OUT"), error = function(e) NULL)
  if (is.null(s) || !"HWAM" %in% names(s)) return(NULL)
  data.frame(treatment = as.integer(s$TRNO), HWAM_orig = suppressWarnings(as.numeric(s$HWAM)))
}

# Adapt a raw DSSAT FileX into a pipeline template: insert the placeholders the
# pipeline substitutes per point — WSTA & ID_FIELD -> "00000000" (weather/point
# id), ID_SOIL -> "SOIL_ID" (generated soil), coords -> LATITUDE/LONGITUDE. The
# pipeline replaces the soil and weather; XCRD/YCRD are only metadata.
adapt_template <- function(meta, out_path) {
  ln <- readLines(meta$filex, warn = FALSE)
  field_hdr <- coord_hdr <- NA
  for (i in seq_along(ln)) {
    if (grepl("^@L", ln[i]) && grepl("ID_FIELD", ln[i])) field_hdr <- i
    if (grepl("^@L", ln[i]) && grepl("XCRD",     ln[i])) coord_hdr <- i
  }
  next_data <- function(h) { for (i in (h + 1):length(ln)) if (grepl("^\\s*[0-9]", ln[i])) return(i); NA }
  if (!is.na(field_hdr)) {
    di <- next_data(field_hdr); tk <- strsplit(trimws(ln[di]), "\\s+")[[1]]
    soil_pos <- if (!is.na(meta$soil_id) && meta$soil_id %in% tk) which(tk == meta$soil_id)[1] else 12L
    tk[2] <- "00000000"; tk[3] <- "00000000"                    # ID_FIELD, WSTA placeholders
    if (!is.na(soil_pos) && soil_pos <= length(tk)) tk[soil_pos] <- "SOIL_ID"
    ln[di] <- paste0(" ", paste(tk, collapse = " "))
  }
  if (!is.na(coord_hdr)) {
    di <- next_data(coord_hdr); tk <- strsplit(trimws(ln[di]), "\\s+")[[1]]
    tk[2] <- "LATITUDE"; tk[3] <- "LONGITUDE"
    ln[di] <- paste0(" ", paste(tk, collapse = " "))
  }
  writeLines(ln, out_path)
}

# ---- Pathway 3: pipeline sweep at the experiment's own coordinates ---------
run_pipeline_sweep <- function(meta, srcs) {
  if (!requireNamespace("sf", quietly = TRUE)) { message("  [pipeline] sf missing — skipped."); return(NULL) }
  suppressMessages(library(sf))
  # 1-point shapefile at the site
  pdir <- file.path(OUT_DIR, "points"); dir.create(pdir, showWarnings = FALSE, recursive = TRUE)
  shp  <- file.path(pdir, paste0(meta$exp, ".shp"))
  pt   <- st_as_sf(data.frame(ID = 1L), geometry = st_sfc(st_point(c(meta$lon, meta$lat)), crs = 4326))
  st_write(pt, shp, append = FALSE, delete_layer = TRUE, quiet = TRUE)
  # make the experiment available to the pipeline as a template. If a template
  # of this name already exists (e.g. the bundled, already-adapted demo), reuse
  # it untouched; otherwise write an adapted copy with the pipeline placeholders.
  tdir <- file.path(PROJECT_ROOT, "dssat_templates")
  tmpl <- file.path(tdir, basename(meta$filex))
  if (!file.exists(tmpl)) adapt_template(meta, tmpl)
  for (e in c("CUL","ECO","SPE")) {                       # custom genotype, if any ships with the experiment
    g <- file.path(CROP_DIR, paste0(meta$exp, ".", e)); if (file.exists(g)) file.copy(g, tdir, overwrite = TRUE)
  }
  # validation experiment.yml (MODE B single point, site-appropriate factors)
  combined <- file.path(OUT_DIR, paste0(meta$exp, "_pipeline_combined.csv"))
  yml <- list(
    experiment_name = paste0("val_", meta$exp),
    base = list(project_name = paste0("val_", meta$exp), grid_spacing_meters = 1000,
                crop_extension = CROP_EXT, use_existing_point_shapefile = TRUE,
                existing_point_shapefile_path = shp, enable_boundary_filter = FALSE,
                template_file_name = basename(meta$filex), run_mode = "experiment",
                treatment_start = 1, treatment_end = meta$n_treatments,
                weather_start_year = meta$planting_year, weather_end_year = meta$planting_year),
    factors = list(weather_source = as.list(srcs$weather), soil_source = as.list(srcs$soil)),
    options = list(stop_on_error = FALSE, reuse_existing = TRUE, validate = TRUE,
                   max_parallel = MAX_PARALLEL, response_vars = list("final_grain_kg_ha")),
    output  = list(combined_csv = combined,
                   summary_csv  = file.path(OUT_DIR, paste0(meta$exp, "_pipeline_summary.csv")),
                   plot_png     = file.path(OUT_DIR, paste0(meta$exp, "_pipeline_boxplot.png")),
                   variance_csv = file.path(OUT_DIR, paste0(meta$exp, "_pipeline_variance.csv"))))
  yml_path <- file.path(OUT_DIR, paste0("experiment_val_", meta$exp, ".yml"))
  yaml::write_yaml(yml, yml_path)
  status <- system2("Rscript", args = c(shQuote("run_experiment.R"), shQuote(yml_path)),
                    stdout = file.path(OUT_DIR, paste0(meta$exp, "_pipeline.log")),
                    stderr = file.path(OUT_DIR, paste0(meta$exp, "_pipeline.log")))
  if (!file.exists(combined)) { message("  [pipeline] no combined output — see ", meta$exp, "_pipeline.log"); return(NULL) }
  d <- read.csv(combined, stringsAsFactors = FALSE)
  if (!all(c("treatment","final_grain_kg_ha","weather_source","soil_source") %in% names(d))) return(NULL)
  agg <- aggregate(final_grain_kg_ha ~ treatment + weather_source + soil_source, data = d, FUN = mean)
  names(agg)[names(agg) == "final_grain_kg_ha"] <- "HWAM_pipe"
  agg
}

# ---- Resolve experiment metadata ------------------------------------------
metas <- lapply(EXPERIMENTS, function(e) tryCatch(parse_filex(e), error = function(err) NULL))
metas <- Filter(Negate(is.null), metas)
if (!length(metas)) stop("No experiments resolved — check EXPERIMENTS / CROP_DIR.")

# ---- Preflight: which sources can actually run here? ----------------------
# Deterministic dependency/key/file checks per source (cannot catch transient
# online/API failures — e.g. a SoilGrids REST timeout — only missing prerequisites).
need_pkgs <- function(...) { p <- c(...); m <- p[!vapply(p, function(x) requireNamespace(x, quietly = TRUE), logical(1))]
                             if (length(m)) paste0("install R pkg(s): ", paste(m, collapse = ", ")) else "" }
source_status <- function(src) {
  r <- switch(src,
    DAYMET            = need_pkgs("daymetr"),
    GRIDMET           = need_pkgs("terra","ncdf4","httr"),
    OPEN_METEO        = need_pkgs("httr","jsonlite"),
    NASA_POWER        = need_pkgs("httr","jsonlite"),
    NASA_POWER_CHIRPS = need_pkgs("terra","ncdf4"),
    NASA_POWER_CHIRPS_V3 = need_pkgs("terra","ncdf4"),
    AGERA5            = { m <- need_pkgs("ecmwfr","terra")
                         k <- if (!HAS_CDS) "run setup_cds_credentials() or set CDSAPI_KEY" else ""
                         paste(Filter(nzchar, c(m, k)), collapse = "; ") },
    SSURGO            = need_pkgs("sf","httr"),
    SOILGRIDS_10K     = if (file.exists(file.path(PROJECT_ROOT, "SoilGrids", "US.SOL"))) "" else "missing SoilGrids/US.SOL",
    SOILGRIDS_ONLINE  = need_pkgs("httr","sf"),       # REST mode; VRT mode also needs python rasterio
    HWSD              = paste(Filter(nzchar, c(need_pkgs("terra","RSQLite"), if (!HAS_HWSD) "missing HWSD raster/db files" else "")), collapse = "; "),
    "")
  list(ok = !nzchar(r), reason = r)
}
all_srcs <- sort(unique(unlist(lapply(metas, function(m) {
  s <- select_sources(region_of(m$lat, m$lon), m$lat, m$lon, m$planting_year); c(s$weather, s$soil) }))))
st <- lapply(all_srcs, source_status); names(st) <- all_srcs
RUNNABLE <- all_srcs[vapply(st, `[[`, logical(1), "ok")]

cat("\n=== Source preflight (this machine) ===\n")
pf <- data.frame(source = all_srcs,
                 runnable = ifelse(vapply(st, `[[`, logical(1), "ok"), "yes", "NO"),
                 note = vapply(st, `[[`, character(1), "reason"), stringsAsFactors = FALSE)
print(pf, row.names = FALSE, right = FALSE)
write.csv(pf, file.path(OUT_DIR, "source_preflight.csv"), row.names = FALSE)
if (!PRUNE_UNRUNNABLE) RUNNABLE <- all_srcs    # report only, don't drop
prune <- function(s) {
  if (is.null(s)) return(NULL)
  w <- intersect(s$weather, RUNNABLE); so <- intersect(s$soil, RUNNABLE)
  if (!length(w) || !length(so)) return(NULL)
  list(weather = w, soil = so)
}
if (CHECK_ONLY) { cat("\n[--check] preflight only; no runs executed.\n"); quit(save = "no") }

# ---- Drive all experiments ------------------------------------------------
cat("\n=== Validation plan ===\n")
plan <- do.call(rbind, lapply(metas, function(m) {
  s <- prune(select_sources(region_of(m$lat, m$lon), m$lat, m$lon, m$planting_year))
  data.frame(exp = m$exp, lat = m$lat, lon = m$lon, region = region_of(m$lat, m$lon) %||% NA,
             year = m$planting_year, n_trt = m$n_treatments,
             weather = if (is.null(s)) NA else paste(s$weather, collapse = "+"),
             soil = if (is.null(s)) NA else paste(s$soil, collapse = "+"),
             stringsAsFactors = FALSE)
}))
print(plan, row.names = FALSE)
write.csv(plan, file.path(OUT_DIR, "validation_plan.csv"), row.names = FALSE)
if (DRY_RUN) { cat("\n[dry-run] metadata resolved; no runs executed.\n"); quit(save = "no") }


long <- list()
for (m in metas) {
  cat(sprintf("\n----- %s  (%.3f, %.3f) -----\n", m$exp, m$lat %||% NA, m$lon %||% NA))
  reg <- region_of(m$lat, m$lon); srcs <- prune(select_sources(reg, m$lat, m$lon, m$planting_year))
  if (is.null(srcs)) { cat("  no runnable coordinates/sources — skipping.\n"); next }

  obs <- tryCatch(parse_filea(file.path(CROP_DIR, paste0(m$exp, ".", CROP_EXT, "A"))), error = function(e) NULL)
  if (!is.null(obs)) long[[length(long)+1]] <- data.frame(exp = m$exp, treatment = obs$treatment,
        pathway = "observed", source = "observed", HWAM = obs$HWAM_obs)

  cat("  [1/2] original DSSAT (local weather/soil)...\n")
  orig <- run_original(m)
  if (!is.null(orig)) long[[length(long)+1]] <- data.frame(exp = m$exp, treatment = orig$treatment,
        pathway = "original", source = "DSSAT-local", HWAM = orig$HWAM_orig)

  cat(sprintf("  [2/2] pipeline sweep (%s x %s)...\n", paste(srcs$weather, collapse=","), paste(srcs$soil, collapse=",")))
  pipe <- run_pipeline_sweep(m, srcs)
  if (!is.null(pipe)) long[[length(long)+1]] <- data.frame(exp = m$exp, treatment = pipe$treatment,
        pathway = "pipeline", source = paste(pipe$weather_source, pipe$soil_source, sep = " | "), HWAM = pipe$HWAM_pipe)
}

if (!length(long)) { cat("\nNo results produced.\n"); quit(save = "no", status = 1) }
LONG <- do.call(rbind, long)
write.csv(LONG, file.path(OUT_DIR, "validation_long.csv"), row.names = FALSE)
cat(sprintf("\nLong table (%d rows) -> validation/validation_long.csv\n", nrow(LONG)))

# ---- Metrics: every simulated source vs observed --------------------------
obs_tab <- LONG[LONG$pathway == "observed", c("exp","treatment","HWAM")]; names(obs_tab)[3] <- "obs"
sims    <- LONG[LONG$pathway != "observed", ]
met <- do.call(rbind, lapply(split(sims, list(sims$exp, sims$source), drop = TRUE), function(g) {
  mg <- merge(g, obs_tab, by = c("exp","treatment"))
  if (!nrow(mg)) return(NULL)
  cbind(data.frame(exp = g$exp[1], pathway = g$pathway[1], source = g$source[1]),
        metrics(mg$obs, mg$HWAM))
}))
met <- met[order(met$exp, -met$d), ]
write.csv(met, file.path(OUT_DIR, "validation_metrics.csv"), row.names = FALSE)
cat("Metrics -> validation/validation_metrics.csv\n\n"); print(met, row.names = FALSE)

# ---- Figures --------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressMessages(library(ggplot2))
  sc <- merge(sims, obs_tab, by = c("exp","treatment"))
  if (nrow(sc)) {
    rng <- range(c(sc$obs, sc$HWAM), na.rm = TRUE)
    p1 <- ggplot(sc, aes(obs, HWAM, colour = pathway, shape = pathway)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
      geom_point(size = 2, alpha = 0.8) + facet_wrap(~exp, scales = "free") +
      labs(title = "Observed vs. simulated grain yield (HWAM)",
        subtitle = "dashed = 1:1 line; each point a treatment",
        x = "Observed HWAM (kg/ha)", y = "Simulated HWAM (kg/ha)") + theme_minimal(base_size = 11)
    ggsave(file.path(OUT_DIR, "fig_obs_vs_sim.png"), p1, width = 10, height = 7, dpi = 150)
    cat("Figure -> validation/fig_obs_vs_sim.png\n")
  }
  if (nrow(met)) {
    p2 <- ggplot(met, aes(x = reorder(source, RMSE), y = RMSE, fill = pathway)) +
      geom_col() + facet_wrap(~exp, scales = "free") + coord_flip() +
      labs(title = "Yield error by input source (RMSE vs observed)", x = NULL, y = "RMSE (kg/ha)") +
      theme_minimal(base_size = 11)
    ggsave(file.path(OUT_DIR, "fig_metrics_rmse.png"), p2, width = 10, height = 7, dpi = 150)
    cat("Figure -> validation/fig_metrics_rmse.png\n")
  }
} else cat("(ggplot2 not installed — figures skipped)\n")

cat("\nDone.\n")
