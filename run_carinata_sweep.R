# =============================================================================
# run_carinata_sweep.R — weather x soil x RESOLUTION sweep for carinata.
#
# Runs dssat_main_pipeline.R once per (weather, soil, grid-resolution) combination
# over the southern-US region, parses the carinata DSSAT Summary.OUT files, tags
# every row with weather_source / soil_source / resolution_km, and assembles one
# combined CSV — the input to compare_model_inputs.py (which now treats grid
# resolution as a third axis of model-input uncertainty).
#
#   Field yields come from DSSAT (BC crop, DSSAT48 executable -> Summary.OUT),
#   so this parses Summary.OUT directly rather than the pipeline's summary.csv.
#
# USAGE
#   DSSAT_EXE=/path/to/DSSAT48/dscsm048 Rscript run_carinata_sweep.R
# Resumable: combos whose run folders already hold Summary.OUT are re-parsed,
# not re-simulated.
# =============================================================================
suppressWarnings(suppressMessages({ library(DSSAT); library(sf); library(yaml) }))

# ---- Configuration --------------------------------------------------------
PROJECT        <- "carinata_sweep"
STATES         <- c("Texas","Louisiana","Mississippi","Alabama","Georgia","Florida","South Carolina")
RESOLUTIONS_KM <- c(100, 50)                       # test set; full study adds 25, 10, 5
WEATHER        <- c("DAYMET", "NASA_POWER")
SOIL           <- c("SSURGO", "SOILGRIDS_10K")
TEMPLATE       <- "UFJA1803.BCX"                   # carinata N-rate FileX
YEAR_START     <- 2018; YEAR_END <- 2019
TREATMENTS     <- c(4, 4)                           # single production N rate (134 kg N/ha) — fast test
DSSAT_EXE      <- Sys.getenv("DSSAT_EXE", unset = file.path("..","DSSAT48","dscsm048"))
OUT_CSV        <- file.path("results", "EXPERIMENT_carinata_resolution_combined.csv")
PIPELINE       <- "dssat_main_pipeline.R"
dir.create("results", showWarnings = FALSE)

# ---- Parse one run's Summary.OUT files into tagged rows --------------------
parse_run <- function(run_dir, gridpts_shp, weather, soil, res_km) {
  if (!dir.exists(run_dir) || !file.exists(gridpts_shp)) return(NULL)
  pts <- st_read(gridpts_shp, quiet = TRUE)
  cc  <- st_coordinates(st_transform(pts, 4326)); idc <- intersect(c("ID","id"), names(pts))[1]
  coords <- data.frame(point_id = sprintf("%08d", as.integer(pts[[idc]])), latitude = cc[,2], longitude = cc[,1])
  rows <- list()
  for (f in list.dirs(run_dir, recursive = FALSE)) {
    sm <- file.path(f, "Summary.OUT"); if (!file.exists(sm)) next
    d <- tryCatch(read_output(sm), error = function(e) NULL); if (is.null(d) || nrow(d) == 0) next
    hy <- if ("HYEAR" %in% names(d)) d$HYEAR else if ("HDAT" %in% names(d)) as.integer(substr(d$HDAT,1,4)) else NA
    rows[[basename(f)]] <- data.frame(point_id = basename(f), treatment = d$TRNO, crop_code = "BC",
      year_harvest = hy, harvest_date = if ("HDAT" %in% names(d)) d$HDAT else NA,
      final_grain_kg_ha = d$HWAM, top_weight_kg_ha = if ("CWAM" %in% names(d)) d$CWAM else NA,
      inorganic_n_applied_kg_ha = if ("NICM" %in% names(d)) d$NICM else NA,
      nitrate_leaching_kg_ha = if ("NLCM" %in% names(d)) d$NLCM else NA,
      cumulative_net_co2_emissions_kg_CO2_ha = if ("CO2EM" %in% names(d)) d$CO2EM else NA,
      cumulative_n2o_emissions_kg_N_ha = if ("N2OEM" %in% names(d)) d$N2OEM else NA,
      final_irrigation_amount_mm = 0)
  }
  if (!length(rows)) return(NULL)
  res <- merge(coords, do.call(rbind, rows), by = "point_id")
  res <- res[!is.na(res$final_grain_kg_ha) & res$final_grain_kg_ha >= 0, ]
  if (!nrow(res)) return(NULL)
  res$weather_source <- weather; res$soil_source <- soil; res$resolution_km <- res_km
  res
}

# ---- Sweep ----------------------------------------------------------------
all_rows <- list(); n_ok <- 0
for (res_km in RESOLUTIONS_KM) for (w in WEATHER) for (s in SOIL) {
  tag <- sprintf("%dkm | %s x %s", res_km, w, s)
  scenario <- sprintf("%s_%dkm_%s_%s", PROJECT, res_km, w, s)
  run_dir  <- file.path("dssat_runs", scenario)
  gridpts  <- file.path("gridpoints", sprintf("%s_%dkm.shp", PROJECT, res_km))

  have <- length(list.files(run_dir, "Summary.OUT", recursive = TRUE)) > 0
  if (!have) {
    cat(sprintf("[run ] %s ...\n", tag)); t0 <- Sys.time()
    cfg <- list(project_name = PROJECT, grid_spacing_meters = res_km * 1000, crop_extension = "BC",
                use_existing_point_shapefile = FALSE, boundary_shapefile_name = "tl_2024_us_state.shp",
                enable_boundary_filter = TRUE, boundary_filter_column = "NAME", state_name_filter = STATES,
                weather_source = w, weather_start_year = YEAR_START, weather_end_year = YEAR_END,
                soil_source = s, template_file_name = TEMPLATE, run_mode = "experiment",
                treatment_start = TREATMENTS[1], treatment_end = TREATMENTS[2], bundle_genotype_files = FALSE)
    cfgf <- file.path(tempdir(), paste0(scenario, ".yml")); write_yaml(cfg, cfgf)
    logf <- file.path("results", paste0(scenario, ".log"))
    st <- system2("Rscript", shQuote(PIPELINE),
                  env = c(paste0("DSSAT_CONFIG_FILE=", cfgf), paste0("DSSAT_EXE=", DSSAT_EXE)),
                  stdout = logf, stderr = logf)
    cat(sprintf("       (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  } else cat(sprintf("[skip] %s (Summary.OUT present)\n", tag))

  r <- parse_run(run_dir, gridpts, w, s, res_km)
  if (!is.null(r)) { all_rows[[scenario]] <- r; n_ok <- n_ok + 1
    cat(sprintf("       parsed %d points\n", length(unique(r$point_id)))) }
  else cat(sprintf("       NO results parsed for %s\n", tag))
}

if (!length(all_rows)) stop("No combinations produced results.")
combined <- do.call(rbind, all_rows)
write.csv(combined, OUT_CSV, row.names = FALSE)
cat(sprintf("\n%d/%d combinations OK -> %s (%d rows)\n",
            n_ok, length(RESOLUTIONS_KM)*length(WEATHER)*length(SOIL), OUT_CSV, nrow(combined)))
cat("points by resolution x source:\n")
print(table(combined$resolution_km, paste(combined$weather_source, combined$soil_source)))
