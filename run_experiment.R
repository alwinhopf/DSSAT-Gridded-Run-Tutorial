# =============================================================================
# run_experiment.R — orchestrate a WEATHER x SOIL [x PERIOD] sensitivity study.
#
# WHAT IT DOES
#   Reads experiment.yml, builds the full factorial of weather_source x
#   soil_source x period (period optional, minus any exclusions), and runs the
#   EXISTING pipeline (dssat_main_pipeline.R) once per combination. Afterwards it
#   stitches every combination's results into one tidy table (tagged with
#   weather_source / soil_source / period), and for each response variable writes
#   a per-combination summary, an ANOVA variance decomposition, and a boxplot —
#   so you can see how much each input choice moves the outcome.
#   For richer post-hoc analysis (spatial sensitivity maps, rank stability, RMSD),
#   run analyze_experiment.R on the combined CSV produced by this script.
#
# HOW TO RUN
#       Rscript run_experiment.R                 # uses experiment.yml
#       Rscript run_experiment.R my_other.yml    # use a different experiment file
#   or in RStudio: open this file and click Source.
#
# INJECTION (non-invasive, parallel-safe)
#   The pipeline reads its settings from config.yml — OR, if the environment
#   variable DSSAT_CONFIG_FILE is set, from that file (see config_loader.R/.py).
#   For each combination this script writes a private merged config to a temp
#   file and runs the pipeline as a fresh Rscript subprocess with
#   DSSAT_CONFIG_FILE pointing at it. Your real config.yml is NEVER modified, and
#   parallel workers never share a config — so nothing can clobber anything.
#
# PARALLELISM (options.max_parallel > 1)
#   Number of combinations to run concurrently. 1 = serial. >1 uses:
#   - parallel::mclapply (fork) on Linux/macOS
#   - parallel::makeCluster + parLapplyLB (PSOCK sockets) on Windows
#   All platforms get true multi-core parallelism. The orchestrator first
#   runs a serial "warm-up" pass that touches each weather/soil source once
#   to populate the shared download caches, then parallelises the rest.
#   Each combination gets a unique project_name so its weather/soil/run
#   folders never collide with a concurrent worker's.
# =============================================================================

suppressWarnings(suppressMessages({
  if (!requireNamespace("yaml", quietly = TRUE))
    install.packages("yaml", repos = "https://cloud.r-project.org")
}))

# --- 0. Locate project root & experiment file ------------------------------
get_script_dir <- function() {
  if (requireNamespace("this.path", quietly = TRUE)) {
    d <- tryCatch(dirname(this.path::this.path()), error = function(e) NA)
    if (!is.na(d)) return(d)
  }
  getwd()
}
PROJECT_ROOT <- normalizePath(get_script_dir(), mustWork = FALSE)
setwd(PROJECT_ROOT)

args <- commandArgs(trailingOnly = TRUE)
EXPERIMENT_FILE <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "experiment.yml"
if (!file.exists(EXPERIMENT_FILE))
  stop(sprintf("Experiment file not found: %s (looked in %s)", EXPERIMENT_FILE, PROJECT_ROOT))

# The pipeline must understand DSSAT_CONFIG_FILE, or per-combination configs are
# silently ignored and every run would read the shared config.yml instead.
loader_path <- file.path(PROJECT_ROOT, "config_loader.R")
if (file.exists(loader_path) &&
    !any(grepl("DSSAT_CONFIG_FILE", readLines(loader_path, warn = FALSE)))) {
  stop("config_loader.R does not support DSSAT_CONFIG_FILE. Update it (this repo's ",
       "version does) before running parallel/period experiments.")
}

exp_cfg   <- yaml::read_yaml(EXPERIMENT_FILE)
base_cfg  <- if (is.null(exp_cfg$base)) list() else exp_cfg$base
factors   <- exp_cfg$factors
opts      <- if (is.null(exp_cfg$options)) list() else exp_cfg$options
out_cfg   <- if (is.null(exp_cfg$output))  list() else exp_cfg$output

STOP_ON_ERROR  <- isTRUE(opts$stop_on_error)
REUSE_EXISTING <- isTRUE(opts$reuse_existing)
DRY_RUN        <- isTRUE(opts$dry_run)
VALIDATE       <- if (is.null(opts$validate)) TRUE else isTRUE(opts$validate)
MAX_PARALLEL   <- max(1L, as.integer(if (is.null(opts$max_parallel)) 1 else opts$max_parallel))
# Note: Windows does not support fork(), but we use makeCluster/parLapply on
# that platform so MAX_PARALLEL is honoured on all operating systems.
EXP_NAME       <- if (is.null(exp_cfg$experiment_name)) "experiment" else exp_cfg$experiment_name

# Response variable(s): accept response_vars (list) or response_var (single).
RESPONSE_VARS <- if (!is.null(opts$response_vars)) {
  as.character(unlist(opts$response_vars))
} else if (!is.null(opts$response_var)) {
  as.character(opts$response_var)
} else "final_grain_kg_ha"

weather_list <- as.character(unlist(factors$weather_source))
soil_list    <- as.character(unlist(factors$soil_source))
if (length(weather_list) == 0 || length(soil_list) == 0)
  stop("experiment.yml must define factors$weather_source and factors$soil_source (each a non-empty list).")

sanitize <- function(x) gsub("[^A-Za-z0-9_\\-]", "_", x)

# --- 1. Periods ------------------------------------------------------------
# Each "START-END" entry overrides the weather years. With no period factor we
# use one window taken from base/config (label kept for the output table).
parse_period <- function(p) {
  m <- regmatches(p, regexec("^\\s*(\\d{4})\\s*-\\s*(\\d{4})\\s*$", p))[[1]]
  if (length(m) != 3) stop(sprintf("Bad period '%s' — expected \"YYYY-YYYY\".", p))
  list(label = paste0(m[2], "_", m[3]), start = as.integer(m[2]), end = as.integer(m[3]))
}
base_y0 <- if (!is.null(base_cfg$weather_start_year)) base_cfg$weather_start_year else 1982
base_y1 <- if (!is.null(base_cfg$weather_end_year))   base_cfg$weather_end_year   else 1983
if (!is.null(factors$period)) {
  periods <- lapply(as.character(unlist(factors$period)), parse_period)
} else {
  periods <- list(list(label = paste0(base_y0, "_", base_y1), start = base_y0, end = base_y1))
}
PERIODS_SWEPT <- length(periods) > 1
PARALLEL_ISO  <- MAX_PARALLEL > 1L

# --- 2. Build the factorial (minus exclusions) & predict output paths ------
grid <- expand.grid(weather_source = weather_list, soil_source = soil_list,
                    period_idx = seq_along(periods),
                    stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

excluded <- function(w, s) {
  ex <- exp_cfg$exclude
  if (is.null(ex)) return(FALSE)
  for (pair in ex)
    if (!is.null(pair$weather_source) && !is.null(pair$soil_source) &&
        pair$weather_source == w && pair$soil_source == s) return(TRUE)
  FALSE
}
grid <- grid[!mapply(excluded, grid$weather_source, grid$soil_source), , drop = FALSE]
rownames(grid) <- NULL

base_project <- if (!is.null(base_cfg$project_name)) base_cfg$project_name else "dssat_spatial_demo"
grid_m       <- if (!is.null(base_cfg$grid_spacing_meters)) base_cfg$grid_spacing_meters else 40000
res_tag      <- if (grid_m < 1000) paste0(grid_m, "m") else paste0(grid_m / 1000, "km")

combos <- data.frame(
  weather_source = grid$weather_source, soil_source = grid$soil_source,
  period_label = vapply(grid$period_idx, function(k) periods[[k]]$label, ""),
  year_start   = vapply(grid$period_idx, function(k) periods[[k]]$start, 0L),
  year_end     = vapply(grid$period_idx, function(k) periods[[k]]$end,   0L),
  stringsAsFactors = FALSE)

# Per-combination project_name. ALWAYS namespace by the year window — the
# pipeline's folder/results names don't encode years, so without this, changing
# the period would silently reuse stale weather and results. In parallel mode
# also append weather+soil so concurrent workers never share a weather/soil/run
# folder (in serial mode those folders are safely reused within a period).
combos$project_name <- sanitize(paste0(
  base_project, "_", combos$period_label,
  if (PARALLEL_ISO) paste0("_", combos$weather_source, "_", combos$soil_source) else ""))
combos$scenario_id  <- sanitize(paste0(combos$project_name, "_", res_tag, "_",
                                       combos$weather_source, "_", combos$soil_source))
combos$results_path <- file.path(PROJECT_ROOT, "results", paste0(combos$scenario_id, "_results.csv"))

# --- 3. Pre-flight validation (coverage guardrails) ------------------------
combos$valid <- TRUE
combos$skip_reason <- ""
soft_warnings <- character()
if (VALIDATE) {
  cdsrc <- file.exists(path.expand("~/.cdsapirc")) || nzchar(Sys.getenv("CDSAPI_RC")) ||
           nzchar(Sys.getenv("CDSAPI_KEY"))
  warned_us <- FALSE; warned_rest <- FALSE
  for (i in seq_len(nrow(combos))) {
    w <- combos$weather_source[i]; s <- combos$soil_source[i]
    # Hard-invalid: drop the combination.
    if (w == "AGERA5" && !cdsrc) {
      combos$valid[i] <- FALSE; combos$skip_reason[i] <- "AGERA5 needs a Copernicus CDS key (~/.cdsapirc)"
    } else if (w %in% c("NASA_POWER", "NASA_POWER_CHIRPS") && combos$year_start[i] < 1984) {
      combos$valid[i] <- FALSE
      combos$skip_reason[i] <- sprintf("%s has no data before 1984 (period starts %d)", w, combos$year_start[i])
    }
    # Soft warnings: keep the combination, flag a risk once.
    if (combos$valid[i] && (w == "GRIDMET" || s == "SSURGO") && !warned_us) {
      soft_warnings <- c(soft_warnings, "GRIDMET/SSURGO are US-only — ensure your region is within the contiguous US.")
      warned_us <- TRUE
    }
    if (combos$valid[i] && s == "SOILGRIDS_ONLINE" && MAX_PARALLEL > 1 && !warned_rest) {
      soft_warnings <- c(soft_warnings, "SOILGRIDS_ONLINE (REST) is rate-limited; high max_parallel may trigger throttling — consider soilgrids_mode: VRT.")
      warned_rest <- TRUE
    }
  }
}

# --- 4. Print the plan -----------------------------------------------------
cat(sprintf("\n=== Experiment: %s ===\n", EXP_NAME))
cat(sprintf("Region/base      : %s_%s  (%s)\n", base_project, res_tag,
            paste(unlist(base_cfg$state_name_filter), collapse = ", ")))
cat(sprintf("Weather sources  : %s\n", paste(weather_list, collapse = ", ")))
cat(sprintf("Soil sources     : %s\n", paste(soil_list, collapse = ", ")))
if (PERIODS_SWEPT)
  cat(sprintf("Periods          : %s\n", paste(vapply(periods, `[[`, "", "label"), collapse = ", ")))
cat(sprintf("Combinations     : %d valid / %d total\n", sum(combos$valid), nrow(combos)))
cat(sprintf("Response var(s)  : %s\n", paste(RESPONSE_VARS, collapse = ", ")))
cat(sprintf("Parallelism      : max_parallel = %d%s\n", MAX_PARALLEL,
            if (PARALLEL_ISO) "  (parallel: isolated per-combo folders)" else "  (serial)"))
if (length(soft_warnings)) { cat("\nWarnings:\n"); for (m in soft_warnings) cat("  ! ", m, "\n", sep = "") }
cat("\n")
for (i in seq_len(nrow(combos))) {
  lab <- if (PERIODS_SWEPT) sprintf("%-12s x %-14s x %s", combos$weather_source[i], combos$soil_source[i], combos$period_label[i])
         else sprintf("%-12s x %-14s", combos$weather_source[i], combos$soil_source[i])
  if (combos$valid[i]) cat(sprintf("  [%02d] %s -> %s\n", i, lab, basename(combos$results_path[i])))
  else                 cat(sprintf("  [--] %s -> SKIPPED: %s\n", lab, combos$skip_reason[i]))
}

if (DRY_RUN) { cat("\n[dry_run] Plan printed above. No DSSAT runs executed.\n"); quit(save = "no") }

run_idx <- which(combos$valid)
if (length(run_idx) == 0) { cat("\nNo valid combinations to run.\n"); quit(save = "no", status = 1) }

# --- 5. Merge config + run a single combination ----------------------------
existing_cfg <- if (file.exists("config.yml")) yaml::read_yaml("config.yml") else list()
merge_over   <- function(a, b) { for (k in names(b)) a[[k]] <- b[[k]]; a }
LOG_DIR      <- file.path(PROJECT_ROOT, "results", "experiment_logs")
CFG_DIR      <- file.path(PROJECT_ROOT, ".experiment_configs")
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CFG_DIR, recursive = TRUE, showWarnings = FALSE)

run_one <- function(i) {
  w <- combos$weather_source[i]; s <- combos$soil_source[i]
  tag <- sprintf("[%02d] %s x %s%s", i, w, s,
                 if (PERIODS_SWEPT) paste0(" x ", combos$period_label[i]) else "")
  log_file <- file.path(LOG_DIR, paste0(combos$scenario_id[i], ".log"))

  if (REUSE_EXISTING && file.exists(combos$results_path[i])) {
    cat(sprintf("%s  -> reusing existing results\n", tag)); return(TRUE)
  }

  # Private, fully-merged config for this combination.
  cfg <- merge_over(existing_cfg, base_cfg)
  cfg$project_name       <- combos$project_name[i]
  cfg$weather_source     <- w
  cfg$soil_source        <- s
  cfg$weather_start_year <- combos$year_start[i]
  cfg$weather_end_year   <- combos$year_end[i]
  cfg$run_tag            <- ""        # default run naming => DSSAT_RUN_NAME == scenario_id
  cfg$run_name_override  <- ""
  cfg$run_name_style     <- "grid"
  cfg_file <- file.path(CFG_DIR, paste0(combos$scenario_id[i], ".yml"))
  yaml::write_yaml(cfg, cfg_file)

  cat(sprintf("%s  -> running... (log: %s)\n", tag, basename(log_file)))
  t0 <- Sys.time()
  status <- system2("Rscript", args = shQuote("dssat_main_pipeline.R"),
                    env = paste0("DSSAT_CONFIG_FILE=", shQuote(cfg_file)),
                    stdout = log_file, stderr = log_file)
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

  ok <- (status == 0) && file.exists(combos$results_path[i])
  if (ok) {
    cat(sprintf("%s  -> OK (%s min)\n", tag, dt))
  } else {
    cat(sprintf("%s  -> FAILED (%s min) — see log: %s\n", tag, dt, log_file))
    if (STOP_ON_ERROR)
      stop(sprintf("Combination failed (%s); stop_on_error=true. See %s", tag, log_file))
  }
  ok
}

# --- 6. Execute: serial warm-up of each source, then parallel remainder ----
# Helper: run a list of combination indices in parallel.
# On Linux/macOS we fork with mclapply (no serialisation overhead).
# On Windows we use a PSOCK cluster so all cores are used (Windows has no fork).
run_parallel <- function(indices, n_workers) {
  if (length(indices) == 0L) return(list())
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(n_workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # Export every object that run_one() closes over, including run_one itself.
    parallel::clusterExport(cl, varlist = c(
      "combos", "REUSE_EXISTING", "PERIODS_SWEPT",
      "STOP_ON_ERROR", "existing_cfg", "base_cfg",
      "CFG_DIR", "LOG_DIR", "merge_over", "run_one"
    ), envir = globalenv())
    parallel::clusterEvalQ(cl, {
      suppressWarnings(suppressMessages({
        if (!requireNamespace("yaml", quietly = TRUE))
          install.packages("yaml", repos = "https://cloud.r-project.org")
        library(yaml)
      }))
    })
    res <- parallel::parLapplyLB(cl, indices, function(i) isTRUE(run_one(i)))
  } else {
    res <- parallel::mclapply(indices, function(i) isTRUE(run_one(i)),
                              mc.cores = n_workers, mc.preschedule = FALSE)
  }
  res
}

results_ok <- setNames(logical(length(run_idx)), run_idx)
if (MAX_PARALLEL <= 1L) {
  for (i in run_idx) results_ok[as.character(i)] <- isTRUE(run_one(i))
} else {
  # Warm-up = the minimal set of combinations that touches every weather (by
  # source+years) and every soil source once, so caches are populated serially.
  seen_w <- character(); seen_s <- character(); warm <- integer()
  for (i in run_idx) {
    wkey <- paste(combos$weather_source[i], combos$year_start[i], combos$year_end[i])
    skey <- combos$soil_source[i]
    if (!(wkey %in% seen_w) || !(skey %in% seen_s)) {
      warm <- c(warm, i); seen_w <- union(seen_w, wkey); seen_s <- union(seen_s, skey)
    }
  }
  rest <- setdiff(run_idx, warm)
  cat(sprintf("\n-- Warm-up pass (serial, %d combos) --\n", length(warm)))
  for (i in warm) results_ok[as.character(i)] <- isTRUE(run_one(i))
  if (length(rest)) {
    cat(sprintf("\n-- Parallel pass (%d combos, %d workers) [%s] --\n",
                length(rest), MAX_PARALLEL,
                if (.Platform$OS.type == "windows") "PSOCK/Windows" else "fork/Unix"))
    res <- run_parallel(rest, MAX_PARALLEL)
    for (k in seq_along(rest)) results_ok[as.character(rest[k])] <- isTRUE(res[[k]])
  }
}

# --- 7. Aggregate every combination's results into one tidy table ----------
read_tagged <- function(i) {
  p <- combos$results_path[i]
  if (!file.exists(p)) return(NULL)
  d <- tryCatch(read.csv(p, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d$weather_source <- combos$weather_source[i]
  d$soil_source    <- combos$soil_source[i]
  d$period         <- combos$period_label[i]
  d$scenario_id    <- combos$scenario_id[i]
  d
}
ok_idx   <- run_idx[results_ok[as.character(run_idx)]]
all_list <- Filter(Negate(is.null), lapply(ok_idx, read_tagged))

if (length(all_list) == 0) {
  cat("\nNo results were produced. Check the logs in results/experiment_logs/.\n")
  quit(save = "no", status = 1)
}
all_cols <- Reduce(union, lapply(all_list, names))
all_list <- lapply(all_list, function(d) { d[setdiff(all_cols, names(d))] <- NA; d[all_cols] })
combined <- do.call(rbind, all_list)

combined_csv <- if (is.null(out_cfg$combined_csv)) file.path("results", paste0("EXPERIMENT_", EXP_NAME, "_combined.csv")) else out_cfg$combined_csv
dir.create(dirname(combined_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(combined, combined_csv, row.names = FALSE)
cat(sprintf("\nCombined results (%d rows) -> %s\n", nrow(combined), combined_csv))

# Per-variable output paths: insert the variable name only when >1 variable, so
# single-variable runs keep the clean names from experiment.yml:output.
out_path <- function(key, default_suffix, var) {
  base <- if (is.null(out_cfg[[key]])) file.path("results", paste0("EXPERIMENT_", EXP_NAME, "_", default_suffix)) else out_cfg[[key]]
  if (length(RESPONSE_VARS) == 1) return(base)
  ext <- tools::file_ext(base); stem <- tools::file_path_sans_ext(base)
  paste0(stem, "_", var, ".", ext)
}

have_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
if (have_ggplot) suppressMessages(library(ggplot2))

# --- 8. Per response variable: summary + ANOVA variance + boxplot ----------
for (RV in RESPONSE_VARS) {
  cat(sprintf("\n----- Response variable: %s -----\n", RV))
  if (!RV %in% names(combined)) {
    cat(sprintf("  not found in results; skipping. Available numeric columns:\n  %s\n",
                paste(names(combined)[vapply(combined, is.numeric, logical(1))], collapse = ", ")))
    next
  }
  y <- suppressWarnings(as.numeric(combined[[RV]]))

  # Summary by combination (weather x soil x period).
  grp <- interaction(combined$weather_source, combined$soil_source, combined$period,
                     sep = " | ", drop = TRUE)
  summ <- do.call(rbind, lapply(split(seq_along(y), grp), function(idx) {
    yi <- y[idx]; yi <- yi[is.finite(yi)]
    data.frame(weather_source = combined$weather_source[idx[1]],
               soil_source = combined$soil_source[idx[1]], period = combined$period[idx[1]],
               n = length(yi), mean = if (length(yi)) mean(yi) else NA,
               sd = if (length(yi)) sd(yi) else NA, min = if (length(yi)) min(yi) else NA,
               max = if (length(yi)) max(yi) else NA, stringsAsFactors = FALSE)
  }))
  summ <- summ[order(-summ$mean), ]; rownames(summ) <- NULL
  if (!PERIODS_SWEPT) summ$period <- NULL
  sp <- out_path("summary_csv", "summary.csv", RV)
  write.csv(summ, sp, row.names = FALSE)
  cat(sprintf("  summary  -> %s\n", sp)); print(summ, digits = 5)

  # ANOVA variance decomposition: share of variance attributable to each factor.
  ok <- is.finite(y)
  terms <- c()
  if (length(unique(combined$weather_source[ok])) > 1) terms <- c(terms, "weather_source")
  if (length(unique(combined$soil_source[ok]))    > 1) terms <- c(terms, "soil_source")
  if (PERIODS_SWEPT && length(unique(combined$period[ok])) > 1) terms <- c(terms, "period")
  if (sum(ok) > length(terms) + 1 && length(terms) >= 1) {
    df <- combined[ok, ]; df$.y <- y[ok]
    fit <- tryCatch(aov(as.formula(paste(".y ~", paste(terms, collapse = " + "))), data = df),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      a <- summary(fit)[[1]]; ss <- a[["Sum Sq"]]
      labs <- sub("Residuals", "residual", trimws(rownames(a)))
      vd <- data.frame(factor = labs, sum_sq = ss,
                       variance_pct = round(100 * ss / sum(ss), 2),
                       response_var = RV, stringsAsFactors = FALSE)
      vp <- out_path("variance_csv", "variance.csv", RV)
      write.csv(vd, vp, row.names = FALSE)
      cat(sprintf("  variance -> %s\n", vp))
      cat(sprintf("  share of variance in %s explained:\n", RV))
      for (k in seq_len(nrow(vd))) cat(sprintf("    %-14s %6.2f%%\n", vd$factor[k], vd$variance_pct[k]))
    }
  } else cat("  (need >1 level on >=1 factor for ANOVA; skipped)\n")

  # Boxplot.
  if (have_ggplot) {
    pdf_df <- combined[ok, ]; pdf_df$.y <- y[ok]
    p <- ggplot(pdf_df, aes(x = weather_source, y = .y, fill = soil_source)) +
      geom_boxplot(position = position_dodge(0.8), outlier.size = 0.6, alpha = 0.85) +
      labs(title = sprintf("Impact of inputs on %s", RV),
           subtitle = sprintf("%s — %s", EXP_NAME, paste(unlist(base_cfg$state_name_filter), collapse = ", ")),
           x = "Weather source", y = RV, fill = "Soil source") +
      theme_minimal(base_size = 12)
    if (PERIODS_SWEPT) p <- p + facet_wrap(~period)
    pp <- out_path("plot_png", "boxplot.png", RV)
    dir.create(dirname(pp), recursive = TRUE, showWarnings = FALSE)
    ggsave(pp, p, width = 9, height = 5.5, dpi = 150)
    cat(sprintf("  boxplot  -> %s\n", pp))
  } else cat("  (ggplot2 not installed — skipped boxplot)\n")
}

cat(sprintf("\nDone. %d/%d valid combinations produced results.\n", length(ok_idx), length(run_idx)))
