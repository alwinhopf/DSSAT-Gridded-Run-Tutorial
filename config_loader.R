# config_loader.R
# ---------------------------------------------------------------------------
# Lightweight, dependency-tolerant loader for the shared config.yml.
#
# Sourced near the top of SECTION 0 in dssat_main_pipeline.R. Exposes one
# function, cfg_get(key, default), which returns the value from config.yml if
# present, otherwise the supplied default. If config.yml or the `yaml` package
# is missing, every call simply returns the default — so the pipeline still
# runs exactly as before.
# ---------------------------------------------------------------------------

.CONFIG <- list()

local({
  # Highest priority: an explicit path in the DSSAT_CONFIG_FILE env var. This
  # lets an orchestrator (e.g. run_experiment.R) point each parallel run at its
  # own private config file without touching the shared config.yml.
  env_path <- Sys.getenv("DSSAT_CONFIG_FILE", "")

  # Otherwise find config.yml: search this script's directory, then cwd.
  candidates <- c()
  this_dir <- tryCatch({
    if (requireNamespace("this.path", quietly = TRUE)) dirname(this.path::this.path()) else NA
  }, error = function(e) NA)
  if (!is.na(this_dir)) candidates <- c(candidates, file.path(this_dir, "config.yml"))
  candidates <- c(candidates, file.path(getwd(), "config.yml"), "config.yml")

  if (nzchar(env_path) && file.exists(env_path)) {
    path <- env_path
  } else {
    if (nzchar(env_path))
      message(sprintf("[config_loader] DSSAT_CONFIG_FILE set but not found (%s) — falling back to config.yml search.", env_path))
    path <- candidates[file.exists(candidates)][1]
  }
  if (is.na(path) || length(path) == 0) {
    message("[config_loader] No config.yml found — using in-script defaults.")
    return(invisible(NULL))
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    message("[config_loader] Package 'yaml' not installed — using in-script defaults.")
    return(invisible(NULL))
  }
  cfg <- tryCatch(yaml::read_yaml(path), error = function(e) {
    message(sprintf("[config_loader] Failed to parse %s (%s) — using defaults.", path, e$message))
    NULL
  })
  if (!is.null(cfg)) {
    assign(".CONFIG", cfg, envir = .GlobalEnv)
    message(sprintf("[config_loader] Loaded %d settings from %s", length(cfg), path))
  }
})

# cfg_get("key", default):
#   - returns default when key is absent or blank ("")
#   - coerces YAML lists to plain vectors (e.g. state_name_filter)
cfg_get <- function(key, default) {
  if (!exists(".CONFIG", envir = .GlobalEnv)) return(default)
  cfg <- get(".CONFIG", envir = .GlobalEnv)
  if (is.null(cfg[[key]])) return(default)
  val <- cfg[[key]]
  if (is.list(val)) val <- unlist(val, use.names = FALSE)
  if (is.character(val) && length(val) == 1 && !nzchar(val)) return(default)  # blank => default
  val
}
