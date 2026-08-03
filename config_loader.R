# config_loader.R
# ---------------------------------------------------------------------------
# Loader for the shared config.yml.
#
# Sourced near the top of SECTION 0 in dssat_main_pipeline.R. Exposes one
# config.yml beside this file is the canonical default mapping. A file named by
# DSSAT_CONFIG_FILE is merged over it as a partial study-specific override.
# This is the R twin of config_loader.py; keep precedence and blank handling
# aligned.
# ---------------------------------------------------------------------------

.CONFIG <- list()

local({
  deep_merge <- function(base, override) {
    if (!is.list(base)) base <- list()
    if (!is.list(override)) return(base)
    for (key in names(override)) {
      if (is.list(base[[key]]) && is.list(override[[key]])) {
        base[[key]] <- deep_merge(base[[key]], override[[key]])
      } else {
        base[[key]] <- override[[key]]
      }
    }
    base
  }

  validate_overlay <- function(base, override, prefix = "") {
    for (key in names(override)) {
      path <- if (nzchar(prefix)) paste(prefix, key, sep = ".") else key
      if (!(key %in% names(base))) stop("Unknown configuration key: ", path, call. = FALSE)
      if (is.list(override[[key]]) && is.list(base[[key]]) &&
          !is.null(names(override[[key]]))) {
        validate_overlay(base[[key]], override[[key]], path)
      } else if (is.list(override[[key]]) != is.list(base[[key]]) &&
                 length(override[[key]]) != 0L) {
        stop("Configuration type mismatch at ", path, call. = FALSE)
      }
    }
  }

  env_path <- Sys.getenv("DSSAT_CONFIG_FILE", "")
  this_dir <- tryCatch({
    if (requireNamespace("this.path", quietly = TRUE)) dirname(this.path::this.path()) else NA
  }, error = function(e) NA)
  if (is.na(this_dir)) this_dir <- getwd()
  base_path <- file.path(this_dir, "config.yml")
  if (!file.exists(base_path)) stop("Canonical pipeline config not found: ", base_path)
  if (nzchar(env_path) && !file.exists(env_path)) {
    stop("DSSAT_CONFIG_FILE points to a missing file: ", env_path)
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to load the central config.yml.")
  }
  cfg <- yaml::read_yaml(base_path)
  if (!is.list(cfg)) stop("Top level must be a YAML mapping: ", base_path)
  loaded <- base_path
  if (nzchar(env_path) && normalizePath(env_path) != normalizePath(base_path)) {
    override <- yaml::read_yaml(env_path)
    if (!is.list(override)) stop("Top level must be a YAML mapping: ", env_path)
    validate_overlay(cfg, override)
    cfg <- deep_merge(cfg, override)
    loaded <- c(loaded, env_path)
  }
  assign(".CONFIG", cfg, envir = .GlobalEnv)
  message(sprintf("[config_loader] Loaded %d settings from %s",
                  length(cfg), paste(loaded, collapse = " + ")))
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
