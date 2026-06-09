# setup_renv.R
# ---------------------------------------------------------------------------
# Reproducible R environment for the DSSAT gridded pipeline.
#
# Run ONCE from the repo root (RStudio: Source; or `Rscript setup_renv.R`).
# It initializes renv, installs the pinned package set, and writes renv.lock,
# which you then COMMIT so collaborators get identical versions via
# `renv::restore()`.
#
# Why a generator instead of a committed lock? An renv.lock records exact
# versions AND content hashes resolved against your repositories/platform, so
# it must be produced on a real machine. This script makes that one command.
# ---------------------------------------------------------------------------

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# Use CRAN cloud to get precompiled Windows binaries for R 4.6.0.
options(repos = c(CRAN = "https://cloud.r-project.org"))

renv::init(bare = TRUE, restart = FALSE)

pkgs <- c(
  # Core
  "sf", "dplyr", "tidyr", "stringr", "lubridate",
  "foreach", "doParallel", "zoo", "R.utils", "processx",
  "ggplot2", "readr", "tibble", "rstudioapi", "this.path",
  # DSSAT interface
  "DSSAT",
  # Soil
  "soilDB", "terra", "httr", "jsonlite",
  # Weather
  "daymetr", "nasapower", "ncdf4",
  # Config + progress
  "yaml", "pbapply"
)

renv::install(pkgs)
renv::snapshot(packages = pkgs, prompt = FALSE)

# Find the script's directory and construct a relative path to dssatutils
.script_dir <- function() {
  argv <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", argv, value = TRUE)
  if (length(file_arg)) {
    path <- sub("^--file=", "", file_arg[1])
    return(dirname(normalizePath(path, mustWork = FALSE)))
  }
  for (i in seq_len(sys.nframe())) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile) && nzchar(ofile)) {
      return(dirname(normalizePath(ofile, mustWork = FALSE)))
    }
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) {
      return(dirname(normalizePath(ctx$path, mustWork = FALSE)))
    }
  }
  getwd()
}

renv::install("git::https://github.com/alwinhopf/dssatutils.git@v0.1.0")
renv::install("../dssatengine")


