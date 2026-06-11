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

# --- Shared Git packages: dssatutils + dssatengine -------------------------
# These live in their own GitHub repos and are reused across projects, so we
# install them FROM GITHUB (not a local "../path") and let renv pin the exact
# commit in renv.lock via the snapshot below. On a fresh machine
# `renv::restore()` then reinstalls the identical commit automatically — no
# devtools::install_local(), and no "did I reinstall after editing?" version
# skew.
#
# Pin style: "@main" installs the current tip and the resolved SHA is frozen
# into renv.lock by the snapshot. Once you cut release tags you can swap
# "@main" for "@<tag>" (e.g. "@v0.1.0") for human-readable pins.
#
# PREREQUISITES (one-time, on the machine that GENERATES the lock):
#   1. Commit + push any local edits to dssatutils / dssatengine FIRST,
#      otherwise the GitHub install will not include them.
#   2. If either repo is private, configure a GitHub PAT so renv can clone it:
#        usethis::create_github_token(); gitcreds::gitcreds_set()
renv::install(c(
  "git::https://github.com/alwinhopf/dssatutils.git@main",
  "git::https://github.com/alwinhopf/dssatengine.git@main"
))

# Snapshot the FULL project library (CRAN packages + both Git packages) so
# renv.lock records every dependency, including the Git RemoteSha values that
# make dssatutils/dssatengine reproducible. COMMIT renv.lock after this runs.
renv::snapshot(type = "all", prompt = FALSE)

message("\nrenv.lock written. Commit it so a fresh machine can `renv::restore()`.")


