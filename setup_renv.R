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

# Use CRAN cloud for precompiled Windows/macOS binaries, and force binary
# installs so system-library packages (openssl, sf, terra, curl) are NEVER built
# from source — source builds need dev headers/Rtools that a fresh laptop lacks
# and are the usual cause of the openssl -> httr -> daymetr install cascade.
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (.Platform$OS.type == "windows" || Sys.info()[["sysname"]] == "Darwin") {
  options(pkgType = "binary", install.packages.check.source = "no")
}

renv::init(bare = TRUE, restart = FALSE)

pkgs <- c(
  # Core
  "sf", "dplyr", "tidyr", "stringr", "lubridate",
  "foreach", "doParallel", "zoo", "R.utils", "processx",
  "ggplot2", "readr", "tibble", "rstudioapi", "this.path",
  # DSSAT interface
  "DSSAT",
  # Soil (SSURGO/gNATSGO via soilDB; iSDAsoil/HWSD rasters via terra; SoilGrids)
  "soilDB", "terra", "httr", "jsonlite",
  # Soil attribute DB readers — HWSD: SQLite (DBI/RSQLite) + FAO Access .mdb
  # (odbc, needs the OS Microsoft Access ODBC driver); LUCAS: optional .xlsx.
  "DBI", "RSQLite", "odbc", "readxl",
  # Weather (Daymet/NASA-POWER point APIs; DWD stations; E-OBS/Xavier/CMFD NetCDF)
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
# We use the "user/repo" shorthand (NOT a "git::https://...@main" URL) on
# purpose: renv resolves this through the GitHub API and downloads a tarball, so
# it works on a bare laptop with NO local git client installed. The resolved SHA
# is still frozen into renv.lock by the snapshot below. The repos are public, so
# no PAT is required; GitHub's unauthenticated API limit (60 req/hr) is plenty
# for two installs. If you later make them private, set a PAT once first:
#   usethis::create_github_token(); gitcreds::gitcreds_set()
#
# PREREQUISITE: commit + push any local edits to dssatutils / dssatengine FIRST,
# otherwise the GitHub install won't include them.
renv::install(c(
  "alwinhopf/dssatutils",
  "alwinhopf/dssatengine"
))

# Snapshot the FULL project library (CRAN packages + both Git packages) so
# renv.lock records every dependency, including the Git RemoteSha values that
# make dssatutils/dssatengine reproducible. COMMIT renv.lock after this runs.
renv::snapshot(type = "all", prompt = FALSE)

message("\nrenv.lock written. Commit it so a fresh machine can `renv::restore()`.")


