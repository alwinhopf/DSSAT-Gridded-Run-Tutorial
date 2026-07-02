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
# from source; source builds need dev headers/Rtools that a fresh laptop lacks
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
  # Soil attribute DB readers: HWSD uses SQLite (DBI/RSQLite) + FAO Access .mdb
  # (odbc, needs the OS Microsoft Access ODBC driver); LUCAS: optional .xlsx.
  "DBI", "RSQLite", "odbc", "readxl",
  # Weather (Daymet/NASA-POWER point APIs; DWD stations; E-OBS/Xavier/CMFD NetCDF)
  "daymetr", "nasapower", "ncdf4",
  # Config + progress
  "yaml", "pbapply"
)

renv::install(pkgs)

# --- Shared Git packages: dssatutils + dssatengine -------------------------
# These live in their own GitHub repos and are reused across projects. During
# local development, prefer sibling checkouts if they exist; on a fresh machine,
# install from GitHub main and let renv pin the exact commit in renv.lock via the
# snapshot below.
#
# The GitHub fallback uses "user/repo" rather than a tag because the package
# versions can move ahead before release tags are cut. The resolved SHA is still
# frozen into renv.lock by the snapshot below. If GitHub prompts for
# authentication (private repo access or rate limits), set a PAT once first:
#   usethis::create_github_token(); gitcreds::gitcreds_set()
#
# PREREQUISITE for fresh machines: commit + push any local edits to dssatutils /
# dssatengine FIRST, otherwise the GitHub install won't include them.
workspace_root <- normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
for (pkg in c("dssatutils", "dssatengine")) {
  local_dir <- file.path(workspace_root, pkg)
  target <- if (file.exists(file.path(local_dir, "DESCRIPTION"))) {
    local_dir
  } else {
    paste0("alwinhopf/", pkg)
  }
  message(sprintf("Installing shared package '%s' from %s", pkg, target))
  renv::install(target, prompt = FALSE)
}

# Snapshot the FULL project library (CRAN packages + both Git packages) so
# renv.lock records every dependency, including the Git RemoteSha values that
# make dssatutils/dssatengine reproducible. COMMIT renv.lock after this runs.
renv::snapshot(type = "all", prompt = FALSE)

message("\nrenv.lock written. Commit it so a fresh machine can `renv::restore()`.")
