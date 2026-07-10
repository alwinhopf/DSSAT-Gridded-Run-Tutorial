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
# These live in their own GitHub repos and are reused across projects. Fresh
# installs use release tags, not branch names, so the environment is reproducible
# without relying on whatever happens to be on `main`.
#
# For active shared-package development only, set USE_LOCAL_SHARED_PACKAGES=1 to
# install from sibling checkouts instead of the pinned tags below, then snapshot
# the resulting lockfile deliberately.
#
# If the repos are private, set a PAT once first:
#   usethis::create_github_token(); gitcreds::gitcreds_set()
workspace_root <- normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
shared_refs <- c(
  dssatutils = "v0.4.0",
  # Release-tag blocker: the dssatengine package is version 0.4.0, but the
  # remote currently has no v0.4.0 tag. Pin the verified release commit until
  # that tag is pushed.
  dssatengine = "84b6e50895e7e2e2a4b02553d2705f4d879d269b"
)
use_local_shared <- identical(Sys.getenv("USE_LOCAL_SHARED_PACKAGES"), "1")
for (pkg in names(shared_refs)) {
  local_dir <- file.path(workspace_root, pkg)
  target <- if (use_local_shared && file.exists(file.path(local_dir, "DESCRIPTION"))) {
    local_dir
  } else {
    paste0("alwinhopf/", pkg, "@", shared_refs[[pkg]])
  }
  message(sprintf("Installing shared package '%s' from %s", pkg, target))
  renv::install(target, prompt = FALSE)
}

# Snapshot the FULL project library (CRAN packages + both Git packages) so
# renv.lock records every dependency, including the Git RemoteSha values that
# make dssatutils/dssatengine reproducible. COMMIT renv.lock after this runs.
renv::snapshot(type = "all", prompt = FALSE)

message("\nrenv.lock written. Commit it so a fresh machine can `renv::restore()`.")
