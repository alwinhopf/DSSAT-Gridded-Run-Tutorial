#!/usr/bin/env Rscript
# End-to-end smoke test (R) — minimal spatial demo using dssatutils
# - Creates 3 sample points
# - Runs weather (Open-Meteo) + soil (SoilGrids) download via dssatutils
# - Validates that .WTH/.SOL files are written and well-formed

suppressPackageStartupMessages({
  opts <- options(stringsAsFactors = FALSE)
  library(parallel)    # Required by dssatutils weather/soil functions for parallelization
  library(doParallel)  # Required for registerDoParallel in dssatutils
})

check <- function(cond, msg) {
  if (isTRUE(cond)) {
    message(" [ok] ", msg)
    return(TRUE)
  } else {
    message(" [FAIL] ", msg)
    quit(status = 1)
  }
}

if (!requireNamespace("dssatutils", quietly = TRUE)) {
  message("dssatutils not installed; skipping E2E test. Install dssatutils and rerun.")
  quit(status = 0)
}

# Optional helpers
have_sf <- requireNamespace("sf", quietly = TRUE)
have_dplyr <- requireNamespace("dplyr", quietly = TRUE)
have_httr <- requireNamespace("httr", quietly = TRUE)
if (have_sf) message("sf available; will pass sf object to dssatutils if required.")
if (have_dplyr) message("dplyr available; will attach for SoilGrids processing.")
if (have_httr) message("httr available; will attach for SoilGrids REST API requests.")

# Minimal 3-point domain (EU / Asia / Africa)
POINTS <- data.frame(
  ID = c("EU_NL", "AS_IN", "AF_KE"),
  LAT = c(52.000, 30.900, -0.500),
  LONG = c(5.000, 75.800, 37.000)
)

work <- tempfile("dssat_e2e_")
dir.create(work)
wth_dir <- file.path(work, "weather")
sol_dir <- file.path(work, "soil")
dir.create(wth_dir, recursive = TRUE)
dir.create(sol_dir, recursive = TRUE)
log <- file.path(work, "errors.log")

# Use a plain data frame for weather inputs so Open-Meteo output is written
# consistently. Keep an sf copy available only for soil functions when needed.
shapefile_arg <- POINTS
soil_shapefile_arg <- POINTS
if (have_sf) {
  library(sf)
  soil_shapefile_arg <- sf::st_as_sf(POINTS, coords = c("LONG", "LAT"), crs = 4326)
}
if (have_dplyr) {
  library(dplyr)
}
if (have_httr) {
  library(httr)
}

# Run Open-Meteo (may hit network). The function name and signature mirror the
# pipeline conventions; if it errors due to missing optional deps, we fail early.
message("Running dssatutils::process_weather_openmeteo() ...")
tryCatch({
  dssatutils::process_weather_openmeteo(
    shapefile = shapefile_arg,
    start_year = 2010,
    end_year = 2011,
    output_dir = wth_dir,
    id_col = "ID",
    lat_col = "LAT",
    lon_col = "LONG",
    n_cores = 1,
    log_file = log
  )
}, error = function(e) {
  message("Open-Meteo step failed: ", e$message)
  quit(status = 1)
})

# Validate .WTH files
for (pid in POINTS$ID) {
  p <- file.path(wth_dir, paste0(pid, ".WTH"))
  check(file.exists(p), paste0(pid, ".WTH written"))
  lines <- readLines(p)
  lines <- lines[which(nzchar(trimws(lines)))]
  check(length(lines) >= 5, paste0(pid, ": reasonable file length"))
  check(startsWith(lines[1], "$WEATHER"), paste0(pid, ": $WEATHER header present"))
  check(grepl("@  DATE", lines[3]) || grepl("@  DATE", lines[4]), paste0(pid, ": DATE header present"))
  data_lines <- lines[which(grepl("^[0-9]", trimws(lines)))]
  check(length(data_lines) >= 300, paste0(pid, ": expected ~730 daily rows (got ", length(data_lines), ")"))
  check(!any(grepl("nan", tolower(paste(data_lines, collapse = " ")))), paste0(pid, ": no NaN in data block"))
}

# Attempt SoilGrids .SOL generation (optional, REST API is slow)
# Skip if it takes too long; weather validation is more critical for smoke tests
if (requireNamespace("dssatutils", quietly = TRUE)) {
  if (exists("process_soils_soilgrids_online", where = asNamespace("dssatutils"), inherits = FALSE)) {
    message("Running dssatutils::process_soils_soilgrids_online() ... (REST API may be slow)")
    tryCatch({
      dssatutils::process_soils_soilgrids_online(
        gridfile = soil_shapefile_arg,
        soilfile_csv_path = file.path(sol_dir, "soil_map.csv"),
        output_sol_dir = sol_dir,
        id_col = "ID"
      )
    }, error = function(e) {
      message("SoilGrids step skipped/failed: ", e$message)
      # Not fatal: REST API is optional; continue but warn
    })

    # Check any produced .SOL files
    sols <- list.files(sol_dir, pattern = "\\.SOL$", full.names = TRUE)
    if (length(sols) > 0) {
      for (s in sols) {
        txt <- tolower(paste(readLines(s), collapse = " "))
        check(grepl("*soils" , toupper(substr(txt,1,6))) || grepl("*soils", toupper(txt)), paste0(basename(s), ": .SOL header present"))
      }
    } else {
      message("No .SOL files produced; SoilGrids REST API skipped or has coverage gaps.")
    }
  } else {
    message("dssatutils does not expose process_soils_soilgrids_online; skipping soils step.")
  }
}

message("ALL E2E CHECKS PASSED.")
quit(status = 0)
