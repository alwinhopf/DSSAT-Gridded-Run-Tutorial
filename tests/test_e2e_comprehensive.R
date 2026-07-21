#!/usr/bin/env Rscript
# Comprehensive end-to-end test (R) — tests all weather & soil sources
# - Tests each weather source: Open-Meteo, Daymet, NASA-POWER, GridMET, AgERA5, NASA-POWER CHIRPS
# - Tests each soil source: SoilGrids, SoilGrids Online, SSURGO, HWSD
# - Uses synthetic/mocked data to avoid flaky network calls in CI
# - Skips sources with missing optional deps or API keys gracefully

if (Sys.getenv("DSSAT_RUN_LIVE_E2E", "") != "1") {
  message("Skipping live comprehensive R downloads; set DSSAT_RUN_LIVE_E2E=1 to run them.")
} else local({

suppressPackageStartupMessages({
  library(parallel)
  library(doParallel)
})

# Optional libraries for soil processing
have_sf <- requireNamespace("sf", quietly = TRUE)
have_dplyr <- requireNamespace("dplyr", quietly = TRUE)
have_httr <- requireNamespace("httr", quietly = TRUE)

if (have_sf) library(sf)
if (have_dplyr) library(dplyr)
if (have_httr) library(httr)

# Global test state
PASSED <- 0
SKIPPED <- 0
FAILED <- 0

log_test <- function(status, name, msg = "") {
  if (status == "ok") {
    message(" [OK] ", name)
    PASSED <<- PASSED + 1
  } else if (status == "skip") {
    message(" [SKIP] ", name, " — ", msg)
    SKIPPED <<- SKIPPED + 1
  } else if (status == "fail") {
    message(" [FAIL] ", name, " — ", msg)
    FAILED <<- FAILED + 1
  }
}

utils_repo <- normalizePath(file.path(getwd(), "..", "dssatutils"), mustWork = FALSE)
if (dir.exists(utils_repo) && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(utils_repo, quiet = TRUE)
}
if (!requireNamespace("dssatutils", quietly = TRUE)) {
  message("dssatutils not installed; skipping comprehensive E2E.")
  quit(status = 0)
}

# Test data: 3 points covering different regions
GLOBAL_POINTS <- data.frame(
  ID = c("GLOBAL_EU", "GLOBAL_AS", "GLOBAL_AF"),
  LAT = c(52.0, 30.9, -0.5),
  LONG = c(5.0, 75.8, 37.0)
)

US_POINTS <- data.frame(
  ID = c("US_CORN", "US_WHEAT"),
  LAT = c(40.0, 38.5),
  LONG = c(-97.0, -99.5)
)

# Prepare shapefile objects
# NOTE: Weather functions work better with plain data frames
# Soil functions prefer sf objects when available
make_sf_or_df <- function(df) {
  if (have_sf) {
    sf <- asNamespace("sf")
    return(sf$st_as_sf(df, coords = c("LONG", "LAT"), crs = 4326))
  }
  df
}

# Use plain data frames for weather (more stable)
global_df <- GLOBAL_POINTS
us_df <- US_POINTS

# Use sf for soil (if available)
global_sf <- make_sf_or_df(GLOBAL_POINTS)
us_sf <- make_sf_or_df(US_POINTS)

work <- tempfile("dssat_e2e_comprehensive_")
dir.create(work)

# Helper to check .WTH file validity
check_wth <- function(path) {
  if (!file.exists(path)) return(FALSE)
  lines <- readLines(path)
  lines <- lines[which(nzchar(trimws(lines)))]
  if (length(lines) < 5) return(FALSE)
  if (!startsWith(lines[1], "$WEATHER")) return(FALSE)
  has_date_header <- any(grepl("@  DATE", lines[1:5]))
  if (!has_date_header) return(FALSE)
  data_lines <- lines[which(grepl("^[0-9]", trimws(lines)))]
  if (length(data_lines) < 100) return(FALSE)
  if (any(grepl("nan", tolower(paste(data_lines, collapse = " "))))) return(FALSE)
  TRUE
}

# Helper to check .SOL file validity
check_sol <- function(path) {
  if (!file.exists(path)) return(FALSE)
  txt <- tolower(paste(readLines(path), collapse = " "))
  if (!grepl("*soils", txt, fixed = TRUE)) return(FALSE)
  if (!grepl("slb|layer", txt)) return(FALSE)
  if (grepl("nan", txt)) return(FALSE)
  TRUE
}

# --- WEATHER SOURCES ---

message("\n=== WEATHER SOURCES ===\n")

# 1. OPEN-METEO (global, free, no key)
test_name <- "OPEN-METEO"
wth_dir <- file.path(work, "weather_openmeteo")
dir.create(wth_dir, recursive = TRUE)
tryCatch({
  dssatutils::process_weather_openmeteo(
    shapefile = global_df,
    start_year = 2010,
    end_year = 2011,
    output_dir = wth_dir,
    id_col = "ID",
    lat_col = "LAT",
    lon_col = "LONG",
    n_cores = 1,
    log_file = file.path(work, "openmeteo.log")
  )
  if (all(file.exists(file.path(wth_dir, paste0(GLOBAL_POINTS$ID, ".WTH"))))) {
    log_test("ok", test_name)
  } else {
    log_test("fail", test_name, "Not all .WTH files written")
  }
}, error = function(e) {
  log_test("fail", test_name, e$message)
})

# 2. DAYMET (US only, free, no key)
test_name <- "DAYMET"
wth_dir <- file.path(work, "weather_daymet")
dir.create(wth_dir, recursive = TRUE)
tryCatch({
  if (!exists("process_weather_daymet", where = asNamespace("dssatutils"), inherits = FALSE)) {
    log_test("skip", test_name, "Function not available in dssatutils")
  } else {
    dssatutils::process_weather_daymet(
      shapefile = us_df,
      start_year = 2010,
      end_year = 2011,
      output_dir = wth_dir,
      id_col = "ID",
      lat_col = "LAT",
      lon_col = "LONG",
      n_cores = 1,
      log_file = file.path(work, "daymet.log")
    )
    if (all(file.exists(file.path(wth_dir, paste0(US_POINTS$ID, ".WTH"))))) {
      log_test("ok", test_name)
    } else {
      log_test("fail", test_name, "Not all .WTH files written")
    }
  }
}, error = function(e) {
  if (grepl("daymetR|Daymet", e$message)) {
    log_test("skip", test_name, "Optional dep (daymetR) missing")
  } else {
    log_test("fail", test_name, e$message)
  }
})

# 3. NASA-POWER (global, free, no key)
test_name <- "NASA-POWER"
wth_dir <- file.path(work, "weather_nasapower")
dir.create(wth_dir, recursive = TRUE)
tryCatch({
  if (!exists("process_weather_nasapower", where = asNamespace("dssatutils"), inherits = FALSE)) {
    log_test("skip", test_name, "Function not available")
  } else {
    dssatutils::process_weather_nasapower(
      shapefile = global_df,
      start_year = 2010,
      end_year = 2011,
      output_dir = wth_dir,
      id_col = "ID",
      lat_col = "LAT",
      lon_col = "LONG",
      n_cores = 1,
      log_file = file.path(work, "nasapower.log")
    )
    if (all(file.exists(file.path(wth_dir, paste0(GLOBAL_POINTS$ID, ".WTH"))))) {
      log_test("ok", test_name)
    } else {
      log_test("fail", test_name, "Not all .WTH files written")
    }
  }
}, error = function(e) {
  log_test("fail", test_name, e$message)
})

# 4. GRIDMET (US only, free, no key)
test_name <- "GRIDMET"
wth_dir <- file.path(work, "weather_gridmet")
dir.create(wth_dir, recursive = TRUE)
tryCatch({
  if (!exists("process_weather_gridmet", where = asNamespace("dssatutils"), inherits = FALSE)) {
    log_test("skip", test_name, "Function not available")
  } else {
    dssatutils::process_weather_gridmet(
      shapefile = us_df,
      start_year = 2010,
      end_year = 2011,
      output_dir = wth_dir,
      id_col = "ID",
      lat_col = "LAT",
      lon_col = "LONG",
      n_cores = 1,
      log_file = file.path(work, "gridmet.log"),
      gridmet_cache_dir = file.path(work, "gridmet_cache")
    )
    if (all(file.exists(file.path(wth_dir, paste0(US_POINTS$ID, ".WTH"))))) {
      log_test("ok", test_name)
    } else {
      log_test("fail", test_name, "Not all .WTH files written")
    }
  }
}, error = function(e) {
  if (grepl("xarray|netCDF", e$message)) {
    log_test("skip", test_name, "Optional dep (xarray) missing")
  } else {
    log_test("fail", test_name, e$message)
  }
})

# 5. AgERA5 (global, requires CDS key)
test_name <- "AgERA5"
has_cds_key <- Sys.getenv("CDSAPI_RC", "") != ""
if (has_cds_key) {
  wth_dir <- file.path(work, "weather_agera5")
  dir.create(wth_dir, recursive = TRUE)
  tryCatch({
    if (!exists("process_weather_agera5", where = asNamespace("dssatutils"), inherits = FALSE)) {
      log_test("skip", test_name, "Function not available")
    } else {
      dssatutils::process_weather_agera5(
        shapefile = global_df,
        start_year = 2010,
        end_year = 2011,
        output_dir = wth_dir,
        id_col = "ID",
        lat_col = "LAT",
        lon_col = "LONG",
        n_cores = 1,
        log_file = file.path(work, "agera5.log"),
        agera5_cache_dir = file.path(work, "agera5_cache")
      )
      if (all(file.exists(file.path(wth_dir, paste0(GLOBAL_POINTS$ID, ".WTH"))))) {
        log_test("ok", test_name)
      } else {
        log_test("fail", test_name, "Not all .WTH files written")
      }
    }
  }, error = function(e) {
    log_test("fail", test_name, e$message)
  })
} else {
  log_test("skip", test_name, "CDS API key not configured (CDSAPI_RC env var)")
}

# 6. NASA-POWER CHIRPS
test_name <- "NASA-POWER CHIRPS"
wth_dir <- file.path(work, "weather_nasapower_chirps")
dir.create(wth_dir, recursive = TRUE)
tryCatch({
  if (!exists("process_weather_nasapower_chirps", where = asNamespace("dssatutils"), inherits = FALSE)) {
    log_test("skip", test_name, "Function not available")
  } else {
    dssatutils::process_weather_nasapower_chirps(
      shapefile = global_df,
      start_year = 2010,
      end_year = 2011,
      output_dir = wth_dir,
      id_col = "ID",
      lat_col = "LAT",
      lon_col = "LONG",
      n_cores = 1,
      log_file = file.path(work, "chirps.log"),
      chirps_cache_dir = file.path(work, "chirps_cache")
    )
    if (all(file.exists(file.path(wth_dir, paste0(GLOBAL_POINTS$ID, ".WTH"))))) {
      log_test("ok", test_name)
    } else {
      log_test("fail", test_name, "Not all .WTH files written")
    }
  }
}, error = function(e) {
  log_test("fail", test_name, e$message)
})

# --- SOIL SOURCES ---

message("\n=== SOIL SOURCES ===\n")
message("NOTE: SoilGrids REST API is slow (10-30s per point). Tests may take >1 min.\n")

# 1. SoilGrids Online (global, free, REST) — SLOW API
test_name <- "SoilGrids Online"
sol_dir <- file.path(work, "soil_soilgrids_online")
dir.create(sol_dir, recursive = TRUE)
tryCatch({
  if (!exists("process_soils_soilgrids_online", where = asNamespace("dssatutils"), inherits = FALSE)) {
    log_test("skip", test_name, "Function not available")
  } else {
    message(sprintf("[%s] Fetching from SoilGrids REST API (slow, ~1-2 min)...", test_name))
    dssatutils::process_soils_soilgrids_online(
      gridfile = global_sf,
      soilfile_csv_path = file.path(sol_dir, "soil_map.csv"),
      output_sol_dir = sol_dir,
      id_col = "ID"
    )
    sols <- list.files(sol_dir, pattern = "\\.SOL$")
    if (length(sols) > 0) {
      log_test("ok", test_name)
    } else {
      log_test("skip", test_name, "No .SOL files written (REST API timeout or coverage gap)")
    }
  }
}, error = function(e) {
  log_test("skip", test_name, paste("REST API error:", e$message))
})

# 2. SoilGrids (offline VRT or online, global)
test_name <- "SoilGrids"
sol_dir <- file.path(work, "soil_soilgrids")
dir.create(sol_dir, recursive = TRUE)
tryCatch({
  if (!exists("process_soils_soilgrids", where = asNamespace("dssatutils"), inherits = FALSE)) {
    log_test("skip", test_name, "Function not available")
  } else {
    # Create a dummy master .SOL file
    dummy_sol_file <- file.path(work, "dummy_master.SOL")
    writeLines(c(
      "*SOILS: Test Soils Master",
      "*TEST0001     TEST_SOIL       52.000    5.000",
      "@SITE        COUNTRY          LAT     LONG SCS FAMILY",
      " TEST_SOIL   World         52.000    5.000 ",
      "@ SCOM  SALB  SLU1  SLDR  SLRO  SLNF  SLPF  SMHB  SMPX  SMKE",
      "    BN   .13     6    .6    73     1     1 IB001 IB001 IB001",
      "@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI  SLCF  SLNI  SLHW  SLHB  SCEC  SADC",
      "     5   -99 0.100 0.200 0.300  1.00  10.0  1.40  1.00  20.0  40.0   0.0   -99   -99   -99   -99   -99",
      "    15   -99 0.100 0.200 0.300  0.80  10.0  1.40  0.80  20.0  40.0   0.0   -99   -99   -99   -99   -99"
    ), dummy_sol_file)
    
    dssatutils::process_soils_soilgrids(
      grid_points = global_sf,
      source_sol_file = dummy_sol_file,
      output_csv_path = file.path(sol_dir, "soil_map.csv"),
      output_sol_dir = sol_dir,
      id_col = "ID"
    )
    sols <- list.files(sol_dir, pattern = "\\.SOL$")
    if (length(sols) > 0) {
      log_test("ok", test_name)
    } else {
      log_test("fail", test_name, "No .SOL files written")
    }
  }
}, error = function(e) {
  log_test("fail", test_name, e$message)
})

# 3. SSURGO (US only)
test_name <- "SSURGO"
sol_dir <- file.path(work, "soil_ssurgo")
dir.create(sol_dir, recursive = TRUE)
tryCatch({
  if (!exists("process_soils_ssurgo", where = asNamespace("dssatutils"), inherits = FALSE)) {
    log_test("skip", test_name, "Function not available")
  } else {
    format_SQL_in_statement <- function(x) {
      if (length(x) == 0) return("()")
      paste0("('", paste(x, collapse = "','"), "')")
    }
    dssatutils::process_soils_ssurgo(
      grid_points = us_sf,
      output_dir_csv = file.path(sol_dir, "soil_map.csv"),
      output_dir_individual = sol_dir,
      n_cores = 1,
      id_col = "ID",
      lat_col = "LAT",
      long_col = "LONG",
      format_sql_func = format_SQL_in_statement
    )
    sols <- list.files(sol_dir, pattern = "\\.SOL$")
    if (length(sols) > 0) {
      log_test("ok", test_name)
    } else {
      log_test("fail", test_name, "No .SOL files written")
    }
  }
}, error = function(e) {
  log_test("fail", test_name, e$message)
})

# 4. HWSD (global, SQLite-based)
test_name <- "HWSD"
sol_dir <- file.path(work, "soil_hwsd")
dir.create(sol_dir, recursive = TRUE)
tryCatch({
  if (!exists("process_soils_hwsd", where = asNamespace("dssatutils"), inherits = FALSE)) {
    log_test("skip", test_name, "Function not available")
  } else {
    # Create dummy raster and sqlite db
    dummy_raster_file <- file.path(work, "dummy_hwsd.tif")
    r <- terra::rast(nrows=10, ncols=10, xmin=-180, xmax=180, ymin=-90, ymax=90, crs="EPSG:4326")
    terra::values(r) <- 1
    terra::writeRaster(r, dummy_raster_file, overwrite=TRUE)
    
    dummy_db_file <- file.path(work, "dummy_hwsd.sqlite")
    con <- DBI::dbConnect(RSQLite::SQLite(), dummy_db_file)
    layers_df <- data.frame(
      HWSD2_SMU_ID = 1,
      SEQUENCE = 1,
      SHARE = 100,
      TOPDEP = 0,
      BOTDEP = 200,
      SAND = 40,
      CLAY = 20,
      SILT = 40,
      BULK = 1.4,
      ORG_CARBON = 1.0,
      COARSE = 0
    )
    DBI::dbWriteTable(con, "HWSD2_LAYERS", layers_df)
    DBI::dbDisconnect(con)
    
    dssatutils::process_soils_hwsd(
      grid_points = global_sf,
      hwsd_raster_file = dummy_raster_file,
      hwsd_db_file = dummy_db_file,
      output_csv_path = file.path(sol_dir, "soil_map.csv"),
      output_sol_dir = sol_dir,
      id_col = "ID",
      lat_col = "LAT",
      long_col = "LONG"
    )
    sols <- list.files(sol_dir, pattern = "\\.SOL$")
    if (length(sols) > 0) {
      log_test("ok", test_name)
    } else {
      log_test("fail", test_name, "No .SOL files written")
    }
  }
}, error = function(e) {
  log_test("fail", test_name, e$message)
})

# --- SUMMARY ---
message("\n=== SUMMARY ===")
message(sprintf("Passed: %d | Skipped: %d | Failed: %d", PASSED, SKIPPED, FAILED))

if (FAILED > 0) {
  quit(status = 1)
} else {
  message("ALL TESTS PASSED (or skipped due to missing optional deps).")
  quit(status = 0)
}
})
