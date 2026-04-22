#=======================================================================
#   DSSAT PIPELINE SCRIPT — Spatial / Gridded Crop Modeling
#   Modular soil & weather download, folder building, local or HPC run
#=======================================================================
#
# HOW TO GET STARTED (beginners):
#   1. Open this script in RStudio and click SOURCE.
#   2. Set DSSAT_EXE in your R console first (one-time):
#        Sys.setenv(DSSAT_EXE = "/path/to/dscsm048")           # Linux/macOS
#        Sys.setenv(DSSAT_EXE = "C:/DSSAT48/DSCSM048.exe")    # Windows
#   3. Download the demo boundary shapefile (one-time, see README):
#        mkdir -p shapefile
#        curl -L -o shapefile/tl_2024_us_state.zip \
#          "https://www2.census.gov/geo/tiger/TIGER2024/STATE/tl_2024_us_state.zip"
#        unzip -o shapefile/tl_2024_us_state.zip -d shapefile
#   4. Place a DSSAT template file in dssat_templates/ and set TEMPLATE_FILE_NAME below.
#
# --- THREE WAYS TO DEFINE YOUR SPATIAL DOMAIN ---
#
#  MODE A — Regular grid (DEFAULT for demo): the pipeline creates a regular
#            grid of points at GRID_SPACING_METERS spacing that falls inside
#            the boundary you specify in BOUNDARY_SHAPEFILE_NAME / STATE_NAME_FILTER.
#            Set USE_EXISTING_POINT_SHAPEFILE <- FALSE  (default)
#
#  MODE B — Custom shapefile (your own field/farm/study area): supply any
#            point or polygon shapefile — polygons are auto-converted to
#            centroids. Set:
#              USE_EXISTING_POINT_SHAPEFILE <- TRUE
#              EXISTING_POINT_SHAPEFILE_PATH <- "path/to/your/points.shp"
#            Useful for: field boundaries, research station coordinates,
#            administrative units, or any pre-defined point set.
#
#  MODE C — Cropland-only points from CDL/NLCD (run on farmland only):
#            First run the two landcover helper scripts to build a cropland
#            point shapefile, then feed it to this pipeline via MODE B.
#            See README section "Optional: run only on cropland" for the
#            full step-by-step workflow.
#
# All three modes produce the same pipeline inputs (a standardised point
# shapefile with ID / LAT / LONG columns), so Steps 1-4 are identical
# regardless of how you defined your spatial domain.
#=======================================================================

#-----------------------------------------------------------------------
# SECTION 0: MASTER CONFIGURATION
#-----------------------------------------------------------------------

# --- 1. Path & Platform Detection ---
# FUNCTION: Robustly detect the directory where this script is located
detect_project_dir <- function() {
  # A. Try RStudio API (Interactive Mode)
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    path <- tryCatch(dirname(rstudioapi::getSourceEditorContext()$path), error = function(e) NA)
    if (!is.na(path) && path != "") return(path)
  }
  
  # B. Try Command Line Arguments (Rscript / HPC Mode)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  
  # C. Fallback to Current Working Directory
  return(getwd())
}

MAIN_PROJECT_DIR <- detect_project_dir()

# --- Project-root aware paths (portable across machines/HPC) ---
# This script is designed to live in the project root of the repository.
# Helper R scripts (soil/weather) should live in: <repo_root>/r_scripts/
PROJECT_ROOT <- MAIN_PROJECT_DIR
R_SCRIPTS_DIR <- file.path(PROJECT_ROOT, "r_scripts")
SHAPEFILE_DIR <- file.path(PROJECT_ROOT, "shapefile")
GRIDPOINTS_DIR <- file.path(PROJECT_ROOT, "gridpoints")
WEATHER_ROOT_DIR <- file.path(PROJECT_ROOT, "weather")
SOIL_ROOT_DIR <- file.path(PROJECT_ROOT, "soil")
RUNS_ROOT_DIR <- file.path(PROJECT_ROOT, "dssat_runs")
RESULTS_ROOT_DIR <- file.path(PROJECT_ROOT, "results")

message(sprintf("Running in Project Directory: %s", MAIN_PROJECT_DIR))

# Detect Operating System for DSSAT Executable
os_system <- Sys.info()["sysname"]

if (os_system == "Windows") {
  # Standard Windows install location
  DSSAT_BASE    <- "C:/DSSAT48"
  DSSAT_EXE_NAME <- "DSCSM048.exe"
} else {
  # macOS / Linux: DSSAT must be compiled from source. No official installer.
  # See README.md → "Installing DSSAT on macOS (Apple Silicon)" for the full tutorial.
  # Quick summary for macOS Apple Silicon (M1/M2/M3/M4/M5):
  #
  #   brew install gcc cmake git
  #   git clone https://github.com/DSSAT/dssat-csm-os.git
  #   git clone https://github.com/DSSAT/dssat-csm-data.git
  #   mkdir -p ~/Documents/GitHub/DSSAT48
  #   cp -r ~/Documents/GitHub/dssat-csm-data/. ~/Documents/GitHub/DSSAT48/
  #   cd ~/Documents/GitHub/dssat-csm-os && mkdir build && cd build
  #   cmake .. -DCMAKE_BUILD_TYPE=RELEASE \
  #            -DCMAKE_Fortran_COMPILER=$(which gfortran) -G "Unix Makefiles"
  #   make -j$(sysctl -n hw.logicalcpu)
  #   cp bin/dscsm048 ~/Documents/GitHub/DSSAT48/dscsm048
  #   chmod +x ~/Documents/GitHub/DSSAT48/dscsm048
  #   cp ~/Documents/GitHub/DSSAT48/{MODEL.ERR,OUTPUT.CDE,DATA.CDE} \
  #      ~/Documents/GitHub/DSSAT48/StandardData/
  #   sed 's|/usr/local|/Users/YOUR_USERNAME/Documents/GitHub/DSSAT48|g' \
  #       ~/Documents/GitHub/dssat-csm-os/Data/DSSATPRO.L48 \
  #       > ~/Documents/GitHub/DSSAT48/DSSATPRO.L48
  #   echo 'export PATH="/Users/YOUR_USERNAME/Documents/GitHub/DSSAT48:$PATH"' >> ~/.zprofile
  #   source ~/.zprofile
  #
  # Override path without editing this file:
  #   Sys.setenv(DSSAT_EXE = "/Users/YOUR_USERNAME/Documents/GitHub/DSSAT48/dscsm048")
  #
  # Linux system-wide: change DSSAT_BASE to "/opt/DSSAT48"
  DSSAT_BASE    <- "~/Documents/GitHub/DSSAT48"
  DSSAT_EXE_NAME <- "dscsm048"
}

# Environment-variable override (recommended — keeps paths out of version control):
#   Sys.setenv(DSSAT_EXE  = "/full/path/to/dscsm048")   # overrides both BASE + NAME
#   Sys.setenv(DSSAT_BASE = "/path/to/DSSAT48")          # overrides base folder only
DSSAT_BASE <- Sys.getenv("DSSAT_BASE", unset = DSSAT_BASE)

# --- 2. Project Settings ---
# PROJECT_NAME: short label used in folder/file naming. No spaces.
PROJECT_NAME      <- "dssat_spatial_demo"
# Grid spacing in meters. Larger = fewer points = faster demo.
# Suggested: 50000 (50 km) for a first test; 5000–10000 for production runs.
GRID_SPACING_METERS <- 50000
# Crop module extension — "MZ" = maize, "WH" = wheat, "SB" = soybean, etc.
CROP_EXTENSION    <- "MZ"

# --- 2b. Optional: Run folder naming (decouple inputs vs runs) ---
# Keep weather/soil/gridpoint folder names tied to PROJECT_NAME/RESOLUTION/SOURCE,
# but allow the DSSAT run folder name (under dssat_runs/) to be anything you want.
RUN_TAG <- ""            # e.g. "run1", "calibA"; set "" to keep default naming
RUN_NAME_STYLE <- "grid"     # "grid"     => <GRID_BASE_NAME>_<RUN_TAG>
                             #              e.g., dssat_spatial_demo_50km_run1
                             # "scenario" => <GRID_BASE_NAME>_<WEATHER>_<SOIL>_<RUN_TAG>
                             #              e.g., dssat_spatial_demo_50km_DAYMET_SSURGO_run1
RUN_NAME_OVERRIDE <- ""      # if non-empty, this exact name is used for the run folder


# --- 2a. Spatial domain: choose ONE of the three modes (A / B / C) ---
#
# MODE A (DEFAULT): generate a regular grid from a boundary polygon.
#   Set USE_EXISTING_POINT_SHAPEFILE <- FALSE and configure the boundary
#   settings in section 2b below.
#
# MODE B: supply your own point or polygon shapefile.
#   Set USE_EXISTING_POINT_SHAPEFILE <- TRUE and point
#   EXISTING_POINT_SHAPEFILE_PATH at your file.
#   The pipeline auto-converts polygons to centroids, re-projects to
#   WGS84, and normalises IDs — so any shapefile format works.
#   Examples:
#     "gridpoints/my_field_centroids.shp"     # farm/field polygons
#     "gridpoints/study_sites.shp"            # research station points
#     "gridpoints/my_state_cropland_5k.shp"   # CDL-derived cropland grid
#                                             # (built with the landcover helpers)
#
# MODE C (CDL/NLCD cropland only):
#   Run r_scripts/landcover_raster.R then r_scripts/landcover_raster_to_gridpoints.R
#   to build a cropland point shapefile, then use MODE B to feed it here.
#   See README → "Optional: run only on cropland" for the full walkthrough.
USE_EXISTING_POINT_SHAPEFILE  <- FALSE   # TRUE = MODE B/C; FALSE = MODE A (demo default)
EXISTING_POINT_SHAPEFILE_PATH <- file.path(MAIN_PROJECT_DIR, "gridpoints", "my_points.shp")  # MODE B/C only

# --- 2b. Boundary settings (MODE A only — ignored when USE_EXISTING_POINT_SHAPEFILE = TRUE) ---
# BOUNDARY_SHAPEFILE_NAME: path relative to shapefile/ folder.
#   Default uses the US Census TIGER/Line state boundaries (download instructions in README).
BOUNDARY_SHAPEFILE_NAME  <- "tl_2024_us_state.shp"
# Set ENABLE_BOUNDARY_FILTER = TRUE to restrict the grid to one or more states/regions.
# BOUNDARY_FILTER_COLUMN: the attribute column to filter on (e.g., "NAME" for state name,
#                         "STUSPS" for two-letter abbreviation).
ENABLE_BOUNDARY_FILTER   <- TRUE
BOUNDARY_FILTER_COLUMN   <- "NAME"
# STATE_NAME_FILTER: one or more values to keep. Must match BOUNDARY_FILTER_COLUMN exactly.
# Demo: Montana at 50 km spacing yields ~60 grid points — fast for a first test.
STATE_NAME_FILTER        <- c("Iowa")

# --- 2c. Auto-Names & Naming Convention ---
if (GRID_SPACING_METERS < 1000) { 
  RESOLUTION_TAG <- paste0(GRID_SPACING_METERS, "m") 
} else { 
  RESOLUTION_TAG <- paste0(GRID_SPACING_METERS / 1000, "km") 
}

# 1. Base Grid Identity (Location + Resolution)
GRID_BASE_NAME <- paste0(PROJECT_NAME, "_", RESOLUTION_TAG)
BOUNDARY_FILTER_VALUE <- STATE_NAME_FILTER

# 2. Weather Settings
# WEATHER_SOURCE: "DAYMET" (North America, best for US), "NASA_POWER" (global),
#                 "GRIDMET" (US, high-res; requires extra packages)
WEATHER_SOURCE     <- "GRIDMET"
# Keep the date range short for a first test (longer ranges = more downloads).
WEATHER_START_YEAR <- 1982 #note: nasa_power not suitable/available before 1984
WEATHER_END_YEAR   <- 1983

# 3. Soil Settings
# SOIL_SOURCE: "SSURGO"          — US only, queries USDA SDA web service per point
#              "SOILGRIDS_10K"   — global, reads a pre-downloaded master .SOL file
#              "SOILGRIDS_ONLINE"— global, queries SoilGrids REST API per point
SOIL_SOURCE        <- "SSURGO"
# EXTERNAL_SOIL_FILE: only needed when SOIL_SOURCE is "SOILGRIDS_10K".
# Pre-formatted DSSAT-ready .SOL files at 10 km resolution (by country):
#   Harvard Dataverse: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/1PEEY0
# Based on: Folberth et al. (2019) Environ. Model. Softw. 111:218-228
#   https://www.sciencedirect.com/science/article/pii/S1364815218313033
# Download the country file you need, place it under SoilGrids/, and adjust the path below.
EXTERNAL_SOIL_FILE <- file.path(MAIN_PROJECT_DIR, "SoilGrids", "US.SOL")

# 4. Construct Dynamic Folder Names
# > Soil & Weather folders: Named by [Location]_[Resolution]_[Source]
SOIL_BASENAME <- paste0(GRID_BASE_NAME, "_", SOIL_SOURCE)
WEATHER_DIR_NAME <- paste0(GRID_BASE_NAME, "_", WEATHER_SOURCE)

# > Scenario ID (used for tracing back to input folders)
SCENARIO_ID <- paste0(GRID_BASE_NAME, "_", WEATHER_SOURCE, "_", SOIL_SOURCE)

# > DSSAT Run Folder: configurable name under dssat_runs/
#   Default (if RUN_TAG and RUN_NAME_OVERRIDE are empty): <SCENARIO_ID>
#   If RUN_TAG set: name becomes either <GRID_BASE_NAME>_<RUN_TAG> or <SCENARIO_ID>_<RUN_TAG> depending on RUN_NAME_STYLE
DEFAULT_RUN_NAME <- SCENARIO_ID

if (nzchar(RUN_NAME_OVERRIDE)) {
  DSSAT_RUN_NAME <- RUN_NAME_OVERRIDE
} else if (nzchar(RUN_TAG)) {
  if (RUN_NAME_STYLE == "scenario") {
    DSSAT_RUN_NAME <- paste0(SCENARIO_ID, "_", RUN_TAG)
  } else {
    DSSAT_RUN_NAME <- paste0(GRID_BASE_NAME, "_", RUN_TAG)
  }
} else {
  DSSAT_RUN_NAME <- DEFAULT_RUN_NAME
}

# Make sure the folder name is filesystem-friendly
DSSAT_RUN_NAME <- gsub("[^A-Za-z0-9_\\-]", "_", DSSAT_RUN_NAME)

# --- Paths (Dynamic) ---
GRIDPOINTS_SUBDIR <- "gridpoints"
CENTRAL_SOIL_DIR <- SOIL_ROOT_DIR
CENTRAL_WEATHER_DIR <- WEATHER_ROOT_DIR
RESULTS_SUBDIR <- "results"
GRIDMET_CACHE_DIR <- file.path(MAIN_PROJECT_DIR, "gridmet_netcdf_cache")

# Input Paths
GRIDPOINTS_OUTPUT_DIR <- file.path(MAIN_PROJECT_DIR, GRIDPOINTS_SUBDIR)
POINT_SHAPEFILE_NAME <- paste0(GRID_BASE_NAME, ".shp") # Shapefile is tied to Location/Res only
POINT_SHAPEFILE_PATH <- file.path(GRIDPOINTS_OUTPUT_DIR, POINT_SHAPEFILE_NAME)

# Run & Output Paths
DSSAT_RUN_DIR <- file.path(MAIN_PROJECT_DIR, "dssat_runs", DSSAT_RUN_NAME)
FINAL_OUTPUT_DIR <- file.path(MAIN_PROJECT_DIR, RESULTS_SUBDIR)
FINAL_RESULTS_PATH <- file.path(FINAL_OUTPUT_DIR, paste0(DSSAT_RUN_NAME, "_results.csv"))
FINAL_PLOT_PATH <- file.path(FINAL_OUTPUT_DIR, paste0(DSSAT_RUN_NAME, "_yield_map.png"))

POINT_ID_COLUMN <- "ID"
LAT_COLUMN <- "LAT"
LONG_COLUMN <- "LONG"

# --- 4a. Weather Extension Settings ---
EXTEND_WEATHER_DATA <- FALSE     # Set to FALSE to disable filling partial years
WEATHER_REFERENCE_YEAR <- 2025  # The historic year used to clone data for filling gaps

# Ensure your HMX template has "ID_SOIL" or "SOIL_ID" in the SLNO column!
TEMPLATE_SOIL_ID_PLACEHOLDER <- "SOIL_ID" 

# --- 5. DSSAT Settings ---
DSSAT_EXE_PATH <- file.path(DSSAT_BASE, DSSAT_EXE_NAME)
DSSAT_EXE_PATH <- Sys.getenv("DSSAT_EXE", unset = DSSAT_EXE_PATH)
TEMPLATE_DIR <- file.path(MAIN_PROJECT_DIR, "dssat_templates")
TEMPLATE_FILE_NAME <- "UFGA8201.MZX"  # DEMO PLACEHOLDER — replace with your own experiment file
                                       # UFGA8201.MZX is a maize file bundled with DSSAT for testing only.
                                       # Any valid DSSAT experiment file works (.MZX, .WHX, .SBX, etc.)
                                       # Match the extension to your CROP_EXTENSION setting above.
TEMPLATE_FILE_PATH <- file.path(TEMPLATE_DIR, TEMPLATE_FILE_NAME)

# --- 6. Run Mode ---
RUN_MODE <- "experiment" #experiment or sequence
TREATMENT_START <- 1
TREATMENT_END <- 4
SEQUENCE_START <- 1
SEQUENCE_END <- 1

# --- 6b. HPC Settings ---
# Set to TRUE to zip the run folder and delete the original (for cloud upload)
ZIP_FOR_HPC <- FALSE 

# --- 7. Switches ---
RUN_STEP_1_SOILS <- TRUE # Set to FALSE to only use already downloaded soil files
RUN_STEP_2_WEATHER <- TRUE  # Set to FALSE to only use already downloaded weather files
RUN_DSSAT_EXECUTION <- TRUE # Set to FALSE for HPC preparation
CLEANUP_RUN_FOLDERS <- FALSE # Set to TRUE to delete all simulation folders after run
RESUME_DSSAT_RUNS <- FALSE # Set to TRUE to skip creation of new folders

# --- 8. Parallel ---
N_CORES_TO_USE <- max(1, parallel::detectCores() - 4)
SOIL_CORES <- N_CORES_TO_USE
WEATHER_CORES <- N_CORES_TO_USE
DSSAT_CORES <- N_CORES_TO_USE


#-----------------------------------------------------------------------
# SECTION 1: LOAD LIBRARIES & HELPERS
#-----------------------------------------------------------------------
message("Loading libraries...")
packages_needed <- c("sf", "lubridate", "foreach", "doParallel", "parallel", "DSSAT", "stringr", "dplyr", "tidyverse", "R.utils", "processx", "soilDB", "pbapply", "tools", "rstudioapi", "zoo", "nasapower")

if (WEATHER_SOURCE == "DAYMET") packages_needed <- c(packages_needed, "daymetr")
if (WEATHER_SOURCE == "NASA_POWER") packages_needed <- c(packages_needed, "nasapower")
if (WEATHER_SOURCE == "GRIDMET") packages_needed <- c(packages_needed, "terra", "ncdf4", "httr")

for (pkg in packages_needed) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg == "climateR") remotes::install_github("mikejohnson51/climateR") else install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

options(DSSAT.CSM = DSSAT_EXE_PATH)
sf_use_s2(FALSE)

if (WEATHER_SOURCE == "GRIDMET") {
  dir.create(GRIDMET_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
  options(climateR.cache_dir = GRIDMET_CACHE_DIR)
}

# --- Load Helper Scripts (Robust Search) ---
possible_script_dirs <- c(
  R_SCRIPTS_DIR,
  tryCatch(dirname(rstudioapi::getSourceEditorContext()$path), error = function(e) NA),
  getwd(),
  MAIN_PROJECT_DIR
)

SCRIPT_DIR <- NA
for (dir in possible_script_dirs) {
  if (!is.na(dir) && dir != "" && file.exists(file.path(dir, "soil_ssurgo.R"))) {
    SCRIPT_DIR <- dir
    break
  }
}

if (is.na(SCRIPT_DIR)) {
  stop("Helper scripts (soil_ssurgo.R, etc.) not found. Expected them under <repo_root>/r_scripts/. Please check your folder layout and paths.")
}
message(sprintf("Sourcing helper scripts from: %s", SCRIPT_DIR))

source(file.path(SCRIPT_DIR, "weather_daymet.R"))
source(file.path(SCRIPT_DIR, "weather_nasapower.R"))
source(file.path(SCRIPT_DIR, "weather_gridmet.R"))
source(file.path(SCRIPT_DIR, "soil_ssurgo.R"))
source(file.path(SCRIPT_DIR, "soil_soilgrids.R"))
source(file.path(SCRIPT_DIR, "soil_soilgrids_online.R"))

# --- Helper: format_SQL_in_statement ---
format_SQL_in_statement <- function(x) {
  if (is.null(x) || length(x) == 0) return("('')")
  cleaned <- unique(na.omit(x))
  if (length(cleaned) == 0) return("('')")
  paste0("(", paste(paste0("'", cleaned, "'"), collapse = ","), ")")
}

# --- Helper: Delete Numbered Folders ---
delete_numbered_folders <- function(ids) {
  dirs <- list.dirs(full.names = TRUE, recursive = FALSE)
  numbered_dirs <- dirs[basename(dirs) %in% ids]
  if(length(numbered_dirs) > 0) {
    message(sprintf("Deleting %d old simulation folders...", length(numbered_dirs)))
    sapply(numbered_dirs, unlink, recursive = TRUE)
  }
}

# --- VALIDATION BLOCK (Place after Section 1) ---
message("Running Pre-flight Checks...")

# 1. Check Boundary File
#if (!file.exists(file.path(MAIN_PROJECT_DIR, BOUNDARY_SHAPEFILE_NAME))) {
#  stop(sprintf("CRITICAL: Boundary shapefile not found: %s", BOUNDARY_SHAPEFILE_NAME))
#}
# 1. Check Boundary File (only needed when we are generating a grid from a boundary)
if (!USE_EXISTING_POINT_SHAPEFILE) {
  if (!file.exists(file.path(SHAPEFILE_DIR, BOUNDARY_SHAPEFILE_NAME))) {
    stop(sprintf("CRITICAL: Boundary shapefile not found: %s", BOUNDARY_SHAPEFILE_NAME))
  }
}

# 2. Check Template File
if (!file.exists(TEMPLATE_FILE_PATH)) {
  stop(sprintf("CRITICAL: DSSAT Template file not found: %s", TEMPLATE_FILE_PATH))
}

# 3. Check External Soil File (if not SSURGO)
if (SOIL_SOURCE != "SSURGO" && SOIL_SOURCE != "SOILGRIDS_ONLINE" && !file.exists(EXTERNAL_SOIL_FILE)) {
  stop(sprintf("CRITICAL: External soil file needed for %s but not found at: %s", SOIL_SOURCE, EXTERNAL_SOIL_FILE))
}

message("All checks passed. Starting pipeline...")

#-----------------------------------------------------------------------
# SECTION 2: HELPER FUNCTIONS (Grid Creation & Weather Extension)
#-----------------------------------------------------------------------
create_grid_points <- function(boundary_shape, spacing_meters, output_path) {
  boundary_projected <- st_transform(boundary_shape, 5070)
  bbox <- st_bbox(boundary_projected)
  x_coords <- seq(floor(bbox$xmin), ceiling(bbox$xmax), by = spacing_meters)
  y_coords <- seq(floor(bbox$ymin), ceiling(bbox$ymax), by = spacing_meters)
  full_grid_df <- expand.grid(X = x_coords, Y = y_coords)
  grid_points_projected <- st_as_sf(full_grid_df, coords = c("X", "Y"), crs = 5070)
  points_in_boundary_projected <- st_filter(grid_points_projected, boundary_projected)
  
  if (nrow(points_in_boundary_projected) == 0) stop("STEP 0 FAILED: No grid points created.")
  
  points_with_coords <- st_transform(points_in_boundary_projected, 4326) %>%
    mutate(!!LAT_COLUMN := round(st_coordinates(.)[,2], 6),
           !!LONG_COLUMN := round(st_coordinates(.)[,1], 6),
           !!POINT_ID_COLUMN := sprintf("%08d", row_number()))
  
  st_write(points_with_coords, output_path, append = FALSE, delete_layer = TRUE)
  return(points_with_coords)
}

# Optional helper: load an existing point shapefile and standardize it to the pipeline schema
load_existing_points <- function(input_path, output_path,
                                 id_col = POINT_ID_COLUMN,
                                 lat_col = LAT_COLUMN,
                                 lon_col = LONG_COLUMN) {
  if (!file.exists(input_path)) stop(sprintf("Existing point shapefile not found at: %s", input_path))
  pts <- st_read(input_path, quiet = TRUE)
  
  # If it's not points, try converting polygons/lines to points (centroids) as a fallback
  gtype <- unique(as.character(st_geometry_type(pts)))
  if (!all(gtype %in% c("POINT", "MULTIPOINT"))) {
    message(sprintf("Existing shapefile geometry is [%s]; converting to centroids.", paste(gtype, collapse = ", ")))
    pts <- st_centroid(pts)
  }
  
  # Standardize to WGS84 and ensure LAT/LONG columns exist
  pts_ll <- st_transform(pts, 4326)
  coords <- st_coordinates(pts_ll)
  pts_ll[[lat_col]] <- round(coords[, 2], 6)
  pts_ll[[lon_col]] <- round(coords[, 1], 6)
  
  # Ensure an ID column exists and is unique; coerce to 8-digit character IDs used elsewhere in the pipeline
  if (!(id_col %in% names(pts_ll))) {
    pts_ll[[id_col]] <- sprintf("%08d", seq_len(nrow(pts_ll)))
  } else {
    ids <- as.character(pts_ll[[id_col]])
    bad <- is.na(ids) | ids == "" | duplicated(ids)
    if (any(bad)) {
      message("ID column exists but has NA/blank/duplicates; regenerating sequential IDs.")
      pts_ll[[id_col]] <- sprintf("%08d", seq_len(nrow(pts_ll)))
    } else {
      # normalize but keep stable ordering
      pts_ll[[id_col]] <- sprintf("%08d", match(ids, unique(ids)))
    }
  }
  
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  st_write(pts_ll, output_path, append = FALSE, delete_layer = TRUE, quiet = TRUE)
  return(pts_ll)
}


# --- NEW: Repeat historic block to fill many future years ---
# --- DEPRECATED (legacy): extend_weather_repeat_single ---
# Superseded by extend_weather_repeat_single_ignore_partial, which handles
# partial/incomplete years more robustly. Retained for reference only.
# Do NOT call directly.
extend_weather_repeat_single <- function(f, ref_start_year, ref_end_year, target_end_year) {
  # f: full path to .WTH file
  # ref_start_year/ref_end_year: inclusive historic block to repeat (e.g., 1984, 2025)
  # target_end_year: final year to reach (e.g., 2083)
  lines <- readLines(f, warn = FALSE)
  data_start_idx <- grep("^\\s*[0-9]+", lines)[1]
  if (is.na(data_start_idx)) return(NULL)
  header_lines <- lines[1:(data_start_idx - 1)]
  data_lines_raw <- lines[data_start_idx:length(lines)]
  
  # sanitize and read as numeric table (preserve first col as numeric date)
  data_lines_clean <- gsub("NA", " -99 ", data_lines_raw, fixed = TRUE)
  data_lines_clean <- gsub("NaN", " -99 ", data_lines_clean, fixed = TRUE)
  d <- tryCatch({
    read.table(text = data_lines_clean, header = FALSE, fill = TRUE,
               colClasses = "numeric", na.strings = c("-99", "-99.0", "-99.00"))
  }, error = function(e) return(NULL))
  if (is.null(d) || nrow(d) == 0) return(NULL)
  
  # determine date format (YYDDD or YYYYDDD)
  sample_date <- as.integer(d[1,1])
  if (nchar(as.character(sample_date)) <= 5) {
    year_format <- "YYDDD"
    get_year <- function(x) {
      yy <- floor(x / 1000)
      ifelse(yy < 80, 2000 + yy, 1900 + yy)
    }
    fmt_date <- function(y, doy) sprintf("%05d", (y %% 100) * 1000 + doy)
  } else {
    year_format <- "YYYYDDD"
    get_year <- function(x) floor(x / 1000)
    fmt_date <- function(y, doy) sprintf("%07d", y * 1000 + doy)
  }
  
  # helper: is leap
  is_leap <- function(yr) (yr %% 4 == 0 & yr %% 100 != 0) | (yr %% 400 == 0)
  
  # split existing data into year-blocks by extracting year and doy
  raw_dates <- as.integer(d[,1])
  years_present <- get_year(raw_dates)
  doys_present <- raw_dates %% 1000
  d$YEAR <- years_present
  d$DOY  <- doys_present
  
  last_year_in_file <- max(d$YEAR, na.rm = TRUE)
  if (last_year_in_file >= target_end_year) return(TRUE) # already good
  
  # Extract reference block(s) from the *file* (prefer exact ref years if present)
  ref_years <- seq(ref_start_year, ref_end_year)
  ref_blocks <- lapply(ref_years, function(yr) {
    rb <- d[d$YEAR == yr, , drop = FALSE]
    # if not present for this yr, skip
    if (nrow(rb) == 0) return(NULL)
    rb[order(rb$DOY), , drop = FALSE]
  })
  # if none found inside file, fall back to first 365 rows as single block
  ref_blocks <- Filter(Negate(is.null), ref_blocks)
  if (length(ref_blocks) == 0) {
    # fallback: take first 365/366 rows grouped by the first year in file
    fallback_year <- min(d$YEAR, na.rm = TRUE)
    ref_blocks <- list(d[d$YEAR == fallback_year, ])
  }
  
  # ensure each block has DOY column and numeric fields only (except YEAR/DOY)
  # prepare a data.frame of only numeric columns except YEAR/DOY in case of extra cols
  data_cols <- setdiff(names(d), c("YEAR", "DOY"))
  # Build sequence of target years to add
  years_to_add <- seq(last_year_in_file + 1, target_end_year)
  
  added_blocks <- list()
  ref_idx <- 1
  for (tgt_year in years_to_add) {
    # select next ref block (cycle through ref_blocks)
    ref_block <- ref_blocks[[ref_idx]]
    ref_idx <- ref_idx + 1
    if (ref_idx > length(ref_blocks)) ref_idx <- 1
    
    ref_is_leap <- is_leap(unique(ref_block$YEAR))
    tgt_is_leap <- is_leap(tgt_year)
    
    # work on a copy of numeric columns (preserve ordering)
    temp_blk <- ref_block
    # adapt DOY rows if leap mismatch
    if (ref_is_leap && !tgt_is_leap) {
      # drop DOY == 60 (Feb29) from ref before mapping
      temp_blk <- temp_blk[temp_blk$DOY != 60, , drop = FALSE]
      # renumber DOY to 1..365 (they already will be)
    } else if (!ref_is_leap && tgt_is_leap) {
      # insert a duplicate of DOY 59 (Feb28) at position 60
      if (any(temp_blk$DOY == 59)) {
        row_59 <- temp_blk[temp_blk$DOY == 59, , drop = FALSE]
        row_59$DOY <- 60
        # increment DOY >= 60 by 1
        temp_blk[temp_blk$DOY >= 60, "DOY"] <- temp_blk[temp_blk$DOY >= 60, "DOY"] + 1
        # insert row_59
        temp_blk <- rbind(temp_blk, row_59)
        temp_blk <- temp_blk[order(temp_blk$DOY), , drop = FALSE]
      } else {
        # if DOY 59 missing, just duplicate first row of March (doy 60) back to 60
        # fallback: duplicate first row and set DOY to 60
        tmp_ins <- temp_blk[1, , drop = FALSE]
        tmp_ins$DOY <- 60
        temp_blk$DOY[temp_blk$DOY >= 60] <- temp_blk$DOY[temp_blk$DOY >= 60] + 1
        temp_blk <- rbind(temp_blk, tmp_ins)
        temp_blk <- temp_blk[order(temp_blk$DOY), , drop = FALSE]
      }
    }
    
    # Now set YEAR to target year and set the date column value appropriately
    temp_blk$YEAR <- tgt_year
    new_date_nums <- as.integer(sprintf("%d", 0)) # placeholder
    # Build date number depending on year_format
    new_date_cols <- mapply(function(y, doy) {
      as.integer(if (year_format == "YYDDD") {
        (y %% 100) * 1000 + doy
      } else {
        y * 1000 + doy
      })
    }, y = temp_blk$YEAR, doy = temp_blk$DOY)
    
    temp_blk[,1] <- new_date_cols  # replace first column with new dates
    # ensure NA replacement consistent
    temp_blk[is.na(temp_blk)] <- -99.0
    
    # strip YEAR/DOY helper columns for binding
    temp_blk_out <- temp_blk[, setdiff(names(temp_blk), c("YEAR", "DOY")), drop = FALSE]
    added_blocks[[length(added_blocks) + 1]] <- temp_blk_out
  }
  
  # combine original data (remove trailing rows if they partially overlap?) + added blocks
  final_df <- rbind(d[, setdiff(names(d), c("YEAR", "DOY")), drop = FALSE],
                    do.call(rbind, added_blocks))
  
  # format output lines like original: first column date formatted, other numeric columns with width
  # detect date width
  if (year_format == "YYDDD") date_fmt <- "%05d" else date_fmt <- "%07d"
  date_col_str <- sprintf(date_fmt, as.integer(final_df[,1]))
  val_cols_list <- lapply(final_df[, -1, drop = FALSE], function(col) {
    col[is.na(col)] <- -99.0
    # keep decimals like existing files (use one decimal like the pipeline)
    sprintf("%6.1f", as.numeric(col))
  })
  formatted_body <- do.call(paste, c(list(date_col_str), val_cols_list, list(sep = "")))
  
  writeLines(c(header_lines, formatted_body), f)
  return(TRUE)
}

extend_weather_repeat_single_ignore_partial <- function(f,
                                                        ref_start_year,
                                                        ref_end_year,
                                                        target_end_year,
                                                        verbose = TRUE) {
  # v4: same as v3 but enforces canonical column NAMES (not just counts)
  lines <- readLines(f, warn = FALSE)
  data_start_idx <- grep("^\\s*[0-9]+", lines)[1]
  if (is.na(data_start_idx)) return(NULL)
  header_lines <- lines[1:(data_start_idx - 1)]
  data_lines_raw <- lines[data_start_idx:length(lines)]
  
  # sanitize; ensure we replace NA/NaN with "-99.0" text so read.table sees numeric -99
  data_lines_clean <- gsub("NA", " -99.0 ", data_lines_raw, fixed = TRUE)
  data_lines_clean <- gsub("NaN", " -99.0 ", data_lines_clean, fixed = TRUE)
  d <- tryCatch({
    read.table(text = data_lines_clean, header = FALSE, fill = TRUE,
               colClasses = "numeric", na.strings = c("-99", "-99.0", "-99.00"))
  }, error = function(e) {
    if (verbose) message("Failed to read file: ", f, " : ", conditionMessage(e))
    return(NULL)
  })
  if (is.null(d) || nrow(d) == 0) return(NULL)
  
  # detect date format
  sample_date <- as.integer(d[1,1])
  if (nchar(as.character(sample_date)) <= 5) {
    year_format <- "YYDDD"
    get_year_from_code <- function(x) { yy <- floor(x / 1000); ifelse(yy < 80, 2000 + yy, 1900 + yy) }
    make_date_code <- function(y, doy) as.integer(sprintf("%05d", (y %% 100) * 1000 + doy))
    date_width_fmt <- "%05d"
  } else {
    year_format <- "YYYYDDD"
    get_year_from_code <- function(x) floor(x / 1000)
    make_date_code <- function(y, doy) as.integer(sprintf("%07d", y * 1000 + doy))
    date_width_fmt <- "%07d"
  }
  
  is_leap <- function(yr) (yr %% 4 == 0 & yr %% 100 != 0) | (yr %% 400 == 0)
  
  # attach helper YEAR/DOY
  raw_dates <- as.integer(d[,1])
  years_present <- get_year_from_code(raw_dates)
  doys_present  <- raw_dates %% 1000
  d$YEAR <- years_present
  d$DOY  <- doys_present
  
  # find complete years
  years_unique <- sort(unique(d$YEAR))
  complete_years <- c()
  for (yr in years_unique) {
    expected_days <- if (is_leap(yr)) 366 else 365
    actual_days <- sum(d$YEAR == yr, na.rm = TRUE)
    if (actual_days == expected_days) complete_years <- c(complete_years, yr)
  }
  
  # choose chosen_ref_year and truncate to last full year
  if (length(complete_years) > 0) {
    if (ref_end_year %in% complete_years) {
      chosen_ref_year <- ref_end_year
    } else {
      prior_candidates <- complete_years[complete_years <= ref_end_year]
      chosen_ref_year <- if (length(prior_candidates) > 0) max(prior_candidates) else max(complete_years)
    }
    last_full_year <- max(complete_years)
    d_trunc <- d[d$YEAR <= last_full_year, , drop = FALSE]
  } else {
    warning(sprintf("No complete years found in %s; falling back to first 365 rows.", f))
    chosen_ref_year <- ref_end_year
    n_take <- min(365, nrow(d))
    d_trunc <- d[1:n_take, , drop = FALSE]
    last_full_year <- get_year_from_code(as.integer(d_trunc[nrow(d_trunc), 1]))
  }
  
  # canonical column names come from the truncated base (preserves original names)
  canonical_colnames <- names(d_trunc)                   # includes date-code column and any others (and YEAR/DOY if present)
  # but downstream we remove YEAR/DOY helper columns for base_df, so canonical for writing should match that subset:
  canonical_body_colnames <- setdiff(canonical_colnames, c("YEAR", "DOY"))
  
  canonical_ncol <- length(canonical_body_colnames)      # number of columns written out per row
  
  # coerce function that honors canonical names
  coerce_block_to_canonical <- function(block_df, src_name = "<block>") {
    # block_df expected to contain at least the date-code column + numeric columns. It *may* include YEAR/DOY helpers.
    # First, drop helper YEAR/DOY if present
    if ("YEAR" %in% names(block_df)) block_df$YEAR <- NULL
    if ("DOY"  %in% names(block_df)) block_df$DOY  <- NULL
    
    # now ensure column count matches canonical_ncol and names match canonical_body_colnames
    # If block_df has different number of columns, add filler columns or trim right-most columns
    if (ncol(block_df) < canonical_ncol) {
      n_missing <- canonical_ncol - ncol(block_df)
      filler <- as.data.frame(matrix(-99.0, nrow = nrow(block_df), ncol = n_missing))
      # name the filler columns using the canonical names for the missing positions
      names_existing <- names(block_df)
      missing_names <- canonical_body_colnames[(length(names_existing) + 1):canonical_ncol]
      names(filler) <- missing_names
      out <- cbind(block_df, filler)
      names(out) <- canonical_body_colnames  # enforce canonical ordering/names
      if (verbose) message(sprintf("Coerced block '%s': added %d filler cols", src_name, n_missing))
      return(out)
    } else if (ncol(block_df) > canonical_ncol) {
      # trim extra columns — keep left-most columns and then enforce canonical names on those kept
      out <- block_df[, seq_len(canonical_ncol), drop = FALSE]
      names(out) <- canonical_body_colnames
      if (verbose) message(sprintf("Coerced block '%s': trimmed from %d -> %d cols", src_name, ncol(block_df), canonical_ncol))
      return(out)
    } else {
      # equal count — but names may differ; rename to canonical order
      out <- block_df
      names(out) <- canonical_body_colnames
      return(out)
    }
  }
  
  # Build ref_block for chosen_ref_year (use d, which still has YEAR/DOY)
  if (year_format == "YYDDD") {
    yy_short <- chosen_ref_year %% 100
    ref_start_code <- yy_short * 1000
    ref_end_code   <- (yy_short + 1) * 1000
  } else {
    ref_start_code <- chosen_ref_year * 1000
    ref_end_code   <- (chosen_ref_year + 1) * 1000
  }
  ref_block <- d[d[,1] > ref_start_code & d[,1] < ref_end_code, , drop = FALSE]
  
  # fallback to last_full_year block if necessary
  if (nrow(ref_block) == 0 && last_full_year %in% years_unique) {
    if (year_format == "YYDDD") {
      yy2 <- last_full_year %% 100
      ref_start_code <- yy2 * 1000
      ref_end_code   <- (yy2 + 1) * 1000
    } else {
      ref_start_code <- last_full_year * 1000
      ref_end_code   <- (last_full_year + 1) * 1000
    }
    ref_block <- d[d[,1] > ref_start_code & d[,1] < ref_end_code, , drop = FALSE]
  }
  if (nrow(ref_block) == 0) {
    ref_block <- d[1:min(365, nrow(d)), , drop = FALSE]
    tmp_first_year <- get_year_from_code(as.integer(ref_block[1,1]))
    chosen_ref_year <- tmp_first_year
  }
  
  # coerce ref_block to canonical body columns
  ref_block_body <- coerce_block_to_canonical(ref_block, src_name = "ref_block")
  
  # prepare canonical base_df from truncated data (drop YEAR/DOY)
  base_df <- d_trunc[, canonical_body_colnames, drop = FALSE]
  base_df <- coerce_block_to_canonical(base_df, src_name = "base_df")  # defensive rename
  
  # if nothing to add
  if (last_full_year >= target_end_year) {
    final_df <- base_df
  } else {
    years_to_add <- seq(last_full_year + 1, target_end_year)
    
    # prepare MM-DD mapping from ref_block (we need ref DOY and the ref year)
    ref_DOY <- as.integer(ref_block[,1]) %% 1000
    ref_year_in_block <- get_year_from_code(as.integer(ref_block[1,1]))
    ref_dates_vec <- as.Date(paste0(ref_year_in_block, "-01-01")) + (ref_DOY - 1)
    ref_mmdd <- format(ref_dates_vec, "%m-%d")
    ref_vals_allcols <- ref_block_body  # already canonical and contains first column as date-code (we'll overwrite)
    
    added_list <- list()
    for (tgt_year in years_to_add) {
      tgt_expected_days <- if (is_leap(tgt_year)) 366 else 365
      tgt_dates <- as.Date(paste0(tgt_year, "-01-01")) + (0:(tgt_expected_days - 1))
      tgt_mmdd <- format(tgt_dates, "%m-%d")
      
      idx_rows <- integer(length(tgt_mmdd))
      for (i in seq_along(tgt_mmdd)) {
        mm <- tgt_mmdd[i]
        candidates <- which(ref_mmdd == mm)
        if (length(candidates) > 0) {
          idx_rows[i] <- candidates[1]
        } else if (mm == "02-29") {
          c2 <- which(ref_mmdd == "02-28")
          if (length(c2) > 0) idx_rows[i] <- c2[1]
          else idx_rows[i] <- 1
        } else {
          found <- FALSE
          for (k in 1:5) {
            prev_mm <- format(tgt_dates[i] - k, "%m-%d")
            pidx <- which(ref_mmdd == prev_mm)
            if (length(pidx) > 0) { idx_rows[i] <- pidx[1]; found <- TRUE; break }
          }
          if (!found) idx_rows[i] <- 1
        }
      }
      
      # select rows and coerce to canonical column names (will assign canonical_body_colnames)
      temp_out <- ref_vals_allcols[idx_rows, , drop = FALSE]
      temp_out <- coerce_block_to_canonical(temp_out, src_name = paste0("temp_out_", tgt_year))
      
      # set correct date codes into first column (which we named canonical_body_colnames[1])
      new_date_codes <- vapply(1:tgt_expected_days, function(doy) make_date_code(tgt_year, doy), integer(1))
      temp_out[, 1] <- as.integer(new_date_codes)
      
      temp_out[is.na(temp_out)] <- -99.0
      added_list[[length(added_list) + 1]] <- temp_out
    }
    
    final_df <- rbind(base_df, do.call(rbind, added_list))
  }
  
  # final formatting: first column date code, others numeric with one decimal
  final_df[,1] <- as.integer(final_df[,1])
  date_col_str <- sprintf(date_width_fmt, as.integer(final_df[,1]))
  if (ncol(final_df) >= 2) {
    val_cols_list <- lapply(final_df[, -1, drop = FALSE], function(col) {
      col[is.na(col)] <- -99.0
      sprintf("%6.1f", as.numeric(col))
    })
    formatted_body <- do.call(paste, c(list(date_col_str), val_cols_list, list(sep = "")))
  } else {
    formatted_body <- date_col_str
  }
  
  writeLines(c(header_lines, formatted_body), f)
  if (verbose) message("Wrote extended file: ", f, " up to ", target_end_year)
  return(TRUE)
}


# --- DEPRECATED (legacy): extend_weather_smart_single ---
# Superseded by extend_weather_repeat_single_ignore_partial.
# Retained for reference only. Do NOT call directly.
extend_weather_smart_single <- function(f, reference_year) {
  
  # --- STAGE 1: FAST READ ---
  lines <- readLines(f)
  
  # 1. Identify data start
  data_start_idx <- grep("^\\s*[0-9]+", lines)[1]
  if(is.na(data_start_idx)) return(NULL) # Skip empty
  
  header_lines <- lines[1:(data_start_idx - 1)]
  data_lines_raw <- lines[data_start_idx:length(lines)]
  
  # 2. Fast text sanitization
  data_lines_clean <- gsub("NA", " -99 ", data_lines_raw, fixed = TRUE)
  data_lines_clean <- gsub("NaN", " -99 ", data_lines_clean, fixed = TRUE)
  
  # 3. Read data
  d <- tryCatch({
    read.table(text = data_lines_clean, header = FALSE, fill = TRUE, 
               colClasses = "numeric", na.strings = c("-99", "-99.0", "-99.00"))
  }, error = function(e) return(NULL))
  
  if(is.null(d) || nrow(d) == 0) return(NULL)
  
  # --- STAGE 2: VECTORIZED INTERPOLATION ---
  if (any(is.na(d))) {
    d[, -1] <- lapply(d[, -1], function(x) {
      if (all(is.na(x))) return(rep(-99.0, length(x))) 
      if (!any(is.na(x))) return(x) 
      zoo::na.approx(x, rule = 2) 
    })
  }
  
  # --- STAGE 3: LOGIC & EXTENSION ---
  last_val <- tail(d[,1], 1)
  if(is.na(last_val)) return(NULL)
  date_str <- as.character(as.integer(last_val))
  
  if (nchar(date_str) == 5) {
    current_year_full <- ifelse(as.numeric(substr(date_str, 1, 2)) < 80, 
                                2000 + as.numeric(substr(date_str, 1, 2)), 
                                1900 + as.numeric(substr(date_str, 1, 2)))
    last_ddd <- as.numeric(substr(date_str, 3, 5))
    year_format <- "YYDDD"
    date_fmt <- "%05d"
  } else {
    current_year_full <- as.numeric(substr(date_str, 1, 4))
    last_ddd <- as.numeric(substr(date_str, 5, 7))
    year_format <- "YYYYDDD"
    date_fmt <- "%07d"
  }
  
  # Extract Reference Data
  #ref_yy_short <- reference_year %% 100
  #if (year_format == "YYDDD") {
  #  ref_start_year <- ref_yy_short * 1000
  #  ref_end_year <- (ref_yy_short + 1) * 1000
  #} else {
  #  ref_start_year <- reference_year * 1000
  #  ref_end_year <- (reference_year + 1) * 1000
  #}
  
  #ref_data <- d[d[,1] > ref_start_year & d[,1] < ref_end_year, ]
  #if(nrow(ref_data) == 0) ref_data <- d[1:min(365, nrow(d)), ]
  
  # --- NEW: Pick reference data only from *complete* years (365/366 days) ---
  # Build year / doy vectors depending on date format
  if (year_format == "YYDDD") {
    yrs <- floor(d[,1] / 1000)
    yrs_resolved <- ifelse(yrs < 80, 2000 + yrs, 1900 + yrs)
  } else {
    yrs_resolved <- floor(d[,1] / 1000)
  }
  doy_vals <- d[,1] %% 1000
  
  # compute counts per year and expected days (handle leap-year)
  is_leap <- function(yr) (yr %% 4 == 0 & yr %% 100 != 0) | (yr %% 400 == 0)
  years_unique <- sort(unique(yrs_resolved))
  complete_years <- c()
  for (yr in years_unique) {
    expected_days <- if (is_leap(yr)) 366 else 365
    actual_days <- sum(yrs_resolved == yr, na.rm = TRUE)
    if (actual_days == expected_days) complete_years <- c(complete_years, yr)
  }
  
  # Choose ref_data from only complete years.
  # Prefer user-specified reference_year if it's complete; else choose the most
  # recent complete year before or equal to reference_year; else choose last complete year.
  chosen_ref_year <- NULL
  if (length(complete_years) > 0) {
    if (reference_year %in% complete_years) {
      chosen_ref_year <- reference_year
    } else {
      # pick the latest complete year <= reference_year, if any
      prior_candidates <- complete_years[complete_years <= reference_year]
      if (length(prior_candidates) > 0) chosen_ref_year <- max(prior_candidates)
      else chosen_ref_year <- max(complete_years)  # fallback to latest available complete year
    }
    # slice ref_data for chosen_ref_year
    if (year_format == "YYDDD") {
      yy_short <- chosen_ref_year %% 100
      ref_start <- yy_short * 1000
      ref_end   <- (yy_short + 1) * 1000
    } else {
      ref_start <- chosen_ref_year * 1000
      ref_end   <- (chosen_ref_year + 1) * 1000
    }
    ref_data <- d[d[,1] > ref_start & d[,1] < ref_end, , drop = FALSE]
  } else {
    # NO complete years found — fall back to previous behavior (first 365 rows)
    ref_data <- d[1:min(365, nrow(d)), , drop = FALSE]
  }
  
  # Prepare Extension Rows
  rows_to_add <- list()
  is_leap <- function(yr) (yr %% 4 == 0 & yr %% 100 != 0) | (yr %% 400 == 0)
  days_in_current_year <- if(is_leap(current_year_full)) 366 else 365
  
  if (last_ddd < days_in_current_year) {
    rows_needed <- days_in_current_year - last_ddd
    filler <- ref_data[rep(seq_len(nrow(ref_data)), length.out = rows_needed), ]
    start_fill_doy <- last_ddd + 1
    doy_seq <- start_fill_doy + (1:nrow(filler)) - 1
    new_dates <- if(year_format=="YYDDD") (current_year_full %% 100)*1000 + doy_seq else current_year_full*1000 + doy_seq
    filler[, 1] <- new_dates
    rows_to_add[[1]] <- filler
  }
  
  # Phase B: Add Buffer Year
  next_year_full <- current_year_full + 1
  days_in_next <- if(is_leap(next_year_full)) 366 else 365
  buffer <- ref_data[rep(seq_len(nrow(ref_data)), length.out = days_in_next), ]
  doy_seq <- 1:nrow(buffer)
  new_dates <- if(year_format=="YYDDD") (next_year_full %% 100)*1000 + doy_seq else next_year_full*1000 + doy_seq
  buffer[, 1] <- new_dates
  rows_to_add[[2]] <- buffer
  
  final_df <- rbind(d, do.call(rbind, rows_to_add))
  
  # --- STAGE 4: VECTORIZED FORMATTING ---
  date_col_str <- sprintf(date_fmt, as.integer(final_df[,1]))
  val_cols_list <- lapply(final_df[, -1], function(col) {
    col[is.na(col)] <- -99.0
    sprintf("%6.1f", col)
  })
  formatted_body <- do.call(paste, c(list(date_col_str), val_cols_list, list(sep = "")))
  
  writeLines(c(header_lines, formatted_body), f)
  return(TRUE)
}

#-----------------------------------------------------------------------
# STEP 0: CREATE GRIDFILE
#-----------------------------------------------------------------------
# STEP 0: CREATE GRIDFILE
#-----------------------------------------------------------------------
message("STEP 0: PREPARING GRIDFILE / POINTS")
setwd(MAIN_PROJECT_DIR) 
dir.create(GRIDPOINTS_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (USE_EXISTING_POINT_SHAPEFILE) {
  message(sprintf("Using existing point shapefile: %s", EXISTING_POINT_SHAPEFILE_PATH))
  gridfile <- load_existing_points(EXISTING_POINT_SHAPEFILE_PATH, POINT_SHAPEFILE_PATH,
                                   id_col = POINT_ID_COLUMN, lat_col = LAT_COLUMN, lon_col = LONG_COLUMN)
} else {
  message(sprintf("Generating grid points at %sm spacing from boundary: %s", GRID_SPACING_METERS, BOUNDARY_SHAPEFILE_NAME))
  BOUNDARY_SHAPEFILE_PATH <- file.path(SHAPEFILE_DIR, BOUNDARY_SHAPEFILE_NAME)
  
  if (!file.exists(BOUNDARY_SHAPEFILE_PATH)) stop(paste("Shapefile not found at:", BOUNDARY_SHAPEFILE_PATH))
  
  boundary_sf <- st_read(BOUNDARY_SHAPEFILE_PATH)
  
  if (ENABLE_BOUNDARY_FILTER) {
    filter_string <- paste(BOUNDARY_FILTER_VALUE, collapse = ", ")
    message(sprintf("Filtering boundary where %s is in: [%s]", BOUNDARY_FILTER_COLUMN, filter_string))
    boundary_sf <- boundary_sf[boundary_sf[[BOUNDARY_FILTER_COLUMN]] %in% BOUNDARY_FILTER_VALUE, ]
    if (nrow(boundary_sf) == 0) stop("Filter resulted in 0 features.")
  }
  gridfile <- create_grid_points(boundary_sf, GRID_SPACING_METERS, POINT_SHAPEFILE_PATH)
}

#-----------------------------------------------------------------------
# STEP 1: SOIL DATA
#-----------------------------------------------------------------------
message("STEP 1: PROCESSING SOIL DATA")
if (RUN_STEP_1_SOILS) {
  dir.create(CENTRAL_SOIL_DIR, showWarnings = FALSE, recursive = TRUE)
  # Uses new naming: [Location]_[Res]_[Source]
  soilfile_path_prefix <- file.path(CENTRAL_SOIL_DIR, SOIL_BASENAME)
  
  soilfile_DSSAT <- paste0(soilfile_path_prefix, ".SOL")
  soilfile_CSV <- paste0(soilfile_path_prefix, ".CSV")
  individual_sol_output_folder <- paste0(soilfile_path_prefix, "_individual_SOL")
  dir.create(individual_sol_output_folder, recursive = TRUE, showWarnings = FALSE)
  
  if (SOIL_SOURCE == "SSURGO") {
    process_soils_ssurgo(gridfile, soilfile_CSV, individual_sol_output_folder, SOIL_CORES, 
                         POINT_ID_COLUMN, LAT_COLUMN, LONG_COLUMN, format_SQL_in_statement)
    
    # Optional Combine for Master File
    combine_sol_files_local <- function(input_folder, output_file_path) {
      sol_files <- list.files(path = input_folder, pattern = "\\.SOL$", full.names = TRUE)
      out_con <- file(output_file_path, open = "wt")
      cat("*SOILS: Combined\n", file = out_con)
      for (f in sol_files) {
        lines <- readLines(f, warn = FALSE)
        # Profile header lines start with "*" regardless of source (SSURGO, SoilGrids, custom)
        start <- grep("^\\*", lines)[1]
        if(!is.na(start)) writeLines(lines[start:length(lines)], out_con)
      }
      close(out_con)
    }
    combine_sol_files_local(individual_sol_output_folder, soilfile_DSSAT)
    
  } else if (SOIL_SOURCE == "SOILGRIDS_10K") {
    if (!file.exists(EXTERNAL_SOIL_FILE)) stop("External soil file not found.")
    process_soils_soilgrids(gridfile, EXTERNAL_SOIL_FILE, soilfile_CSV, 
                            output_sol_dir = individual_sol_output_folder,
                            id_col = POINT_ID_COLUMN)
  } else if (SOIL_SOURCE == "SOILGRIDS_ONLINE") {
    process_soils_soilgrids_online(gridfile, soilfile_CSV, 
                                   output_sol_dir = individual_sol_output_folder,
                                   id_col = POINT_ID_COLUMN)
  } else { stop("Unknown SOIL_SOURCE") }
}

#-----------------------------------------------------------------------
# STEP 2: WEATHER DATA
#-----------------------------------------------------------------------
message("STEP 2: DOWNLOADING WEATHER DATA")
if (RUN_STEP_2_WEATHER) {
  dir.create(CENTRAL_WEATHER_DIR, showWarnings = FALSE, recursive = TRUE)
  
  # Uses new naming: [Location]_[Res]_[Source]
  output_dir <- file.path(CENTRAL_WEATHER_DIR, WEATHER_DIR_NAME)
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # --- SMART RESUME BLOCK ---
  existing_files <- tools::file_path_sans_ext(list.files(output_dir, pattern = "\\.WTH$"))
  missing_mask <- ! (as.character(gridfile[[POINT_ID_COLUMN]]) %in% existing_files)
  points_to_process <- gridfile[missing_mask, ]
  
  if (nrow(points_to_process) == 0) {
    message("All weather files already exist. Skipping processing.")
  } else {
    message(sprintf("Resuming: Processing %d remaining points...", nrow(points_to_process)))
    
    log_file <- file.path(output_dir, "download_errors.log")
    
    common_args <- list(shapefile = points_to_process, 
                        start_year = WEATHER_START_YEAR, 
                        end_year = WEATHER_END_YEAR, 
                        output_dir = output_dir, 
                        id_col = POINT_ID_COLUMN, 
                        lat_col = LAT_COLUMN, 
                        lon_col = LONG_COLUMN, 
                        n_cores = WEATHER_CORES, 
                        log_file = log_file)
    
    if (WEATHER_SOURCE == "DAYMET") do.call(process_weather_daymet, common_args)
    else if (WEATHER_SOURCE == "NASA_POWER") do.call(process_weather_nasapower, common_args)
    else if (WEATHER_SOURCE == "GRIDMET") do.call(process_weather_gridmet, c(common_args, list(gridmet_cache_dir = GRIDMET_CACHE_DIR)))
  }
  
  # === WEATHER EXTENSION LOGIC (PARALLEL) ===
  if (EXTEND_WEATHER_DATA) {
    all_wth_files <- list.files(output_dir, pattern = "\\.WTH$", full.names = TRUE)
    message("Checking which files need extension...")
    
    files_to_extend <- Filter(function(f) {
      tryCatch({
        last_line <- tail(readLines(f, n = -1), 1)
        if (length(last_line) == 0) return(TRUE) 
        last_date_str <- regmatches(last_line, regexpr("^\\s*\\d+", last_line))
        if (length(last_date_str) == 0) return(TRUE)
        last_date_num <- as.numeric(last_date_str)
        if (last_date_num > 99999) last_year <- floor(last_date_num / 1000) 
        else {
          yr_short <- floor(last_date_num / 1000)
          last_year <- ifelse(yr_short < 80, 2000 + yr_short, 1900 + yr_short)
        }
        return(last_year <= WEATHER_END_YEAR)
      }, error = function(e) TRUE) 
    }, all_wth_files)
    
    if(length(files_to_extend) > 0) {
      message(sprintf("Running Fast Extension on %d files...", length(files_to_extend)))
      cl <- makeCluster(WEATHER_CORES)
      #clusterExport(cl, c("extend_weather_smart_single", "WEATHER_REFERENCE_YEAR"), envir = environment())
      #clusterEvalQ(cl, library(zoo)) 
      #parLapply(cl, files_to_extend, function(f) {
      #  extend_weather_smart_single(f, WEATHER_REFERENCE_YEAR)
      #})
      # Export the new function and parameters into the cluster
      clusterExport(cl, c("extend_weather_repeat_single_ignore_partial", 
                          "WEATHER_REFERENCE_YEAR",
                          "WEATHER_END_YEAR", 
                          "WEATHER_START_YEAR"), 
                    envir = environment())
      clusterEvalQ(cl, { library(zoo) }) 
      
      # Choose reference start/end years (assumes WEATHER_START_YEAR..WEATHER_REFERENCE_YEAR inclusive is your historic block)
      #ref_start_year <- WEATHER_START_YEAR
      #ref_end_year   <- WEATHER_REFERENCE_YEAR   # e.g., 1984..2025 if WEATHER_REFERENCE_YEAR == 2025
      
      #parLapply(cl, files_to_extend, function(f) {
      # call the new repeater: f, ref_start, ref_end, target_end_year
      #  extend_weather_repeat_single(f, ref_start_year, ref_end_year, WEATHER_END_YEAR)
      #})
      parLapply(cl, files_to_extend, function(f) {
        extend_weather_repeat_single_ignore_partial(
          f,
          ref_start_year = WEATHER_START_YEAR,
          ref_end_year   = WEATHER_REFERENCE_YEAR,
          target_end_year = WEATHER_END_YEAR
        )
      })
      
      stopCluster(cl)
    } else {
      message("All files appear to be already extended.")
    }
  }
}

#-----------------------------------------------------------------------
# STEP 3: RUN SIMULATIONS (RENAMING TO SOIL.SOL)
#-----------------------------------------------------------------------
message("STEP 3: RUNNING DSSAT SIMULATIONS")

# Check config
if (!exists("TEMPLATE_SOIL_ID_PLACEHOLDER")) TEMPLATE_SOIL_ID_PLACEHOLDER <- "SOIL_ID"

if (!dir.exists(DSSAT_RUN_DIR)) dir.create(DSSAT_RUN_DIR, recursive = TRUE)
setwd(DSSAT_RUN_DIR)

# Ensure points ID is character
points <- gridfile 
points[[POINT_ID_COLUMN]] <- as.character(points[[POINT_ID_COLUMN]])

# --- 3.1. Load Soil Mapping ---
soil_mapping_file <- file.path(CENTRAL_SOIL_DIR, paste0(SOIL_BASENAME, ".CSV"))

if(file.exists(soil_mapping_file)) {
  message("Loading soil mapping file...")
  soil_map <- read.csv(soil_mapping_file, colClasses = "character")
  soil_map[[POINT_ID_COLUMN]] <- as.character(soil_map[[POINT_ID_COLUMN]])
  common_cols <- setdiff(intersect(names(points), names(soil_map)), POINT_ID_COLUMN)
  if(length(common_cols) > 0) soil_map <- soil_map[, !names(soil_map) %in% common_cols]
  points <- left_join(points, soil_map, by = POINT_ID_COLUMN)
} else {
  message("WARNING: Soil mapping CSV not found. Attempting legacy logic.")
  if (SOIL_SOURCE == "SSURGO") points$SOIL_ID <- points[[POINT_ID_COLUMN]]
}

# --- 3.1b. CRITICAL FIX FOR SSURGO / MISSING COLUMNS ---
if (!("SOIL_ID" %in% names(points))) {
  if (SOIL_SOURCE == "SSURGO") {
    message("SSURGO Mode: Defaulting SOIL_ID to Grid Point ID.")
    points$SOIL_ID <- points[[POINT_ID_COLUMN]]
  } else {
    message("WARNING: SOIL_ID column missing. Initializing as NA.")
    points$SOIL_ID <- NA_character_
  }
}

individual_soil_folder <- file.path(CENTRAL_SOIL_DIR, paste0(SOIL_BASENAME, "_individual_SOL"))
weather_repo <- file.path(CENTRAL_WEATHER_DIR, WEATHER_DIR_NAME) 

# --- 3.2. Clean folders if not resuming ---
if (!RESUME_DSSAT_RUNS) {
  delete_numbered_folders(points[[POINT_ID_COLUMN]])
}

# --- 3.3. Create Folders Function ---
create_folders_and_files <- function(i) {
  ID <- points[[POINT_ID_COLUMN]][i]       
  assigned_soil_id <- points$SOIL_ID[i]    
  
  dir.create(ID, showWarnings = FALSE)
  
  if (is.na(assigned_soil_id)) {
    writeLines(c("*SOIL ERROR", "No Soil ID assigned"), file.path(ID, "SOIL.SOL"))
    return(NULL)
  }
  
  if (SOIL_SOURCE == "SSURGO") {
    source_filename <- paste0(ID, ".SOL")
    hmx_replacement_id <- ID
  } else {
    source_filename <- paste0(assigned_soil_id, ".SOL")
    hmx_replacement_id <- assigned_soil_id
  }
  
  # 1. HANDLE SOIL FILE 
  src_path <- file.path(individual_soil_folder, source_filename)
  dest_path <- file.path(ID, "SOIL.SOL")
  if (file.exists(src_path)) file.copy(src_path, dest_path, overwrite = TRUE)
  else writeLines(c("*SOIL ERROR", paste("Source missing:", source_filename)), dest_path)
  
  # 2. HANDLE EXPERIMENT FILE
  tryCatch({
    content <- readLines(TEMPLATE_FILE_PATH)
    if (!any(grepl(TEMPLATE_SOIL_ID_PLACEHOLDER, content, fixed=TRUE))) {
      if(any(grepl("ID_SOIL", content, fixed=TRUE))) {
        content <- gsub("ID_SOIL", hmx_replacement_id, content, fixed = TRUE)
      }
    } else {
      content <- gsub(TEMPLATE_SOIL_ID_PLACEHOLDER, hmx_replacement_id, content, fixed = TRUE)
    }
    content <- gsub("00000000", ID, content, fixed = TRUE)
    content <- gsub(tools::file_path_sans_ext(TEMPLATE_FILE_NAME), ID, content, fixed = TRUE)
    writeLines(content, file.path(ID, paste0(ID, ".", tools::file_ext(TEMPLATE_FILE_NAME))))
  }, error = function(e) return(NULL))
  
  # 3. HANDLE WEATHER FILE
  wth <- file.path(weather_repo, paste0(ID, ".WTH"))
  if (file.exists(wth)) file.copy(wth, file.path(ID, basename(wth)))
  
  # 4. EXTRAS
  files_to_copy <- list.files(TEMPLATE_DIR, pattern = "\\.(WDA|SDA|CUL|ECO|SPE)$", full.names = TRUE)
  if (length(files_to_copy) > 0) file.copy(files_to_copy, ID)
}

# --- 3.4. Execute Folder Creation ---
if (!RESUME_DSSAT_RUNS) {
  message("Creating simulation folders...")
  cl <- makeCluster(DSSAT_CORES)
  clusterExport(cl, varlist = c("points", "SOIL_SOURCE", "individual_soil_folder", "TEMPLATE_FILE_PATH", "create_folders_and_files", "TEMPLATE_FILE_NAME", "TEMPLATE_SOIL_ID_PLACEHOLDER", "POINT_ID_COLUMN", "weather_repo", "TEMPLATE_DIR"), envir = environment())
  parLapply(cl, 1:nrow(points), create_folders_and_files)
  stopCluster(cl)
}

# --- 3.5. Run Simulations ---
run_simulation <- function(ID) {
  options(DSSAT.CSM = DSSAT_EXE_PATH)
  tryCatch({ setwd(ID) }, error = function(e) return(NULL))
  on.exit(setwd(".."))
  
  template_ext <- tools::file_ext(TEMPLATE_FILE_NAME)
  experiment_file <- list.files(pattern = paste0("\\.", template_ext, "$"))[1]
  if (is.na(experiment_file)) return(NULL)
  
  results_template <- data.frame(
    point_id = character(), run_number = numeric(), treatment = numeric(), crop_code = character(),
    latitude = numeric(), longitude = numeric(), weather_station_id = character(), soil_profile_id = character(),
    dssat_file_id = character(), dssat_description = character(), planting_date = numeric(), emergence_date = numeric(),
    harvest_date = numeric(), year_planting = numeric(), year_harvest = numeric(), top_weight_kg_ha = numeric(),
    final_grain_kg_ha = numeric(), removed_residue_kg_ha = numeric(), soil_organic_carbon_start_kg_C_ha = numeric(),
    soil_organic_carbon_end_kg_C_ha = numeric(), soil_organic_carbon_delta_kg_C_ha = numeric(),
    final_irrigation_applications_count = numeric(), final_irrigation_amount_mm= numeric(),       
    inorganic_n_applied_count = numeric(), inorganic_n_applied_kg_ha = numeric(), nitrate_leaching_kg_ha = numeric(),
    cumulative_net_co2_emissions_kg_CO2_ha = numeric(), cumulative_n2o_emissions_kg_N_ha = numeric(),
    stringsAsFactors = FALSE
  )
  results <- results_template
  
  read_supp_file <- function(fname) {
    if(file.exists(fname)) {
      d <- try(suppressWarnings(read_csv(fname, show_col_types = FALSE)), silent = TRUE)
      if(inherits(d, "try-error") || is.null(d) || nrow(d) == 0) return(NULL)
      return(d)
    }
    return(NULL)
  }
  
  tryCatch({
    
    # MODE A: EXPERIMENT 
    if (RUN_MODE == "experiment") {
      batch_file_path <- file.path(getwd(), 'DSSBatch.V48')
      write_dssbatch(x = experiment_file, trtno = TREATMENT_START:TREATMENT_END, file_name = batch_file_path)
      run_dssat()
      
      summary <- suppressWarnings(read_csv('summary.csv', show_col_types = FALSE))
      
      if(nrow(summary) == 0) {
        treatments_vec <- TREATMENT_START:TREATMENT_END
        n_years <- (WEATHER_END_YEAR - WEATHER_START_YEAR)
        n_runs <- length(treatments_vec) * n_years
        summary <- tibble(
          RUNNO = 1:n_runs, TRNO = rep(treatments_vec, each = n_years),
          PYEAR = rep(WEATHER_START_YEAR:(WEATHER_END_YEAR - 1), times = length(treatments_vec)),
          CR = NA_character_, LAT = NA_real_, LONG = NA_real_, WSTA = NA_character_,
          SOIL_ID = NA_character_, EXNAME = NA_character_, TNAM = NA_character_,
          PDAT = NA_real_, EDAT = NA_real_, HDAT = NA_real_, HYEAR = NA_real_,
          CWAM = NA_real_, HWAM = NA_real_, BWAH = NA_real_, CO2EM = NA_real_, N2OEM = NA_real_
        )
      } else {
        summary$PYEAR = substr(summary$PDAT, 1, 4)
      }
      
      master_runs <- tibble(RUNNO = summary$RUNNO)
      soil_org <- read_supp_file('soilorg.csv')
      if(!is.null(soil_org)) {
        soil_org_sum <- soil_org %>% group_by(RUN) %>% summarise(SOMCT_start = head(SOMCT, 1), SOMCT_end = tail(SOMCT, 1)) %>% rename(RUNNO = RUN)
        soil_organic_summarized <- left_join(master_runs, soil_org_sum, by = "RUNNO")
      } else { soil_organic_summarized <- master_runs %>% mutate(SOMCT_start=NA, SOMCT_end=NA) }
      
      soil_ni <- read_supp_file('soilni.csv')
      if(!is.null(soil_ni)) {
        soil_ni_sum <- soil_ni %>% group_by(RUN) %>% summarise(NAPC = tail(NAPC, 1), NLCC = tail(NLCC, 1), `NI#M` = tail(`NI#M`,1)) %>% rename(RUNNO = RUN)
        soilnitrogen_summarized <- left_join(master_runs, soil_ni_sum, by = "RUNNO")
      } else { soilnitrogen_summarized <- master_runs %>% mutate(NAPC=NA, NLCC=NA, `NI#M`=NA) }
      
      soil_wat <- read_supp_file('soilwat.csv')
      if(!is.null(soil_wat)) {
        soil_wat_sum <- soil_wat %>% group_by(RUN) %>% summarise(`IR#C` = tail(`IR#C`, 1), IRRC = tail(IRRC, 1)) %>% rename(RUNNO = RUN)
        irrigation_summarized <- left_join(master_runs, soil_wat_sum, by = "RUNNO")
      } else { irrigation_summarized <- master_runs %>% mutate(`IR#C`=NA, IRRC=NA) }
      
      run_results <- data.frame(
        point_id = ID, run_number = summary$RUNNO, treatment = summary$TRNO, crop_code = summary$CR,
        latitude = summary$LAT, longitude = summary$LONG, weather_station_id = summary$WSTA,
        soil_profile_id = summary$SOIL_ID, dssat_file_id = summary$EXNAME, dssat_description = summary$TNAM,
        planting_date = summary$PDAT, emergence_date = summary$EDAT, harvest_date = summary$HDAT,
        year_planting = as.integer(summary$PYEAR), year_harvest = summary$HYEAR, top_weight_kg_ha = summary$CWAM,
        final_grain_kg_ha = summary$HWAM, removed_residue_kg_ha = summary$BWAH,
        soil_organic_carbon_start_kg_C_ha = soil_organic_summarized$SOMCT_start,
        soil_organic_carbon_end_kg_C_ha = soil_organic_summarized$SOMCT_end,
        soil_organic_carbon_delta_kg_C_ha = soil_organic_summarized$SOMCT_end - soil_organic_summarized$SOMCT_start,
        final_irrigation_applications_count = irrigation_summarized$`IR#C`, final_irrigation_amount_mm= irrigation_summarized$IRRC,
        inorganic_n_applied_count = soilnitrogen_summarized$`NI#M`, inorganic_n_applied_kg_ha = soilnitrogen_summarized$NAPC,
        nitrate_leaching_kg_ha = soilnitrogen_summarized$NLCC, cumulative_net_co2_emissions_kg_CO2_ha = summary$CO2EM,
        cumulative_n2o_emissions_kg_N_ha = summary$N2OEM
      )
      results <- rbind(results, run_results)
      
      # MODE B: SEQUENCE 
    } else if (RUN_MODE == "sequence") {
      for(trt in TREATMENT_START:TREATMENT_END) {
        seq_vec <- SEQUENCE_START:SEQUENCE_END
        n_seq <- length(seq_vec)
        batch_data <- tibble(FILEX = experiment_file, TRTNO = rep(trt, n_seq), RP = 1, SQ = seq_vec, OP = 1, CO = 0)
        write_dssbatch(batch_data)
        run_dssat(run_mode = "Q", suppress_output = TRUE)
        summary <- suppressWarnings(read_csv('summary.csv', show_col_types = FALSE))
        
        if (nrow(summary) > 0) {
          summary$PYEAR = substr(summary$PDAT, 1, 4)
          master_runs <- tibble(RUNNO = summary$RUNNO)
          
          soil_org <- read_supp_file('soilorg.csv')
          if(!is.null(soil_org)) {
            soil_org_sum <- soil_org %>% group_by(RUN) %>% summarise(SOMCT_start = head(SOMCT, 1), SOMCT_end = tail(SOMCT, 1)) %>% rename(RUNNO = RUN)
            soil_organic_summarized <- left_join(master_runs, soil_org_sum, by = "RUNNO")
          } else { soil_organic_summarized <- master_runs %>% mutate(SOMCT_start=NA, SOMCT_end=NA) }
          
          soil_ni <- read_supp_file('soilni.csv')
          if(!is.null(soil_ni)) {
            soil_ni_sum <- soil_ni %>% group_by(RUN) %>% summarise(NAPC = tail(NAPC, 1), NLCC = tail(NLCC, 1), `NI#M` = tail(`NI#M`,1)) %>% rename(RUNNO = RUN)
            soilnitrogen_summarized <- left_join(master_runs, soil_ni_sum, by = "RUNNO")
          } else { soilnitrogen_summarized <- master_runs %>% mutate(NAPC=NA, NLCC=NA, `NI#M`=NA) }
          
          soil_wat <- read_supp_file('soilwat.csv')
          if(!is.null(soil_wat)) {
            soil_wat_sum <- soil_wat %>% group_by(RUN) %>% summarise(`IR#C` = tail(`IR#C`, 1), IRRC = tail(IRRC, 1)) %>% rename(RUNNO = RUN)
            irrigation_summarized <- left_join(master_runs, soil_wat_sum, by = "RUNNO")
          } else { irrigation_summarized <- master_runs %>% mutate(`IR#C`=NA, IRRC=NA) }
          
          seq_results <- data.frame(
            point_id = ID, run_number = summary$RUNNO, treatment = summary$TRNO, crop_code = summary$CR,
            latitude = summary$LAT, longitude = summary$LONG, weather_station_id = summary$WSTA,
            soil_profile_id = summary$SOIL_ID, dssat_file_id = summary$EXNAME, dssat_description = summary$TNAM,
            planting_date = summary$PDAT, emergence_date = summary$EDAT, harvest_date = summary$HDAT,
            year_planting = as.integer(summary$PYEAR), year_harvest = summary$HYEAR, top_weight_kg_ha = summary$CWAM,
            final_grain_kg_ha = summary$HWAM, removed_residue_kg_ha = summary$BWAH,
            soil_organic_carbon_start_kg_C_ha = soil_organic_summarized$SOMCT_start,
            soil_organic_carbon_end_kg_C_ha = soil_organic_summarized$SOMCT_end,
            soil_organic_carbon_delta_kg_C_ha = soil_organic_summarized$SOMCT_end - soil_organic_summarized$SOMCT_start,
            final_irrigation_applications_count = irrigation_summarized$`IR#C`, final_irrigation_amount_mm= irrigation_summarized$IRRC,
            inorganic_n_applied_count = soilnitrogen_summarized$`NI#M`, inorganic_n_applied_kg_ha = soilnitrogen_summarized$NAPC,
            nitrate_leaching_kg_ha = soilnitrogen_summarized$NLCC, cumulative_net_co2_emissions_kg_CO2_ha = summary$CO2EM,
            cumulative_n2o_emissions_kg_N_ha = summary$N2OEM
          )
          results <- rbind(results, seq_results)
        }
      } 
    }
    
    if(nrow(results) > 0) write_csv(results, paste0("results_", ID, ".csv"), na = "")
    return(results) 
    
  }, error = function(e) { 
    message(sprintf("--- FATAL ERROR processing ID %s: %s ---", ID, e$message))
    return(NULL) 
  })
}

# --- 3.6. Parallel Execution ---
all_ids <- points[[POINT_ID_COLUMN]]

if (RUN_DSSAT_EXECUTION) {
  message("Starting DSSAT execution...")
  ids_to_run <- all_ids
  
  if (RESUME_DSSAT_RUNS) {
    check_file_exists <- function(ID) file.exists(file.path(ID, paste0("results_", ID, ".csv")))
    completed_mask <- unlist(pblapply(all_ids, check_file_exists))
    ids_to_run <- all_ids[!completed_mask]
  }
  
  if(length(ids_to_run) > 0) {
    cl <- makeCluster(DSSAT_CORES)
    clusterEvalQ(cl, { library(DSSAT); library(dplyr); library(readr); library(tibble); library(stringr) })
    clusterExport(cl, c("DSSAT_EXE_PATH", "RUN_MODE", "TREATMENT_START", "TREATMENT_END", 
                        "SEQUENCE_START", "SEQUENCE_END", "TEMPLATE_FILE_NAME", 
                        "run_simulation", "WEATHER_START_YEAR", "WEATHER_END_YEAR"), 
                  envir = .GlobalEnv)
    parLapply(cl, ids_to_run, run_simulation)
    stopCluster(cl)
  }
  
  # --- 3.7. Combine Results ---
  if (!dir.exists(FINAL_OUTPUT_DIR)) dir.create(FINAL_OUTPUT_DIR, recursive = TRUE)
  
  combine_results <- function(folder_ids) {
    all_results <- list()
    for (ID in folder_ids) {
      f <- file.path(ID, paste0("results_", ID, ".csv"))
      if (file.exists(f)) {
        all_results[[ID]] <- read_csv(f, show_col_types = FALSE,
                                      col_types = cols(point_id = col_character(), soil_profile_id = col_character(),
                                                       weather_station_id = col_character(), dssat_file_id = col_character(),
                                                       .default = col_guess()))
      }
    }
    bind_rows(all_results)
  }
  
  final_data <- combine_results(all_ids)
  if(!is.null(final_data) && nrow(final_data) > 0) {
    write_csv(final_data, FINAL_RESULTS_PATH)
    message("Results combined.")
  }
  
} else {
  # ============================================================================
  # HPC PREP MODE (Folders Created, No Execution)
  # ============================================================================
  message(paste(rep("=", 60), collapse = ""))
  message(" HPC PREP MODE COMPLETE")
  message(paste(rep("=", 60), collapse = ""))
  
  # --- WRITE RUN METADATA ---
  metadata_file <- file.path(DSSAT_RUN_DIR, "README_CONFIG.txt")
  metadata_content <- c(
    "============================================================",
    sprintf("DSSAT RUN CONFIGURATION: %s", Sys.time()),
    "============================================================",
    sprintf("Project Name:      %s", PROJECT_NAME),
    sprintf("Resolution:        %s (%d meters)", RESOLUTION_TAG, GRID_SPACING_METERS),
    sprintf("Grid/Shapefile:    %s", POINT_SHAPEFILE_NAME),
    sprintf("Scenario ID:       %s", SCENARIO_ID),
    sprintf("Run Folder Name:   %s", DSSAT_RUN_NAME),
    sprintf("Weather Folder:    %s", WEATHER_DIR_NAME),
    sprintf("Soil Folder:       %s", SOIL_BASENAME),
    "",
    "--- DATA SOURCES ---",
    sprintf("Weather Source:    %s", WEATHER_SOURCE),
    sprintf("Years:             %d - %d", WEATHER_START_YEAR, WEATHER_END_YEAR),
    sprintf("Weather Extended?: %s (Ref Year: %s)", EXTEND_WEATHER_DATA, WEATHER_REFERENCE_YEAR),
    sprintf("Soil Source:       %s", SOIL_SOURCE),
    "",
    "--- SIMULATION SETTINGS ---",
    sprintf("Mode:              %s", RUN_MODE),
    sprintf("Treatments:        %d - %d", TREATMENT_START, TREATMENT_END),
    sprintf("Sequences:         %d - %d", SEQUENCE_START, SEQUENCE_END),
    sprintf("Template File:     %s", TEMPLATE_FILE_NAME),
    "",
    "--- PATHS ---",
    sprintf("Soil Map CSV:      %s", soil_mapping_file)
  )
  writeLines(metadata_content, metadata_file)
  message("Metadata file created.")
  
  if (ZIP_FOR_HPC) {
    message(sprintf("Zipping run directory for HPC upload: %s", DSSAT_RUN_NAME))
    
    # 1. Setup paths
    parent_dir <- dirname(DSSAT_RUN_DIR)
    folder_name <- basename(DSSAT_RUN_DIR)
    zip_file_name <- paste0(folder_name, ".zip")
    zip_full_path <- file.path(parent_dir, zip_file_name)
    
    # 2. Change dir to parent to avoid absolute paths in the zip
    # We do NOT capture getwd() here because it is the folder we are about to delete.
    setwd(parent_dir)
    
    tryCatch({
      # 3. Zip the folder
      utils::zip(zip_file_name, files = folder_name, flags = "-r")
      
      # 4. Verify and Delete Original
      if (file.exists(zip_file_name)) {
        message("Zip created successfully. Deleting uncompressed folder to save space...")
        unlink(folder_name, recursive = TRUE)
        message(sprintf("READY: %s", zip_full_path))
      } else {
        message("ERROR: Zip file was not created. Original folder preserved.")
      }
      
    }, error = function(e) {
      message(sprintf("Error during zipping: %s", e$message))
    }, finally = {
      # 5. Return to Project Root (Safest location)
      setwd(MAIN_PROJECT_DIR)
    })
    
  } else {
    message(sprintf("Folders and Input Files created in: %s", DSSAT_RUN_DIR))
    message("DSSAT execution skipped. You can now transfer this directory to your HPC.")
  }
  
  message("Note: DSSBatch.V48 files are usually generated at runtime. Ensure your HPC script handles batch file creation if needed.")
}

#-----------------------------------------------------------------------
# STEP 4: VISUALIZE RESULTS
#-----------------------------------------------------------------------
if (RUN_DSSAT_EXECUTION) {
  message(sprintf("\n%s\nSTEP 4: VISUALIZING RESULTS\n%s", paste(rep("=", 60), collapse = ""), paste(rep("=", 60), collapse = "")))
  
  library(readr); library(dplyr); library(ggplot2); library(sf); library(tools)
  file_path <- FINAL_RESULTS_PATH
  
  if (file.exists(file_path)) {
    message("Loading results for final plot...")
    sim_data <- read_csv(file_path, show_col_types = FALSE) 
    
    avg_yield_data_by_treatment <- sim_data %>%
      filter(!is.na(longitude) & !is.na(latitude)) %>%
      group_by(point_id, latitude, longitude, treatment) %>% 
      summarise(avg_grain_yield = mean(final_grain_kg_ha, na.rm = TRUE), .groups = 'drop') %>%
      filter(!is.na(avg_grain_yield))
    
    unique_treatments <- unique(avg_yield_data_by_treatment$treatment)
    
    if(exists("boundary_sf")) {
      boundary_sf_4326 <- st_transform(boundary_sf, 4326)
      
      for (trt in unique_treatments) {
        message(sprintf("--- Generating map for Treatment %d ---", trt))
        data_for_this_plot <- avg_yield_data_by_treatment %>% filter(treatment == trt)
        avg_yield_data_sf <- st_as_sf(data_for_this_plot, coords = c("longitude", "latitude"), crs = 4326)
        
        yield_map <- ggplot() +
          geom_sf(data = boundary_sf_4326, fill = "grey90", color = "black", linewidth = 0.5) +
          geom_sf(data = avg_yield_data_sf, aes(color = avg_grain_yield), alpha = 0.8, size = 2.5) +
          scale_color_viridis_c(option = "plasma", name = "Avg. Yield (kg/ha)") +
          coord_sf(crs = 4326) + 
          labs(title = sprintf("Average Simulated Grain Yield (Treatment %d)", trt),
               subtitle = paste(sprintf("Weather Data: %s (%d-%d)", WEATHER_SOURCE, WEATHER_START_YEAR, WEATHER_END_YEAR), sep = "\n")) +
          theme_minimal() + theme(panel.background = element_rect(fill = "aliceblue", color = NA))
        
        base_plot_path <- tools::file_path_sans_ext(FINAL_PLOT_PATH)
        new_plot_path <- sprintf("%s_treatment%d.png", base_plot_path, trt)
        ggsave(new_plot_path, yield_map, width = 10, height = 8, dpi = 300)
        message(paste("Plot saved to:", new_plot_path))
      } 
    }
  } else {
    message(paste("Plotting skipped: '", file_path, "' not found. Did simulations run correctly?"))
  }
  message(sprintf("\n%s\nPIPELINE FINISHED!\n%s", paste(rep("*", 60), collapse = ""), paste(rep("*", 60), collapse = "")))
} else {
  message("Skipping Visualization (HPC Prep Mode)")
}
