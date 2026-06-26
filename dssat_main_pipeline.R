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
CODE_ROOT_DIR <- MAIN_PROJECT_DIR

message(sprintf("Running engine code from: %s", CODE_ROOT_DIR))

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

# --- Shared config overlay (config.yml) ---
# Loads config.yml (if present) so R and Python share one set of settings.
# cfg_get(key, default) returns the YAML value or the default below.
local({
  loader_dir <- tryCatch(dirname(this.path::this.path()), error = function(e) getwd())
  for (cand in c(file.path(CODE_ROOT_DIR, "config_loader.R"),
                 file.path(loader_dir, "config_loader.R"),
                 "config_loader.R")) {
    if (file.exists(cand)) { source(cand); break }
  }
})
if (!exists("cfg_get")) cfg_get <- function(key, default) default

resolve_config_path <- function(path, base = CODE_ROOT_DIR) {
  path <- trimws(as.character(if (is.null(path)) "" else path))
  if (!nzchar(path)) return("")
  if (grepl("^([A-Za-z]:|/|\\\\\\\\|~)", path)) return(normalizePath(path, mustWork = FALSE))
  normalizePath(file.path(base, path), mustWork = FALSE)
}

# --- Project-root aware paths (portable across machines/HPC) ---
# CODE_ROOT_DIR holds code/templates/static resources. OUTPUT_ROOT_DIR holds
# generated model artifacts (gridpoints, weather, soils, run folders, results).
# This lets consuming projects use the engine without filling the engine repo.
INPUT_ROOT_DIR  <- resolve_config_path(cfg_get("input_root_dir", CODE_ROOT_DIR), CODE_ROOT_DIR)
OUTPUT_ROOT_DIR <- resolve_config_path(cfg_get("output_root_dir", CODE_ROOT_DIR), CODE_ROOT_DIR)
PROJECT_ROOT <- OUTPUT_ROOT_DIR
R_SCRIPTS_DIR <- file.path(CODE_ROOT_DIR, "r_scripts")
SHAPEFILE_DIR <- resolve_config_path(cfg_get("shapefile_dir", file.path(INPUT_ROOT_DIR, "shapefile")), CODE_ROOT_DIR)
GRIDPOINTS_DIR <- resolve_config_path(cfg_get("gridpoints_dir", file.path(OUTPUT_ROOT_DIR, "gridpoints")), CODE_ROOT_DIR)
WEATHER_ROOT_DIR <- resolve_config_path(cfg_get("weather_root_dir", file.path(OUTPUT_ROOT_DIR, "weather")), CODE_ROOT_DIR)
SOIL_ROOT_DIR <- resolve_config_path(cfg_get("soil_root_dir", file.path(OUTPUT_ROOT_DIR, "soil")), CODE_ROOT_DIR)
RUNS_ROOT_DIR <- resolve_config_path(cfg_get("runs_root_dir", file.path(OUTPUT_ROOT_DIR, "dssat_runs")), CODE_ROOT_DIR)
RESULTS_ROOT_DIR <- resolve_config_path(cfg_get("results_root_dir", file.path(OUTPUT_ROOT_DIR, "results")), CODE_ROOT_DIR)

message(sprintf("Engine input/static root: %s", INPUT_ROOT_DIR))
message(sprintf("Engine generated-output root: %s", OUTPUT_ROOT_DIR))

# --- 2. Project Settings ---
# PROJECT_NAME: short label used in folder/file naming. No spaces.
PROJECT_NAME      <- cfg_get("project_name", "dssat_spatial_demo")
# Grid spacing in meters. Larger = fewer points = faster demo.
# Suggested: 50000 (50 km) for a first test; 5000–10000 for production runs.
GRID_SPACING_METERS <- cfg_get("grid_spacing_meters", 50000)
# Crop module extension — "MZ" = maize, "WH" = wheat, "SB" = soybean, etc.
CROP_EXTENSION    <- cfg_get("crop_extension", "MZ")

# --- 2b. Optional: Run folder naming (decouple inputs vs runs) ---
# Keep weather/soil/gridpoint folder names tied to PROJECT_NAME/RESOLUTION/SOURCE,
# but allow the DSSAT run folder name (under dssat_runs/) to be anything you want.
RUN_TAG <- cfg_get("run_tag", "")            # e.g. "run1", "calibA"; set "" to keep default naming
RUN_NAME_STYLE <- cfg_get("run_name_style", "grid")  # "grid"     => <GRID_BASE_NAME>_<RUN_TAG>
                             #              e.g., dssat_spatial_demo_50km_run1
                             # "scenario" => <GRID_BASE_NAME>_<WEATHER>_<SOIL>_<RUN_TAG>
                             #              e.g., dssat_spatial_demo_50km_DAYMET_SSURGO_run1
RUN_NAME_OVERRIDE <- cfg_get("run_name_override", "")  # if non-empty, this exact name is used for the run folder


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
USE_EXISTING_POINT_SHAPEFILE  <- cfg_get("use_existing_point_shapefile", FALSE)   # TRUE = MODE B/C; FALSE = MODE A (demo default)
EXISTING_POINT_SHAPEFILE_PATH <- cfg_get("existing_point_shapefile_path",
                                         file.path(GRIDPOINTS_DIR, "my_points.shp"))  # MODE B/C only

# --- 2b. Boundary settings (MODE A only — ignored when USE_EXISTING_POINT_SHAPEFILE = TRUE) ---
# BOUNDARY_SHAPEFILE_NAME: path relative to shapefile/ folder.
#   Default uses the US Census TIGER/Line state boundaries (download instructions in README).
BOUNDARY_SHAPEFILE_NAME  <- cfg_get("boundary_shapefile_name", "tl_2024_us_state.shp")
# Set ENABLE_BOUNDARY_FILTER = TRUE to restrict the grid to one or more states/regions.
# BOUNDARY_FILTER_COLUMN: the attribute column to filter on (e.g., "NAME" for state name,
#                         "STUSPS" for two-letter abbreviation).
ENABLE_BOUNDARY_FILTER   <- cfg_get("enable_boundary_filter", TRUE)
BOUNDARY_FILTER_COLUMN   <- cfg_get("boundary_filter_column", "NAME")
# STATE_NAME_FILTER: one or more values to keep. Must match BOUNDARY_FILTER_COLUMN exactly.
# Demo: Montana at 50 km spacing yields ~60 grid points — fast for a first test.
STATE_NAME_FILTER        <- as.character(cfg_get("state_name_filter", c("Iowa")))

# --- 2d. Optional cropland mask ---
# If enabled, grid cells are filtered to cropland-bearing cells before soil,
# weather, and DSSAT execution. The shapefile carries short field names because
# ESRI Shapefile truncates names longer than 10 characters:
#   crop_frac = cropland fraction [0-1]
#   crop_pct  = cropland percent [0-100]
#   crop_ha   = cropland hectares in the grid cell
#   cell_ha   = full grid-cell hectares
USE_CROPLAND_MASK <- isTRUE(as.logical(cfg_get("use_cropland_mask", FALSE)))
CROPLAND_RASTER_FILE <- as.character(cfg_get("cropland_raster_file", ""))
CROPLAND_CLASSES <- as.integer(unlist(cfg_get("cropland_classes", c(82))))
CROPLAND_MIN_FRACTION <- as.numeric(cfg_get("cropland_min_fraction", 0))
CROPLAND_STRICT <- isTRUE(as.logical(cfg_get("cropland_strict", FALSE)))
REUSE_CROPLAND_GRID <- isTRUE(as.logical(cfg_get("reuse_cropland_grid", TRUE)))

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
WEATHER_SOURCE     <- cfg_get("weather_source", "DAYMET")
# Keep the date range short for a first test (longer ranges = more downloads).
WEATHER_START_YEAR <- cfg_get("weather_start_year", 1982) #note: nasa_power not suitable/available before 1984
WEATHER_END_YEAR   <- cfg_get("weather_end_year", 1983)
# NASA_POWER_CHIRPS only: CHIRPS rainfall resolution "p05" (~5.5km) or "p25" (~28km).
CHIRPS_RESOLUTION  <- as.character(cfg_get("chirps_resolution", "p05"))

# 3. Soil Settings
# SOIL_SOURCE: "SSURGO"          — US only, queries USDA SDA web service per point
#              "SOILGRIDS_10K"   — global, reads a pre-downloaded master .SOL file
#              "AGMIP"           — global AgMIP/Han DSSAT-ready 5 arc-min .SOL file
#              "SOILGRIDS_ONLINE"— global, queries SoilGrids REST API per point
#              "POLARIS"         — US 30 m, streams POLARIS GeoTIFF tiles per point
SOIL_SOURCE        <- cfg_get("soil_source", "SOILGRIDS_10K")
STATSGO            <- isTRUE(as.logical(cfg_get("statsgo", FALSE)))
STANDARDIZE_LAYERS <- isTRUE(as.logical(cfg_get("standardize_layers", FALSE)))
# EXTERNAL_SOIL_FILE: needed when SOIL_SOURCE is "SOILGRIDS_10K" or "AGMIP".
# Pre-formatted DSSAT-ready .SOL files at 5 arc-min / ~10 km (by country):
#   Harvard Dataverse: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/1PEEY0
# Based on: Han et al. (2019) Environ. Model. Softw. 119:70-83
#   https://doi.org/10.1016/j.envsoft.2019.05.012
# Download the country file you need, place it under SoilGrids/, and adjust the path below.
EXTERNAL_SOIL_FILE <- cfg_get("external_soil_file",
                              file.path(INPUT_ROOT_DIR, "SoilGrids", "US.SOL"))
# SOILGRIDS_ONLINE only: "REST" (JSON API, rate-limited) or "VRT" (GDAL virtual
# rasters via terra; batch-friendly, better coverage). Drives USE_REST_API.
SOILGRIDS_MODE     <- toupper(cfg_get("soilgrids_mode", "REST"))
# POLARIS only (CONUS 30 m): statistic layer to build the profile from. "p50"
# (median) is the deterministic drop-in; p5/p95 are reserved for a future soil-
# uncertainty ensemble. Optional polaris_cache_dir caches the GeoTIFF tiles.
POLARIS_STAT       <- cfg_get("polaris_stat", "p50")
POLARIS_CACHE_DIR  <- { .v <- cfg_get("polaris_cache_dir", ""); if (nzchar(.v)) .v else NULL }
# HWSD only: paths to the FAO HWSD v2.0 raster (SMU IDs) + SQLite database,
# downloaded once from FAO (blank = script defaults under HWSD/).
HWSD_RASTER_FILE   <- cfg_get("hwsd_raster_file",
                              file.path(INPUT_ROOT_DIR, "HWSD", "HWSD2.bil"))
HWSD_DB_FILE       <- cfg_get("hwsd_db_file",
                              file.path(INPUT_ROOT_DIR, "HWSD", "HWSD2.sqlite"))
# E-OBS only: folder of pre-downloaded E-OBS NetCDFs (tx/tn/rr/qq...). Set
# eobs_use_cds: true to fetch an area subset via the Copernicus CDS instead.
EOBS_NC_DIR        <- cfg_get("eobs_nc_dir",
                              file.path(INPUT_ROOT_DIR, "eobs_netcdf"))
EOBS_USE_CDS       <- isTRUE(as.logical(cfg_get("eobs_use_cds", FALSE)))
# Xavier (Brazil) / CMFD (China) only: folders of pre-downloaded NetCDFs.
XAVIER_NC_DIR      <- cfg_get("xavier_nc_dir", file.path(INPUT_ROOT_DIR, "xavier_netcdf"))
CMFD_NC_DIR        <- cfg_get("cmfd_nc_dir", file.path(INPUT_ROOT_DIR, "cmfd_netcdf"))
# Large local/cache-backed weather products.
CHELSA_NC_DIR      <- cfg_get("chelsa_nc_dir", file.path(INPUT_ROOT_DIR, "chelsa_w5e5_netcdf"))
AGMERRA_NC_DIR     <- cfg_get("agmerra_nc_dir", file.path(INPUT_ROOT_DIR, "agmerra_netcdf"))
AGCFSR_NC_DIR      <- cfg_get("agcfsr_nc_dir", file.path(INPUT_ROOT_DIR, "agcfsr_netcdf"))
SILO_NC_DIR        <- cfg_get("silo_nc_dir", file.path(INPUT_ROOT_DIR, "silo_netcdf"))
PRISM_CACHE_DIR    <- cfg_get("prism_cache_dir", file.path(INPUT_ROOT_DIR, "prism_cache"))
MSWX_NC_DIR        <- cfg_get("mswx_nc_dir", file.path(INPUT_ROOT_DIR, "mswx_netcdf"))
MSWEP_NC_DIR       <- cfg_get("mswep_nc_dir", file.path(INPUT_ROOT_DIR, "mswep_netcdf"))
CRUJRA_NC_DIR      <- cfg_get("crujra_nc_dir", file.path(INPUT_ROOT_DIR, "crujra_netcdf"))
TERRACLIMATE_NC_DIR <- cfg_get("terraclimate_nc_dir", file.path(INPUT_ROOT_DIR, "terraclimate_netcdf"))
# LUCAS (Europe) only: the downloaded ESDAC LUCAS topsoil table (CSV/XLSX).
LUCAS_CSV          <- cfg_get("lucas_csv", file.path(INPUT_ROOT_DIR, "LUCAS", "lucas_topsoil.csv"))
# Local/cache-backed soil products.
HIHYDROSOIL_RASTER_DIR <- cfg_get("hihydrosoil_raster_dir", file.path(INPUT_ROOT_DIR, "HiHydroSoil"))
SLGA_RASTER_DIR        <- cfg_get("slga_raster_dir", file.path(INPUT_ROOT_DIR, "SLGA"))
WISE30SEC_RASTER_DIR   <- cfg_get("wise30sec_raster_dir", file.path(INPUT_ROOT_DIR, "WISE30sec"))
WOSIS_PROFILE_CSV      <- cfg_get("wosis_profile_csv", file.path(INPUT_ROOT_DIR, "WoSIS", "wosis_processed_profiles.csv"))

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
# Downloaded source-data caches (reanalysis grids, station pulls) live under the
# INPUT/engine root, NOT the per-study OUTPUT root. The raw grids are independent
# of which study consumes them, so keeping them with the engine lets every output
# project reuse one cache instead of re-downloading (e.g. GridMET 1984-2025) per
# study. These dirs are gitignored — they hold large generated downloads.
GRIDMET_CACHE_DIR <- file.path(INPUT_ROOT_DIR, "gridmet_netcdf_cache")
CHIRPS_CACHE_DIR <- file.path(INPUT_ROOT_DIR, "chirps_netcdf_cache")
AGERA5_CACHE_DIR <- file.path(INPUT_ROOT_DIR, "agera5_netcdf_cache")
AGERA5_MAX_CONCURRENT_REQUESTS <- as.integer(cfg_get("agera5_max_concurrent_requests", 4))
DWD_CACHE_DIR    <- file.path(INPUT_ROOT_DIR, "dwd_station_cache")
EOBS_CACHE_DIR   <- file.path(INPUT_ROOT_DIR, "eobs_cds_cache")

# Input Paths
GRIDPOINTS_OUTPUT_DIR <- GRIDPOINTS_DIR
ALL_LAND_POINT_SHAPEFILE_NAME <- paste0(GRID_BASE_NAME, ".shp")
CROPLAND_GRID_TAG <- ""
if (USE_CROPLAND_MASK) {
  class_tag <- paste(CROPLAND_CLASSES, collapse = "-")
  min_tag <- gsub("\\.", "p", format(CROPLAND_MIN_FRACTION, trim = TRUE, scientific = FALSE))
  CROPLAND_GRID_TAG <- gsub("[^A-Za-z0-9_\\-]", "_", paste0("_cropland_", class_tag, "_min", min_tag))
}
POINT_SHAPEFILE_NAME <- paste0(GRID_BASE_NAME, CROPLAND_GRID_TAG, ".shp")
ALL_LAND_POINT_SHAPEFILE_PATH <- file.path(GRIDPOINTS_OUTPUT_DIR, ALL_LAND_POINT_SHAPEFILE_NAME)
POINT_SHAPEFILE_PATH <- file.path(GRIDPOINTS_OUTPUT_DIR, POINT_SHAPEFILE_NAME)

# Run & Output Paths
DSSAT_RUN_DIR <- file.path(RUNS_ROOT_DIR, DSSAT_RUN_NAME)
FINAL_OUTPUT_DIR <- RESULTS_ROOT_DIR
FINAL_RESULTS_PATH <- file.path(FINAL_OUTPUT_DIR, paste0(DSSAT_RUN_NAME, "_results.csv"))
FINAL_PLOT_PATH <- file.path(FINAL_OUTPUT_DIR, paste0(DSSAT_RUN_NAME, "_yield_map.png"))

POINT_ID_COLUMN <- "ID"
LAT_COLUMN <- "LAT"
LONG_COLUMN <- "LONG"

# --- 4a. Weather Extension Settings ---
EXTEND_WEATHER_DATA <- isTRUE(as.logical(cfg_get("extend_weather_data", FALSE)))     # Set to FALSE to disable filling partial years
WEATHER_REFERENCE_YEAR <- as.integer(cfg_get("weather_reference_year", 2025))  # The historic year used to clone data for filling gaps
REPAIR_WEATHER_MISSING_VALUES <- isTRUE(as.logical(cfg_get("repair_weather_missing_values", FALSE)))
REPAIR_WEATHER_DATE_GAPS <- isTRUE(as.logical(cfg_get("repair_weather_date_gaps", FALSE)))
REPAIR_WEATHER_TEMPERATURE_INVERSIONS <- isTRUE(as.logical(cfg_get("repair_weather_temperature_inversions", FALSE)))
AUDIT_WEATHER_QUALITY <- isTRUE(as.logical(cfg_get("audit_weather_quality", FALSE)))
WEATHER_REPAIR_MAX_GAP_DAYS <- as.integer(cfg_get("weather_repair_max_gap_days", 3))
WEATHER_REPAIR_WINDOW_DAYS <- as.integer(cfg_get("weather_repair_window_days", 2))
WEATHER_QUALITY_FLATLINE_DAYS <- as.integer(cfg_get("weather_quality_flatline_days", 10))
WEATHER_REPAIR_VARIABLES <- as.character(unlist(cfg_get(
  "weather_repair_variables",
  c("SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND")
), use.names = FALSE))

# Soil/weather placeholders in the *FIELDS data row of the FileX template:
#   SID00000 -> per-point soil ID   (preferred; 8-char, ID_SOIL column only)
#   WID00000 -> per-point weather/WSTA ID (preferred; 8-char, WSTA column only)
# Legacy templates may instead use SOIL_ID / ID_SOIL for soil and 00000000 for
# WSTA; both engines still handle those as a fallback (see HANDLE EXPERIMENT FILE).
TEMPLATE_SOIL_ID_PLACEHOLDER <- "SOIL_ID"

# --- 5. DSSAT Settings ---
DSSAT_EXE_PATH <- file.path(DSSAT_BASE, DSSAT_EXE_NAME)
DSSAT_EXE_PATH <- Sys.getenv("DSSAT_EXE", unset = DSSAT_EXE_PATH)
TEMPLATE_DIR <- resolve_config_path(cfg_get("template_dir", file.path(INPUT_ROOT_DIR, "dssat_templates")), CODE_ROOT_DIR)
TEMPLATE_FILE_NAME <- cfg_get("template_file_name", "UFGA8201.MZX")  # DEMO PLACEHOLDER — replace with your own experiment file
                                       # UFGA8201.MZX is a maize file bundled with DSSAT for testing only.
                                       # Any valid DSSAT experiment file works (.MZX, .WHX, .SBX, etc.)
                                       # Match the extension to your CROP_EXTENSION setting above.
TEMPLATE_FILE_PATH <- file.path(TEMPLATE_DIR, TEMPLATE_FILE_NAME)

# --- 6. Run Mode ---
RUN_MODE <- cfg_get("run_mode", "experiment") #experiment or sequence
TREATMENT_START <- as.integer(cfg_get("treatment_start", 1))
TREATMENT_END <- as.integer(cfg_get("treatment_end", 4))
legacy_treatments_key <- cfg_get("treatments", NULL)
if (!is.null(legacy_treatments_key)) {
  stop("Config key 'treatments' is ambiguous in the gridded engine. ",
       "Use treatment_start/treatment_end for a contiguous range, or ",
       "treatment_list for explicit non-contiguous treatment IDs.",
       call. = FALSE)
}
TREATMENT_LIST <- cfg_get("treatment_list", NULL)
if (!is.null(TREATMENT_LIST)) {
  TREATMENT_LIST <- suppressWarnings(as.integer(unlist(TREATMENT_LIST, use.names = FALSE)))
  TREATMENT_LIST <- unique(TREATMENT_LIST[!is.na(TREATMENT_LIST)])
}
if (is.null(TREATMENT_LIST)) TREATMENT_LIST <- integer(0)
if (is.na(TREATMENT_START) || is.na(TREATMENT_END)) {
  stop("treatment_start and treatment_end must be valid integers.", call. = FALSE)
}
if (TREATMENT_END < TREATMENT_START) {
  stop(sprintf("treatment_end (%d) must be >= treatment_start (%d).",
               TREATMENT_END, TREATMENT_START), call. = FALSE)
}
if (length(TREATMENT_LIST) && any(TREATMENT_LIST < 1L)) {
  stop("treatment_list must contain positive integer treatment IDs.", call. = FALSE)
}
TREATMENT_RUN_LABEL <- if (length(TREATMENT_LIST)) {
  paste(TREATMENT_LIST, collapse = ",")
} else {
  sprintf("%d-%d", TREATMENT_START, TREATMENT_END)
}
SEQUENCE_START <- cfg_get("sequence_start", 1)
SEQUENCE_END <- cfg_get("sequence_end", 1)

# --- 6b. HPC Settings ---
# Set to TRUE to zip the run folder and delete the original (for cloud upload)
ZIP_FOR_HPC <- FALSE 

# --- 7. Switches ---
RUN_STEP_1_SOILS <- isTRUE(as.logical(cfg_get("run_step_1_soils", TRUE))) # Set to FALSE to only use already downloaded soil files
RUN_STEP_2_WEATHER <- isTRUE(as.logical(cfg_get("run_step_2_weather", TRUE)))  # Set to FALSE to only use already downloaded weather files
RUN_DSSAT_EXECUTION <- isTRUE(as.logical(cfg_get("run_dssat_execution", TRUE))) # Set to FALSE for HPC preparation
CLEANUP_RUN_FOLDERS <- isTRUE(as.logical(cfg_get("cleanup_run_folders", FALSE))) # Set to TRUE to delete all simulation folders after run
RESUME_DSSAT_RUNS <- isTRUE(as.logical(cfg_get("resume_dssat_runs", FALSE))) # Set to TRUE to skip creation of new folders

# Validation and retry controls
CHECK_WEATHER_DOWNLOADS <- isTRUE(as.logical(cfg_get("check_weather_downloads", FALSE)))
WEATHER_DOWNLOAD_RETRIES <- as.integer(cfg_get("weather_download_retries", 3))
CHECK_SOIL_DOWNLOADS <- isTRUE(as.logical(cfg_get("check_soil_downloads", FALSE)))
SOIL_DOWNLOAD_RETRIES <- as.integer(cfg_get("soil_download_retries", 3))


# Performance: when DSSATPRO.V48 sits next to the executable, DSSAT resolves
# genotype/SDA/CO2 files from the install dir, so they need not be copied into
# every run folder (big metadata saving at scale). Copy them only when shipping
# folders elsewhere (ZIP_FOR_HPC), when forced, or when no DSSATPRO is present.
BUNDLE_GENOTYPE_FILES <- isTRUE(as.logical(cfg_get("bundle_genotype_files", FALSE)))
COPY_SUPPORT_FILES <- BUNDLE_GENOTYPE_FILES || ZIP_FOR_HPC ||
  !file.exists(file.path(dirname(DSSAT_EXE_PATH), "DSSATPRO.V48"))

is_custom_or_missing <- function(fname, template_dir, dssat_dir) {
  ext <- toupper(tools::file_ext(fname))
  stock_folder <- if (ext %in% c("CUL", "ECO", "SPE")) "Genotype" else "StandardData"
  stock_path <- file.path(dssat_dir, stock_folder, fname)
  if (!file.exists(stock_path)) {
    return(TRUE)
  }
  local_path <- file.path(template_dir, fname)
  if (!file.exists(local_path)) {
    return(FALSE)
  }
  if (file.info(local_path)$size != file.info(stock_path)$size) {
    return(TRUE)
  }
  
  # Compare MD5
  get_hash <- function(p) {
    tryCatch({
      as.character(tools::md5sum(p))
    }, error = function(e) "")
  }
  
  if (get_hash(local_path) != get_hash(stock_path)) {
    return(TRUE)
  }
  return(FALSE)
}

# Precompute support files to copy
all_templates <- list.files(TEMPLATE_DIR, pattern = "\\.(WDA|SDA|CUL|ECO|SPE)$", full.names = FALSE)
if (COPY_SUPPORT_FILES) {
  SUPPORT_FILES <- all_templates
  message("Genotype files resolution: BUNDLE_GENOTYPE_FILES/ZIP_FOR_HPC or no DSSATPRO.V48 (copying all support files).")
} else {
  SUPPORT_FILES <- character()
  for (f in all_templates) {
    if (is_custom_or_missing(f, TEMPLATE_DIR, DSSAT_BASE)) {
      SUPPORT_FILES <- c(SUPPORT_FILES, f)
    }
  }
  message(sprintf("Genotype files resolved via DSSATPRO.V48. Copying %d custom/missing support file(s) per point.", length(SUPPORT_FILES)))
}

# --- 8. Parallel ---
# Core counts are read from config.yml (soil_cores / weather_cores / dssat_cores).
# "auto" (the default) = all logical CPUs minus 1, which is appropriate for most
# workloads. Set an explicit integer in config.yml to tune per stage:
#   - soil/weather steps are API/IO-bound  → can benefit from MORE than n_physical_cores
#   - DSSAT step is CPU+disk-bound          → stay at or below n_physical_cores
# makeCluster/parLapply (PSOCK, not fork) is used throughout the main pipeline,
# so all core counts work correctly on Windows, Linux, and macOS.
resolve_cores <- function(key, default_n) {
  v <- cfg_get(key, "auto")
  if (identical(v, "auto") || is.null(v) || !nzchar(as.character(v)))
    return(max(1L, default_n))
  max(1L, as.integer(v))
}
default_cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
SOIL_CORES    <- resolve_cores("soil_cores",    default_cores)
WEATHER_CORES <- resolve_cores("weather_cores", default_cores)
DSSAT_CORES   <- resolve_cores("dssat_cores",   default_cores)
message(sprintf("Parallelism: soil=%d  weather=%d  dssat=%d core(s)",
                SOIL_CORES, WEATHER_CORES, DSSAT_CORES))


#-----------------------------------------------------------------------
# SECTION 1: LOAD LIBRARIES & HELPERS
#-----------------------------------------------------------------------
message("Loading libraries...")

# Every CRAN package the pipeline attaches when it sources its helper modules.
# All helper modules are sourced at startup (some share helper functions), so
# the FULL set must be present for `source()` to succeed — that is why we check
# them all up front rather than only the selected weather/soil source's deps.
packages_needed <- c(
  # --- core pipeline ---
  "sf", "lubridate", "foreach", "doParallel", "parallel", "DSSAT", "stringr",
  "dplyr", "ggplot2", "tibble", "R.utils", "processx", "pbapply", "tools",
  "rstudioapi", "zoo",
  # --- weather / soil helper modules (attached at source time) ---
  "nasapower", "daymetr", "terra", "ncdf4", "httr", "jsonlite", "tidyr",
  "readr", "soilDB", "dssatutils", "dssatengine"
)
# Source-specific packages that are loaded LAZILY inside their module (so
# sourcing the file never needs them) — require them only when that source is
# actually selected, so e.g. a Daymet user is never asked to install RSQLite.
if (WEATHER_SOURCE == "AGERA5") packages_needed <- c(packages_needed, "ecmwfr")
if (SOIL_SOURCE == "HWSD")      packages_needed <- c(packages_needed, "DBI", "RSQLite")

# Packages installed from GitHub rather than CRAN (name -> repo). None are
# currently required beyond the shared DSSAT helper packages — the GRIDMET
# module downloads netCDF directly via httr + terra, so no climateR dependency.
# Add entries here only if a helper module genuinely `library()`s another
# GitHub-only package.
github_pkgs <- list(
  dssatutils = "alwinhopf/dssatutils",
  dssatengine = "alwinhopf/dssatengine"
)

# ---------------------------------------------------------------------------
# ensure_packages(): install any missing packages (prompting first in an
# interactive session, auto-installing in batch/Rscript), then attach them all.
# This makes "Source" in RStudio self-bootstrapping instead of erroring on the
# first missing dependency.
#
# Windows/macOS gotcha this handles: install.packages() defaults to
# pkgType = "both", so for packages whose CRAN *source* is newer than the
# *binary* (sf, terra, soilDB, ...) it tries to COMPILE from source — which
# needs Rtools + system GDAL/GEOS/PROJ and typically fails. We therefore force
# type = "binary" on Windows/macOS (with a source fallback), and route installs
# through renv::install() when an renv project is active so they land in the
# project library rather than fighting the renv shim.
# ---------------------------------------------------------------------------

# Make sure a real CRAN mirror is set (a bare "@CRAN@" makes installs prompt /
# fail in non-interactive Rscript sessions).
.repos <- getOption("repos")
if (is.null(.repos) || is.na(.repos["CRAN"]) || .repos["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

# Prefer CRAN binaries on Windows/macOS; Linux CRAN has no binaries (use source).
.is_binary_platform <- .Platform$OS.type == "windows" ||
                       Sys.info()[["sysname"]] == "Darwin"
.install_type <- if (.is_binary_platform) "binary" else "source"
# renv sets RENV_PROJECT when a project is active; route installs through it.
.renv_active <- nzchar(Sys.getenv("RENV_PROJECT")) &&
                requireNamespace("renv", quietly = TRUE)

# Install a single package (CRAN name or, if github_repo given, a GitHub repo),
# preferring binaries and falling back to source if the binary path fails.
.install_one <- function(pkg, github_repo = NULL) {
  target <- if (!is.null(github_repo)) github_repo else pkg
  attempt <- function(type) {
    if (.renv_active) {
      renv::install(target, type = type, prompt = FALSE)
    } else if (!is.null(github_repo)) {
      remotes::install_github(github_repo, upgrade = "never")
    } else {
      install.packages(pkg, type = type)
    }
  }
  ok <- tryCatch({ attempt(.install_type); TRUE },
                 error = function(e) {
                   message(sprintf("  binary/default install of '%s' failed: %s",
                                   pkg, conditionMessage(e))); FALSE })
  # Source fallback only meaningful where we tried a binary first.
  if (!ok && .is_binary_platform && is.null(github_repo)) {
    message(sprintf("  retrying '%s' from source ...", pkg))
    ok <- tryCatch({ attempt("source"); TRUE },
                   error = function(e) {
                     message(sprintf("  source install of '%s' also failed: %s",
                                     pkg, conditionMessage(e))); FALSE })
  }
  ok
}

ensure_packages <- function(pkgs, github = list()) {
  pkgs    <- unique(pkgs)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) > 0) {
    message("\nThe following required R package(s) are not installed:\n  ",
            paste(missing, collapse = ", "), "\n")
    message(sprintf("Install method: %s%s",
                    if (.renv_active) "renv::install" else "install.packages",
                    sprintf(" (type = '%s')", .install_type)))
    do_install <- TRUE
    if (interactive()) {
      ans <- readline(sprintf(
        "Install these %d package(s) now? [Y/n]: ", length(missing)))
      do_install <- !nzchar(ans) || tolower(substr(ans, 1, 1)) == "y"
    }
    if (!do_install)
      stop("Required packages are missing and installation was declined. ",
           "Install them manually, then re-run.", call. = FALSE)

    if (length(intersect(missing, names(github))) > 0 &&
        !requireNamespace("remotes", quietly = TRUE)) {
      .install_one("remotes")
    }
    for (pkg in missing) {
      message(sprintf("Installing '%s' ...", pkg))
      .install_one(pkg, github_repo = github[[pkg]])
    }
  }

  # Attach everything; collect (don't stop on) any that still won't load.
  failed <- pkgs[!vapply(pkgs, function(p)
    suppressWarnings(suppressMessages(
      require(p, character.only = TRUE, quietly = TRUE))), logical(1))]
  if (length(failed) > 0)
    stop("These package(s) could not be installed/loaded automatically:\n  ",
         paste(failed, collapse = ", "),
         sprintf("\nTry installing them manually as binaries, e.g.:\n  %s",
                 if (.renv_active)
                   sprintf("renv::install(c(%s), type = \"binary\")",
                           paste(sprintf('"%s"', failed), collapse = ", "))
                 else
                   sprintf("install.packages(c(%s), type = \"binary\")",
                           paste(sprintf('"%s"', failed), collapse = ", "))),
         "\nthen re-run.", call. = FALSE)
}

ensure_packages(packages_needed, github = github_pkgs)

options(DSSAT.CSM = DSSAT_EXE_PATH)
sf_use_s2(FALSE)

if (WEATHER_SOURCE == "GRIDMET") {
  dir.create(GRIDMET_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
}

# --- Load Helper Scripts (Robust Search) ---
possible_script_dirs <- c(
  R_SCRIPTS_DIR,
  tryCatch(dirname(rstudioapi::getSourceEditorContext()$path), error = function(e) NA),
  getwd(),
  MAIN_PROJECT_DIR
)

# Weather/soil helpers now come from the dssatutils package (loaded below); only
# the landcover helpers remain in r_scripts/, so detect the dir via one of those.
SCRIPT_DIR <- NA
for (dir in possible_script_dirs) {
  if (!is.na(dir) && dir != "" && file.exists(file.path(dir, "landcover_raster.R"))) {
    SCRIPT_DIR <- dir
    break
  }
}

if (is.na(SCRIPT_DIR)) {
  stop("Landcover helper scripts (landcover_raster.R, etc.) not found. Expected them under <repo_root>/r_scripts/. Please check your folder layout and paths.")
}
message(sprintf("Using landcover helper scripts from: %s", SCRIPT_DIR))

library(dssatutils)  # [dssatutils] shared weather/soil sources
library(dssatengine) # [dssatengine] shared gridded run engine
# Soil sources that write one Saxton-&-Rawls .SOL per grid point named by point
# ID (so SOIL_ID == point ID and the per-point combine logic applies).
PER_POINT_SOIL <- c("SSURGO", "GNATSGO", "ISDASOIL", "LUCAS", "SSURGO_ALDERMAN")
# Soil sources needing no pre-downloaded external file (queried online).
KEYLESS_ONLINE_SOIL <- c("SSURGO", "GNATSGO", "ISDASOIL", "SOILGRIDS_ONLINE", "POLARIS", "SSURGO_ALDERMAN")
# [dssatutils] source(file.path(SCRIPT_DIR, "weather_daymet.R"))
# [dssatutils] source(file.path(SCRIPT_DIR, "weather_nasapower.R"))
# [dssatutils] source(file.path(SCRIPT_DIR, "weather_gridmet.R"))
# [dssatutils] source(file.path(SCRIPT_DIR, "weather_openmeteo.R"))
# [dssatutils] source(file.path(SCRIPT_DIR, "weather_nasapower_chirps.R"))
# [dssatutils] source(file.path(SCRIPT_DIR, "weather_agera5.R"))
# [dssatutils] source(file.path(SCRIPT_DIR, "soil_ssurgo.R"))
# [dssatutils] source(file.path(SCRIPT_DIR, "soil_soilgrids.R"))
# [dssatutils] source(file.path(SCRIPT_DIR, "soil_soilgrids_online.R"))
# [dssatutils] source(file.path(SCRIPT_DIR, "soil_hwsd.R"))

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

# --- Helper: Validate Weather File ---
is_wth_valid <- function(id, dir, end_yr) {
  f <- file.path(dir, paste0(id, ".WTH"))
  if (!file.exists(f)) return(FALSE)
  if (file.info(f)$size == 0) return(FALSE)
  
  tryCatch({
    lines <- tail(readLines(f, warn = FALSE, encoding = "UTF-8"), 50)
    if (length(lines) == 0) return(FALSE)
    data_lines <- grep("^\\s*[0-9]{5,7}\\b", lines, value = TRUE)
    if (length(data_lines) == 0) return(FALSE)
    last_line <- tail(data_lines, 1)
    last_date_str <- regmatches(last_line, regexpr("^\\s*\\d+", last_line))
    if (length(last_date_str) == 0) return(FALSE)
    last_date_num <- as.numeric(last_date_str)
    if (last_date_num > 99999) {
      last_year <- floor(last_date_num / 1000)
    } else {
      yr_short <- floor(last_date_num / 1000)
      last_year <- ifelse(yr_short < 80, 2000 + yr_short, 1900 + yr_short)
    }
    return(last_year >= end_yr)
  }, error = function(e) FALSE)
}

# --- Helper: Clean Invalid Soil Files ---
clean_invalid_soils <- function(dir, ids) {
  if (!dir.exists(dir)) return()
  for (id in ids) {
    f <- file.path(dir, paste0(id, ".SOL"))
    if (file.exists(f)) {
      is_valid <- tryCatch({
        if (file.info(f)$size == 0) FALSE
        else {
          lines <- readLines(f, n = 5, warn = FALSE, encoding = "UTF-8")
          if (length(lines) == 0) FALSE
          else if (any(grepl("ERROR", lines, ignore.case = TRUE))) FALSE
          else if (!any(grepl("^\\s*\\*", lines))) FALSE
          else TRUE
        }
      }, error = function(e) FALSE)
      
      if (!is_valid) {
        message(sprintf("Deleting invalid/empty soil file: %s", f))
        unlink(f)
      }
    }
  }
}

clear_run_diagnostics <- function(ids) {
  artifacts <- c("_run_error.log", "dssat_B_stdout_stderr.log", "dssat_Q_stdout_stderr.log",
                 "ERROR.OUT", "WARNING.OUT", "INFO.OUT", "Summary.OUT", "summary.csv")
  for (id in ids) {
    unlink(file.path(id, artifacts), force = TRUE)
  }
}

write_input_error <- function(id, reason) {
  dir.create(id, showWarnings = FALSE, recursive = TRUE)
  line <- sprintf("[%s] ID %s: INPUT: %s", format(Sys.time()), id, reason)
  con <- file(file.path(id, "_run_error.log"), open = "w", encoding = "UTF-8")
  writeLines(line, con = con)
  close(con)
}

soil_input_issue <- function(id) {
  f <- file.path(id, "SOIL.SOL")
  if (!file.exists(f)) return("SOIL.SOL is missing")
  lines <- tryCatch(readLines(f, warn = FALSE, encoding = "UTF-8"), error = function(e) character())
  if (!length(lines)) return("SOIL.SOL is empty or unreadable")
  if (any(grepl("^\\*SOIL ERROR", lines))) {
    return(paste(trimws(lines[grepl("^\\*SOIL ERROR|^Source missing|^No Soil ID", lines)]),
                 collapse = " | "))
  }

  hdr <- grep("^@\\s+SLB\\b", lines)
  if (!length(hdr)) return("SOIL.SOL has no @ SLB layer table")
  layer_lines <- lines[(hdr[1] + 1):length(lines)]
  layer_lines <- layer_lines[nzchar(trimws(layer_lines))]
  layer_lines <- layer_lines[grepl("^\\s*\\d+\\s+", layer_lines)]
  depths <- suppressWarnings(as.integer(sub("^\\s*(\\d+).*", "\\1", layer_lines)))
  depths <- depths[is.finite(depths)]
  if (!length(depths)) return("SOIL.SOL has no parseable SLB layer depths")
  if (length(depths) > 19) return(sprintf("SOIL.SOL has %d layers; DSSAT accepts at most 19", length(depths)))
  if (any(diff(depths) <= 0)) return(sprintf("SOIL.SOL layer depths are not strictly increasing: %s",
                                             paste(depths, collapse = ",")))
  NULL
}

weather_input_issue <- function(id) {
  f <- file.path(id, paste0(id, ".WTH"))
  if (!file.exists(f)) return(paste0(basename(f), " is missing"))
  if (file.info(f)$size == 0) return(paste0(basename(f), " is empty"))
  NULL
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

# 3. Check soil input files
if (SOIL_SOURCE == "HWSD") {
  for (f in c(HWSD_RASTER_FILE, HWSD_DB_FILE)) {
    if (!file.exists(f))
      stop(sprintf("CRITICAL: HWSD file not found: %s\nDownload HWSD v2.0 from FAO and set hwsd_raster_file / hwsd_db_file in config.yml.", f))
  }
} else if (SOIL_SOURCE == "LUCAS") {
  if (!file.exists(LUCAS_CSV))
    stop(sprintf("CRITICAL: LUCAS topsoil table not found: %s\nRequest it free from ESDAC (esdac.jrc.ec.europa.eu) and set lucas_csv in config.yml.", LUCAS_CSV))
} else if (!(SOIL_SOURCE %in% KEYLESS_ONLINE_SOIL) && !file.exists(EXTERNAL_SOIL_FILE)) {
  stop(sprintf("CRITICAL: External soil file needed for %s but not found at: %s", SOIL_SOURCE, EXTERNAL_SOIL_FILE))
}

message("All checks passed. Starting pipeline...")

#-----------------------------------------------------------------------
# --- DEPRECATED (legacy): extend_weather_smart_single ---
# Superseded by extend_weather_repeat_single_ignore_partial.
# Retained for reference only. Do NOT call directly.
extend_weather_smart_single <- function(f, reference_year) {
  
  # --- STAGE 1: FAST READ ---
  lines <- readLines(f, encoding = "UTF-8")
  
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
  
  con <- file(f, open = "w", encoding = "UTF-8")
  writeLines(c(header_lines, formatted_body), con = con)
  close(con)
  return(TRUE)
}

resolve_optional_path <- function(path) {
  path <- trimws(as.character(if (is.null(path)) "" else path))
  if (!nzchar(path)) return("")
  if (grepl("^([A-Za-z]:|/|\\\\\\\\|~)", path)) return(normalizePath(path, mustWork = FALSE))
  normalizePath(file.path(INPUT_ROOT_DIR, path), mustWork = FALSE)
}

apply_cropland_mask <- function(points_sf, raster_file, grid_spacing_m, classes,
                                min_fraction = 0, strict = FALSE) {
  fallback <- function(text) {
    if (strict) stop(text, call. = FALSE)
    message("WARNING: ", text, " Continuing with all grid cells.")
    points_sf
  }

  raster_file <- resolve_optional_path(raster_file)
  if (!nzchar(raster_file) || !file.exists(raster_file)) {
    return(fallback(sprintf("Cropland mask enabled, but cropland_raster_file is missing or not found: '%s'", raster_file)))
  }
  if (!requireNamespace("terra", quietly = TRUE)) {
    return(fallback("Cropland mask enabled, but package 'terra' is not installed."))
  }
  if (length(classes) == 0 || any(is.na(classes))) {
    return(fallback("Cropland mask enabled, but cropland_classes is empty or invalid."))
  }

  r <- tryCatch(terra::rast(raster_file), error = function(e) e)
  if (inherits(r, "error")) {
    return(fallback(sprintf("Could not read cropland raster '%s': %s", raster_file, r$message)))
  }
  if (terra::is.lonlat(r)) {
    return(fallback("Cropland raster is longitude/latitude. Reproject it to a meter-based CRS before grid-cell fraction extraction."))
  }

  message(sprintf("Applying cropland mask from %s (classes: %s)", raster_file, paste(classes, collapse = ", ")))
  pts_r <- st_transform(points_sf, terra::crs(r))
  xy <- st_coordinates(pts_r)
  half <- grid_spacing_m / 2
  frac <- vapply(seq_len(nrow(pts_r)), function(i) {
    x <- xy[i, 1]; y <- xy[i, 2]
    cell_ext <- terra::ext(x - half, x + half, y - half, y + half)
    cell_r <- tryCatch(terra::crop(r, cell_ext, snap = "out"), error = function(e) NULL)
    if (is.null(cell_r) || terra::ncell(cell_r) == 0) return(NA_real_)
    vals <- terra::values(cell_r, mat = FALSE, na.rm = TRUE)
    if (length(vals) == 0) return(NA_real_)
    mean(vals %in% classes)
  }, numeric(1))
  frac[is.nan(frac)] <- NA_real_
  frac <- pmax(0, pmin(1, frac))

  out <- points_sf
  out$crop_frac <- frac
  out$crop_pct <- round(100 * frac, 4)
  out$cell_ha <- (grid_spacing_m^2) / 10000
  out$crop_ha <- round(out$cell_ha * frac, 4)

  keep <- !is.na(out$crop_frac)
  if (min_fraction <= 0) keep <- keep & out$crop_frac > 0 else keep <- keep & out$crop_frac >= min_fraction
  kept <- sum(keep)
  message(sprintf("Cropland mask retained %d of %d grid cells (%.1f%%).",
                  kept, nrow(out), if (nrow(out) > 0) 100 * kept / nrow(out) else 0))
  if (kept == 0) stop("Cropland mask removed all grid cells. Lower cropland_min_fraction or check cropland raster/classes.", call. = FALSE)

  out <- out[keep, ]
  tryCatch(
    sf::st_write(out, POINT_SHAPEFILE_PATH, append = FALSE, delete_layer = TRUE, quiet = TRUE),
    error = function(e) message("WARNING: Could not rewrite cropland-filtered grid shapefile: ", e$message)
  )
  out
}

#-----------------------------------------------------------------------
# STEP 0: CREATE GRIDFILE
#-----------------------------------------------------------------------
# STEP 0: CREATE GRIDFILE
#-----------------------------------------------------------------------
message("STEP 0: PREPARING GRIDFILE / POINTS")
setwd(OUTPUT_ROOT_DIR)
dir.create(GRIDPOINTS_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!USE_EXISTING_POINT_SHAPEFILE && USE_CROPLAND_MASK && REUSE_CROPLAND_GRID && file.exists(POINT_SHAPEFILE_PATH)) {
  message(sprintf("Reusing existing cropland grid shapefile: %s", POINT_SHAPEFILE_PATH))
  gridfile <- st_read(POINT_SHAPEFILE_PATH, quiet = TRUE)
} else if (USE_EXISTING_POINT_SHAPEFILE) {
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
  raw_grid_path <- if (USE_CROPLAND_MASK) ALL_LAND_POINT_SHAPEFILE_PATH else POINT_SHAPEFILE_PATH
  gridfile <- create_grid_points(boundary_sf, GRID_SPACING_METERS, raw_grid_path)
}

if (USE_CROPLAND_MASK && !("crop_pct" %in% names(gridfile))) {
  gridfile <- apply_cropland_mask(
    points_sf = gridfile,
    raster_file = CROPLAND_RASTER_FILE,
    grid_spacing_m = GRID_SPACING_METERS,
    classes = CROPLAND_CLASSES,
    min_fraction = CROPLAND_MIN_FRACTION,
    strict = CROPLAND_STRICT
  )
}

message(sprintf("Grid points ready: %d point(s)", nrow(gridfile)))

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
  
  if (SOIL_SOURCE %in% PER_POINT_SOIL) {
    # All of these write one Saxton-&-Rawls .SOL per point named by point ID, so
    # the combine/retry logic is shared. gNATSGO fills SSURGO's US gaps;
    # iSDAsoil = Africa 30 m; LUCAS = EU topsoil (needs lucas_csv).
    .soil_fn <- switch(SOIL_SOURCE,
      SSURGO   = process_soils_ssurgo,
      SSURGO_ALDERMAN = function(g, c, i, n, id, la, lo, fm)
                          process_soils_ssurgo_alderman(g, c, i, n, id, la, lo, fm, STATSGO = STATSGO, standardize_layers = STANDARDIZE_LAYERS),
      GNATSGO  = process_soils_gnatsgo,
      ISDASOIL = process_soils_isdasoil,
      LUCAS    = function(gridfile, csv, indiv, cores, idc, latc, lonc, fmt)
                   process_soils_lucas(gridfile, csv, indiv, cores, idc, latc, lonc, fmt, lucas_csv = LUCAS_CSV))
    if (CHECK_SOIL_DOWNLOADS) {
      ids <- as.character(gridfile[[POINT_ID_COLUMN]])
      clean_invalid_soils(individual_sol_output_folder, ids)
    }
    
    max_retries <- if (CHECK_SOIL_DOWNLOADS) SOIL_DOWNLOAD_RETRIES else 1
    retry_count <- 0
    ids <- as.character(gridfile[[POINT_ID_COLUMN]])
    
    # Define combine helper
    combine_sol_files_local <- function(input_folder, output_file_path) {
      sol_files <- list.files(path = input_folder, pattern = "\\.SOL$", full.names = TRUE)
      out_con <- file(output_file_path, open = "wt", encoding = "UTF-8")
      cat("*SOILS: Combined\n", file = out_con)
      for (f in sol_files) {
        lines <- readLines(f, warn = FALSE, encoding = "UTF-8")
        start <- grep("^\\*", lines)[1]
        if(!is.na(start)) writeLines(lines[start:length(lines)], out_con)
      }
      close(out_con)
    }

    repeat {
      .soil_fn(gridfile, soilfile_CSV, individual_sol_output_folder, SOIL_CORES,
               POINT_ID_COLUMN, LAT_COLUMN, LONG_COLUMN, format_SQL_in_statement)

      combine_sol_files_local(individual_sol_output_folder, soilfile_DSSAT)
      
      retry_count <- retry_count + 1
      
      if (CHECK_SOIL_DOWNLOADS && retry_count < max_retries) {
        clean_invalid_soils(individual_sol_output_folder, ids)
        existing_files <- tools::file_path_sans_ext(list.files(individual_sol_output_folder, pattern = "\\.SOL$"))
        missing_count <- sum(!ids %in% existing_files)
        if (missing_count == 0) {
          break
        } else {
          message(sprintf("%s: %d profiles missing/invalid. Retrying...", SOIL_SOURCE, missing_count))
        }
      } else {
        break
      }
    }
    
  } else if (SOIL_SOURCE %in% c("SOILGRIDS_10K", "AGMIP")) {
    if (!file.exists(EXTERNAL_SOIL_FILE)) stop("External soil file not found.")
    .external_sol_fn <- if (SOIL_SOURCE == "AGMIP") process_soils_agmip else process_soils_soilgrids
    .external_sol_fn(gridfile, EXTERNAL_SOIL_FILE, soilfile_CSV,
                     output_sol_dir = individual_sol_output_folder,
                     id_col = POINT_ID_COLUMN)
  } else if (SOIL_SOURCE == "SOILGRIDS_ONLINE") {
    ids <- as.character(gridfile[[POINT_ID_COLUMN]])
    if (CHECK_SOIL_DOWNLOADS) {
      clean_invalid_soils(individual_sol_output_folder, ids)
    }
    existing_files <- tools::file_path_sans_ext(list.files(individual_sol_output_folder, pattern = "\\.SOL$"))
    missing_mask <- ! (ids %in% existing_files)
    points_to_process <- gridfile[missing_mask, ]
    
    if (nrow(points_to_process) == 0) {
      message("All online soil profiles already exist. Skipping SoilGrids processing.")
    } else {
      USE_REST_API <<- (SOILGRIDS_MODE != "VRT")
      message(sprintf("SoilGrids online mode: %s. Processing %d missing points...",
                      if (SOILGRIDS_MODE == "VRT") "VRT" else "REST API", nrow(points_to_process)))
                      
      max_retries <- if (CHECK_SOIL_DOWNLOADS) SOIL_DOWNLOAD_RETRIES else 1
      retry_count <- 0
      
      while (nrow(points_to_process) > 0 && retry_count < max_retries) {
        if (retry_count > 0) {
          message(sprintf("Retry %d/%d: Fetching %d failed online soil points...",
                          retry_count, max_retries, nrow(points_to_process)))
        }
        
        ok <- tryCatch({
          process_soils_soilgrids_online(points_to_process, soilfile_CSV,
                                         output_sol_dir = individual_sol_output_folder,
                                         id_col = POINT_ID_COLUMN)
          TRUE
        }, error = function(e) {
          warning(sprintf("SoilGrids online extraction failed for this retry batch (%d point(s)): %s",
                          nrow(points_to_process), conditionMessage(e)))
          FALSE
        })
        
        retry_count <- retry_count + 1
        
        if (CHECK_SOIL_DOWNLOADS) {
          ids_left <- as.character(points_to_process[[POINT_ID_COLUMN]])
          clean_invalid_soils(individual_sol_output_folder, ids_left)
          if (!ok) break
          existing_left <- tools::file_path_sans_ext(list.files(individual_sol_output_folder, pattern = "\\.SOL$"))
          missing_left_mask <- ! (ids_left %in% existing_left)
          points_to_process <- points_to_process[missing_left_mask, ]
        } else {
          # No deep validity check requested — still confirm the .SOL files now
          # exist on disk so successfully-fetched points are cleared and we don't
          # emit a spurious "failed" warning for them.
          existing_left <- tools::file_path_sans_ext(list.files(individual_sol_output_folder, pattern = "\\.SOL$"))
          ids_left <- as.character(points_to_process[[POINT_ID_COLUMN]])
          points_to_process <- points_to_process[!(ids_left %in% existing_left), ]
          break
        }
      }
      
      if (nrow(points_to_process) > 0) {
        warning(sprintf("Failed to successfully download soil data for %d points after %d retries.",
                        nrow(points_to_process), max_retries))
      }
    }

    valid_online_ids <- tools::file_path_sans_ext(list.files(individual_sol_output_folder, pattern = "\\.SOL$"))
    valid_online_ids <- intersect(ids, valid_online_ids)
    missing_online_ids <- setdiff(ids, valid_online_ids)
    if (length(valid_online_ids) == 0) {
      stop("SOILGRIDS_ONLINE produced no valid .SOL files. Check SoilGrids connectivity, coordinates, or try SOILGRIDS_10K / SSURGO.",
           call. = FALSE)
    }
    if (length(missing_online_ids) > 0) {
      message(sprintf("SOILGRIDS_ONLINE: %d/%d valid profile(s) available; missing IDs will be skipped during DSSAT execution: %s",
                      length(valid_online_ids), length(ids),
                      paste(head(missing_online_ids, 20), collapse = ", ")))
    }
    
    # Always write/rebuild the complete CSV mapping since we might have run a subset
    mapping_df <- data.frame(ID = gridfile[[POINT_ID_COLUMN]], SOIL_ID = gridfile[[POINT_ID_COLUMN]])
    colnames(mapping_df)[1] <- POINT_ID_COLUMN
    write.csv(mapping_df, soilfile_CSV, row.names = FALSE)

  } else if (SOIL_SOURCE == "POLARIS") {
    # POLARIS = 30 m probabilistic disaggregation of SSURGO (CONUS). Same output
    # contract as SOILGRIDS_ONLINE: one .SOL per point named by point ID + a CSV,
    # streamed via GDAL /vsicurl. Water limits come from POLARIS's van Genuchten
    # curve (stat=p50 deterministic). Resume by skipping points already on disk.
    ids <- as.character(gridfile[[POINT_ID_COLUMN]])
    existing_files <- tools::file_path_sans_ext(
      list.files(individual_sol_output_folder, pattern = "\\.SOL$"))
    points_to_process <- gridfile[!(ids %in% existing_files), ]

    if (nrow(points_to_process) == 0) {
      message("All POLARIS soil profiles already exist. Skipping POLARIS processing.")
    } else {
      message(sprintf("POLARIS (statistic=%s): processing %d missing point(s)...",
                      POLARIS_STAT, nrow(points_to_process)))
      tryCatch({
        process_soils_polaris(points_to_process, soilfile_CSV,
                              output_sol_dir = individual_sol_output_folder,
                              id_col = POINT_ID_COLUMN,
                              stat = POLARIS_STAT, cache_dir = POLARIS_CACHE_DIR)
      }, error = function(e) {
        warning(sprintf("POLARIS extraction failed for this batch (%d point(s)): %s",
                        nrow(points_to_process), conditionMessage(e)))
      })
    }

    valid_ids <- intersect(ids, tools::file_path_sans_ext(
      list.files(individual_sol_output_folder, pattern = "\\.SOL$")))
    if (length(valid_ids) == 0) {
      stop("POLARIS produced no valid .SOL files. Check connectivity, coordinates (CONUS only), or try SSURGO.",
           call. = FALSE)
    }
    missing_ids <- setdiff(ids, valid_ids)
    if (length(missing_ids) > 0) {
      message(sprintf("POLARIS: %d/%d valid profile(s) available; missing IDs will be skipped: %s",
                      length(valid_ids), length(ids),
                      paste(head(missing_ids, 20), collapse = ", ")))
    }
    mapping_df <- data.frame(ID = gridfile[[POINT_ID_COLUMN]], SOIL_ID = gridfile[[POINT_ID_COLUMN]])
    colnames(mapping_df)[1] <- POINT_ID_COLUMN
    write.csv(mapping_df, soilfile_CSV, row.names = FALSE)

  } else if (SOIL_SOURCE == "HWSD") {
    process_soils_hwsd(gridfile, HWSD_RASTER_FILE, HWSD_DB_FILE,
                       output_csv_path = soilfile_CSV,
                       output_sol_dir = individual_sol_output_folder,
                       id_col = POINT_ID_COLUMN,
                       lat_col = LAT_COLUMN, long_col = LONG_COLUMN)
  } else if (SOIL_SOURCE == "HIHYDROSOIL") {
    process_soils_hihydrosoil(gridfile, HIHYDROSOIL_RASTER_DIR,
                              output_csv_path = soilfile_CSV,
                              output_sol_dir = individual_sol_output_folder,
                              id_col = POINT_ID_COLUMN,
                              lat_col = LAT_COLUMN, long_col = LONG_COLUMN)
  } else if (SOIL_SOURCE == "SLGA") {
    process_soils_slga(gridfile, SLGA_RASTER_DIR,
                       output_csv_path = soilfile_CSV,
                       output_sol_dir = individual_sol_output_folder,
                       id_col = POINT_ID_COLUMN,
                       lat_col = LAT_COLUMN, long_col = LONG_COLUMN)
  } else if (SOIL_SOURCE == "WISE30SEC") {
    process_soils_wise30sec(gridfile, WISE30SEC_RASTER_DIR,
                            output_csv_path = soilfile_CSV,
                            output_sol_dir = individual_sol_output_folder,
                            id_col = POINT_ID_COLUMN,
                            lat_col = LAT_COLUMN, long_col = LONG_COLUMN)
  } else if (SOIL_SOURCE == "WOSIS") {
    process_soils_wosis(gridfile, WOSIS_PROFILE_CSV,
                        output_csv_path = soilfile_CSV,
                        output_sol_dir = individual_sol_output_folder,
                        id_col = POINT_ID_COLUMN,
                        lat_col = LAT_COLUMN, long_col = LONG_COLUMN)
  } else { stop("Unknown SOIL_SOURCE") }
}

#-----------------------------------------------------------------------
# STEP 2: WEATHER DATA
#-----------------------------------------------------------------------
message("STEP 2: DOWNLOADING WEATHER DATA")
output_dir <- file.path(CENTRAL_WEATHER_DIR, WEATHER_DIR_NAME)
if (RUN_STEP_2_WEATHER) {
  dir.create(CENTRAL_WEATHER_DIR, showWarnings = FALSE, recursive = TRUE)
  
  # Uses new naming: [Location]_[Res]_[Source]
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # --- SMART RESUME BLOCK ---
  ids <- as.character(gridfile[[POINT_ID_COLUMN]])
  if (CHECK_WEATHER_DOWNLOADS) {
    message("Verifying existing weather files for validity...")
    valid_mask <- vapply(ids, function(id) is_wth_valid(id, output_dir, WEATHER_END_YEAR), logical(1))
    missing_mask <- !valid_mask
  } else {
    existing_files <- tools::file_path_sans_ext(list.files(output_dir, pattern = "\\.WTH$"))
    missing_mask <- ! (ids %in% existing_files)
  }
  points_to_process <- gridfile[missing_mask, ]
  
  if (nrow(points_to_process) == 0) {
    message("All weather files already exist and are valid. Skipping processing.")
  } else {
    log_file <- file.path(output_dir, "download_errors.log")
    
    max_retries <- if (CHECK_WEATHER_DOWNLOADS) WEATHER_DOWNLOAD_RETRIES else 1
    retry_count <- 0
    
    while (nrow(points_to_process) > 0 && retry_count < max_retries) {
      if (retry_count > 0) {
        message(sprintf("Retry %d/%d: Downloading %d failed/incomplete weather points...", 
                        retry_count, max_retries, nrow(points_to_process)))
      } else {
        message(sprintf("Processing %d weather points...", nrow(points_to_process)))
      }
      
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
      else if (WEATHER_SOURCE == "OPEN_METEO") do.call(process_weather_openmeteo, common_args)
      else if (WEATHER_SOURCE == "NASA_POWER_CHIRPS") {
        CHIRPS_RESOLUTION <<- if (exists("CHIRPS_RESOLUTION")) CHIRPS_RESOLUTION else "p05"
        do.call(process_weather_nasapower_chirps,
                c(common_args, list(chirps_cache_dir = CHIRPS_CACHE_DIR)))
      }
      else if (WEATHER_SOURCE == "AGERA5")
        {
          agera5_args <- common_args
          agera5_args$n_cores <- AGERA5_MAX_CONCURRENT_REQUESTS
          do.call(process_weather_agera5, c(agera5_args, list(agera5_cache_dir = AGERA5_CACHE_DIR)))
        }
      else if (WEATHER_SOURCE == "DWD")
        do.call(process_weather_dwd, c(common_args, list(dwd_cache_dir = DWD_CACHE_DIR)))
      else if (WEATHER_SOURCE == "EOBS")
        do.call(process_weather_eobs, c(common_args, list(eobs_nc_dir = EOBS_NC_DIR,
                                                          eobs_cache_dir = EOBS_CACHE_DIR,
                                                          eobs_use_cds = EOBS_USE_CDS)))
      else if (WEATHER_SOURCE == "XAVIER")
        do.call(process_weather_xavier, c(common_args, list(xavier_nc_dir = XAVIER_NC_DIR)))
      else if (WEATHER_SOURCE == "CMFD")
        do.call(process_weather_cmfd, c(common_args, list(cmfd_nc_dir = CMFD_NC_DIR)))
      else if (WEATHER_SOURCE == "CHELSA_W5E5")
        do.call(process_weather_chelsa_w5e5, c(common_args, list(chelsa_nc_dir = CHELSA_NC_DIR)))
      else if (WEATHER_SOURCE == "AGMERRA")
        do.call(process_weather_agmerra, c(common_args, list(agmerra_nc_dir = AGMERRA_NC_DIR)))
      else if (WEATHER_SOURCE == "AGCFSR")
        do.call(process_weather_agcfsr, c(common_args, list(agcfsr_nc_dir = AGCFSR_NC_DIR)))
      else if (WEATHER_SOURCE == "SILO")
        do.call(process_weather_silo, c(common_args, list(silo_nc_dir = SILO_NC_DIR)))
      else if (WEATHER_SOURCE == "PRISM")
        do.call(process_weather_prism, c(common_args, list(prism_cache_dir = PRISM_CACHE_DIR)))
      else if (WEATHER_SOURCE == "MSWX")
        do.call(process_weather_mswx, c(common_args, list(mswx_nc_dir = MSWX_NC_DIR)))
      else if (WEATHER_SOURCE == "MSWEP")
        do.call(process_weather_mswep, c(common_args, list(mswep_nc_dir = MSWEP_NC_DIR)))
      else if (WEATHER_SOURCE == "CRUJRA")
        do.call(process_weather_crujra, c(common_args, list(crujra_nc_dir = CRUJRA_NC_DIR)))
      else if (WEATHER_SOURCE == "TERRACLIMATE")
        do.call(process_weather_terraclimate, c(common_args, list(terraclimate_nc_dir = TERRACLIMATE_NC_DIR)))

      retry_count <- retry_count + 1
      
      if (CHECK_WEATHER_DOWNLOADS) {
        ids_left <- as.character(points_to_process[[POINT_ID_COLUMN]])
        valid_mask <- vapply(ids_left, function(id) is_wth_valid(id, output_dir, WEATHER_END_YEAR), logical(1))
        points_to_process <- points_to_process[!valid_mask, ]
      } else {
        # No deep validity check requested — still confirm the .WTH files now
        # exist on disk, so points that downloaded fine are cleared from the
        # pending list and we don't emit a spurious "failed" warning for them.
        existing_files <- tools::file_path_sans_ext(list.files(output_dir, pattern = "\\.WTH$"))
        ids_left <- as.character(points_to_process[[POINT_ID_COLUMN]])
        points_to_process <- points_to_process[!(ids_left %in% existing_files), ]
        break
      }
    }
    
    if (nrow(points_to_process) > 0) {
      warning(sprintf("Failed to successfully download weather data for %d points after %d retries.", 
                      nrow(points_to_process), max_retries))
    }
  }
}

if (REPAIR_WEATHER_MISSING_VALUES) {
  if (dir.exists(output_dir)) {
    ids <- as.character(gridfile[[POINT_ID_COLUMN]])
    repair_log <- file.path(output_dir, "weather_repair.log")
    repair_summary <- dssatutils::repair_weather_missing_values(
      weather_dir = output_dir,
      ids = ids,
      max_gap_days = WEATHER_REPAIR_MAX_GAP_DAYS,
      window_days = WEATHER_REPAIR_WINDOW_DAYS,
      variables = WEATHER_REPAIR_VARIABLES,
      log_file = repair_log
    )
    repaired_values <- sum(repair_summary$repaired_count, na.rm = TRUE)
    unrepaired_values <- sum(repair_summary$unrepaired_count, na.rm = TRUE)
    message(sprintf(
      "Weather missing-value repair complete: %d value(s) repaired; %d missing value(s) left unrepaired. Log: %s",
      repaired_values, unrepaired_values, repair_log
    ))
  } else {
    warning(sprintf(
      "repair_weather_missing_values is TRUE, but weather directory does not exist: %s",
      output_dir
    ))
  }
}

if (REPAIR_WEATHER_DATE_GAPS) {
  if (dir.exists(output_dir)) {
    ids <- as.character(gridfile[[POINT_ID_COLUMN]])
    repair_log <- file.path(output_dir, "weather_repair.log")
    repair_summary <- dssatutils::repair_weather_date_gaps(
      weather_dir = output_dir,
      ids = ids,
      max_gap_days = WEATHER_REPAIR_MAX_GAP_DAYS,
      window_days = WEATHER_REPAIR_WINDOW_DAYS,
      variables = WEATHER_REPAIR_VARIABLES,
      log_file = repair_log
    )
    repaired_values <- sum(repair_summary$repaired_count, na.rm = TRUE)
    unrepaired_values <- sum(repair_summary$unrepaired_count, na.rm = TRUE)
    message(sprintf(
      "Weather date-gap repair complete: %d missing day row(s) inserted; %d missing day row(s) left unrepaired. Log: %s",
      repaired_values, unrepaired_values, repair_log
    ))
  } else {
    warning(sprintf(
      "repair_weather_date_gaps is TRUE, but weather directory does not exist: %s",
      output_dir
    ))
  }
}

if (REPAIR_WEATHER_TEMPERATURE_INVERSIONS) {
  if (dir.exists(output_dir)) {
    ids <- as.character(gridfile[[POINT_ID_COLUMN]])
    repair_log <- file.path(output_dir, "weather_repair.log")
    repair_summary <- dssatutils::repair_weather_temperature_inversions(
      weather_dir = output_dir,
      ids = ids,
      max_gap_days = WEATHER_REPAIR_MAX_GAP_DAYS,
      window_days = WEATHER_REPAIR_WINDOW_DAYS,
      log_file = repair_log
    )
    repaired_values <- sum(repair_summary$repaired_count, na.rm = TRUE)
    unrepaired_values <- sum(repair_summary$unrepaired_count, na.rm = TRUE)
    message(sprintf(
      "Weather Tmax/Tmin inversion repair complete: %d day(s) repaired; %d inversion day(s) left unrepaired. Log: %s",
      repaired_values, unrepaired_values, repair_log
    ))
  } else {
    warning(sprintf(
      "repair_weather_temperature_inversions is TRUE, but weather directory does not exist: %s",
      output_dir
    ))
  }
}

if (AUDIT_WEATHER_QUALITY) {
  if (dir.exists(output_dir)) {
    ids <- as.character(gridfile[[POINT_ID_COLUMN]])
    repair_log <- file.path(output_dir, "weather_repair.log")
    audit_csv <- file.path(output_dir, "weather_quality_audit.csv")
    audit_summary <- dssatutils::audit_weather_quality(
      weather_dir = output_dir,
      ids = ids,
      audit_csv = audit_csv,
      flatline_days = WEATHER_QUALITY_FLATLINE_DAYS,
      log_file = repair_log
    )
    message(sprintf(
      "Weather QA audit complete: %d finding row(s). Audit CSV: %s",
      nrow(audit_summary), audit_csv
    ))
  } else {
    warning(sprintf(
      "audit_weather_quality is TRUE, but weather directory does not exist: %s",
      output_dir
    ))
  }
}

# === WEATHER EXTENSION LOGIC (PARALLEL) ===
if (EXTEND_WEATHER_DATA) {
  if (dir.exists(output_dir)) {
    all_wth_files <- list.files(output_dir, pattern = "\\.WTH$", full.names = TRUE)
    message("Checking which files need extension...")
    
    files_to_extend <- Filter(function(f) {
      tryCatch({
        last_line <- tail(readLines(f, n = -1, encoding = "UTF-8"), 1)
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
  } else {
    warning(sprintf(
      "extend_weather_data is TRUE, but weather directory does not exist: %s",
      output_dir
    ))
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
  # The soil mapping CSV is a per-LAYER table (many rows per point id). We only
  # need point-level attributes (SOIL_ID, coords) here — the full profile lives
  # in each point's .SOL file — so collapse to ONE row per id before joining.
  # Without this, left_join() does a one-to-many expansion that multiplies
  # `points` into one-row-per-soil-layer (e.g. 23 points -> 317 rows), which
  # then runs DSSAT many redundant times per folder in parallel (racing on the
  # same dir) and breaks the per-point coordinate lookup.
  if (anyDuplicated(soil_map[[POINT_ID_COLUMN]])) {
    n_before <- nrow(soil_map)
    soil_map <- soil_map[!duplicated(soil_map[[POINT_ID_COLUMN]]), , drop = FALSE]
    message(sprintf("Soil mapping has %d layer-rows for %d points; collapsed to one row per point.",
                    n_before, nrow(soil_map)))
  }
  points <- left_join(points, soil_map, by = POINT_ID_COLUMN)
} else {
  message("WARNING: Soil mapping CSV not found. Attempting legacy logic.")
  if (SOIL_SOURCE %in% PER_POINT_SOIL) points$SOIL_ID <- points[[POINT_ID_COLUMN]]
}

# --- 3.1b. CRITICAL FIX FOR SSURGO / MISSING COLUMNS ---
if (!("SOIL_ID" %in% names(points))) {
  if (SOIL_SOURCE %in% PER_POINT_SOIL) {
    message(sprintf("%s Mode: Defaulting SOIL_ID to Grid Point ID.", SOIL_SOURCE))
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
    con <- file(file.path(ID, "SOIL.SOL"), open = "w", encoding = "UTF-8")
    writeLines(c("*SOIL ERROR", "No Soil ID assigned"), con = con)
    close(con)
    return(NULL)
  }
  
  if (SOIL_SOURCE %in% PER_POINT_SOIL) {
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
  else {
    con <- file(dest_path, open = "w", encoding = "UTF-8")
    writeLines(c("*SOIL ERROR", paste("Source missing:", source_filename)), con = con)
    close(con)
  }
  
  # 2. HANDLE EXPERIMENT FILE
  tryCatch({
    content <- readLines(TEMPLATE_FILE_PATH, encoding = "UTF-8")
    # --- Soil ID substitution ---
    # Preferred: the unambiguous 8-char SID00000 placeholder that occupies ONLY the
    # *FIELDS ID_SOIL column. Legacy fallback: the older SOIL_ID / ID_SOIL tokens.
    # (Mirrors the Python engine so R/Python output stays in parity - AGENTS.md S5.)
    if (any(grepl("SID00000", content, fixed=TRUE))) {
      content <- gsub("SID00000", hmx_replacement_id, content, fixed = TRUE)
    } else if (!any(grepl(TEMPLATE_SOIL_ID_PLACEHOLDER, content, fixed=TRUE))) {
      if(any(grepl("ID_SOIL", content, fixed=TRUE))) {
        if (any(grepl("   ID_SOIL", content, fixed=TRUE))) {
          content <- gsub("   ID_SOIL", sprintf("%-10s", hmx_replacement_id), content, fixed = TRUE)
        } else {
          content <- gsub("ID_SOIL", hmx_replacement_id, content, fixed = TRUE)
        }
      }
    } else {
      if (any(grepl("   SOIL_ID", content, fixed=TRUE))) {
        content <- gsub("   SOIL_ID", sprintf("%-10s", hmx_replacement_id), content, fixed = TRUE)
      } else {
        content <- gsub(TEMPLATE_SOIL_ID_PLACEHOLDER, hmx_replacement_id, content, fixed = TRUE)
      }
    }
    # --- Weather station (WSTA) substitution ---
    # Preferred: the unambiguous WID00000 placeholder. Legacy fallback: the blanket
    # 00000000 -> point ID (kept for un-migrated templates). With WID00000 present
    # the blanket is skipped, so ID_FIELD's own 00000000 is preserved exactly as in
    # the Python engine.
    if (any(grepl("WID00000", content, fixed=TRUE))) {
      content <- gsub("WID00000", ID, content, fixed = TRUE)
    } else {
      content <- gsub("00000000", ID, content, fixed = TRUE)
    }
    content <- gsub(tools::file_path_sans_ext(TEMPLATE_FILE_NAME), ID, content, fixed = TRUE)

    # Substitute real per-point coordinates into the *FIELDS tier-2 line so DSSAT
    # does not log "Error reading latitude/longitude/elevation". The template
    # placeholders (LATITUDE/LONGITUDE/ELEV) are right-justified at their field
    # edges, so each is swapped for a SAME-LENGTH number to preserve the
    # fixed-width column alignment. Only the single data line that carries the
    # LATITUDE/LONGITUDE placeholders is edited — "ELEV" also appears in the @L
    # header, so a global replace would corrupt it. Anything unexpected (bad
    # coord, width overflow) leaves the placeholder untouched: never fatal.
    content <- tryCatch({
      fld_idx <- which(grepl("LATITUDE", content, fixed = TRUE) &
                       grepl("LONGITUDE", content, fixed = TRUE))
      if (length(fld_idx) == 1) {
        lat  <- suppressWarnings(as.numeric(points[[LAT_COLUMN]][i]))
        lon  <- suppressWarnings(as.numeric(points[[LONG_COLUMN]][i]))
        elev <- -99  # DSSAT "missing"; valid number so no read error
        wth_src <- file.path(weather_repo, paste0(ID, ".WTH"))
        if (file.exists(wth_src)) {
          hdr <- readLines(wth_src, n = 4, warn = FALSE, encoding = "UTF-8")
          h_i <- grep("\\bELEV\\b", hdr)
          if (length(h_i) >= 1 && h_i[1] < length(hdr)) {
            nm <- strsplit(trimws(sub("^@", "", hdr[h_i[1]])), "\\s+")[[1]]
            vl <- strsplit(trimws(hdr[h_i[1] + 1]), "\\s+")[[1]]
            ec <- which(nm == "ELEV")
            if (length(ec) == 1 && ec <= length(vl)) {
              e <- suppressWarnings(as.numeric(vl[ec]))
              if (is.finite(e) && e > -90) elev <- e
            }
          }
        }
        # format to EXACTLY the placeholder width, else skip (keep alignment)
        fit <- function(x, width, digits) {
          s <- formatC(x, format = "f", digits = digits, width = width)
          if (nchar(s) == width) s else NULL
        }
        line  <- content[fld_idx]
        s_lat <- if (is.finite(lat)) fit(lat, 8, 3) else NULL  # "LATITUDE"  = 8
        s_lon <- if (is.finite(lon)) fit(lon, 9, 3) else NULL  # "LONGITUDE" = 9
        s_ele <- fit(elev, 4, 0)                               # "ELEV"      = 4
        if (!is.null(s_lat)) line <- sub("LATITUDE",  s_lat, line, fixed = TRUE)
        if (!is.null(s_lon)) line <- sub("LONGITUDE", s_lon, line, fixed = TRUE)
        if (!is.null(s_ele)) line <- sub("ELEV",      s_ele, line, fixed = TRUE)
        content[fld_idx] <- line
      }
      content
    }, error = function(e) content)

    con <- file(file.path(ID, paste0(ID, ".", tools::file_ext(TEMPLATE_FILE_NAME))), open = "w", encoding = "UTF-8")
    writeLines(content, con = con)
    close(con)
  }, error = function(e) return(NULL))
  
  # 3. HANDLE WEATHER FILE
  wth <- file.path(weather_repo, paste0(ID, ".WTH"))
  if (file.exists(wth)) file.copy(wth, file.path(ID, basename(wth)))
  
  # 4. EXTRAS — genotype/SDA/CO2 files (only copy files listed in SUPPORT_FILES).
  if (length(SUPPORT_FILES) > 0) {
    files_to_copy <- file.path(TEMPLATE_DIR, SUPPORT_FILES)
    file.copy(files_to_copy, ID, overwrite = TRUE)
  }
}

# --- 3.4. Execute Folder Creation ---
message("Creating/verifying simulation folders...")
ids_to_create_indices <- 1:nrow(points)
if (RESUME_DSSAT_RUNS) {
  completed_mask <- file.exists(file.path(points[[POINT_ID_COLUMN]], paste0("results_", points[[POINT_ID_COLUMN]], ".csv")))
  ids_to_create_indices <- which(!completed_mask)
}
if (length(ids_to_create_indices) > 0) {
  cl <- makeCluster(SOIL_CORES)
  on.exit(tryCatch(stopCluster(cl), error = function(e) NULL), add = TRUE)
  clusterExport(cl, varlist = c("points", "SOIL_SOURCE", "PER_POINT_SOIL", "individual_soil_folder", "TEMPLATE_FILE_PATH", "create_folders_and_files", "TEMPLATE_FILE_NAME", "TEMPLATE_SOIL_ID_PLACEHOLDER", "POINT_ID_COLUMN", "weather_repo", "TEMPLATE_DIR", "SUPPORT_FILES", "LAT_COLUMN", "LONG_COLUMN"), envir = environment())
  parLapply(cl, ids_to_create_indices, create_folders_and_files)
  stopCluster(cl)
  on.exit()  # clear the guard — cluster stopped cleanly
}

# --- 3.5. Run Simulations ---
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
    clear_run_diagnostics(ids_to_run)

    input_ok <- vapply(ids_to_run, function(id) {
      issue <- soil_input_issue(id)
      if (is.null(issue)) issue <- weather_input_issue(id)
      if (!is.null(issue)) {
        write_input_error(id, issue)
        return(FALSE)
      }
      TRUE
    }, logical(1))

    skipped_input_ids <- ids_to_run[!input_ok]
    if (length(skipped_input_ids) > 0) {
      message(sprintf("Skipping %d point(s) with invalid DSSAT inputs; see _run_error.log in each folder.",
                      length(skipped_input_ids)))
      message("  Invalid input IDs: ", paste(head(skipped_input_ids, 20), collapse = ", "),
              if (length(skipped_input_ids) > 20) sprintf(" ... (+%d more)", length(skipped_input_ids) - 20) else "")
    }
    ids_to_run <- ids_to_run[input_ok]
  }

  if(length(ids_to_run) > 0) {
    cl <- makeCluster(DSSAT_CORES)
    on.exit(tryCatch(stopCluster(cl), error = function(e) NULL), add = TRUE)
    clusterEvalQ(cl, { library(dssatengine) })
    SUPPORTS_TREATMENT_LIST <- "treatment_list" %in% names(formals(dssatengine::run_simulation))
    if (length(TREATMENT_LIST) && !SUPPORTS_TREATMENT_LIST) {
      stop("This config sets treatment_list, but the installed dssatengine package ",
           "does not support explicit treatment lists. Reinstall/update dssatengine ",
           "before running; falling back to a contiguous range would change the run.",
           call. = FALSE)
    }
    run_simulation_wrapper <- function(ID) {
      args <- list(
        ID = ID,
        dssat_run_dir = DSSAT_RUN_DIR,
        crop_extension = CROP_EXTENSION,
        template_file_name = TEMPLATE_FILE_NAME,
        template_file_path = TEMPLATE_FILE_PATH,
        run_mode = RUN_MODE,
        treatment_start = TREATMENT_START,
        treatment_end = TREATMENT_END,
        sequence_start = SEQUENCE_START,
        sequence_end = SEQUENCE_END,
        weather_start_year = WEATHER_START_YEAR,
        weather_end_year = WEATHER_END_YEAR,
        dssat_exe_path = DSSAT_EXE_PATH,
        cleanup_run_folders = CLEANUP_RUN_FOLDERS,
        points_df = points
      )
      if (length(TREATMENT_LIST)) args$treatment_list <- TREATMENT_LIST
      do.call(dssatengine::run_simulation, args)
    }
    clusterExport(cl, c("DSSAT_RUN_DIR", "CROP_EXTENSION", "TEMPLATE_FILE_NAME",
                        "TEMPLATE_FILE_PATH", "RUN_MODE", "TREATMENT_START",
                        "TREATMENT_END", "TREATMENT_LIST", "SUPPORTS_TREATMENT_LIST",
                        "SEQUENCE_START", "SEQUENCE_END",
                        "WEATHER_START_YEAR", "WEATHER_END_YEAR", "DSSAT_EXE_PATH",
                        "CLEANUP_RUN_FOLDERS", "points"),
                  envir = globalenv())
    parLapply(cl, ids_to_run, run_simulation_wrapper)
    stopCluster(cl)
    on.exit()  # clear guard — cluster stopped cleanly

    # --- Loud per-point failure report ---
    # run_simulation() runs in parallel workers whose message()/error output is
    # invisible in the parent console, so a failed point leaves no trace except a
    # missing results_<ID>.csv. Surface that here instead of silently combining
    # nothing. Each failed folder also contains a '_run_error.log' with the cause.
    produced  <- file.exists(file.path(ids_to_run, paste0("results_", ids_to_run, ".csv")))
    failed_ids <- ids_to_run[!produced]
    if (length(failed_ids) > 0) {
      message(sprintf("Retrying %d DSSAT point(s) serially after parallel failure: %s",
                      length(failed_ids), paste(head(failed_ids, 20), collapse = ", ")),
              if (length(failed_ids) > 20) sprintf(" ... (+%d more)", length(failed_ids) - 20) else "")
      invisible(lapply(failed_ids, run_simulation_wrapper))
      produced <- file.exists(file.path(ids_to_run, paste0("results_", ids_to_run, ".csv")))
    }
    failed_ids <- ids_to_run[!produced]
    if (length(failed_ids) > 0) {
      message(sprintf("WARNING: %d of %d point(s) produced NO results_<ID>.csv.",
                      length(failed_ids), length(ids_to_run)))
      message("  Failed IDs: ", paste(head(failed_ids, 20), collapse = ", "),
              if (length(failed_ids) > 20) sprintf(" ... (+%d more)", length(failed_ids) - 20) else "")
      message("  Cause is logged in '_run_error.log' inside each failed folder under:\n    ", DSSAT_RUN_DIR)
    } else {
      message(sprintf("All %d point(s) produced results.", length(ids_to_run)))
    }
  }

  # --- 3.7. Combine Results ---
  if (!dir.exists(FINAL_OUTPUT_DIR)) dir.create(FINAL_OUTPUT_DIR, recursive = TRUE)

  combine_results <- function(folder_ids) {
    col_spec <- cols(point_id = col_character(), soil_profile_id = col_character(),
                     weather_station_id = col_character(), dssat_file_id = col_character(),
                     .default = col_guess())
    paths <- file.path(folder_ids, paste0("results_", folder_ids, ".csv"))
    exists_mask <- file.exists(paths)
    if (!any(exists_mask)) return(tibble::tibble())
    bind_rows(lapply(paths[exists_mask], function(p)
      tryCatch(read_csv(p, show_col_types = FALSE, col_types = col_spec),
               error = function(e) { message("WARNING: could not read ", p, ": ", e$message); NULL })
    ))
  }

  final_data <- combine_results(all_ids)
  if (!is.null(final_data) && nrow(final_data) > 0) {
    point_attrs <- st_drop_geometry(gridfile)
    keep_cols <- intersect(c(POINT_ID_COLUMN, "crop_frac", "crop_pct", "crop_ha", "cell_ha"), names(point_attrs))
    if (length(keep_cols) > 1) {
      point_attrs <- point_attrs[, keep_cols, drop = FALSE]
      names(point_attrs)[names(point_attrs) == POINT_ID_COLUMN] <- "point_id"
      point_attrs$point_id <- as.character(point_attrs$point_id)
      final_data <- dplyr::left_join(final_data, point_attrs, by = "point_id")
      if ("cell_ha" %in% names(final_data)) final_data$gridcell_area_ha <- final_data$cell_ha
      if ("crop_ha" %in% names(final_data)) final_data$cropland_ha <- final_data$crop_ha
      if (all(c("final_grain_kg_ha", "cropland_ha") %in% names(final_data))) {
        final_data$final_grain_production_kg <- final_data$final_grain_kg_ha * final_data$cropland_ha
      }
      if (all(c("top_weight_kg_ha", "cropland_ha") %in% names(final_data))) {
        final_data$top_weight_production_kg <- final_data$top_weight_kg_ha * final_data$cropland_ha
      }
    }
    write_csv(final_data, FINAL_RESULTS_PATH)
    message(sprintf("Results combined: %d rows -> %s", nrow(final_data), FINAL_RESULTS_PATH))
  } else {
    message("WARNING: No results were produced. Check DSSAT logs in the per-point run folders.")
  }

  # --- 3.8. Optional cleanup of per-point run folders --------------------------
  # When CLEANUP_RUN_FOLDERS is TRUE, delete the per-point simulation subfolders
  # (and all of their DSSAT I/O) now that the combined summary CSV has been
  # written — so a large sweep does not leave thousands of folders cluttering the
  # drive. Guarded on the combined CSV actually existing on disk: if the combine
  # step produced nothing, the folders are KEPT so the failure can be diagnosed.
  # The per-scenario results CSV (FINAL_RESULTS_PATH) lives under
  # RESULTS_ROOT_DIR, outside DSSAT_RUN_DIR, so it survives the cleanup.
  if (CLEANUP_RUN_FOLDERS) {
    if (file.exists(FINAL_RESULTS_PATH)) {
      point_dirs <- file.path(DSSAT_RUN_DIR, all_ids)
      point_dirs <- point_dirs[dir.exists(point_dirs)]
      if (length(point_dirs) > 0) {
        message(sprintf("Cleanup: deleting %d per-point run folder(s) under %s",
                        length(point_dirs), DSSAT_RUN_DIR))
        unlink(point_dirs, recursive = TRUE, force = TRUE)
      }
    } else {
      message("Cleanup requested, but the combined results CSV was not written — keeping run folders for troubleshooting.")
    }
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
    sprintf("Treatments:        %s", TREATMENT_RUN_LABEL),
    sprintf("Sequences:         %d - %d", SEQUENCE_START, SEQUENCE_END),
    sprintf("Template File:     %s", TEMPLATE_FILE_NAME),
    "",
    "--- PATHS ---",
    sprintf("Soil Map CSV:      %s", soil_mapping_file)
  )
  metadata_con <- file(metadata_file, open = "w", encoding = "UTF-8")
  writeLines(metadata_content, con = metadata_con)
  close(metadata_con)
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
      setwd(OUTPUT_ROOT_DIR)
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

  # Visualization is best-effort: the results CSV is already written and combined
  # by this point, so a plotting failure (bad geometry, missing boundary_sf, a
  # treatment with all-NA yields, ...) must NEVER abort the run or the sweep loop.
  # (Parity: the same guard wraps STEP 4 in dssat_main_pipeline.py - AGENTS S2.)
  tryCatch({
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
  }, error = function(e) {
    message(sprintf("Visualization skipped (non-fatal): %s", conditionMessage(e)))
  })
  message(sprintf("\n%s\nPIPELINE FINISHED!\n%s", paste(rep("*", 60), collapse = ""), paste(rep("*", 60), collapse = "")))
} else {
  message("Skipping Visualization (HPC Prep Mode)")
}
