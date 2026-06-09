#!/usr/bin/env python3
"""
dssat_main_pipeline.py
======================
Python port of dssat_main_pipeline.R

End-to-end spatial / gridded DSSAT crop modeling pipeline:
Step 0 — Create or load a grid / point shapefile
Step 1 — Download and format soil data
Step 2 — Download and format weather data
Step 3 — Build DSSAT run folders and (optionally) run simulations
Step 4 — Combine results and visualize

HOW TO USE (beginners)
1. Edit SECTION 0 below (paths, crop, weather/soil sources, years).
2. Set DSSAT_EXE in your shell first (once):
   export DSSAT_EXE=/full/path/to/dscsm048   # Linux / macOS
   set DSSAT_EXE=C:\\DSSAT48\\DSCSM048.EXE    # Windows
3. Run:
   python dssat_main_pipeline.py

THREE SPATIAL DOMAIN MODES (same as R pipeline)
MODE A — Regular grid from a boundary polygon (default demo)
MODE B — Your own point / polygon shapefile
MODE C — Cropland-only points from CDL / NLCD (build with landcover scripts,
         then feed shapefile here via MODE B)
"""

# =============================================================================
# STANDARD LIBRARY
# =============================================================================
import csv
import gc
import math
import multiprocessing
import os
import platform
import re
import shutil
import subprocess
import sys
import zipfile
from datetime import date, datetime
from pathlib import Path
from typing import Optional

# =============================================================================
# THIRD-PARTY (install: pip install geopandas pandas numpy matplotlib pyproj)
# =============================================================================
try:
    import geopandas as gpd
    import numpy as np
    import pandas as pd
    from pyproj import CRS, Transformer
    from shapely.geometry import Point
    from shapely.ops import transform as shp_transform
except ImportError as exc:
    sys.exit(
        f"Required package missing: {exc}\n"
        "Install with: pip install geopandas pandas numpy matplotlib pyproj shapely"
    )

# =============================================================================
# SECTION 0 — MASTER CONFIGURATION
# =============================================================================

# --- 0.1 Path & platform detection -----------------------------------------

def _detect_project_dir() -> str:
    """Return the directory that contains this script."""
    return str(Path(__file__).resolve().parent)

MAIN_PROJECT_DIR = _detect_project_dir()
PROJECT_ROOT     = MAIN_PROJECT_DIR

R_SCRIPTS_DIR    = os.path.join(PROJECT_ROOT, "r_scripts")   # unused in Python port
PY_SCRIPTS_DIR   = os.path.join(PROJECT_ROOT, "py_scripts")
SHAPEFILE_DIR    = os.path.join(PROJECT_ROOT, "shapefile")
GRIDPOINTS_DIR   = os.path.join(PROJECT_ROOT, "gridpoints")
WEATHER_ROOT_DIR = os.path.join(PROJECT_ROOT, "weather")
SOIL_ROOT_DIR    = os.path.join(PROJECT_ROOT, "soil")
RUNS_ROOT_DIR    = os.path.join(PROJECT_ROOT, "dssat_runs")
RESULTS_ROOT_DIR = os.path.join(PROJECT_ROOT, "results")

print(f"Running in Project Directory: {MAIN_PROJECT_DIR}")

# --- 0.2 DSSAT executable ---------------------------------------------------
_os = platform.system()
if _os == "Windows":
    DSSAT_BASE     = r"C:\DSSAT48"
    DSSAT_EXE_NAME = "DSCSM048.EXE"
else:
    # macOS / Linux: compile from source — see README
    DSSAT_BASE     = os.path.expanduser("~/Documents/GitHub/DSSAT48")
    DSSAT_EXE_NAME = "dscsm048"

DSSAT_BASE = os.environ.get("DSSAT_BASE", DSSAT_BASE)

# --- Shared config overlay (config.yml) -------------------------------------
# Loads config.yml (if present) so R and Python share one set of settings.
# cfg_get(key, default) returns the YAML value or the default below.
try:
    import sys as _sys
    _sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from config_loader import cfg_get
except Exception:  # noqa: BLE001
    def cfg_get(key, default):
        return default

# --- 0.3 Project settings ---------------------------------------------------
PROJECT_NAME        = cfg_get("project_name", "dssat_spatial_demo")
GRID_SPACING_METERS = cfg_get("grid_spacing_meters", 50_000)   # 50 km test; 5–10 km production
CROP_EXTENSION      = cfg_get("crop_extension", "MZ")     # MZ=maize  WH=wheat  SB=soybean …

# --- 0.3b Optional run-folder naming ----------------------------------------
RUN_TAG           = cfg_get("run_tag", "")        # e.g. "run1" / "calibA"
RUN_NAME_STYLE    = cfg_get("run_name_style", "grid")    # "grid" | "scenario"
RUN_NAME_OVERRIDE = cfg_get("run_name_override", "")        # If set, used verbatim as run-folder name

# --- 0.4 Spatial domain (choose MODE A / B / C) ----------------------------
USE_EXISTING_POINT_SHAPEFILE  = cfg_get("use_existing_point_shapefile", False)
EXISTING_POINT_SHAPEFILE_PATH = cfg_get("existing_point_shapefile_path",
                                        os.path.join(MAIN_PROJECT_DIR, "gridpoints", "my_points.shp"))

# MODE A boundary settings (ignored when USE_EXISTING_POINT_SHAPEFILE = True)
BOUNDARY_SHAPEFILE_NAME = cfg_get("boundary_shapefile_name", "tl_2024_us_state.shp")
ENABLE_BOUNDARY_FILTER  = cfg_get("enable_boundary_filter", True)
BOUNDARY_FILTER_COLUMN  = cfg_get("boundary_filter_column", "NAME")
STATE_NAME_FILTER       = list(cfg_get("state_name_filter", ["Iowa"]))

# --- 0.5 Auto-naming convention ---------------------------------------------
if GRID_SPACING_METERS < 1000:
    RESOLUTION_TAG = f"{GRID_SPACING_METERS}m"
else:
    RESOLUTION_TAG = f"{GRID_SPACING_METERS // 1000}km"

GRID_BASE_NAME        = f"{PROJECT_NAME}_{RESOLUTION_TAG}"
BOUNDARY_FILTER_VALUE = STATE_NAME_FILTER

# Weather settings
WEATHER_SOURCE     = cfg_get("weather_source", "GRIDMET")   # DAYMET | NASA_POWER | GRIDMET | OPEN_METEO | NASA_POWER_CHIRPS
WEATHER_START_YEAR = cfg_get("weather_start_year", 1982)
WEATHER_END_YEAR   = cfg_get("weather_end_year", 1983)
# NASA_POWER_CHIRPS only: CHIRPS rainfall resolution "p05" (~5.5 km, recommended)
# or "p25" (~28 km, lighter download).
CHIRPS_RESOLUTION  = str(cfg_get("chirps_resolution", "p05"))

# Soil settings
SOIL_SOURCE        = cfg_get("soil_source", "SSURGO")   # "SSURGO" | "SOILGRIDS_10K" | "SOILGRIDS_ONLINE"
EXTERNAL_SOIL_FILE = cfg_get("external_soil_file",
                             os.path.join(MAIN_PROJECT_DIR, "SoilGrids", "US.SOL"))
# SOILGRIDS_ONLINE only: "REST" (JSON API, rate-limited, no extra deps) or
# "VRT" (GDAL virtual rasters via rasterio; batch-friendly, better coverage).
SOILGRIDS_MODE     = str(cfg_get("soilgrids_mode", "REST")).upper()
# HWSD only: paths to the FAO HWSD v2.0 raster (SMU IDs) + SQLite database,
# downloaded once from FAO (blank = script defaults under HWSD/).
HWSD_RASTER_FILE   = cfg_get("hwsd_raster_file",
                             os.path.join(MAIN_PROJECT_DIR, "HWSD", "HWSD2.bil"))
HWSD_DB_FILE       = cfg_get("hwsd_db_file",
                             os.path.join(MAIN_PROJECT_DIR, "HWSD", "HWSD2.sqlite"))

# Constructed names
SOIL_BASENAME    = f"{GRID_BASE_NAME}_{SOIL_SOURCE}"
WEATHER_DIR_NAME = f"{GRID_BASE_NAME}_{WEATHER_SOURCE}"
SCENARIO_ID      = f"{GRID_BASE_NAME}_{WEATHER_SOURCE}_{SOIL_SOURCE}"

if RUN_NAME_OVERRIDE:
    DSSAT_RUN_NAME = RUN_NAME_OVERRIDE
elif RUN_TAG:
    if RUN_NAME_STYLE == "scenario":
        DSSAT_RUN_NAME = f"{SCENARIO_ID}_{RUN_TAG}"
    else:
        DSSAT_RUN_NAME = f"{GRID_BASE_NAME}_{RUN_TAG}"
else:
    DSSAT_RUN_NAME = SCENARIO_ID

DSSAT_RUN_NAME = re.sub(r"[^A-Za-z0-9_\-]", "_", DSSAT_RUN_NAME)

# --- 0.6 Dynamic paths ------------------------------------------------------
GRIDMET_CACHE_DIR     = os.path.join(MAIN_PROJECT_DIR, "gridmet_netcdf_cache")
CHIRPS_CACHE_DIR      = os.path.join(MAIN_PROJECT_DIR, "chirps_netcdf_cache")
AGERA5_CACHE_DIR      = os.path.join(MAIN_PROJECT_DIR, "agera5_netcdf_cache")
GRIDPOINTS_OUTPUT_DIR = os.path.join(MAIN_PROJECT_DIR, "gridpoints")
POINT_SHAPEFILE_NAME  = f"{GRID_BASE_NAME}.shp"
POINT_SHAPEFILE_PATH  = os.path.join(GRIDPOINTS_OUTPUT_DIR, POINT_SHAPEFILE_NAME)
DSSAT_RUN_DIR         = os.path.join(MAIN_PROJECT_DIR, "dssat_runs", DSSAT_RUN_NAME)
FINAL_OUTPUT_DIR      = os.path.join(MAIN_PROJECT_DIR, "results")
FINAL_RESULTS_PATH    = os.path.join(FINAL_OUTPUT_DIR, f"{DSSAT_RUN_NAME}_results.csv")
FINAL_PLOT_PATH       = os.path.join(FINAL_OUTPUT_DIR, f"{DSSAT_RUN_NAME}_yield_map.png")

POINT_ID_COLUMN = "ID"
LAT_COLUMN      = "LAT"
LONG_COLUMN     = "LONG"

# --- 0.7 Weather extension --------------------------------------------------
EXTEND_WEATHER_DATA    = False
WEATHER_REFERENCE_YEAR = 2025

TEMPLATE_SOIL_ID_PLACEHOLDER = "SOIL_ID"

# --- 0.8 DSSAT settings -----------------------------------------------------
DSSAT_EXE_PATH     = os.environ.get("DSSAT_EXE",
                                    os.path.join(DSSAT_BASE, DSSAT_EXE_NAME))
TEMPLATE_DIR       = os.path.join(MAIN_PROJECT_DIR, "dssat_templates")
TEMPLATE_FILE_NAME = cfg_get("template_file_name", "UFGA8201.MZX")   # DEMO PLACEHOLDER — replace with your own
TEMPLATE_FILE_PATH = os.path.join(TEMPLATE_DIR, TEMPLATE_FILE_NAME)

# --- 0.9 Run mode -----------------------------------------------------------
RUN_MODE        = cfg_get("run_mode", "experiment")   # "experiment" | "sequence"
TREATMENT_START = cfg_get("treatment_start", 1)
TREATMENT_END   = cfg_get("treatment_end", 4)
SEQUENCE_START  = cfg_get("sequence_start", 1)
SEQUENCE_END    = cfg_get("sequence_end", 1)

# --- 0.10 HPC & switches ----------------------------------------------------
ZIP_FOR_HPC         = False
RUN_STEP_1_SOILS    = True
RUN_STEP_2_WEATHER  = True
RUN_DSSAT_EXECUTION = True
CLEANUP_RUN_FOLDERS = False
RESUME_DSSAT_RUNS   = False
# When a DSSATPRO.V48 is available next to the executable, DSSAT resolves
# genotype/species/SDA/CO2 support files from the install directory, so they do
# NOT need copying into every run folder — a big metadata saving at scale
# (thousands of points × ~27 files). Set bundle_genotype_files: true to force
# copying them anyway, for self-contained folders you zip and ship to a machine
# whose DSSATPRO does not point at a matching install.
BUNDLE_GENOTYPE_FILES = bool(cfg_get("bundle_genotype_files", False))

# --- 0.11 Parallelism -------------------------------------------------------
N_CORES_TO_USE = max(1, multiprocessing.cpu_count() - 4)
SOIL_CORES     = N_CORES_TO_USE
WEATHER_CORES  = N_CORES_TO_USE
DSSAT_CORES    = N_CORES_TO_USE

# =============================================================================
# SECTION 1 — LOAD HELPER MODULES
# =============================================================================
sys.path.insert(0, PY_SCRIPTS_DIR)

# Try to resolve local dssatutils and dssatengine package paths if not already installed in the environment
try:
    import dssatutils
except ImportError:
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dssatutils_path = os.path.join(parent_dir, "dssatutils", "python")
    if os.path.isdir(dssatutils_path):
        sys.path.insert(0, dssatutils_path)

try:
    import dssatengine
except ImportError:
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dssatengine_path = os.path.join(parent_dir, "dssatengine", "python")
    if os.path.isdir(dssatengine_path):
        sys.path.insert(0, dssatengine_path)

from dssatutils.weather_daymet        import process_weather_daymet
from dssatutils.weather_nasapower     import process_weather_nasapower
from dssatutils.weather_gridmet       import process_weather_gridmet
from dssatutils.weather_openmeteo     import process_weather_openmeteo
from dssatutils.weather_nasapower_chirps import process_weather_nasapower_chirps
from dssatutils.weather_agera5        import process_weather_agera5
from dssatutils.soil_ssurgo           import process_soils_ssurgo
from dssatutils.soil_soilgrids        import process_soils_soilgrids
from dssatutils.soil_soilgrids_online import process_soils_soilgrids_online
from dssatutils.soil_hwsd             import process_soils_hwsd

print(f"Sourcing helper modules from: {PY_SCRIPTS_DIR}")

# =============================================================================
# SECTION 1b — PRE-FLIGHT CHECKS
# =============================================================================
print("Running Pre-flight Checks...")

if not USE_EXISTING_POINT_SHAPEFILE:
    shp = os.path.join(SHAPEFILE_DIR, BOUNDARY_SHAPEFILE_NAME)
    if not os.path.exists(shp):
        sys.exit(f"CRITICAL: Boundary shapefile not found: {shp}")

if not os.path.exists(TEMPLATE_FILE_PATH):
    sys.exit(f"CRITICAL: DSSAT Template file not found: {TEMPLATE_FILE_PATH}")

if SOIL_SOURCE == "HWSD":
    for _f, _what in [(HWSD_RASTER_FILE, "raster"), (HWSD_DB_FILE, "database")]:
        if not os.path.exists(_f):
            sys.exit(
                f"CRITICAL: HWSD {_what} file not found: {_f}\n"
                f"Download HWSD v2.0 from FAO and set hwsd_raster_file / "
                f"hwsd_db_file in config.yml."
            )
elif (SOIL_SOURCE not in ("SSURGO", "SOILGRIDS_ONLINE")
        and not os.path.exists(EXTERNAL_SOIL_FILE)):
    sys.exit(
        f"CRITICAL: External soil file needed for {SOIL_SOURCE} "
        f"but not found at: {EXTERNAL_SOIL_FILE}"
    )

print("All checks passed. Starting pipeline...")

# =============================================================================
# SECTION 2 — HELPER FUNCTIONS
# =============================================================================

from dssatengine import (
    create_grid_points,
    load_existing_points,
    extend_weather_repeat_single_ignore_partial,
    _write_dssbatch,
    _write_dssbatch_sequence,
    _run_dssat,
    _read_csv_safe,
    _merge_supplemental,
    _build_result_rows,
    _run_simulation,
    _run_one_point
)


if __name__ == '__main__':
    # =============================================================================
    # STEP 0 — CREATE / LOAD GRID POINTS
    # =============================================================================
    print("=" * 60)
    print("STEP 0: PREPARING GRIDFILE / POINTS")
    print("=" * 60)

    os.makedirs(GRIDPOINTS_OUTPUT_DIR, exist_ok=True)

    if USE_EXISTING_POINT_SHAPEFILE:
        print(f"Using existing point shapefile: {EXISTING_POINT_SHAPEFILE_PATH}")
        gridfile = load_existing_points(EXISTING_POINT_SHAPEFILE_PATH, POINT_SHAPEFILE_PATH)
    else:
        print(f"Generating grid at {GRID_SPACING_METERS} m from: {BOUNDARY_SHAPEFILE_NAME}")
        BOUNDARY_SHAPEFILE_PATH = os.path.join(SHAPEFILE_DIR, BOUNDARY_SHAPEFILE_NAME)
        if not os.path.exists(BOUNDARY_SHAPEFILE_PATH):
            sys.exit(f"Shapefile not found: {BOUNDARY_SHAPEFILE_PATH}")

        boundary_sf = gpd.read_file(BOUNDARY_SHAPEFILE_PATH)
        if ENABLE_BOUNDARY_FILTER:
            fstr = ", ".join(map(str, BOUNDARY_FILTER_VALUE))
            print(f"Filtering boundary where {BOUNDARY_FILTER_COLUMN} in [{fstr}]")
            boundary_sf = boundary_sf[
                boundary_sf[BOUNDARY_FILTER_COLUMN].isin(BOUNDARY_FILTER_VALUE)
            ]

        if boundary_sf.empty:
            sys.exit("Filter resulted in 0 features.")

        gridfile = create_grid_points(boundary_sf, GRID_SPACING_METERS, POINT_SHAPEFILE_PATH)

    print(f"Grid points ready: {len(gridfile)} points")

    # =============================================================================
    # STEP 1 — SOIL DATA
    # =============================================================================
    print("=" * 60)
    print("STEP 1: PROCESSING SOIL DATA")
    print("=" * 60)

    if RUN_STEP_1_SOILS:
        os.makedirs(SOIL_ROOT_DIR, exist_ok=True)
        soilfile_prefix    = os.path.join(SOIL_ROOT_DIR, SOIL_BASENAME)
        soilfile_DSSAT     = soilfile_prefix + ".SOL"
        soilfile_CSV       = soilfile_prefix + ".CSV"
        individual_sol_dir = soilfile_prefix + "_individual_SOL"
        os.makedirs(individual_sol_dir, exist_ok=True)

        if SOIL_SOURCE == "SSURGO":
            process_soils_ssurgo(
                grid_points=gridfile,
                output_dir_csv=soilfile_CSV,
                output_dir_individual=individual_sol_dir,
                n_cores=SOIL_CORES,
                id_col=POINT_ID_COLUMN,
                lat_col=LAT_COLUMN,
                long_col=LONG_COLUMN,
            )
            sol_files = sorted([f for f in os.listdir(individual_sol_dir) if f.endswith(".SOL")])
            with open(soilfile_DSSAT, "w") as out_fh:
                out_fh.write("*SOILS: Combined\n")
                for sf in sol_files:
                    with open(os.path.join(individual_sol_dir, sf)) as in_fh:
                        data       = in_fh.read()
                        first_star = data.find("\n*")
                        if first_star == -1:
                            first_star = 0 if data.startswith("*") else len(data)
                        else:
                            first_star += 1
                        out_fh.write(data[first_star:])

        elif SOIL_SOURCE == "SOILGRIDS_10K":
            if not os.path.exists(EXTERNAL_SOIL_FILE):
                sys.exit(f"External soil file not found: {EXTERNAL_SOIL_FILE}")
            process_soils_soilgrids(
                grid_points=gridfile,
                source_sol_file=EXTERNAL_SOIL_FILE,
                output_csv_path=soilfile_CSV,
                output_sol_dir=individual_sol_dir,
                id_col=POINT_ID_COLUMN,
            )

        elif SOIL_SOURCE == "SOILGRIDS_ONLINE":
            # Choose REST (one HTTP request per point, rate-limited) or VRT (GDAL
            # virtual rasters streamed from ISRIC; reads each global raster once and
            # samples all points, so it scales better and avoids rate limits, but
            # needs rasterio). Set via config key `soilgrids_mode`.
            import dssatutils.soil_soilgrids_online as _sg_mod
            _sg_mod.USE_REST_API = (SOILGRIDS_MODE != "VRT")
            print(f"SoilGrids online mode: {'VRT' if SOILGRIDS_MODE == 'VRT' else 'REST API'}")
            process_soils_soilgrids_online(
                gridfile=gridfile,
                soilfile_csv_path=soilfile_CSV,
                output_sol_dir=individual_sol_dir,
                id_col=POINT_ID_COLUMN,
            )

        elif SOIL_SOURCE == "HWSD":
            process_soils_hwsd(
                grid_points=gridfile,
                hwsd_raster_file=HWSD_RASTER_FILE,
                hwsd_db_file=HWSD_DB_FILE,
                output_csv_path=soilfile_CSV,
                output_sol_dir=individual_sol_dir,
                id_col=POINT_ID_COLUMN,
                lat_col=LAT_COLUMN,
                long_col=LONG_COLUMN,
            )

        else:
            sys.exit(f"Unknown SOIL_SOURCE: {SOIL_SOURCE}")

    # =============================================================================
    # STEP 2 — WEATHER DATA
    # =============================================================================
    print("=" * 60)
    print("STEP 2: DOWNLOADING WEATHER DATA")
    print("=" * 60)

    if RUN_STEP_2_WEATHER:
        os.makedirs(WEATHER_ROOT_DIR, exist_ok=True)
        weather_dir = os.path.join(WEATHER_ROOT_DIR, WEATHER_DIR_NAME)
        os.makedirs(weather_dir, exist_ok=True)

        existing_wth      = {os.path.splitext(f)[0]
                             for f in os.listdir(weather_dir) if f.endswith(".WTH")}
        missing_mask      = ~gridfile[POINT_ID_COLUMN].astype(str).isin(existing_wth)
        points_to_process = gridfile[missing_mask].reset_index(drop=True)

        if points_to_process.empty:
            print("All weather files already exist. Skipping processing.")
        else:
            print(f"Resuming: Processing {len(points_to_process)} remaining points...")
            log_file    = os.path.join(weather_dir, "download_errors.log")
            common_args = dict(
                shapefile=points_to_process,
                start_year=WEATHER_START_YEAR,
                end_year=WEATHER_END_YEAR,
                output_dir=weather_dir,
                id_col=POINT_ID_COLUMN,
                lat_col=LAT_COLUMN,
                lon_col=LONG_COLUMN,
                n_cores=WEATHER_CORES,
                log_file=log_file,
            )
            if WEATHER_SOURCE == "DAYMET":
                process_weather_daymet(**common_args)
            elif WEATHER_SOURCE == "NASA_POWER":
                process_weather_nasapower(**common_args)
            elif WEATHER_SOURCE == "GRIDMET":
                process_weather_gridmet(**common_args, gridmet_cache_dir=GRIDMET_CACHE_DIR)
            elif WEATHER_SOURCE == "OPEN_METEO":
                process_weather_openmeteo(**common_args)
            elif WEATHER_SOURCE == "NASA_POWER_CHIRPS":
                import dssatutils.weather_nasapower_chirps as _wc
                _wc.CHIRPS_RESOLUTION = CHIRPS_RESOLUTION
                process_weather_nasapower_chirps(
                    **common_args, chirps_cache_dir=CHIRPS_CACHE_DIR)
            elif WEATHER_SOURCE == "AGERA5":
                process_weather_agera5(**common_args, agera5_cache_dir=AGERA5_CACHE_DIR)
            else:
                sys.exit(f"Unknown WEATHER_SOURCE: {WEATHER_SOURCE}")

        if EXTEND_WEATHER_DATA:
            all_wth = [os.path.join(weather_dir, f)
                       for f in os.listdir(weather_dir) if f.endswith(".WTH")]
            print("Checking which files need extension...")

            def _needs_extension(f: str) -> bool:
                try:
                    with open(f) as fh:
                        last = fh.readlines()[-1]
                    m = re.match(r"^\s*(\d+)", last)
                    if not m:
                        return True
                    dc   = int(m.group(1))
                    year = (dc // 1000 if len(str(dc)) > 5
                            else (2000 + dc // 1000 if dc // 1000 < 80 else 1900 + dc // 1000))
                    return year <= WEATHER_END_YEAR
                except Exception:
                    return True

            to_extend = [f for f in all_wth if _needs_extension(f)]
            if to_extend:
                print(f"Extending {len(to_extend)} file(s)...")
                for wth_f in to_extend:
                    extend_weather_repeat_single_ignore_partial(
                        wth_f,
                        ref_start_year=WEATHER_START_YEAR,
                        ref_end_year=WEATHER_REFERENCE_YEAR,
                        target_end_year=WEATHER_END_YEAR,
                    )
            else:
                print("All files already extended.")

    # =============================================================================
    # STEP 3 — BUILD DSSAT FOLDERS AND RUN SIMULATIONS
    # =============================================================================
    print("=" * 60)
    print("STEP 3: RUNNING DSSAT SIMULATIONS")
    print("=" * 60)

    os.makedirs(DSSAT_RUN_DIR, exist_ok=True)

    points = gridfile.copy()
    points[POINT_ID_COLUMN] = points[POINT_ID_COLUMN].astype(str)

    soil_mapping_file      = os.path.join(SOIL_ROOT_DIR, f"{SOIL_BASENAME}.CSV")
    individual_soil_folder = os.path.join(SOIL_ROOT_DIR, f"{SOIL_BASENAME}_individual_SOL")
    weather_repo           = os.path.join(WEATHER_ROOT_DIR, WEATHER_DIR_NAME)

    if os.path.exists(soil_mapping_file):
        print("Loading soil mapping file...")
        soil_map    = pd.read_csv(soil_mapping_file, dtype=str)
        if POINT_ID_COLUMN in soil_map.columns:
            soil_map = soil_map.drop_duplicates(subset=[POINT_ID_COLUMN])
        common_cols = [c for c in set(points.columns) & set(soil_map.columns)
                       if c != POINT_ID_COLUMN]
        if common_cols:
            soil_map = soil_map.drop(columns=common_cols)
        if POINT_ID_COLUMN in soil_map.columns:
            points = points.merge(soil_map, on=POINT_ID_COLUMN, how="left")
    else:
        print("WARNING: Soil mapping CSV not found. Attempting legacy logic.")
        if SOIL_SOURCE == "SSURGO":
            points["SOIL_ID"] = points[POINT_ID_COLUMN]

    if "SOIL_ID" not in points.columns:
        points["SOIL_ID"] = points[POINT_ID_COLUMN] if SOIL_SOURCE == "SSURGO" else None

    # --- Precompute run-folder support files ONCE (not per point) ---------------
    # DSSATPRO.V48 lets DSSAT find genotype/SDA/CO2 files in the install dir, so we
    # only copy those ~27 files per point when DSSATPRO is unavailable or the user
    # explicitly wants self-contained folders (BUNDLE_GENOTYPE_FILES).
    _DSSATPRO_SRC      = os.path.join(os.path.dirname(DSSAT_EXE_PATH), "DSSATPRO.V48")
    _DSSATPRO_OK       = os.path.exists(_DSSATPRO_SRC)
    _SUPPORT_EXTS      = {".CUL", ".ECO", ".SPE", ".SDA", ".WDA", ".CDE"}
    # Bundle support files when shipping folders elsewhere (ZIP_FOR_HPC) — the local
    # DSSATPRO paths won't be valid on the target — or when explicitly requested, or
    # as a fallback when no DSSATPRO is available next to the executable.
    _COPY_SUPPORT      = BUNDLE_GENOTYPE_FILES or ZIP_FOR_HPC or not _DSSATPRO_OK
    _SUPPORT_FILES     = (
        [f for f in os.listdir(TEMPLATE_DIR)
         if os.path.splitext(f)[1].upper() in _SUPPORT_EXTS]
        if _COPY_SUPPORT else []
    )
    if _DSSATPRO_OK and not BUNDLE_GENOTYPE_FILES:
        print("Genotype files resolved via DSSATPRO.V48 (not copied per point — "
              "faster, fewer files).")
    elif not _DSSATPRO_OK:
        print(f"DSSATPRO.V48 not found next to executable; bundling "
              f"{len(_SUPPORT_FILES)} support files into each run folder.")


    def _create_folders_and_files(i: int) -> None:
        """Build one DSSAT run folder for row i of the points GeoDataFrame."""
        row              = points.iloc[i]
        ID               = str(row[POINT_ID_COLUMN])
        assigned_soil_id = str(row["SOIL_ID"]) if pd.notna(row.get("SOIL_ID")) else None

        point_dir = os.path.join(DSSAT_RUN_DIR, ID)
        os.makedirs(point_dir, exist_ok=True)

        if assigned_soil_id is None:
            with open(os.path.join(point_dir, "SOIL.SOL"), "w") as fh:
                fh.write("*SOIL ERROR\nNo Soil ID assigned\n")
            return

        if SOIL_SOURCE == "SSURGO":
            src_fname          = f"{ID}.SOL"
            hmx_replacement_id = ID
        else:
            src_fname          = f"{assigned_soil_id}.SOL"
            hmx_replacement_id = assigned_soil_id

        src_sol = os.path.join(individual_soil_folder, src_fname)
        dst_sol = os.path.join(point_dir, "SOIL.SOL")
        if os.path.exists(src_sol):
            shutil.copy2(src_sol, dst_sol)
        else:
            with open(dst_sol, "w") as fh:
                fh.write(f"*SOIL ERROR\nSource missing: {src_fname}\n")

        try:
            with open(TEMPLATE_FILE_PATH) as fh:
                content = fh.read()
            if TEMPLATE_SOIL_ID_PLACEHOLDER in content:
                content = content.replace(TEMPLATE_SOIL_ID_PLACEHOLDER, hmx_replacement_id)
            elif "ID_SOIL" in content:
                content = content.replace("ID_SOIL", hmx_replacement_id)
            # Patch WSTA in *FIELDS section: DSSAT opens <WSTA>.WTH from the run
            # folder. The template uses "00000000"; replace with the 8-char point
            # ID so it resolves to the per-point weather file (e.g. 00000001.WTH).
            # The *FIELDS header is "@L ID_FIELD WSTA...."; WSTA is the token right
            # after the ID_FIELD token, so we replace only the WSTA occurrence by
            # targeting the specific pattern in the fields data line.
            content = re.sub(
                r'(?m)^( \d+\s+\S{8}\s+)00000000(\s)',
                lambda m: m.group(1) + ID[:8].ljust(8) + m.group(2),
                content,
            )
            # NOTE: Do NOT patch the LATITUDE/LONGITUDE/ELEV text placeholders in the
            # *FIELDS coordinate line. When DSSAT encounters non-numeric text there it
            # falls back to reading coordinates from the .WTH file header (@ INSI LAT
            # LONG …), which already has the correct per-point coordinates written by
            # the weather-download step. Patching with numeric values causes DSSAT to
            # use those values instead and then the values end up wrong in summary.csv.
            with open(os.path.join(point_dir, TEMPLATE_FILE_NAME), "w") as fh:
                fh.write(content)
        except Exception as exc:
            print(f"Warning: could not write experiment file for {ID}: {exc}")

        # Copy DSSATPRO.V48 into the run folder so DSSAT resolves the genotype /
        # support directory paths (which point to the DSSAT48 installation root).
        if _DSSATPRO_OK:
            dssat_pro_dst = os.path.join(point_dir, "DSSATPRO.V48")
            if not os.path.exists(dssat_pro_dst):
                try:
                    os.link(_DSSATPRO_SRC, dssat_pro_dst)
                except OSError:
                    shutil.copy2(_DSSATPRO_SRC, dssat_pro_dst)

        wth_src = os.path.join(weather_repo, f"{ID}.WTH")
        wth_dst = os.path.join(point_dir, f"{ID}.WTH")
        if os.path.exists(wth_src) and not os.path.exists(wth_dst):
            try:
                os.symlink(wth_src, wth_dst)
            except (OSError, NotImplementedError):
                shutil.copy2(wth_src, wth_dst)

        # Genotype / support files (cultivar, ecotype, species, SDA, CO2, …). Only
        # bundled when DSSATPRO can't resolve them or the user forces it; otherwise
        # skipped entirely (precomputed _SUPPORT_FILES is empty). Hard-link to save
        # disk + metadata; fall back to copy across filesystems.
        for fname in _SUPPORT_FILES:
            src = os.path.join(TEMPLATE_DIR, fname)
            dst = os.path.join(point_dir, fname)
            if not os.path.exists(dst):
                try:
                    os.link(src, dst)
                except OSError:
                    shutil.copy2(src, dst)


    if not RESUME_DSSAT_RUNS:
        all_ids       = points[POINT_ID_COLUMN].tolist()
        existing_dirs = [
            d for d in os.listdir(DSSAT_RUN_DIR)
            if os.path.isdir(os.path.join(DSSAT_RUN_DIR, d)) and d in all_ids
        ]
        if existing_dirs:
            print(f"Deleting {len(existing_dirs)} old simulation folders...")
            for d in existing_dirs:
                shutil.rmtree(os.path.join(DSSAT_RUN_DIR, d), ignore_errors=True)

    print("Building DSSAT run folders...")
    for i in range(len(points)):
        _create_folders_and_files(i)

    # =============================================================================
    # STEP 3 (continued) — EXECUTE DSSAT
    # =============================================================================

    all_ids = points[POINT_ID_COLUMN].tolist()

    if RUN_DSSAT_EXECUTION:
        ids_to_run = all_ids
        if RESUME_DSSAT_RUNS:
            ids_to_run = [
                ID for ID in all_ids
                if not os.path.exists(
                    os.path.join(DSSAT_RUN_DIR, ID, f"results_{ID}.csv")
                )
            ]

        print(f"Starting DSSAT execution for {len(ids_to_run)} point(s)...")

        from concurrent.futures import ProcessPoolExecutor, as_completed

        # On macOS/Windows the default start method is 'spawn', which re-imports
        # this module in each worker and triggers the pool submission code again.
        # Using 'fork' avoids that for a pure command-line pipeline.
        import multiprocessing as _mp
        import platform
        _mp_ctx = None
        if hasattr(_mp, "get_context"):
            if platform.system() == "Windows":
                _mp_ctx = _mp.get_context("spawn")
            else:
                try:
                    _mp_ctx = _mp.get_context("fork")
                except ValueError:
                    _mp_ctx = _mp.get_context("spawn")

        # Package the task arguments for the top-level helper
        tasks = []
        for ID in ids_to_run:
            row = points[points[POINT_ID_COLUMN] == ID].iloc[0]
            row_dict = row.to_dict()
            tasks.append({
                "ID": ID,
                "row_dict": row_dict,
                "dssat_run_dir": DSSAT_RUN_DIR,
                "crop_extension": CROP_EXTENSION,
                "template_file_name": TEMPLATE_FILE_NAME,
                "template_file_path": TEMPLATE_FILE_PATH,
                "run_mode": RUN_MODE,
                "treatment_start": TREATMENT_START,
                "treatment_end": TREATMENT_END,
                "sequence_start": SEQUENCE_START,
                "sequence_end": SEQUENCE_END,
                "weather_start_year": WEATHER_START_YEAR,
                "weather_end_year": WEATHER_END_YEAR,
                "dssat_exe_path": DSSAT_EXE_PATH,
            })

        all_results = []
        with ProcessPoolExecutor(max_workers=DSSAT_CORES,
                                 mp_context=_mp_ctx) as pool:
            futures = {pool.submit(_run_one_point, t): t["ID"] for t in tasks}
            for fut in as_completed(futures):
                ID = futures[fut]
                try:
                    res = fut.result()
                    if res is not None and not res.empty:
                        all_results.append(res)
                except Exception as exc:
                    print(f"ERROR on ID {ID}: {exc}")

        os.makedirs(FINAL_OUTPUT_DIR, exist_ok=True)
        if all_results:
            final_data = pd.concat(all_results, ignore_index=True)
            final_data.to_csv(FINAL_RESULTS_PATH, index=False, na_rep="")
            print(f"Results combined -> {FINAL_RESULTS_PATH}")
        else:
            print("WARNING: No results produced.")

        if CLEANUP_RUN_FOLDERS:
            print("Cleaning up run folders...")
            for ID in all_ids:
                d = os.path.join(DSSAT_RUN_DIR, ID)
                if os.path.isdir(d):
                    shutil.rmtree(d, ignore_errors=True)

    else:
        # ==========================================================================
        # HPC PREP MODE — folders built, no execution
        # ==========================================================================
        print("=" * 60)
        print("  HPC PREP MODE COMPLETE")
        print("=" * 60)

        metadata_path = os.path.join(DSSAT_RUN_DIR, "README_CONFIG.txt")
        with open(metadata_path, "w") as fh:
            fh.write("=" * 60 + "\n")
            fh.write(f"DSSAT RUN CONFIGURATION: {datetime.now()}\n")
            fh.write("=" * 60 + "\n")
            fh.write(f"Project Name:    {PROJECT_NAME}\n")
            fh.write(f"Resolution:      {RESOLUTION_TAG} ({GRID_SPACING_METERS} meters)\n")
            fh.write(f"Grid/Shapefile:  {POINT_SHAPEFILE_NAME}\n")
            fh.write(f"Scenario ID:     {SCENARIO_ID}\n")
            fh.write(f"Run Folder Name: {DSSAT_RUN_NAME}\n")
            fh.write(f"Weather Folder:  {WEATHER_DIR_NAME}\n")
            fh.write(f"Soil Folder:     {SOIL_BASENAME}\n")
            fh.write("\n--- DATA SOURCES ---\n")
            fh.write(f"Weather Source:  {WEATHER_SOURCE}\n")
            fh.write(f"Years:           {WEATHER_START_YEAR} - {WEATHER_END_YEAR}\n")
            fh.write(f"Weather Extended?: {EXTEND_WEATHER_DATA} "
                     f"(Ref Year: {WEATHER_REFERENCE_YEAR})\n")
            fh.write(f"Soil Source:     {SOIL_SOURCE}\n")
            fh.write("\n--- SIMULATION SETTINGS ---\n")
            fh.write(f"Mode:            {RUN_MODE}\n")
            fh.write(f"Treatments:      {TREATMENT_START} - {TREATMENT_END}\n")
            fh.write(f"Sequences:       {SEQUENCE_START} - {SEQUENCE_END}\n")
            fh.write(f"Template File:   {TEMPLATE_FILE_NAME}\n")
            fh.write("\n--- PATHS ---\n")
            fh.write(f"Soil Map CSV:    {soil_mapping_file}\n")
        print(f"Metadata written -> {metadata_path}")

        if ZIP_FOR_HPC:
            print(f"Zipping run directory: {DSSAT_RUN_NAME}")
            zip_path = DSSAT_RUN_DIR + ".zip"
            with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
                for root, dirs, files in os.walk(DSSAT_RUN_DIR):
                    for file in files:
                        abs_path = os.path.join(root, file)
                        arc_path = os.path.relpath(abs_path, os.path.dirname(DSSAT_RUN_DIR))
                        zf.write(abs_path, arc_path)
            if os.path.exists(zip_path):
                print("Zip created successfully. Deleting uncompressed folder...")
                shutil.rmtree(DSSAT_RUN_DIR, ignore_errors=True)
                print(f"READY: {zip_path}")
            else:
                print("ERROR: Zip file not created. Original folder preserved.")
        else:
            print(f"Folders created in: {DSSAT_RUN_DIR}")
            print("DSSAT execution skipped. Transfer this directory to your HPC.")

    # =============================================================================
    # STEP 4 — VISUALIZE RESULTS
    # =============================================================================
    if RUN_DSSAT_EXECUTION and os.path.exists(FINAL_RESULTS_PATH):
        print("=" * 60)
        print("STEP 4: VISUALIZING RESULTS")
        print("=" * 60)

        try:
            import matplotlib.pyplot as plt
            import matplotlib.cm as cm

            sim_data = pd.read_csv(FINAL_RESULTS_PATH)

            avg_yield = (
                sim_data
                .dropna(subset=["longitude", "latitude", "final_grain_kg_ha"])
                .groupby(["point_id", "latitude", "longitude"], as_index=False)
                ["final_grain_kg_ha"].mean()
            )

            if avg_yield.empty:
                print("No yield data available for mapping.")
            else:
                fig, ax = plt.subplots(figsize=(10, 8))
                sc = ax.scatter(
                    avg_yield["longitude"],
                    avg_yield["latitude"],
                    c=avg_yield["final_grain_kg_ha"],
                    cmap="YlGn",
                    s=60,
                    edgecolors="none",
                    alpha=0.85,
                )
                plt.colorbar(sc, ax=ax, label="Mean Grain Yield (kg/ha)")
                ax.set_title(f"Simulated Yield — {DSSAT_RUN_NAME}", fontsize=13)
                ax.set_xlabel("Longitude")
                ax.set_ylabel("Latitude")
                ax.set_aspect("equal")
                plt.tight_layout()
                os.makedirs(FINAL_OUTPUT_DIR, exist_ok=True)
                fig.savefig(FINAL_PLOT_PATH, dpi=150)
                plt.close(fig)
                print(f"Yield map saved -> {FINAL_PLOT_PATH}")

        except Exception as exc:
            print(f"Step 4 visualization failed: {exc}")