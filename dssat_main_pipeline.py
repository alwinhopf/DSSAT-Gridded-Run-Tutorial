#!/usr/bin/env python3
"""
dssat_main_pipeline.py
======================
Python port of dssat_main_pipeline.R

End-to-end spatial / gridded DSSAT crop modeling pipeline:
Step 0 - Create or load a grid / point shapefile
Step 1 - Download and format soil data
Step 2 - Download and format weather data
Step 3 - Build DSSAT run folders and (optionally) run simulations
Step 4 - Combine results and visualize

HOW TO USE (beginners)
1. Edit SECTION 0 below (paths, crop, weather/soil sources, years).
2. Set DSSAT_EXE in your shell first (once):
   export DSSAT_EXE=/full/path/to/dscsm048   # Linux / macOS
   set DSSAT_EXE=C:\\DSSAT48\\DSCSM048.EXE    # Windows
3. Run:
   python dssat_main_pipeline.py

THREE SPATIAL DOMAIN MODES (same as R pipeline)
MODE A - Regular grid from a boundary polygon (default demo)
MODE B - Your own point / polygon shapefile
MODE C - Cropland-only points from CDL / NLCD (build with landcover scripts,
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
    from shapely.geometry import Point, box, mapping
    from shapely.ops import transform as shp_transform
except ImportError as exc:
    sys.exit(
        f"Required package missing: {exc}\n"
        "Install with: pip install geopandas pandas numpy matplotlib pyproj shapely"
    )

# =============================================================================
# SECTION 0 - MASTER CONFIGURATION
# =============================================================================

# --- 0.1 Path & platform detection -----------------------------------------

def _detect_project_dir() -> str:
    """Return the directory that contains this script."""
    return str(Path(__file__).resolve().parent)

MAIN_PROJECT_DIR = _detect_project_dir()
CODE_ROOT_DIR    = MAIN_PROJECT_DIR

print(f"Running engine code from: {CODE_ROOT_DIR}")

# --- 0.2 DSSAT executable ---------------------------------------------------
_os = platform.system()
if _os == "Windows":
    DSSAT_BASE     = r"C:\DSSAT48"
    DSSAT_EXE_NAME = "DSCSM048.EXE"
else:
    # macOS / Linux: compile from source - see README
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


def resolve_config_path(path, base=CODE_ROOT_DIR) -> str:
    path = str(path or "").strip()
    if not path:
        return ""
    p = Path(path).expanduser()
    if not p.is_absolute():
        p = Path(base) / p
    return str(p.resolve())


# CODE_ROOT_DIR holds code/templates/static resources. OUTPUT_ROOT_DIR holds
# generated model artifacts (gridpoints, weather, soils, run folders, results).
INPUT_ROOT_DIR  = resolve_config_path(cfg_get("input_root_dir", CODE_ROOT_DIR), CODE_ROOT_DIR)
OUTPUT_ROOT_DIR = resolve_config_path(cfg_get("output_root_dir", CODE_ROOT_DIR), CODE_ROOT_DIR)
PROJECT_ROOT     = OUTPUT_ROOT_DIR
R_SCRIPTS_DIR    = os.path.join(CODE_ROOT_DIR, "r_scripts")   # unused in Python port
PY_SCRIPTS_DIR   = os.path.join(CODE_ROOT_DIR, "py_scripts")
SHAPEFILE_DIR    = resolve_config_path(cfg_get("shapefile_dir", os.path.join(INPUT_ROOT_DIR, "shapefile")), CODE_ROOT_DIR)
GRIDPOINTS_DIR   = resolve_config_path(cfg_get("gridpoints_dir", os.path.join(OUTPUT_ROOT_DIR, "gridpoints")), CODE_ROOT_DIR)
WEATHER_ROOT_DIR = resolve_config_path(cfg_get("weather_root_dir", os.path.join(OUTPUT_ROOT_DIR, "weather")), CODE_ROOT_DIR)
SOIL_ROOT_DIR    = resolve_config_path(cfg_get("soil_root_dir", os.path.join(OUTPUT_ROOT_DIR, "soil")), CODE_ROOT_DIR)
RUNS_ROOT_DIR    = resolve_config_path(cfg_get("runs_root_dir", os.path.join(OUTPUT_ROOT_DIR, "dssat_runs")), CODE_ROOT_DIR)
RESULTS_ROOT_DIR = resolve_config_path(cfg_get("results_root_dir", os.path.join(OUTPUT_ROOT_DIR, "results")), CODE_ROOT_DIR)

print(f"Engine input/static root: {INPUT_ROOT_DIR}")
print(f"Engine generated-output root: {OUTPUT_ROOT_DIR}")

# --- 0.3 Project settings ---------------------------------------------------
PROJECT_NAME        = cfg_get("project_name", "dssat_spatial_demo")
GRID_SPACING_METERS = cfg_get("grid_spacing_meters", 50_000)   # 50 km test; 5-10 km production
CROP_EXTENSION      = cfg_get("crop_extension", "MZ")     # MZ=maize  WH=wheat  SB=soybean ...

# --- 0.3b Optional run-folder naming ----------------------------------------
RUN_TAG           = cfg_get("run_tag", "")        # e.g. "run1" / "calibA"
RUN_NAME_STYLE    = cfg_get("run_name_style", "grid")    # "grid" | "scenario"
RUN_NAME_OVERRIDE = cfg_get("run_name_override", "")        # If set, used verbatim as run-folder name

# --- 0.4 Spatial domain (choose MODE A / B / C) ----------------------------
USE_EXISTING_POINT_SHAPEFILE  = cfg_get("use_existing_point_shapefile", False)
EXISTING_POINT_SHAPEFILE_PATH = cfg_get("existing_point_shapefile_path",
                                        os.path.join(GRIDPOINTS_DIR, "my_points.shp"))

# MODE A boundary settings (ignored when USE_EXISTING_POINT_SHAPEFILE = True)
BOUNDARY_SHAPEFILE_NAME = cfg_get("boundary_shapefile_name", "tl_2024_us_state.shp")
ENABLE_BOUNDARY_FILTER  = cfg_get("enable_boundary_filter", True)
BOUNDARY_FILTER_COLUMN  = cfg_get("boundary_filter_column", "NAME")
STATE_NAME_FILTER       = list(cfg_get("state_name_filter", ["Iowa"]))

# Optional cropland mask. Shapefile fields stay <=10 chars:
# crop_frac [0-1], crop_pct [0-100], crop_ha, cell_ha.
USE_CROPLAND_MASK     = bool(cfg_get("use_cropland_mask", False))
CROPLAND_RASTER_FILE  = str(cfg_get("cropland_raster_file", ""))
CROPLAND_CLASSES      = [int(x) for x in cfg_get("cropland_classes", [82])]
CROPLAND_MIN_FRACTION = float(cfg_get("cropland_min_fraction", 0))
CROPLAND_STRICT       = bool(cfg_get("cropland_strict", False))
REUSE_CROPLAND_GRID   = bool(cfg_get("reuse_cropland_grid", True))

# --- 0.5 Auto-naming convention ---------------------------------------------
if GRID_SPACING_METERS < 1000:
    RESOLUTION_TAG = f"{GRID_SPACING_METERS}m"
else:
    RESOLUTION_TAG = f"{GRID_SPACING_METERS // 1000}km"

GRID_BASE_NAME        = f"{PROJECT_NAME}_{RESOLUTION_TAG}"
BOUNDARY_FILTER_VALUE = STATE_NAME_FILTER

# Weather settings
WEATHER_SOURCE     = cfg_get("weather_source", "GRIDMET")   # DAYMET | NASA_POWER | GRIDMET | OPEN_METEO | NASA_POWER_CHIRPS
WEATHER_START_YEAR = int(cfg_get("weather_start_year", 1982))
WEATHER_END_YEAR   = int(cfg_get("weather_end_year", 1983))
# NASA_POWER_CHIRPS only: CHIRPS rainfall resolution "p05" (~5.5 km, recommended)
# or "p25" (~28 km, lighter download).
CHIRPS_RESOLUTION  = str(cfg_get("chirps_resolution", "p05"))

# Soil settings
SOIL_SOURCE        = cfg_get("soil_source", "SSURGO")   # "SSURGO" | "SOILGRIDS_10K" | "SOILGRIDS_ONLINE" | "POLARIS"
EXTERNAL_SOIL_FILE = cfg_get("external_soil_file",
                             os.path.join(INPUT_ROOT_DIR, "SoilGrids", "US.SOL"))
# SOILGRIDS_ONLINE only: "REST" (JSON API, rate-limited, no extra deps) or
# "VRT" (GDAL virtual rasters via rasterio; batch-friendly, better coverage).
SOILGRIDS_MODE     = str(cfg_get("soilgrids_mode", "REST")).upper()
# POLARIS only (CONUS 30 m): which statistic layer to build the profile from.
# "p50" (median) is the deterministic drop-in; p5/p95 are reserved for a future
# soil-uncertainty ensemble. Optional polaris_cache_dir caches the GeoTIFF tiles.
POLARIS_STAT       = str(cfg_get("polaris_stat", "p50"))
POLARIS_CACHE_DIR  = cfg_get("polaris_cache_dir", "") or None
# HWSD only: paths to the FAO HWSD v2.0 raster (SMU IDs) + SQLite database,
# downloaded once from FAO (blank = script defaults under HWSD/).
HWSD_RASTER_FILE   = cfg_get("hwsd_raster_file",
                             os.path.join(INPUT_ROOT_DIR, "HWSD", "HWSD2.bil"))
HWSD_DB_FILE       = cfg_get("hwsd_db_file",
                             os.path.join(INPUT_ROOT_DIR, "HWSD", "HWSD2.sqlite"))
# E-OBS only: folder of pre-downloaded E-OBS NetCDFs (tx/tn/rr/qq...). Set
# eobs_use_cds: true to fetch an area subset via the Copernicus CDS instead
# (needs the same ~/.cdsapirc key as AgERA5).
EOBS_NC_DIR        = cfg_get("eobs_nc_dir",
                             os.path.join(INPUT_ROOT_DIR, "eobs_netcdf"))
EOBS_USE_CDS       = bool(cfg_get("eobs_use_cds", False))
# Xavier (Brazil) / CMFD (China) only: folders of pre-downloaded NetCDFs.
XAVIER_NC_DIR      = cfg_get("xavier_nc_dir", os.path.join(INPUT_ROOT_DIR, "xavier_netcdf"))
CMFD_NC_DIR        = cfg_get("cmfd_nc_dir", os.path.join(INPUT_ROOT_DIR, "cmfd_netcdf"))
# LUCAS (Europe) only: the downloaded ESDAC LUCAS topsoil table (CSV/XLSX).
LUCAS_CSV          = cfg_get("lucas_csv", os.path.join(INPUT_ROOT_DIR, "LUCAS", "lucas_topsoil.csv"))

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
GRIDMET_CACHE_DIR     = os.path.join(OUTPUT_ROOT_DIR, "gridmet_netcdf_cache")
CHIRPS_CACHE_DIR      = os.path.join(OUTPUT_ROOT_DIR, "chirps_netcdf_cache")
AGERA5_CACHE_DIR      = os.path.join(OUTPUT_ROOT_DIR, "agera5_netcdf_cache")
DWD_CACHE_DIR         = os.path.join(OUTPUT_ROOT_DIR, "dwd_station_cache")
EOBS_CACHE_DIR        = os.path.join(OUTPUT_ROOT_DIR, "eobs_cds_cache")
GRIDPOINTS_OUTPUT_DIR = GRIDPOINTS_DIR
ALL_LAND_POINT_SHAPEFILE_NAME = f"{GRID_BASE_NAME}.shp"
CROPLAND_GRID_TAG = ""
if USE_CROPLAND_MASK:
    class_tag = "-".join(map(str, CROPLAND_CLASSES))
    min_tag = f"{CROPLAND_MIN_FRACTION:g}".replace(".", "p")
    CROPLAND_GRID_TAG = re.sub(r"[^A-Za-z0-9_\-]", "_", f"_cropland_{class_tag}_min{min_tag}")
POINT_SHAPEFILE_NAME  = f"{GRID_BASE_NAME}{CROPLAND_GRID_TAG}.shp"
ALL_LAND_POINT_SHAPEFILE_PATH = os.path.join(GRIDPOINTS_OUTPUT_DIR, ALL_LAND_POINT_SHAPEFILE_NAME)
POINT_SHAPEFILE_PATH  = os.path.join(GRIDPOINTS_OUTPUT_DIR, POINT_SHAPEFILE_NAME)
DSSAT_RUN_DIR         = os.path.join(RUNS_ROOT_DIR, DSSAT_RUN_NAME)
FINAL_OUTPUT_DIR      = RESULTS_ROOT_DIR
FINAL_RESULTS_PATH    = os.path.join(FINAL_OUTPUT_DIR, f"{DSSAT_RUN_NAME}_results.csv")
FINAL_PLOT_PATH       = os.path.join(FINAL_OUTPUT_DIR, f"{DSSAT_RUN_NAME}_yield_map.png")

POINT_ID_COLUMN = "ID"
LAT_COLUMN      = "LAT"
LONG_COLUMN     = "LONG"

# --- 0.7 Weather extension --------------------------------------------------
EXTEND_WEATHER_DATA    = bool(cfg_get("extend_weather_data", False))
WEATHER_REFERENCE_YEAR = int(cfg_get("weather_reference_year", WEATHER_END_YEAR))

TEMPLATE_SOIL_ID_PLACEHOLDER = "SOIL_ID"

# --- 0.8 DSSAT settings -----------------------------------------------------
DSSAT_EXE_PATH     = os.environ.get("DSSAT_EXE",
                                    os.path.join(DSSAT_BASE, DSSAT_EXE_NAME))
# Single canonical template dir shared by ALL workflows. Genotype + FileX
# templates live ONLY in DSSAT_Gridded_Run_Tutorial/dssat_templates (this repo) -
# runs do NOT fall back to the DSSAT48 install for these, so copy any new
# .CUL/.ECO/.SPE/.SDA/.WDA/FileX into that folder. Override with the
# `template_dir` config key or the DSSAT_TEMPLATE_DIR env var.
_DEFAULT_TEMPLATE_DIR = os.path.join(INPUT_ROOT_DIR, "dssat_templates")
TEMPLATE_DIR       = os.environ.get("DSSAT_TEMPLATE_DIR",
                                    cfg_get("template_dir", _DEFAULT_TEMPLATE_DIR))
TEMPLATE_FILE_NAME = cfg_get("template_file_name", "UFGA8201.MZX")   # DEMO PLACEHOLDER - replace with your own
TEMPLATE_FILE_PATH = os.path.join(TEMPLATE_DIR, TEMPLATE_FILE_NAME)

# --- 0.9 Run mode -----------------------------------------------------------
RUN_MODE        = cfg_get("run_mode", "experiment")   # "experiment" | "sequence"
TREATMENT_START = int(cfg_get("treatment_start", 1))
TREATMENT_END   = int(cfg_get("treatment_end", 4))
SEQUENCE_START  = int(cfg_get("sequence_start", 1))
SEQUENCE_END    = int(cfg_get("sequence_end", 1))

# --- 0.10 HPC & switches ----------------------------------------------------
ZIP_FOR_HPC         = False
RUN_STEP_1_SOILS    = bool(cfg_get("run_step_1_soils", True))
RUN_STEP_2_WEATHER  = bool(cfg_get("run_step_2_weather", True))
RUN_DSSAT_EXECUTION = bool(cfg_get("run_dssat_execution", True))
CLEANUP_RUN_FOLDERS = bool(cfg_get("cleanup_run_folders", False))
RESUME_DSSAT_RUNS   = bool(cfg_get("resume_dssat_runs", False))
# When a DSSATPRO.V48 is available next to the executable, DSSAT resolves
# genotype/species/SDA/CO2 support files from the install directory, so they do
# NOT need copying into every run folder - a big metadata saving at scale
# (thousands of points x ~27 files). Set bundle_genotype_files: true to force
# copying them anyway, for self-contained folders you zip and ship to a machine
# whose DSSATPRO does not point at a matching install.
BUNDLE_GENOTYPE_FILES = bool(cfg_get("bundle_genotype_files", False))
# Link genotype/FileX support files into each run folder by symlink (cheap, no
# data copied) instead of copying. Falls back to a real copy automatically when
# the OS/filesystem does not support symlinks. Set use_symlinks: false to always
# copy (e.g. for folders you zip and ship elsewhere).
USE_SYMLINKS = bool(cfg_get("use_symlinks", True))


def _link_or_copy(src: str, dst: str, use_symlinks: bool = True) -> None:
    """Provision a file into a run folder: symlink (cheap) when use_symlinks, with
    an automatic copy fallback if the OS/filesystem rejects symlinks; otherwise
    copy. No-op if src is missing or dst already exists."""
    if not os.path.exists(src) or os.path.exists(dst):
        return
    if use_symlinks:
        try:
            os.symlink(src, dst)
            return
        except (OSError, NotImplementedError):
            pass
    try:
        shutil.copy2(src, dst)
    except OSError:
        pass


# --- 0.11 Parallelism -------------------------------------------------------
# Core counts are read from config.yml (soil_cores / weather_cores / dssat_cores).
# "auto" (the default) = all logical CPUs minus 1.
# Soil/weather steps are API/IO-bound; dssat step is CPU+disk-bound.
def _resolve_cores(key: str) -> int:
    v = cfg_get(key, "auto")
    if str(v).lower() == "auto" or not str(v).strip():
        return max(1, multiprocessing.cpu_count() - 1)
    return max(1, int(v))


def _resolve_optional_path(path: str) -> str:
    path = str(path or "").strip()
    if not path:
        return ""
    p = Path(path).expanduser()
    if not p.is_absolute():
        p = Path(INPUT_ROOT_DIR) / p
    return str(p)


def apply_cropland_mask(points_gdf, raster_file: str, grid_spacing_m: float,
                        classes, min_fraction: float = 0,
                        strict: bool = False):
    """Attach cropland fractions and keep only cropland-bearing grid cells."""
    def fallback(message: str):
        if strict:
            raise RuntimeError(message)
        print(f"WARNING: {message} Continuing with all grid cells.")
        return points_gdf

    raster_file = _resolve_optional_path(raster_file)
    if not raster_file or not os.path.exists(raster_file):
        return fallback(f"Cropland mask enabled, but cropland_raster_file is missing or not found: '{raster_file}'")

    try:
        import rasterio
        from rasterio.mask import mask as rio_mask
    except Exception as exc:  # noqa: BLE001
        return fallback(f"Cropland mask enabled, but rasterio is not available: {exc}")

    if not classes:
        return fallback("Cropland mask enabled, but cropland_classes is empty.")

    with rasterio.open(raster_file) as src:
        if src.crs is None:
            return fallback(f"Cropland raster has no CRS: {raster_file}")
        if src.crs.is_geographic:
            return fallback("Cropland raster is longitude/latitude. Reproject it to a meter-based CRS before grid-cell fraction extraction.")

        print(f"Applying cropland mask from {raster_file} (classes: {', '.join(map(str, classes))})")
        pts_r = points_gdf.to_crs(src.crs)
        half = grid_spacing_m / 2.0
        fractions = []
        for geom in pts_r.geometry:
            cell = box(geom.x - half, geom.y - half, geom.x + half, geom.y + half)
            try:
                arr, _ = rio_mask(src, [mapping(cell)], crop=True, filled=False)
                data = arr[0]
                valid = ~np.ma.getmaskarray(data)
                if src.nodata is not None:
                    valid &= np.asarray(data) != src.nodata
                if not np.any(valid):
                    fractions.append(np.nan)
                    continue
                fractions.append(float(np.isin(np.asarray(data)[valid], classes).mean()))
            except Exception:  # noqa: BLE001
                fractions.append(np.nan)

    out = points_gdf.copy()
    frac = np.clip(np.asarray(fractions, dtype=float), 0, 1)
    out["crop_frac"] = frac
    out["crop_pct"] = np.round(frac * 100, 4)
    out["cell_ha"] = (grid_spacing_m ** 2) / 10000.0
    out["crop_ha"] = np.round(out["cell_ha"] * out["crop_frac"], 4)

    if min_fraction <= 0:
        keep = out["crop_frac"].notna() & (out["crop_frac"] > 0)
    else:
        keep = out["crop_frac"].notna() & (out["crop_frac"] >= min_fraction)
    kept = int(keep.sum())
    pct = 100 * kept / len(out) if len(out) else 0
    print(f"Cropland mask retained {kept} of {len(out)} grid cells ({pct:.1f}%).")
    if kept == 0:
        raise RuntimeError("Cropland mask removed all grid cells. Lower cropland_min_fraction or check cropland raster/classes.")

    out = out.loc[keep].copy()
    try:
        out.to_file(POINT_SHAPEFILE_PATH)
    except Exception as exc:  # noqa: BLE001
        print(f"WARNING: Could not rewrite cropland-filtered grid shapefile: {exc}")
    return out


SOIL_CORES     = _resolve_cores("soil_cores")
WEATHER_CORES  = _resolve_cores("weather_cores")
DSSAT_CORES    = _resolve_cores("dssat_cores")
print(f"Parallelism: soil={SOIL_CORES}  weather={WEATHER_CORES}  dssat={DSSAT_CORES} core(s)")

# =============================================================================
# SECTION 1 - LOAD HELPER MODULES
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
from dssatutils.weather_dwd           import process_weather_dwd
from dssatutils.weather_eobs          import process_weather_eobs
from dssatutils.weather_xavier        import process_weather_xavier
from dssatutils.weather_cmfd          import process_weather_cmfd

# Soil sources that write one .SOL per grid point named by the point ID (so
# SOIL_ID == point ID and the per-point combine logic below applies). The other
# sources (SoilGrids / HWSD) map points to shared profile IDs instead.
_PER_POINT_SOIL = ("SSURGO", "GNATSGO", "ISDASOIL", "LUCAS")
# Soil sources that need no pre-downloaded external file (queried online).
_KEYLESS_ONLINE_SOIL = ("SSURGO", "GNATSGO", "ISDASOIL", "SOILGRIDS_ONLINE", "POLARIS")
from dssatutils.soil_ssurgo           import process_soils_ssurgo
from dssatutils.soil_gnatsgo          import process_soils_gnatsgo
from dssatutils.soil_isdasoil         import process_soils_isdasoil
from dssatutils.soil_lucas            import process_soils_lucas
from dssatutils.soil_soilgrids        import process_soils_soilgrids
from dssatutils.soil_soilgrids_online import process_soils_soilgrids_online
from dssatutils.soil_polaris          import process_soils_polaris
from dssatutils.soil_hwsd             import process_soils_hwsd

print(f"Sourcing helper modules from: {PY_SCRIPTS_DIR}")

# =============================================================================
# SECTION 1b - PRE-FLIGHT CHECKS
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
elif SOIL_SOURCE == "LUCAS":
    if not os.path.exists(LUCAS_CSV):
        sys.exit(
            f"CRITICAL: LUCAS topsoil table not found: {LUCAS_CSV}\n"
            f"Request it free from ESDAC (esdac.jrc.ec.europa.eu) and set "
            f"lucas_csv in config.yml."
        )
elif (SOIL_SOURCE not in _KEYLESS_ONLINE_SOIL
        and not os.path.exists(EXTERNAL_SOIL_FILE)):
    sys.exit(
        f"CRITICAL: External soil file needed for {SOIL_SOURCE} "
        f"but not found at: {EXTERNAL_SOIL_FILE}"
    )

print("All checks passed. Starting pipeline...")

# =============================================================================
# SECTION 2 - HELPER FUNCTIONS
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
    # STEP 0 - CREATE / LOAD GRID POINTS
    # =============================================================================
    print("=" * 60)
    print("STEP 0: PREPARING GRIDFILE / POINTS")
    print("=" * 60)

    os.makedirs(GRIDPOINTS_OUTPUT_DIR, exist_ok=True)

    if (not USE_EXISTING_POINT_SHAPEFILE and USE_CROPLAND_MASK
            and REUSE_CROPLAND_GRID and os.path.exists(POINT_SHAPEFILE_PATH)):
        print(f"Reusing existing cropland grid shapefile: {POINT_SHAPEFILE_PATH}")
        gridfile = gpd.read_file(POINT_SHAPEFILE_PATH)
    elif USE_EXISTING_POINT_SHAPEFILE:
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

        raw_grid_path = ALL_LAND_POINT_SHAPEFILE_PATH if USE_CROPLAND_MASK else POINT_SHAPEFILE_PATH
        gridfile = create_grid_points(boundary_sf, GRID_SPACING_METERS, raw_grid_path)

    if USE_CROPLAND_MASK and "crop_pct" not in gridfile.columns:
        gridfile = apply_cropland_mask(
            gridfile,
            CROPLAND_RASTER_FILE,
            GRID_SPACING_METERS,
            CROPLAND_CLASSES,
            CROPLAND_MIN_FRACTION,
            CROPLAND_STRICT,
        )

    print(f"Grid points ready: {len(gridfile)} points")

    # =============================================================================
    # STEP 1 - SOIL DATA
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

        if SOIL_SOURCE in _PER_POINT_SOIL:
            # These sources all write one Saxton-&-Rawls .SOL per point named by
            # the point ID, so the combine step below is shared. gNATSGO fills
            # SSURGO's US gaps; iSDAsoil = Africa 30 m; LUCAS = EU topsoil.
            _soil_kwargs = dict(
                grid_points=gridfile,
                output_dir_csv=soilfile_CSV,
                output_dir_individual=individual_sol_dir,
                n_cores=SOIL_CORES,
                id_col=POINT_ID_COLUMN,
                lat_col=LAT_COLUMN,
                long_col=LONG_COLUMN,
            )
            if SOIL_SOURCE == "SSURGO":
                process_soils_ssurgo(**_soil_kwargs)
            elif SOIL_SOURCE == "GNATSGO":
                process_soils_gnatsgo(**_soil_kwargs)
            elif SOIL_SOURCE == "ISDASOIL":
                process_soils_isdasoil(**_soil_kwargs)
            else:  # LUCAS
                process_soils_lucas(**_soil_kwargs, lucas_csv=LUCAS_CSV)
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
            try:
                process_soils_soilgrids_online(
                    gridfile=gridfile,
                    soilfile_csv_path=soilfile_CSV,
                    output_sol_dir=individual_sol_dir,
                    id_col=POINT_ID_COLUMN,
                )
            except Exception as exc:
                valid_ids = {
                    os.path.splitext(f)[0]
                    for f in os.listdir(individual_sol_dir)
                    if f.upper().endswith(".SOL")
                }
                valid_ids &= {str(x) for x in gridfile[POINT_ID_COLUMN].tolist()}
                if not valid_ids:
                    raise
                print(
                    "WARNING: SoilGrids online extraction failed, but "
                    f"{len(valid_ids)}/{len(gridfile)} existing valid profile(s) "
                    f"are available; missing IDs will be skipped. Cause: {exc}"
                )
                if not os.path.exists(soilfile_CSV):
                    pd.DataFrame({
                        POINT_ID_COLUMN: gridfile[POINT_ID_COLUMN].astype(str),
                        "SOIL_ID": gridfile[POINT_ID_COLUMN].astype(str),
                    }).to_csv(soilfile_CSV, index=False)

        elif SOIL_SOURCE == "POLARIS":
            # POLARIS = 30 m probabilistic disaggregation of SSURGO (CONUS).
            # Same output contract as SOILGRIDS_ONLINE: one .SOL per point named
            # by point ID + a CSV, streamed via GDAL /vsicurl. Water limits come
            # from POLARIS's van Genuchten curve (stat=p50 deterministic).
            try:
                process_soils_polaris(
                    gridfile=gridfile,
                    soilfile_csv_path=soilfile_CSV,
                    output_sol_dir=individual_sol_dir,
                    id_col=POINT_ID_COLUMN,
                    stat=POLARIS_STAT,
                    cache_dir=POLARIS_CACHE_DIR,
                )
            except Exception as exc:
                valid_ids = {
                    os.path.splitext(f)[0]
                    for f in os.listdir(individual_sol_dir)
                    if f.upper().endswith(".SOL")
                }
                valid_ids &= {str(x) for x in gridfile[POINT_ID_COLUMN].tolist()}
                if not valid_ids:
                    raise
                print(
                    "WARNING: POLARIS extraction failed, but "
                    f"{len(valid_ids)}/{len(gridfile)} existing valid profile(s) "
                    f"are available; missing IDs will be skipped. Cause: {exc}"
                )
                if not os.path.exists(soilfile_CSV):
                    pd.DataFrame({
                        POINT_ID_COLUMN: gridfile[POINT_ID_COLUMN].astype(str),
                        "SOIL_ID": gridfile[POINT_ID_COLUMN].astype(str),
                    }).to_csv(soilfile_CSV, index=False)

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
    # STEP 2 - WEATHER DATA
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
            elif WEATHER_SOURCE == "DWD":
                process_weather_dwd(**common_args, dwd_cache_dir=DWD_CACHE_DIR)
            elif WEATHER_SOURCE == "EOBS":
                process_weather_eobs(
                    **common_args,
                    eobs_nc_dir=EOBS_NC_DIR,
                    eobs_cache_dir=EOBS_CACHE_DIR,
                    eobs_use_cds=EOBS_USE_CDS,
                )
            elif WEATHER_SOURCE == "XAVIER":
                process_weather_xavier(**common_args, xavier_nc_dir=XAVIER_NC_DIR)
            elif WEATHER_SOURCE == "CMFD":
                process_weather_cmfd(**common_args, cmfd_nc_dir=CMFD_NC_DIR)
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
    # STEP 3 - BUILD DSSAT FOLDERS AND RUN SIMULATIONS
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
        if SOIL_SOURCE in _PER_POINT_SOIL:
            points["SOIL_ID"] = points[POINT_ID_COLUMN]

    if "SOIL_ID" not in points.columns:
        points["SOIL_ID"] = points[POINT_ID_COLUMN] if SOIL_SOURCE in _PER_POINT_SOIL else None

    # --- Precompute run-folder support files ONCE (not per point) ---------------
    # DSSATPRO.V48 lets DSSAT find genotype/SDA/CO2 files in the install dir, so we
    # only copy those ~27 files per point when DSSATPRO is unavailable or the user
    # explicitly wants self-contained folders (BUNDLE_GENOTYPE_FILES).
    _DSSATPRO_SRC      = os.path.join(os.path.dirname(DSSAT_EXE_PATH), "DSSATPRO.V48")
    _DSSATPRO_OK       = os.path.exists(_DSSATPRO_SRC)
    _SUPPORT_EXTS      = {".CUL", ".ECO", ".SPE", ".SDA", ".WDA", ".CDE"}
    # Genotype/support files ALWAYS come from the shared template dir, never the
    # DSSAT48 install: every run folder gets the .CUL/.ECO/.SPE/.SDA/.WDA/.CDE that
    # live in TEMPLATE_DIR, so the local copy overrides whatever DSSATPRO would
    # otherwise resolve from DSSAT48. Linked by symlink (USE_SYMLINKS, cheap) with
    # an automatic copy fallback. To add a crop, drop its files into TEMPLATE_DIR.
    _SUPPORT_FILES     = [f for f in os.listdir(TEMPLATE_DIR)
                          if os.path.splitext(f)[1].upper() in _SUPPORT_EXTS]
    # Symlink for local runs; force a real copy when building a portable bundle
    # (ZIP_FOR_HPC / bundle_genotype_files) - symlinks don't survive the move.
    _SUPPORT_SYMLINK   = USE_SYMLINKS and not ZIP_FOR_HPC and not BUNDLE_GENOTYPE_FILES
    print(f"Provisioning {len(_SUPPORT_FILES)} genotype/support file(s) per run "
          f"from {TEMPLATE_DIR} ({'symlink' if _SUPPORT_SYMLINK else 'copy'}).")


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

        if SOIL_SOURCE in _PER_POINT_SOIL:
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
                if f"   {TEMPLATE_SOIL_ID_PLACEHOLDER}" in content:
                    content = content.replace(f"   {TEMPLATE_SOIL_ID_PLACEHOLDER}", hmx_replacement_id.ljust(10))
                else:
                    content = content.replace(TEMPLATE_SOIL_ID_PLACEHOLDER, hmx_replacement_id)
            elif "ID_SOIL" in content:
                if "   ID_SOIL" in content:
                    content = content.replace("   ID_SOIL", hmx_replacement_id.ljust(10))
                else:
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
            # Substitute real per-point coordinates into the *FIELDS tier-2 line.
            # Preserve placeholder widths so DSSAT's fixed-column parser stays aligned.
            def _fit_field(value, width, digits):
                try:
                    value = float(value)
                except Exception:
                    return None
                if not math.isfinite(value):
                    return None
                out = f"{value:{width}.0f}" if digits == 0 else f"{value:{width}.{digits}f}"
                return out if len(out) == width else None

            def _wth_elev_default():
                wth = os.path.join(weather_repo, f"{ID}.WTH")
                try:
                    with open(wth, encoding="utf-8", errors="ignore") as wh:
                        lines = [next(wh, "") for _ in range(3)]
                    if len(lines) >= 3 and lines[1].lstrip().startswith("@"):
                        header = lines[1].split()
                        values = lines[2].split()
                        if "ELEV" in header:
                            idx = header.index("ELEV")
                            if idx < len(values):
                                elev = float(values[idx])
                                if math.isfinite(elev):
                                    return elev
                except Exception:
                    pass
                return -99.0

            coord_match = re.search(r"(?m)^.*LATITUDE.*LONGITUDE.*$", content)
            if coord_match:
                line = coord_match.group(0)
                s_lat = _fit_field(row.get(LAT_COLUMN), 8, 3)
                s_lon = _fit_field(row.get(LONG_COLUMN), 9, 3)
                s_elev = _fit_field(_wth_elev_default(), 4, 0)
                if s_lat is not None:
                    line = line.replace("LATITUDE", s_lat, 1)
                if s_lon is not None:
                    line = line.replace("LONGITUDE", s_lon, 1)
                if s_elev is not None:
                    line = line.replace("ELEV", s_elev, 1)
                content = content[:coord_match.start()] + line + content[coord_match.end():]

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

        # Genotype / support files (cultivar, ecotype, species, SDA, CO2, ...). Only
        # bundled when DSSATPRO can't resolve them or the user forces it; otherwise
        # skipped entirely (precomputed _SUPPORT_FILES is empty). Hard-link to save
        # disk + metadata; fall back to copy across filesystems.
        for fname in _SUPPORT_FILES:
            _link_or_copy(os.path.join(TEMPLATE_DIR, fname),
                          os.path.join(point_dir, fname), _SUPPORT_SYMLINK)


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
    # STEP 3 (continued) - EXECUTE DSSAT
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
            crop_cols = [c for c in ["crop_frac", "crop_pct", "crop_ha", "cell_ha"] if c in gridfile.columns]
            if crop_cols and POINT_ID_COLUMN in gridfile.columns:
                point_attrs = gridfile[[POINT_ID_COLUMN] + crop_cols].copy()
                point_attrs = point_attrs.rename(columns={POINT_ID_COLUMN: "point_id"})
                point_attrs["point_id"] = point_attrs["point_id"].astype(str)
                final_data["point_id"] = final_data["point_id"].astype(str)
                final_data = final_data.merge(point_attrs, on="point_id", how="left")
                if "cell_ha" in final_data.columns:
                    final_data["gridcell_area_ha"] = final_data["cell_ha"]
                if "crop_ha" in final_data.columns:
                    final_data["cropland_ha"] = final_data["crop_ha"]
                if {"final_grain_kg_ha", "cropland_ha"}.issubset(final_data.columns):
                    final_data["final_grain_production_kg"] = (
                        final_data["final_grain_kg_ha"] * final_data["cropland_ha"]
                    )
                if {"top_weight_kg_ha", "cropland_ha"}.issubset(final_data.columns):
                    final_data["top_weight_production_kg"] = (
                        final_data["top_weight_kg_ha"] * final_data["cropland_ha"]
                    )
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
        # HPC PREP MODE - folders built, no execution
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
    # STEP 4 - VISUALIZE RESULTS
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
                ax.set_title(f"Simulated Yield - {DSSAT_RUN_NAME}", fontsize=13)
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
