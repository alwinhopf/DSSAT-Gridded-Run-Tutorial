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

def _is_leap(year: int) -> bool:
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)


def create_grid_points(boundary_shape: gpd.GeoDataFrame,
                       spacing_m: int,
                       output_path: str) -> gpd.GeoDataFrame:
    """
    Create a regular grid of points inside *boundary_shape* at *spacing_m*
    metre spacing, re-project to WGS84, add LAT / LONG / ID columns, and
    write a shapefile. Mirrors ``create_grid_points`` in R.
    """
    ALBERS_CRS = "EPSG:5070"   # USA Contiguous Albers Equal Area (metres)

    projected = boundary_shape.to_crs(ALBERS_CRS)
    minx, miny, maxx, maxy = projected.total_bounds

    xs = np.arange(math.floor(minx), math.ceil(maxx) + spacing_m, spacing_m)
    ys = np.arange(math.floor(miny), math.ceil(maxy) + spacing_m, spacing_m)
    grid_pts = gpd.GeoDataFrame(
        geometry=[Point(x, y) for y in ys for x in xs],
        crs=ALBERS_CRS,
    )

    inside = gpd.sjoin(grid_pts, projected[["geometry"]], how="inner",
                       predicate="within").drop_duplicates("geometry")

    if inside.empty:
        raise RuntimeError("STEP 0 FAILED: No grid points created inside boundary.")

    inside = inside.to_crs("EPSG:4326").reset_index(drop=True)
    inside[LAT_COLUMN]      = inside.geometry.y.round(6)
    inside[LONG_COLUMN]     = inside.geometry.x.round(6)
    inside[POINT_ID_COLUMN] = [f"{i+1:08d}" for i in range(len(inside))]

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    inside.to_file(output_path)
    return inside


def load_existing_points(input_path: str, output_path: str) -> gpd.GeoDataFrame:
    """
    Load and normalise an existing point (or polygon) shapefile to the
    pipeline schema (WGS84, ID / LAT / LONG columns). Polygons are
    converted to centroids. Mirrors ``load_existing_points`` in R.
    """
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Existing point shapefile not found: {input_path}")

    pts = gpd.read_file(input_path)

    if not all(pts.geometry.geom_type.isin(["Point", "MultiPoint"])):
        print(f"Geometry is [{pts.geometry.geom_type.unique()}]; converting to centroids.")
        pts = pts.copy()
        pts["geometry"] = pts.geometry.centroid

    pts = pts.to_crs("EPSG:4326").reset_index(drop=True)
    pts[LAT_COLUMN]  = pts.geometry.y.round(6)
    pts[LONG_COLUMN] = pts.geometry.x.round(6)

    if POINT_ID_COLUMN not in pts.columns:
        pts[POINT_ID_COLUMN] = [f"{i+1:08d}" for i in range(len(pts))]
    else:
        ids = pts[POINT_ID_COLUMN].astype(str)
        bad = ids.isna() | (ids == "") | ids.duplicated()
        if bad.any():
            print("ID column has NA/blank/duplicates; regenerating sequential IDs.")
            pts[POINT_ID_COLUMN] = [f"{i+1:08d}" for i in range(len(pts))]
        else:
            unique_order = dict(zip(ids.unique(), range(len(ids.unique()))))
            pts[POINT_ID_COLUMN] = ids.map(lambda x: f"{unique_order[x]+1:08d}")

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    pts.to_file(output_path)
    return pts


# ---------------------------------------------------------------------------
# Weather extension helpers
# ---------------------------------------------------------------------------

def _read_wth_file(path: str):
    """
    Read a DSSAT .WTH file; return (header_lines, data_df) where data_df has
    a first column of integer date codes plus numeric data columns.
    """
    with open(path, "r") as fh:
        lines = fh.readlines()

    data_start = next(
        (i for i, ln in enumerate(lines) if re.match(r"^\s*\d", ln)),
        None
    )

    if data_start is None:
        return lines, None

    header_lines = [ln.rstrip() for ln in lines[:data_start]]
    data_lines   = [ln.rstrip() for ln in lines[data_start:]]

    clean  = [re.sub(r"\bNA\b|\bNaN\b", " -99.0 ", ln) for ln in data_lines]
    rows   = [ln.split() for ln in clean if ln.strip()]

    if not rows:
        return header_lines, None

    ncols  = max(len(r) for r in rows)
    padded = [r + ["-99.0"] * (ncols - len(r)) for r in rows]
    arr    = np.array(padded, dtype=float)
    df     = pd.DataFrame(arr)
    return header_lines, df


def _get_year_doy(date_code: float, year_format: str):
    dc  = int(date_code)
    if year_format == "YYDDD":
        yy   = dc // 1000
        doy  = dc % 1000
        year = (2000 + yy) if yy < 80 else (1900 + yy)
    else:
        year = dc // 1000
        doy  = dc % 1000
    return year, doy


def _make_date_code(year: int, doy: int, year_format: str) -> int:
    if year_format == "YYDDD":
        return (year % 100) * 1000 + doy
    return year * 1000 + doy


def _format_wth_row(date_code: int, vals: np.ndarray, year_format: str) -> str:
    fmt_date  = f"{date_code:07d}" if year_format == "YYYYDDD" else f"{date_code:05d}"
    num_parts = "".join(f"{v:6.1f}" for v in vals)
    return f"{fmt_date:>7s}{num_parts}".replace(" -99.0", "  -99")


def extend_weather_repeat_single_ignore_partial(
    path: str,
    ref_start_year: int,
    ref_end_year: int,
    target_end_year: int,
    verbose: bool = True,
) -> bool:
    """
    Extend a .WTH file to *target_end_year* by repeating a reference historic
    block (ref_start_year … ref_end_year), skipping partial/incomplete years.
    Handles leap-year day-count mismatches.
    Mirrors ``extend_weather_repeat_single_ignore_partial`` in R.
    """
    header_lines, df = _read_wth_file(path)
    if df is None or df.empty:
        return False

    sample      = int(df.iloc[0, 0])
    year_format = "YYDDD" if len(str(sample)) <= 5 else "YYYYDDD"

    df["_year"] = df.iloc[:, 0].apply(lambda x: _get_year_doy(x, year_format)[0])
    df["_doy"]  = df.iloc[:, 0].apply(lambda x: _get_year_doy(x, year_format)[1])

    year_counts    = df.groupby("_year").size()
    complete_years = [
        yr for yr, cnt in year_counts.items()
        if cnt == (366 if _is_leap(yr) else 365)
    ]

    if not complete_years:
        import warnings
        warnings.warn(f"No complete years in {path}; using first 365 rows as fallback.")
        chosen_ref = ref_end_year
        last_full  = int(df["_year"].iloc[0])
        df_trunc   = df.iloc[:365].copy()
    else:
        prior      = [y for y in complete_years if y <= ref_end_year]
        chosen_ref = max(prior) if prior else max(complete_years)
        last_full  = max(complete_years)
        df_trunc   = df[df["_year"] <= last_full].copy()

    if last_full >= target_end_year:
        return True   # already covers target

    body_cols = [c for c in df_trunc.columns if c not in ("_year", "_doy")]

    ref_block = df[df["_year"] == chosen_ref][body_cols].copy().reset_index(drop=True)
    if ref_block.empty:
        ref_block = df_trunc[body_cols].iloc[:365].copy().reset_index(drop=True)

    ref_doys        = ref_block.iloc[:, 0].apply(lambda x: _get_year_doy(x, year_format)[1])
    ref_yr_in_block = _get_year_doy(ref_block.iloc[0, 0], year_format)[0]
    ref_dates       = pd.to_datetime(ref_yr_in_block * 1000 + ref_doys.values, format="%Y%j")
    ref_mmdd        = ref_dates.strftime("%m-%d")

    base_df = df_trunc[body_cols].copy().reset_index(drop=True)

    added = []
    for tgt_year in range(last_full + 1, target_end_year + 1):
        n_days    = 366 if _is_leap(tgt_year) else 365
        tgt_dates = pd.date_range(f"{tgt_year}-01-01", periods=n_days, freq="D")
        tgt_mmdd  = tgt_dates.strftime("%m-%d")

        idx_rows = []
        for mm in tgt_mmdd:
            cands = np.where(ref_mmdd == mm)[0]
            if len(cands):
                idx_rows.append(int(cands[0]))
            elif mm == "02-29":
                c228 = np.where(ref_mmdd == "02-28")[0]
                idx_rows.append(int(c228[0]) if len(c228) else 0)
            else:
                fallback = 0
                for k in range(1, 6):
                    prev = (tgt_dates[tgt_mmdd.tolist().index(mm)] -
                            pd.Timedelta(days=k)).strftime("%m-%d")
                    pc   = np.where(ref_mmdd == prev)[0]
                    if len(pc):
                        fallback = int(pc[0])
                        break
                idx_rows.append(fallback)

        blk            = ref_block.iloc[idx_rows].copy().reset_index(drop=True)
        new_codes      = [_make_date_code(tgt_year, doy + 1, year_format) for doy in range(n_days)]
        blk.iloc[:, 0] = new_codes
        blk            = blk.fillna(-99.0)
        added.append(blk)

    final_df = pd.concat([base_df] + added, ignore_index=True) if added else base_df

    date_codes = final_df.iloc[:, 0].astype(int)
    val_cols   = final_df.iloc[:, 1:]
    date_strs  = date_codes.apply(
        lambda dc: (f"{dc:07d}" if year_format == "YYYYDDD" else f"{dc:05d}")
    )
    val_strs = val_cols.apply(
        lambda col: col.fillna(-99.0).map(lambda v: f"{v:6.1f}"), axis=0
    )
    rows_out = date_strs.str.rjust(7) + val_strs.apply(lambda r: "".join(r.values), axis=1)
    rows_out = rows_out.str.replace(" -99.0", "  -99")

    with open(path, "w") as fh:
        fh.write("\n".join(header_lines) + "\n")
        fh.write("\n".join(rows_out.tolist()) + "\n")

    if verbose:
        print(f"Wrote extended file: {path} up to {target_end_year}")
    return True


# ---------------------------------------------------------------------------
# DSSAT batch / execution helpers
# ---------------------------------------------------------------------------

def _write_dssbatch(experiment_file: str, trtno_list: list,
                    batch_path: str, run_mode: str = "experiment") -> None:
    """
    Write a DSSBatch.V48 file for experiment mode.
    Mirrors R DSSAT::write_dssbatch().
    """
    mode_tag = "EXPERIMENT" if run_mode == "experiment" else "SEQUENCE"
    header   = (
        f"$BATCH({mode_tag})\n"
        "!\n"
        "@ FILEX                                                                                        "
        "TRTNO RP SQ OP CO\n"
    )
    fname = os.path.basename(experiment_file)
    lines = []
    for trt in trtno_list:
        filex_padded = f" {fname:<93s}"
        lines.append(f"{filex_padded}{trt:6d}  1  0  1  0")

    with open(batch_path, "w") as fh:
        fh.write(header)
        fh.write("\n".join(lines) + "\n")


def _write_dssbatch_sequence(experiment_file: str, trt: int,
                              seq_start: int, seq_end: int,
                              batch_path: str) -> None:
    """Write a DSSBatch.V48 for sequence mode."""
    fname  = os.path.basename(experiment_file)
    header = (
        "$BATCH(SEQUENCE)\n"
        "!\n"
        "@ FILEX                                                                                        "
        "TRTNO RP SQ OP CO\n"
    )
    lines = []
    for sq in range(seq_start, seq_end + 1):
        filex_padded = f" {fname:<93s}"
        lines.append(f"{filex_padded}{trt:6d}  1{sq:6d}  1  0")

    with open(batch_path, "w") as fh:
        fh.write(header)
        fh.write("\n".join(lines) + "\n")


def _run_dssat(run_dir: str, exe: str, run_mode_flag: str = "A",
               filex: str = "") -> None:
    """Execute DSSAT in *run_dir* using *exe*.

    Run-mode mapping (from DSSAT help):
      A  FileX      – Run all treatments in the specified FileX.
      B  BatchFile  – Batch: batchfile lists experiments and treatments.
      Q  BatchFile  – Sequence analysis.

    Mode A passes the FileX directly (e.g. UFGA8201.MZX); mode B/Q passes
    DSSBatch.V48.  Note: passing DSSBatch.V48 to mode A silently fails
    (DSSAT treats it as a FileX and can't find the *TREATMENTS section).
    """
    if run_mode_flag in ("B", "Q", "N", "S"):
        arg = "DSSBatch.V48"
    else:
        # mode A (and C, G, …): pass the FileX filename directly
        arg = filex if filex else "DSSBatch.V48"
    cmd = [exe, run_mode_flag, arg]
    subprocess.run(cmd, cwd=run_dir, check=False,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def _read_csv_safe(path: str) -> Optional[pd.DataFrame]:
    """Read a CSV, return None if missing or empty."""
    if not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path, index_col=False)
        return df if not df.empty else None
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Supplemental-file merge helper (shared by both experiment and sequence mode)
# ---------------------------------------------------------------------------

def _merge_supplemental(point_dir: str, master_runs: pd.DataFrame) -> pd.DataFrame:
    """
    Read soilorg.csv, soilni.csv, soilwat.csv from *point_dir* and
    left-join their summarised columns onto *master_runs* (keyed on RUNNO).
    Returns an augmented copy of master_runs.

    Column mapping mirrors the R pipeline exactly:
      soilorg  → SOMCT_start (first), SOMCT_end (last)  per RUN
      soilni   → NAPC (last), NLCC (last), NI#M (last)  per RUN
      soilwat  → IR#C (last), IRRC (last)               per RUN
    """
    mr = master_runs.copy()

    # --- soilorg ---
    soilorg = _read_csv_safe(os.path.join(point_dir, "soilorg.csv"))
    if soilorg is not None and "RUN" in soilorg.columns and "SOMCT" in soilorg.columns:
        so = (soilorg.groupby("RUN")
                     .agg(SOMCT_start=("SOMCT", "first"),
                          SOMCT_end=("SOMCT", "last"))
                     .reset_index()
                     .rename(columns={"RUN": "RUNNO"}))
        mr = mr.merge(so, on="RUNNO", how="left")
    else:
        mr["SOMCT_start"] = None
        mr["SOMCT_end"]   = None

    # --- soilni ---
    soilni = _read_csv_safe(os.path.join(point_dir, "soilni.csv"))
    if soilni is not None and "RUN" in soilni.columns:
        agg_dict = {}
        if "NAPC" in soilni.columns:
            agg_dict["NAPC"] = ("NAPC", "last")
        if "NLCC" in soilni.columns:
            agg_dict["NLCC"] = ("NLCC", "last")
        if "NI#M" in soilni.columns:
            agg_dict["NIM"]  = ("NI#M", "last")
        if agg_dict:
            sn = (soilni.groupby("RUN")
                        .agg(**agg_dict)
                        .reset_index()
                        .rename(columns={"RUN": "RUNNO"}))
            mr = mr.merge(sn, on="RUNNO", how="left")
    if "NAPC" not in mr.columns: mr["NAPC"] = None
    if "NLCC" not in mr.columns: mr["NLCC"] = None
    if "NIM"  not in mr.columns: mr["NIM"]  = None

    # --- soilwat ---
    soilwat = _read_csv_safe(os.path.join(point_dir, "soilwat.csv"))
    if soilwat is not None and "RUN" in soilwat.columns:
        agg_dict = {}
        if "IR#C" in soilwat.columns:
            agg_dict["IRC"]  = ("IR#C", "last")
        if "IRRC" in soilwat.columns:
            agg_dict["IRRC"] = ("IRRC", "last")
        if agg_dict:
            sw = (soilwat.groupby("RUN")
                         .agg(**agg_dict)
                         .reset_index()
                         .rename(columns={"RUN": "RUNNO"}))
            mr = mr.merge(sw, on="RUNNO", how="left")
    if "IRC"  not in mr.columns: mr["IRC"]  = None
    if "IRRC" not in mr.columns: mr["IRRC"] = None

    return mr


def _build_result_rows(ID: str, summary: pd.DataFrame,
                       mr: pd.DataFrame) -> pd.DataFrame:
    """
    Assemble the full 27-column result DataFrame from a summary table and the
    augmented master_runs table (mr). Shared by both experiment and sequence mode.
    """
    mr_idx = mr.set_index("RUNNO") if "RUNNO" in mr.columns else mr

    def _reindex_col(col_name: str) -> np.ndarray:
        if col_name in mr_idx.columns:
            return mr_idx[col_name].reindex(summary["RUNNO"]).values
        return np.full(len(summary), None)

    def _reindex_numeric(col_name: str) -> np.ndarray:
        """Reindex *col_name* from mr_idx, coercing None/object → float NaN."""
        if col_name not in mr_idx.columns:
            return np.full(len(summary), np.nan)
        return (pd.to_numeric(mr_idx[col_name].reindex(summary["RUNNO"]),
                              errors="coerce")
                  .to_numpy(dtype=float, na_value=np.nan))

    somct_start = _reindex_numeric("SOMCT_start")
    somct_end   = _reindex_numeric("SOMCT_end")
    somct_delta = somct_end - somct_start   # NaN − NaN = NaN, never raises

    # In DSSAT mode-A, TRNO is not reliably populated (always reports 1).
    # Use RUNNO as the treatment identifier — it correctly increments 1-N.
    _trno = summary.get("TRNO")
    _trno_all_one = (
        _trno is not None
        and hasattr(_trno, "nunique")
        and _trno.nunique() == 1
    )
    treatment_col = summary["RUNNO"] if _trno_all_one else _trno

    # DSSAT maps XCRD→LONG and YCRD→LAT in its summary output. After our
    # coordinate patch the summary LAT/LONG columns are now correct.
    return pd.DataFrame({
        "point_id":                               ID,
        "run_number":                             summary["RUNNO"],
        "treatment":                              treatment_col,
        "crop_code":                              summary.get("CR"),
        "latitude":                               summary.get("LAT"),
        "longitude":                              summary.get("LONG"),
        "weather_station_id":                     summary.get("WSTA"),
        "soil_profile_id":                        summary.get("SOIL_ID"),
        "dssat_file_id":                          summary.get("EXNAME"),
        "dssat_description":                      summary.get("TNAM"),
        "planting_date":                          summary.get("PDAT"),
        "emergence_date":                         summary.get("EDAT"),
        "harvest_date":                           summary.get("HDAT"),
        "year_planting":                          pd.to_numeric(summary.get("PYEAR"), errors="coerce"),
        "year_harvest":                           summary.get("HYEAR"),
        "top_weight_kg_ha":                       summary.get("CWAM"),
        "final_grain_kg_ha":                      summary.get("HWAM"),
        "removed_residue_kg_ha":                  summary.get("BWAH"),
        "soil_organic_carbon_start_kg_C_ha":      somct_start,
        "soil_organic_carbon_end_kg_C_ha":        somct_end,
        "soil_organic_carbon_delta_kg_C_ha":      somct_delta,
        "final_irrigation_applications_count":    _reindex_col("IRC"),
        "final_irrigation_amount_mm":             _reindex_col("IRRC"),
        "inorganic_n_applied_count":              _reindex_col("NIM"),
        "inorganic_n_applied_kg_ha":              _reindex_col("NAPC"),
        "nitrate_leaching_kg_ha":                 _reindex_col("NLCC"),
        "cumulative_net_co2_emissions_kg_CO2_ha": summary.get("CO2EM"),
        "cumulative_n2o_emissions_kg_N_ha":       summary.get("N2OEM"),
    })


# ---------------------------------------------------------------------------
# Per-point simulation
# ---------------------------------------------------------------------------

def _run_simulation(ID: str,
                    points_row: pd.Series,
                    dssat_run_dir: str) -> Optional[pd.DataFrame]:
    """
    Build (if needed) and run a DSSAT simulation for one grid point.
    Returns a DataFrame of results or None on error.
    Handles both RUN_MODE == "experiment" and RUN_MODE == "sequence".
    """
    point_dir = os.path.join(dssat_run_dir, ID)
    os.makedirs(point_dir, exist_ok=True)

    ext_map   = {"MZ": "MZX", "WH": "WHX", "SB": "SBX", "SC": "SCX",
                 "BA": "BAX", "SG": "SGX", "RI": "RIX"}
    ext       = ext_map.get(CROP_EXTENSION, TEMPLATE_FILE_NAME.rsplit(".", 1)[-1])
    exp_fname = f"{TEMPLATE_FILE_NAME.rsplit('.', 1)[0]}.{ext}"
    exp_path  = os.path.join(point_dir, exp_fname)

    results_template = {
        "point_id": [], "run_number": [], "treatment": [], "crop_code": [],
        "latitude": [], "longitude": [], "weather_station_id": [],
        "soil_profile_id": [], "dssat_file_id": [], "dssat_description": [],
        "planting_date": [], "emergence_date": [], "harvest_date": [],
        "year_planting": [], "year_harvest": [],
        "top_weight_kg_ha": [], "final_grain_kg_ha": [], "removed_residue_kg_ha": [],
        "soil_organic_carbon_start_kg_C_ha": [], "soil_organic_carbon_end_kg_C_ha": [],
        "soil_organic_carbon_delta_kg_C_ha": [],
        "final_irrigation_applications_count": [], "final_irrigation_amount_mm": [],
        "inorganic_n_applied_count": [], "inorganic_n_applied_kg_ha": [],
        "nitrate_leaching_kg_ha": [],
        "cumulative_net_co2_emissions_kg_CO2_ha": [],
        "cumulative_n2o_emissions_kg_N_ha": [],
    }
    results = pd.DataFrame(results_template)

    try:
        # ------------------------------------------------------------------ #
        # EXPERIMENT MODE                                                      #
        # ------------------------------------------------------------------ #
        if RUN_MODE == "experiment":
            batch_path = os.path.join(point_dir, "DSSBatch.V48")
            _write_dssbatch(exp_path,
                            list(range(TREATMENT_START, TREATMENT_END + 1)),
                            batch_path, run_mode="experiment")
            _run_dssat(point_dir, DSSAT_EXE_PATH, "A", filex=exp_fname)

            summary = _read_csv_safe(os.path.join(point_dir, "summary.csv"))
            if summary is None or summary.empty:
                trts  = list(range(TREATMENT_START, TREATMENT_END + 1))
                n_yr  = max(1, WEATHER_END_YEAR - WEATHER_START_YEAR)
                summary = pd.DataFrame({
                    "RUNNO":   range(1, len(trts) * n_yr + 1),
                    "TRNO":    [t for t in trts for _ in range(n_yr)],
                    "PYEAR":   [y for _ in trts
                                  for y in range(WEATHER_START_YEAR, WEATHER_END_YEAR)],
                    "CR": None, "LAT": None, "LONG": None, "WSTA": None,
                    "SOIL_ID": None, "EXNAME": None, "TNAM": None,
                    "PDAT": None, "EDAT": None, "HDAT": None, "HYEAR": None,
                    "CWAM": None, "HWAM": None, "BWAH": None,
                    "CO2EM": None, "N2OEM": None,
                })
            else:
                summary["PYEAR"] = summary["PDAT"].astype(str).str[:4]

            master_runs = summary[["RUNNO"]].copy()
            master_runs = _merge_supplemental(point_dir, master_runs)
            run_results = _build_result_rows(ID, summary, master_runs)
            results     = pd.concat([results, run_results], ignore_index=True)

        # ------------------------------------------------------------------ #
        # SEQUENCE MODE                                                        #
        # ------------------------------------------------------------------ #
        elif RUN_MODE == "sequence":
            all_seq_results = []

            for trt in range(TREATMENT_START, TREATMENT_END + 1):
                batch_path = os.path.join(point_dir, "DSSBatch.V48")
                _write_dssbatch_sequence(exp_path, trt,
                                         SEQUENCE_START, SEQUENCE_END,
                                         batch_path)
                _run_dssat(point_dir, DSSAT_EXE_PATH, "Q")

                summary = _read_csv_safe(os.path.join(point_dir, "summary.csv"))
                if summary is None or summary.empty:
                    n_seq   = SEQUENCE_END - SEQUENCE_START + 1
                    summary = pd.DataFrame({
                        "RUNNO":   range(1, n_seq + 1),
                        "TRNO":    trt,
                        "PYEAR":   [y for y in range(WEATHER_START_YEAR,
                                                      WEATHER_START_YEAR + n_seq)],
                        "CR": None, "LAT": None, "LONG": None, "WSTA": None,
                        "SOIL_ID": None, "EXNAME": None, "TNAM": None,
                        "PDAT": None, "EDAT": None, "HDAT": None, "HYEAR": None,
                        "CWAM": None, "HWAM": None, "BWAH": None,
                        "CO2EM": None, "N2OEM": None,
                    })
                else:
                    summary["PYEAR"] = summary["PDAT"].astype(str).str[:4]
                    if "TRNO" not in summary.columns or summary["TRNO"].isna().all():
                        summary["TRNO"] = trt

                master_runs = summary[["RUNNO"]].copy()
                master_runs = _merge_supplemental(point_dir, master_runs)
                trt_results = _build_result_rows(ID, summary, master_runs)
                all_seq_results.append(trt_results)

            if all_seq_results:
                results = pd.concat([results] + all_seq_results, ignore_index=True)

        # ------------------------------------------------------------------ #
        # Patch coordinates and save per-point CSV                            #
        # ------------------------------------------------------------------ #
        # DSSAT's summary.csv LAT/LONG columns are unreliable (DSSAT stores
        # the WTH LONG value in its summary LAT column due to a coordinate
        # convention quirk). Override with the authoritative gridpoint
        # coordinates from points_row, which come directly from our shapefile.
        if not results.empty:
            try:
                results["latitude"]  = float(points_row.get("LAT",  np.nan))
                results["longitude"] = float(points_row.get("LONG", np.nan))
            except Exception:
                pass
            out_csv = os.path.join(point_dir, f"results_{ID}.csv")
            results.to_csv(out_csv, index=False, na_rep="")
        return results

    except Exception as exc:
        print(f"--- FATAL ERROR processing ID {ID}: {exc} ---")
        return None


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

        def _run_one(ID: str) -> Optional[pd.DataFrame]:
            row = points[points[POINT_ID_COLUMN] == ID].iloc[0]
            return _run_simulation(ID, row, DSSAT_RUN_DIR)

        from concurrent.futures import ProcessPoolExecutor, as_completed

        # On macOS/Windows the default start method is 'spawn', which re-imports
        # this module in each worker and triggers the pool submission code again.
        # Using 'fork' avoids that for a pure command-line pipeline.
        import multiprocessing as _mp
        _mp_ctx = _mp.get_context("fork") if hasattr(_mp, "get_context") else None

        all_results = []
        with ProcessPoolExecutor(max_workers=DSSAT_CORES,
                                 mp_context=_mp_ctx) as pool:
            futures = {pool.submit(_run_one, ID): ID for ID in ids_to_run}
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
            print(f"Results combined → {FINAL_RESULTS_PATH}")
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
        print(f"Metadata written → {metadata_path}")

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
                print(f"Yield map saved → {FINAL_PLOT_PATH}")

        except Exception as exc:
            print(f"Step 4 visualization failed: {exc}")