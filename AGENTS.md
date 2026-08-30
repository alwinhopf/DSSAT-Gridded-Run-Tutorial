# AGENTS.md — DSSAT_Gridded_Run_Tutorial

> **Workspace context:** Read the root [`../AGENTS.md`](../AGENTS.md) first. This document
> holds rules and guidance specific to the `DSSAT_Gridded_Run_Tutorial` repository.

## 1. Role in the Workspace

`DSSAT_Gridded_Run_Tutorial` is the **canonical gridded crop modeling engine driver**:
- It consumes the shared foundation libraries **`dssatutils`** (weather and soil download) and **`dssatengine`** (grid generation, batch writing, DSSAT execution, output parsing).
- It provides the reference end-to-end gridded workflows (`dssat_main_pipeline.py` and `dssat_main_pipeline.R`).
- Downstream study repositories (`Bioenergy_Model_Input_Comparison`, `DSSAT_SubField_MILP_Analysis`, etc.) drive this engine via `ENGINE_DIR` or sibling imports — they do **not** fork it.

## 2. Critical Domain Rules & Common Traps

### A. 8-Character FileX Basename Rule
DSSAT Fortran (`CSM.for`) expects FileX experiment file basenames to be **exactly 8 characters** (plus a 3-character extension ending in `X`, e.g. `UFJA1803.WHX` or `00010001.MZX`).
- Longer base names crash DSSAT with "File not found" or corrupted character buffer errors.
- Both pipeline twins auto-map grid point IDs to canonical 8-character experiment file names.

### B. Artifact Location Discipline
- The code repository holds code, templates (`dssat_templates/`), and static shapefiles (`shapefile/`).
- **Never commit generated artifacts:** `gridpoints/`, `weather/`, `soil/`, `dssat_runs/`, `results/`, `*_cache/`, `*.OUT`, or `.migration_backup_*`.
- Set `output_root_dir` (or `*_ROOT` / `*_dir` keys in `config.yml`) to route generated simulation artifacts into the calling study's folder.
- For large spatial sweeps, enable `cleanup_run_folders: true` to prevent disk exhaustion.

### C. Coordinate Reference Systems (CRS)
- Use **EPSG:4326 (WGS84)** for geographic API queries (e.g. Daymet, NASA POWER) and point coordinate exchange.
- Use an **equal-area projected CRS** (e.g. `EPSG:6933` equal-area cylindrical, or `EPSG:5070` CONUS Albers) when generating regular metric grids (`GRID_SPACING_METERS`) and computing cropland fractions or cell areas.

### D. Cropland Masks & Spatial Modes
- **Mode A:** Regular grid clipped to a boundary polygon.
- **Mode B:** User-supplied point/polygon shapefile (`USE_EXISTING_POINT_SHAPEFILE`).
- **Mode C:** Cropland-filtered points built using `r_scripts/landcover_raster.R` / `python_scripts/landcover_raster.py` (CDL / NLCD class 82).

## 3. R ↔ Python Parity Contract

`dssat_main_pipeline.py` and `dssat_main_pipeline.R` are twins:
- Any configuration key, command-line argument, or step logic added to one must be added to the other.
- Configuration loading is handled by `config_loader.py` and `config_loader.R`. Both read the same `config.yml` (and optional overlay config).
- Test cross-language invariants with `pytest tests/test_cross_language_parity.py`.

## 4. Verification & Testing

Before completing any changes, run the test suites:

```bash
# Python smoke & cross-language parity tests
pytest tests/test_smoke.py
pytest tests/test_cross_language_parity.py

# Comprehensive E2E tests (when DSSAT executable is configured)
pytest tests/test_e2e_comprehensive.py
Rscript -e "testthat::test_file('tests/test_e2e_comprehensive.R')"
```
