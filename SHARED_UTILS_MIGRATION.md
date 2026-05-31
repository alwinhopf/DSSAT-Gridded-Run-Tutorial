# Shared Utilities Migration — Playbook & State

**Goal:** Extract the duplicated weather + soil download functions out of
`DSSAT Gridded Run Tutorial` (source of truth) and `DSSAT ML Phenology Prediction`
into a dedicated **private** shared package `dssat-spatial-utils`, then have both
repos depend on a pinned version.

This file is a resumable handoff: it captures the plan, the confirmed file/function
inventory, the exact commands, and the code-review findings so far. Safe to delete
once migration is complete.

---

## Decisions locked in
- **Source of truth:** `DSSAT Gridded Run Tutorial` (it's the superset — more
  sources, R+Python parity). Copy FROM here.
- **Shared repo:** `dssat-spatial-utils`, **private**, owner `alwinhopf`.
- **Scope:** weather + soil only (NOT landcover — see open question below).
- **Architecture:** dedicated package repo (NOT submodule, NOT monorepo). The
  Gridded repo is too heavy with data dirs (SoilGrids/, *_cache/, results/) to be
  installed directly, which is exactly why a lean package repo is correct.

## Environment facts (this machine)
- `brew` ✅ `/opt/homebrew/bin/brew` · `git` ✅ · `R`/`Rscript` ✅ `/usr/local/bin`
  · `python3` ✅ miniconda
- `gh` ❌ NOT installed. Either `brew install gh` (then `gh auth login`) to script
  the private-repo creation, OR create the empty private repo in the GitHub UI and
  add it as a remote manually.
- Both existing repos use HTTPS remotes under `github.com/alwinhopf/...` and have
  **unrelated uncommitted changes** (.DS_Store, README, pipeline tweaks) — do NOT
  commit those as part of this migration.

---

## Confirmed inventory (from source repo)

### R — public entry points to EXPORT (NAMESPACE)
| File (r_scripts/) | Exported function |
|---|---|
| weather_daymet.R | `process_weather_daymet` |
| weather_gridmet.R | `process_weather_gridmet` |
| weather_nasapower.R | `process_weather_nasapower` |
| weather_agera5.R | `process_weather_agera5` |
| weather_openmeteo.R | `process_weather_openmeteo` |
| weather_nasapower_chirps.R | `process_weather_nasapower_chirps` |
| soil_soilgrids.R | `process_soils_soilgrids` |
| soil_soilgrids_online.R | `process_soils_soilgrids_online` |
| soil_ssurgo.R | `process_soils_ssurgo` |
| soil_hwsd.R | `process_soils_hwsd` |

R internal helpers (keep in package namespace, do NOT export): `robust_SDA_query`,
`robust_SDA_spatialQuery`, `calculate_soil_properties`, `format_dssat_soil_single`
(soil_ssurgo.R); `calculate_soil_physics`, `format_dssat_sol_file`,
`fetch_soilgrids_rest`, `fetch_soilgrids_vrt` (soil_soilgrids_online.R).

### Python — public API for `__init__.py` (process_* only; `_`-prefixed stay private)
Same 10 names as above: `process_weather_{daymet,gridmet,nasapower,agera5,openmeteo,nasapower_chirps}`,
`process_soils_{soilgrids,soilgrids_online,ssurgo,hwsd}`.

### R library deps to put in DESCRIPTION Imports
daymetr, nasapower, soilDB, sf, terra, ncdf4, httr, jsonlite, dplyr, tidyr,
lubridate, doParallel, foreach, pbapply, DSSAT
Suggests (optional sources): ecmwfr (AgERA5)

### Python deps (already in source repo requirements.txt — reuse)
numpy, pandas, geopandas, shapely, pyproj, fiona, rasterio, xarray, netCDF4,
requests, tqdm, PyYAML; optional: cdsapi (AgERA5)

---

## Step-by-step

### Phase 1 — Build `dssat-spatial-utils` (additive, low risk)
```
ROOT=/Users/alwinhopf/Documents/GitHub
SRC="$ROOT/DSSAT Gridded Run Tutorial"
PKG="$ROOT/dssat-spatial-utils"

mkdir -p "$PKG/R" "$PKG/python/dssatutils" "$PKG/tests" "$PKG/.github/workflows"

# Copy R sources (weather+soil only)
cp "$SRC"/r_scripts/weather_*.R "$SRC"/r_scripts/soil_*.R "$PKG/R/"
# Copy Python sources
cp "$SRC"/python_scripts/weather_*.py "$SRC"/python_scripts/soil_*.py "$PKG/python/dssatutils/"
# Reuse infra
cp "$SRC"/tests/test_*.py "$PKG/tests/"
cp "$SRC"/.github/workflows/smoke.yml "$PKG/.github/workflows/"
cp "$SRC"/requirements.txt "$PKG/requirements.txt"
```
Then create (contents below): `$PKG/DESCRIPTION`, `$PKG/NAMESPACE`,
`$PKG/python/dssatutils/__init__.py`, `$PKG/pyproject.toml`, `$PKG/README.md`,
`$PKG/.gitignore`, `$PKG/LICENSE` (MIT).

`git init`, commit, create **private** GitHub repo, push, then `git tag v0.1.0 && git push --tags`.

### Phase 2 — Repo 1 (ML Phenology) consumes the package
- Delete duplicated `r_scripts/weather_*.R` and `r_scripts/soil_*.R`.
- In `pipeline/config.R` add a guarded `remotes::install_github("alwinhopf/dssat-spatial-utils@v0.1.0")` + `library(dssatutils)`.
- Replace any `source("r_scripts/weather_*.R")` with `library(dssatutils)`.
- `renv::init(); renv::snapshot()` (repo currently has NO lock files).

### Phase 3 — Repo 2 (Gridded, this repo) consumes the package
- Delete `r_scripts/` + `python_scripts/` weather/soil files (keep landcover unless moved).
- `dssat_main_pipeline.R`: replace `source("r_scripts/...")` with `library(dssatutils)`.
- `dssat_main_pipeline.py`: replace `from python_scripts.weather_* import ...` with `from dssatutils import ...`.
- R: add install line to `setup_renv.R`. Python: add to `requirements.txt`:
  `dssatutils @ git+https://github.com/alwinhopf/dssat-spatial-utils.git@v0.1.0`

### Phase 4 — Versioning going forward
Branch in shared repo → CI smoke tests → merge → tag `vX.Y.Z` → bump pin in both
consumer repos. Both pin to a TAG, never `main`, so neither breaks unexpectedly.

---

## Code review findings (weather R files read so far)

Reviewed: weather_daymet.R, weather_gridmet.R, weather_nasapower.R,
weather_agera5.R, weather_openmeteo.R. **Still to review:** all soil R files,
weather_nasapower_chirps.R, and all Python files.

1. **GridMET RH2M / TDEW are crude proxies** (`weather_gridmet.R:188-189`):
   `TDEW = TMIN - 2.5` and `RH2M = 100 - (TMAX-TMIN)*2` clamped to [20,100]. These
   are rough estimates, not measured — fine for DSSAT runs but should be documented
   so downstream users don't treat them as observations.
2. **Open-Meteo emits TDEW=-99 and RH2M=-99** (`weather_openmeteo.R:89-90`): no
   daily dewpoint/RH from that API. DSSAT tolerates -99, but ET methods needing RH
   will silently degrade. Worth a note in the source's docstring (already partly noted).
3. **`weather_nasapower.R` has leftover "← FIX 1..4" comments** marking past
   bug fixes (unclosed blocks). Harmless but clutter — clean up during the move.
4. **Daymet leap-year handling duplicates DOY 365 as 366** (`weather_daymet.R:73-82`)
   rather than interpolating — acceptable approximation, keep the comment.
5. **Inconsistent TAV/AMP computation:** gridmet uses `DSSAT::calc_TAV/ calc_AMP`,
   while daymet/nasapower/openmeteo/agera5 hand-roll the monthly-mean amplitude.
   Results are close but not identical. Consider a single shared helper in the
   package so all sources agree. (Python side already has `_calc_tav/_calc_amp`
   duplicated per-module — same consolidation opportunity.)
6. **`library(...)` calls at top of each R file** must become `@importFrom` /
   DESCRIPTION Imports when packaged (R CMD check rejects bare `library()` in a
   package). Minor mechanical change during Phase 1.

(Complete the soil + Python review in the clean implementation pass.)

---

## Open questions for the user
1. **Landcover** (`landcover_raster.*`, `landcover_raster_to_gridpoints.*`) is ALSO
   duplicated across both repos but is out of the stated "weather + soil" scope.
   Move it into the shared package too, or leave it?
2. **Package name:** `dssatutils` (R+Python importable name) under repo
   `dssat-spatial-utils` — OK, or prefer a different name?
