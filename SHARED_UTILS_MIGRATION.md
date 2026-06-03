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

### ⚠️ CRITICAL — Python cross-module imports MUST be fixed when packaging
Two modules import private helpers from sibling modules using FLAT imports that
only work when all files sit in one directory on `sys.path`. Inside a package
these break with `ModuleNotFoundError`. Convert to RELATIVE imports during the copy:
- `soil_hwsd.py`:
  `from soil_soilgrids_online import _calculate_soil_physics, _format_dssat_sol_file`
  → `from .soil_soilgrids_online import _calculate_soil_physics, _format_dssat_sol_file`
- `weather_nasapower_chirps.py`:
  `from weather_nasapower import _fetch_nasa_power, _calc_tav, _calc_amp`
  → `from .weather_nasapower import _fetch_nasa_power, _calc_tav, _calc_amp`
(Implication: CHIRPS reuses NASA-POWER internals; HWSD reuses SoilGrids-online
internals — so even a "weather only" or "soil only" partial move is not possible;
move all 10 together.) Re-grep for any other `from (weather_|soil_)` flat imports
before finalizing.

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

## Resolved decisions (user-confirmed 2026-05-31)
1. **Landcover stays put** — `landcover_raster.*` and
   `landcover_raster_to_gridpoints.*` remain in the Gridded model pipeline's
   `r_scripts/` + `python_scripts/`. Do NOT move them into the shared package.
   (The R pipeline sources them at `dssat_main_pipeline.R:252-253` — leave those
   two `source()` lines intact; only the weather/soil `source()` lines 242-251
   get replaced by `library(dssatutils)`.)
2. **Repo name = `dssatutils`** (NOT `dssat-spatial-utils`). Importable name in
   both R and Python is also `dssatutils`. Install URLs:
   - R: `remotes::install_github("alwinhopf/dssatutils@v0.1.0")`
   - Python: `dssatutils @ git+https://github.com/alwinhopf/dssatutils.git@v0.1.0`

## Pipeline wiring facts
**Rewire by CONTENT MATCH, not line number** (line numbers reported inconsistently
this session due to corruption — do not trust them; match the actual text).
- `dssat_main_pipeline.R`: there is a block of 10 consecutive
  `source(file.path(SCRIPT_DIR, "weather_*.R" / "soil_*.R"))` calls, immediately
  FOLLOWED by 2 `source(... "landcover_raster.R")` /
  `"landcover_raster_to_gridpoints.R"` calls. Replace ONLY the 10 weather/soil
  `source()` lines with a single `library(dssatutils)`. KEEP the 2 landcover
  `source()` lines (landcover stays local). Also keep the SCRIPT_DIR bootstrap
  block earlier in the file (it uses `source(cand)` for config detection).
- `dssat_main_pipeline.py`: a block of 10 `from <module> import process_*`
  (weather_daymet, weather_nasapower, weather_gridmet, weather_openmeteo,
  weather_nasapower_chirps, weather_agera5, soil_ssurgo, soil_soilgrids,
  soil_soilgrids_online, soil_hwsd). Replace with `from dssatutils import (...)`.
  NOTE: there are ALSO two lazy/local imports deeper in the file —
  `import soil_soilgrids_online as _sg_mod` and
  `import weather_nasapower_chirps as _wc` — update these to
  `from dssatutils import soil_soilgrids_online as _sg_mod` etc.
  (i.e. `import dssatutils.soil_soilgrids_online as _sg_mod`).
- Landcover imports in the Python pipeline (if any) stay as-is.

## ✅ PHASE 2 + 3 EXECUTED (2026-06-01)
Both consumer repos rewired to the local dssatutils package. Staged (NOT committed)
so you can review/adjust. Backups in each repo's `.migration_backup_<ts>/`.

PHASE 2 — DSSAT ML Phenology Prediction:
- 6 dup files deleted from r_scripts/ (weather_{daymet,gridmet,nasapower}.R,
  soil_{soilgrids,soilgrids_online,ssurgo}.R); landcover kept.
- 4 source() sites rewired to suppressMessages(library(dssatutils)): utils.R (x2,
  inside local=TRUE blocks), 01_particle_filter.R, 04_cohesive_calibration.R,
  scratch/run_g29_debug.R.
- Guarded install+library prepended to pipeline/config.R (install_local, USE_LOCAL).
- dssatutils 0.1.0 built & installed (R CMD: * DONE), loads with both soil fns.
- renv::init() => INIT_OK; renv.lock has 131 pkgs incl dssatutils 0.1.0 (Source
  Local). .Rprofile + renv/ created. Separate snapshot() unnecessary (init captured
  everything incl dssatutils).

PHASE 3 — DSSAT Gridded Run Tutorial:
- 20 dup files deleted (10 r_scripts/ + 10 python_scripts/); landcover_* kept.
- R: 10 weather/soil source() lines commented, library(dssatutils) inserted.
- Py: 10 `from <mod> import` -> `from dssatutils.<mod> import`; 2 lazy imports ->
  `import dssatutils.soil_soilgrids_online as _sg_mod` / `...weather_nasapower_chirps as _wc`.
- requirements.txt += `-e /Users/.../dssatutils`; setup_renv.R += install_local.

### ⚠️ TWO BUGS FOUND IN PHASE 3 SCRIPT & FIXED BY HAND (also fixed in the scripts)
1. Regex `(weather|soil)_[a-z_]*\.R` did NOT match the DIGIT in "agera5", so
   `source(... "weather_agera5.R")` stayed ACTIVE while the file was deleted ->
   would crash. Fixed the orphaned line by hand; fixed sed in BOTH 02 & 03 scripts
   to `[a-z0-9_]*`.
2. SCRIPT_DIR auto-detection probed `file.exists(file.path(dir, "soil_ssurgo.R"))`
   — that file is now deleted, so detection would hit stop(). Repointed the
   sentinel to the still-present `landcover_raster.R` (and updated the error msg).

### Verification after fixes (all PASS)
- Rscript parse: dssat_main_pipeline.R + all 5 modified ML files = PARSE_OK.
- python3 -m py_compile dssat_main_pipeline.py = OK.
- All 10 `dssatutils.*` submodule specs resolve via find_spec.
- ML repo: no orphaned source()/file.exists() of deleted files (only backups match).
- Gridded: 20 deletions confirmed; landcover kept; all weather/soil sources commented.

### REMAINING (your call)
- Review diffs, then COMMIT in each repo (changes are staged, not committed).
- When the GitHub remote exists, change the local-path pins to the tag:
  config.R/setup_renv.R install_github("alwinhopf/dssatutils@v0.1.0");
  requirements.txt `dssatutils @ git+https://github.com/alwinhopf/dssatutils.git@v0.1.0`.
- Optional: `pip install -e /…/dssatutils` (or into the project env) so the Python
  pipeline import works at runtime; deep-run renv::restore() on another machine.
- The migrate/ scripts now have the corrected regex for any re-runs.

## PHASE 2/3 ATTEMPT (2026-05-31) — superseded by the EXECUTED section above
Display went blind mid-run; stopped before any destructive edits. NO changes were
made to either consumer repo. Findings that REQUIRED fixing the pre-written scripts:

- ML Phenology does NOT use literal `source("r_scripts/...")`. It sources via a
  variable, from multiple files, some with `local = TRUE`:
    01_particle_filter.R:36   source(file.path(R_SCRIPTS_DIR, "soil_soilgrids_online.R"))
    04_cohesive_calibration.R:23  source(file.path(R_SCRIPTS_DIR, "soil_soilgrids_online.R"))
    pipeline/utils.R:577      source(file.path(R_SCRIPTS_DIR, "soil_soilgrids.R"), local = TRUE)
    pipeline/utils.R:601      source(file.path(R_SCRIPTS_DIR, "soil_soilgrids_online.R"), local = TRUE)
    scratch/run_g29_debug.R:10  source(file.path(PROJECT_ROOT, "r_scripts/soil_soilgrids_online.R"))
  Only SOIL functions are used in the ML pipeline (process_soils_soilgrids,
  process_soils_soilgrids_online). The 3 weather scripts in its r_scripts/ appear
  unused by the pipeline but are duplicates -> still safe to delete in Phase 2.
  ML r_scripts/ has only 6 dup files (weather_{daymet,gridmet,nasapower}.R,
  soil_{soilgrids,soilgrids_online,ssurgo}.R) + landcover (keep).
- `02_phase2_ml_phenology.sh` was REWRITTEN to match these real patterns: it now
  replaces every weather_/soil_ source() line (R_SCRIPTS_DIR or r_scripts/ form,
  with/without `, local = TRUE`) with `suppressMessages(library(dssatutils))`,
  prepends a guarded install+library to pipeline/config.R, and deletes the 6 dup
  files (keeps landcover). Also fixed a portability bug: replaced `mapfile`
  (bash 4+) with a `while read` loop since macOS /bin/bash is 3.2.
- Both repos were essentially clean before stopping (ML had only an untracked
  scratch/dssat-csm-os/, unrelated; Gridded clean).
- `03_phase3_gridded.sh` was NOT changed this session; the Gridded pipeline uses
  the clean literal `source(file.path(SCRIPT_DIR, "weather_*.R"))` form that its
  sed already matches. Still review its diff after running.

### How to RUN phase 2/3 yourself (you can see output; tool display was glitching)
```
cd /Users/alwinhopf/Documents/GitHub/dssatutils
USE_LOCAL=1 ./migrate/02_phase2_ml_phenology.sh
git -C "/Users/alwinhopf/Documents/GitHub/DSSAT ML Phenology Prediction" diff   # review
# then in R, from the ML repo root:
#   renv::init(); renv::snapshot()
USE_LOCAL=1 ./migrate/03_phase3_gridded.sh
git -C "/Users/alwinhopf/Documents/GitHub/DSSAT Gridded Run Tutorial" diff      # review
```
Each script backs up touched files to `<repo>/.migration_backup_<ts>/`. To undo:
restore from there or `git checkout -- <file>` (edits are staged, not committed).

## ✅ DONE so far (local, no remote yet)
- Package built at `/Users/alwinhopf/Documents/GitHub/dssatutils` and VERIFIED while
  display was working: R/ has 10 files, python/dssatutils/ has 10 + __init__.py,
  relative cross-imports confirmed (`from .soil_soilgrids_online`,
  `from .weather_nasapower`), no flat cross-imports remain, smoke.yml present,
  `python3 -c "import dssatutils"` returns all 10 process_* in __all__.
- `git init` + 2 commits made locally (commit 1 = package; commit 2 = migrate/
  scripts). Tag `v0.1.0` created locally. NO remote yet (user chose skip-for-now).
- Ready-to-run migration scripts written under `dssatutils/migrate/`:
  - `00_create_remote.sh` (gh OR manual; creates PRIVATE repo, pushes code+tag)
  - `01_local_install_check.sh` (python import + R sys.source smoke test)
  - `02_phase2_ml_phenology.sh` (idempotent rewire + backups; USE_LOCAL=1 supported)
  - `03_phase3_gridded.sh` (idempotent rewire, keeps landcover, adds pins; USE_LOCAL=1)
  - `migrate/README.md` (run order + rollback)

## Verification results (TIER 1, while display worked)
- Python byte-compile of all 10 modules: PASS
- `import dssatutils` (lazy) exposes all 10 process_* in __all__: PASS
- R `parse()` of all 10 R files: PASS
- TIER 2 (real per-module Python import) needs numpy/pandas/geopandas/rasterio/
  xarray etc. The conda *base* env lacks them, so deep import is opt-in
  (`DEEP=1 ./migrate/01_local_install_check.sh`) — run it inside the project's
  env (the Gridded repo's environment.yml / requirements.txt) to exercise fully.
- Two bugs were found & fixed in the CHECK SCRIPT itself (not the package):
  R regex `\\.R$` -> `[.]R$`; and split syntax-check from heavy real-import.

## TO RUN when ready (needs the remote, or USE_LOCAL=1)
1. `cd dssatutils && ./migrate/01_local_install_check.sh`  (sanity)
2. `./migrate/00_create_remote.sh`  (after `brew install gh && gh auth login`, or make empty private repo in UI)
3. `./migrate/02_phase2_ml_phenology.sh` then review `git -C "../DSSAT ML Phenology Prediction" diff`; then `renv::init(); renv::snapshot()`
4. `./migrate/03_phase3_gridded.sh` then review `git -C "../DSSAT Gridded Run Tutorial" diff`
Each script backs up touched files to `<repo>/.migration_backup_<ts>/`.

## Code review — COMPLETE for R; Python partial
All 10 R files reviewed clean. Final two:
- `weather_nasapower_chirps.R`: solid. NASA-POWER fetched per point, RAIN replaced
  by CHIRPS within |lat|<=50 with per-day match + graceful fallback to NASA rain
  outside coverage / no-data cells. Notes: (a) serial loop — the `n_cores`
  argument is accepted for interface parity but UNUSED (get_power is the
  bottleneck & self-throttles); fine, just document. (b) relies on global
  `CHIRPS_RESOLUTION` via exists() — same pattern as soilgrids USE_REST_API;
  candidate to convert to a function arg later.
- `soil_hwsd.R`: solid. Tolerant SQLite column matching, dominant-component pick,
  lazy DBI/RSQLite load, points over no-data skipped with warning. IMPORTANT
  packaging note: it calls `calculate_soil_physics()` + `format_dssat_sol_file()`
  defined in soil_soilgrids_online.R — in an R PACKAGE all R/ files share one
  namespace so this resolves automatically (NO import fix needed, unlike the
  Python side). Minor: mapping CSV writes a row for every ID incl. skipped ones
  (SOIL_ID==ID) though skipped points have no .SOL — could mislead downstream.

Python: byte-compile PASS for all 10 (TIER 1a) and lazy `import dssatutils` PASS;
a deeper line-by-line Python read is the only remaining nice-to-have (blocked by
intermittent blind-tool display this session). Known structural facts already
captured: CHIRPS py reuses `_fetch_nasa_power/_calc_tav/_calc_amp` and HWSD py
reuses `_calculate_soil_physics/_format_dssat_sol_file` — both converted to
relative imports and byte-compile-verified.

## BUILD STATE — Phase 1 executed (2026-05-31)
The package was created at `/Users/alwinhopf/Documents/GitHub/dssatutils`.
Operations were issued while the tool DISPLAY was intermittently blind, so the
fresh session MUST verify the following actually landed (the file ops themselves
succeed even when output doesn't render):

Created / expected to exist:
- `dssatutils/R/` — 10 files (weather_{daymet,gridmet,nasapower,agera5,openmeteo,
  nasapower_chirps}.R, soil_{soilgrids,soilgrids_online,ssurgo,hwsd}.R)
- `dssatutils/python/dssatutils/` — same 10 as .py + `__init__.py` (lazy PEP-562
  re-export of the 10 process_* functions)
- `dssatutils/tests/` — test_global_sources.py, test_smoke.py
- `dssatutils/` root — DESCRIPTION, NAMESPACE (10 explicit export()),
  pyproject.toml (packages.find where=["python"]; extras: agera5/all/dev),
  requirements.txt, README.md, LICENSE (MIT), .gitignore
- `dssatutils/.github/workflows/` — dir created; smoke.yml copy NOT yet confirmed
  (re-copy from Gridded repo `.github/workflows/smoke.yml` if missing).

Python cross-import fixes applied via `sed -i ''`:
- soil_hwsd.py: `from soil_soilgrids_online import` -> `from .soil_soilgrids_online import`
- weather_nasapower_chirps.py: `from weather_nasapower import` -> `from .weather_nasapower import`
VERIFY both took (grep '^from \.'); if not, apply manually.

### VERIFY checklist (run first in fresh session)
```
PKG=/Users/alwinhopf/Documents/GitHub/dssatutils
ls "$PKG" "$PKG/R" "$PKG/python/dssatutils" "$PKG/tests" "$PKG/.github/workflows"
grep -n '^from \.' "$PKG/python/dssatutils/soil_hwsd.py" "$PKG/python/dssatutils/weather_nasapower_chirps.py"
python3 -c "import sys; sys.path.insert(0,'$PKG/python'); import dssatutils; print(dssatutils.__all__)"
Rscript -e 'cat(readLines(file.path("'"$PKG"'","NAMESPACE")), sep="\n")'
```
Also re-copy smoke.yml if missing; consider `pip install -e "$PKG"` and
`R CMD INSTALL "$PKG"` (or devtools::load_all) as smoke tests.

### REMAINING WORK (needs working display; do NOT do blind)
1. smoke.yml present? else copy it.
2. `cd dssatutils && git init && git add -A && git commit` (first commit).
3. Create PRIVATE GitHub repo `alwinhopf/dssatutils`. `gh` NOT installed:
   `brew install gh && gh auth login && gh repo create alwinhopf/dssatutils --private --source=. --push`
   OR create empty private repo in UI then `git remote add origin ... && git push -u origin main`.
4. `git tag v0.1.0 && git push --tags`.
5. Phase 2 (ML Phenology repo): delete its duplicated r_scripts/weather_*.R &
   soil_*.R; add guarded remotes::install_github(".../dssatutils@v0.1.0")+library;
   replace any source() of those; `renv::init(); renv::snapshot()`.
6. Phase 3 (Gridded repo): replace the 10 weather/soil `source()` lines in
   dssat_main_pipeline.R with `library(dssatutils)` (KEEP 2 landcover source()
   lines); in dssat_main_pipeline.py replace the 10 `from <mod> import process_*`
   with `from dssatutils import (...)` AND fix the 2 lazy imports
   (`import dssatutils.soil_soilgrids_online as _sg_mod`,
   `import dssatutils.weather_nasapower_chirps as _wc`); delete the moved
   weather/soil files from r_scripts/ & python_scripts/ (KEEP landcover_*);
   add the pip pin to requirements.txt + install line to setup_renv.R.
   NOTE: Gridded pipeline sets `USE_REST_API` in global env before calling
   soilgrids_online — that still works with the packaged function.
7. Finish code review: weather_nasapower_chirps.R, soil_hwsd.R, all Python.
8. Optional polish: convert bare library() in R files to @importFrom; unify
   TAV/AMP into one shared helper (R + Python); make soilgrids mode an arg.

## ⚠️ Session status / WARNING for whoever resumes
Implementation was NOT executed. The tool layer (Read AND Bash) returned corrupted/
truncated results in this session (the harness itself flagged results as unreliable;
e.g. `wc -l` duplicated a filename with two different counts; Read line numbers
jumped non-monotonically). Do the real work in a FRESH session.

**FALSE ALARMS — RESOLVED.** Earlier suspected soil "bugs" were display artifacts
from the corruption, now DISPROVEN by two consistent channels (clean Read + Bash
grep). The soil R files are fine:
- `soil_ssurgo.R`: correctly calls `calculate_soil_properties()` (defined L35,
  called L210). NO undefined function. SSURGO SQL is REAL (full queries L167,
  L196-200). Smart-resume + Saxton-Rawls physics all present.
- `soil_soilgrids_online.R`: `format_dssat_sol_file()` (L78) and
  `fetch_soilgrids_vrt()` (L242) are FULLY implemented; REST mode writes .SOL
  files in a tryCatch loop with error logging (L380-394). Not stubs.

**Code review status:** cleanly reviewed = all weather R EXCEPT
weather_nasapower_chirps.R, plus soil_soilgrids.R, soil_ssurgo.R,
soil_soilgrids_online.R. STILL TO REVIEW in a clean session:
weather_nasapower_chirps.R, soil_hwsd.R, and ALL Python modules.

Confirmed weather findings stand (see weather section above): GridMET RH2M/TDEW
proxies, Open-Meteo -99 TDEW/RH, NASA-POWER leftover "FIX" comments, inconsistent
TAV/AMP (gridmet uses DSSAT::calc_TAV/AMP; others hand-roll), bare library() calls
to convert to Imports.

Verified-clean soil findings (minor, optional):
- `soil_soilgrids_online.R` hardcodes `USE_REST_API <- TRUE` at top and reads it
  via `exists()`. As a package this becomes a package-level binding; the pipeline
  overrides by assigning `USE_REST_API` in the global env before the call. Works,
  but a cleaner API would make it a function argument (`mode = "REST"|"VRT"`).
- SSURGO `process_soils_ssurgo` only parallelizes on Windows (`makeCluster`);
  on macOS/Linux `cl <- n_cores` is passed to `pblapply`, which treats an integer
  as fork-based cores — fine, just asymmetric. Worth a comment.
