# DSSAT Crop-Modeling Workspace — Architecture

> **Scope.** This document describes how the sibling DSSAT repositories in this
> workspace fit together, the role of each, the data they exchange, and a prioritized
> set of suggested structural improvements.
>
> **Last verified:** 2026-06-09 · **Vantage point:** `DSSAT_Gridded_Run_Tutorial`
> (the canonical engine).

## Contents

1. [The stack at a glance](#1-the-stack-at-a-glance)
2. [Glossary — DSSAT artifacts](#2-glossary--dssat-artifacts)
3. [Layer 1 — Foundation](#3-layer-1--foundation-read-only-dependencies)
4. [Layer 2 — The engine](#4-layer-2--the-engine)
5. [Layer 3 — Applications](#5-layer-3--applications)
6. [Data flow](#6-data-flow)
7. [Dependency matrix](#7-dependency-matrix)
8. [Suggested improvements](#8-suggested-improvements)

---

## 1. The stack at a glance

Each folder is an independent git repository. They layer into three tiers: a shared
**foundation** (compiled model + data-download library), a reusable **gridded-simulation
engine**, and **specialized applications** that consume them.

```
┌─────────────────────────────────────────────────────────────────┐
│  APPLICATIONS (specialized studies)                              │
│  Bioenergy Input Comparison · dssat_lca_tea · DSSAT_Calibration  │
│  ML Phenology Prediction · SubField MILP · DSSAT_acceleration   │
│  (DSSAT_LAI_Assimilation — empty placeholder)                   │
└───────────────┬─────────────────────────────────────────────────┘
                │ fork / reuse
┌───────────────▼─────────────────────────────────────────────────┐
│  ENGINE:  DSSAT_Gridded_Run_Tutorial                             │
│  dssat_main_pipeline.{R,py} — points→weather+soil→FileX→run→parse│
│  + HPC MPI runner (SLURM)                                        │
└───────────────┬─────────────────────────────────────────────────┘
                │ imports (pinned @vX.Y.Z)        │ invokes binary
┌───────────────▼──────────────────┐  ┌───────────▼─────────────────┐
│  dssatutils  (shared library)    │  │  DSSAT48  (model install)    │
│  R + Python, weather/soil → .WTH │  │  dscsm048 Fortran binary,    │
│  /.SOL  (private GitHub package) │  │  crop genotype files, .CDE   │
└──────────────────────────────────┘  └──────────────────────────────┘
```

---

## 2. Glossary — DSSAT artifacts

Cross-referenced throughout this document so non-DSSAT readers can follow the data flow.

| Artifact | Extension | Meaning |
|---|---|---|
| Weather file | `.WTH` | Daily weather for one point (one file per point-year set). |
| Soil profile | `.SOL` | Layered soil properties for one point or `mukey`. |
| FileX (experiment) | `.??X` (e.g. `.WHX`, `.MZX`, `.SQX`) | Experiment definition: fields, planting, management, treatments. |
| Cultivar | `.CUL` | Genotype coefficients for a named cultivar. |
| Ecotype | `.ECO` | Coefficients shared by a group of cultivars. |
| Species | `.SPE` | Crop-wide physiological constants. |
| Code definitions | `.CDE` | Controlled vocabularies the model reads. |
| Output | `Summary.OUT`, `*.csv` | Per-run results; pipelines typically parse only the summary. |
| HRU | — | Homogeneous Response Unit: a soil+climate(+topo) class simulated once, mapped to many pixels. |

---

## 3. Layer 1 — Foundation (read-only dependencies)

### DSSAT48
A compiled installation of **DSSAT-CSM v4.8**, not project source code:

- The serial, single-point Fortran binary (`dscsm048` / `DSCSM048.EXE`, `DSSATDebug.exe`).
- ~50 per-crop coefficient folders (Wheat, Soybean, Carinata, Peanut, Cotton, …).
- Genotype files (`.CUL`, `.ECO`, `.SPE`), `.CDE` files, `MODEL.ERR`, and the install
  registry `DSSATPRO.V48`.

Every simulation pipeline ultimately spawns this binary. Treated as **immutable** — no
Fortran source is present to recompile.

### dssatutils
A **dual-language shared package** (`github.com/alwinhopf/dssatutils`, private) with
parallel **R** (`R/`) and **Python** (`python/dssatutils/`) implementations of the same
API. Its sole responsibility: fetch public weather/soil data and write DSSAT-format
`.WTH` / `.SOL` files for a set of grid points.

| Domain | Functions (identical names in R and Python) |
|---|---|
| Weather | `process_weather_daymet`, `process_weather_gridmet`, `process_weather_nasapower`, `process_weather_openmeteo`, `process_weather_agera5`, `process_weather_nasapower_chirps` |
| Soil | `process_soils_ssurgo`, `process_soils_soilgrids`, `process_soils_soilgrids_online`, `process_soils_hwsd` |

Coverage: Daymet = North America; GridMET/SSURGO = USA; NASA POWER / Open-Meteo /
AgERA5 / SoilGrids / HWSD2 = global. AgERA5 requires a Copernicus CDS key; CHIRPS fuses
NASA POWER with high-resolution rainfall (50S–50N).

- Function names are identical across R and Python by design.
- Versioned with git tags; **consumers pin to a tag** (`@v0.1.0`), never `main`, so
  upstream changes never silently break a pipeline until the pin is deliberately bumped.
- Product of the extraction migration that de-duplicated download code formerly copied
  between the Gridded Run Tutorial and ML Phenology repos (see `SHARED_UTILS_MIGRATION.md`).

---

## 4. Layer 2 — The engine

### DSSAT_Gridded_Run_Tutorial
The "source of truth" gridded workflow and the most complete repo (renv-locked R, conda
environment, `tests/`, `.github/`). The core is `dssat_main_pipeline.{R,py}`
(~1,900 / ~1,400 lines), a documented 3-step pipeline:

1. **STEP 1 — Soils:** acquire one `.SOL` per grid point (via `dssatutils`).
2. **STEP 2 — Weather:** acquire one `.WTH` per grid point.
3. **STEP 3 — Run:** build a DSSAT-ready folder per point, spawn `dscsm048`, parse outputs.

Key capabilities:

- **Spatial-domain modes:** regular grid from a boundary polygon, a user shapefile, or
  cropland-only points derived from CDL/NLCD rasters (`r_scripts/`, `python_scripts/`).
- **Sensitivity sweep:** the weather × soil × resolution sweep driver now lives in the
  `Bioenergy_Model_Input_Comparison` repo (`run_carinata_sweep.{R,py}`) and references this
  engine via `ENGINE_DIR`.
- **Scale-out:** `hpc/dssat_mpi_runner.py` (SLURM + MPI) with dynamic work-stealing,
  streaming per-rank writes, and node-local scratch staging.
- **Caching:** `*_netcdf_cache/`, `soil_cache/`, `weather_cache/` hold acquired inputs.

---

## 5. Layer 3 — Applications

| Repo | Language | What it adds on top of the engine |
|---|---|---|
| **Bioenergy_Model_Input_Comparison** | Python + R | Compares how weather × soil data-source choices (DAYMET / NASA_POWER / GridMET × SSURGO / SoilGrids) change carinata outputs. Outputs per-combination run dirs and comparison heatmaps. |
| **dssat_lca_tea** | Python **and** R (full parity) | Numbered **LCA/TEA pipeline** (`pipeline_01`…`07`, `pipeline_99_run_master`) on DSSAT yields: ISO-14044 LCA + bottom-up HEFA techno-economic analysis for camelina/carinata SAF. `test_excel_parity` validates against client Excel workbooks; see `PIPELINE_REVIEW.md`. |
| **DSSAT_Calibration** | R (`dssatcal`) | AgMIP-protocol genotype calibration: Morris→Sobol screening → stepwise AICc/BIC selection → optional DREAM-zs Bayesian UQ. CroptimizR-compatible wrapper (`wrapper.R`); `.CUL` free, `.ECO` gated, `.SPE` blocked by default. |
| **DSSAT_ML_Phenology_Prediction** | R (`torch`, xgboost, lightgbm) | **22-model hybrid phenology pipeline** (`pipeline/00`…`13`) — DSSAT physics, SMC-MCMC assimilation, GBMs, sequence DL (LSTM/CNN/Transformer/PINN), Worrall model-guided LSTMs. `features.R` is the single feature source; `MODELS_AND_FEATURES.md` is authoritative. |
| **DSSAT_SubField_MILP_Analysis** | Python (+R companion) | Subfield (50 m) cover-crop bioenergy modeling + **MILP optimization** for spatial placement. Two tracks (single-season vs corn-soy rotation). `dssat_main_pipeline.py` is a thin wrapper that imports `dssatengine` (no local engine defs); `README_bioenergy.md` documents FileX/soil column-alignment fixes. |
| **DSSAT_acceleration** | Python | **New standalone** performance pipeline (30 m / southern-US / 30 yr / 10 mgmt). Treats engine, `dssatutils`, `DSSAT48` as read-only; replaces the per-point loop with HRU dedup + output suppression + RAM-disk batching. Spec in `ACCELERATION_PLAN.md`. |
| **DSSAT_LAI_Assimilation** | — | **Empty** placeholder (likely planned LAI remote-sensing assimilation). |

---

## 6. Data flow

```
shapefile / boundary / CDL-NLCD raster
            │  (gridpoints generation)
            ▼
       grid points  ──►  dssatutils.process_weather_*  ──►  *.WTH
            │       └─►  dssatutils.process_soils_*     ──►  *.SOL
            ▼
   FileX template (.??X) + genotype (.CUL/.ECO/.SPE)
            │  (STEP 3: build per-point run folder)
            ▼
       DSSAT48/dscsm048   ──►  Summary.OUT + ~30 output files
            │  (parser keeps the summary)
            ▼
   merged results.csv  ──►  application layer
                            (maps · LCA/TEA · calibration · ML · MILP)
```

---

## 7. Dependency matrix

Engine logic now lives in the shared **`dssatengine`** package (R + Python, v0.1.0); the
"Engine" column records how each repo obtains it — import, `ENGINE_DIR` reference, or not
at all. **No repo carries a hand-copied engine fork anymore.**

| Repo | DSSAT48 | dssatutils | Engine (`dssatengine`) |
|---|:--:|:--:|:--:|
| DSSAT_Gridded_Run_Tutorial | ✓ | ✓ | imports (canonical wrapper) |
| Bioenergy_Model_Input_Comparison | ✓ | ✓ | via `ENGINE_DIR` → Gridded |
| dssat_lca_tea | ✓ (via yields) | — (reads DSSAT output CSVs) | — (consumes CSVs) |
| DSSAT_Calibration | ✓ | — | — (own CroptimizR wrapper) |
| DSSAT_ML_Phenology_Prediction | ✓ | ✓ | — (own pipeline) |
| DSSAT_SubField_MILP_Analysis | ✓ | ✓ | imports |
| DSSAT_acceleration | ✓ | ✓ (read-only) | imports (transitively via SubField) |
| DSSAT_LAI_Assimilation | — | — | — (empty) |

---

## 8. Suggested improvements

Prioritized by leverage. **Status (2026-06):** the highest-priority item — #1, extract
the engine — is **done**, so the original *silent drift between hand-copied forks* risk
is largely closed; #4 (naming) and #5 (workspace index) are also done. Remaining items
stand.

| # | Improvement | Problem it solves | Effort |
|---|---|---|---|
| 1 | ✅ **Done — engine extracted into the shared `dssatengine` package** (R + Python, v0.1.0), the way `dssatutils` was. `DSSAT_Gridded_Run_Tutorial` and `DSSAT_SubField_MILP_Analysis` now import it (zero local engine defs); `Bioenergy` references it via `ENGINE_DIR`; `DSSAT_acceleration` uses it transitively. No hand-copied forks remain — a fix lands once and reaches every consumer (e.g. the leading-space `DSSBatch` fix). | A bug fixed in one fork never reaching the others — now resolved. | High (done) |
| 2 | **Pin and record the dependency versions each consumer uses** — one line per repo (`dssatutils@v0.1.0`, engine tag, DSSAT48 build). Add a `DEPENDENCIES.md` or a row in each README. | Today you can't tell which repo runs which code without diffing. Makes results reproducible and upgrades deliberate. | Low |
| 3 | **Centralize `dssat_templates/` genotype files.** The same `.CUL/.ECO/.SPE/.CDE` are duplicated across most repos. Ship them from `DSSAT48` (resolved via `DSSATPRO.V48`) or a small shared `dssat-templates` package; keep only project-specific FileX locally. | Genotype edits (e.g. the cereal-rye `TKFH = -25°C` ecotype) must currently be hand-synced. | Medium |
| 4 | ✅ **Done — repo naming standardized** to `snake_case`; the spaced folder names (`DSSAT Gridded Run Tutorial`, …) were renamed with underscores and all references updated. | A whole class of shell-glob / `cd` quoting bugs — now removed. | Low–Medium (done) |
| 5 | ✅ **Done — workspace index README** exists (`dssatengine/README.md`), with a repo table (tier / language / purpose). | New contributors no longer have to open each folder to learn what it is. | Low (done) |
| 6 | **Decide an R↔Python parity policy.** `dssatutils`, `dssat_lca_tea`, and the engine all maintain mirrored implementations by hand. Either (a) declare one language canonical and generate/wrap the other, or (b) add cross-language parity tests (camelina already has `test_excel_parity` — extend the idea). | Halves maintenance and prevents the two implementations from quietly diverging. | Medium |
| 7 | **Repository hygiene.** Commit-ignore and purge transient state: `*_cache/`, `dssat_runs/`, `.RData`, `.Rhistory`, `.DS_Store`, and leftover `.migration_backup_*/` dirs. Verify each repo's `.gitignore` covers them. | Shrinks repos, removes machine-specific noise from diffs, prevents accidental cache commits. | Low |
| 8 | **Scaffold or remove `DSSAT_LAI_Assimilation`.** It is an empty directory. Either add a README + skeleton describing the intended LAI assimilation design, or drop it until work starts. | An empty repo in the stack is ambiguous — planned vs abandoned is unclear. | Low |
| 9 | **Promote shared parsing/output logic.** The DSSAT output parser (keeps `Summary.OUT`, merges to one CSV) is reimplemented per fork. Fold it into the extracted engine package (#1) so the I/O-suppression work in `DSSAT_acceleration` can be shared rather than re-derived. | Stops every consumer from re-solving the same I/O-bound bottleneck. | Medium |
| 10 | **Keep this document live.** Add a short note in each repo README pointing here, and re-verify the dependency matrix whenever a `dssatutils` tag is bumped or a fork is updated. | A snapshot architecture doc decays; a linked, dated one stays trustworthy. | Low |

### Recommended sequence

1. **Quick wins** — #5 and #4 are done; remaining quick wins: #2, #7, #8.
2. ~~Naming standardization (#4)~~ — done.
3. **The structural refactor** (~~#1~~ → #3 → #9 → #6) — #1 (extract the engine) is **done** (`dssatengine`); next centralize templates (#3), promote shared parsing (#9), and settle the parity policy (#6). This mirrors the already-successful `dssatutils` extraction.
