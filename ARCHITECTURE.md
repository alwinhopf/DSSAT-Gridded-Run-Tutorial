# DSSAT Crop-Modeling Workspace — Architecture

> **Scope.** This document describes how the sibling DSSAT repositories in this
> workspace fit together, the role of each, the data they exchange, and a prioritized
> set of suggested structural improvements.
>
> **Last verified:** 2026-07-10 · **Vantage point:** `DSSAT_Gridded_Run_Tutorial`
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
│  Bioenergy Input Comparison · dssat_lca_tea · dssatcalibrator    │
│  ML Phenology Prediction · SubField MILP · pythia (3rd-party)    │
│  (see the application table in section 5)                        │
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
| Weather | `process_weather_daymet`, `process_weather_gridmet`, `process_weather_nasapower`, `process_weather_openmeteo`, `process_weather_agera5`, `process_weather_era5_land`, `process_weather_nasapower_chirps`, `process_weather_nasapower_chirps_v3`, `process_weather_dwd`, `process_weather_eobs`, `process_weather_xavier`, `process_weather_cmfd`, plus local/raster products such as CHELSA-W5E5, AgMERRA/AgCFSR, SILO, PRISM, MSWX/MSWEP, CRU-JRA, TerraClimate, APHRODITE, ANUSPLIN, TAMSAT, GHCN, PGF, and MERRA-2. |
| Soil | `process_soils_ssurgo`, `process_soils_ssurgo_alderman`, `process_soils_gnatsgo`, `process_soils_isdasoil`, `process_soils_lucas`, `process_soils_polaris`, `process_soils_soilgrids`, `process_soils_soilgrids_online`, `process_soils_hwsd`, plus AgMIP/Han, HiHydroSoil, SLGA, WISE30sec, WoSIS, GSDE, China BNU, FEBR, SLC, ESDB, and OpenLandMap. |
| Credentials | `setup_cds_credentials` and compatibility alias `era5land_set_cds_key` are public in both R and Python. They configure CDS-backed sources from `CDSAPI_KEY` / `CDSAPI_URL`, an existing `~/.cdsapirc`, or an interactive prompt. |

Coverage: Daymet = North America; GridMET/SSURGO/gNATSGO = USA; DWD = Germany;
E-OBS / LUCAS = Europe; Xavier (BR-DWGD) = Brazil; CMFD = China; iSDAsoil = Africa;
NASA POWER / Open-Meteo / AgERA5 / ERA5-Land / SoilGrids / HWSD2 = global. AgERA5,
ERA5-Land, and optional E-OBS CDS mode require a Copernicus CDS token; run
`setup_cds_credentials()` once before using them. CHIRPS fuses NASA POWER with high-resolution
rainfall (50S–50N). gNATSGO is the gap-free 30 m US grid; POLARIS is a 30 m probabilistic
disaggregation of SSURGO (USA); iSDAsoil is 30 m Africa; LUCAS is measured EU topsoil
(extrapolated below 20 cm); DWD estimates SRAD from sunshine duration; E-OBS and Xavier
carry daily global radiation directly; CMFD is aggregated from 3-hourly.

- Function names are identical across R and Python by design.
- Provider defaults that are not study-specific live in the shared package
  `config.yml` and are merged with consumer config files. Current package-level
  defaults include Open-Meteo throttling, CHIRPS v3 settings, SoilGrids
  REST/VRT behavior, and the CDS API URL.
- Public exports are hand-maintained in `dssatutils/NAMESPACE` and
  `python/dssatutils/__init__.py`; keep them in parity and run the export audit
  whenever a source or helper is added.
- Versioned with git tags; **consumers pin to a tag** (`@vX.Y.Z`), never `main`, so
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
| **dssatcalibrator** | Python **and** R (full parity) | General **config-driven** framework to calibrate DSSAT-CSM for a new environment / crop / scenario via Monte-Carlo sampling + Bayesian techniques: spawns many perturbed runs in parallel and fits parameter set(s) to a central observations store (LAI, biomass, grain yield, phenology) across one or many experiments (one `config.yaml` read by both languages; mirrored `run_calibration.{py,R}`). **Now the workspace's single calibration repo:** the former R-only `DSSAT_Calibration` (CroptimizR / AgMIP-stepwise) is superseded by this package's R twin, and in-season LAI assimilation (the former `DSSAT_LAI_Assimilation` scaffold) is a **mode** here (satellite LAI source → coupled recalibration / `nowcast` forecast), not a separate repo. |
| **DSSAT_ML_Phenology_Prediction** | R (`torch`, xgboost, lightgbm) | **22-model hybrid phenology pipeline** (`pipeline/00`…`13`) — DSSAT physics, SMC-MCMC assimilation, GBMs, sequence DL (LSTM/CNN/Transformer/PINN), Worrall model-guided LSTMs. `features.R` is the single feature source; `MODELS_AND_FEATURES.md` is authoritative. |
| **DSSAT_SubField_MILP_Analysis** | Python (+R companion) | Subfield (50 m) cover-crop bioenergy modeling + **MILP optimization** for spatial placement. Two tracks (single-season vs corn-soy rotation). Application orchestrators write engine overlay configs and shell out to this sibling repo's `dssat_main_pipeline.{py,R}` through `DSSAT_ENGINE_DIR`; `README_bioenergy.md` documents FileX/soil column-alignment fixes. |
| **pythia** | Python | Vendored **third-party** DSSAT-Pythia tool (independent gridded runner). Not part of the layered stack — treated as read-only upstream; see its own install guide. |

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

Engine logic now lives in the shared **`dssatengine`** package (R + Python, 0.4.0); the
"Engine" column records how each repo obtains it — import, `ENGINE_DIR` reference, or not
at all. **No repo carries a hand-copied engine fork anymore.**

| Repo | DSSAT48 | dssatutils | Engine (`dssatengine`) |
|---|:--:|:--:|:--:|
| DSSAT_Gridded_Run_Tutorial | ✓ | ✓ | imports (canonical wrapper) |
| Bioenergy_Model_Input_Comparison | ✓ | ✓ | via `ENGINE_DIR` → Gridded |
| dssat_lca_tea | ✓ (via yields) | — (reads DSSAT output CSVs) | — (consumes CSVs) |
| dssatcalibrator | ✓ | optional (`[acquire]`) | optional executor (`[shared]`) |
| DSSAT_ML_Phenology_Prediction | ✓ | ✓ | — (own pipeline) |
| DSSAT_SubField_MILP_Analysis | ✓ | ✓ | imports |
| pythia | — | — | — (independent third-party tool) |

---

## 8. Suggested improvements

Prioritized by leverage. **Status (2026-06):** the highest-priority item — #1, extract
the engine — is **done**, so the original *silent drift between hand-copied forks* risk
is largely closed; #4 (naming) and #5 (workspace index) are also done. Remaining items
stand.

| # | Improvement | Problem it solves | Effort |
|---|---|---|---|
| 1 | ✅ **Done — engine extracted into the shared `dssatengine` package** (R + Python, 0.4.0), the way `dssatutils` was. `DSSAT_Gridded_Run_Tutorial` and `DSSAT_SubField_MILP_Analysis` now import it (zero local engine defs); `Bioenergy` references it via `ENGINE_DIR`; `dssatcalibrator` can use its public executor behind `execution.backend: dssatengine`. No hand-copied forks remain — a fix lands once and reaches every consumer (e.g. the leading-space `DSSBatch` fix, explicit `treatment_list` support, and fail-loud DSSAT execution). | A bug fixed in one fork never reaching the others — now resolved. | High (done) |
| 2 | ✅ **Done — dependency versions are recorded** in `dssatengine/DEPENDENCIES.md` and committed manifests. Consumers pin the immutable workspace baselines `dssatutils@e9c859fa1d915623df23e2eb13084cb085dbfe3e` and `dssatengine@31085c7eac1628db949e3ad9fdb16947a65d0834`; `USE_LOCAL_SHARED_PACKAGES=1` opts into sibling checkouts for development. | Makes results reproducible and upgrades deliberate. | Low (done) |
| 3 | **Centralize `dssat_templates/` genotype files.** The same `.CUL/.ECO/.SPE/.CDE` are duplicated across most repos. Ship them from `DSSAT48` (resolved via `DSSATPRO.V48`) or a small shared `dssat-templates` package; keep only project-specific FileX locally. | Genotype edits (e.g. the cereal-rye `TKFH = -25°C` ecotype) must currently be hand-synced. | Medium |
| 4 | ✅ **Done — spaced and hyphenated local folder names removed**; current local repo names use underscores while preserving existing capitalization. New repos should use `snake_case`. | A whole class of shell-glob / `cd` quoting bugs — now removed. | Low–Medium (done) |
| 5 | ✅ **Done — workspace index README** exists at the workspace root (`README.md`), with a repo table (tier / language / purpose). | New contributors no longer have to open each folder to learn what it is. | Low (done) |
| 6 | **Decide an R↔Python parity policy.** `dssatutils`, `dssat_lca_tea`, and the engine all maintain mirrored implementations by hand. Either (a) declare one language canonical and generate/wrap the other, or (b) add cross-language parity tests (camelina already has `test_excel_parity` — extend the idea). | Halves maintenance and prevents the two implementations from quietly diverging. | Medium |
| 7 | **Repository hygiene.** Commit-ignore and purge transient state: `*_cache/`, `dssat_runs/`, `.RData`, `.Rhistory`, `.DS_Store`, and leftover `.migration_backup_*/` dirs. Verify each repo's `.gitignore` covers them. | Shrinks repos, removes machine-specific noise from diffs, prevents accidental cache commits. | Low |
| 8 | ✅ **Done — the empty `DSSAT_LAI_Assimilation` placeholder is no longer in the workspace checkout.** | An empty repo in the stack was ambiguous — now removed. | Low (done) |
| 9 | **Promote shared parsing/output logic.** The DSSAT output parser (keeps `Summary.OUT`, merges to one CSV) is reimplemented per fork. Fold it into the extracted engine package (#1) so the I/O-suppression work explored for acceleration can be shared rather than re-derived. | Stops every consumer from re-solving the same I/O-bound bottleneck. | Medium |
| 10 | **Keep this document live.** Add a short note in each repo README pointing here, and re-verify the dependency matrix whenever a `dssatutils` tag is bumped or a fork is updated. | A snapshot architecture doc decays; a linked, dated one stays trustworthy. | Low |

### Recommended sequence

1. **Quick wins** — #5 and #4 are done; remaining quick wins: #2, #7, #8.
2. ~~Naming standardization (#4)~~ — done.
3. **The structural refactor** (~~#1~~ → #3 → #9 → #6) — #1 (extract the engine) is **done** (`dssatengine`); next centralize templates (#3), promote shared parsing (#9), and settle the parity policy (#6). This mirrors the already-successful `dssatutils` extraction.
