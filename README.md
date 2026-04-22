# Spatial Gridded Crop Modeling with DSSAT

This repository provides an end-to-end, beginner-friendly workflow for **spatial (gridded / point-based) crop modeling in DSSAT**. The pipeline downloads weather and soil data for a set of geographic points, builds a DSSAT-ready input folder for every point, runs the simulations locally or prepares them for HPC, and maps the results.

**Start here:** open `dssat_main_pipeline.R` in **RStudio** and click **Source** to run a small spatial demo locally (no command line needed). Once that works you can scale up on HPC/cloud using **SLURM + MPI**.

> Scripts and examples are a work in progress; several aspects are not yet optimised or finalised.

---

## Table of Contents

- [Status, limitations, and future work](#status-limitations-and-future-work)
- [What is DSSAT and why spatial modeling?](#what-is-dssat-and-why-spatial-modeling)
- [Prerequisites](#prerequisites)
- [Quick start: run the demo in RStudio](#quick-start-run-the-demo-in-rstudio)
- [Choosing your spatial input: shapefile types explained](#choosing-your-spatial-input-shapefile-types-explained)
  - [Option 1 — US state/county boundary (TIGER/Line)](#option-1--us-statecounty-boundary-tigerline)
  - [Option 2 — Cropland Data Layer (USDA CDL)](#option-2--cropland-data-layer-usda-cdl)
  - [Option 3 — NLCD land-cover raster](#option-3--nlcd-land-cover-raster)
  - [Option 4 — Your own shapefile](#option-4--your-own-shapefile)
  - [Which option should I use?](#which-option-should-i-use)
- [Defining your spatial domain (pipeline modes)](#defining-your-spatial-domain-pipeline-modes)
  - [Mode A — Regular grid from a boundary polygon](#mode-a--regular-grid-from-a-boundary-polygon-default-demo)
  - [Mode B — Your own shapefile](#mode-b--your-own-shapefile)
  - [Mode C — Cropland-only points (CDL / NLCD)](#mode-c--cropland-only-points-cdl--nlcd)
- [Soil data sources](#soil-data-sources)
- [Weather data sources](#weather-data-sources)
- [Repository layout](#repository-layout)
- [Background: DSSAT inputs and run-folder anatomy](#background-dssat-inputs-and-run-folder-anatomy)
  - [Regular vs crop-sequence runs](#regular-vs-crop-sequence-runs)
  - [Simulation Controls explained](#simulation-controls-explained)
- [Concepts and conventions](#concepts-and-conventions)
- [Step 0 — Prepare DSSAT templates](#step-0--prepare-dssat-templates)
- [Step 1 — Configure and run the R pipeline](#step-1--configure-and-run-the-r-pipeline)
- [Validate one point folder before scaling](#validate-one-point-folder-before-scaling)
- [Option A — Run locally](#option-a--run-locally)
- [Option B — Run on HPC / cloud (SLURM + MPI)](#option-b--run-on-hpc--cloud-slurm--mpi)
- [Outputs and what they mean](#outputs-and-what-they-mean)
- [Post-processing: join results back to maps](#post-processing-join-results-back-to-maps)
- [How the MPI Python runner works](#how-the-mpi-python-runner-works)
- [Troubleshooting](#troubleshooting)
- [Performance tips](#performance-tips)
- [Reproducibility checklist](#reproducibility-checklist)
- [Publishing checklist](#publishing-checklist)
- [Installing DSSAT](#installing-dssat)
  - [Windows](#windows)
  - [macOS (Apple Silicon M1/M2/M3/M4/M5)](#macos-apple-silicon-m1m2m3m4m5)
  - [Linux](#linux)
- [R package requirements](#r-package-requirements)
- [Advanced: HPC Python environment setup](#advanced-hpc-python-environment-setup)
- [Related tools and gridded crop modeling ecosystems](#related-tools-and-gridded-crop-modeling-ecosystems)
- [References](#references)

---

## Status, limitations, and future work

This repository is a **work in progress**. The core workflow (R folder builder → local or HPC execution → parsed results) is usable today, but scaling, portability, and robustness are still being improved.

### What works today

- **Folder generation in R:** builds one DSSAT-ready folder per grid point (weather + soil + FileX/SQX template).
- **Local runs:** useful for debugging templates and small studies.
- **HPC runs with SLURM + MPI:** distributes points across ranks, runs DSSAT, parses CSV outputs, and merges into a single results file.
- **Output parsing:** merges yield/biomass/management/emissions metrics plus seasonal-average stress indicators.

### Known limitations

- **I/O-heavy at scale:** DSSAT creates many small files per run. On shared filesystems (Lustre/GPFS/NFS), performance can become dominated by metadata overhead once you run many hundreds to thousands of concurrent ranks.
- **Static work allocation:** MPI ranks are assigned points via simple slicing. Some points may run slower, causing load imbalance.
- **Memory accumulation:** each rank accumulates results in memory before writing `temp_results_rank_<rank>.csv`. For very large point counts this can become memory-heavy.
- **Single-process merge:** rank 0 merges all rank files; this is a bottleneck for huge outputs.
- **Assumes CSV outputs:** the parser expects DSSAT CSV output files (`summary.csv`, `soilwat.csv`, etc.). If your DSSAT setup produces only `*.OUT` text files, you must enable CSV outputs or adapt the parser.

### Roadmap

- Streaming rank output to disk (write rows as you go).
- Dynamic scheduling (work queue) for better load balancing.
- Node-local scratch execution to reduce shared filesystem pressure.
- Checkpoint/restart: skip completed points, log failures to `failures.csv`.
- Parquet outputs for faster downstream analytics.
- `environment.yml` / Apptainer container for reproducible HPC installs.

---

## What is DSSAT and why spatial modeling?

**DSSAT (Decision Support System for Agrotechnology Transfer)** is a crop modeling framework that simulates crop growth, development, and yield as a function of soil–plant–atmosphere dynamics. It uses daily weather, soil physical and chemical properties, crop genetics, and management as inputs, and supports many crops (maize, wheat, soybean, sorghum, rice, and many more).

A **crop model** is a mathematical representation of how crops respond to weather (temperature, rainfall, solar radiation, humidity, wind), soil conditions (water holding, texture, nutrients, organic matter), and management decisions (planting date, cultivar choice, irrigation, fertilisation, tillage, harvest). Why model spatially?

- To scale results from plot → farm → region.
- To explore "what-if" scenarios (irrigation, fertiliser rates, cultivar selection) across many locations simultaneously.
- To examine long-term effects (rotations, cover crops, soil C/N changes) at landscape scale.
- To identify spatial heterogeneity in yield gaps, water stress, or nutrient limitations.

> Models are only as good as their **inputs** and **calibration**. Always validate against field observations before applying at scale.

Official overview: https://dssat.net/

---

## Prerequisites

### DSSAT (required)

You need a DSSAT installation with a command-line executable. See [Installing DSSAT](#installing-dssat) for platform-specific instructions.

Typical executable names:
- `dscsm048` on Linux/macOS
- `DSCSM048.EXE` on Windows

### R (required for folder building and local runs)

- R 4.x recommended
- RStudio (recommended for beginners)
- Required packages: `sf`, `dplyr`, `tidyverse`, plus soil/weather-specific packages — see [R package requirements](#r-package-requirements).

### Python + MPI (required for HPC execution; optional locally)

- Python 3.9+
- `mpi4py`
- Standard library modules used by the runner: `argparse`, `csv`, `glob`, `os`, `re`, `shutil`, `subprocess`

### SLURM (optional)

If you have a SLURM-managed HPC cluster, you can run the MPI runner with `sbatch`. Any other scheduler (PBS, LSF, SGE) works too — translate `sbatch`/`srun` to your scheduler's equivalents.

---

## Quick start: run the demo in RStudio

The easiest path: open the R pipeline in RStudio, set the DSSAT path once, and click **Source**.

### 1 — Open the repo in RStudio

Open the repository as an RStudio Project (recommended), or set your working directory to the repo root, then open `dssat_main_pipeline.R`.

### 2 — Tell the script where DSSAT lives

Run one of these **once per R session** in the RStudio Console (no code edits required):

```r
# Linux / macOS
Sys.setenv(DSSAT_EXE = "/full/path/to/dscsm048")

# Windows
Sys.setenv(DSSAT_EXE = "C:/DSSAT48/DSCSM048.exe")
```

### 3 — Download the demo boundary shapefile (one-time)

The default demo uses US Census TIGER/Line state boundaries:

```bash
mkdir -p shapefile
curl -L -o shapefile/tl_2024_us_state.zip \
  "https://www2.census.gov/geo/tiger/TIGER2024/STATE/tl_2024_us_state.zip"
unzip -o shapefile/tl_2024_us_state.zip -d shapefile
```

### 4 — Place a DSSAT template file

Copy a known-good experiment file (e.g., `UFGA8201.MZX` from your DSSAT `Maize/` folder) into `dssat_templates/`, then set `TEMPLATE_FILE_NAME` in Section 0 of the pipeline script. `UFGA8201.MZX` is a maize file bundled with DSSAT for testing. Any valid DSSAT experiment file works (`.MZX`, `.WHX`, `.SBX`, `.SQX`, etc.) — match the extension to your `CROP_EXTENSION` and `RUN_MODE` settings.

### 5 — Source the pipeline

Click **Source** in RStudio. With default demo settings (Montana at 50 km spacing, 3 years of Daymet weather, SSURGO soil, local DSSAT run), the script will:

- Generate ~60 grid points across Montana
- Download soil and weather data for each point
- Build one DSSAT run folder per point
- Run DSSAT locally
- Parse and merge results
- Write a merged results CSV and a yield map under `results/`

First run can take 10–30 minutes depending on your internet connection and machine speed.

### 6 — Confirm you got results

After the run completes, look for:

- **Merged results:** `results/<RUN_NAME>_results.csv`
- **Maps:** `results/<RUN_NAME>_yield_map_treatment<k>.png`
- **Per-point folders:** `dssat_runs/<RUN_NAME>/<POINT_ID>/` (includes `SOIL.SOL`, `*.WTH`, the per-point FileX/SQX, and DSSAT outputs)

### 7 — If the first run fails

- DSSAT path not found → confirm `Sys.setenv(DSSAT_EXE = ...)` is set correctly.
- Missing boundary shapefile → download the TIGER/Line file (Step 3 above) or point the script to your own boundary.
- Template/run-mode mismatch (e.g., `RUN_MODE = "sequence"` but template is `*.MZX`) → see [Background: DSSAT inputs](#background-dssat-inputs-and-run-folder-anatomy) and update `RUN_MODE` or `TEMPLATE_FILE_NAME` in Section 0.

Once the demo works, customise settings in **SECTION 0: MASTER CONFIGURATION** of `dssat_main_pipeline.R`.

---

## Choosing your spatial input: shapefile types explained

Before running the pipeline, you need to decide **which locations to simulate**. This choice determines how densely and how accurately your simulations represent actual agricultural conditions. Below is a plain-language explanation of the most common options.

---

### Option 1 — US state/county boundary (TIGER/Line)

**What it is:** The US Census Bureau TIGER/Line dataset provides administrative boundary polygons for every US state, county, congressional district, and more. The `tl_2024_us_state.shp` file used in the demo contains one polygon per state.

**What the pipeline does with it:** The pipeline reads the state polygon(s) you select (e.g., Montana), places a regular grid of points at your chosen spacing inside that polygon, and simulates at every grid point — regardless of what land cover is actually there. Some points may fall on forests, lakes, or cities.

**Best for:**
- Quick demos and regional studies where you want comprehensive spatial coverage.
- Climate impact assessments across all land types.
- Any area where you do not have a more specific point set.

**Download:**
```bash
mkdir -p shapefile
curl -L -o shapefile/tl_2024_us_state.zip \
  "https://www2.census.gov/geo/tiger/TIGER2024/STATE/tl_2024_us_state.zip"
unzip -o shapefile/tl_2024_us_state.zip -d shapefile
```

For counties, substitute `STATE` → `COUNTY` in the URL:
```bash
curl -L -o shapefile/tl_2024_us_county.zip \
  "https://www2.census.gov/geo/tiger/TIGER2024/COUNTY/tl_2024_us_county.zip"
unzip -o shapefile/tl_2024_us_county.zip -d shapefile
```

TIGER/Line products are free to use and not copyrighted. If you publish results, acknowledge the US Census Bureau as the data source.

---

### Option 2 — Cropland Data Layer (USDA CDL)

**What it is:** The USDA National Agricultural Statistics Service (NASS) Cropland Data Layer (CDL) is an annual, crop-specific land-cover raster for the contiguous United States at ~30 m resolution. Each pixel is classified as a specific crop type (corn, soybeans, winter wheat, cotton, etc.) or a non-agricultural class (forest, water, urban, etc.). CDL is updated every year and is the gold standard for identifying where specific crops are grown in the US.

**Why use it for DSSAT:** Running DSSAT on every grid cell in a state wastes compute resources on forests, cities, and water bodies. By masking to CDL cropland pixels only, you ensure that every simulation point corresponds to actual (or likely) agricultural land. You can also filter to a **specific crop** — for example, simulate only on corn pixels if you are calibrating or validating a maize template.

**Key CDL class codes:**
```r
cropland_values <- c(1)        # corn only
cropland_values <- c(5)        # soybeans only
cropland_values <- c(1, 5)     # corn + soybeans
cropland_values <- c(2)        # cotton
cropland_values <- c(24)       # winter wheat
cropland_values <- c(1, 5, 24) # corn + soybeans + winter wheat
```

Full CDL class code list: https://www.nass.usda.gov/Research_and_Science/Cropland/sarsfaqs2.php

**Data sources:**
- Interactive map and state-level downloads: https://croplandcros.scinet.usda.gov/
- National annual GeoTIFF download by state and year: https://www.nass.usda.gov/Research_and_Science/Cropland/Release/index.php
- Google Earth Engine (scripted access): https://developers.google.com/earth-engine/datasets/catalog/USDA_NASS_CDL

**Workflow:** CDL is a raster — the pipeline cannot use it directly as a boundary polygon. You must first run the two landcover helper scripts to convert it to a point shapefile. See [Mode C](#mode-c--cropland-only-points-cdl--nlcd) below.

---

### Option 3 — NLCD land-cover raster

**What it is:** The National Land Cover Database (NLCD) is a 30 m land-cover classification for the contiguous US, updated every 2–3 years. NLCD class `82` ("Cultivated Crops") is a broad cropland mask that captures all agricultural land without distinguishing specific crop types.

**Compared to CDL:**

| | CDL | NLCD |
|--|-----|------|
| Crop specificity | Per-crop codes (corn, soy, wheat, etc.) | Class 82 = all cultivated crops |
| Update frequency | Annual | Every 2–3 years |
| Best use | Run on a specific crop's footprint | General cropland mask for any crop |

**Download:** https://www.mrlc.gov/data (Annual NLCD GeoTIFF)
**Class legend:** https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description

**Workflow:** Same as CDL — run the two landcover helper scripts to convert to a point shapefile, then use Mode B/C.

---

### Option 4 — Your own shapefile

If you already have a set of locations — farm field boundaries, experimental station coordinates, watershed polygons, research site networks, or any other vector dataset — supply it directly.

- Polygons are automatically converted to centroids by the pipeline.
- Any coordinate system is accepted; the pipeline reprojects to WGS84.
- No specific attribute schema is required — the pipeline only needs geometry.

Useful for precision-agriculture field-scale simulations, research station networks, or any study area for which you have an existing digital boundary.

---

### Which option should I use?

| Goal | Recommended input |
|------|-------------------|
| Quick first test / debugging | Mode A — TIGER/Line state boundary, 50 km spacing |
| Regional study across all land types | Mode A — TIGER/Line, production spacing (5–10 km) |
| Simulate only on cropland (all crop types) | Mode C — NLCD class 82 raster |
| Simulate only on corn or soybean footprints | Mode C — CDL with crop-specific codes |
| Known field boundaries or research sites | Mode B — your own shapefile |
| Non-US region (global simulation) | Mode A with a custom boundary polygon, or Mode B |

---

## Defining your spatial domain (pipeline modes)

All three modes produce the same standardised point shapefile (with `ID` / `LAT` / `LONG` columns) that the rest of the pipeline consumes. You can switch modes without changing anything after Section 0.

The controlling flag in `dssat_main_pipeline.R` Section 0:

```r
USE_EXISTING_POINT_SHAPEFILE <- FALSE   # MODE A — generate a regular grid
USE_EXISTING_POINT_SHAPEFILE <- TRUE    # MODE B / C — supply your own points
```

---

### Mode A — Regular grid from a boundary polygon (default demo)

The pipeline reads a boundary shapefile, filters to the region(s) you want, and places points on a regular grid at `GRID_SPACING_METERS` spacing.

```r
# SECTION 0 of dssat_main_pipeline.R
USE_EXISTING_POINT_SHAPEFILE <- FALSE

BOUNDARY_SHAPEFILE_NAME  <- "tl_2024_us_state.shp"
ENABLE_BOUNDARY_FILTER   <- TRUE
BOUNDARY_FILTER_COLUMN   <- "NAME"
STATE_NAME_FILTER        <- c("Montana")      # one or more state names

GRID_SPACING_METERS      <- 50000            # 50 km → ~60 points (fast demo)
                                              # 5000  → ~6 000 points (production)
```

---

### Mode B — Your own shapefile

Supply any point or polygon shapefile. Polygons are auto-converted to centroids.

```r
USE_EXISTING_POINT_SHAPEFILE  <- TRUE
EXISTING_POINT_SHAPEFILE_PATH <- file.path(MAIN_PROJECT_DIR,
                                            "gridpoints", "my_study_sites.shp")
```

The `load_existing_points()` function handles the rest automatically:
- Polygons / lines → converted to centroids
- Any CRS → re-projected to WGS84 (EPSG:4326)
- ID column → standardised to zero-padded 8-digit strings
- Missing / duplicate IDs → regenerated sequentially

---

### Mode C — Cropland-only points (CDL / NLCD)

This is Mode B with a cropland-masked point shapefile that you build first using two helper scripts.

#### Step 1 — Download a landcover raster

| Dataset | Class code | Coverage | URL |
|---------|-----------|----------|-----|
| NLCD | `82` = Cultivated Crops (all crops) | Contiguous US, ~30 m | https://www.mrlc.gov/data |
| CDL | Crop-specific codes (e.g., `1` = corn) | Contiguous US, ~30 m | https://croplandcros.scinet.usda.gov/ |
| ESA WorldCover | `40` = Cropland | Global, 10 m | https://worldcover2021.esa.int/ |

Place the downloaded `.tif` under `data/landcover/`.

#### Step 2 — Create a binary cropland mask (`r_scripts/landcover_raster.R`)

Edit the USER SETTINGS block at the top of the script:

```r
input_raster    <- file.path("data", "landcover",
                              "Annual_NLCD_LndCov_2024_CU_C1V1.tif")
boundary_vector <- file.path("shapefile", "tl_2024_us_state.shp")
state_names     <- c("Montana", "North Dakota", "South Dakota")

# NLCD class 82 = all cultivated crops (broad)
# CDL: use crop-specific codes, e.g. c(1, 5) for corn + soybeans
cropland_values <- c(82)

output_dir      <- file.path("data", "landcover", "derived")
write_per_state <- TRUE
```

Run it:
```bash
Rscript r_scripts/landcover_raster.R
```

Output: `data/landcover/derived/cropland_mask_<state>.tif` — a binary raster (1 = cropland).

#### Step 3 — Aggregate to grid points (`r_scripts/landcover_raster_to_gridpoints.R`)

This aggregates the high-resolution (~30 m) mask to your target DSSAT grid spacing and writes a point shapefile with a `crop_pct` attribute (fraction of the grid cell covered by cropland).

```r
crop_raster_file    <- file.path("data", "landcover", "derived",
                                  "cropland_mask_montana.tif")
output_dir          <- "gridpoints"
grid_resolution_m   <- 5000   # must match GRID_SPACING_METERS in dssat_main_pipeline.R
cropland_threshold  <- 0      # keep any cell with any cropland
                               # raise to e.g. 0.5 to keep only cells ≥50% cropland
```

Run it:
```bash
Rscript r_scripts/landcover_raster_to_gridpoints.R
```

Outputs in `gridpoints/`:
- `montana_cropland_5k.shp` — all aggregated cells with a `crop_pct` attribute
- `montana_cropland_5k_above_threshold.shp` — filtered to cells above `cropland_threshold`

#### Step 4 — Feed into Mode B

```r
USE_EXISTING_POINT_SHAPEFILE  <- TRUE
EXISTING_POINT_SHAPEFILE_PATH <- file.path(MAIN_PROJECT_DIR,
                                 "gridpoints", "montana_cropland_5k_above_threshold.shp")
```

---

## Soil data sources

Set `SOIL_SOURCE` in Section 0 of `dssat_main_pipeline.R`.

| Value | Coverage | Method | Notes |
|-------|----------|--------|-------|
| `SSURGO` | United States only | Queries USDA Soil Data Access (SDA) web service per point | Most detailed US data; requires internet; queries one point at a time |
| `SOILGRIDS_10K` | Global | Reads a pre-downloaded master `.SOL` file at 10 km resolution | Fastest for large domains; download the country file once (see below) |
| `SOILGRIDS_ONLINE` | Global | Queries ISRIC SoilGrids 2.0 REST API or VRT/GDAL per point | Flexible; use `USE_REST_API` flag to switch between REST (interactive) and VRT (HPC batch) mode |

### SoilGrids 10K pre-formatted files

Pre-formatted DSSAT-ready `.SOL` files at 10 km resolution, organised by country:

- **Download (Harvard Dataverse):** https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/1PEEY0
- **Citation:** Folberth et al. (2019). *Environmental Modelling & Software* 111:218–228. https://www.sciencedirect.com/science/article/pii/S1364815218313033

Download the country file you need (e.g. `US.SOL`), place it under `SoilGrids/`, and set:
```r
EXTERNAL_SOIL_FILE <- file.path(MAIN_PROJECT_DIR, "SoilGrids", "US.SOL")
```

### SoilGrids online: VRT vs REST API

The `USE_REST_API` flag in `soil_soilgrids_online.R` controls how data is fetched:

```r
USE_REST_API <- TRUE    # JSON REST API — best for interactive / local use;
                         # rate-limited; includes exponential back-off retry (up to 5 attempts)
USE_REST_API <- FALSE   # VRT/GDAL virtual rasters — best for HPC / large batch jobs;
                         # avoids per-point HTTP overhead; benefits from GDAL block caching
```

---

## Weather data sources

Set `WEATHER_SOURCE` in Section 0 of `dssat_main_pipeline.R`.

| Value | Coverage | Approx. resolution | Notes |
|-------|----------|--------------------|-------|
| `DAYMET` | North America | ~1 km daily | Best choice for US simulations; uses the `daymetr` package |
| `NASA_POWER` | Global | ~0.5° daily | Recommended for international runs; outputs real wind speed (`WS2M`) at `WNDHT = 2.0` m |
| `GRIDMET` | Contiguous US | ~4 km daily | High spatial resolution for US; requires `terra`, `ncdf4`, and `httr` |

> **NASA-POWER wind note:** Previous versions of `weather_nasapower.R` wrote `WIND = -99` (missing) and set `WNDHT = -99` in the `.WTH` header. The current version writes the real `WS2M` value and correctly sets `REFHT = 2.0` and `WNDHT = 2.0`, matching the NASA-POWER AG community product specification.

---

## Repository layout

```text
.
├── dssat_main_pipeline.R              # MAIN entrypoint — start here
├── r_scripts/
│   ├── soil_ssurgo.R
│   ├── soil_soilgrids.R
│   ├── soil_soilgrids_online.R
│   ├── weather_daymet.R
│   ├── weather_nasapower.R
│   ├── weather_gridmet.R
│   ├── landcover_raster.R
│   └── landcover_raster_to_gridpoints.R
├── dssat_templates/                   # template FileX/SQX + cultivar/ecotype/species files
├── shapefile/                         # boundary polygons (TIGER/Line or custom)
├── gridpoints/                        # generated or user-supplied point shapefiles
├── data/landcover/                    # landcover rasters and derived masks (Mode C)
├── soil/                              # downloaded/cached soil products
├── weather/                           # downloaded/cached weather products
├── SoilGrids/                         # pre-downloaded master .SOL files (SOILGRIDS_10K)
├── dssat_runs/                        # generated per-point DSSAT run folders
├── results/                           # merged CSVs and yield maps
├── dssat_mpi_runner.py                # Advanced: MPI runner for HPC
└── run_dssat_python.slurm             # Advanced: SLURM submit script
```

> **`.gitignore` note:** the provided `.gitignore` prevents accidentally committing generated run folders (`dssat_runs/`), weather/soil inputs, and merged outputs (`results/`). Review and adjust for your team's data-sharing preferences.

---

## Background: DSSAT inputs and run-folder anatomy

A **DSSAT run folder** is a directory of text files that the DSSAT executable reads from the command line. In this repo, every grid point gets its **own run folder** so the same template scenario is repeated across space.

Each point folder typically needs:

- **A scenario file:** a crop-specific FileX (e.g., `*.MZX`, `*.WHX`, `*.SBX`) for a single-season experiment, or a `*.SQX` for a multi-year cropping sequence.
- **Weather:** `*.WTH` (daily weather) referenced by station ID in the scenario file.
- **Soil:** `SOIL.SOL` containing the soil profile ID referenced by the scenario file.
- **Support files:** cultivar (`*.CUL`), ecotype (`*.ECO`), species (`*.SPE`) files required by your crop module.
- **DSSBatch:** `DSSBatch.V48` tells DSSAT which scenario file to run and which treatments/sequences to execute. This is written automatically by the pipeline.

---

### Regular vs crop-sequence runs

DSSAT supports two fundamentally different run types. Your choice determines the template file extension and the `RUN_MODE` setting in the pipeline.

| Concept | Regular experiment (e.g., `UFGA8201.MZX`) | Crop sequence (e.g., `UFGA7804.SQX`) |
|---|---|---|
| Primary use | Single-season management comparison | Multi-year rotation / long-term soil dynamics |
| Years simulated | Usually `NYERS = 1` | Often `NYERS > 1` (5, 10, 30 years) |
| Treatments | Each row is an independent treatment | One sequence has ordered steps (Bean → Fallow → Soybean) |
| Soil state between seasons | Reset to initial conditions each run | Carried forward (water, nitrate, residue, SOM) |
| Best for | Yield response to N, irrigation, planting date | Rotation impacts, cover crop effects, long-term soil C/N |
| Pipeline `RUN_MODE` setting | `"experiment"` | `"sequence"` |

**Use `experiment`** when each treatment is independent (typical single-season management comparisons).
**Use `sequence`** when you want multi-year carryover (rotation / cover crop / soil C-N dynamics).

> **MPI runner note:** the current `dssat_mpi_runner.py` always looks for a `*.SQX` file and runs DSSAT with the `Q` switch. Treat it as **sequence-mode only** unless you extend it to detect and handle seasonal FileX types.

---

### Simulation Controls explained

The `*SIMULATION CONTROLS` section inside a FileX/SQX is the most important section to understand — it controls which processes DSSAT simulates and how management is applied.

Example from `UFGA8201.MZX`:

```text
*SIMULATION CONTROLS
@N GENERAL     NYERS NREPS START SDATE RSEED SNAME
 1 GE          1     1     S 82056  2150 RAINFED LOW N

@N OPTIONS     WATER NITRO SYMBI PHOSP POTAS DISES CHEM TILL CO2
 1 OP          Y     Y     N     N     N     N     N    Y    M

@N METHODS     WTHER INCON LIGHT EVAPO INFIL PHOTO HYDRO
 1 ME          M     M     E     R     S     R     R

@N MANAGEMENT  PLANT IRRIG FERTI RESID HARVS
 1 MA          R     R     R     N     M
```

What these mean for a beginner:

- `NYERS = 1`: one year simulated (set higher for sequence/rotation runs).
- `WATER = Y`: soil water balance is on.
- `NITRO = Y`: nitrogen balance is on.
- `SYMBI = N`: symbiotic N fixation off (correct for maize; set `Y` for legumes like soybean or bean).
- `TILL = Y`: tillage and residue handling is enabled.
- `CO2 = M`: CO₂ is taken from a measured values file. See https://dssat.net/weather-module/ for CO₂ options.
- `PLANT = R`, `IRRIG = R`, `FERTI = R`: management events are read from the explicit schedule in the experiment file sections.
- `HARVS = M`: harvest is managed/automatic (maturity-driven).

For a crop sequence run with `UFGA7804.SQX` (bean–fallow–soybean rotation), you will see `NYERS = 10` and `SYMBI = Y` for legume steps. The `*TREATMENTS` section in a sequence file contains ordered steps (not independent treatments) — each step has the same `N` (treatment number) but increments `R` (sequence step).

---

## Concepts and conventions

### Point IDs — the canonical key

When the pipeline generates a grid, each point is assigned a **zero-padded 8-digit string ID** (e.g., `00000001`). This ID is used consistently for folder names, weather file names, and soil profile IDs.

### Folder structure

Step 3 creates one folder per point:
```
dssat_runs/<RUN_NAME>/<POINT_ID>/
```

Inside each folder:
- `SOIL.SOL` — the soil profile for that point
- `<POINT_ID>.<ext>` — a copy of your template FileX/SQX with placeholders replaced
- `<POINT_ID>.WTH` — the weather file for that point
- Support files copied from `dssat_templates/` (`*.CUL`, `*.ECO`, `*.SPE`, etc.)

### Template placeholders

When building each point folder, the pipeline performs simple string replacement on your template:

- Soil: if `TEMPLATE_SOIL_ID_PLACEHOLDER` (default `"SOIL_ID"`) is present in the template, it is replaced with the resolved soil profile ID.
- Point: `00000000` is replaced with the point ID (use this for `WSTA`, field IDs, etc.).
- Template filename base: occurrences of the template filename (without extension) are replaced with the point ID.

### Run folder naming

`DSSAT_RUN_NAME` controls the run directory name under `dssat_runs/`. Configure in Section 0:

```r
RUN_TAG        <- "run1"      # appended to the base name; leave "" for default
RUN_NAME_STYLE <- "grid"      # "grid"     => <GRID_BASE_NAME>_<RUN_TAG>
                               #              e.g., dssat_spatial_demo_50km_run1
                               # "scenario" => <GRID_BASE_NAME>_<WEATHER>_<SOIL>_<RUN_TAG>
                               #              e.g., dssat_spatial_demo_50km_DAYMET_SSURGO_run1
RUN_NAME_OVERRIDE <- ""       # if non-empty, this exact string is used as the run folder name
```

---

## Step 0 — Prepare DSSAT templates

At minimum, each point folder needs:

- **Experiment / sequence file** (`*.MZX`, `*.SQX`, etc. depending on crop and `RUN_MODE`)
- **Cultivar / ecotype / species files** required by your crop module

Copy a working DSSAT example into `dssat_templates/`. For a quick test, copy `UFGA8201.MZX` from your DSSAT installation's `Maize/` folder plus the supporting `*.CUL`, `*.ECO`, `*.SPE` files for the CERES-Maize module.

> **Key requirement:** your template must contain the **placeholder strings** that the pipeline replaces per point. Open the template in a text editor and confirm that the `*FIELDS` section has `WSTA` set to `00000000` and `ID_SOIL` set to `SOIL_ID`. The pipeline replaces these with the actual point ID and soil profile ID for each point.

---

## Step 1 — Configure and run the R pipeline

The main entrypoint is `dssat_main_pipeline.R`. The only section you need to edit for most runs is **SECTION 0: MASTER CONFIGURATION**.

Key settings:

```r
# --- Spatial domain ---
USE_EXISTING_POINT_SHAPEFILE  <- FALSE       # TRUE = Mode B/C; FALSE = Mode A
BOUNDARY_SHAPEFILE_NAME       <- "tl_2024_us_state.shp"
STATE_NAME_FILTER             <- c("Montana")
GRID_SPACING_METERS           <- 50000

# --- Weather ---
WEATHER_SOURCE     <- "DAYMET"               # "DAYMET", "NASA_POWER", or "GRIDMET"
WEATHER_START_YEAR <- 1999
WEATHER_END_YEAR   <- 2001

# --- Soil ---
SOIL_SOURCE        <- "SSURGO"               # "SSURGO", "SOILGRIDS_10K", or "SOILGRIDS_ONLINE"

# --- DSSAT template ---
TEMPLATE_FILE_NAME <- "UFGA8201.MZX"        # DEMO PLACEHOLDER — replace with your file
RUN_MODE           <- "experiment"           # "experiment" or "sequence"
TREATMENT_START    <- 1
TREATMENT_END      <- 4

# --- Execution mode ---
RUN_DSSAT_EXECUTION <- TRUE                  # FALSE = HPC prep only (no local execution)
ZIP_FOR_HPC         <- FALSE                 # TRUE = zip run folder for HPC upload
```

### Common first customisations (after the demo works)

- **Spatial scope:** change `GRID_SPACING_METERS`, `STATE_NAME_FILTER`, or switch to `USE_EXISTING_POINT_SHAPEFILE <- TRUE`.
- **Weather:** change `WEATHER_SOURCE` and extend `WEATHER_START_YEAR` / `WEATHER_END_YEAR`.
- **Soil:** change `SOIL_SOURCE`; if using `SOILGRIDS_10K`, set `EXTERNAL_SOIL_FILE`.
- **Crop:** update `CROP_EXTENSION` and `TEMPLATE_FILE_NAME` to match your crop module.
- **HPC prep:** set `RUN_DSSAT_EXECUTION <- FALSE` and `ZIP_FOR_HPC <- TRUE`.

### Running non-interactively (command line)

```bash
Rscript dssat_main_pipeline.R
```

---

## Validate one point folder before scaling

Before running hundreds or thousands of points, validate a **tiny run** end-to-end. This saves significant time.

### Recommended validation flow

1. Reduce to 5–20 points: increase `GRID_SPACING_METERS`, tighten `STATE_NAME_FILTER`, or use a small existing shapefile.
2. Keep the run minimal: `TREATMENT_START == TREATMENT_END`, 2–3 years of weather.
3. Set `CLEANUP_RUN_FOLDERS <- FALSE` so you can inspect DSSAT logs afterwards.
4. Click **Source** in RStudio.

### What a successful run looks like

Inspect `dssat_runs/<RUN_NAME>/00000001/` and confirm:

- `SOIL.SOL` exists and contains the soil profile ID referenced in the FileX/SQX.
- `00000001.WTH` exists and covers your simulation period.
- The per-point FileX/SQX exists (e.g., `00000001.MZX`) with placeholders replaced.
- DSSAT produced outputs (`summary.csv`, `LUN.LST`, etc.).
- `results_00000001.csv` exists.

### If it fails

- Inspect `LUN.LST` (and any `*.LST` / `*.OUT`) inside the point folder — DSSAT errors are described there.
- Set `DSSAT_CORES <- 1` to make debugging simpler.
- Confirm your template is configured to write **CSV outputs** (`summary.csv`), since the parser expects them.

---

## Option A — Run locally

### Option A1: Run from R (default)

If `RUN_DSSAT_EXECUTION <- TRUE` in Section 0, the pipeline runs DSSAT per point in parallel using `doParallel`, parses outputs, and writes the merged results CSV and yield map automatically.

### Option A2: Run the MPI runner locally

This is the best way to test the HPC execution path on a workstation before submitting to a cluster:

```bash
mpirun -n 4 python dssat_mpi_runner.py \
  --base_dir "dssat_runs/<RUN_NAME>" \
  --summary_dir "results" \
  --exe_path "/path/to/dscsm048" \
  --run_mode sequence \
  --trt_start 1 --trt_end 5 \
  --seq_start 1 --seq_end 2 \
  --cleanup_mode never \
  --archive_outputs
```

---

## Option B — Run on HPC / cloud (SLURM + MPI)

### 1 — HPC smoke test: verify your environment

Before submitting a large gridded job, confirm that MPI, mpi4py, and DSSAT all work on your cluster.

**Quick interactive checks:**

```bash
# MPI + mpi4py
python -c "import mpi4py; print('mpi4py OK')"
mpirun -n 2 python -c "from mpi4py import MPI; c=MPI.COMM_WORLD; print(f'rank={c.Get_rank()} size={c.Get_size()}')"

# DSSAT
dscsm048        # should print a usage banner or version
```

**Tiny SLURM test job** — create `hpc_smoke_test.slurm`:

```bash
#!/bin/bash
#SBATCH --job-name=dssat_smoke
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --output=slurm-smoke-%j.out

# Load site modules and activate environment (adapt to your cluster)
# module load python openmpi
# source activate dssat_env

export OMP_NUM_THREADS=1
export PYTHONUNBUFFERED=1
export DSSAT_EXE="/path/to/dscsm048"

echo "=== MPI test ==="
mpirun -n ${SLURM_NTASKS} python -c \
  "from mpi4py import MPI; c=MPI.COMM_WORLD; print(f'rank={c.Get_rank()} size={c.Get_size()}')"

echo "=== DSSAT test ==="
"$DSSAT_EXE" || true

echo "Smoke test complete."
```

```bash
sbatch hpc_smoke_test.slurm
squeue -u $USER
```

### 2 — Build run folders and zip (R)

In Section 0 set `RUN_DSSAT_EXECUTION <- FALSE` and `ZIP_FOR_HPC <- TRUE`, then:

```bash
Rscript dssat_main_pipeline.R
```

This builds all point folders and produces `dssat_runs/<RUN_NAME>.zip`.

### 3 — Transfer to HPC scratch

```bash
scp dssat_runs/<RUN_NAME>.zip user@cluster:/scratch/user/
ssh user@cluster
cd /scratch/user && unzip <RUN_NAME>.zip
```

### 4 — Submit the SLURM job

Edit `run_dssat_python.slurm` to set `DSSAT_EXE`, `BASE_DIR`, `SUMMARY_DIR`, node count, and walltime. Then:

```bash
sbatch run_dssat_python.slurm
squeue -u $USER
```

> Keep the SLURM script **generic**: avoid committing personal email addresses, allocation keys, or user-specific paths. Use placeholders or environment variables.

### 5 — End-to-end mini run (strongly recommended before scaling)

```bash
mpirun -n 2 python dssat_mpi_runner.py \
  --base_dir "dssat_runs/<RUN_NAME>" \
  --summary_dir "results" \
  --exe_path "$DSSAT_EXE" \
  --run_mode sequence \
  --trt_start 1 --trt_end 1 \
  --seq_start 1 --seq_end 1 \
  --cleanup_mode success
```

If this succeeds, you are ready to scale up.

---

## Outputs and what they mean

**Per-point:** `dssat_runs/<RUN_NAME>/<point_id>/results_<point_id>.csv`

**Merged (rank 0):** `<summary_dir>/results_<RUN_NAME>.csv`

Columns include:
- IDs and coordinates
- Planting, emergence, and harvest dates
- Grain yield and above-ground biomass (kg/ha)
- Irrigation totals (mm)
- Nitrogen application and leaching (kg/ha)
- CO₂ and N₂O emission metrics (if present in DSSAT outputs)
- Seasonal average water and nitrogen stress indices (from PlantGro)

---

## Post-processing: join results back to maps

### Option 1 (recommended): join by `point_id` in R

```r
library(sf)
library(dplyr)

pts <- st_read("gridpoints/<YOUR_GRIDPOINTS>.shp")
res <- read.csv("results/results_<RUN_NAME>.csv", stringsAsFactors = FALSE)

pts <- pts %>% mutate(point_id = sprintf("%08d", as.integer(ID)))
pts_res <- pts %>% left_join(res, by = "point_id")

plot(pts_res["final_grain_kg_ha"])
st_write(pts_res, "results/results_<RUN_NAME>.gpkg", delete_dsn = TRUE)
```

### Option 2: rebuild points from latitude / longitude

```r
library(sf)
res <- read.csv("results/results_<RUN_NAME>.csv", stringsAsFactors = FALSE)
pts_res <- st_as_sf(res, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
plot(pts_res["final_grain_kg_ha"])
```

### Mapping tips

- **Units:** most DSSAT outputs are in **kg/ha** and **mm**. Confirm assumptions before plotting.
- **Sequence runs:** decide whether you want all seasons, a specific crop step, or per-point mean yield across rotation cycles.
- **Large runs:** stream or chunk instead of loading the full results table into memory.

---

## How the MPI Python runner works

1. Rank 0 scans folders under `--base_dir` and identifies point folders (digit-only names).
2. MPI ranks split the folder list using simple rank slicing (`folders[rank::size]`).
3. Each rank runs DSSAT per point and per treatment, and writes `temp_results_rank_<rank>.csv`.
4. Rank 0 merges all temp files into the final `results_<RUN_NAME>.csv`.

The runner supports optional output cleanup (`--cleanup_mode never/success/always`) and output archiving (`--archive_outputs`).

---

## Troubleshooting

| Problem | Likely cause | Solution |
|---------|-------------|----------|
| `CRITICAL: Boundary shapefile not found` | Shapefile not downloaded or path wrong | Run the `curl` + `unzip` commands in Quick Start; confirm `BOUNDARY_SHAPEFILE_NAME` matches the `.shp` filename |
| `CRITICAL: DSSAT Template file not found` | Template missing from `dssat_templates/` | Copy `UFGA8201.MZX` from `DSSAT48/Maize/` for the demo; replace with your own file for production |
| DSSAT executable not found at runtime | `DSSAT_EXE` env var not set | `Sys.setenv(DSSAT_EXE = "/full/path/to/dscsm048")` before sourcing |
| `No .SQX file found` (MPI runner) | MPI runner expects `.SQX`; you have a `.MZX` | Switch to a `.SQX` template or extend the runner to support seasonal FileX |
| `DSSAT completed but no parsable outputs` | `summary.csv` missing or empty | Rerun with `--cleanup_mode never`; inspect `*.LST` / `LUN.LST`; ensure template writes CSV outputs |
| `No data returned from NASA-POWER` | Point is over ocean or invalid coordinates | Filter out water pixels before running or use Mode C |
| Weather download stalls or rate-limit errors | Too many parallel requests | Set `N_CORES_TO_USE <- 1`; SoilGrids REST back-off retry handles 429 errors automatically |
| Soil processing skips many points | Critical data (sand/clay/BD) is NA | Check `soil_processing_errors.log` in the soil output folder; consider switching soil source |
| Soil ID mismatch in DSSAT log | Placeholder `SOIL_ID` not replaced | Check `TEMPLATE_SOIL_ID_PLACEHOLDER` in Section 0; confirm the placeholder string is in your template |
| Weather does not cover full simulation period | Simulation years extend beyond downloaded range | For sequence runs: a 10-year rotation needs at least 10 + `NYERS` years of weather; extend `WEATHER_END_YEAR` or use the weather extension functions |
| `mpi4py` import error on HPC | Built against wrong MPI library | Rebuild: `pip install --no-cache-dir --no-binary :all: mpi4py` after loading the correct MPI module |
| Results CSV is empty after run | Parsing failed; DSSAT wrote `.OUT` only | Check per-point `LUN.LST`; enable CSV output in the `*OUTPUTS` section of your template |

---

## Performance tips

- Use **scratch storage** (fast local SSD or parallel filesystem) for `base_dir` — DSSAT creates many small files per run.
- Use `--cleanup_mode always` in production to remove transient DSSAT files after each point.
- Always **debug with a single point** first (`DSSAT_CORES <- 1`).
- For very large domains, set `ZIP_FOR_HPC <- TRUE` and run on HPC rather than locally.

---

## Reproducibility checklist

- Record DSSAT version and executable path used.
- Record the template file(s) and all cultivar/ecotype/species files.
- Record soil and weather data sources and year ranges.
- Save the `README_CONFIG.txt` written by the pipeline to the run folder (it captures key Section 0 settings automatically).
- Commit Section 0 changes to version control so every run can be traced back to a specific configuration.

---

## Publishing checklist

**1 — DSSAT executables and licensed content**
- Official DSSAT downloads (via dssat.net) are tied to a licence. **Do not** redistribute proprietary DSSAT installers or binaries.
- The open-source DSSAT-CSM-OS codebase (BSD 3-Clause) may be redistributed as compiled binaries if you keep the licence and document provenance (commit hash, build flags, platform).
- Only include template files (`*.SQX`, `*.MZX`, etc.) if you have the right to redistribute them.

**2 — Sanitise cluster-specific values**
Avoid committing personal email addresses (`#SBATCH --mail-user=...`), allocation keys (`--account`), or user-specific paths. Use placeholders or environment variables.

**3 — Keep large generated data out of Git**
`dssat_runs/`, `weather/`, `soil/`, and `results/` can be massive. The provided `.gitignore` is configured to ignore them.

**4 — Provide small example data**
For a public release, include a tiny boundary or points file (5–20 points) and a short example weather period so users can reproduce an end-to-end run quickly.

**5 — Cite data sources**
If your pipeline downloads weather/soil data, include citations and terms of use for gridMET / Daymet / NASA POWER, SSURGO / SoilGrids, CDL / NLCD, and any other sources used.

---

## Installing DSSAT

### Windows

1. Download the official installer from https://dssat.net/download/
2. Run the installer and accept the default path (`C:\DSSAT48\`).
3. The executable is `C:\DSSAT48\DSCSM048.exe`.
4. In R:
   ```r
   Sys.setenv(DSSAT_EXE = "C:/DSSAT48/DSCSM048.exe")
   ```

---

### macOS (Apple Silicon M1/M2/M3/M4/M5)

DSSAT has no official macOS installer. You must compile from source using Homebrew GCC and CMake. Tested on Apple Silicon running macOS Ventura 14+.

#### Part 1 — Install build tools

```bash
xcode-select --install          # Xcode Command Line Tools
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
brew install gcc cmake git
which gfortran && gfortran --version
```

#### Part 2 — Clone repositories

```bash
cd ~/Documents/GitHub
git clone https://github.com/DSSAT/dssat-csm-os.git
git clone https://github.com/DSSAT/dssat-csm-data.git
mkdir -p ~/Documents/GitHub/DSSAT48
cp -r ~/Documents/GitHub/dssat-csm-data/. ~/Documents/GitHub/DSSAT48/
```

#### Part 3 — Compile

```bash
cd ~/Documents/GitHub/dssat-csm-os
mkdir build && cd build
cmake .. \
  -DCMAKE_BUILD_TYPE=RELEASE \
  -DCMAKE_Fortran_COMPILER=$(which gfortran) \
  -G "Unix Makefiles"
make -j$(sysctl -n hw.logicalcpu)
```

A successful build produces `build/bin/dscsm048`.

#### Part 4 — Install executable and support files

```bash
cp build/bin/dscsm048 ~/Documents/GitHub/DSSAT48/dscsm048
chmod +x ~/Documents/GitHub/DSSAT48/dscsm048
cp ~/Documents/GitHub/DSSAT48/{MODEL.ERR,OUTPUT.CDE,DATA.CDE} \
   ~/Documents/GitHub/DSSAT48/StandardData/
```

#### Part 5 — Add to PATH

```bash
echo 'export PATH="/Users/YOUR_USERNAME/Documents/GitHub/DSSAT48:$PATH"' >> ~/.zprofile
source ~/.zprofile
which dscsm048
```

Replace `YOUR_USERNAME` with your macOS username throughout.

#### Part 6 — Generate DSSATPRO.L48

```bash
sed 's|/usr/local|/Users/YOUR_USERNAME/Documents/GitHub/DSSAT48|g' \
    ~/Documents/GitHub/dssat-csm-os/Data/DSSATPRO.L48 \
    > ~/Documents/GitHub/DSSAT48/DSSATPRO.L48
grep "StandardData" ~/Documents/GitHub/DSSAT48/DSSATPRO.L48
```

Rules: file must be named `.L48` (not `.v48`); all paths must use forward slashes; no `.EXE` suffix.

#### Part 7 — Test

```bash
cd ~/Documents/GitHub/DSSAT48/Maize
dscsm048 A UFGA8201.MZX 1
```

A successful run produces `Summary.OUT`.

#### Directory layout after install

```
~/Documents/GitHub/DSSAT48/
├── dscsm048
├── DSSATPRO.L48
├── MODEL.ERR / OUTPUT.CDE / DATA.CDE
├── StandardData/
├── Maize/
├── Wheat/
└── ...
```

#### Troubleshooting (macOS)

| Symptom | Fix |
|---------|-----|
| `cmake: No Fortran compiler found` | Pass `-DCMAKE_Fortran_COMPILER=$(which gfortran)` explicitly |
| `make` fails with linker errors | Set `export CC=$(brew --prefix gcc)/bin/gcc-14` before cmake |
| `dscsm048: command not found` | Run `source ~/.zprofile` or open a new terminal tab |
| `ERROR: DSSATPRO.L48 not found` | Confirm file is in the same directory as `dscsm048` |
| `ERROR: Cannot find StandardData` | Re-run the `sed` command; verify with `grep StandardData DSSATPRO.L48` |
| Test run exits with no output | `cd` into the crop subdirectory (e.g., `Maize/`) before running |
| R `DSSAT` package cannot find executable | Add `options(DSSAT.CSM = Sys.getenv("DSSAT_EXE"))` before sourcing the pipeline |

Resources: https://github.com/DSSAT/dssat-csm-os · https://dssat.net/forum/ · https://brew.sh

---

### Linux

```bash
sudo apt-get install gcc gfortran cmake git   # Debian / Ubuntu
# or: sudo yum install gcc gcc-gfortran cmake git   # RHEL / Rocky / CentOS

git clone https://github.com/DSSAT/dssat-csm-os.git
git clone https://github.com/DSSAT/dssat-csm-data.git
sudo mkdir -p /opt/DSSAT48
sudo cp -r dssat-csm-data/. /opt/DSSAT48/

cd dssat-csm-os && mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE \
         -DCMAKE_Fortran_COMPILER=$(which gfortran) \
         -G "Unix Makefiles"
make -j$(nproc)

sudo cp bin/dscsm048 /opt/DSSAT48/ && sudo chmod +x /opt/DSSAT48/dscsm048
sudo cp /opt/DSSAT48/{MODEL.ERR,OUTPUT.CDE,DATA.CDE} /opt/DSSAT48/StandardData/
sed 's|/usr/local|/opt/DSSAT48|g' ../Data/DSSATPRO.L48 > /opt/DSSAT48/DSSATPRO.L48
echo 'export PATH="/opt/DSSAT48:$PATH"' >> ~/.bashrc && source ~/.bashrc
```

For a user-level install replace `/opt/DSSAT48` with `~/DSSAT48` and omit `sudo`.

In R:
```r
Sys.setenv(DSSAT_EXE = "/opt/DSSAT48/dscsm048")
```

Official compile guide: https://dssat.net/source-code/

---

## R package requirements

```r
# Core (always required)
install.packages(c(
  "sf", "dplyr", "tidyr", "stringr", "lubridate",
  "foreach", "doParallel", "parallel",
  "zoo", "R.utils", "processx", "tools",
  "ggplot2", "readr", "tibble", "rstudioapi"
))

# DSSAT R interface
install.packages("DSSAT")
# If not on CRAN: remotes::install_github("palderman/DSSAT")

# Soil packages
install.packages("soilDB")                           # SSURGO
install.packages(c("terra", "httr", "jsonlite"))     # SoilGrids online / VRT

# Weather — install the package(s) matching your WEATHER_SOURCE
install.packages("daymetr")                          # DAYMET
install.packages("nasapower")                        # NASA_POWER
install.packages(c("terra", "ncdf4", "httr"))        # GRIDMET

# Optional
install.packages("pbapply")                          # progress bars
```

---

## Advanced: HPC Python environment setup

If you need a dedicated Python environment for the MPI runner, here are two common patterns. Adapt module names to your cluster.

### Option 1: venv + pip

```bash
module load python gcc openmpi   # adapt to your cluster's module names

python -m venv ~/envs/dssat_env
source ~/envs/dssat_env/bin/activate
pip install --upgrade pip
pip install numpy pandas tqdm shapely pyproj pyogrio fiona geopandas

# Build mpi4py against the loaded MPI (important — do not use a pre-built wheel)
pip install --no-cache-dir --no-binary :all: mpi4py

# Register as Jupyter kernel (optional)
pip install ipykernel
python -m ipykernel install --user --name=dssat_env --display-name "DSSAT env"
```

### Option 2: conda / mamba

```bash
module load python gcc openmpi

mamba create -y --prefix ~/envs/dssat_env python=3.11 ipykernel numpy pandas tqdm geopandas
conda activate ~/envs/dssat_env
pip install --no-cache-dir --no-binary :all: mpi4py
python -m ipykernel install --user --name=dssat_env --display-name "DSSAT env"
```

### Validation

```bash
python -c "import numpy, pandas; print('python OK')"
python -c "import geopandas; print('geopandas OK')"
python -c "from mpi4py import MPI; c=MPI.COMM_WORLD; print('rank', c.Get_rank(), 'size', c.Get_size())"
```

---

## Related tools and gridded crop modeling ecosystems

### Other crop models commonly used for multi-site simulation

- **APSIM Next Generation:** widely used cropping systems framework; supports complex rotations.
  - https://www.apsim.info/apsim-next-generation/ and https://github.com/APSIMInitiative/ApsimX
  - R interface (`apsimx`): https://github.com/femiguez/apsimx
- **AquaCrop-OSPy:** crop-water model for water-limited yield estimation and irrigation studies.
  - https://aquacropos.github.io/aquacrop/
- **PCSE (Python Crop Simulation Environment):** includes WOFOST and other models.
  - https://pcse.readthedocs.io/

### Frameworks for gridded / ensemble impact modeling

- **pSIMS:** HPC-oriented framework for running site-based impact models over a geospatial grid.
  - Paper: https://www.sciencedirect.com/science/article/pii/S1364815214001121
- **AgMIP Ag-GRID:** protocols for gridded crop modeling and multi-model ensembles.
  - https://agmip.org/ag-grid-2/
- **ISIMIP:** global framework for intercomparing sectoral impact models.
  - https://www.isimip.org/

### DSSAT-specific automation tools

- **DSSATTools (Python):** builds DSSAT inputs and runs DSSAT from Python.
  - https://py-dssattools.readthedocs.io/ and https://github.com/daquinterop/Py_DSSATTools
- **Pythia:** DSSAT automation and spatial pipeline framework.
  - https://github.com/DSSAT/pythia
- **Alderman (2021):** parallel gridded DSSAT-CSM using MPI and NetCDF.
  - https://doi.org/10.5194/gmd-14-6541-2021

> These tools have different licensing models, assumptions, and supported processes. Always validate model outputs against field observations before applying at scale.

---

## References

**DSSAT**
- Official site: https://dssat.net/
- Command-line / power-user tools: https://dssat.net/tools/tools-for-power-users/
- Weather module (CO₂ specification options): https://dssat.net/weather-module/
- Automatic irrigation guide: https://dssat.net/wp-content/uploads/2020/03/Users-Guide-to-Automatic-Irrigation-using-Growth-Stage-Controls_v7.pdf
- Source code: https://github.com/DSSAT/dssat-csm-os
- Data repository: https://github.com/DSSAT/dssat-csm-data

**Weather datasets**
- gridMET: https://www.climatologylab.org/gridmet.html
- Daymet: https://daymet.ornl.gov/
- NASA POWER daily API: https://power.larc.nasa.gov/docs/services/api/temporal/daily/

**Soil datasets**
- SSURGO (NRCS): https://www.nrcs.usda.gov/resources/data-and-reports/soil-survey-geographic-database-ssurgo
- soilDB R package: https://ncss-tech.github.io/soilDB/
- SoilGrids portal: https://soilgrids.org/
- SoilGrids docs: https://docs.isric.org/globaldata/soilgrids/
- SoilGrids 10K pre-formatted DSSAT files (Harvard Dataverse): https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/1PEEY0
- Folberth et al. (2019) *Environ. Model. Softw.* 111:218–228: https://www.sciencedirect.com/science/article/pii/S1364815218313033

**Landcover**
- NLCD download: https://www.mrlc.gov/data
- NLCD class legend: https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description
- USDA CDL (CroplandCROS): https://croplandcros.scinet.usda.gov/
- CDL national annual download: https://www.nass.usda.gov/Research_and_Science/Cropland/Release/index.php
- CDL class codes: https://www.nass.usda.gov/Research_and_Science/Cropland/sarsfaqs2.php
- CDL via Google Earth Engine: https://developers.google.com/earth-engine/datasets/catalog/USDA_NASS_CDL
- ESA WorldCover (global, 10 m): https://worldcover2021.esa.int/

**Spatial boundaries**
- US Census TIGER/Line states: https://www2.census.gov/geo/tiger/TIGER2024/STATE/tl_2024_us_state.zip
- US Census TIGER/Line counties: https://www2.census.gov/geo/tiger/TIGER2024/COUNTY/tl_2024_us_county.zip
