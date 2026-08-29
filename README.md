# Spatial Gridded Crop Modeling with DSSAT

> **AI agents & maintainers:** read [`../AGENTS.md`](../AGENTS.md) before editing this repo.

This repository provides an end-to-end, beginner-friendly workflow for **spatial (gridded / point-based) crop modeling in DSSAT**. The pipeline downloads weather and soil data for a set of geographic points, builds a DSSAT-ready input folder for every point, runs the simulations locally or prepares them for HPC, and maps the results.

**Start here:** open `dssat_main_pipeline.R` in **RStudio** and click **Source** to run a small spatial demo locally (no command line needed). Once that works you can scale up on HPC/cloud using **SLURM + MPI**.

> Scripts and examples are a work in progress; several aspects are not yet optimised or finalised.

---

## Table of Contents

- [Status, limitations, and future work](#status-limitations-and-future-work)
- [What is DSSAT and why spatial modeling?](#what-is-dssat-and-why-spatial-modeling)
- [Prerequisites](#prerequisites)
- [How to start on a fresh Windows laptop](#how-to-start-on-a-fresh-windows-laptop)
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
- [Sensitivity experiments: sweep weather × soil combinations](#sensitivity-experiments-sweep-weather--soil-combinations)
- [Validation against observed data](#validation-against-observed-data)
- [Repository layout](#repository-layout)
- [Background: DSSAT inputs and run-folder anatomy](#background-dssat-inputs-and-run-folder-anatomy)
  - [Regular vs crop-sequence runs](#regular-vs-crop-sequence-runs)
  - [Simulation Controls explained](#simulation-controls-explained)
- [Concepts and conventions](#concepts-and-conventions)
- [Step 0 — Prepare DSSAT templates](#step-0--prepare-dssat-templates)
- [Step 1 — Configure and run the R pipeline](#step-1--configure-and-run-the-r-pipeline)
- [Validate one point folder before scaling](#validate-one-point-folder-before-scaling)
- [Testing](#testing)
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

### Scaling features (HPC MPI runner)

The MPI runner (`hpc/dssat_mpi_runner.py`) was built to scale to large grids:

- **Dynamic work-stealing:** a manager hands points to workers on demand (not static `folders[rank::size]` slicing), so slow points don't stall a rank.
- **Streaming per-point writes:** each rank streams results straight to its part file instead of accumulating them in memory.
- **Node-local scratch staging** (`--scratch_dir`, defaults to `$TMPDIR`): each rank stages its small files on fast local disk to avoid a metadata storm on Lustre/GPFS, then copies results back.
- **Optional merge skip** (`--merge_mode none`): skip the rank-0 serial merge on very large grids and merge the per-rank parts later.
- **Tiny run folders:** with `DSSATPRO.V48` next to the executable, genotype/support files are resolved from the install instead of copied per point (see [Performance tips](#performance-tips)) — far fewer files for the filesystem to track.

### Known limitations

- **Assumes CSV outputs:** the parser expects DSSAT CSV output files (`summary.csv`, `soilwat.csv`, etc.). If your DSSAT setup produces only `*.OUT` text files, you must enable CSV outputs or adapt the parser.
- **Per-point network sources:** Daymet / NASA POWER / Open-Meteo make one HTTP request per point (parallelised, cached). For very large grids prefer a gridded source (GridMET, CHIRPS, AgERA5) that downloads once and extracts all points locally.
- **AgERA5 licence + queue:** AgERA5 needs a Copernicus CDS key and a one-time licence acceptance, and its requests are queued server-side (first run for a new area/period can take minutes).

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

## How to start on a fresh Windows laptop

End-to-end setup on a clean Windows machine. The R environment is pinned with
[`renv`](https://rstudio.github.io/renv/), so `renv::restore()` reproduces the
exact package versions — including the shared `dssatengine` / `dssatutils`
packages, which install straight from GitHub (no side-by-side clones needed).

### 1 — Install the non-R prerequisites

1. **R** — install the version recorded in `renv.lock` (currently R 4.6.x) from
   [CRAN](https://cran.r-project.org/bin/windows/base/).
2. **RStudio** — [Posit Desktop](https://posit.co/download/rstudio-desktop/);
   the pipeline uses `rstudioapi` to locate itself.
3. **Rtools** — [matching your R version](https://cran.r-project.org/bin/windows/Rtools/).
   Most packages install as binaries, but Rtools is a safe fallback for any
   source build.
4. **Git** — [git-scm.com](https://git-scm.com/download/win). The shared
   packages `dssatutils` and `dssatengine` install from GitHub. If Git prompts
   for authentication (for private repos, rate limits, or account-specific access),
   authenticate once with your GitHub account:
   * **For R / `renv`:**
     Generate a Personal Access Token (PAT) by running `usethis::create_github_token()` in R.
     Store the token by running `gitcreds::gitcreds_set()`. Alternatively, you can add `GITHUB_PAT=your_token_here` directly to your local `~/.Renviron` file.
   * **For Python / `pip`:**
     Ensure your Git Credential Manager is active (it will prompt for authentication if needed), or install using your PAT directly if access requires it:
     `pip install "git+https://github.com/alwinhopf/dssatutils.git@e9c859fa1d915623df23e2eb13084cb085dbfe3e"`
     PAT fallback:
     `pip install "git+https://<PAT>@github.com/alwinhopf/dssatutils.git@e9c859fa1d915623df23e2eb13084cb085dbfe3e"`
     Use the CDS extra for Copernicus-backed weather sources:
     `pip install "dssatutils[cds] @ git+https://github.com/alwinhopf/dssatutils.git@e9c859fa1d915623df23e2eb13084cb085dbfe3e"`
     If using SSH, verify your SSH keys are added to your GitHub account: `ssh -T git@github.com`.
5. **DSSAT 4.8** — install from [dssat.net](https://dssat.net) to the default
   **`C:\DSSAT48`**. The pipeline auto-detects `C:\DSSAT48\DSCSM048.EXE` on
   Windows. This is the crop-model executable; `renv` does **not** provide it. If installed elsewhere, see [Troubleshooting](#troubleshooting).

### 2 — Clone the project and restore the R environment

```powershell
git clone https://github.com/alwinhopf/DSSAT-Gridded-Run-Tutorial.git DSSAT_Gridded_Run_Tutorial
cd DSSAT_Gridded_Run_Tutorial
```

Open the folder in RStudio (the project's `.Rprofile` auto-bootstraps `renv`),
then install every pinned package — CRAN dependencies **and** the
`dssatengine` / `dssatutils` packages from GitHub — in one step:

```r
renv::restore()
```

Optional for Copernicus CDS weather sources (`AGERA5`, `ERA5_LAND`, and E-OBS
CDS mode): with `dssatutils` v0.4.0 or newer, configure the CDS Personal Access
Token once after restore:

```r
library(dssatutils)
setup_cds_credentials()
```

Python users can run the same setup helper:

```python
from dssatutils import setup_cds_credentials
setup_cds_credentials()
```

The helper uses `CDSAPI_KEY` / `CDSAPI_URL`, imports an existing `~/.cdsapirc`,
or prompts in an interactive session. It writes a `cdsapi`-compatible
`.cdsapirc`; the R helper also stores the token for `ecmwfr`.

This replaces the old manual `install.packages(..., type = "binary")` /
`devtools::install_local(...)` dance. The project `.Rprofile` sets
`options(pkgType = "binary")` on Windows/macOS **before** `renv` installs
anything, so every package — including transitive system-library dependencies
like `openssl`, `curl`, `sf`, and `terra` — is pulled as a prebuilt binary and
never compiled from source. (Source builds of those packages are the usual cause
of a cascade failure such as `openssl → httr → daymetr → dssatutils` on a clean
laptop, because they need system dev headers that a fresh machine lacks.)

> **If `renv::restore()` still tries to build something from source** (rare —
> usually a package whose binary isn't yet on CRAN for your R version), install
> just that one as a binary and re-snapshot:
> ```r
> renv::install("openssl", type = "binary")   # name the offending package
> renv::snapshot()                             # update renv.lock so it sticks
> ```

> **Generating the lock (one-time, maintainers only):** if `renv.lock` does not
> yet pin `dssatengine` / `dssatutils`, run `Rscript setup_renv.R` once on a
> machine where those repos are committed and pushed, then commit the
> regenerated `renv.lock`.

### 3 — Get the demo boundary shapefile (one-time)

Follow [Quick start → step 3](#3--download-the-demo-boundary-shapefile-one-time)
to download `tl_2024_us_state.zip` into `shapefile/`.

### 4 — Run it

Open `dssat_main_pipeline.R` and click **Source**. DSSAT is auto-detected on
Windows; override only if you installed it elsewhere:

```r
Sys.setenv(DSSAT_EXE = "C:/DSSAT48/DSCSM048.exe")
```

A successful run writes a combined results CSV and a yield map under `results/`.
If any soil or weather point cannot be downloaded, the pipeline now reports it
explicitly and writes a `*_download_failures.csv` log next to the data, so a
point landing on water or outside coverage is auditable rather than a silent gap
(see [Troubleshooting](#troubleshooting)).

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

The repository default is the included `dssat_templates/UFGA9201.MZX`. To use
another experiment, copy a known-good FileX into `dssat_templates/`, then update
`template_file_name`, `crop_extension`, and the weather years together. Any
valid DSSAT experiment file works (`.MZX`, `.WHX`, `.SBX`, `.SQX`, etc.) when
its crop extension, planting dates, and `run_mode` agree with the configuration.

### 5 — Source the pipeline

Click **Source** in RStudio. With the committed default settings (Iowa at 20 km
spacing, 1992 NASA POWER weather, SSURGO_ALDERMAN soil, `UFGA9201.MZX`, and a
local DSSAT run), the script will:

- Generate a few hundred grid points across Iowa
- Download soil and weather data for each point
- Build one DSSAT run folder per point
- Run DSSAT locally
- Parse and merge results
- Write a merged results CSV and a yield map under `results/`

The first run can take from tens of minutes to several hours because soil and
weather services are live and rate-limited. For the first validation, use a
coarser grid such as 300 km; this normally yields a single Iowa point.

### 6 — Confirm you got results

After the run completes, look for:

- **Merged results:** `results/<RUN_NAME>_results.csv`
- **Maps:** `results/<RUN_NAME>_yield_map_treatment<k>.png`
- **Per-point folders:** `dssat_runs/<RUN_NAME>/<POINT_ID>/` (includes `SOIL.SOL`, `*.WTH`, the per-point FileX/SQX, and DSSAT outputs)

### 7 — If the first run fails

- DSSAT path not found → confirm `Sys.setenv(DSSAT_EXE = ...)` is set correctly.
- Missing boundary shapefile → download the TIGER/Line file (Step 3 above) or point the script to your own boundary.
- Template/run-mode mismatch (e.g., `run_mode: "sequence"` but the template is `*.MZX`) → see [Background: DSSAT inputs](#background-dssat-inputs-and-run-folder-anatomy) and update `run_mode` or `template_file_name` in `config.yml`.

Once the demo works, customise settings in the central **`config.yml`** file.

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

All three modes produce the same standardised point shapefile (with `ID` / `LAT` / `LONG` columns) that the rest of the pipeline consumes.

Set the controlling flag in `config.yml`:

```yaml
use_existing_point_shapefile: false   # MODE A — generate a regular grid
# use_existing_point_shapefile: true  # MODE B / C — supply your own points
```

---

### Mode A — Regular grid from a boundary polygon (default demo)

The pipeline reads a boundary shapefile, filters to the region(s) you want, and places points on a regular grid at `grid_spacing_meters` spacing.

```yaml
use_existing_point_shapefile: false
boundary_shapefile_name: "tl_2024_us_state.shp"
enable_boundary_filter: true
boundary_filter_column: "NAME"
state_name_filter: ["Montana"]
grid_spacing_meters: 50000  # 50 km fast demo; 5000 for production
```

For resolution experiments, Mode A can optionally use a persistent nested
master lattice instead of independently placing every grid:

```yaml
use_existing_point_shapefile: false
use_master_grid: true
master_grid_spacing_meters: 1000
master_grid_crs: "EPSG:5070"  # CONUS; EPSG:6933 is the global default
master_grid_origin_x: 0
master_grid_origin_y: 0
master_grid_phase_row: 0
master_grid_phase_col: 0
master_grid_path: ""           # automatic GeoPackage under gridpoints/
reuse_master_grid: true
grid_spacing_meters: 100000    # must be an integer multiple of 1000
```

The master GeoPackage is generated once. Every target grid selects rows and
columns at an integer stride, so coarse point IDs and coordinates are retained
in all compatible finer grids. Generated point files include `MROW`, `MCOL`,
`MSPACE_M`, `SAMP_M`, and `NEST_F` for auditing; these fields are also carried
into combined simulation results. Set `use_master_grid: false`
to restore the original independent placement. `use_master_grid` and
`use_existing_point_shapefile` are mutually exclusive.

The nested lattice controls sampling locations; it does not change the native
resolution of weather or soil products. When cropland filtering is active, the
pipeline still calculates cropland support at each target cell spacing after
the nested points have been selected. Two requested levels are strict subsets
of one another only when the coarser spacing is an integer multiple of the
finer spacing; otherwise they remain separate branches of the same master.
For strict post-filter nesting, set `cropland_filter_basis: "point"`; membership
then depends on the land-cover pixel at the shared point while target-cell
fractions and hectares remain available as weights. The legacy
`"cell_fraction"` basis can retain a coarse cell because it contains cropland
even when its anchor is not cropland, so filtered membership may differ by level.

---

### Mode B — Your own shapefile

Supply any point or polygon shapefile. Polygons are auto-converted to centroids.

```yaml
use_existing_point_shapefile: true
existing_point_shapefile_path: "gridpoints/my_study_sites.shp"
```

The `load_existing_points()` function handles the rest automatically:
- Polygons / lines → converted to centroids
- Any CRS → re-projected to WGS84 (EPSG:4326)
- ID column → standardised to zero-padded 8-digit strings
- Missing / duplicate IDs → regenerated sequentially

---

### Mode C — Cropland-only points (CDL / NLCD)

The preferred workflow is now integrated into `dssat_main_pipeline.R` and
`dssat_main_pipeline.py`: keep Mode A/Mode B as usual, then enable the optional
cropland mask in `config.yml`. The pipeline computes cropland percentage for
each grid cell, keeps only cells with cropland, and carries cropland area into
the final results CSV.

#### Step 1 — Download a landcover raster

| Dataset | Class code | Coverage | URL |
|---------|-----------|----------|-----|
| Annual NLCD Land Cover CONUS | `82` = Cultivated Crops (all crops) | Contiguous US, ~30 m | https://www.mrlc.gov/data |
| CDL | Crop-specific codes (e.g., `1` = corn) | Contiguous US, ~30 m | https://croplandcros.scinet.usda.gov/ |
| ESA WorldCover | `40` = Cropland | Global, 10 m | https://worldcover2021.esa.int/ |

For Annual NLCD, use the MRLC data portal: open
https://www.mrlc.gov/data, choose **Annual NLCD → Land Cover → Land Cover
(CONUS) 2024**, download the bundle, unzip it, and place the `.tif` under
`data/`. The current direct bundle URL is
https://www.mrlc.gov/downloads/sciweb1/shared/mrlc/data-bundles/Annual_NLCD_LndCov_2024_CU_C1V1.zip,
but MRLC may change this URL in future years.

#### Step 2 — Enable cropland filtering in `config.yml`

For all cultivated cropland with NLCD, use class `82`:

```yaml
use_cropland_mask:      true
cropland_raster_file:   "data/Annual_NLCD_LndCov_2024_CU_C1V1.tif"
cropland_classes:       [82]
cropland_min_fraction:  0      # keep grid cells with >0 cropland
cropland_strict:        false  # warn and continue all-land if raster is missing
reuse_cropland_grid:    true
```

The output grid shapefile stores:

- `crop_frac` — cropland fraction from 0 to 1
- `crop_pct` — cropland percentage from 0 to 100
- `crop_ha` — cropland hectares inside the grid cell
- `cell_ha` — total grid-cell hectares

Cropland grids are saved separately from all-land grids and reused on later
runs with the same project, resolution, class list, and threshold. For example,
`carinata_sweep_100km.shp` is the all-land grid, while
`carinata_sweep_100km_cropland_82_min0.shp` is the NLCD class 82 cropland grid.
Set `reuse_cropland_grid: false` to force regeneration.

The final results CSV also includes `cropland_ha`, `gridcell_area_ha`,
`final_grain_production_kg`, and `top_weight_production_kg` when crop area is
available. Downstream bioenergy comparison and LCA/TEA scripts can use these
columns directly.

Set `use_cropland_mask: false` to run all grid cells. If `cropland_strict:
false`, a missing raster only prints a warning and the pipeline falls back to
all-land mode.

#### Legacy helper path — build a cropland shapefile first

The older helper scripts are still useful if you want to inspect or reuse a
pre-built cropland-only point shapefile outside the main pipeline.

##### Create a binary cropland mask (`r_scripts/landcover_raster.R`)

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

##### Aggregate to grid points (`r_scripts/landcover_raster_to_gridpoints.R`)

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

##### Feed into Mode B

```yaml
use_existing_point_shapefile: true
existing_point_shapefile_path: "gridpoints/montana_cropland_5k_above_threshold.shp"
```

---

## Soil data sources

Set `soil_source` in `config.yml`.

| Value | Coverage | Method | Notes |
|-------|----------|--------|-------|
| `SSURGO` | United States only | Queries USDA Soil Data Access (SDA) web service per point | Most detailed US data; requires internet; queries one point at a time |
| `SOILGRIDS_10K` | Global | Reads a pre-downloaded master `.SOL` file at 10 km resolution | Fastest for large domains; download the country file once (see below) |
| `AGMIP` | Global cropland | Reads the AgMIP/Han DSSAT-ready country `.SOL` files at 5 arc-min (~10 km) | Semantic alias for the Han et al. global DSSAT soil profile database; use this when you want the source named explicitly in scenario IDs |
| `SOILGRIDS_ONLINE` | Global | Queries ISRIC SoilGrids 2.0 via REST API or VRT/GDAL virtual rasters | Flexible; set `soilgrids_mode: REST` (interactive) or `VRT` (HPC batch) in `config.yml` — see [VRT vs REST](#soilgrids-online-vrt-vs-rest-api) |
| `HWSD` | Global | Reads the FAO Harmonized World Soil Database v2.0 (raster of mapping-unit IDs + SQLite attribute DB) | The FAO "official" ~1 km global product; the long-standing reference for global gridded crop-model studies. **One-time manual download** from FAO (not a streaming API): set `hwsd_raster_file` + `hwsd_db_file` in `config.yml`. Samples the dominant soil per mapping unit and computes DSSAT physics. Requires `rasterio` (Python) / `terra`+`RSQLite` (R). |
| `HIHYDROSOIL` | Global | Reads local HiHydroSoil v2.0 hydraulic GeoTIFF/VRT rasters | Uses pF4.2, pF2/pF2.5, and saturated-water rasters directly for LL/DUL/SAT; set `hihydrosoil_raster_dir`. If sand/silt/clay rasters are absent it uses the product's USDA texture-class raster as a documented approximation. |
| `SLGA` | Australia | Reads local Soil and Landscape Grid of Australia rasters | Australian regional soil source; set `slga_raster_dir`. Computes DSSAT hydraulics from local sand/clay/silt/BD/OC rasters using the same Saxton-Rawls path as SoilGrids. |
| `WISE30SEC` | Global | Reads local WISE30sec GeoTIFF/VRT rasters | Global soil fallback using the same texture/BD/OC to Saxton-Rawls path as SLGA; set `wise30sec_raster_dir`. |
| `WOSIS` | Global point profiles | Reads a processed WoSIS layer CSV and assigns nearest profile | For calibration/validation or sparse measured-profile runs. Requires a harmonized CSV (`wosis_profile_csv`) with sand/clay/silt/BD/OC and layer depths. |

### AgMIP / SoilGrids 10K pre-formatted files

Pre-formatted DSSAT-ready `.SOL` files at 5 arc-min (~10 km) resolution, organised by country:

- **Download (Harvard Dataverse):** https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/1PEEY0
- **Citation:** Han, Ines, and Koo (2019). *Environmental Modelling & Software* 119:70–83. https://doi.org/10.1016/j.envsoft.2019.05.012

Download the country file you need (e.g. `US.SOL`), place it under `SoilGrids/`, and set either source name:
```yaml
soil_source: "AGMIP"          # or "SOILGRIDS_10K" for the legacy label
external_soil_file: "SoilGrids/US.SOL"
```

### HWSD v2.0 (FAO Harmonized World Soil Database)

`HWSD` reads two files you download **once** from FAO (it is a curated dataset, not a streaming API):

- The HWSD2 **raster** of mapping-unit (SMU) IDs (GeoTIFF/BIL).
- The HWSD2 **attribute database** (SQLite) with per-layer soil properties.
- **Source:** FAO HWSD v2.0 — https://www.fao.org/soils-portal/data-hub/soil-maps-and-databases/harmonized-world-soil-database-v2-0/

Place them anywhere and point `config.yml` at them:
```yaml
soil_source:       "HWSD"
hwsd_raster_file:  "HWSD/HWSD2.bil"       # or .tif
hwsd_db_file:      "HWSD/HWSD2.sqlite"
```
The module samples the raster at each point to get the SMU ID, selects the **dominant** soil component (largest area share) from the SQLite DB, computes DSSAT physics (Saxton & Rawls), and writes per-point `.SOL` files. Column names are matched tolerantly; on first use, confirm one output `.SOL` looks sane against your particular HWSD2 release. Points over no-data cells are skipped with a warning.

### SoilGrids online: VRT vs REST API

`SOILGRIDS_ONLINE` can fetch data two ways. Select with the **`soilgrids_mode`** key in `config.yml` (both pipelines honour it):

```yaml
soilgrids_mode: "REST"   # JSON REST API (default)
soilgrids_mode: "VRT"    # GDAL virtual rasters
```

| Mode | Best for | How it works | Dependencies |
|------|----------|--------------|--------------|
| `REST` | Interactive / small local runs | One HTTPS request **per point**; rate-limited (~5 req/min) with exponential back-off (up to 5 attempts) | none beyond `requests` / `httr` |
| `VRT` | HPC / large batch jobs | Streams each global SoilGrids raster **once** (via `/vsicurl/`) and samples all points from it — no per-point HTTP, no rate limit, and often better coverage at the cell level | Python: `rasterio` (GDAL ≥ 3); R: `terra` |

Both modes project points to SoilGrids' native Interrupted Goode's Homolosine CRS, compute DSSAT physics with Saxton & Rawls (2006), and write identical `.SOL` output (the file header records which mode produced it). Points that fall on a SoilGrids no-data cell (e.g. over water) are masked and skipped with a warning rather than producing invalid soil. For grids beyond a handful of points, prefer `VRT`: it reads 36 rasters (6 properties × 6 depths) total regardless of point count, whereas `REST` makes one rate-limited request per point.

Direct package calls use the merged `dssatutils` `config.yml` default
(`soil.soilgrids_online.use_rest_api`; `false` = VRT, `true` = REST). The
pipeline-level `soilgrids_mode` key remains the recommended control for normal
tutorial runs.

---

## Weather data sources

Use `weather_source: "NASA_POWER_CHIRPS_V3"` for NASA POWER with CHIRPS v3
daily rainfall. v3 options are `chirps_v3_product` (`rnl` gauge-adjusted final
or `sat` satellite-only), `chirps_v3_stream` (`final` or `prelim`), and
`chirps_v3_fetch_mode` (`monthly_netcdf` recommended because yearly v3 files are
large).

Set `weather_source` in `config.yml`.

| Value | Coverage | Approx. resolution | Notes |
|-------|----------|--------------------|-------|
| `DAYMET` | North America | ~1 km daily | Best choice for US simulations; uses the `daymetr` package |
| `NASA_POWER` | Global | ~0.5° daily | Recommended for international runs; outputs real wind speed (`WS2M`) at `WNDHT = 2.0` m |
| `GRIDMET` | Contiguous US | ~4 km daily | High spatial resolution for US; requires `terra`, `ncdf4`, and `httr` |
| `OPEN_METEO` | Global (1940–present) | ~9–11 km temperature/humidity plus ERA5 forcing (ERA5-Seamless) | Keyless, no registration; higher-resolution alternative to NASA-POWER for Europe / Asia / Africa / Oceania / South America. Writes API daily mean dewpoint and relative humidity to `TDEW` / `RH2M`; ERA5 supplies radiation, precipitation, and wind. Requires `httr` + `jsonlite` (R) / `requests` (Python). Data is CC-BY 4.0 (Copernicus/ECMWF) — cite when publishing. |
| `NASA_POWER_CHIRPS` | Global, **rainfall 50°S–50°N** | NASA POWER ~0.5° + CHIRPS rain ~0.05° (~5.5 km) | **Hybrid**: all variables from NASA POWER, but rainfall replaced with high-resolution, station-blended CHIRPS — markedly better precipitation for the tropics / semi-arid (Africa, India). Outside 50°S–50°N (and over CHIRPS no-data) it falls back to NASA POWER rain, so output stays global. Keyless. Downloads CHIRPS yearly netCDF (cached); resolution via `chirps_resolution` (`p05` default / `p25`). Requires `xarray`+`netCDF4` (Python) / `terra` (R). Cite Funk et al. 2015. |
| `AGERA5` | Global (1979–present), incl. poles | ~0.1° (~10 km) daily | ECMWF agrometeorological reanalysis (ERA5 reprocessed for agriculture); all variables, higher-res than NASA POWER, covers high latitudes (unlike CHIRPS). The default `timeseries` backend caches compact all-variable CSVs on globally anchored AgERA5 cells; `gridded` retains the legacy daily-NetCDF ZIP path. **Not keyless**: run `setup_cds_credentials()` or provide `CDSAPI_KEY` / `~/.cdsapirc`, then accept the dataset licence once (see below). Requires `cdsapi` (Python) / `ecmwfr` (R). |
| `CHELSA_W5E5` | Global (1979–2016) | 30 arcsec (~1 km) daily | Topographically downscaled daily climate forcing. Requires local NetCDFs (`chelsa_nc_dir`) with `tasmax`, `tasmin`, `pr`, and `rsds`. |
| `AGMERRA` / `AGCFSR` | Global (1980–2010) | 0.25° daily | AgMIP-standard historical climate forcing for benchmark/intercomparison studies. Requires local NetCDFs (`agmerra_nc_dir` / `agcfsr_nc_dir`). |
| `SILO` | Australia (1889–present) | Australian gridded daily | Australian regional weather source. Requires local SILO NetCDFs via `silo_nc_dir`. |
| `PRISM` | Contiguous US (daily 1981–present) | 4 km daily | Downloads/caches PRISM public 4 km daily grids (`prism_cache_dir`) for precipitation, Tmax, Tmin, and dewpoint. SRAD/wind/RH are written missing (`-99`). |
| `MSWX` | Global | ~0.1° daily | Full local NetCDF weather forcing; set `mswx_nc_dir`. |
| `MSWEP` | Global precipitation | ~0.1° daily | Hybrid source: NASA POWER supplies non-rain variables, local MSWEP NetCDFs replace rainfall; set `mswep_nc_dir`. |
| `CRUJRA` | Global | 0.5° daily | Local CRU-JRA NetCDF weather forcing; set `crujra_nc_dir`. |
| `TERRACLIMATE` | Global | ~4 km monthly | Monthly product written only as screening/climatology input; set `terraclimate_nc_dir`. Do not treat as true daily weather. |

Between them, `NASA_POWER`, `OPEN_METEO`, `NASA_POWER_CHIRPS`, and `AGERA5` give full global daily coverage, so Europe, Asia, and Africa runs need no US-only source. For **rainfed crops in the tropics/semi-arid**, `NASA_POWER_CHIRPS` gives the best rainfall; for **all-variable global incl. high latitudes**, `AGERA5` (key required) is the highest-resolution option. Soil coverage is likewise global via `SOILGRIDS_10K` / `SOILGRIDS_ONLINE` (ISRIC SoilGrids) or `HWSD` (FAO).

> **NASA-POWER wind note:** Previous versions of `weather_nasapower.R` wrote `WIND = -99` (missing) and set `WNDHT = -99` in the `.WTH` header. The current version writes the real `WS2M` value and correctly sets `REFHT = 2.0` and `WNDHT = 2.0`, matching the NASA-POWER AG community product specification.

> **GridMET variables (Wind, Humidity, PET) note:** Standard DSSAT runs require only solar radiation, temperature, and precipitation. However, if using advanced evapotranspiration (such as FAO-56 Penman-Monteith, `MEEVP` = 'F' or 'G' in DSSAT simulation controls), the model consumes daily wind speed (`WIND` in `.WTH` files) and relative humidity/dewpoint (`TDEW`). GridMET's wind speed (`vs`) and specific humidity (`sph`) can be converted to dewpoint and wind speed to make Penman-Monteith crop-water simulations highly accurate, bypassing empirical temperature-based approximations (e.g. `Tmin - 2.5` for dewpoint). The GridMET reference evapotranspiration (`pet`) is not directly consumed by standard DSSAT configurations since DSSAT simulates crop transpiration and potential/actual ET dynamically based on crop leaf area index, canopy cover, and water balance.

### Copernicus CDS one-time setup

`AGERA5`, `ERA5_LAND`, and E-OBS CDS mode need a free Copernicus CDS account and
**three** one-time steps:

1. **Register** at https://cds.climate.copernicus.eu/ and copy your *Personal Access Token* from your profile.
2. **Run the dssatutils setup helper**, or set `CDSAPI_KEY` / `CDSAPI_URL`
   yourself:
   ```r
   library(dssatutils)
   setup_cds_credentials()
   ```
   The helper writes `~/.cdsapirc` (Linux/macOS) or `%USERPROFILE%\.cdsapirc`
   (Windows) and configures `ecmwfr` for R.
3. **Accept the dataset licence** (otherwise the first request fails with `403 … required licences not accepted`): open
   https://cds.climate.copernicus.eu/datasets/sis-agrometeorological-indicators?tab=download and accept the licence(s) under *Manage licences*.

Then install the CDS client dependency (`pip install "dssatutils[cds] ..."` for
Python, `install.packages("ecmwfr")` for R) and set the desired weather source.
Requests are queued server-side, so the first run for a new area/period can take
minutes; downloads are cached under the source cache directory.

AgERA5 cache behavior is configured in this repository's `config.yml`:

```yaml
agera5_backend: "timeseries"         # recommended; use "gridded" only for the legacy ZIP workflow
agera5_data_format: "csv"
agera5_timeseries_chunk_degrees: 0.1  # one globally reusable AgERA5 cell per cache key
agera5_max_concurrent_requests: 4
```

`0.1` minimizes transfer and disk use for sparse crop-model points. A larger
chunk groups neighboring cells into fixed globally aligned tiles, reducing the
number of CDS requests while downloading more cells. Both the R and Python
pipelines use the shared cache at `agera5_netcdf_cache/` under the engine input
root.

---

## Sensitivity experiments: sweep weather × soil combinations

`config.yml` runs **one** `weather_source` × `soil_source` at a time. To study how
much the choice of input data drives your results, you often want to run the *same*
simulation under **every** combination of weather and soil and compare. Two files
on top of the pipeline do this — you never edit the pipeline itself:

| File | Role |
|------|------|
| **`experiment.yml`** | Single source of truth for the study: the constant *base* simulation (region, crop, template, grid) plus the *factor lists* to sweep — `weather_source`, `soil_source`, and an optional `period` (weather-year window). |
| **`run_experiment.R`** | Orchestrator: builds the full factorial, runs `dssat_main_pipeline.R` once per combination (optionally in parallel), then aggregates everything into one tidy table, plus a summary, ANOVA variance decomposition, and boxplot **per response variable**. |

### Run the whole workflow in one go (`run_all.R`)

The stages below (sweep → analysis → validation) can be chained with a single
command:

```bash
Rscript run_all.R                       # sweep (experiment.yml) + analysis
Rscript run_all.R --experiment my.yml   # use a different experiment file
Rscript run_all.R --validate            # also run the observed-data validation
Rscript run_all.R --validate-only       # just the validation stage
Rscript run_all.R --no-analysis         # sweep only
```

`run_all.R` reads the experiment file to locate the sweep's `combined.csv` and
response variable(s), runs `analyze_experiment.R` on each, and (with
`--validate`) finishes with `validate_against_observed.R`. It stops at the first
failing stage. The individual scripts below are still there if you want to run a
single stage or re-run just the analysis on an existing sweep.

### How it works

Each combination needs its own settings without disturbing your `config.yml`.
The loaders (`config_loader.R` / `config_loader.py`) honour a **`DSSAT_CONFIG_FILE`**
environment variable that, when set, overrides the usual `config.yml` lookup. So
for every combination the orchestrator:

1. merges settings — `config.yml < experiment.yml:base < the combination's factor values`,
2. writes that to a **private temp config file** under `.experiment_configs/`,
3. runs the pipeline as a fresh `Rscript` subprocess with `DSSAT_CONFIG_FILE`
   pointing at it.

Your real `config.yml` is **never modified**, and parallel workers never share a
config — so nothing can clobber anything. Each combination's results land in its
own `results/<scenario>_results.csv`, tagged on read with the factor values.

### Define the experiment (`experiment.yml`)

```yaml
experiment_name: "wx_soil_sensitivity"

base:                         # held constant across every combination (config.yml keys)
  project_name:        "dssat_spatial_demo"
  grid_spacing_meters: 40000
  crop_extension:      "MZ"
  state_name_filter:   ["Iowa"]    # region of interest
  weather_start_year:  1982        # must cover the template's planting year (UFGA8201 = 1982)
  weather_end_year:    1983
  template_file_name:  "UFGA8201.MZX"
  treatment_start:     1
  treatment_end:       4
  treatment_list:      []            # optional explicit IDs, e.g. [1, 5, 10]

factors:                      # the full factorial of these lists is run
  weather_source: ["DAYMET", "OPEN_METEO"]   # keyless, both cover 1982
  soil_source:    ["SSURGO", "SOILGRIDS_10K"]
  # period:       ["1982-1983", "1984-1985"]   # OPTIONAL 3rd factor: weather-year windows

exclude: []                   # optional: drop specific (weather, soil) pairs

options:
  stop_on_error:  false       # continue past a failed combination
  reuse_existing: true        # skip combos already computed (resumable)
  dry_run:        false       # true = print the plan and exit, run nothing
  validate:       true        # pre-flight checks: drop impossible combos, warn on risky ones
  max_parallel:   1           # >1 runs combinations concurrently (see below)
  response_vars:              # one or many; each gets its own summary/variance/plot
    - "final_grain_kg_ha"
    # - "soil_organic_carbon_delta_kg_C_ha"
    # - "nitrate_leaching_kg_ha"
```

**Period as a third factor.** Inter-annual weather variability often rivals the
choice of dataset. Adding `period:` reruns each weather × soil pair over several
year windows so the ANOVA can separate "*which dataset*" from "*which years*".
Note that the pipeline does **not** rebase the template's planting year, so each
window must cover the planting year in your FileX — the `UFGA8201` demo template
is fixed to 1982. Use a multi-year experiment file (or per-window templates) to
sweep periods meaningfully.

**Multiple response variables.** Inputs that barely move yield can still strongly
move the water or nitrogen balance, so `response_vars` accepts a list — each is
summarised, decomposed, and plotted separately (`response_var` still works for a
single one). Any numeric column of the results CSV is valid (e.g.
`final_grain_kg_ha`, `soil_organic_carbon_delta_kg_C_ha`, `nitrate_leaching_kg_ha`,
`final_irrigation_amount_mm`).

### Run it

```bash
Rscript run_experiment.R                     # uses experiment.yml
Rscript run_experiment.R my_other.yml        # or any experiment file
```

(or open `run_experiment.R` in RStudio and click **Source**.) Tip: set
`dry_run: true` first to preview the combinations, validation results, and output
paths without running DSSAT.

### Pre-flight validation (`validate: true`)

Before running, the orchestrator screens every combination and prints the plan:

- **Drops impossible combinations** — `AGERA5` without a configured Copernicus
  CDS token (`setup_cds_credentials()`, `CDSAPI_KEY`, or `~/.cdsapirc`), or
  `NASA_POWER` with a period starting before 1984.
- **Warns on risky ones** — `GRIDMET`/`SSURGO` are US-only; `SOILGRIDS_ONLINE`
  (REST) is rate-limited and may throttle under high `max_parallel`.

### Parallel execution (`max_parallel > 1`)

Combinations are independent, so they can run concurrently
(`parallel::mclapply`; Windows falls back to serial). Because combinations that
share a weather or soil source also share the framework's download caches, the
orchestrator first runs a short **serial warm-up pass** that touches each
weather/soil source once to populate those caches, then runs the rest in
parallel so workers only *read* them — avoiding cold-download races. Parallel (or
period-swept) runs also give each combination a unique `project_name` so their
weather/soil/run folders never collide.

### Outputs

| Output | Contents |
|--------|----------|
| `results/EXPERIMENT_<name>_combined.csv` | Every combination's point-level results stacked into one tidy long table, tagged with `weather_source`, `soil_source`, `period`, `scenario_id`. |
| `results/EXPERIMENT_<name>[_<var>]_summary.csv` | Mean / sd / min / max of the response variable per combination. |
| `results/EXPERIMENT_<name>[_<var>]_variance.csv` | ANOVA share of variance attributable to each factor. |
| `results/EXPERIMENT_<name>[_<var>]_boxplot.png` | Response distribution, weather on the x-axis, soil as colour (faceted by period when swept). |
| `results/experiment_logs/<scenario>.log` | Full pipeline stdout/stderr for each combination (where to look if one fails). |

(The `_<var>` infix appears only when you request more than one response variable.)
The `*_variance.csv` / console output is a **Type-I ANOVA variance decomposition** —
the share of variance in the response attributable to each factor — a quick read
on which input matters more, e.g.:

```
share of variance in final_grain_kg_ha explained:
  weather_source  81.30%
  soil_source      7.00%
  residual        11.70%
```

> **Notes.** The default factors (`DAYMET`, `OPEN_METEO` × `SSURGO`,
> `SOILGRIDS_10K`) are all keyless and valid for the 1982 US demo, so the
> experiment runs as-is. With `reuse_existing: true` the sweep is fully
> resumable — rerun after a crash and it picks up where it left off.

### Deeper analysis of a sweep (`analyze_experiment.R`)

`run_experiment.R` emits a boxplot + variance CSV; for the richer views, run the
companion analyser on any sweep's `combined.csv`:

```bash
Rscript analyze_experiment.R results/EXPERIMENT_<name>_combined.csv
Rscript analyze_experiment.R <combined.csv> final_grain_kg_ha --treatment 1 --out analysis/run1
Rscript analyze_experiment.R <combined.csv> --boundary shapefile/world.shp   # custom outline
```

The maps draw a **boundary outline underneath the points** (state/country
borders, via `coord_sf`) — by default the bundled US states
(`shapefile/tl_2024_us_state.shp`), cropped to your grid. Pass `--boundary
<file.shp>` for another region (e.g. a world/country layer for non-US studies),
or `--boundary none` to disable it.

It averages each grid point over treatments/years per combination (or restrict to
one treatment with `--treatment`) and writes, under `analysis/<name>/`:

| Output | What it shows |
|--------|---------------|
| `fig_yield_maps.png` | Small-multiple maps, one panel per combination — eyeball where combos agree/diverge. |
| `fig_sensitivity_cv_map.png` + `sensitivity_by_point.csv` | Per-point **coefficient of variation across combinations** — *where* the input choice moves the result most. |
| `fig_variance_decomposition.png` + `variance_decomposition.csv` | Share of variance attributable to weather / soil / period / residual. |
| `fig_pairwise_rmsd.png` + `pairwise_rmsd.csv` | RMSD between every pair of combinations — how interchangeable two inputs are. |
| `fig_rank_stability.png` + `rank_stability.csv` | Spearman correlation of the spatial yield pattern between combinations ("does good land stay good?"). |
| `summary_by_combo.csv` | mean / sd / CV / range per combination. |

The two complementary signals: the **ANOVA** shares how much of the *total* spread
each factor explains (spatial variation usually dominates the residual), while the
**per-point CV** isolates input sensitivity after removing the spatial gradient.

---

## Validation against observed data

Many DSSAT example experiments ship with **measured** end-of-season data (a
"FileA", e.g. `DSSAT48/Maize/UFGA8201.MZA` holds observed grain yield `HWAM` per
treatment). `validate_against_observed.R` uses them to benchmark the pipeline by
comparing three pathways of yield **at each experiment's own location**:

| Pathway | Inputs | Meaning |
|---------|--------|---------|
| **observed** | — | Measured `HWAM` per treatment, read from the FileA. |
| **original** | the experiment's *local* `.WTH` weather + its soil | The calibrated "textbook" DSSAT run (run via the `DSSAT` R package). |
| **pipeline** | this pipeline's *gridded* weather × soil sweep | The same experiment, but with downloaded inputs — so you also see *which gridded input best reproduces the observed / locally-simulated yield*. |

```bash
Rscript validate_against_observed.R               # all configured experiments
Rscript validate_against_observed.R UFGA8201      # one or more by name
Rscript validate_against_observed.R --dry-run     # resolve metadata only, run nothing
```

For each experiment the script parses the FileX (treatments, planting year,
weather station, soil id, **coordinates from the weather-file header**), then:
runs the original experiment with its native inputs; runs the pipeline at the
experiment's coordinates via a single-point MODE B shapefile, **delegating the
weather × soil sweep to `run_experiment.R`**; and reads the observed FileA.

It is **location- and coverage-aware**. With `ALL_SOURCES <- TRUE` (the default,
in the config block) every site is run against *all* gridded sources that are
geographically and temporally feasible for it — e.g. a contiguous-US 2002 site
gets `DAYMET`+`GRIDMET`+`OPEN_METEO`+`NASA_POWER`+`NASA_POWER_CHIRPS`(+`AGERA5`
if a CDS key is present) × `SSURGO`+`SOILGRIDS_10K`+`SOILGRIDS_ONLINE`, while a
global pre-1984 site drops the US-only and post-1984-only sources automatically.
Set `ALL_SOURCES <- FALSE` for a single recommended pair per site (far fewer
runs). Combinations that can't run (missing key/files, transient source errors)
are skipped, and whatever succeeds is reported.

**Preflight your sources first.** Before launching the long sweep, check which
sources are actually runnable on your machine (R packages installed, keys/files
present):

```bash
Rscript validate_against_observed.R --check
```

```
=== Source preflight (this machine) ===
 source            runnable note
 AGERA5            NO       install R pkg(s): ecmwfr; then run setup_cds_credentials()
 DAYMET            yes
 GRIDMET           yes
 ...
```

With `PRUNE_UNRUNNABLE <- TRUE` (default) un-runnable sources are dropped from
every site's matrix automatically, so you never waste runs on guaranteed
failures. The check is deterministic (deps/keys/files) — it can't predict a
transient API timeout, but it catches the common "source X never works here"
case up front. Results are also written to `validation/source_preflight.csv`.

**Outputs** (under `validation/`): `validation_long.csv` (tidy obs/original/
pipeline yields per treatment), `validation_metrics.csv` (RMSE, nRMSE, mean bias,
Willmott's *d*, modelling efficiency `EF`, `R²` of every simulated source vs.
observed), and figures `fig_obs_vs_sim.png` (observed-vs-simulated scatter with a
1:1 line) and `fig_metrics_rmse.png` (error by input source).

A typical reading: the **original** local run tracks observed best; among gridded
inputs the closest combination quantifies how much accuracy you trade for global
coverage, and a poor soil/weather source stands out immediately (negative `EF`).

> **Notes.** Configure paths / crop / experiment list in the block at the top of
> the script (defaults: maize, `../DSSAT48`). Raw DSSAT FileX templates are
> auto-adapted with the pipeline's placeholders (`SOIL_ID`, `00000000`); the
> bundled demo template is reused untouched. The full 10-site run with global
> SoilGrids/Open-Meteo fetches is long but resumable (`reuse_existing`).

---

## Repository layout

```text
.
├── config.yml                         # CENTRAL CONFIG — shared by R and Python pipelines
├── config_loader.R                    # loads config.yml into the R pipeline
├── config_loader.py                   # loads config.yml into the Python pipeline
├── dssat_main_pipeline.R              # MAIN entrypoint (R) — start here
├── dssat_main_pipeline.py             # Python port of the R pipeline (serial/multiprocessing)
├── experiment.yml                     # weather × soil sensitivity experiment definition
├── run_all.R                          # one-command runner: sweep → analysis → (validation)
├── run_experiment.R                   # orchestrator: sweeps combinations, aggregates results
├── analyze_experiment.R               # sweep analysis: maps, variance, pairwise/stability tables
├── validate_against_observed.R        # observed vs original-DSSAT vs pipeline yield validation
├── r_scripts/
│   ├── soil_ssurgo.R
│   ├── soil_soilgrids.R
│   ├── soil_soilgrids_online.R
│   ├── soil_hwsd.R                    # FAO HWSD v2.0 global soil (external files)
│   ├── weather_daymet.R
│   ├── weather_nasapower.R
│   ├── weather_gridmet.R
│   ├── weather_openmeteo.R            # global, keyless ERA5 (EU/Asia/Africa/…)
│   ├── weather_nasapower_chirps.R     # NASA POWER + high-res CHIRPS rainfall
│   ├── weather_agera5.R               # AgERA5 (CDS key required)
│   ├── landcover_raster.R
│   └── landcover_raster_to_gridpoints.R
├── python_scripts/                    # Python twins of the r_scripts/ modules
├── dssat_templates/                   # template FileX/SQX + cultivar/ecotype/species files
├── shapefile/                         # boundary polygons (TIGER/Line or custom)
├── gridpoints/                        # generated or user-supplied point shapefiles
├── data/landcover/                    # landcover rasters and derived masks (Mode C)
├── soil/                              # downloaded/cached soil products
├── weather/                           # downloaded/cached weather products
├── SoilGrids/                         # pre-downloaded master .SOL files (SOILGRIDS_10K)
├── dssat_runs/                        # generated per-point DSSAT run folders
├── results/                           # merged CSVs and yield maps
├── tests/
│   ├── test_smoke.py                  # offline CI test (config + imports + .WTH writer)
│   ├── test_cross_language_parity.py  # R/Python config and dependency contracts
│   ├── test_dssat_model_e2e.py        # opt-in real one-point DSSAT execution
│   ├── test_dssat_model_e2e.R         # R twin of the real one-point check
│   ├── test_global_sources.py         # live EU/Asia/Africa test (Open-Meteo + SoilGrids)
│   ├── test_e2e_comprehensive.py      # live Python provider matrix
│   ├── test_e2e.R                     # compact live R provider smoke test
│   └── test_e2e_comprehensive.R       # live R provider matrix
├── .github/workflows/smoke.yml        # CI: cross-platform offline checks
├── .github/workflows/e2e.yml          # CI: live Python/R provider checks
├── environment.yml / requirements.txt # Python dependency pins (conda / pip)
├── setup_renv.R                       # R dependency setup (renv)
└── hpc/
    ├── dssat_mpi_runner.py            # Advanced: MPI runner for HPC
    └── run_dssat_python.slurm         # Advanced: SLURM submit script
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

> **MPI runner note:** `hpc/dssat_mpi_runner.py` supports both run modes via `--run_mode`. In `sequence` mode it runs the per-point `*.SQX` with the `Q` switch; in `experiment` mode it runs the seasonal FileX (`*.MZX`, `*.WHX`, etc.) with the `A`/`B` switch. Match `--run_mode` (and `--trt_*` / `--seq_*`) to the template you built the folders with.

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

The run-name keys control the directory under `dssat_runs/`:

```yaml
run_tag: "run1"          # appended to the generated base name
run_name_style: "grid"   # grid | scenario
run_name_override: ""    # non-empty value is used verbatim
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

### Central configuration: `config.yml` (recommended)

The repository ships a single **`config.yml`** at the project root that is the **shared source of truth for both pipelines** — the R pipeline (`dssat_main_pipeline.R`) and the Python port (`dssat_main_pipeline.py`) each load it at startup via a small loader (`config_loader.R` / `config_loader.py`). Edit settings **once** in `config.yml` and both pipelines pick them up, so the R and Python configurations cannot drift apart.

```yaml
# config.yml — edit these; both R and Python read them.
project_name:         "dssat_spatial_demo"
grid_spacing_meters:  50000          # 50 km demo; 5000–10000 for production
crop_extension:       "MZ"           # MZ=maize, WH=wheat, SB=soybean, ...

weather_source:       "DAYMET"       # DAYMET | NASA_POWER | GRIDMET | OPEN_METEO | NASA_POWER_CHIRPS | NASA_POWER_CHIRPS_V3
weather_start_year:   1982
weather_end_year:     1983

soil_source:          "SOILGRIDS_10K"  # SSURGO | SOILGRIDS_10K | SOILGRIDS_ONLINE
external_soil_file:   ""                # blank = script default (SoilGrids/US.SOL)

template_file_name:   "UFGA8201.MZX"
run_mode:             "experiment"   # experiment | sequence
treatment_start:      1
treatment_end:        4
treatment_list:       []             # optional explicit non-contiguous IDs
```

**How precedence works:** the repository `config.yml` is the complete shared
default for R and Python. Set `DSSAT_CONFIG_FILE` to a study-specific YAML file
to merge only that study's overrides over the repository defaults. `DSSAT_EXE`
and `DSSAT_BASE` environment variables have highest priority for machine-local
paths. Missing or malformed configuration fails immediately instead of silently
running with hidden defaults. PyYAML (Python) or `yaml` (R) is therefore a
required dependency.

Do not edit SECTION 0 in either pipeline for normal customization. It contains
configuration plumbing and derived values; all user-facing knobs belong in
`config.yml`. The HPC MPI runner remains configured by command-line flags (see
[HPC execution](#hpc-execution-mpi--slurm)).

### Common first customisations (after the demo works)

- **Spatial scope:** change `grid_spacing_meters`, `state_name_filter`, or set `use_existing_point_shapefile: true`.
- **Weather:** change `weather_source`, `weather_start_year`, and `weather_end_year`.
- **Soil:** change `soil_source`; for `SOILGRIDS_10K`, set `external_soil_file`.
- **Crop:** update `crop_extension` and `template_file_name` together.
- **HPC prep:** set `run_dssat_execution: false` and `zip_for_hpc: true`.

### R/Python parity note

Both pipelines consume the same merged YAML and run the same engine contract.
One internal optimization remains implementation-specific: Python can opt into
symlinked support files with `use_symlinks`, while R currently copies them.
Download validation/retry and per-point DSSAT failure isolation are matched in
both drivers and do not change scientific output schemas. Keep
`use_symlinks: false` for identical, portable run-folder contents.

### Running non-interactively (command line)

```bash
Rscript dssat_main_pipeline.R
```

---

## Validate one point folder before scaling

Before running hundreds or thousands of points, validate a **tiny run** end-to-end. This saves significant time.

### Recommended validation flow

1. Reduce to 5–20 points: increase `grid_spacing_meters`, tighten `state_name_filter`, or use a small existing shapefile.
2. Keep the run minimal: set equal `treatment_start` / `treatment_end` and use 2–3 weather years.
3. Set `cleanup_run_folders: false` so you can inspect DSSAT logs afterwards.
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
- Set `dssat_cores: 1` to make debugging simpler.
- Confirm your template is configured to write **CSV outputs** (`summary.csv`), since the parser expects them.

---

## Testing

The repository ships five Python and three R test scripts under `tests/`, plus
separate smoke and live-provider GitHub Actions workflows.

### Offline smoke test (no internet, no DSSAT)

`tests/test_smoke.py` validates the parts of the pipeline that don't need live
APIs or a DSSAT install: `config.yml` loading and fallback semantics, helper
imports, MPI CLI help without a cluster runtime, and deterministic `.WTH`
formatting. `tests/test_cross_language_parity.py` checks central configuration
and selected R/Python contracts. GitHub checks out the exact dependency source
revisions required by those parity assertions; they may skip in a standalone
clone that has no sibling or `.ci-deps` source checkout.

```bash
# from the repo root, in your Python environment
python -m pytest -m "not integration and not slow" -v
#   or as a plain script:
python tests/test_smoke.py
```

The provider-heavy comprehensive module is deliberately marked `integration`.
Run it explicitly when you want real downloads (several minutes and potentially
large cached NetCDFs). It requires every expected output and retries only
recognised transient HTTP/transport failures. Set
`DSSAT_USE_WORKSPACE_SIBLINGS=1` only when intentionally testing un-released
sibling `dssatutils` source instead of the installed pinned package:

```bash
python -m pytest tests/test_e2e_comprehensive.py -m integration -v -s
DSSAT_RUN_LIVE_E2E=1 Rscript --vanilla tests/test_e2e_comprehensive.R
```

### Live global-source test (internet required, no DSSAT)

`tests/test_global_sources.py` exercises the **global** data sources end-to-end on three points (Europe / Asia / Africa): it downloads Open-Meteo weather and SoilGrids soil and checks that the emitted `.WTH` / `.SOL` files are well-formed. It accounts for every input point individually — a point that hits a genuine SoilGrids coverage gap is reported as `[skip-no-data]` rather than silently dropped. Because it hits live, rate-limited APIs it is meant to be run on your own machine, not in CI.

```bash
python tests/test_global_sources.py
python tests/test_global_sources.py --keep   # keep the output dir for inspection
```

### Continuous integration

`.github/workflows/smoke.yml` runs the offline Python tests, byte-compiles the
Python modules, and validates R configuration/syntax on Ubuntu, macOS, and
Windows. This catches platform-specific path/config regressions on every push
and pull request. The workflow installs the R `yaml` dependency with
`Rscript --vanilla`, matching the cross-language parity subprocess rather than
the repository's `renv`-activated interactive environment.

`.github/workflows/e2e.yml` runs live Python and R provider checks on Ubuntu on
pushes, pull requests, a weekly schedule, and manual dispatch. Live provider
outages can therefore produce explicit skips after bounded retries; test-code,
schema, and missing-output defects still fail. JUnit and R logs are retained
even after failures. The R job disables project `renv` auto-activation and uses
its explicitly installed integration-test packages, preventing misleading
partial-lock warnings in the main process and parallel workers.

These provider checks do **not** execute the licensed DSSAT model. A true model
validation requires an installed DSSAT executable. The committed one-point
fixture provides an explicit real-model check:

```bash
DSSAT_EXE=/verified/path/to/dscsm048 \
  python -m pytest tests/test_dssat_model_e2e.py -v
Rscript --vanilla tests/test_dssat_model_e2e.R
```

It skips when DSSAT or its adjacent `DSSATPRO` support profile is absent. The
one-point validation flow above and MPI mini run below remain the production
checks. Do not describe the provider matrix as proof that DSSAT itself runs on
every hosted operating system.

---

## Option A — Run locally

### Option A1: Run from R (default)

With `run_dssat_execution: true` in `config.yml`, the pipeline runs DSSAT per point in parallel, parses outputs, and writes the merged results CSV and yield map automatically.

### Option A2: Run the MPI runner locally

This is the best way to test the HPC execution path on a workstation before submitting to a cluster:

```bash
mpirun -n 4 python hpc/dssat_mpi_runner.py \
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

Set `run_dssat_execution: false` and `zip_for_hpc: true` in `config.yml`, then:

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

Edit `hpc/run_dssat_python.slurm` to set `DSSAT_EXE`, the `RUN_DIRS` array (one or more `--base_dir` targets under `/scratch/$USER/`), `SUMMARY_DIR`, node count, and walltime. Supply the conda environment name or full prefix through `CONDA_ENV`, then:

```bash
CONDA_ENV=dssat_env sbatch hpc/run_dssat_python.slurm
squeue -u $USER
```

> Keep the SLURM script **generic**: avoid committing personal email addresses, allocation keys, or user-specific paths. Use placeholders or environment variables.

### 5 — End-to-end mini run (strongly recommended before scaling)

```bash
mpirun -n 2 python hpc/dssat_mpi_runner.py \
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

1. Rank 0 scans folders under `--base_dir` and identifies point folders (digit-only names), then broadcasts the task list.
2. Work is distributed by **dynamic work-stealing**: rank 0 acts as a manager and hands out one folder at a time to each worker on demand, so slow points can't stall a rank's static share while other ranks sit idle. (With `-n 1` it falls back to a single-process loop.)
3. Each worker runs DSSAT per point/treatment and **streams** rows straight to its `temp_results_rank_<rank>.csv` as each point finishes — no rank-wide in-memory accumulation, and progress survives a late rank failure.
4. Part files are written to **node-local scratch** (`--scratch_dir`, defaults to `$TMPDIR` in the SLURM script) to avoid a small-file metadata storm on shared Lustre/GPFS, then each rank stages its own part back to the shared summary dir in parallel.
5. Final merge is controlled by `--merge_mode`: `concat` (default) streams all parts into one `results_<RUN_NAME>.csv`; `none` leaves the per-rank parts in place (each a complete, valid CSV) to skip the rank-0 serial read on very large grids.

The runner also supports optional output cleanup (`--cleanup_mode never/success/always`) and output archiving (`--archive_outputs`).

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
| NASA POWER reports a missing weather day/value | The provider has a missing forcing value for that cell/date | Keep `check_weather_downloads: true`; invalid cached files are removed before each genuine retry. Before DSSAT starts, core forcing (`SRAD`, `TMAX`, `TMIN`, `RAIN`) is checked strictly; if the provider remains incomplete, that point is logged and omitted without aborting the others. Enable the optional short-gap repair only when interpolation is appropriate for the study. |
| Weather download stalls or rate-limit errors | Too many parallel requests | Set `weather_cores: 1`; SoilGrids REST back-off retry handles 429 errors automatically |
| Soil processing skips many points | Critical data (sand/clay/BD) is NA | Check `soil_processing_errors.log` in the soil output folder; consider switching soil source |
| Soil ID mismatch in DSSAT log | Placeholder `SOIL_ID` not replaced | Check `template_soil_id_placeholder` in `config.yml`; confirm the placeholder string is in your template |
| Weather does not cover full simulation period | Simulation years extend beyond downloaded range | For sequence runs: a 10-year rotation needs at least 10 + `NYERS` years of weather; extend `WEATHER_END_YEAR` or use the weather extension functions |
| `mpi4py` import error on HPC | Built against wrong MPI library | Rebuild: `pip install --no-cache-dir --no-binary :all: mpi4py` after loading the correct MPI module |
| Results CSV is empty after run | Parsing failed; DSSAT wrote `.OUT` only | Check per-point `LUN.LST`; enable CSV output in the `*OUTPUTS` section of your template |
| Custom DSSAT installation not found | Default assumptions (e.g. `C:/DSSAT48`) do not match your setup | For the gridded pipeline, specify `DSSAT_EXE` and `standard_dssat_dir` in `config.yml` / `config.yaml` or set the `DSSAT_EXE` environment variable |
| Shared package clone failures | Authentication mismatch on `dssatutils` / `dssatengine`, missing access, or GitHub rate limits | Verify that Git Credential Manager, a GitHub Personal Access Token (PAT) via `gitcreds` / `GITHUB_PAT`, or SSH keys are configured for the install path you are using |

---

## Performance tips

- **Tiny run folders (default):** when `DSSATPRO.V48` sits next to the DSSAT executable, the pipeline lets DSSAT resolve genotype/species/SDA/CO₂ files from the install directory instead of copying ~27 files into every point folder. Each folder then holds just the 4–5 essential files (`.WTH`, `SOIL.SOL`, FileX, `DSSBatch.V48`, `DSSATPRO.V48`). On a 10,000-point grid that's ~270,000 fewer files — a large saving on shared filesystems (Lustre/GPFS/NFS) where metadata ops dominate. Set `bundle_genotype_files: true` (auto-forced by `ZIP_FOR_HPC`) only when you need self-contained folders to ship to a host whose `DSSATPRO` doesn't match.
- **Weather downloads are parallel/cached:** per-point sources (Daymet, NASA POWER, Open-Meteo) download across `WEATHER_CORES`; gridded sources (GridMET, CHIRPS, AgERA5) download once into a cache and extract all points locally. AgERA5 submits its CDS variable requests **concurrently** so the server-side queue waits overlap and keys cache entries by geographic area. With `check_weather_downloads: true`, both drivers validate fixed-width DSSAT rows, consecutive dates, and plausible forcing before reuse; an invalid output is removed before each of up to `weather_download_retries` genuine regeneration attempts.
- **Point failures stay local:** a DSSAT error at one grid cell is recorded in that point's `_run_error.log`; the remaining cells continue and are combined. The parent summary reports every omitted point ID so partial spatial coverage remains explicit.
- Use **scratch storage** (fast local SSD or parallel filesystem) for `base_dir` — DSSAT creates many small files per run.
- Use `--cleanup_mode always` in production to remove transient DSSAT files after each point.
- Always **debug with a single point** first (`dssat_cores: 1`).
- For very large domains, set `zip_for_hpc: true` and run on HPC rather than locally.

---

## Reproducibility checklist

- Record DSSAT version and executable path used.
- Record the template file(s) and all cultivar/ecotype/species files.
- Record soil and weather data sources and year ranges.
- Save the `README_CONFIG.txt` written by the pipeline to the run folder (it captures key configuration automatically).
- Commit the study YAML so every run can be traced back to a specific configuration.

---

## Publishing checklist

**1 — DSSAT executables and licensed content**
- Official DSSAT downloads (via dssat.net) are tied to a licence. **Do not** redistribute proprietary DSSAT installers or binaries.
- The open-source DSSAT-CSM-OS codebase (BSD 3-Clause) may be redistributed as compiled binaries if you keep the licence and document provenance (commit hash, build flags, platform).
- The legacy files under `dssat_executable/` have unknown build provenance and are not self-contained; see `dssat_executable/README.md`. Do not publish them as release assets until provenance is established or they are replaced.
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
git clone https://github.com/DSSAT/dssat-csm-os.git dssat_csm_os
git clone https://github.com/DSSAT/dssat-csm-data.git dssat_csm_data
mkdir -p ~/Documents/GitHub/DSSAT48
cp -r ~/Documents/GitHub/dssat_csm_data/. ~/Documents/GitHub/DSSAT48/
```

#### Part 3 — Compile

```bash
cd ~/Documents/GitHub/dssat_csm_os
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
    ~/Documents/GitHub/dssat_csm_os/Data/DSSATPRO.L48 \
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

git clone https://github.com/DSSAT/dssat-csm-os.git dssat_csm_os
git clone https://github.com/DSSAT/dssat-csm-data.git dssat_csm_data
sudo mkdir -p /opt/DSSAT48
sudo cp -r dssat_csm_data/. /opt/DSSAT48/

cd dssat_csm_os && mkdir build && cd build
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

Run `renv::restore()` once before sourcing `dssat_main_pipeline.R`. The pipeline
checks every required package and stops with the missing package names; it never
changes the environment during a scientific run. This keeps Windows, Linux,
macOS, workstation, and CI executions on the recorded dependency graph. It
stops with a clear message naming exactly which package is missing,
instead of a cryptic `there is no package called …` error mid-source.

> **Using RStudio?** RStudio may point at a *different* R installation (and
> package library) than your command-line `Rscript`. If a package looks
> "installed" in a terminal but the pipeline reports it missing, install it in
> the R library used by RStudio and rerun. The pipeline never installs packages.

The full set recorded by the lockfile includes:

```r
# Core (always required)
install.packages(c(
  "sf", "dplyr", "tidyr", "stringr", "lubridate",
  "foreach", "doParallel", "parallel",
  "zoo", "R.utils", "processx", "tools",
  "ggplot2", "readr", "tibble", "rstudioapi", "pbapply"
))

# DSSAT R interface
install.packages("DSSAT")
# If not on CRAN: remotes::install_github("palderman/DSSAT")

# Soil
install.packages("soilDB")                           # SSURGO
install.packages(c("terra", "httr", "jsonlite"))     # SoilGrids online / VRT
install.packages(c("terra", "DBI", "RSQLite"))       # HWSD

# Weather (match your WEATHER_SOURCE)
install.packages("daymetr")                          # DAYMET
install.packages("nasapower")                        # NASA_POWER / NASA_POWER_CHIRPS / NASA_POWER_CHIRPS_V3
install.packages(c("terra", "ncdf4", "httr"))        # GRIDMET
install.packages(c("httr", "jsonlite"))              # OPEN_METEO
install.packages(c("nasapower", "terra"))            # NASA_POWER_CHIRPS / NASA_POWER_CHIRPS_V3
install.packages(c("ecmwfr", "terra"))               # CDS weather; then run setup_cds_credentials()
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
