# Weather and Soil Data Sources Reference

Provider implementations are installed from `cropmodel/packages/cropmodel_data`
(Python `cropmodel_data`, R `cropmodeldata`). Existing provider caches stay in
the configured input roots; migration does not move or regenerate them.

This document outlines all weather and soil data sources available in `cropmodel_data`, their coverage, requirements, and limitations.

## Package APIs versus tutorial configuration

This catalog describes the shared package, not a promise that every API is a
`weather_source` or `soil_source` value in the tutorial. The tutorial currently
routes weather sources DAYMET, NASA_POWER, GRIDMET, OPEN_METEO, NASA_POWER_CHIRPS,
NASA_POWER_CHIRPS_V3, AGERA5, DWD, EOBS, XAVIER, CMFD, CHELSA_W5E5, AGMERRA,
AGCFSR, SILO, PRISM, MSWX, MSWEP, CRUJRA and TERRACLIMATE. `gridpoints_weather`
uses existing files. ERA5-Land and the additional package sources listed below
must be called through their package APIs until both tutorial dispatchers support
them. See config.yml for source-specific options and required local datasets.

## Weather Data Sources

### 1. **Open-Meteo** (Global, Free, ERA5-Seamless)
- **Coverage**: Global (all land areas)
- **Spatial Resolution**: ~11 km (0.1°)
- **Temporal Resolution**: Daily (1979–present)
- **API Key Required**: No
- **Optional Dependencies**: None (requests only)
- **Status**: ✅ Stable, recommended for global studies
- **Implementation**: HTTP requests to open-meteo.com
- **Notes**: 
  - No API key or registration required
  - Fast (~1-2s per point)
  - ERA5-Land temperature/humidity combined with ERA5 forcing variables
  - Daily mean dewpoint and relative humidity are written to `TDEW` and `RH2M`
  - Free tier suitable for research

### 2. **Daymet** (US only, Free)
- **Coverage**: CONUS + Alaska (US only)
- **Spatial Resolution**: 1 km grid
- **Temporal Resolution**: Daily (1980–present)
- **API Key Required**: No
- **Optional Dependencies**: `daymetr` R package
- **Status**: ✅ Implemented; requires the R `daymetr` package on the R path
- **Implementation**: Raster download via ORNL DAAC
- **Notes**:
  - Requires `daymetr` package: `install.packages("daymetr")`
  - Only covers US territory (CONUS, Alaska, Hawaii)
  - Very high quality: ground station interpolation
  - Can be slow for large regions (~30s for single point)

### 3. **NASA-POWER** (Global, Free)
- **Coverage**: Global (land + water)
- **Spatial Resolution**: 0.5° (~56 km)
- **Temporal Resolution**: Daily (1981–present)
- **API Key Required**: No
- **Optional Dependencies**: None
- **Status**: ✅ Stable
- **Implementation**: HTTP requests to POWER API
- **Notes**:
  - Fast and reliable (~1-3s per point)
  - Coarser than Daymet/GridMET but better than nothing outside US
  - Used by many global crop models
  - Free tier allows unlimited requests

### 4. **GridMET** (US only, Free)
- **Coverage**: CONUS + parts of Mexico/Canada
- **Spatial Resolution**: 4 km grid
- **Temporal Resolution**: Daily (1979–present)
- **API Key Required**: No
- **Optional Dependencies**: `xarray`, `netCDF4` Python packages
- **Status**: ✅ Implemented; requires NetCDF-capable dependencies (`xarray`/`netCDF4` in Python, `httr`/`ncdf4` in R)
- **Implementation**: NetCDF remote files via Thredds server
- **Notes**:
  - Requires optional deps: `pip install xarray netcdf4`
  - High quality gridded data
  - Slower than Open-Meteo (~10-20s per point) due to NetCDF parsing
  - Can fail if server is down or network issues occur
  - **Variable Relevance for DSSAT:**
    * **Baseline variables (`tmmn`, `tmmx`, `pr`, `srad`):** Essential. Always required.
    * **Wind Speed (`vs`):** Wind speed at 10m. Crucial if using the **FAO-56 Penman-Monteith method** for potential evapotranspiration (PET) in DSSAT (mapped to the `WIND` column in `.WTH` files). If omitted, DSSAT falls back to simpler Priestley-Taylor ET or default wind constants.
    * **Specific Humidity (`sph`):** Specific humidity. Used to calculate actual Relative Humidity (`RHUM`) or Dewpoint Temperature (`TDEW`, written to `.WTH` files) for Penman-Monteith ET calculations, avoiding the need for empirical temperature-based approximations (e.g. `Tmin - 2.5`).
    * **Reference ET (`pet`):** Alfalfa reference ET. Not directly consumed by standard DSSAT since DSSAT computes potential and actual ET dynamically based on crop leaf area and water balance. It is valuable as a baseline for model validation or regional water policy comparisons.

### 5. **AgERA5** (Global, Requires CDS Key)
- **Coverage**: Global land areas
- **Spatial Resolution**: 0.1° (~11 km)
- **Temporal Resolution**: Daily (1979–present, 10-day lag)
- **API Key Required**: Yes (Copernicus Data Store CDS account)
- **Optional Dependencies**: Python `cdsapi` via `cropmodel_data[cds]`; R `ecmwfr`
- **Status**: ⚠️ Requires CDS API credentials
- **Implementation**: CDS API requests
- **Notes**:
  - Requires free CDS account setup: https://cds.climate.copernicus.eu/
  - Run `cropmodeldata::setup_cds_credentials()` in R or `cropmodel_data.setup_cds_credentials()` in Python, or create `~/.cdsapirc`
  - Accept the AgERA5 dataset licence once in the CDS web UI
  - High-quality fusion of ERA5 + station data
  - Daily agrometeorological indicators downloaded once per variable-year and cached
  - The adapter writes provider values without a separate physical-quality gate. The shared engine validator then checks dates, ranges, temperature ordering, and all seven AgERA5 forcing columns consistently with the other weather sources.
  - Time-series backend validates cached annual CSVs per grid cell, acquires locks across workers, and enforces exact calendar coverage (refusing partial weather histories).
  - Both R and Python assemble one tile across all requested years, publish complete point weather files, then release the tile data before advancing. Concurrent requests are limited to years within the current tile; memory scales with its points and history length rather than the whole study. Existing annual caches remain reusable.
  - Optional weather repairs (e.g. bounded temperature inversion swap/neighbor repair) run before simulation validation, and newly valid points are pruned from `unresolvable_points.json`.

### 6. **NASA-POWER CHIRPS Hybrid** (Global, Free)
- **Coverage**: Global (focus on tropics/subtropics for CHIRPS)
- **Spatial Resolution**: 0.05° (~5.5 km for CHIRPS)
- **Temporal Resolution**: Daily (1981–present)
- **API Key Required**: No
- **Optional Dependencies**: None (requests only)
- **Status**: ✅ Implemented as a NASA-POWER + CHIRPS rainfall hybrid
- **Implementation**: Hybrid NASA-POWER (temp/humidity) + CHIRPS (rainfall)
- **Notes**:
  - Combines NASA-POWER weather with high-res rainfall from CHIRPS
  - CHIRPS coverage: |Lat| < 50° (tropics/subtropics mainly)
  - Falls back to NASA-POWER rain if outside CHIRPS coverage
  - Good for tropical/subtropical regions with intense rainfall
  - Slower than NASA-POWER alone (~5-10s per point)

> The v0.3.0 backends (PRISM, CHELSA-W5E5, AgMERRA/AgCFSR, SILO, MSWX, MSWEP,
> CRU-JRA, TerraClimate) are catalogued in the **priority table** under the Soil
> section below. 10 m winds are converted to 2 m (FAO-56, ×0.748) in the shared
> extractor, and SILO vapour pressure is converted to dewpoint (inverse Magnus).

---

## Soil Data Sources

### 1. **SoilGrids Online** (Global, Free, REST API or VRT)
- **Coverage**: Global
- **Spatial Resolution**: 250 m (finer resolution available)
- **Data Type**: Soil properties (sand/silt/clay, organic matter, pH, etc.)
- **API Key Required**: No
- **Optional Dependencies**: `sf`, `dplyr`, `httr`, `terra` for VRT (R); `geopandas`, `requests`, `rasterio` for VRT (Python)
- **Status**: ✅ Implemented; REST is slow/rate-limited, VRT is preferred for batch runs
- **Implementation**: ISRIC SoilGrids REST API or GDAL VRT raster sampling
- **Notes**:
  - Select via `soilgrids_mode` in the tutorial `config.yml`; direct package calls use merged `cropmodel_data` `config.yml` (`soil.soilgrids_online.use_rest_api`)
  - REST mode makes one HTTP request per point
  - VRT mode streams global SoilGrids rasters and samples all points locally
  - Coverage gaps over water, deserts, some urban areas
  - Expected to produce partial results (some points may fail)
  - REST API is slow but does not require local data; VRT needs GDAL-backed raster support
  - **CI Recommendation**: Skip or set long timeout (>1 min for 3 points)

### 2. **SoilGrids (Offline)** (Global, Free, Local VRT/GeoTIFF)
- **Coverage**: Global
- **Spatial Resolution**: 250 m
- **Data Type**: Soil properties (same as SoilGrids Online)
- **API Key Required**: No
- **Optional Dependencies**: `sf`, `terra`/`raster` (R); `geopandas`, `rasterio` (Python)
- **Status**: ✅ Implemented for local SoilGrids raster/VRT inputs
- **Implementation**: Local GeoTIFF files via GDAL
- **Notes**:
  - Requires pre-downloaded SoilGrids GeoTIFF files (~30 GB for global)
  - Much faster than REST API (~1-2s per point) once cached
  - Can be set up with `.vrt` (Virtual Raster) for efficient access
  - Recommended for production: cache locally, use offline mode

### 3. **SSURGO** (US only, Free, USDA Soil Data Access)
- **Coverage**: CONUS only (some gaps in Alaska/Hawaii)
- **Spatial Resolution**: 30 m
- **Data Type**: Detailed soil properties (USDA Soil Taxonomy)
- **API Key Required**: No
- **Optional Dependencies**: `sf`, `terra` (R); `geopandas`, `shapely` (Python)
- **Status**: ✅ Implemented through the shared `cropmodel_data` adapter
- **Implementation**: USDA SDA spatial/tabular queries, cached layer mapping CSVs
  and derived per-point DSSAT `.SOL` files. gNATSGO uses a WCS map-unit lookup
  followed by SDA tabular queries; Alderman uses the full-profile conversion.
- **Notes**:
  - Highest quality soil data for US (30 m resolution)
  - Can have coverage gaps in mountainous/remote areas
  - Historical SSURGO/gNATSGO and Alderman layer-width errors require rebuilding
    derived `.SOL` files, not deleting the raw data. See
    [offline recovery](../cropmodel_data/README.md#input-cache-integrity-fixes-2026-08-31).
  - R GRIDMET's historical invalid-cell row shift likewise requires re-extracting
    affected `.WTH` files from the existing native NetCDF cache. Valid-looking
    weather dates/values alone cannot detect weather assigned to the wrong ID.

### 4. **HWSD** (Global, Free, Local Database)
- **Coverage**: Global
- **Spatial Resolution**: 1 km
- **Data Type**: Coarser soil properties (FAO Harmonized World Soil Database)
- **API Key Required**: No
- **Optional Dependencies**: `sf`, `terra`/`raster` (R); `geopandas`, `rasterio` (Python)
- **Status**: ✅ Implemented in `cropmodel_data` for local HWSD raster + database inputs
- **Implementation**: Local SQLite database
- **Notes**:
  - Coarser than SoilGrids (1 km) and SSURGO (30 m) but global coverage
  - Good fallback for regions without SSURGO
  - Faster than SoilGrids REST API (~2-3s per point)
  - Requires local HWSD database (~2-3 GB)
  - Database: https://www.fao.org/documents/card/en/c/CA12305EN/

### 5. **AgMIP/Han DSSAT Soil Profile Database** (Global cropland, Free, Local `.SOL`)
- **Coverage**: Global cropland, country files
- **Spatial Resolution**: 5 arc-min (~10 km)
- **Data Type**: DSSAT-ready soil profiles to 2 m depth
- **API Key Required**: No
- **Optional Dependencies**: Same as the external `.SOL` mapper (`sf`/`geopandas` useful; haversine fallback in Python)
- **Status**: ✅ Implemented as `AGMIP` / `process_soils_agmip`
- **Implementation**: Reads a downloaded country `.SOL` master file and selects the nearest source profile for each grid point
- **Notes**:
  - Data DOI: https://doi.org/10.7910/DVN/1PEEY0
  - Source paper: Han, Ines, and Koo (2019), *Environmental Modelling & Software* 119:70-83
  - This is the most DSSAT-specific global soil source: hydraulic properties are already present in DSSAT format, so the pipeline does not recompute PTFs

---

## Candidate Sources Worth Implementing Next

The following sources were reviewed online in June 2026 and have now been implemented as local/cache-backed backends.

| Priority | Source | Type | Implementation | Main caveat |
|----------|--------|------|------------------------|-------------|
| 1 | **[HiHydroSoil v2.0](https://www.futurewater.eu/tools/hihydrosoil/)** | Soil hydraulics, global 250 m | `HIHYDROSOIL` / `process_soils_hihydrosoil`, reading local GeoTIFF/VRT rasters. | Requires local rasters; hydraulic layers are often integer-scaled by 10000, so `integer_scale` defaults to 0.0001. |
| 2 | **[CHELSA-W5E5 daily](https://www.chelsa-climate.org/datasets/chelsaw5e5)** | Weather, global 30 arcsec (~1 km), 1979-2016 | `CHELSA_W5E5` / `process_weather_chelsa_w5e5`, reading local NetCDF files. | Historical only through 2016; humidity/wind are not required and are written missing unless provided by companion files. |
| 3 | **[SILO](https://www.longpaddock.qld.gov.au/silo/) + [SLGA](https://www.csiro.au/en/research/natural-environment/land/soil-and-landscape-grid-of-australia)** | Weather + soil, Australia | `SILO` / `process_weather_silo` for local NetCDF weather; `SLGA` / `process_soils_slga` for local soil rasters. | Regional Australia-only pair; data access/download is managed outside the pipeline cache. |
| 4 | **[AgMERRA / AgCFSR](https://data.giss.nasa.gov/impacts/agmipcf/)** | Weather, global 0.25 degree, 1980-2010 | `AGMERRA` and `AGCFSR` / local NetCDF readers. | Ends in 2010, so it is for historical baselines and AgMIP intercomparison. |
| 5 | **[PRISM](https://prism.oregonstate.edu/)** | Weather, US 4 km | `PRISM` / `process_weather_prism`: default `prism_backend: "acis"` queries daily point histories; `"nacse"` downloads CONUS rasters. Both use `prism_cache_dir`. | SRAD is estimated with Bristow-Campbell by default; ACIS dewpoint/wind/RH are -99. Versioned point caches validate coordinates, dates, variables and units; incomplete core forcing is retried and not published. Old ID-only caches are preserved but not trusted. |
| 6 | **WISE30sec** | Soil, global 30 arcsec | `WISE30SEC` / `process_soils_wise30sec`, reading local GeoTIFF/VRT rasters. | Uses WISE30sec's own 7 depth layers — rasters must carry a property token (sand/clay/silt/BD/OC) **and** a depth token `d1`–`d7` in the filename. WISE is natively a map-unit raster + attribute table; rasterize that join to per-property/per-depth grids first. |
| 7 | **MSWX / MSWEP** | Weather, global ~0.1 degree | `MSWX` local full-weather NetCDF reader; `MSWEP` NASA POWER + local MSWEP rainfall hybrid. | MSWEP is precipitation-only, so it is not a standalone full-weather source. |
| 8 | **CRU-JRA** | Weather, global 0.5 degree | `CRUJRA` / local NetCDF reader. | Coarser than the other gridded products. |
| 9 | **TerraClimate** | Weather/climatology, global monthly | `TERRACLIMATE` / local NetCDF reader; each month is disaggregated to continuous daily records (constant T/SRAD/wind, monthly precip spread evenly) so the `.WTH` is runnable. | Carries **no day-to-day variability** (rain intensity, dry spells, heat-stress days) — screening/climatology only, not production. |
| 10 | **WoSIS** | Soil point profiles | `WOSIS` / processed layer CSV nearest-profile ingestion. | Raw WoSIS exports must be harmonized before use; best for calibration/validation. |

### Regional coverage retained in cropmodel_data / cropmodeldata

Targeted at regions previously served only by coarse global reanalysis. Two are
**live downloads** (no local data); the rest read local NetCDF / rasters.

| Source | Function | Type / Coverage | Access | Caveat |
|--------|----------|-----------------|--------|--------|
| **APHRODITE** | `process_weather_aphrodite` | Weather, monsoon Asia (0.25-0.5°, gauge precip) | Local NetCDF (rain) + NASA-POWER | Rainfall-only hybrid: T/SRAD/RH/wind from NASA-POWER. |
| **ANUSPLIN** | `process_weather_anusplin` | Temperature/precipitation, Canada 10 km (1950-2015) | Local NetCDF | Core ANUSPLIN has no SRAD; standalone DSSAT WTH creation is rejected until a solar-radiation layer is supplied. |
| **TAMSAT** | `process_weather_tamsat` | Weather, Africa 4 km rainfall | Local NetCDF (rain) + NASA-POWER | Rainfall-only hybrid; pairs with iSDAsoil. |
| **GHCN-Daily** | `process_weather_ghcn` | Weather, global station obs | **Live (NOAA)** | Nearest-station snap; per-station record ~11 MB; SRAD/RH/wind = -99. Best as obs/validation. |
| **Princeton PGF** | `process_weather_pgf` | Weather, global 0.25° | Local NetCDF | Full-variable reanalysis alt. to AgMERRA/CRU-JRA. |
| **MERRA-2** | `process_weather_merra2` | Weather, global ~0.5° (1980-) | Local NetCDF | Kelvin temps auto-convert. |
| **GSDE** | `process_soils_gsde` | Soil, global 1 km, 8 layers to 2.3 m | Local rasters (`l1`-`l8` tokens) | Richer vertical profile than HWSD. |
| **China BNU** | `process_soils_china` | Soil, China 1 km, 8 layers | Local rasters (`l1`-`l8`) | Pairs with CMFD weather. |
| **FEBR / Embrapa** | `process_soils_febr` | Soil, Brazil point profiles | Harmonized CSV | Nearest-profile (WoSIS mechanism); pairs with Xavier. |
| **SLC** | `process_soils_slc` | Soil, Canada (top/sub) | Local rasters (`top`/`sub`) | Rasterize the SLC polygon+component join first; pairs with ANUSPLIN. |
| **ESDB** | `process_soils_esdb` | Soil, Europe full profile (top/sub) | Local rasters (`top`/`sub`) | Full profile vs topsoil-only LUCAS. |
| **OpenLandMap** | `process_soils_openlandmap` | Soil, global 250 m (0-30/30-60/60-100) | **Live (STAC COG)** | No local data; texture 120 m + BD/OC 250 m sampled over HTTP via the OpenLandMap STAC catalog. |

---

## Regional Selection Guide

### **For Global Studies (Use One):**
| Source | Coverage | Speed | Quality | Notes |
|--------|----------|-------|---------|-------|
| **Open-Meteo** | ✅ Global | Fast | Good | **Recommended**: Free, fast, no setup |
| NASA-POWER | ✅ Global | Fast | Fair | Coarser resolution (0.5°) |
| AgERA5 | ✅ Global | Slow | Excellent | Requires CDS API key + setup |
| SoilGrids (Online) | ✅ Global | Very Slow | Good | REST API: ~30s/point, has gaps |
| SoilGrids (Offline) | ✅ Global | Fast | Good | Requires ~30 GB local cache |
| HWSD | ✅ Global | Fast | Fair | Coarser (1 km) but fully coverage |

### **For US Studies (Use Regional + Soil):**
| Component | Source | Speed | Quality |
|-----------|--------|-------|---------|
| **Weather** | Daymet | Slow | Excellent | 1 km, ground-based interpolation |
| | GridMET | Slow | Excellent | 4 km, comprehensive variables |
| | Open-Meteo | Fast | Good | Fallback if others fail |
| **Soil** | SSURGO | Fast | Excellent | 30 m, detailed USDA data |
| | SoilGrids (Offline) | Fast | Good | Faster alternative if SSURGO unavailable |

### **For Tropical/Subtropical Studies:**
| Component | Source | Notes |
|-----------|--------|-------|
| **Weather** | NASA-POWER CHIRPS | Combines NASA-POWER + high-res rainfall (CHIRPS) |
| | Open-Meteo | Simple alternative |
| **Soil** | SoilGrids (Online/Offline) | Global coverage, REST API slow but workable |

---

## Testing & CI Integration

### **Smoke Test (Fast, Offline):**
- ✅ Config loading
- ✅ Module imports
- ✅ .WTH/.SOL file format validation
- ⏱️ Runtime: ~30 seconds

### **E2E Test (Weather Only, Network):**
- ✅ Open-Meteo: 1-2 minutes for 3 global points
- ✅ NASA-POWER: 1-2 minutes for 3 global points
- ⚠️ Daymet: Slow (~30s/point, skip on slow networks)
- ⚠️ GridMET: Slow (~20s/point, can timeout)
- ❌ SoilGrids REST: Very slow (10-30s/point, skip in CI)
- ⏱️ Total runtime: 2-5 minutes (weather only)

### **Recommended CI Strategy:**
1. **Smoke test** (always run): Fast validation of imports/formats
2. **Weather E2E** (always run): Open-Meteo + NASA-POWER (fast + reliable)
3. **Optional E2E** (scheduled/manual): Full test with all sources + soil (slow)
4. **Cache frequently-used soil data** (SoilGrids offline, SSURGO local)

---

## Setup Instructions

### **Open-Meteo (No Setup)**
```bash
# Works out of the box, no API key needed
```

### **Daymet (R)**
```r
install.packages("daymetR")
# Also requires: sf, terra, dplyr
```

### **GridMET (R)**
```r
install.packages("httr")  # For HTTP GET requests
```

### **AgERA5 (Setup Required)**
```r
# 1. Create free CDS account at https://cds.climate.copernicus.eu/
# 2. Copy your Personal Access Token from your CDS profile.
# 3. Configure credentials for R and Python clients:
library(cropmodeldata)
setup_cds_credentials()
```

For Python-only environments, install the CDS extra and run:

```python
from cropmodel_data import setup_cds_credentials
setup_cds_credentials()
```

Accept the AgERA5 dataset licence in the CDS web UI before the first run.

### **SoilGrids Offline (Advanced)**
```bash
# Download SoilGrids GeoTIFF files (~30 GB)
# Create .vrt file for efficient access
# Configure cropmodel_data to use local path
```

### **SSURGO (US Only)**
```bash
# 1. Visit https://websoilsurvey.nrcs.usda.gov/
# 2. Download SSURGO database for your region
# 3. Configure cropmodel_data SSURGO path
```

---

## Troubleshooting

| Error | Source | Solution |
|-------|--------|----------|
| `could not find function "leap_year"` | Daymet | Install `daymetR` package: `install.packages("daymetR")` |
| `could not find function "GET"` | GridMET | Load `httr` library: `library(httr)` |
| `could not find function "st_transform"` | SoilGrids REST | Load `sf` library: `library(sf)` |
| `No soil data extracted` | SoilGrids REST | REST API timeout or coverage gap (normal, graceful skip) |
| Timeout after 60s | GridMET/Daymet | These sources are slow; increase timeout or use Open-Meteo |
| HTTP 429 (rate limit) | Any source | Add delay between points or use async requests |
| SSL certificate error | Any source | Update `curl`/`openssl`: `brew install curl openssl` (macOS) |

---

## Performance Benchmarks (Single Point, Local Network)

| Source | Time | Notes |
|--------|------|-------|
| Open-Meteo | ~1 s | Fastest, reliable |
| NASA-POWER | ~1-2 s | Fast, reliable |
| NASA-POWER CHIRPS | ~5-10 s | Hybrid processing |
| Daymet | ~30 s | Ground-station interpolation |
| GridMET | ~10-20 s | NetCDF parsing |
| AgERA5 | ~5-10 s | CDS API overhead |
| SoilGrids (Offline) | ~1-2 s | Fast local access |
| SoilGrids (REST) | ~10-30 s | HTTP requests + processing |
| SSURGO | ~1-2 s | Fast local DB |
| HWSD | ~2-3 s | SQLite queries |

---

## Reference Links

- **Open-Meteo**: https://open-meteo.com/
- **Daymet**: https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1328
- **NASA-POWER**: https://power.larc.nasa.gov/
- **GridMET**: https://www.climatologylab.org/gridmet.html
- **AgERA5**: https://www.copernicus.eu/en/access-data/copernicus-data-space-ecosystem
- **SoilGrids**: https://soilgrids.org/
- **SSURGO**: https://websoilsurvey.nrcs.usda.gov/
- **HWSD**: https://www.fao.org/documents/card/en/c/CA12305EN/

AgERA5 recovery integrity (2026-09-06): the mirrored downloaders validate annual
per-cell calendars and all seven forcing variables, serialize cache publication
with compatible OS locks, and preserve failed-cache evidence. The R dependency
`filelock` is required. Full boundary-date validation and retryable weather
failure records are described in the engine README; source-data inversions still
require the explicitly configured shared repair policy.

Generic daily NetCDF sources now require explicit physical units and geographic coordinates.
Supply wind `height_m` metadata; ambiguous units, projected coordinates, subdaily or duplicate dates,
and outside-domain points must be resolved before acquisition. See cropmodel_data README for the contract.

### Wind serialization for the active Carinata comparison adapters

NASA POWER, NASA POWER + CHIRPS, and AgERA5 wind is supplied in m/s and must
be multiplied by 86.4 for DSSAT v4.8 WTH `WIND` (km/day; see `WEATHER.CDE`).
The shared R/Python writers perform this conversion at serialization, retaining
missing markers and raw provider cache units. Older generated WTH caches require
writer-provenance verification and regeneration or a backed-up conversion;
do not guess units from magnitude or rescale already-correct files.

Daymet station metadata must preserve DSSAT header column boundaries; the R
writer now matches Python rather than shifting its coordinates one column.
GRIDMET AMP uses the full mean annual range of monthly temperatures, matching
the other active adapters. DSSAT STEMP divides that value by two internally.
Existing weather files require deliberate header regeneration for these fixes.
The five active comparison adapters also reserve nine station-longitude
columns in both R and Python. Three-digit negative longitudes must not shift
ELEV, TAV or AMP relative to the header labels; whitespace token checks alone
do not establish correct DSSAT parsing.
AgERA5 writes two temperature decimals for pairs within 0.05°C of zero so
rounding does not unnecessarily create a DSSAT-rejected zero pair. True source
zeros remain zero. Existing affected caches need repair from the raw provider
values, not an arbitrary temperature offset.


Download recovery: transient provider requests use five bounded attempts with
5/10/20/40-second default waits (configurable helper cap: 60 seconds). The R and
Python weather drivers catch exhausted connectivity errors, recheck completed
files, and retry only missing or invalid points within the configured batch
budget. Non-transient configuration errors remain fatal. Exhausted weather
points stay retryable; download-only runs report partial coverage with exit 2.
SoilGrids batch retries preserve existing profiles and do not treat historical
SKIPPED log entries as permanent geographic exclusions. Required raster layers
must succeed before any new profiles are published. The Bioenergy scheduler
continues independent source/resolution tasks after a child failure and reports
all failures in its final summary.

AgERA5 time-series caches allow fully masked cells alongside complete cells in a coastal tile. Partial missing data and incomplete calendars still fail validation; each site keeps its nearest cell and requires complete forcing.
# Soil cache recovery

Existing derived SOL profiles must pass the shared fixed-column and hydraulic
validator before reuse. Missing optional SoilGrids silt/coarse fragments are
serialized as DSSAT `-99`; missing required physical properties remain failures.
POLARIS uses consistent integer tile ownership and a 1e-9-degree numerical
tolerance for raster-edge roundoff, without relocating interior sampling points.



## Optional bounded land sampling

When disabled, extraction uses the configured point. Both R and Python can
optionally retry a missing soil or weather input at land-mask pixel centres,
ordered by geodesic distance. Every accepted candidate must produce one complete,
validated soil profile or the entire requested weather period at the same location.
The model point and ID remain unchanged. Values from different locations are never
mixed. Transport, authentication and corrupt-raster errors remain retryable errors;
they do not justify a geographic substitution.

```yaml
input_fallback_enabled: false
input_fallback_land_mask: ""
input_fallback_land_values: []
input_fallback_max_distance_m: 500
input_fallback_max_candidates: 10000
input_fallback_within_cell: true
input_fallback_cell_crs: "EPSG:5070"
```

To enable land mode, provide a readable categorical raster and its explicit valid land
classes. For example, NLCD class 82 means cultivated crops; codes from another
land-cover product must be chosen from that product's legend. The search is
bounded by distance, candidate count and, by default, the original square study
cell. A count limit can stop the search before the full radius is exhausted.
Candidates are centres of the supplied mask, rather than necessarily centres of
the source raster. Station and prebuilt-profile sources retain their native
matching rules. This feature therefore finds the nearest *tested mask location*
that returns a valid input; it does not establish equivalent soil or climate.

The policy creates distinct input-cache directories. Shifts and failures are
written to run logs; no soil/weather provenance sidecars are added. Changing the
contents of a mask at the same path requires clearing the affected input cache.
NoData alone does not establish that the original point is water or a quarry.
A cropland centroid can also fall outside cultivated land.

SOILGRIDS_ONLINE samples the native product. SOILGRIDS_10K selects prebuilt source
profiles; it does not calculate an area-weighted mean over the current study grid
cell. Optional relocation preserves these source-specific extraction rules.


## Sampling mode selection (11 October 2026)

`input_fallback_enabled: false` disables geographic fallback. When enabled,
`input_fallback_mode` selects one of:

- `nearest_valid`: no land mask required. Test locations on a local metric lattice
  (`input_fallback_search_step_m: 50`) in increasing distance order, bounded by
  radius, candidate cap and study cell. This finds the nearest valid *tested
  location*, not necessarily the exact nearest source pixel. It does not establish
  that the replacement is land or cropland.
- `nearest_valid_land`: require `input_fallback_land_mask` and explicit
  `input_fallback_land_values`; test eligible mask pixel centres as previously.

Both modes require complete validated inputs at one replacement location and
retain the same network/error safeguards. The Bioenergy study configuration now
activates `nearest_valid` with a 500 m radius and 50 m step. Change the mode to
`nearest_valid_land` and supply mask/classes to enable filtering; set enabled to
false to disable. Policy keys separate caches for the modes and search spacing.
Existing files remain untouched; running processes need a relaunch to adopt it.


Converted soil/weather inputs no longer require conversion-provenance sidecars.
Existing markers are ignored and retained; file-content validation still applies.
The shared R fallback wrapper binds each provider when registered, so regional
and point-validation loops retain the correct provider signature.


Fallback recovery treats explicit per-point NoData/outside-coverage outcomes
independently of connection failures reported for other points in the same batch.
Unknown failures remain retryable. The exact SoilGrids all-null response permits
spatial fallback; it does not establish permanent geographic exclusion.


### Input recovery and independent relocation switches

Regional input preparation and sweeps use the top-level `input_fallback_enabled`
setting. Set it to `false` to download at the configured original grid points
without searching nearby coordinates. Point validation independently uses
`validation_pipeline.input_fallback_enabled`; set that to `false` for direct
site/selected-point downloads or `true` for the configured nearest-valid mode.
An explicit `INPUT_FALLBACK_ENABLED` environment override takes precedence.
Disabling relocation does not request entire provider raster archives: the
configured sources, points, years and provider cache mechanisms still apply.
Relocated and direct-point inputs use separate cache identities.

Configured weather repairs are applied before deciding to search for another
coordinate. DAYMET multi-year HTTP 5xx failures are retried as individual years
at the same coordinate. POLARIS no-data and Alderman query failures retain
distinct missing-data versus retryable outcomes. The Alderman field-capacity
fallback is unchanged. R and Python use a shared operating-system file lock
per input directory while preparing inputs; it is released when the process
exits. All concurrent processes must use the updated package and drivers.


## Shared SoilGrids raster cache

The R/Python soil stages accept `soilgrids_cache_dir`, `soilgrids_cache_download`,
`soilgrids_cache_scope` (region/world), `soilgrids_cache_bbox`, and
`soilgrids_cache_buffer_degrees`. Blank directory keeps remote access; a relative
directory resolves under `provider_cache_root` (default `input_root_dir`). The
Bioenergy downloader enables `soilgrids_raster_cache` here, beside the weather
caches. Region acquisition uses the full prepared grid's WGS84 bounds with a
configurable margin; world acquisition is explicitly opt-in. Both use the shared
`cropmodel_data` tile downloader; R acquisition uses its Python bridge, while
R/Python extraction stays native. Failed acquisition is an explicit resumable
soil-task error. Local tile reads never silently fall back to the network.
The updated shared packages must be installed; historical remote pins do not yet
include this unreleased feature. Existing caches and generated profiles are preserved.


Resolution-relative fallback: `input_fallback_resolution_fraction: 0.10` sets the maximum distance to 10% of the grid width (300 km → 30 km; 200 → 20; 100 → 10; 50 → 5; 25 → 2.5; 10 → 1; 5 → 0.5 km). Zero retains `input_fallback_max_distance_m`. `input_fallback_enabled` independently enables relocation. Exact-site validation uses the fixed metre limit; regional validation supplies its actual resolution through `input_fallback_resolution_m` and preserves the original grid-cell centre. Keep `input_fallback_within_cell: true` to constrain candidates to that cell.

Relative nearest-valid searches retain the configured fine spacing in the first 20 steps and double the spacing in successive outward bands, keeping large searches bounded. Candidates are globally ordered by distance and limited by `input_fallback_max_candidates`; the first valid tested coordinate wins. This does not guarantee the nearest valid native source pixel. Network/server failures remain retryable and do not trigger relocation. The radius policy and resolution are part of the fallback cache key.

## Sampling provenance

Explicit sampling uses the declared crop mask, independent of the compared soil product. Provider recovery records per-layer requested and recovered coordinates and input hashes. Legacy caches without records are unknown. Cell/stratum supports and crop eligibility constrain recovery. A request coordinate is not a guarantee of native raster-pixel location. See the sibling study SAMPLING_EXPERIMENT.md.
