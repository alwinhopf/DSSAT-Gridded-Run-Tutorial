# Weather and Soil Data Sources Reference

This document outlines all weather and soil data sources available in `dssatutils`, their coverage, requirements, and limitations.

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
- **Optional Dependencies**: Python `cdsapi` via `dssatutils[cds]`; R `ecmwfr`
- **Status**: ⚠️ Requires CDS API credentials
- **Implementation**: CDS API requests
- **Notes**:
  - Requires free CDS account setup: https://cds.climate.copernicus.eu/
  - Run `dssatutils::setup_cds_credentials()` in R or `dssatutils.setup_cds_credentials()` in Python, or create `~/.cdsapirc`
  - Accept the AgERA5 dataset licence once in the CDS web UI
  - High-quality fusion of ERA5 + station data
  - Daily agrometeorological indicators downloaded once per variable-year and cached
  - Requests are queued server-side; transient CDS/DNS failures should be retried with the same cache directory
  - The adapter writes provider values without a separate physical-quality gate. The shared engine validator then checks dates, ranges, temperature ordering, and all seven AgERA5 forcing columns consistently with the other weather sources.

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
  - Select via `soilgrids_mode` in the tutorial `config.yml`; direct package calls use merged `dssatutils` `config.yml` (`soil.soilgrids_online.use_rest_api`)
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

### 3. **SSURGO** (US only, Free, Local Database)
- **Coverage**: CONUS only (some gaps in Alaska/Hawaii)
- **Spatial Resolution**: 30 m
- **Data Type**: Detailed soil properties (USDA Soil Taxonomy)
- **API Key Required**: No
- **Optional Dependencies**: `sf`, `terra` (R); `geopandas`, `shapely` (Python)
- **Status**: ✅ Implemented for local SSURGO / gSSURGO-style inputs
- **Implementation**: Local SQLite database via USDA NRCS Web Soil Survey
- **Notes**:
  - Highest quality soil data for US (30 m resolution)
  - Requires local SSURGO database setup
  - Database file: ~5-10 GB for CONUS
  - Can have coverage gaps in mountainous/remote areas
  - Must download via NRCS WSS: https://websoilsurvey.nrcs.usda.gov/

### 4. **HWSD** (Global, Free, Local Database)
- **Coverage**: Global
- **Spatial Resolution**: 1 km
- **Data Type**: Coarser soil properties (FAO Harmonized World Soil Database)
- **API Key Required**: No
- **Optional Dependencies**: `sf`, `terra`/`raster` (R); `geopandas`, `rasterio` (Python)
- **Status**: ✅ Implemented in `dssatutils` for local HWSD raster + database inputs
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
| 5 | **[PRISM](https://prism.oregonstate.edu/)** | Weather, US 4 km | `PRISM` / `process_weather_prism`, downloading public daily grids into `prism_cache_dir`. | No daily SRAD/wind/RH in this path → those columns are `-99`. NACSE throttles rapid requests, so the backend adds a 1 s inter-request delay and validates each zip (throttled days degrade to a skip, not a crash). |
| 6 | **WISE30sec** | Soil, global 30 arcsec | `WISE30SEC` / `process_soils_wise30sec`, reading local GeoTIFF/VRT rasters. | Uses WISE30sec's own 7 depth layers — rasters must carry a property token (sand/clay/silt/BD/OC) **and** a depth token `d1`–`d7` in the filename. WISE is natively a map-unit raster + attribute table; rasterize that join to per-property/per-depth grids first. |
| 7 | **MSWX / MSWEP** | Weather, global ~0.1 degree | `MSWX` local full-weather NetCDF reader; `MSWEP` NASA POWER + local MSWEP rainfall hybrid. | MSWEP is precipitation-only, so it is not a standalone full-weather source. |
| 8 | **CRU-JRA** | Weather, global 0.5 degree | `CRUJRA` / local NetCDF reader. | Coarser than the other gridded products. |
| 9 | **TerraClimate** | Weather/climatology, global monthly | `TERRACLIMATE` / local NetCDF reader; each month is disaggregated to continuous daily records (constant T/SRAD/wind, monthly precip spread evenly) so the `.WTH` is runnable. | Carries **no day-to-day variability** (rain intensity, dry spells, heat-stress days) — screening/climatology only, not production. |
| 10 | **WoSIS** | Soil point profiles | `WOSIS` / processed layer CSV nearest-profile ingestion. | Raw WoSIS exports must be harmonized before use; best for calibration/validation. |

### Regional-coverage fill (dssatutils v0.4.0)

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
library(dssatutils)
setup_cds_credentials()
```

For Python-only environments, install the CDS extra and run:

```python
from dssatutils import setup_cds_credentials
setup_cds_credentials()
```

Accept the AgERA5 dataset licence in the CDS web UI before the first run.

### **SoilGrids Offline (Advanced)**
```bash
# Download SoilGrids GeoTIFF files (~30 GB)
# Create .vrt file for efficient access
# Configure dssatutils to use local path
```

### **SSURGO (US Only)**
```bash
# 1. Visit https://websoilsurvey.nrcs.usda.gov/
# 2. Download SSURGO database for your region
# 3. Configure dssatutils SSURGO path
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
