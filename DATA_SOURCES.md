# Weather and Soil Data Sources Reference

This document outlines all weather and soil data sources available in `dssatutils`, their coverage, requirements, and limitations.

## Weather Data Sources

### 1. **Open-Meteo** (Global, Free, ERA5)
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
  - Backed by ERA5 reanalysis data
  - Free tier suitable for research

### 2. **Daymet** (US only, Free)
- **Coverage**: CONUS + Alaska (US only)
- **Spatial Resolution**: 1 km grid
- **Temporal Resolution**: Daily (1980–present)
- **API Key Required**: No
- **Optional Dependencies**: `daymetR` R package (must install separately)
- **Status**: ⚠️ Optional dep missing in dssatutils (leap_year function)
- **Implementation**: Raster download via ORNL DAAC
- **Notes**:
  - Requires `daymetR` package: `install.packages("daymetR")`
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
- **Status**: ⚠️ Requires httr for HTTP GET requests (missing in R env)
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
- **Optional Dependencies**: `cdsapi` Python package
- **Status**: ⚠️ Requires CDS API credentials
- **Implementation**: CDS API requests
- **Notes**:
  - Requires free CDS account setup: https://cds.climate.copernicus.eu/
  - Create `~/.cdsapi` config file with credentials
  - High-quality fusion of ERA5 + station data
  - Slower (~5-10s per point due to API overhead)
  - Monthly data then interpolated to daily
  - 10-day lag (data published ~mid-month)

### 6. **NASA-POWER CHIRPS Hybrid** (Global, Free)
- **Coverage**: Global (focus on tropics/subtropics for CHIRPS)
- **Spatial Resolution**: 0.05° (~5.5 km for CHIRPS)
- **Temporal Resolution**: Daily (1981–present)
- **API Key Required**: No
- **Optional Dependencies**: None (requests only)
- **Status**: ⚠️ Coordinate handling issue in R (needs debugging)
- **Implementation**: Hybrid NASA-POWER (temp/humidity) + CHIRPS (rainfall)
- **Notes**:
  - Combines NASA-POWER weather with high-res rainfall from CHIRPS
  - CHIRPS coverage: |Lat| < 50° (tropics/subtropics mainly)
  - Falls back to NASA-POWER rain if outside CHIRPS coverage
  - Good for tropical/subtropical regions with intense rainfall
  - Slower than NASA-POWER alone (~5-10s per point)

---

## Soil Data Sources

### 1. **SoilGrids Online** (Global, Free, REST API)
- **Coverage**: Global
- **Spatial Resolution**: 250 m (finer resolution available)
- **Data Type**: Soil properties (sand/silt/clay, organic matter, pH, etc.)
- **API Key Required**: No
- **Optional Dependencies**: `sf`, `dplyr`, `httr` (R); `geopandas`, `requests` (Python)
- **Status**: ⚠️ Slow REST API (~10-30s per point) with coverage gaps
- **Implementation**: REST API calls to ISRIC SoilGrids service
- **Notes**:
  - One HTTP request per point (no batch mode)
  - Coverage gaps over water, deserts, some urban areas
  - Expected to produce partial results (some points may fail)
  - REST API is slow but does not require local data
  - Python implementation faster than R for parallel requests
  - **CI Recommendation**: Skip or set long timeout (>1 min for 3 points)

### 2. **SoilGrids (Offline)** (Global, Free, Local VRT/GeoTIFF)
- **Coverage**: Global
- **Spatial Resolution**: 250 m
- **Data Type**: Soil properties (same as SoilGrids Online)
- **API Key Required**: No
- **Optional Dependencies**: `sf`, `terra`/`raster` (R); `geopandas`, `rasterio` (Python)
- **Status**: ⚠️ Function signature issue in comprehensive test
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
- **Status**: ⚠️ Function signature issue in comprehensive test
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
- **Status**: ⚠️ Function signature issue in comprehensive test
- **Implementation**: Local SQLite database
- **Notes**:
  - Coarser than SoilGrids (1 km) and SSURGO (30 m) but global coverage
  - Good fallback for regions without SSURGO
  - Faster than SoilGrids REST API (~2-3s per point)
  - Requires local HWSD database (~2-3 GB)
  - Database: https://www.fao.org/documents/card/en/c/CA12305EN/

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
```bash
# 1. Create free CDS account at https://cds.climate.copernicus.eu/
# 2. Create ~/.cdsapi file with credentials:
cat > ~/.cdsapi << EOF
url: https://cds.climate.copernicus.eu/api/v2
key: XXXX:YYYYYYYYYYYYYY
EOF
chmod 600 ~/.cdsapi

# 3. Install cdsapi
pip install cdsapi
```

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
