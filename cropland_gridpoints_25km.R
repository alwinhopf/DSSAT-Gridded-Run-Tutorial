# =============================================================================
# cropland_gridpoints.R — build CROPLAND-ONLY grid points for the pipeline.
#
# Generates (or samples) a grid over a set of US states and keeps only the cells
# that contain agricultural cropland (NLCD class 82 = Cultivated Crops), storing
# the per-cell cropland FRACTION. The output point shapefile feeds the pipeline
# in MODE B (use_existing_point_shapefile: true) so DSSAT runs only on cropland;
# the fraction CSV feeds dssat_to_lca.py (--crop-fraction-csv) so the LCA's
# cropland_ha reflects the real cropland area per cell.
#
# LANDCOVER SOURCE (NLCD via the MRLC WCS, auto-downloaded per cell — no manual
# raster needed). The whole-region 30 m raster is too large for one request, so
# each cell's small window is fetched on demand (cached, parallel). A local NLCD
# GeoTIFF (NLCD_RASTER) is used instead when provided (fully offline).
#
# USAGE
#   Rscript cropland_gridpoints.R
# Configure the block below (states, spacing, cropland classes, WCS year).
# =============================================================================
suppressWarnings(suppressMessages({ library(sf); library(terra) }))

# ---- Configuration --------------------------------------------------------
STATES          <- c("Texas","Louisiana","Mississippi","Alabama","Georgia",
                     "Florida","South Carolina")          # whole southern US incl. TX
GRID_SPACING_M  <- 25000
BOUNDARY_SHP    <- file.path("shapefile", "tl_2024_us_state.shp")
BOUNDARY_COL    <- "NAME"
CROPLAND_CLASSES<- c(82)                                   # NLCD 82 = Cultivated Crops
MIN_CROP_FRAC   <- 0.01      # keep cells with >= 1% cropland
PROJECT_TAG     <- "carinata_southus"
# If set, sample cropland fraction AT these existing points (e.g. a pipeline grid)
# instead of generating a fresh grid — lets you mask an already-built grid.
EXISTING_POINTS <- ""        # e.g. "gridpoints/carinata_southus_full_80km.shp"
# Landcover source: "" => auto-download NLCD per cell from the MRLC WCS;
# otherwise a path to a local NLCD GeoTIFF (offline, fast).
NLCD_RASTER     <- ""
WCS_BASE        <- "https://www.mrlc.gov/geoserver/ows"
WCS_COVERAGE    <- "mrlc_download__NLCD_2021_Land_Cover_L48"   # see WCS GetCapabilities for years
CACHE_DIR       <- file.path("data", "landcover", "cropland_cache")
N_CORES         <- max(1, parallel::detectCores() - 1)

RES_TAG <- if (GRID_SPACING_M < 1000) paste0(GRID_SPACING_M, "m") else paste0(GRID_SPACING_M/1000, "km")
OUT_SHP <- file.path("gridpoints", paste0("cropland_", PROJECT_TAG, "_", RES_TAG, ".shp"))
OUT_CSV <- file.path("gridpoints", paste0("cropland_", PROJECT_TAG, "_", RES_TAG, "_fractions.csv"))
dir.create("gridpoints", showWarnings = FALSE); dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Candidate cell centres (existing points OR a fresh grid) ----------
if (nzchar(EXISTING_POINTS) && file.exists(EXISTING_POINTS)) {
  pts0 <- st_transform(st_read(EXISTING_POINTS, quiet = TRUE), 4326)
  cc   <- st_coordinates(st_centroid(st_geometry(pts0)))
  cand <- data.frame(cell_id = sprintf("%08d", seq_len(nrow(cc))), lat = cc[, 2], lon = cc[, 1])
  cat(sprintf("Sampling cropland at %d existing points from %s\n", nrow(cand), basename(EXISTING_POINTS)))
} else {
  bnd <- st_read(BOUNDARY_SHP, quiet = TRUE); bnd <- bnd[bnd[[BOUNDARY_COL]] %in% STATES, ]
  if (!nrow(bnd)) stop("No states matched in boundary shapefile.")
  albers <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +datum=NAD83 +units=m"
  region_m <- st_transform(st_union(st_transform(bnd, 4326)), albers)
  grid <- st_make_grid(region_m, cellsize = GRID_SPACING_M, what = "centers")
  grid <- grid[st_intersects(st_sf(geometry = grid), region_m, sparse = FALSE)[, 1]]
  cc <- st_coordinates(st_transform(st_sf(geometry = grid), 4326))
  cand <- data.frame(cell_id = sprintf("%08d", seq_len(nrow(cc))), lat = cc[, 2], lon = cc[, 1])
  cat(sprintf("Region: %s | %d candidate cells at %s spacing\n", paste(STATES, collapse=", "), nrow(cand), RES_TAG))
}

# half-cell window (degrees) around each point for the cropland fraction
half_lat <- (GRID_SPACING_M / 2) / 111320
half_lon <- (GRID_SPACING_M / 2) / (111320 * cos(mean(cand$lat) * pi / 180))

# ---- 2. Per-cell cropland fraction ---------------------------------------
wcs_url <- function(latmin, latmax, lonmin, lonmax) sprintf(
  "%s?service=WCS&version=2.0.1&request=GetCoverage&coverageId=%s&subsettingCrs=http://www.opengis.net/def/crs/EPSG/0/4326&subset=Lat(%f,%f)&subset=Long(%f,%f)&format=image/tiff",
  WCS_BASE, WCS_COVERAGE, latmin, latmax, lonmin, lonmax)

frac_one <- function(row) {
  cf <- file.path(CACHE_DIR, paste0(row$cell_id, ".rds"))
  if (file.exists(cf)) return(readRDS(cf))
  tif <- file.path(CACHE_DIR, paste0(row$cell_id, ".tif"))
  val <- tryCatch({
    if (!file.exists(tif)) {
      u <- wcs_url(row$lat - half_lat, row$lat + half_lat, row$lon - half_lon, row$lon + half_lon)
      utils::download.file(u, tif, mode = "wb", quiet = TRUE, method = "libcurl")
    }
    v <- terra::values(terra::rast(tif)); v <- v[!is.na(v)]
    if (length(v)) mean(v %in% CROPLAND_CLASSES) else NA_real_
  }, error = function(e) NA_real_)
  if (file.exists(tif)) unlink(tif)         # keep only the cached fraction, not the tile
  saveRDS(val, cf); val
}

if (nzchar(NLCD_RASTER) && file.exists(NLCD_RASTER)) {
  cat(sprintf("Using local NLCD raster: %s\n", NLCD_RASTER))
  r <- terra::rast(NLCD_RASTER)
  cv <- terra::project(terra::vect(cand[, c("lon","lat")], geom = c("lon","lat"), crs = "EPSG:4326"), terra::crs(r))
  buf <- terra::buffer(cv, GRID_SPACING_M / 2)
  fr <- terra::extract(r, buf, fun = function(x) mean(x %in% CROPLAND_CLASSES, na.rm = TRUE))
  cand$crop_frac <- as.numeric(fr[, 2])
} else {
  cat(sprintf("Auto-downloading NLCD per cell from MRLC WCS [%s] (%d cores, cached)...\n", WCS_COVERAGE, N_CORES))
  rows <- split(cand, seq_len(nrow(cand)))
  cl <- parallel::makeCluster(N_CORES)
  parallel::clusterExport(cl, c("CROPLAND_CLASSES","CACHE_DIR","WCS_BASE","WCS_COVERAGE",
                                "half_lat","half_lon","wcs_url"), envir = environment())
  parallel::clusterEvalQ(cl, suppressMessages(library(terra)))
  cand$crop_frac <- as.numeric(unlist(parallel::parLapply(cl, rows, frac_one)))
  parallel::stopCluster(cl)
}

# ---- 3. Filter to cropland + write ---------------------------------------
crop <- cand[!is.na(cand$crop_frac) & cand$crop_frac >= MIN_CROP_FRAC, ]
pts  <- st_as_sf(data.frame(ID = sprintf("%08d", seq_len(nrow(crop))),
                            LAT = round(crop$lat, 6), LONG = round(crop$lon, 6),
                            crop_pct = round(100 * crop$crop_frac, 2)),
                 coords = c("LONG", "LAT"), crs = 4326, remove = FALSE)
st_write(pts, OUT_SHP, append = FALSE, delete_layer = TRUE, quiet = TRUE)
write.csv(st_drop_geometry(pts), OUT_CSV, row.names = FALSE)
cat(sprintf("\nKept %d/%d cells as cropland (>=%.0f%%). Cropland %%: median %.1f, max %.1f\n",
            nrow(crop), nrow(cand), 100*MIN_CROP_FRAC,
            ifelse(nrow(crop)>0, median(pts$crop_pct), NA), ifelse(nrow(crop)>0, max(pts$crop_pct), NA)))
cat(sprintf("Wrote %s and %s\n", OUT_SHP, OUT_CSV))
