# =============================================================================
# cropland_gridpoints.R — build CROPLAND-ONLY grid points for the pipeline.
#
# Generates a regular grid over a set of US states, then keeps only the cells
# that contain agricultural cropland (NLCD class 82 = Cultivated Crops), storing
# the per-cell cropland FRACTION. The output point shapefile feeds the pipeline
# in MODE B (use_existing_point_shapefile: true), so DSSAT runs only on cropland;
# the fraction CSV feeds dssat_to_lca.py (--crop-fraction-csv) so the LCA's
# cropland_ha reflects the real cropland area per cell (not the whole cell).
#
# LANDCOVER SOURCE (in priority order):
#   1. NLCD_RASTER  — a local national/clipped NLCD GeoTIFF (fast, offline).
#                     Download once: https://www.mrlc.gov/data  (Annual NLCD).
#   2. FedData live — get_nlcd() per cell (no download, but slow; parallel+cached).
#
# USAGE
#   Rscript cropland_gridpoints.R
# Configure the block below (states, spacing, raster path, cropland classes).
# =============================================================================
suppressWarnings(suppressMessages({ library(sf); library(terra) }))

# ---- Configuration --------------------------------------------------------
STATES          <- c("Texas","Louisiana","Mississippi","Alabama","Georgia",
                     "Florida","South Carolina")          # whole southern US incl. TX
GRID_SPACING_M  <- 55000
BOUNDARY_SHP    <- file.path("shapefile", "tl_2024_us_state.shp")
BOUNDARY_COL    <- "NAME"
CROPLAND_CLASSES<- c(82)                                   # NLCD 82 = Cultivated Crops
NLCD_RASTER     <- ""        # path to a local NLCD GeoTIFF; "" => live FedData per cell
NLCD_YEAR       <- 2021
MIN_CROP_FRAC   <- 0.01      # keep cells with >= 1% cropland
PROJECT_TAG     <- "carinata_southus"
CACHE_DIR       <- file.path("data", "landcover", "cropland_cache")
N_CORES         <- max(1, parallel::detectCores() - 1)

RES_TAG <- if (GRID_SPACING_M < 1000) paste0(GRID_SPACING_M, "m") else paste0(GRID_SPACING_M/1000, "km")
OUT_SHP <- file.path("gridpoints", paste0("cropland_", PROJECT_TAG, "_", RES_TAG, ".shp"))
OUT_CSV <- file.path("gridpoints", paste0("cropland_", PROJECT_TAG, "_", RES_TAG, "_fractions.csv"))
dir.create("gridpoints", showWarnings = FALSE); dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Region + candidate grid ------------------------------------------
bnd <- st_read(BOUNDARY_SHP, quiet = TRUE)
bnd <- bnd[bnd[[BOUNDARY_COL]] %in% STATES, ]
if (!nrow(bnd)) stop("No states matched in boundary shapefile.")
region <- st_union(st_transform(bnd, 4326))
# project to an equal-distance CRS for a regular metric grid (US Albers)
albers <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +datum=NAD83 +units=m"
region_m <- st_transform(region, albers)
grid <- st_make_grid(region_m, cellsize = GRID_SPACING_M, what = "polygons")
grid <- grid[st_intersects(grid, region_m, sparse = FALSE)[, 1]]
cells <- st_transform(st_sf(geometry = grid), 4326)
cells$cell_id <- sprintf("%08d", seq_len(nrow(cells)))
cat(sprintf("Region: %s | %d candidate cells at %s spacing\n",
            paste(STATES, collapse = ", "), nrow(cells), RES_TAG))

# ---- 2. Per-cell cropland fraction ---------------------------------------
crop_fraction_local <- function(cells_sf, rast) {
  r <- terra::rast(rast)
  cv <- terra::project(terra::vect(cells_sf), terra::crs(r))
  fr <- terra::extract(r, cv, fun = function(x) mean(x %in% CROPLAND_CLASSES, na.rm = TRUE))
  as.numeric(fr[, 2])
}
crop_fraction_fed <- function(cell_sf) {
  cid <- cell_sf$cell_id
  cf <- file.path(CACHE_DIR, paste0(cid, ".rds"))
  if (file.exists(cf)) return(readRDS(cf))
  frac <- tryCatch({
    suppressMessages(nl <- FedData::get_nlcd(template = cell_sf, label = cid, year = NLCD_YEAR,
                                             extraction.dir = CACHE_DIR, force.redo = FALSE))
    v <- terra::values(nl); v <- v[!is.na(v)]
    if (length(v)) mean(v %in% CROPLAND_CLASSES) else NA_real_
  }, error = function(e) NA_real_)
  saveRDS(frac, cf); frac
}

if (nzchar(NLCD_RASTER) && file.exists(NLCD_RASTER)) {
  cat(sprintf("Using local NLCD raster: %s\n", NLCD_RASTER))
  cells$crop_frac <- crop_fraction_local(cells, NLCD_RASTER)
} else {
  if (!requireNamespace("FedData", quietly = TRUE)) stop("Set NLCD_RASTER, or install FedData for live fetch.")
  cat(sprintf("No local raster — fetching NLCD per cell via FedData (%d cores, cached). This can be slow.\n", N_CORES))
  cl <- parallel::makeCluster(N_CORES)
  parallel::clusterExport(cl, c("CROPLAND_CLASSES","NLCD_YEAR","CACHE_DIR"), envir = environment())
  fracs <- parallel::parLapply(cl, split(cells, seq_len(nrow(cells))), crop_fraction_fed)
  parallel::stopCluster(cl)
  cells$crop_frac <- as.numeric(unlist(fracs))
}

# ---- 3. Filter to cropland + write ---------------------------------------
crop <- cells[!is.na(cells$crop_frac) & cells$crop_frac >= MIN_CROP_FRAC, ]
ctr  <- st_centroid(st_geometry(crop))
cc   <- st_coordinates(ctr)
pts  <- st_sf(ID = sprintf("%08d", seq_len(nrow(crop))),
              LAT = round(cc[, 2], 6), LONG = round(cc[, 1], 6),
              crop_pct = round(100 * crop$crop_frac, 2),
              geometry = ctr)
st_write(pts, OUT_SHP, append = FALSE, delete_layer = TRUE, quiet = TRUE)
write.csv(st_drop_geometry(pts), OUT_CSV, row.names = FALSE)
cat(sprintf("\nKept %d/%d cells as cropland (>=%.0f%%). Cropland fraction: median %.1f%%, max %.1f%%\n",
            nrow(crop), nrow(cells), 100*MIN_CROP_FRAC, median(pts$crop_pct), max(pts$crop_pct)))
cat(sprintf("Wrote %s and %s\n", OUT_SHP, OUT_CSV))
cat("Next: set use_existing_point_shapefile: true and existing_point_shapefile_path to the shapefile,\n",
    "     then run the pipeline; pass", basename(OUT_CSV), "to dssat_to_lca.py --crop-fraction-csv.\n")
