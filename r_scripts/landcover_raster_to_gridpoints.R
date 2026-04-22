# create_shapefile_from_landcover_raster.R  (repo-ready version)
# ------------------------------------------------------------
# PURPOSE
#   Convert a binary cropland raster (0/1) into a gridded point shapefile
#   where each point represents an aggregated grid cell with a cropland
#   fraction.  The output can be fed directly into dssat_main_pipeline.R
#   via MODE B (USE_EXISTING_POINT_SHAPEFILE = TRUE).
#
# TYPICAL WORKFLOW
#   Step 1: landcover_raster.R
#           Create a binary cropland mask raster for your region.
#           Output: data/landcover/derived/cropland_mask_<region>.tif
#
#   Step 2: THIS script
#           Aggregate the mask to your desired DSSAT grid spacing and
#           export point shapefiles.
#           Output: gridpoints/<region>_cropland_<res>.shp
#
#   Step 3: dssat_main_pipeline.R
#           Set USE_EXISTING_POINT_SHAPEFILE <- TRUE and point
#           EXISTING_POINT_SHAPEFILE_PATH at the filtered shapefile
#           produced by this script.
#
# NOTE
#   This script must be run from the project root (or sourced via
#   dssat_main_pipeline.R which sets the working directory).  All paths
#   below are relative to the project root.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
})

# -----------------------
# USER SETTINGS (edit these)
# -----------------------

# Binary cropland raster produced by landcover_raster.R (0 = non-cropland, 1 = cropland).
# Path is relative to the project root.
crop_raster_file <- file.path("data", "landcover", "derived", "cropland_mask_montana.tif")

# Output directory for the point shapefiles (relative to project root).
output_dir <- "gridpoints"

# Output shapefile names — change the region prefix to match your study area.
output_shapefile_all      <- file.path(output_dir, "montana_cropland_5k.shp")
output_shapefile_filtered <- file.path(output_dir, "montana_cropland_filtered_5k.shp")

# Grid cell size in meters — must match the target GRID_SPACING_METERS in
# dssat_main_pipeline.R so the point density is consistent.
grid_resolution_m <- 5000

# Minimum cropland proportion to keep a point (0 = keep any cell with > 0 % cropland;
# 0.05 = keep only cells where at least 5 % of the area is cropland, etc.).
cropland_threshold <- 0

# -----------------------
# END USER SETTINGS
# -----------------------

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading cropland raster...\n")
crop_raster <- rast(crop_raster_file)

# Sanity check: aggregation by meters only makes sense in projected CRS
if (is.lonlat(crop_raster)) {
  stop("The cropland raster appears to be lon/lat (degrees). Reproject to a projected CRS in meters before aggregating by grid_resolution_m.")
}

# Compute aggregation factor from raster resolution
res_m <- res(crop_raster)[1]
aggregation_factor <- round(grid_resolution_m / res_m)
if (aggregation_factor < 1) aggregation_factor <- 1

cat("Aggregating raster (factor =", aggregation_factor, ")...\n")
aggregated_raster <- aggregate(
  crop_raster,
  fact = aggregation_factor,
  fun = "mean",
  na.rm = TRUE
)

cat("Converting aggregated cells to points...\n")
grid_points <- as.points(aggregated_raster)

# Format attributes
names(grid_points)[1] <- "crop_prop"
grid_points$crop_pct <- grid_points$crop_prop * 100

# Add ID + LAT/LONG
cat("Adding ID, latitude, and longitude columns...\n")
grid_points$ID <- sprintf("%08d", 1:nrow(grid_points))

# Reproject a copy to WGS84 to extract lon/lat
points_latlon <- project(grid_points, "EPSG:4326")
coords <- crds(points_latlon)
grid_points$LONG <- coords[, 1]
grid_points$LAT  <- coords[, 2]

# Reorder columns
grid_points <- grid_points[, c("ID", "crop_prop", "crop_pct", "LONG", "LAT")]

# Write ALL points
writeVector(grid_points, output_shapefile_all, overwrite = TRUE)
cat("Saved ALL grid points to:", output_shapefile_all, "\n")

# Filter points
cat("Filtering points by cropland_threshold >", cropland_threshold, "...\n")
filtered_points <- grid_points[grid_points$crop_prop > cropland_threshold, ]
writeVector(filtered_points, output_shapefile_filtered, overwrite = TRUE)
cat("Saved FILTERED points to:", output_shapefile_filtered, "\n")
cat("Done.\n")
