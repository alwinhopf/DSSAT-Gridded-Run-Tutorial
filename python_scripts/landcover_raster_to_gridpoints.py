# File: landcover_raster_to_gridpoints.py
# Python port of landcover_raster_to_gridpoints.R
#
# Converts a binary cropland raster (0/1) into a gridded point shapefile
# where each point represents one aggregated grid cell with a cropland
# fraction.  The output feeds directly into dssat_main_pipeline.py via
# MODE B (USE_EXISTING_POINT_SHAPEFILE = True).
#
# TYPICAL WORKFLOW
#   Step 1: landcover_raster.py  →  data/landcover/derived/cropland_mask_*.tif
#   Step 2: THIS script           →  gridpoints/<region>_cropland_<res>.shp
#   Step 3: dssat_main_pipeline.py (MODE B, EXISTING_POINT_SHAPEFILE_PATH = …)

import os
import sys

import numpy as np

try:
    import rasterio
    from rasterio.transform import from_origin
    from rasterio.crs import CRS
    import geopandas as gpd
    from shapely.geometry import Point
    from pyproj import Transformer
except ImportError as exc:
    sys.exit(
        f"Missing required packages: {exc}\n"
        "Install with: pip install rasterio geopandas shapely pyproj"
    )


def raster_to_cropland_points(
    crop_raster_file: str,
    output_dir: str,
    output_shapefile_all: str,
    output_shapefile_filtered: str,
    grid_resolution_m: int = 5000,
    cropland_threshold: float = 0.0,
) -> None:
    """
    Aggregate a binary cropland raster to *grid_resolution_m* spacing, convert
    cells to centroid points, attach cropland fraction, and export two
    shapefiles (all points and filtered points above *cropland_threshold*).
    Mirrors the R ``landcover_raster_to_gridpoints.R`` script exactly.

    Parameters
    ----------
    crop_raster_file : str
        Path to the binary cropland raster (0 = non-cropland, 1 = cropland).
        Must be in a projected (metre-based) CRS.
    output_dir : str
        Directory where output shapefiles are written.
    output_shapefile_all : str
        Path for the shapefile containing ALL grid cells (even 0 % cropland).
    output_shapefile_filtered : str
        Path for the shapefile containing only cells above *cropland_threshold*.
    grid_resolution_m : int
        Target grid spacing in metres.  Must match GRID_SPACING_METERS in the
        main pipeline for consistent point density.
    cropland_threshold : float
        Minimum cropland proportion [0, 1) to retain a point in the filtered
        output.  0 = keep any cell with > 0 cropland pixels.
    """
    os.makedirs(output_dir, exist_ok=True)

    print("Loading cropland raster...")
    with rasterio.open(crop_raster_file) as src:
        crs_obj = src.crs
        transform = src.transform
        data = src.read(1).astype(float)
        nodata = src.nodata
        res_x = abs(transform.a)  # pixel width in CRS units
        res_y = abs(transform.e)  # pixel height in CRS units

        # Check that CRS is projected (not geographic degrees)
        if crs_obj.is_geographic:
            raise ValueError(
                "The cropland raster appears to be lon/lat (degrees). "
                "Reproject to a projected CRS in metres before calling this function."
            )

    # Replace nodata with NaN
    if nodata is not None:
        data[data == nodata] = np.nan

    # -----------------------------------------------------------------------
    # Aggregation factor (round to nearest integer)
    # -----------------------------------------------------------------------
    agg_factor = max(1, round(grid_resolution_m / res_x))
    print(f"Aggregating raster (factor = {agg_factor})...")

    n_rows, n_cols = data.shape
    new_rows = n_rows // agg_factor
    new_cols = n_cols // agg_factor

    # Trim to exact multiple
    data_trim = data[:new_rows * agg_factor, :new_cols * agg_factor]

    # Block-mean aggregation (equivalent to terra::aggregate with fun="mean")
    agg_data = (
        data_trim
        .reshape(new_rows, agg_factor, new_cols, agg_factor)
        .mean(axis=(1, 3))
    )  # shape: (new_rows, new_cols)

    # -----------------------------------------------------------------------
    # Compute cell centre coordinates in the raster CRS
    # -----------------------------------------------------------------------
    # New transform after aggregation
    new_res = res_x * agg_factor
    orig_left = transform.c       # x of left edge
    orig_top  = transform.f       # y of top edge
    agg_transform = rasterio.transform.Affine(
        new_res, 0, orig_left,
        0, -new_res, orig_top
    )

    row_indices, col_indices = np.where(~np.isnan(agg_data))
    # Rasterio's xy() returns (x, y) for row/col centre
    xs, ys = rasterio.transform.xy(agg_transform, row_indices, col_indices)
    xs = np.array(xs)
    ys = np.array(ys)
    crop_props = agg_data[row_indices, col_indices]

    print("Adding ID, latitude, and longitude columns...")

    # Re-project centres to WGS84 to get lon/lat
    transformer = Transformer.from_crs(crs_obj.to_string(), "EPSG:4326", always_xy=True)
    lons, lats = transformer.transform(xs, ys)

    n_pts = len(xs)
    ids = [f"{i+1:08d}" for i in range(n_pts)]

    # -----------------------------------------------------------------------
    # Build GeoDataFrame in native projected CRS
    # -----------------------------------------------------------------------
    print("Converting aggregated cells to points...")
    geometries = [Point(x, y) for x, y in zip(xs, ys)]

    gdf = gpd.GeoDataFrame(
        {
            "ID":        ids,
            "crop_prop": crop_props,
            "crop_pct":  crop_props * 100.0,
            "LONG":      lons,
            "LAT":       lats,
        },
        geometry=geometries,
        crs=crs_obj,
    )
    # Reorder columns to match R output
    gdf = gdf[["ID", "crop_prop", "crop_pct", "LONG", "LAT", "geometry"]]

    # -----------------------------------------------------------------------
    # Write ALL points
    # -----------------------------------------------------------------------
    gdf.to_file(output_shapefile_all)
    print(f"Saved ALL grid points to: {output_shapefile_all}")

    # -----------------------------------------------------------------------
    # Filter and write FILTERED points
    # -----------------------------------------------------------------------
    print(f"Filtering points by cropland_threshold > {cropland_threshold}...")
    filtered = gdf[gdf["crop_prop"] > cropland_threshold].copy()
    filtered.to_file(output_shapefile_filtered)
    print(f"Saved FILTERED points to: {output_shapefile_filtered}")
    print("Done.")


# ---------------------------------------------------------------------------
# Script entry point (run directly from project root)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # -------------------------
    # USER SETTINGS
    # -------------------------
    CROP_RASTER_FILE         = os.path.join("data", "landcover", "derived",
                                            "cropland_mask_montana.tif")
    OUTPUT_DIR               = "gridpoints"
    OUTPUT_SHAPEFILE_ALL     = os.path.join(OUTPUT_DIR, "montana_cropland_5k.shp")
    OUTPUT_SHAPEFILE_FILTERED = os.path.join(OUTPUT_DIR, "montana_cropland_filtered_5k.shp")
    GRID_RESOLUTION_M        = 5000
    CROPLAND_THRESHOLD       = 0.0   # > 0 = any cropland pixel; > 0.05 = ≥5 %
    # -------------------------

    raster_to_cropland_points(
        crop_raster_file=CROP_RASTER_FILE,
        output_dir=OUTPUT_DIR,
        output_shapefile_all=OUTPUT_SHAPEFILE_ALL,
        output_shapefile_filtered=OUTPUT_SHAPEFILE_FILTERED,
        grid_resolution_m=GRID_RESOLUTION_M,
        cropland_threshold=CROPLAND_THRESHOLD,
    )
