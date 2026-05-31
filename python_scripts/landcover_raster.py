# File: landcover_raster.py
# Python port of landcover_raster.R
#
# Clips a national landcover raster (NLCD or CDL) to a set of US states and
# produces a binary cropland mask (1 = cropland, 0 = everything else).
# The mask is the input for landcover_raster_to_gridpoints.py.
#
# Run from the project root:
#   python py_scripts/landcover_raster.py
# Or import and call clip_landcover_to_states() from dssat_main_pipeline.py.
#
# Download NLCD: https://www.mrlc.gov/data
# Download CDL:  https://croplandcros.scinet.usda.gov/

import os
import sys

import numpy as np

try:
    import rasterio
    from rasterio.mask import mask as rio_mask
    from rasterio.transform import from_bounds
    import fiona
    from shapely.geometry import shape, mapping
    from shapely.ops import unary_union
    import geopandas as gpd
except ImportError as exc:
    sys.exit(
        f"Missing required packages: {exc}\n"
        "Install with: pip install rasterio fiona shapely geopandas"
    )

# -----------------------
# STATE LOOKUP TABLE
# -----------------------
_STATE_LOOKUP = {
    "Alabama": "AL", "Alaska": "AK", "Arizona": "AZ", "Arkansas": "AR",
    "California": "CA", "Colorado": "CO", "Connecticut": "CT", "Delaware": "DE",
    "District of Columbia": "DC", "Florida": "FL", "Georgia": "GA", "Hawaii": "HI",
    "Idaho": "ID", "Illinois": "IL", "Indiana": "IN", "Iowa": "IA",
    "Kansas": "KS", "Kentucky": "KY", "Louisiana": "LA", "Maine": "ME",
    "Maryland": "MD", "Massachusetts": "MA", "Michigan": "MI", "Minnesota": "MN",
    "Mississippi": "MS", "Missouri": "MO", "Montana": "MT", "Nebraska": "NE",
    "Nevada": "NV", "New Hampshire": "NH", "New Jersey": "NJ", "New Mexico": "NM",
    "New York": "NY", "North Carolina": "NC", "North Dakota": "ND", "Ohio": "OH",
    "Oklahoma": "OK", "Oregon": "OR", "Pennsylvania": "PA", "Rhode Island": "RI",
    "South Carolina": "SC", "South Dakota": "SD", "Tennessee": "TN", "Texas": "TX",
    "Utah": "UT", "Vermont": "VT", "Virginia": "VA", "Washington": "WA",
    "West Virginia": "WV", "Wisconsin": "WI", "Wyoming": "WY",
}
_ABBR_TO_NAME = {v: k for k, v in _STATE_LOOKUP.items()}


def _resolve_states(state_names: list[str]) -> list[str]:
    """
    Convert a list of state names or 2-letter abbreviations to full names.
    Mirrors ``resolve_states`` in R.
    """
    resolved = []
    for s in state_names:
        s = s.strip()
        if not s:
            continue
        if len(s) == 2:
            full = _ABBR_TO_NAME.get(s.upper())
        else:
            full = s if s in _STATE_LOOKUP else None
        if full is not None:
            resolved.append(full)
    return list(dict.fromkeys(resolved))  # unique, order-preserving


def _project_geometry(geom, src_crs: str, dst_crs: str):
    """Reproject a Shapely geometry from src_crs to dst_crs using pyproj."""
    from pyproj import Transformer
    from shapely.ops import transform as shp_transform
    tr = Transformer.from_crs(src_crs, dst_crs, always_xy=True)
    return shp_transform(tr.transform, geom)


def clip_landcover_to_states(
    input_raster: str,
    boundary_vector: str,
    state_names: list[str],
    cropland_values: list[int],
    output_dir: str,
    output_mask: str,
    write_per_state: bool = True,
) -> None:
    """
    Clip *input_raster* to the union of *state_names*, create a binary
    cropland mask, and write it (plus optional per-state masks) to *output_dir*.
    Mirrors the top-level script logic in the R ``landcover_raster.R``.

    Parameters
    ----------
    input_raster : str
        Path to the raw landcover raster (GeoTIFF).
    boundary_vector : str
        Path to the state boundary shapefile (must have NAME and STUSPS fields).
    state_names : list[str]
        State names or 2-letter abbreviations to process.
    cropland_values : list[int]
        Raster class codes to treat as cropland (all others → 0).
    output_dir : str
        Directory for output masks.
    output_mask : str
        Path for the merged regional cropland mask.
    write_per_state : bool
        If True, also write one mask file per state.
    """
    os.makedirs(output_dir, exist_ok=True)

    allowed_states = _resolve_states(state_names)
    if not allowed_states:
        raise ValueError("No valid states found in state_names.")
    print(f"Processing {len(allowed_states)} state(s): {', '.join(allowed_states)}")

    # Load boundary
    bnd_gdf = gpd.read_file(boundary_vector)
    if not all(c in bnd_gdf.columns for c in ("NAME", "STUSPS")):
        raise ValueError(
            f"Boundary must have NAME and STUSPS columns. "
            f"Found: {list(bnd_gdf.columns)}"
        )

    bnd_gdf = bnd_gdf[bnd_gdf["NAME"].str.strip().isin(allowed_states)].copy()
    if bnd_gdf.empty:
        raise ValueError("No matching features found for the requested state_names.")

    print("Loading landcover raster metadata...")
    with rasterio.open(input_raster) as src:
        raster_crs = src.crs.to_string()
        raster_profile = src.profile.copy()

    # Re-project boundary to raster CRS
    bnd_gdf = bnd_gdf.to_crs(raster_crs)

    def _write_mask(geom_union, out_path: str, tag: str) -> None:
        """Crop, mask, and binarise the raster for a given polygon."""
        geom_json = [mapping(geom_union)]
        with rasterio.open(input_raster) as src:
            try:
                arr, tf = rio_mask(src, geom_json, crop=True, nodata=0,
                                   all_touched=False)
            except Exception as exc:
                print(f"  Skipping {tag}: {exc}")
                return
            meta = src.meta.copy()
            meta.update({
                "driver": "GTiff",
                "height": arr.shape[1],
                "width":  arr.shape[2],
                "transform": tf,
                "dtype": "uint8",
                "nodata": 0,
                "count": 1,
                "compress": "lzw",
            })
            # Binary mask: 1 if pixel value is in cropland_values, else 0
            binary = np.isin(arr[0], cropland_values).astype(np.uint8)
            binary[arr[0] == src.nodata if src.nodata is not None else False] = 0

        with rasterio.open(out_path, "w", **meta) as dst:
            dst.write(binary, 1)
        print(f"  Saved: {out_path}")

    # --- Regional mask (union of all states) ---
    print("Creating regional cropland mask...")
    region_union = unary_union(bnd_gdf.geometry.values)
    _write_mask(region_union, output_mask, "region")
    print(f"Saved cropland mask to: {output_mask}")

    # --- Per-state masks ---
    if write_per_state:
        print("Writing per-state masks...")
        for name in allowed_states:
            subset = bnd_gdf[bnd_gdf["NAME"].str.strip() == name]
            if subset.empty:
                continue
            state_union = unary_union(subset.geometry.values)
            safe_name = name.lower().replace(" ", "_")
            out_state = os.path.join(output_dir, f"cropland_mask_{safe_name}.tif")
            _write_mask(state_union, out_state, name)

    print("Done.")


# ---------------------------------------------------------------------------
# Script entry point (run directly from project root)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # -------------------------
    # USER SETTINGS
    # -------------------------
    INPUT_RASTER = os.path.join(
        "data", "landcover",
        "Annual_NLCD_LndCov_2024_CU_C1V1",
        "Annual_NLCD_LndCov_2024_CU_C1V1.tif"
    )
    BOUNDARY_VECTOR = os.path.join("shapefile", "tl_2024_us_state.shp")
    STATE_NAMES = [
        "Montana", "North Dakota", "South Dakota", "Wyoming",
        "Nebraska", "Kansas", "Colorado", "Oklahoma",
        "Minnesota", "Iowa", "Missouri", "Wisconsin", "Illinois", "Michigan",
        "Indiana", "Ohio", "Idaho", "Utah", "Washington", "Oregon", "Nevada",
    ]
    # NLCD 82 = Cultivated Crops; for CDL use crop-specific codes e.g. [1, 5]
    CROPLAND_VALUES = [82]
    OUTPUT_DIR  = os.path.join("data", "landcover", "derived")
    OUTPUT_MASK = os.path.join(OUTPUT_DIR, "cropland_mask_region.tif")
    WRITE_PER_STATE = True
    # -------------------------

    clip_landcover_to_states(
        input_raster=INPUT_RASTER,
        boundary_vector=BOUNDARY_VECTOR,
        state_names=STATE_NAMES,
        cropland_values=CROPLAND_VALUES,
        output_dir=OUTPUT_DIR,
        output_mask=OUTPUT_MASK,
        write_per_state=WRITE_PER_STATE,
    )
