import ast
import json
import math
import os
import shutil
import subprocess
from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import Point, box, mapping


ROOT = Path(__file__).resolve().parents[1]


def _python_mask_function(tmp_path):
    tree = ast.parse((ROOT / "dssat_main_pipeline.py").read_text(encoding="utf-8"))
    fn = next(
        node for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == "apply_cropland_mask"
    )
    namespace = {
        "os": os, "math": math, "np": np, "gpd": gpd, "Point": Point,
        "box": box, "mapping": mapping,
        "_resolve_optional_path": lambda value: str(value),
        "POINT_SHAPEFILE_PATH": str(tmp_path / "python_mask.shp"),
    }
    exec(compile(ast.Module(body=[fn], type_ignores=[]), str(ROOT / "dssat_main_pipeline.py"), "exec"), namespace)
    return namespace["apply_cropland_mask"]


def _fixture(tmp_path):
    raster = tmp_path / "crop.tif"
    values = np.zeros((4, 4), dtype=np.uint8)
    values[:, 0] = 82
    with rasterio.open(
        raster, "w", driver="GTiff", height=4, width=4, count=1,
        dtype=values.dtype, crs="EPSG:3857", transform=from_origin(0, 4, 1, 1),
        nodata=255,
    ) as dst:
        dst.write(values, 1)
    points = gpd.GeoDataFrame(
        {"ID": ["00000001"], "LAT": [0.0], "LONG": [0.0]},
        geometry=[Point(2, 2)], crs="EPSG:3857",
    )
    boundary = gpd.GeoDataFrame(
        {"name": ["half-cell"]}, geometry=[box(0, 0, 2, 4)], crs="EPSG:3857"
    )
    return raster, points, boundary


def test_boundary_clipping_and_anchor_relocation_python(tmp_path):
    raster, points, boundary = _fixture(tmp_path)
    apply_mask = _python_mask_function(tmp_path)
    clipped = apply_mask(
        points, raster, 4, [82], filter_basis="cell_fraction",
        boundary_gdf=boundary, clip_to_boundary=True, relocate_anchor=True,
        strict=True,
    )
    assert len(clipped) == 1
    assert clipped.iloc[0].crop_frac == pytest.approx(0.5)
    assert clipped.iloc[0].full_ha == pytest.approx(16 / 10000)
    assert clipped.iloc[0].cell_ha == pytest.approx(8 / 10000, abs=1e-4)
    assert clipped.iloc[0].crop_ha == pytest.approx(4 / 10000, abs=1e-4)
    assert clipped.iloc[0].anchor_mv == 1
    assert clipped.geometry.iloc[0].x == pytest.approx(0.5)
    assert 1.5 <= clipped.geometry.iloc[0].y <= 2.5

    complete = apply_mask(
        points, raster, 4, [82], filter_basis="cell_fraction",
        boundary_gdf=boundary, clip_to_boundary=False, relocate_anchor=False,
        strict=True,
    )
    assert complete.iloc[0].crop_frac == pytest.approx(0.25)
    assert complete.iloc[0].cell_ha == pytest.approx(16 / 10000, abs=1e-4)
    assert complete.iloc[0].anchor_mv == 0


def test_boundary_clipping_and_anchor_relocation_r_matches_python(tmp_path):
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is not installed")
    probe = subprocess.run(
        [rscript, "--vanilla", "-e", "quit(status=if(requireNamespace('sf',quietly=TRUE)&&requireNamespace('terra',quietly=TRUE))0 else 1)"],
        check=False,
    )
    if probe.returncode:
        pytest.skip("R sf/terra packages are not installed")

    raster, points, boundary = _fixture(tmp_path)
    apply_mask = _python_mask_function(tmp_path)
    py = apply_mask(
        points, raster, 4, [82], filter_basis="cell_fraction",
        boundary_gdf=boundary, clip_to_boundary=True, relocate_anchor=True,
        strict=True,
    ).iloc[0]
    r_code = f"""
    library(sf)
    library(terra)
    parsed <- parse(file={json.dumps(str(ROOT / 'dssat_main_pipeline.R'))})
    hit <- which(vapply(parsed, function(x) grepl('^apply_cropland_mask <- function', paste(deparse(x), collapse=' ')), logical(1)))
    eval(parsed[[hit]])
    resolve_optional_path <- function(path) path
    POINT_SHAPEFILE_PATH <- {json.dumps(str(tmp_path / 'r_mask.shp'))}
    pts <- st_sf(ID='00000001', LAT=0, LONG=0, geometry=st_sfc(st_point(c(2,2)), crs=3857))
    boundary <- st_sf(name='half-cell', geometry=st_sfc(st_polygon(list(matrix(c(0,0,2,0,2,4,0,4,0,0), ncol=2, byrow=TRUE))), crs=3857))
    out <- apply_cropland_mask(pts, {json.dumps(str(raster))}, 4, c(82), strict=TRUE,
      filter_basis='cell_fraction', boundary_sf=boundary,
      clip_to_boundary=TRUE, relocate_anchor=TRUE)
    xy <- st_coordinates(out)
    cat(jsonlite::toJSON(list(crop_frac=out$crop_frac[[1]], full_ha=out$full_ha[[1]],
      cell_ha=out$cell_ha[[1]], crop_ha=out$crop_ha[[1]], anchor_mv=out$anchor_mv[[1]],
      anchor_km=out$anchor_km[[1]], x=xy[1,1], y=xy[1,2]), auto_unbox=TRUE))
    """
    result = subprocess.run(
        [rscript, "--vanilla", "-e", r_code], cwd=ROOT,
        check=False, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    got = json.loads(result.stdout[result.stdout.index("{"):])
    assert got["crop_frac"] == pytest.approx(float(py.crop_frac))
    assert got["full_ha"] == pytest.approx(float(py.full_ha))
    assert got["cell_ha"] == pytest.approx(float(py.cell_ha), abs=1e-4)
    assert got["crop_ha"] == pytest.approx(float(py.crop_ha), abs=1e-4)
    assert got["anchor_mv"] == int(py.anchor_mv)
    assert got["x"] == pytest.approx(py.geometry.x)
    assert got["y"] == pytest.approx(py.geometry.y)
