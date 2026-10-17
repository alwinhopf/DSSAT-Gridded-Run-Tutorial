"""One-point real DSSAT smoke test for machines with a verified installation."""

import csv
import os
import platform
import shutil
import subprocess
from pathlib import Path

import pytest


pytestmark = [pytest.mark.integration, pytest.mark.slow]

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "dssat_model_smoke"


def _find_dssat_executable():
    configured = os.getenv("DSSAT_EXE", "").strip()
    names = ["DSCSM048.EXE"] if platform.system() == "Windows" else ["dscsm048"]
    candidates = [Path(configured).expanduser()] if configured else []
    candidates.extend(ROOT.parent.joinpath("DSSAT48", name) for name in names)
    for name in names:
        found = shutil.which(name)
        if found:
            candidates.append(Path(found))
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    return None


def test_one_point_real_dssat_run(tmp_path):
    """Launch DSSAT and assert deterministic yield from the committed fixture."""
    executable = _find_dssat_executable()
    if executable is None:
        pytest.skip("verified DSSAT executable not installed; set DSSAT_EXE")

    profile_names = ("DSSATPRO.V48", "DSSATPRO.L48")
    if not any((executable.parent / name).is_file() for name in profile_names):
        pytest.skip(f"DSSAT support profile missing next to {executable}")

    run_dir = tmp_path / "point"
    shutil.copytree(FIXTURE, run_dir)
    result = subprocess.run(
        [str(executable), "B", "DSSBatch.V48"],
        cwd=run_dir,
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, (
        f"DSSAT exited {result.returncode}\nstdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )

    summary = run_dir / "summary.csv"
    assert summary.is_file(), "DSSAT completed without summary.csv"
    with summary.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1
    assert rows[0]["CR"].strip() == "MZ"
    assert int(float(rows[0]["HWAM"])) == 6292
    assert int(float(rows[0]["CWAM"])) == 15971


@pytest.mark.parametrize('language', ['python', 'r'])
def test_complete_cached_driver_and_receipt(tmp_path, language):
    """Run grid loading -> cached forcing -> model -> completion in both twins."""
    executable = _find_dssat_executable()
    if executable is None:
        pytest.skip('set DSSAT_EXE to a verified installation')
    import sys
    import json
    import yaml
    import geopandas as gpd
    from shapely.geometry import Point
    if language == 'r':
        from tests.helpers.discovery import find_rscript
        rscript = find_rscript()
        if not rscript:
            pytest.skip('Rscript unavailable')
    points = tmp_path/'points.gpkg'
    gpd.GeoDataFrame({'ID':['00000001']}, geometry=[Point(-92.9372,43.1027)], crs=4326).to_file(points)
    weather = tmp_path/'weather/audit_20km_gridpoints_weather'
    soil = tmp_path/'soil/audit_20km_SOILGRIDS_10K_individual_SOL'
    weather.mkdir(parents=True); soil.mkdir(parents=True)
    # The direct-model fixture has one missing winter SRAD value. Supply a
    # synthetic value in this temporary QA fixture so strict input checks pass.
    forcing = (FIXTURE/'00000001.WTH').read_text()
    assert '1992357   -99' in forcing
    (weather/'00000001.WTH').write_text(forcing.replace('1992357   -99', '1992357   3.0'))
    import pandas as pd
    from cropmodel_data.soil_soilgrids_online import _format_dssat_sol_file
    _format_dssat_sol_file(pd.DataFrame(dict(ID=['00000001']*3,latitude=[43.1027]*3,
        longitude=[-92.9372]*3,depth_bottom=[20,60,150],depth_center=[10,40,105],
        sand=[40]*3,clay=[20]*3,silt=[40]*3,cfvo=[0]*3,bdod=[1.3]*3,
        soc_pct=[1]*3,SLLL=[.12]*3,SDUL=[.25]*3,SSAT=[.45]*3)),str(soil))
    (soil.parent/'audit_20km_SOILGRIDS_10K.CSV').write_text('ID,SOIL_ID\n00000001,00000001\n')
    cfg = dict(project_name='audit',output_root_dir=tmp_path.as_posix(),
               use_existing_point_shapefile=True,existing_point_shapefile_path=points.as_posix(),
               weather_source='gridpoints_weather',soil_source='SOILGRIDS_10K',
               run_step_1_soils=False,run_step_2_weather=False,template_dir=FIXTURE.as_posix(),
               template_file_name='00000001.MZX',treatment_start=1,treatment_end=1,
               weather_start_year=1992,weather_end_year=1992,dssat_cores=1,
               soil_cores=1,weather_cores=1,resume_dssat_runs=True,cleanup_run_folders=True)
    overlay=tmp_path/'driver.yml';overlay.write_text(yaml.safe_dump(cfg))
    receipt=tmp_path/'receipt.json'
    env=dict(os.environ,DSSAT_CONFIG_FILE=str(overlay),DSSAT_EXE=str(executable),
             DSSAT_RESULT_RECEIPT=str(receipt),PYTHONUTF8='1')
    command = ([sys.executable,str(ROOT/'dssat_main_pipeline.py')] if language=='python' else
               [rscript,'--vanilla',str(ROOT/'dssat_main_pipeline.R')])
    for attempt in range(2):
        if receipt.exists(): receipt.unlink()
        result=subprocess.run(command,cwd=ROOT,env=env,capture_output=True,text=True,timeout=180)
        assert result.returncode==0, result.stdout[-6000:]+result.stderr[-6000:]
        actual=Path(json.loads(receipt.read_text())['result_path'])
        status=json.loads(Path(str(actual)+'.status.json').read_text())
        assert status['status']=='complete' and status['completed_ids']==['00000001']
        with actual.open() as handle:
            rows=list(csv.DictReader(handle))
        assert len(rows)==1
        assert float(rows[0]['final_grain_kg_ha']) > 0
        if not attempt:
            first_yield = rows[0]['final_grain_kg_ha']
        else:
            assert rows[0]['final_grain_kg_ha'] == first_yield
        if attempt:
            assert ('0 point(s)' in result.stdout+result.stderr or
                    'All simulations already completed' in result.stdout+result.stderr)
