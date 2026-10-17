"""Consumer import contracts with legacy support namespaces unavailable."""
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def test_driver_loads_monorepo_helpers_without_legacy_packages(tmp_path):
    overlay = tmp_path / 'config.yml'
    overlay.write_text(
        f'output_root_dir: "{tmp_path.as_posix()}"\n'
        f'input_root_dir: "{tmp_path.as_posix()}"\n'
        f'template_dir: "{(ROOT / "dssat_templates").as_posix()}"\n'
        f'prism_cache_dir: "{tmp_path.as_posix()}/prism"\n'
        'soil_source: SSURGO\nuse_existing_point_shapefile: true\n',
        encoding='utf-8',
    )
    code = r'''
import importlib.abc, runpy, sys
class BlockLegacy(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in ('dssatutils', 'dssatengine', 'dssatcalibrator'):
            raise AssertionError('Retired support dependency: ' + fullname)
sys.meta_path.insert(0, BlockLegacy())
driver = runpy.run_path('dssat_main_pipeline.py', run_name='consumer_contract')
assert driver['process_weather_daymet'].__module__.startswith('cropmodel_data.')
assert driver['soil_file_issue'].__module__.startswith('cropmodel_data.')
assert isinstance(driver['CHECK_SOIL_DOWNLOADS'], bool)
assert callable(driver['clean_invalid_soils'])
assert driver['_run_simulation'].__module__.startswith('cropmodel_runtime.')
assert driver['is_wth_valid'](sys.argv[1], 1992, start_year=1992)
'''
    weather = ROOT / 'tests/fixtures/dssat_model_smoke/00000001.WTH'
    env = dict(os.environ, DSSAT_CONFIG_FILE=str(overlay), PYTHONUTF8='1')
    result = subprocess.run([sys.executable, '-c', code, str(weather)], cwd=ROOT,
                            env=env, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
