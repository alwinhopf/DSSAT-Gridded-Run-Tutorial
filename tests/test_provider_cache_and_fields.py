"""Run the engine's cache binding and FileX substitution blocks offline."""
import ast
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
CACHES = {'GRIDMET_CACHE_DIR': 'gridmet_netcdf_cache',
          'CHIRPS_CACHE_DIR': 'chirps_netcdf_cache',
          'CHIRPS_V3_CACHE_DIR': 'chirps_v3_netcdf_cache',
          'AGERA5_CACHE_DIR': 'agera5_netcdf_cache',
          'DWD_CACHE_DIR': 'dwd_station_cache', 'EOBS_CACHE_DIR': 'eobs_cds_cache'}


@pytest.mark.parametrize('override', ['', 'shared-products'])
def test_provider_cache_paths_match_across_languages(tmp_path, override):
    input_root = tmp_path / 'main-engine'
    output_root = tmp_path / 'study'
    input_root.mkdir(); output_root.mkdir()
    tree = ast.parse((ROOT / 'dssat_main_pipeline.py').read_text())
    selected = [n for n in tree.body if
                (isinstance(n, ast.FunctionDef) and n.name == 'resolve_config_path') or
                (isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and
                 t.id in set(CACHES) | {'PROVIDER_CACHE_ROOT'} for t in n.targets))]
    env = dict(Path=Path, os=os, CODE_ROOT_DIR=str(ROOT), INPUT_ROOT_DIR=str(input_root),
               OUTPUT_ROOT_DIR=str(output_root), cfg_get=lambda k, default: override)
    exec(compile(ast.Module(body=selected, type_ignores=[]), 'cache-binding', 'exec'), env)
    expected_root = input_root / override if override else input_root
    expected = {k: str(expected_root / v) for k, v in CACHES.items()}
    assert {k: env[k] for k in CACHES} == expected
    # A populated cache resolves to the very same file, never a study-local copy.
    for key in CACHES:
        path = Path(expected[key]); path.mkdir(parents=True, exist_ok=True)
        (path / 'cached.nc').write_bytes(b'existing provider product')
        assert (Path(env[key]) / 'cached.nc').read_bytes() == b'existing provider product'
    rscript = shutil.which('Rscript')
    if not rscript:
        pytest.skip('Rscript unavailable')
    code = '''a<-commandArgs(TRUE);CODE_ROOT_DIR<-getwd();INPUT_ROOT_DIR<-a[1];OUTPUT_ROOT_DIR<-a[2]
cfg_get<-function(k,default) a[3]
keys<-c("PROVIDER_CACHE_ROOT","GRIDMET_CACHE_DIR","CHIRPS_CACHE_DIR","CHIRPS_V3_CACHE_DIR","AGERA5_CACHE_DIR","DWD_CACHE_DIR","EOBS_CACHE_DIR")
for(e in parse("dssat_main_pipeline.R")) {
 if(is.call(e)&&identical(e[[1]],as.name("<-"))&&as.character(e[[2]])[1] %in% c(keys,"resolve_config_path")) eval(e)
 if(exists("PROVIDER_CACHE_ROOT")&&!nzchar(PROVIDER_CACHE_ROOT)) PROVIDER_CACHE_ROOT<-INPUT_ROOT_DIR
}
cat(jsonlite::toJSON(mget(keys[-1]),auto_unbox=TRUE))'''
    result = subprocess.run([rscript, '--vanilla', '-e', code, str(input_root), str(output_root), override],
                            cwd=ROOT, capture_output=True, text=True, check=True)
    assert {k: Path(v).resolve() for k, v in json.loads(result.stdout).items()} == {k: Path(v).resolve() for k, v in expected.items()}


def test_actual_coordinate_substitution_and_regional_template(tmp_path):
    template = (ROOT / 'dssat_templates' / 'CARINATA1984.SQX').read_text()
    fields = [line for line in template.splitlines() if 'LATITUDE' in line and 'LONGITUDE' in line]
    assert len(fields) == 2
    for line in fields:
        assert line[2:18].strip() == 'LONGITUDE'
        assert line[18:34].strip() == 'LATITUDE'
    content = '@L ...........XCRD ...........YCRD .....ELEV\n' + '\n'.join(fields)
    tree = ast.parse((ROOT / 'dssat_main_pipeline.py').read_text())
    function = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == '_replace_field_coordinates')
    substitution = next(n for n in ast.walk(tree) if isinstance(n, ast.Assign) and
                        isinstance(n.value, ast.Call) and any(isinstance(a, ast.Name) and
                        a.id == '_replace_field_coordinates' for a in n.value.args))
    env = dict(re=re, content=content, s_lat=f'{31.753:8.3f}', s_lon=f'{-87.508:9.3f}', s_elev='  19')
    exec(compile(ast.Module(body=[function, substitution], type_ignores=[]), 'fields', 'exec'), env)
    for before, after in zip(content.splitlines(), env['content'].splitlines()):
        assert len(before) == len(after)
        if before.startswith('@'):
            assert before == after
        else:
            assert float(after[2:18]) == -87.508
            assert float(after[18:34]) == 31.753
