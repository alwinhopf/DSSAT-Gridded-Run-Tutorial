"""Run the real engine selection blocks without loading or launching the pipeline."""
import ast
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
UTILS = ROOT.parent / 'dssatutils'
if not UTILS.is_dir():
    UTILS = ROOT / '.ci-deps' / 'dssatutils'
if (UTILS / 'python').is_dir():
    sys.path.insert(0, str(UTILS / 'python'))
from tests.helpers.discovery import find_rscript

try:
    from dssatutils.weather_validation import is_wth_valid
except (ImportError, Exception):
    is_wth_valid = None


@pytest.mark.parametrize('language',['python','r'])
def test_real_engine_keeps_exclusions_but_retries_old_failures(tmp_path,language):
    p = tmp_path / 'unresolvable_points.json'
    p.write_text(json.dumps({'00000001':{'reason':'failed_after_1_retries'},
                            '00000002':{'reason':'outside_coverage'},
                            '00000003':{'reason':'failed_after_1_retries'},
                            'unrelated':{'reason':'outside_coverage'}}))
    dates = pd.date_range('2001-01-01','2001-12-31')
    (tmp_path/'00000001.WTH').write_text('\n'.join(f'{d:%Y%j} 12 25 15 2 10 60 3' for d in dates))
    (tmp_path/'00000002.WTH').write_text('invalid excluded evidence')
    if language=='python':
        if is_wth_valid is None:
            pytest.skip("dssatutils.weather_validation is not available")
        tree = ast.parse((ROOT/'dssat_main_pipeline.py').read_text())
        names = {'load_unresolvable_point_ids','load_weather_exclusions',
                 'remove_unresolvable_point_ids','save_unresolvable_point_ids','weather_file_is_valid'}
        nodes=[n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name in names]
        env=dict(os=os,json=json,datetime=datetime,pd=pd,is_wth_valid=is_wth_valid,
                 WEATHER_START_YEAR=2001,WEATHER_END_YEAR=2001,WEATHER_VALIDATION_END_DATE='2001-12-31',
                 weather_required_columns=lambda: ['SRAD','TMAX','TMIN','RAIN','TDEW','RH2M','WIND'])
        exec(compile(ast.Module(body=nodes,type_ignores=[]),'helpers','exec'),env)
        ids=pd.Series(['00000001','00000002','00000003'])
        env.update(weather_ids=ids,gridfile=pd.DataFrame({'ID':ids}),POINT_ID_COLUMN='ID',
                   weather_dir=str(tmp_path),weather_unresolvable_file=str(p),
                   weather_retryable_file=str(tmp_path/'retryable_weather_points.json'),
                   unresolvable_weather_ids=env['load_weather_exclusions'](str(p)),CHECK_WEATHER_DOWNLOADS=True)
        outer=next(n for n in ast.walk(tree) if isinstance(n,ast.If) and ast.unparse(n.test)=='RUN_STEP_2_WEATHER')
        selection=next(n for n in outer.body if isinstance(n,ast.If) and ast.unparse(n.test)=='CHECK_WEATHER_DOWNLOADS')
        exec(compile(ast.Module(body=[selection],type_ignores=[]),'selection','exec'),env)
        assert env['points_to_process'].ID.tolist()==['00000003']
        # Execute actual exhausted-retry handling, preserving unrelated records.
        failed=next(n for n in ast.walk(outer) if isinstance(n,ast.If) and n.body and
                    isinstance(n.body[0],ast.Assign) and isinstance(n.body[0].targets[0],ast.Name)
                    and n.body[0].targets[0].id=='failed_wth_ids')
        env['WEATHER_DOWNLOAD_RETRIES']=1
        exec(compile(ast.Module(body=[failed],type_ignores=[]),'failure','exec'),env)
        assert '00000003' in json.loads((tmp_path/'retryable_weather_points.json').read_text())
    else:
        rscript = find_rscript()
        if not rscript:
            pytest.skip("Rscript unavailable")
        code = '''a<-commandArgs(TRUE)
# Load the tested source revision, even if another version is installed.
pkgload::load_all(a[1], quiet = TRUE)
exprs<-parse(a[2]); env<-new.env()
for(e in exprs) if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]])[1] %in%
 c('load_unresolvable_point_ids','load_weather_exclusions','remove_unresolvable_point_ids','save_unresolvable_point_ids','is_wth_valid','weather_required_columns')) eval(e,env)
env$WEATHER_START_YEAR<-2001;env$WEATHER_END_YEAR<-2001;env$WEATHER_VALIDATION_END_DATE<-as.Date('2001-12-31');env$WEATHER_SOURCE<-'AGERA5'
env$ids<-c('00000001','00000002','00000003');env$gridfile<-data.frame(ID=env$ids);env$POINT_ID_COLUMN<-'ID'
env$output_dir<-a[3];env$weather_unresolvable_file<-file.path(a[3],'unresolvable_points.json');env$weather_retryable_file<-file.path(a[3],'retryable_weather_points.json')
env$unresolvable_weather_ids<-env$load_weather_exclusions(env$weather_unresolvable_file);env$CHECK_WEATHER_DOWNLOADS<-TRUE
outer<-Filter(function(e)is.call(e)&&identical(e[[1]],as.name('if'))&&identical(e[[2]],as.name('RUN_STEP_2_WEATHER')),as.list(exprs))[[1]]
selection<-Filter(function(e)is.call(e)&&identical(e[[1]],as.name('if'))&&identical(e[[2]],as.name('CHECK_WEATHER_DOWNLOADS')),as.list(outer[[3]]))[[1]]
eval(selection,env)
# R selection stores its mask, then selects points immediately afterwards.
selected<-env$ids[env$missing_mask & !(env$ids %in% env$unresolvable_weather_ids)]
stopifnot(identical(selected,'00000003'))
'''
        # Use the resolved executable path. On Windows, setup-r can expose
        # Rscript.exe through a PATH entry that is visible to shutil.which()
        # but not reliably discoverable again from a child process launched
        # with the bare command name.
        proc = subprocess.run([rscript, '--vanilla', '-e', code, str(UTILS), str(ROOT / 'dssat_main_pipeline.R'), str(tmp_path)], capture_output=True, text=True, timeout=60)
        assert proc.returncode == 0, proc.stdout + proc.stderr
    remaining=json.loads(p.read_text())
    assert '00000001' not in remaining
    assert remaining['00000002']['reason']=='outside_coverage'
    assert 'unrelated' in remaining
    assert (tmp_path/'00000002.WTH').read_text()=='invalid excluded evidence'
