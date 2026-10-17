"""Run actual SoilGrids selection/retry code against historical failure records."""
import ast
import os
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest
from tests.helpers.discovery import find_rscript

ROOT = Path(__file__).resolve().parents[1]
from tests.helpers.packages import package_source
DATA = package_source('cropmodel_data')


@pytest.mark.parametrize('language', ['python', 'r'])
@pytest.mark.parametrize('check_downloads', [False, True])
@pytest.mark.parametrize('driver_root', [ROOT])
def test_stale_skipped_log_does_not_suppress_batch_retry(tmp_path, language, check_downloads, driver_root):
    if not driver_root.is_dir():
        pytest.skip('mirrored workflow unavailable')
    from cropmodel_data.soil_soilgrids_online import _format_dssat_sol_file
    from cropmodel_data.soil_validation import soil_file_issue
    _format_dssat_sol_file(pd.DataFrame(dict(ID=['A'], latitude=[32], longitude=[-90],
        depth_bottom=[5], depth_center=[2.5], sand=[40], clay=[20], silt=[40], cfvo=[0],
        bdod=[1.3], soc_pct=[1], SLLL=[.12], SDUL=[.25], SSAT=[.45])), str(tmp_path))
    good = (tmp_path / 'A.SOL').read_text()
    (tmp_path / 'soil_processing_errors.log').write_text('SKIPPED: ID: A | critical physical data contains NAs\n')
    (tmp_path / 'B.SOL').write_text(good.replace('  0.0', '   NA'))
    assert soil_file_issue(tmp_path / 'B.SOL') is not None
    if language == 'python':
        tree = ast.parse((driver_root / 'dssat_main_pipeline.py').read_text(encoding='utf-8'))
        branch = next(n for n in ast.walk(tree) if isinstance(n, ast.If)
                      and ast.unparse(n.test) == "SOIL_SOURCE == 'SOILGRIDS_ONLINE'")
        calls = []
        def provider(gridfile, **kwargs):
            calls.append(gridfile.ID.tolist())
            if len(calls) == 1:
                raise RuntimeError('required layer failed; HTTP 503')
            (tmp_path / 'B.SOL').write_text(good)
        env = dict(os=os, pd=pd, gridfile=pd.DataFrame({'ID': ['A', 'B']}),
                   POINT_ID_COLUMN='ID', individual_sol_dir=str(tmp_path),
                   soilfile_CSV=str(tmp_path / 'map.csv'), soil_unresolvable_file='',
                   unresolvable_soil_ids=set(), SOILGRIDS_MODE='VRT', SOILGRIDS_CACHE_DIR='',
                   cfg_get=lambda key, default=None: check_downloads if key == 'check_soil_downloads' else 2,
                   soil_file_issue=soil_file_issue, process_soils_soilgrids_online=provider,
                   remove_unresolvable_point_ids=lambda *args: None,
                   save_unresolvable_point_ids=lambda *args: pytest.fail('recovered points must not be excluded'))
        settings = [n for n in tree.body if isinstance(n, ast.Assign) and
                    any(isinstance(t, ast.Name) and t.id in ('SOIL_DOWNLOAD_RETRIES', 'CHECK_SOIL_DOWNLOADS') for t in n.targets)]
        cleaner = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'clean_invalid_soils')
        exec(compile(ast.Module(body=settings + [cleaner], type_ignores=[]), 'soil-config', 'exec'), env)
        assert env['CHECK_SOIL_DOWNLOADS'] is check_downloads
        exec(compile(ast.Module(body=branch.body, type_ignores=[]), 'soil-branch', 'exec'), env)
        assert calls == [['B'], ['B']]
    else:
        rscript = find_rscript()
        if not rscript:
            pytest.skip('Rscript unavailable')
        script = tmp_path / 'soil.R'
        script.write_text('''a<-commandArgs(TRUE)
pkgload::load_all(a[3],quiet=TRUE)
find_branch <- function(x) {
  if(missing(x)||!is.call(x)) return(NULL)
  if(identical(x[[1]],as.name('if')) && identical(paste(deparse(x[[2]]),collapse=' '),'SOIL_SOURCE == "SOILGRIDS_ONLINE"')) return(x[[3]])
  for(child in as.list(x)[-1]) {z<-find_branch(child);if(!is.null(z))return(z)}
  NULL
}
branch<-Filter(Negate(is.null),lapply(as.list(parse(a[1])),find_branch))[[1]]
cfg_get<-function(key,default=NULL) if(key=="check_soil_downloads") a[4]=="true" else 2
for(expr in as.list(parse(a[1]))) {
 if(is.call(expr)&&identical(expr[[1]],as.name("<-"))&&is.symbol(expr[[2]])&&as.character(expr[[2]]) %in%
    c("clean_invalid_soils","CHECK_SOIL_DOWNLOADS","SOIL_DOWNLOAD_RETRIES")) eval(expr)
}
gridfile<-data.frame(ID=c('A','B'),LAT=c(32,33));POINT_ID_COLUMN<-'ID'
individual_sol_output_folder<-a[2]
soilfile_CSV<-file.path(a[2],'map.csv');soil_unresolvable_file<-file.path(a[2],'old.json')
unresolvable_soil_ids<-character();SOILGRIDS_MODE<-'VRT';SOILGRIDS_CACHE_DIR<-''
calls<-list()
save_unresolvable_point_ids<-function(...) stop('recovered points must not be excluded')
process_soils_soilgrids_online<-function(gridfile,...) {
 calls[[length(calls)+1]]<<-as.character(gridfile$ID)
 if(length(calls)==1)stop('required layer failed; HTTP 503')
 file.copy(file.path(a[2],'A.SOL'),file.path(a[2],'B.SOL'),overwrite=TRUE)
}
eval(branch)
stopifnot(identical(calls,list('B','B')))
''', encoding='utf-8')
        proc = subprocess.run([rscript, '--vanilla', str(script), str(driver_root / 'dssat_main_pipeline.R'),
                               str(tmp_path), str(DATA / 'R/cropmodeldata'), str(check_downloads).lower()],
                              capture_output=True, text=True, timeout=60)
        assert proc.returncode == 0, proc.stdout + proc.stderr
    assert (tmp_path / 'A.SOL').read_text() == good
    assert soil_file_issue(tmp_path / 'B.SOL') is None
