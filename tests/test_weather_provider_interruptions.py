"""Execute actual pipeline retry loops with partial publication and provider failures."""
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
from cropmodel_data.provider_retry import ProviderConnectivityError


@pytest.mark.parametrize('language', ['python', 'r'])
@pytest.mark.parametrize('mode', ['recover', 'wrapped', 'exhaust', 'fatal'])
def test_interrupted_provider_rechecks_completed_files(tmp_path, language, mode):
    if language == 'python':
        tree = ast.parse((ROOT / 'dssat_main_pipeline.py').read_text(encoding='utf-8'))
        loop = next(n for n in ast.walk(tree) if isinstance(n, ast.While)
                    and 'points_to_process.empty' in ast.unparse(n.test))
        calls = []
        def provider(points):
            calls.append(points.ID.tolist())
            if mode == 'fatal':
                raise ValueError('invalid configuration')
            (tmp_path / 'A.WTH').write_text('completed')
            if mode == 'wrapped' and len(calls) == 1:
                raise RuntimeError('task 1 failed: Provider connectivity unavailable: SSL connection reset')
            if mode == 'exhaust' or len(calls) == 1:
                raise ProviderConnectivityError('SSL connection reset')
            (tmp_path / 'B.WTH').write_text('completed')
        env = dict(os=os, points_to_process=pd.DataFrame({'ID': ['A', 'B']}),
                   attempt=0, WEATHER_DOWNLOAD_RETRIES=2, WEATHER_SOURCE='NASA_POWER',
                   process_weather_points=provider, POINT_ID_COLUMN='ID',
                   weather_dir=str(tmp_path), CHECK_WEATHER_DOWNLOADS=True,
                   weather_file_is_valid=lambda p: Path(p).exists(),
                   apply_weather_repairs=lambda ids: None,
                   remove_unresolvable_point_ids=lambda *a: None,
                   weather_unresolvable_file='', weather_retryable_file='',
                   unresolvable_weather_ids=set())
        execute = lambda: exec(compile(ast.Module(body=[loop], type_ignores=[]), 'actual-loop', 'exec'), env)
        if mode == 'fatal':
            with pytest.raises(ValueError, match='invalid configuration'):
                execute()
            assert calls == [['A', 'B']]
        else:
            execute()
            assert calls == [['A', 'B'], ['B']]
            assert env['attempt'] == 2
            assert env['points_to_process'].ID.tolist() == (['B'] if mode == 'exhaust' else [])
            assert (tmp_path / 'A.WTH').read_text() == 'completed'
            assert 'SSL connection reset' in (tmp_path / 'download_errors.log').read_text()
    else:
        rscript = find_rscript()
        if not rscript:
            pytest.skip('Rscript unavailable')
        script = tmp_path / 'test.R'
        script.write_text('''a <- commandArgs(TRUE)
pkgload::load_all(a[1], quiet=TRUE)
find_loop <- function(x) {
  if (missing(x)) return(NULL)
  if (!is.call(x)) return(NULL)
  if (identical(x[[1]], as.name('while')) && grepl('points_to_process', paste(deparse(x[[2]]), collapse=' '))) return(x)
  for (child in as.list(x)[-1]) { z <- find_loop(child); if (!is.null(z)) return(z) }
  NULL
}
outer <- Filter(function(x) is.call(x) && identical(x[[1]], as.name('if')) &&
  identical(x[[2]], as.name('RUN_STEP_2_WEATHER')), as.list(parse(a[2])))[[1]]
loops <- list(find_loop(outer))
loop <- Filter(Negate(is.null), loops)[[1]]
mode <- a[4]; calls <- list(); points_to_process <- data.frame(ID=c('A','B'), LAT=c(32,33))
WEATHER_SOURCE <- 'DAYMET'; WEATHER_START_YEAR <- 2001; WEATHER_END_YEAR <- 2001
POINT_ID_COLUMN <- 'ID'; LAT_COLUMN <- 'LAT'; LONG_COLUMN <- 'LONG'; WEATHER_CORES <- 1
output_dir <- a[3]; log_file <- file.path(output_dir,'download_errors.log')
retry_count <- 0; max_retries <- 2; CHECK_WEATHER_DOWNLOADS <- TRUE
unresolvable_weather_ids <- character(); weather_unresolvable_file <- ''; weather_retryable_file <- ''
apply_weather_repairs <- function(...) NULL
remove_unresolvable_point_ids <- function(...) NULL
is_wth_valid <- function(id, dir, ...) file.exists(file.path(dir,paste0(id,'.WTH')))
process_weather_daymet <- function(shapefile, ...) {
  calls[[length(calls)+1]] <<- as.character(shapefile$ID)
  if (mode == 'fatal') stop('invalid configuration')
  writeLines('completed',file.path(output_dir,'A.WTH'))
  if (mode == 'wrapped' && length(calls)==1) stop('task 1 failed: Provider connectivity unavailable: SSL connection reset')
  if (mode == 'exhaust' || length(calls)==1) stop(structure(list(message='SSL connection reset',call=NULL), class=c('dssat_connectivity_error','error','condition')))
  writeLines('completed',file.path(output_dir,'B.WTH'))
}
error <- tryCatch({eval(loop);NULL}, error=function(e)e)
if (mode=='fatal') {
  stopifnot(inherits(error,'error'),grepl('invalid configuration',conditionMessage(error)),length(calls)==1)
} else {
  if (!is.null(error)) stop(error)
  stopifnot(is.null(error),identical(calls,list(c('A','B'),'B')),retry_count==2)
  stopifnot(identical(as.character(points_to_process$ID), if(mode=='exhaust') 'B' else character()))
  stopifnot(readLines(file.path(output_dir,'A.WTH'))=='completed')
  stopifnot(any(grepl('SSL connection reset',readLines(log_file))))
}
''', encoding='utf-8')
        proc = subprocess.run([rscript, '--vanilla', str(script), str(DATA / 'R/cropmodeldata'),
                               str(ROOT / 'dssat_main_pipeline.R'), str(tmp_path), mode],
                              capture_output=True, text=True, timeout=90)
        assert proc.returncode == 0, proc.stdout + proc.stderr
