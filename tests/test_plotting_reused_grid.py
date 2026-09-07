"""Exercise real plotting blocks without running downloads or simulations."""
import ast
from pathlib import Path
import subprocess
import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point
ROOT = Path(__file__).resolve().parents[1]

@pytest.mark.parametrize('language', ['R', 'python'])
@pytest.mark.parametrize('backdrop', ['available', 'missing', 'unreadable', 'existing_points'])
def test_reused_grid_renders(language, backdrop, tmp_path):
    data = tmp_path / 'results.csv'
    pd.DataFrame({'point_id':['00000001'], 'latitude':[30.], 'longitude':[-85.],
                  'treatment':[4], 'final_grain_kg_ha':[1200.]}).to_csv(data,index=False)
    boundary = tmp_path / 'boundary.geojson'
    if backdrop in ('available', 'existing_points'):
        gpd.GeoDataFrame({'NAME':['keep','exclude']},geometry=[Point(-85,30),Point(0,0)],crs=4326).to_file(boundary,driver='GeoJSON')
    elif backdrop == 'unreadable': boundary.write_text('invalid geometry file')
    plot = tmp_path / 'map.png'
    env=dict(RUN_DSSAT_EXECUTION=True,FINAL_RESULTS_PATH=str(data),FINAL_OUTPUT_DIR=str(tmp_path),
             FINAL_PLOT_PATH=str(plot),DSSAT_RUN_NAME='test',boundary_sf=None,
             USE_EXISTING_POINT_SHAPEFILE=backdrop=='existing_points',SHAPEFILE_DIR=str(tmp_path),
             BOUNDARY_SHAPEFILE_NAME=boundary.name,ENABLE_BOUNDARY_FILTER=True,
             BOUNDARY_FILTER_COLUMN='NAME',BOUNDARY_FILTER_VALUE=['keep'],
             WEATHER_SOURCE='TEST',WEATHER_START_YEAR=1984,WEATHER_END_YEAR=2025)
    if language == 'python':
        import matplotlib
        matplotlib.use('Agg')
        tree=ast.parse((ROOT/'dssat_main_pipeline.py').read_text())
        block=max((n for n in ast.walk(tree) if isinstance(n,ast.If) and 'STEP 4: VISUALIZING RESULTS' in ast.unparse(n)), key=lambda n:n.lineno)
        env.update(pd=pd,gpd=gpd,os=__import__('os'))
        exec(compile(ast.Module(body=[block],type_ignores=[]),'<plotting>','exec'),env)
        assert (env['plot_boundary'] is not None)==(backdrop=='available')
        if backdrop=='available': assert len(env['plot_boundary'])==1
    else:
        import json, shutil
        if not shutil.which('Rscript'): pytest.skip('Rscript unavailable')
        def literal(x):
            if x is None: return 'NULL'
            if isinstance(x,bool): return 'TRUE' if x else 'FALSE'
            if isinstance(x,list): return 'c('+','.join(map(literal,x))+')'
            return json.dumps(x)
        setup='\n'.join(k+' <- '+literal(v) for k,v in env.items())+'\n'
        block=(ROOT/'dssat_main_pipeline.R').read_text().split('# STEP 4: VISUALIZE RESULTS',1)[1]
        script=tmp_path/'render.R'
        script.write_text(setup+block+'\nstopifnot(is.null(boundary_sf_4326) == '+literal(backdrop!='available')+')\n')
        completed=subprocess.run(['Rscript','--vanilla',str(script)],capture_output=True,text=True,timeout=60)
        assert completed.returncode==0,completed.stderr
        plot=tmp_path/'map_treatment4.png'
    assert plot.exists() and plot.stat().st_size>1000
