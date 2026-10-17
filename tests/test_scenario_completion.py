import importlib.util
import json
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('scenario_completion',ROOT/'python_scripts/scenario_completion.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)


def test_partial_not_complete_and_verified_exclusion_is_explicit(tmp_path):
    path=tmp_path/'results.csv'
    data=pd.DataFrame({'point_id':['00000001','00000001'],'treatment':[1,2]})
    data.to_csv(path,index=False)
    status=c.publish_completion(path,data,['00000001','00000002'],[1,2])
    assert status['status']=='partial' and status['failed_ids']==['00000002']
    assert not c.scenario_is_complete(path)
    c.publish_completion(path,data,['00000001','00000002'],[1,2],['00000002'])
    assert c.scenario_is_complete(path)
    data.iloc[:1].to_csv(path,index=False)
    assert not c.scenario_is_complete(path)


def test_resume_recovers_success_after_cleanup_only_matching_identity(tmp_path):
    path=tmp_path/'results.csv';run=tmp_path/'runs'
    data=pd.DataFrame({'point_id':['00000001','00000001','00000002'],'treatment':[1,2,1],'yield':[100,200,300]})
    data.to_csv(path,index=False)
    Path(str(path)+'.manifest.json').write_text(json.dumps({'config':123}))
    assert c.restore_aggregate_points(path,run,{'config':456},[1,2],['00000001','00000002'])==set()
    assert c.restore_aggregate_points(path,run,{'config':123},[1,2],['00000001','00000002'])=={'00000001'}
    assert len(pd.read_csv(run/'00000001/results_00000001.csv'))==2
    assert not (run/'00000002').exists()


def test_added_missing_inputs_preserve_unchanged_points_after_cleanup(tmp_path):
    path=tmp_path/'results.csv'
    data=pd.DataFrame({'point_id':['00000001','00000001'],'treatment':[1,2]})
    data.to_csv(path,index=False)
    old={'model':'same','resolved_inputs':{'weather_wth_sha256':'oldcollection','soil_sol_sha256':'oldcollection','dssat_executable_sha256':'same'}}
    new={'model':'same','resolved_inputs':{'weather_wth_sha256':'newcollection','soil_sol_sha256':'newcollection','dssat_executable_sha256':'same'}}
    signatures={'00000001':{'soil':'unchanged','weather':'unchanged'}}
    Path(str(path)+'.manifest.json').write_text(json.dumps(old))
    c.publish_completion(path,data,['00000001','00000002'],[1,2],identity=old,input_signatures=signatures)
    c.start_completion(path)
    assert not c.scenario_is_complete(path)
    assert c.restore_aggregate_points(path,tmp_path/'runs',new,[1,2],['00000001','00000002'],signatures)=={'00000001'}
    assert c.restore_aggregate_points(path,tmp_path/'changed',new,[1,2],['00000001'],{'00000001':{'soil':'changed','weather':'unchanged'}})==set()


def test_cached_input_initialization_publishes_partial_in_both_languages(tmp_path):
    import ast,os,textwrap,subprocess,shutil
    source=(ROOT/'dssat_main_pipeline.py').read_text(encoding='utf-8')
    env=dict(os=os,SOIL_ROOT_DIR=str(tmp_path),SOIL_BASENAME='test',weather_dir=str(tmp_path),
             load_soil_exclusions=lambda path:set(),load_weather_exclusions=lambda path:set())
    markers=[('# Completion exclusions exist even when providers are disabled', '    if RUN_STEP_1_SOILS:'),
             ('# Completion exclusions exist even when weather download is disabled.', '    if RUN_STEP_2_WEATHER:')]
    for start_marker,end_marker in markers:
        start=source.rfind("\n",0,source.index(start_marker))+1;end=source.index(end_marker,start)
        exec(textwrap.dedent(source[start:end]),env)
    data=pd.DataFrame({'point_id':['00000001'],'treatment':[1]})
    path=tmp_path/'results.csv';data.to_csv(path,index=False)
    status=c.publish_completion(path,data,['00000001','00000002'],[1],
                                 env['unresolvable_soil_ids']|env['unresolvable_weather_ids'])
    assert status['status']=='partial' and status['failed_ids']==['00000002']
    assert not c.scenario_is_complete(path)
    rscript=shutil.which('Rscript')
    if not rscript:return
    rsource=(ROOT/'dssat_main_pipeline.R').read_text(encoding='utf-8')
    blocks=[]
    for marker,end_marker in [('# Completion exclusions exist even when providers are disabled','if (RUN_STEP_1_SOILS)'),
                              ('# Completion exclusions exist even when weather download is disabled.','if (RUN_STEP_2_WEATHER)')]:
        start=rsource.index(marker);end=rsource.index(end_marker,start);blocks.append(rsource[start:end])
    code='args<-commandArgs(TRUE)\nCENTRAL_SOIL_DIR<-args[1];SOIL_BASENAME<-"test";output_dir<-args[1]\nload_soil_exclusions<-function(path) character();load_weather_exclusions<-function(path) character()\n'+'\n'.join(blocks)+'\nsource(args[2]);path<-file.path(args[1],"native-results.csv")\ndata<-data.frame(point_id="00000001",treatment=1);write.csv(data,path,row.names=FALSE)\nstatus<-publish_completion(path,data,c("00000001","00000002"),1,union(unresolvable_soil_ids,unresolvable_weather_ids))\nstopifnot(status$status=="partial",identical(unlist(status$failed_ids),"00000002"),!scenario_is_complete(path))\n'
    script=tmp_path/'native.R';script.write_text(code,encoding='utf-8')
    result=subprocess.run([rscript,'--vanilla',str(script),str(tmp_path),str(ROOT/'r_scripts/scenario_completion.R')],capture_output=True,text=True,timeout=90)
    assert result.returncode==0,result.stdout+result.stderr
