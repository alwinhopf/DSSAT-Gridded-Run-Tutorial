"""Native R acceptance checks for driver-issued result receipts."""
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def test_receipts_require_complete_results_and_preserve_custom_paths(tmp_path, rscript):
    script = tmp_path/'contract.R'
    script.write_text(r'''
args <- commandArgs(TRUE)
source(file.path(args[1], 'r_scripts/scenario_completion.R'))
source(file.path(args[1], 'r_scripts/experiment_contract.R'))
root <- file.path(args[2], 'custom results'); dir.create(root)
path <- file.path(root, 'study_sampling_custom_results.csv')
receipt <- file.path(args[2], 'receipt.json')
stopifnot(is.null(completed_receipt_path(receipt)))
data <- data.frame(point_id='00000001',treatment=1)
write.csv(data,path,row.names=FALSE)
jsonlite::write_json(list(schema=1L,result_path=path),receipt,auto_unbox=TRUE)
stopifnot(is.null(completed_receipt_path(receipt)))
jsonlite::write_json(list(config='fixture'),paste0(path,'.manifest.json'),auto_unbox=TRUE)
publish_completion(path,data,c('00000001','00000002'),1)
stopifnot(is.null(completed_receipt_path(receipt)))
data <- rbind(data,data.frame(point_id='00000002',treatment=1))
write.csv(data,path,row.names=FALSE)
publish_completion(path,data,c('00000001','00000002'),1)
stopifnot(identical(completed_receipt_path(receipt),path))
# A damaged CSV must invalidate an otherwise complete sidecar.
write.csv(data[1,,drop=FALSE],path,row.names=FALSE)
stopifnot(is.null(completed_receipt_path(receipt)))
''', encoding='utf-8')
    result = subprocess.run([rscript, '--vanilla', str(script), str(ROOT), str(tmp_path)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr


def test_experiment_plan_uses_canonical_defaults_without_running_models(tmp_path, rscript):
    overlay = tmp_path/'experiment.yml'
    overlay.write_text('experiment_name: contract\nbase: {}\nfactors:\n  weather_source: [NASA_POWER]\n'
                       '  soil_source: [SSURGO]\noptions:\n  dry_run: true\n  validate: false\n')
    result = subprocess.run([rscript, '--vanilla', str(ROOT/'run_experiment.R'), str(overlay)],
                            cwd=tmp_path, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert '20km' in result.stdout
    assert '1992_1992' in result.stdout
    assert 'driver-resolved output' in result.stdout


def test_orchestrator_invokes_driver_on_reuse_and_excludes_partial_results(tmp_path, rscript):
    import shutil
    import csv
    project = tmp_path/'consumer'
    project.mkdir()
    (project/'r_scripts').mkdir()
    for name in ['run_experiment.R', 'config.yml']:
        shutil.copyfile(ROOT/name, project/name)
    for name in ['scenario_completion.R', 'experiment_contract.R']:
        shutil.copyfile(ROOT/'r_scripts'/name, project/'r_scripts'/name)
    # A controlled driver exercises the real subprocess/receipt/aggregation path.
    (project/'dssat_main_pipeline.R').write_text(r'''
cfg <- yaml::yaml.load(paste(readLines(Sys.getenv('DSSAT_CONFIG_FILE'),encoding='UTF-8'),collapse='\n'))
stopifnot(isTRUE(cfg$resume_dssat_runs),isTRUE(cfg$check_provenance),cfg$grid_spacing_meters==20000)
source('r_scripts/scenario_completion.R')
cat(cfg$soil_source,'\n',file='invocations.txt',append=TRUE)
dir.create(cfg$results_root_dir,recursive=TRUE,showWarnings=FALSE)
path <- file.path(cfg$results_root_dir,paste0(cfg$soil_source,'_sampling_driver_suffix.csv'))
d <- data.frame(point_id=c('00000001','00000002'),treatment=1,
                latitude=c(40,41),longitude=c(-90,-91),final_grain_kg_ha=c(100,200))
if(cfg$soil_source=='POLARIS') d <- d[1,,drop=FALSE]
write.csv(d,path,row.names=FALSE)
publish_completion(path,d,c('00000001','00000002'),1)
jsonlite::write_json(list(model='fixture'),paste0(path,'.manifest.json'),auto_unbox=TRUE)
jsonlite::write_json(list(schema=1L,result_path=path),Sys.getenv('DSSAT_RESULT_RECEIPT'),auto_unbox=TRUE)
''', encoding='utf-8')
    experiment = project/'experiment.yml'
    experiment.write_text(f'''experiment_name: receipt_contract
base:
  results_root_dir: "{(tmp_path/'custom outputs').as_posix()}"
factors:
  weather_source: [NASA_POWER]
  soil_source: [SSURGO, POLARIS]
options:
  reuse_existing: true
  validate: false
  stop_on_error: false
''')
    for _ in range(2):
        result = subprocess.run([rscript, '--vanilla', str(project/'run_experiment.R'), str(experiment)],
                                cwd=tmp_path, capture_output=True, text=True, timeout=90)
        assert result.returncode == 0, result.stdout + result.stderr
    assert len((project/'invocations.txt').read_text().splitlines()) == 4
    with (project/'results/EXPERIMENT_receipt_contract_combined.csv').open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
    assert {row['soil_source'] for row in rows} == {'SSURGO'}
    assert rows[0]['point_id'] == '00000001'
