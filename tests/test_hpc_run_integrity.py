"""Exercise the actual single-rank runner with isolated model subprocesses."""
import csv
import importlib.util
import json
from pathlib import Path
import runpy
import subprocess
import sys
import types

import pytest

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / 'hpc/dssat_mpi_runner.py'


class SingleRank:
    def Get_rank(self): return 0
    def Get_size(self): return 1
    def Barrier(self): pass
    def bcast(self, value, root=0): return value
    def gather(self, value, root=0): return [value]


@pytest.mark.parametrize('successful', [False, True])
@pytest.mark.parametrize('merge_mode', ['concat', 'none'])
def test_main_isolates_old_parts_and_reports_failures(tmp_path, monkeypatch, successful, merge_mode):
    base = tmp_path / 'runs'
    point = base / '00000001'
    point.mkdir(parents=True)
    (point / '00000001.MZX').write_text('fixture')
    summary = tmp_path / 'summary'
    old = summary / 'runs' / 'temp_results_rank_9.csv'
    old.parent.mkdir(parents=True)
    old.write_text('point_id,final_grain_kg_ha\n99999999,123\n')
    (base / 'hpc_manifest.json').write_text(json.dumps(dict(
        treatments=[1, 3], run_mode='experiment',
        point_metadata={'00000001': {'weight_ha': 2.5}})))
    monkeypatch.setitem(sys.modules, 'mpi4py', types.SimpleNamespace(MPI=types.SimpleNamespace(COMM_WORLD=SingleRank())))
    monkeypatch.setattr(sys, 'argv', [str(RUNNER), '--base_dir', str(base), '--summary_dir', str(summary),
        '--exe_path', sys.executable, '--scratch_dir', str(tmp_path/'scratch'), '--merge_mode', merge_mode])
    def model(*args, **kwargs):
        if not successful:
            raise subprocess.CalledProcessError(1, args[0], stderr='model failed')
        # Batch treatment is the second token on the last line.
        batch = (point/'DSSBatch.V48').read_text().splitlines()[-1].split()
        treatment = int(batch[1])
        (point/'summary.csv').write_text(f'RUNNO,TRNO,PDAT,HWAM\n1,{treatment},1992001,100\n')
        return subprocess.CompletedProcess(args[0], 0, '', '')
    monkeypatch.setattr(subprocess, 'run', model)
    with pytest.raises(SystemExit) as exc:
        runpy.run_path(str(RUNNER), run_name='__main__')
    assert exc.value.code == (0 if successful else 1)
    final = summary/'results_runs.csv'
    status = json.loads(Path(str(final)+'.status.json').read_text())
    assert status['status'] == ('complete' if successful else 'partial')
    assert status['completed_ids'] == (['00000001'] if successful else [])
    assert status['failed_ids'] == ([] if successful else ['00000001'])
    assert status['treatments'] == [1, 3]
    assert str(old) not in status['parts'] and old.exists()
    assert final.exists() == (merge_mode == 'concat')
    paths = [final] if final.exists() else map(Path, status['parts'])
    rows = [r for path in paths for r in csv.DictReader(path.open())]
    assert {r['point_id'] for r in rows} == ({'00000001'} if successful else set())
    if successful:
        assert {float(r['treatment']) for r in rows} == {1, 3}
        assert all(float(r['weight_ha']) == 2.5 for r in rows)


@pytest.mark.parametrize('filename', ['plantgro.csv', 'PlantGro.OUT'])
def test_stress_from_previous_run_cannot_leak(tmp_path, monkeypatch, filename):
    spec = importlib.util.spec_from_file_location('runner_integrity', RUNNER)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    (tmp_path/'00000001.MZX').write_text('fixture')
    old = tmp_path/filename
    old.write_text('RUN,WSPD,WSGD,NSTD\n1,0.8,0.7,0.6\n')
    def model(*args, **kwargs):
        assert not old.exists()
        (tmp_path/'summary.csv').write_text('RUNNO,TRNO,PDAT,HWAM\n1,1,1992001,100\n')
        return subprocess.CompletedProcess(args[0], 0, '', '')
    monkeypatch.setattr(subprocess, 'run', model)
    result = runner.run_dssat_simulation_orchestrator(str(tmp_path), '00000001', 1, dict(
        run_mode='experiment', executable_path=sys.executable, control_file_name='DSSBatch.V48',
        dssat_model_code='MZCER048', sequence_start=1, sequence_end=1, cleanup_mode='never'))
    assert result[3] is None
    assert result[2][0]['water_stress_photosynthesis_avg'] is None


def test_staging_failure_is_not_a_success(tmp_path, monkeypatch):
    base = tmp_path/'runs'; base.mkdir()
    monkeypatch.setitem(sys.modules, 'mpi4py', types.SimpleNamespace(MPI=types.SimpleNamespace(COMM_WORLD=SingleRank())))
    monkeypatch.setattr(sys, 'argv', [str(RUNNER), '--base_dir', str(base), '--summary_dir', str(tmp_path/'out'),
        '--exe_path', sys.executable, '--scratch_dir', str(tmp_path/'scratch')])
    import shutil
    def fail(*args, **kwargs): raise OSError('staging fixture failure')
    monkeypatch.setattr(shutil, 'move', fail)
    with pytest.raises(SystemExit) as exc:
        runpy.run_path(str(RUNNER), run_name='__main__')
    assert exc.value.code == 1
    status = json.loads((tmp_path/'out/results_runs.csv.status.json').read_text())
    assert status['status'] == 'partial'
    assert any('staging fixture failure' in error for error in status['errors'])
