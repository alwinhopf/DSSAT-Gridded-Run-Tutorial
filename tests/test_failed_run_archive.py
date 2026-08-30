import ast
import os
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _load_archive_helper(namespace):
    source = (ROOT / "dssat_main_pipeline.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    helper = next(
        node for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "archive_failed_run_folders"
    )
    module = ast.Module(body=[helper], type_ignores=[])
    ast.fix_missing_locations(module)
    exec(compile(module, str(ROOT / "dssat_main_pipeline.py"), "exec"), namespace)
    return namespace["archive_failed_run_folders"]


def test_failed_point_folders_are_archived_with_diagnostics(tmp_path):
    run_dir = tmp_path / "runs" / "scenario_cache_key"
    results_dir = tmp_path / "results"
    archive_path = results_dir / "scenario_failed_run_folders.zip"
    for point_id in ("00000001", "00000002", "00000003"):
        (run_dir / point_id).mkdir(parents=True)
    (run_dir / "00000001" / "_run_error.log").write_text("fatal\n", encoding="utf-8")
    (run_dir / "00000001" / "WARNING.OUT").write_text("warning\n", encoding="utf-8")
    (run_dir / "00000002" / "SOIL.SOL").write_text("soil\n", encoding="utf-8")
    (run_dir / "00000003" / "results_00000003.csv").write_text("ok\n", encoding="utf-8")

    namespace = {
        "os": os,
        "zipfile": zipfile,
        "DSSAT_RUN_DIR": str(run_dir),
        "FINAL_OUTPUT_DIR": str(results_dir),
        "FAILED_RUN_ARCHIVE_PATH": str(archive_path),
    }
    archive = _load_archive_helper(namespace)

    assert archive(["00000001", "00000002"])
    assert archive_path.is_file()
    with zipfile.ZipFile(archive_path) as zf:
        names = set(zf.namelist())
        assert zf.read("FAILED_IDS.txt").decode().splitlines() == ["00000001", "00000002"]
        assert "00000001/_run_error.log" in names
        assert "00000001/WARNING.OUT" in names
        assert "00000002/SOIL.SOL" in names
        assert not any(name.startswith("00000003/") for name in names)

    # A clean rerun must not leave an obsolete failure ZIP in results/.
    assert archive([])
    assert not archive_path.exists()
