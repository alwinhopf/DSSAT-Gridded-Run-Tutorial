"""Test that clear_run_diagnostics unlinks results_<ID>.csv and run logs in both Python and R."""
import ast
import os
from pathlib import Path
import shutil
import subprocess
import pytest

from tests.helpers.discovery import find_rscript

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("language", ["python", "r"])
def test_clear_run_diagnostics_unlinks_results_and_logs(tmp_path, language):
    point_dir = tmp_path / "00000001"
    point_dir.mkdir(parents=True)

    artifacts = [
        "results_00000001.csv",
        "_run_error.log",
        "dssat_A_stdout_stderr.log",
        "dssat_B_stdout_stderr.log",
        "dssat_Q_stdout_stderr.log",
        "ERROR.OUT",
        "WARNING.OUT",
        "INFO.OUT",
        "Summary.OUT",
        "summary.csv",
    ]
    preserved = ["SOIL.SOL", "00000001.WTH", "SOIL_ID.EXP"]

    for art in artifacts:
        (point_dir / art).write_text("dummy artifact content", encoding="utf-8")
    for pres in preserved:
        (point_dir / pres).write_text("keep this input", encoding="utf-8")

    if language == "python":
        tree = ast.parse((ROOT / "dssat_main_pipeline.py").read_text(encoding="utf-8"))
        node = next(
            n for n in ast.walk(tree)
            if isinstance(n, ast.FunctionDef) and n.name == "clear_run_diagnostics"
        )
        namespace = {"os": os, "DSSAT_RUN_DIR": str(tmp_path)}
        exec(compile(ast.Module(body=[node], type_ignores=[]), "diag", "exec"), namespace)
        namespace["clear_run_diagnostics"](["00000001"])
    else:
        rscript = find_rscript()
        if not rscript:
            pytest.skip("Rscript unavailable")
        code = '''args <- commandArgs(TRUE)
expressions <- parse(args[1])
env <- new.env()
for (expr in expressions) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("clear_run_diagnostics"))) eval(expr, env)
}
env$DSSAT_RUN_DIR <- args[2]
env$clear_run_diagnostics(c("00000001"))
'''
        import tempfile
        with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False) as tf:
            tf.write(code)
            tf_path = tf.name
        try:
            proc = subprocess.run(
                [rscript, "--vanilla", tf_path, str(ROOT / "dssat_main_pipeline.R"), str(tmp_path)],
                cwd=tmp_path, capture_output=True, text=True, timeout=60,
            )
            assert proc.returncode == 0, proc.stderr
        finally:
            if os.path.exists(tf_path):
                os.unlink(tf_path)

    for art in artifacts:
        assert not (point_dir / art).exists(), f"{art} should have been unlinked"
    for pres in preserved:
        assert (point_dir / pres).exists(), f"{pres} should have been preserved"
