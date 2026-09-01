"""Exercise both driver preflight functions without launching the pipeline."""
import ast
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Optional

import pytest

ROOT = Path(__file__).resolve().parents[1]
UTILS = ROOT.parent / "dssatutils"
if not UTILS.is_dir():
    UTILS = ROOT / ".ci-deps" / "dssatutils"


@pytest.mark.parametrize("language", ["python", "r"])
def test_driver_rejects_shifted_soil_cache(tmp_path, language):
    if not (UTILS / "python" / "dssatutils" / "soil_validation.py").is_file():
        pytest.skip("updated shared dssatutils source checkout needed")
    sys.path.insert(0, str(UTILS / "python"))
    from dssatutils import rebuild_soil_files_from_mapping, soil_file_issue
    records = rebuild_soil_files_from_mapping(
        UTILS / "tests/fixtures/soil_mapping_rebuild.csv", tmp_path / "soil", "SSURGO")
    point = tmp_path / "00000001"
    point.mkdir()
    soil = point / "SOIL.SOL"
    shutil.copyfile(records.path.iloc[0], soil)

    if language == "python":
        tree = ast.parse((ROOT / "dssat_main_pipeline.py").read_text())
        node = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "soil_input_issue")
        namespace = dict(os=os, Optional=Optional, soil_file_issue=soil_file_issue, DSSAT_RUN_DIR=str(tmp_path))
        exec(compile(ast.Module(body=[node], type_ignores=[]), "preflight", "exec"), namespace)
        check = lambda: namespace["soil_input_issue"]("00000001") or "OK"
    else:
        rscript = shutil.which("Rscript")
        if not rscript:
            pytest.skip("Rscript unavailable")
        code = '''args <- commandArgs(TRUE)
pkgload::load_all(args[1], quiet=TRUE)
expressions <- parse(args[2])
env <- new.env()
for (expr in expressions) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("soil_input_issue"))) eval(expr, env)
}
result <- env$soil_input_issue("00000001")
cat(if (is.null(result)) "OK" else result)
'''
        def check():
            proc = subprocess.run([rscript, "--vanilla", "-e", code, str(UTILS),
                                   str(ROOT / "dssat_main_pipeline.R")], cwd=tmp_path,
                                  capture_output=True, text=True, timeout=60)
            assert proc.returncode == 0, proc.stderr
            return proc.stdout.strip()
    assert check() == "OK"
    lines = soil.read_text().splitlines()
    h = next(i for i, row in enumerate(lines) if row.startswith("@  SLB"))
    lines[h + 2] = " " + lines[h + 2]
    soil.write_text("\n".join(lines))
    assert "fixed-width" in check()
