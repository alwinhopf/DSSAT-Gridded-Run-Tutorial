"""One-point real DSSAT smoke test for machines with a verified installation."""

import csv
import os
import platform
import shutil
import subprocess
from pathlib import Path

import pytest


pytestmark = [pytest.mark.integration, pytest.mark.slow]

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "dssat_model_smoke"


def _find_dssat_executable():
    configured = os.getenv("DSSAT_EXE", "").strip()
    names = ["DSCSM048.EXE"] if platform.system() == "Windows" else ["dscsm048"]
    candidates = [Path(configured).expanduser()] if configured else []
    candidates.extend(ROOT.parent.joinpath("DSSAT48", name) for name in names)
    for name in names:
        found = shutil.which(name)
        if found:
            candidates.append(Path(found))
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    return None


def test_one_point_real_dssat_run(tmp_path):
    """Launch DSSAT and assert deterministic yield from the committed fixture."""
    executable = _find_dssat_executable()
    if executable is None:
        pytest.skip("verified DSSAT executable not installed; set DSSAT_EXE")

    profile_names = ("DSSATPRO.V48", "DSSATPRO.L48")
    if not any((executable.parent / name).is_file() for name in profile_names):
        pytest.skip(f"DSSAT support profile missing next to {executable}")

    run_dir = tmp_path / "point"
    shutil.copytree(FIXTURE, run_dir)
    result = subprocess.run(
        [str(executable), "B", "DSSBatch.V48"],
        cwd=run_dir,
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, (
        f"DSSAT exited {result.returncode}\nstdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )

    summary = run_dir / "summary.csv"
    assert summary.is_file(), "DSSAT completed without summary.csv"
    with summary.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1
    assert rows[0]["CR"].strip() == "MZ"
    assert int(float(rows[0]["HWAM"])) == 6292
    assert int(float(rows[0]["CWAM"])) == 15971
