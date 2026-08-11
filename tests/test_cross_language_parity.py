import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import yaml
import pytest


ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def _dependency_path(repo: str, *parts: str) -> Path:
    """Resolve workspace siblings locally and exact-ref CI checkouts online."""
    candidates = (
        WORKSPACE / repo,
        ROOT / ".ci-deps" / repo,
    )
    for candidate in candidates:
        path = candidate.joinpath(*parts)
        if path.exists():
            return path
    return candidates[0].joinpath(*parts)


def _require_sources(*paths: Path) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if not missing:
        return
    message = "required parity source checkout(s) missing: " + ", ".join(missing)
    if os.getenv("DSSAT_REQUIRE_PARITY_SOURCES", "") == "1":
        pytest.fail(message)
    pytest.skip(message)


def test_ssurgo_failure_diagnostics_are_present_in_r_and_python():
    r_path = _dependency_path("dssatutils", "R", "soil_ssurgo.R")
    py_path = _dependency_path("dssatutils", "python", "dssatutils", "soil_ssurgo.py")
    _require_sources(r_path, py_path)
    r_src = _read(r_path)
    py_src = _read(py_path)

    for marker in ("no-coverage", "no-soil", "no-layers", "_download_failures.csv"):
        assert marker in r_src
        assert marker in py_src


def test_engine_missing_summary_failure_logging_is_present_in_r_and_python():
    r_path = _dependency_path("dssatengine", "R", "engine.R")
    py_path = _dependency_path("dssatengine", "python", "dssatengine", "engine.py")
    _require_sources(r_path, py_path)
    r_src = _read(r_path)
    py_src = _read(py_path)

    for src in (r_src, py_src):
        assert "_run_error.log" in src
        # Robust to wording ("DSSAT produced no ..." / "DSSAT completed but
        # produced no ...") — assert the functional parity, not the exact phrase.
        assert "produced no 'summary.csv'" in src
        assert "FMOPT" in src
        assert "summary.csv" in src


def test_gridmet_netcdf_validation_is_thread_safe_in_pinned_source():
    """Prevent concurrent netCDF4/HDF5 reads from returning to CI."""
    path = _dependency_path(
        "dssatutils", "python", "dssatutils", "weather_gridmet.py"
    )
    _require_sources(path)
    source = _read(path)
    assert "_GRIDMET_NETCDF_LOCK = threading.RLock()" in source
    assert "with _GRIDMET_NETCDF_LOCK:" in source
    assert "download_workers = max(1, min(download_workers, 8))" in source


def test_live_provider_retry_policy_is_present_in_r_and_python():
    """Keep transient GridMET failures bounded and explicit in both suites."""
    r_source = _read(ROOT / "tests" / "test_e2e_comprehensive.R").lower()
    py_source = _read(ROOT / "tests" / "test_e2e_comprehensive.py").lower()
    for source in (r_source, py_source):
        assert "dssat_provider_retries" in source
        assert "invalid/corrupt netcdf cache file" in source
        assert "remained unavailable" in source


def test_filex_coordinate_substitution_policy_is_aligned():
    r_src = _read(ROOT / "dssat_main_pipeline.R")
    py_src = _read(ROOT / "dssat_main_pipeline.py")

    for src in (r_src, py_src):
        assert "LATITUDE" in src
        assert "LONGITUDE" in src
        assert "ELEV" in src
        assert "fixed-column parser" in src or "fixed = TRUE" in src

    assert "Do NOT patch the LATITUDE/LONGITUDE/ELEV" not in py_src


def test_cropland_mask_config_and_outputs_are_aligned():
    r_src = _read(ROOT / "dssat_main_pipeline.R")
    py_src = _read(ROOT / "dssat_main_pipeline.py")
    cfg = _read(ROOT / "config.yml")

    markers = (
        "use_cropland_mask",
        "cropland_raster_file",
        "cropland_classes",
        "cropland_min_fraction",
        "cropland_strict",
        "crop_frac",
        "crop_pct",
        "crop_ha",
        "cell_ha",
        "cropland_ha",
        "gridcell_area_ha",
        "final_grain_production_kg",
        "top_weight_production_kg",
    )

    for marker in markers:
        assert marker in r_src
        assert marker in py_src

    for marker in markers[:5]:
        assert marker in cfg


def _config_keys_used(path: Path) -> set[str]:
    return set(re.findall(r"cfg_get\([\"']([^\"']+)", _read(path)))


def test_every_runtime_setting_is_declared_in_central_yaml():
    """Prevent new user settings from being hidden only in pipeline code."""
    cfg = yaml.safe_load((ROOT / "config.yml").read_text(encoding="utf-8"))
    declared = set(cfg)
    used = (
        _config_keys_used(ROOT / "dssat_main_pipeline.py")
        | _config_keys_used(ROOT / "dssat_main_pipeline.R")
    )
    # This removed key is read only to emit a migration error when supplied.
    used.discard("treatments")
    assert used <= declared, f"settings missing from config.yml: {sorted(used - declared)}"
    for source in ("dssat_main_pipeline.py", "dssat_main_pipeline.R"):
        text = _read(ROOT / source)
        assert re.search(r"ZIP_FOR_HPC\s*(?:=|<-)\s*.*cfg_get\([\"']zip_for_hpc", text)


def test_conda_lock_covers_supported_desktop_platforms():
    """Keep fresh Windows, Linux, and Apple Silicon installs reproducible."""
    lock = yaml.safe_load((ROOT / "conda_lock.yml").read_text(encoding="utf-8"))
    platforms = set(lock.get("metadata", {}).get("platforms", []))
    assert {"win-64", "linux-64", "osx-arm64"} <= platforms


def test_shared_package_pins_are_consistent_across_install_paths():
    """Keep CI, conda, renv, and documented fresh installs on one revision."""
    expected = {
        "dssatutils": "e9c859fa1d915623df23e2eb13084cb085dbfe3e",
        "dssatengine": "9f5bbde0def31dd74c5f881bf6b3be30f787c6a0",
    }
    e2e = yaml.safe_load(_read(ROOT / ".github" / "workflows" / "e2e.yml"))
    smoke = yaml.safe_load(_read(ROOT / ".github" / "workflows" / "smoke.yml"))
    renv = yaml.safe_load(_read(ROOT / "renv.lock"))
    for package, revision in expected.items():
        env_key = f"{package.upper()}_REF"
        assert e2e["env"][env_key] == revision
        assert smoke["env"][env_key] == revision
        assert renv["Packages"][package]["RemoteSha"] == revision
        for path in ("environment.yml", "conda_lock.yml", "setup_renv.R"):
            assert revision in _read(ROOT / path), f"stale {package} pin in {path}"
    assert expected["dssatutils"] in _read(ROOT / "README.md")


def test_resume_provenance_hashes_resolved_model_inputs_in_both_languages():
    """Resume keys must change with generated inputs and the DSSAT runtime."""
    required = (
        "weather_wth", "soil_sol", "soil_mapping", "dssat_executable",
        "dssatpro", "support_files",
    )
    python_source = _read(ROOT / "dssat_main_pipeline.py")
    r_source = _read(ROOT / "dssat_main_pipeline.R")
    for token in required:
        assert token in python_source
        assert token in r_source


def test_python_override_yaml_is_merged_over_central_defaults(tmp_path):
    override = tmp_path / "study.yml"
    override.write_text("project_name: override_study\nweather_start_year: 2001\n", encoding="utf-8")
    env = os.environ.copy()
    env["DSSAT_CONFIG_FILE"] = str(override)
    code = (
        "import config_loader; "
        "print(config_loader.cfg_get('project_name', 'missing')); "
        "print(config_loader.cfg_get('weather_source', 'missing'))"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], cwd=ROOT, env=env,
        check=True, capture_output=True, text=True,
    )
    assert [line.strip() for line in result.stdout.splitlines()[-2:]] == [
        "override_study", "NASA_POWER"
    ]


def test_r_override_yaml_is_merged_over_central_defaults(tmp_path):
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is not installed")
    override = tmp_path / "study.yml"
    override.write_text("project_name: override_study\nweather_start_year: 2001\n", encoding="utf-8")
    env = os.environ.copy()
    env["DSSAT_CONFIG_FILE"] = str(override)
    code = (
        "source('config_loader.R'); "
        "cat(cfg_get('project_name','missing'),'\\n'); "
        "cat(cfg_get('weather_source','missing'),'\\n')"
    )
    result = subprocess.run(
        [rscript, "--vanilla", "-e", code], cwd=ROOT, env=env,
        check=False, capture_output=True, text=True,
    )
    assert result.returncode == 0, (
        "R config loader failed.\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert [line.strip() for line in result.stdout.splitlines()[-2:]] == [
        "override_study", "NASA_POWER"
    ]
