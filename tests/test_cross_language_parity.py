from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def test_ssurgo_failure_diagnostics_are_present_in_r_and_python():
    r_src = _read(WORKSPACE / "dssatutils" / "R" / "soil_ssurgo.R")
    py_src = _read(WORKSPACE / "dssatutils" / "python" / "dssatutils" / "soil_ssurgo.py")

    for marker in ("no-coverage", "no-soil", "no-layers", "_download_failures.csv"):
        assert marker in r_src
        assert marker in py_src


def test_engine_missing_summary_failure_logging_is_present_in_r_and_python():
    r_src = _read(WORKSPACE / "dssatengine" / "R" / "engine.R")
    py_src = _read(WORKSPACE / "dssatengine" / "python" / "dssatengine" / "engine.py")

    for src in (r_src, py_src):
        assert "_run_error.log" in src
        # Robust to wording ("DSSAT produced no ..." / "DSSAT completed but
        # produced no ...") — assert the functional parity, not the exact phrase.
        assert "produced no 'summary.csv'" in src
        assert "FMOPT" in src
        assert "summary.csv" in src


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
