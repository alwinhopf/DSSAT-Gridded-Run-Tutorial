from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_pipelines_forward_all_agera5_v2_controls():
    for filename in ("dssat_main_pipeline.R", "dssat_main_pipeline.py"):
        source = (ROOT / filename).read_text(encoding="utf-8")
        for config_key, runtime_name in (
            ("agera5_backend", "AGERA5_BACKEND"),
            ("agera5_data_format", "AGERA5_DATA_FORMAT"),
            ("agera5_timeseries_chunk_degrees", "AGERA5_TIMESERIES_CHUNK_DEGREES"),
            ("agera5_max_concurrent_requests", "AGERA5_MAX_CONCURRENT_REQUESTS"),
        ):
            assert config_key in source
            assert runtime_name in source


def test_default_agera5_config_values():
    config = (ROOT / "config.yml").read_text(encoding="utf-8")
    assert 'agera5_backend: "timeseries"' in config
    assert 'agera5_data_format: "csv"' in config
    assert "agera5_timeseries_chunk_degrees: 0.1" in config
    assert "agera5_max_concurrent_requests: 4" in config
