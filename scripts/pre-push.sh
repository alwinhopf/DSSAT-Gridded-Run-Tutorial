#!/usr/bin/env bash
# Exact fast CI lane, including the central R configuration and syntax checks.
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="$(command -v python || command -v python3)"
"$PYTHON_BIN" -m pytest tests/test_smoke.py tests/test_agera5_config.py tests/test_failed_run_archive.py tests/test_provider_cache_and_fields.py -v --tb=short --junit-xml=test-results.xml 2>&1 | tee fast-tests.log
RSCRIPT_BIN="$("$PYTHON_BIN" -c 'from dssatutils import find_rscript; path=find_rscript(); assert path, "Rscript is required for the fast CI lane"; print(path)')"
"$RSCRIPT_BIN" --vanilla -e 'source("config_loader.R"); stopifnot(cfg_get("weather_source", "") != ""); parse(file="dssat_main_pipeline.R"); parse(file="tests/test_e2e.R"); parse(file="tests/test_e2e_comprehensive.R"); parse(file="tests/test_dssat_model_e2e.R")' > r-syntax.log
"$PYTHON_BIN" -m compileall -q dssat_main_pipeline.py config_loader.py python_scripts tests
