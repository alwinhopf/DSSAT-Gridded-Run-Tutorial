#!/usr/bin/env bash
# File: scripts/pre-push.sh
# -----------------------------------------------------------------------------
# Local pre-push validation script. Runs the exact fast unit test suite executed
# by the primary PR CI lane (unit-tests in .github/workflows/smoke.yml).
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

echo "======================================================================"
echo "1. Running fast unit tests (pytest)"
echo "======================================================================"
pytest tests/test_smoke.py tests/test_agera5_config.py tests/test_failed_run_archive.py -v

echo "======================================================================"
echo "2. Running central config validation (Rscript)"
echo "======================================================================"

# Discover Rscript (respects RSCRIPT or PATH, checks common locations)
RSCRIPT_BIN="${RSCRIPT:-$(which Rscript 2>/dev/null || true)}"
if [ -z "${RSCRIPT_BIN}" ]; then
  for candidate in \
    "/usr/local/bin/Rscript" \
    "/opt/homebrew/bin/Rscript" \
    "/usr/bin/Rscript" \
    "/Library/Frameworks/R.framework/Resources/bin/Rscript" \
    "C:/Program Files/R/R-*/bin/Rscript.exe"; do
    if [ -x "${candidate}" ]; then
      RSCRIPT_BIN="${candidate}"
      break
    fi
  done
fi

if [ -n "${RSCRIPT_BIN}" ] && [ -x "${RSCRIPT_BIN}" ]; then
  "${RSCRIPT_BIN}" --vanilla -e 'source("config_loader.R"); stopifnot(cfg_get("weather_source", "") != "")'
  echo "R config check passed."
else
  echo "Warning: Rscript not found on PATH or standard locations; skipping R config check."
fi

echo "======================================================================"
echo "All local pre-push checks passed successfully!"
echo "======================================================================"
