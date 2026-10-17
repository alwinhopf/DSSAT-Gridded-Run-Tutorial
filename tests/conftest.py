"""Pytest configuration and fixtures for DSSAT_Gridded_Run_Tutorial tests."""
from pathlib import Path
import sys

import pytest

# Ensure tests/helpers and repository root are on sys.path
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tests.helpers.packages import configure_workspace_packages
configure_workspace_packages()

from tests.helpers.discovery import find_rscript


@pytest.fixture
def rscript() -> str:
    """Fixture returning resolved absolute path to Rscript, or skips cleanly."""
    path = find_rscript()
    if not path:
        pytest.skip("Rscript executable is not installed or discoverable")
    return path


def pytest_report_header(config):
    import importlib
    return [f"{name}: {importlib.import_module(name).__file__}" for name in ('cropmodel_data','cropmodel_runtime')]
