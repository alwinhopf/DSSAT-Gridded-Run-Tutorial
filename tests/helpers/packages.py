"""Source discovery without silently replacing the installed test packages."""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def monorepo_root():
    explicit = os.environ.get('CROPMODEL_ROOT')
    if explicit:
        path = Path(explicit).expanduser().resolve()
        if not (path / 'packages').is_dir():
            raise RuntimeError(f'CROPMODEL_ROOT has no packages directory: {path}')
        return path
    ci = ROOT / '.ci-deps' / 'cropmodel'
    return ci if ci.is_dir() else ROOT.parent / 'cropmodel'


def package_source(name):
    return monorepo_root() / 'packages' / name


def configure_workspace_packages():
    if os.environ.get('DSSAT_USE_WORKSPACE_SIBLINGS') != '1':
        return
    for name in ('cropmodel_data', 'cropmodel_runtime'):
        source = package_source(name) / 'python'
        if not (source / name / '__init__.py').is_file():
            raise RuntimeError(f'Requested local package is missing: {source}')
        if str(source) not in sys.path:
            sys.path.insert(0, str(source))
