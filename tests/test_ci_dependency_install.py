"""Private dependency credentials must stay in the install subprocess."""
import base64
import importlib.util
import os
from pathlib import Path

import pytest


def installer():
    path = Path(__file__).resolve().parents[1] / 'scripts/install_ci_dependencies.py'
    spec = importlib.util.spec_from_file_location('ci_install', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_private_credentials_are_scoped_to_pip(monkeypatch):
    module = installer()
    monkeypatch.setenv('CROPMODEL_READ_TOKEN', 'fixture-secret')
    monkeypatch.setenv('GIT_CONFIG_COUNT', '1')
    monkeypatch.setenv('GIT_CONFIG_KEY_0', 'core.autocrlf')
    monkeypatch.setenv('GIT_CONFIG_VALUE_0', 'false')
    calls = []
    monkeypatch.setattr(module.subprocess, 'run', lambda args, **kwargs: calls.append((args, kwargs)))
    module.main()
    args, kwargs = calls[0]
    assert 'fixture-secret' not in ' '.join(args)
    assert kwargs['check'] is True
    env = kwargs['env']
    assert env['GIT_CONFIG_COUNT'] == '2'
    assert env['GIT_CONFIG_VALUE_0'] == 'false'
    assert env['GIT_CONFIG_KEY_1'] == 'http.https://github.com/.extraheader'
    assert base64.b64decode(env['GIT_CONFIG_VALUE_1'].split()[-1]).decode() == 'x-access-token:fixture-secret'
    assert os.environ['GIT_CONFIG_COUNT'] == '1'


def test_missing_private_credentials_stop_before_install(monkeypatch):
    module = installer()
    monkeypatch.delenv('CROPMODEL_READ_TOKEN', raising=False)
    with pytest.raises(SystemExit, match='CROPMODEL_READ_TOKEN Actions secret'):
        module.main()
