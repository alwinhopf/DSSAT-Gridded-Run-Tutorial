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
    assert env['GIT_CONFIG_KEY_1'] == 'http.https://github.com/alwinhopf/cropmodel.git.extraheader'
    assert base64.b64decode(env['GIT_CONFIG_VALUE_1'].split()[-1]).decode() == 'x-access-token:fixture-secret'
    assert os.environ['GIT_CONFIG_COUNT'] == '1'


def test_missing_private_credentials_stop_before_install(monkeypatch):
    module = installer()
    monkeypatch.delenv('CROPMODEL_READ_TOKEN', raising=False)
    with pytest.raises(SystemExit, match='CROPMODEL_READ_TOKEN Actions secret'):
        module.main()


@pytest.mark.parametrize('token', ['', '   '])
def test_preflight_rejects_missing_or_blank_token(monkeypatch, token):
    module = installer()
    monkeypatch.setenv('CROPMODEL_READ_TOKEN', token)
    monkeypatch.setattr(module.subprocess, 'run', lambda *a, **k: pytest.fail('must not install'))
    with pytest.raises(SystemExit, match='GITHUB_TOKEN cannot read another private repository'):
        module.main(check_only=True)


def test_preflight_does_not_install_or_change_git_environment(monkeypatch):
    module = installer()
    monkeypatch.setenv('CROPMODEL_READ_TOKEN', 'fixture-secret')
    monkeypatch.setenv('GIT_CONFIG_COUNT', '7')
    monkeypatch.setattr(module.subprocess, 'run', lambda *a, **k: pytest.fail('must not install'))
    before = os.environ.copy()
    module.main(check_only=True)
    assert os.environ == before

def test_git_header_is_not_sent_to_other_repositories(monkeypatch, tmp_path):
    import shutil
    import subprocess
    module = installer()
    git = shutil.which('git')
    if not git:
        pytest.skip('Git unavailable')
    monkeypatch.setenv('CROPMODEL_READ_TOKEN', 'fixture-secret')
    monkeypatch.setenv('GIT_CONFIG_COUNT', '0')
    calls = []
    with monkeypatch.context() as patch:
        patch.setattr(module.subprocess, 'run', lambda args, **kwargs: calls.append(kwargs))
        module.main()
    env = calls[0]['env']
    # Actions checkout configures its own repository token. Test our header in
    # an empty directory, independent of checkout/user/system credentials.
    env['GIT_CONFIG_GLOBAL'] = os.devnull
    env['GIT_CONFIG_NOSYSTEM'] = '1'
    for url, expected in [
        ('https://github.com/alwinhopf/cropmodel.git', 0),
        ('https://github.com/alwinhopf/unrelated.git', 1),
        ('https://example.org/alwinhopf/cropmodel.git', 1),
    ]:
        result = subprocess.run([git, 'config', '--get-urlmatch', 'http.extraheader', url],
                                cwd=tmp_path, env=env, capture_output=True, text=True)
        assert result.returncode == expected
        if expected == 0:
            assert result.stdout.strip() == env['GIT_CONFIG_VALUE_0']
        else:
            assert not result.stdout
