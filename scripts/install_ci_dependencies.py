"""Install pinned dependencies with scoped access to the private monorepo."""
import base64
import argparse
import os
from pathlib import Path
import subprocess
import sys


def main(check_only=False):
    token = os.environ.get('CROPMODEL_READ_TOKEN', '').strip()
    if not token:
        raise SystemExit(
            'Set the CROPMODEL_READ_TOKEN Actions secret to a token with '
            'Contents: read access to alwinhopf/cropmodel. The automatic '
            'GITHUB_TOKEN cannot read another private repository. Fork pull '
            'requests require a trusted run after review.')
    if check_only:
        return
    env = os.environ.copy()
    # Give only pip's Git children the header; never persist credentials or put
    # them in URLs, command arguments, logs, or the global Git configuration.
    count = int(env.get('GIT_CONFIG_COUNT', '0'))
    env[f'GIT_CONFIG_KEY_{count}'] = 'http.https://github.com/alwinhopf/cropmodel.git.extraheader'
    auth = base64.b64encode(f'x-access-token:{token}'.encode()).decode()
    env[f'GIT_CONFIG_VALUE_{count}'] = f'AUTHORIZATION: basic {auth}'
    env['GIT_CONFIG_COUNT'] = str(count + 1)
    env['GIT_TERMINAL_PROMPT'] = '0'
    requirements = Path(__file__).resolve().parents[1] / 'requirements.txt'
    subprocess.run([sys.executable, '-m', 'pip', 'install', '-r', str(requirements)],
                   env=env, check=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only', action='store_true',
                        help='Check required credentials without installing packages.')
    main(parser.parse_args().check_only)
