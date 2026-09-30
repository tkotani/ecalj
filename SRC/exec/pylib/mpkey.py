"""The Materials Project API key (2026-10-01).

The key is personal: it must not be in the repository. Looked for in this order:
  1. the environment variable MP_API_KEY
  2. <ecalj>/MaterialProject.key  -- the top directory of the ecalj tree this file belongs to (ignored by git)
  3. ~/ecalj/MaterialProject.key  -- for a frozen copy of the bindir, which is not inside a tree
The file holds the key on one line; '#' starts a comment. <ecalj>/MaterialProject.key.example is the public template.
Bug fixed 2026-10-01: the key was read from ecalj_auto/config.ini, which is tracked, and auto_jobsubmit.py copied it into
every OUTPUT/<run>/config.ini; the key went into the public repository with them.
"""
import os
import sys
from pathlib import Path

KEYFILE = 'MaterialProject.key'
PLACEHOLDER = 'YourAPIkeyToMaterialProject'


def ecalj_top():
    """<ecalj>: this file is <ecalj>/SRC/exec/pylib/mpkey.py (bin/pylib is a symlink to it; resolve() follows it)."""
    return Path(__file__).resolve().parents[3]


def _read(p):
    try:
        for line in p.read_text().splitlines():
            line = line.split('#')[0].strip()
            if line:
                return line
    except OSError:
        pass
    return None


def mp_api_key(required=True):
    """The key, or None when there is none and required is False."""
    key = os.environ.get('MP_API_KEY', '').strip() or None
    if key is None:
        for p in (ecalj_top() / KEYFILE, Path.home() / 'ecalj' / KEYFILE):
            key = _read(p)
            if key is not None:
                break
    if key == PLACEHOLDER:
        key = None
    if key is None and required:
        sys.exit(f'Materials Project API key not found. Write it (one line) in {ecalj_top() / KEYFILE}\n'
                 f'  (cp {ecalj_top() / KEYFILE}.example {ecalj_top() / KEYFILE}; the file is ignored by git)\n'
                 f'  or set MP_API_KEY. The key is on your dashboard at https://materialsproject.org')
    return key
