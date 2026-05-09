"""Auto-convert legacy GWinput to GWinput.toml when needed.

Designed to be called once at the very start of any driver script
(gwsc, genMLWFx, job_mlo, ...), before any ecalj binary is launched.

Behavior:
  - GWinput.toml present  -> nothing to do
  - GWinput present       -> run gwinput2toml.py to generate GWinput.toml
  - neither               -> error out with clear message

This is the script-side counterpart of m_GWinput::gwinput_init,
which now refuses to start without GWinput.toml.
"""
from __future__ import annotations
import os
import sys
import subprocess
from pathlib import Path


def ensure_gwinput_toml(exec_dir=None, cwd=None):
    """Ensure GWinput.toml exists in cwd. Returns True if created or already present.

    exec_dir : Path | None
        Directory containing gwinput2toml.py. Defaults to dirname of this file.
    cwd : Path | None
        Working directory to operate in. Defaults to os.getcwd().
    """
    if exec_dir is None:
        exec_dir = Path(__file__).resolve().parent
    else:
        exec_dir = Path(exec_dir)
    if cwd is None:
        cwd = Path.cwd()
    else:
        cwd = Path(cwd)

    toml = cwd / 'GWinput.toml'
    legacy = cwd / 'GWinput'

    if toml.is_file():
        return True
    if legacy.is_file():
        converter = exec_dir / 'gwinput2toml.py'
        if not converter.is_file():
            print(f'ERROR: GWinput.toml is missing and converter not found at {converter}', file=sys.stderr)
            return False
        print(f'Auto-generating GWinput.toml from GWinput (via {converter.name})')
        result = subprocess.run([sys.executable, str(converter), 'GWinput'], cwd=str(cwd))
        if result.returncode != 0:
            print('ERROR: gwinput2toml.py failed', file=sys.stderr)
            return False
        return toml.is_file()
    print(f'ERROR: neither GWinput.toml nor GWinput found in {cwd}', file=sys.stderr)
    return False


if __name__ == '__main__':
    ok = ensure_gwinput_toml()
    sys.exit(0 if ok else 1)
