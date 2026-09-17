#!/usr/bin/env python3
"""
ctrlg_absorb.py <sname> — fold retired side files into ctrlg.<sname>.toml.

The Fortran binaries read ctrlg.<sname>.toml and nothing else; when they
find one of these files they abort and name this script:

    esm_input.dat     -> [esm] section            (ESM for slabs)
    PB.<sname>.toml   -> [product_basis] nlx / valence / core
                         (the per-atom product-basis tables, a separate
                          file until 2026-09)

Each original is kept as <name>.bk with a header saying where it went.
ctrlg.<sname>.toml is then re-annotated and tidied (content unchanged).
Nothing happens when neither file is present.
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.ctrlg_absorb import absorb_esm, absorb_pb, finish


def main():
    if len(sys.argv) != 2 or sys.argv[1].startswith('-'):
        sys.exit(__doc__)
    sname = sys.argv[1]
    ctrlg = Path(f'ctrlg.{sname}.toml')
    if not ctrlg.exists():
        sys.exit(f'ctrlg_absorb: {ctrlg} not found in cwd '
                 f'(run Legacy2toml.py {sname} for a legacy ctrl.{sname} + GWinput directory)')
    done = absorb_esm(ctrlg)
    done = absorb_pb(ctrlg, sname) or done
    if done:
        finish(ctrlg)
    else:
        print(f'ctrlg_absorb: nothing to do (no esm_input.dat, no PB.{sname}.toml)')


if __name__ == '__main__':
    main()
