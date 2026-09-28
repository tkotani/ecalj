#!/usr/bin/env python3
"""Bring ctrlg.<sname>.toml files up to the current rules (2026-09-28); the calculation stays the same.

  ctrlg_update.py [--no-backup] <ctrlg file> [...]      (in place; the original is kept as <file>.bak_update)

  1. [struc] nbas / nspec are removed (lmf stops on them): the numbers of sites and species are those of the
     [[site]] / [[spec]] tables. If nbas was smaller than the tables (the old way to use only the first sites), the
     [[site]] blocks beyond it are commented out. To leave a site out from now on, comment out its [[site]] block.
  2. [gw] t_tetrakbt is required (0 = no smearing of chi0, > 0 = finite-T tetrahedron, < 0 = the Gaussian smearing
     of Im chi0 as wide as Fermi-Dirac at |T|). SmearX0 = s becomes t_tetrakbt = -T with the same width
     (0.0057 Ha -> -992.4 K); a missing t_tetrakbt becomes 0, the former default.
The rules themselves are in pylib/ctrlg_rules.py.
"""
import shutil
import sys
from pylib.ctrlg_rules import update

if __name__ == '__main__':
    args = sys.argv[1:]
    backup = '--no-backup' not in args
    files = [a for a in args if a != '--no-backup']
    if not files:
        sys.exit(__doc__)
    for f in files:
        text = open(f).read()
        new, msg = update(text, f)
        if new != text:
            if backup:
                shutil.copy(f, f + '.bak_update')
            open(f, 'w').write(new)
        print(f'{f}: {msg}')
