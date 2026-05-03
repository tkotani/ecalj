#!/usr/bin/env python3
"""
Legacy2toml.py — one-shot converter from legacy ecalj inputs to TOML.

Reads:
    ctrl.<sname>      (legacy ctrl text)
    GWinput           (legacy GWinput; optional, only needed for GW jobs)

Writes:
    ctrlG.<sname>.toml    [ctrl sections + [gw], [product_basis] scalars, [blocks]]
    PB.toml               [per-atom product-basis tables: nlx, valence, core]

Usage:
    Legacy2toml.py <sname> [-vfoo=bar ...]

Internally it invokes the existing ctrl2ctrltoml.py and gwinput2toml.py,
then post-processes their output:
  1. ctrl.<sname>.toml -> ctrlG.<sname>.toml (rename)
  2. GWinput.toml -> append [gw], [product_basis] (slim), [blocks] sections
     to ctrlG.<sname>.toml; move per-atom arrays to PB.toml.

Verbose: each step is echoed with a banner so the user knows exactly
which legacy file was consumed and which TOML files were written.
"""
import os, sys, re, shutil, subprocess
from pathlib import Path

def banner(msg):
    print(f'=== Legacy2toml.py: {msg}', flush=True)

def main():
    if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
        sys.exit(__doc__.strip())
    sname = sys.argv[1]
    extra_args = sys.argv[2:]

    here = Path(os.path.dirname(os.path.realpath(__file__)))

    ctrl_legacy = Path(f'ctrl.{sname}')
    if not ctrl_legacy.exists():
        sys.exit(f'Legacy2toml.py: {ctrl_legacy} not found')

    # ------------------------------------------------------------------
    # 1. ctrl.<sname>  ->  ctrlG.<sname>.toml  (via ctrl2ctrltoml.py)
    # ------------------------------------------------------------------
    banner(f'ctrl.{sname} -> ctrlG.{sname}.toml')
    out_path = Path(f'ctrlG.{sname}.toml')
    with open(out_path, 'wb') as fout:
        rc = subprocess.run(
            [str(here / 'ctrl2ctrltoml.py'), *extra_args],
            stdin=open(ctrl_legacy, 'rb'), stdout=fout, stderr=sys.stderr,
        ).returncode
    if rc != 0:
        sys.exit(f'Legacy2toml.py: ctrl2ctrltoml.py returned {rc}')

    # ------------------------------------------------------------------
    # 2. GWinput  ->  ctrlG.<sname>.toml [gw, product_basis (slim), blocks]
    #              +  PB.toml   [per-atom nlx / valence / core]
    # (skipped if no GWinput in cwd)
    # ------------------------------------------------------------------
    gwinput = Path('GWinput')
    if not gwinput.exists():
        banner(f'no GWinput in cwd; ctrlG.{sname}.toml emitted without GW sections')
        return
    banner(f'GWinput -> append [gw]/[product_basis]/[blocks] to ctrlG.{sname}.toml + PB.toml')
    # Run gwinput2toml.py to get a GWinput.toml in the legacy "all in one" form
    tmp_gwinput_toml = Path('GWinput.toml.tmp.l2t')
    rc = subprocess.run(
        [str(here / 'gwinput2toml.py'), 'GWinput', '-o', str(tmp_gwinput_toml)],
        stdout=subprocess.DEVNULL, stderr=sys.stderr,
    ).returncode
    if rc != 0:
        sys.exit(f'Legacy2toml.py: gwinput2toml.py returned {rc}')

    text = tmp_gwinput_toml.read_text()
    tmp_gwinput_toml.unlink()

    # Split [product_basis]:
    #   - keep tolerance/lcutmx in ctrlG (rename to pb_tolerance/pb_lcutmx)
    #   - move nlx/valence/core to PB.toml
    pb_tol_match = re.search(r'^tolerance\s*=\s*(\[[^\]]*\])', text, re.MULTILINE)
    pb_lcm_match = re.search(r'^lcutmx\s*=\s*(\[[^\]]*\])',    text, re.MULTILINE)
    pb_tol = pb_tol_match.group(1) if pb_tol_match else '[1e-3]'
    pb_lcm = pb_lcm_match.group(1) if pb_lcm_match else '[]'

    def grab_block(name):
        m = re.search(rf'^{name}\s*=\s*\[\n(.*?)^\]', text, re.DOTALL | re.MULTILINE)
        return m.group(0) if m else None
    nlx_block     = grab_block('nlx')
    valence_block = grab_block('valence')
    core_block    = grab_block('core')

    # Build slim [gw] + [product_basis] + [blocks] for ctrlG
    # Strategy: take everything after the first [gw] in tmp output, drop
    # the per-atom arrays, drop the legacy 'tolerance'/'lcutmx' assignment
    # (we'll write pb_tolerance/pb_lcutmx instead), and keep [blocks] as-is.
    gw_section = re.search(r'^\[gw\].*?(?=^\[)', text, re.DOTALL | re.MULTILINE)
    blocks_section = re.search(r'^\[blocks\].*', text, re.DOTALL | re.MULTILINE)

    gw_text = gw_section.group(0).rstrip() if gw_section else '[gw]\n'
    pb_text = (
        '[product_basis]\n'
        f'pb_tolerance = {pb_tol}\n'
        f'pb_lcutmx    = {pb_lcm}\n'
        f'\n# Per-atom product-basis tables (nlx / valence / core) live in PB.toml\n'
    )
    blocks_text = blocks_section.group(0).rstrip() + '\n' if blocks_section else ''
    # Strip vestigial PRODUCT_BASIS raw-text block from [blocks]
    blocks_text = re.sub(
        r'PRODUCT_BASIS\s*=\s*"""\n.*?\n"""\s*\n?',
        '', blocks_text, flags=re.DOTALL,
    )

    with open(out_path, 'a') as f:
        f.write('\n')
        f.write(gw_text)
        f.write('\n')
        f.write(pb_text)
        f.write('\n')
        f.write(blocks_text)

    # Build PB.toml
    with open('PB.toml', 'w') as f:
        f.write('# PB.toml -- per-atom product-basis tables (auto-generated by Legacy2toml.py)\n\n')
        f.write('[product_basis]\n\n')
        for blk in (nlx_block, valence_block, core_block):
            if blk:
                f.write(blk.rstrip() + '\n\n')

    banner(f'wrote ctrlG.{sname}.toml + PB.toml')

if __name__ == '__main__':
    main()
