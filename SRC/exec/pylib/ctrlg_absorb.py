"""Fold retired side files into ctrlg.<sname>.toml.

The Fortran reads ctrlg.<sname>.toml and nothing else. Two files that older
working directories still carry are converted here, never by the binaries
(they abort and name this module's CLI, ctrlg_absorb.py):

  esm_input.dat      (positional, 5 lines)  -> [esm] section
  PB.<sname>.toml    ([product_basis] nlx / valence / core, until 2026-09)
                                            -> appended to [product_basis]

Each original is kept as <name>.bk with a header saying where it went.
Legacy2toml.py calls absorb_esm() after converting ctrl.<sname>.
"""
from __future__ import annotations
import datetime
import re
import sys
import tomllib
from pathlib import Path

ESM_BOUNDARY = {
    0: 'off', 1: 'vac/slab/vac', 2: 'metal/slab/metal', 3: 'vac/slab/metal',
    4: 'metal/slab/vac', 5: 'vac/slab/vac:field', 6: 'metal/slab/metal:v-e',
    7: 'metal/slab/metal:e-v', 10: 'periodic:esm', 11: 'periodic:esm',
}


def _say(msg):
    print(f'=== ctrlg_absorb: {msg}', flush=True)


def _archive(src: Path, ctrlg: Path, note: list[str]):
    when = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    head = [f'# {src.name} is retired; this is the archived original.',
            f'# Moved {when}.'] + note + ['#', '# --- original contents below ---']
    Path(str(src) + '.bk').write_text('\n'.join(head) + '\n' + src.read_text())
    src.unlink()


def absorb_esm(ctrlg: Path) -> bool:
    """esm_input.dat -> [esm] in ctrlg. Returns True when something was done."""
    esm = Path('esm_input.dat')
    if not esm.exists():
        return False
    if re.search(r'^\[esm\]', ctrlg.read_text(), re.M):
        _archive(esm, ctrlg, [
            f'# {ctrlg.name} already had an [esm] section, so these',
            f'# settings were NOT copied anywhere -- the TOML one is what runs.',
            f'# Compare the two if you expected this file to be in effect.'])
        _say(f'{ctrlg.name} already has [esm]; {esm} -> {esm}.bk (unused)')
        return True
    nums = []
    for line in esm.read_text().splitlines():
        line = line.split('#')[0].strip()
        if line:
            nums.extend(line.split())
    if len(nums) < 9:
        sys.exit(f'ctrlg_absorb: {esm} has {len(nums)} values, expected 9')
    jesm = int(float(nums[0]))
    jtresm = int(float(nums[1]))
    tresm, z1, z2, vp, vm, ep, em = (float(x) for x in nums[2:9])
    # the commented-out template ctrlgenToml.py writes is superseded by the real section
    from pylib.toml_comments import ESM_SAMPLE, SECTION_HEADER
    text = ctrlg.read_text()
    sample = '\n'.join(ESM_SAMPLE) + '\n'
    if sample in text:
        ctrlg.write_text(text.replace(sample, ''))
    with open(ctrlg, 'a') as f:
        f.write('\n' + '\n'.join(SECTION_HEADER['esm']) + '\n')
        f.write(f'''# Converted from esm_input.dat.
[esm]
boundary  = "{ESM_BOUNDARY.get(jesm, 'off')}"
            # "off" / "vac/slab/vac" / "metal/slab/metal" / "vac/slab/metal" /
            # "metal/slab/vac" / "vac/slab/vac:field" / "metal/slab/metal:v-e" /
            # "metal/slab/metal:e-v" / "periodic:esm"
origin    = {tresm!r}  # (a.u.) z-translation of the density; code applies -origin
shiftmode = {jtresm}  # 0: origin is absolute, 1: in units of the cell length
zb        = [{z1!r}, {z2!r}]  # (a.u.) boundaries z1esm, z2esm
potential = [{vp!r}, {vm!r}]  # (Ry) vesmp, vesmm on the +z / -z sides
field     = [{ep!r}, {em!r}]  # (Ry/a.u.) eesmp, eesmm on the +z / -z sides
''')
    _archive(esm, ctrlg, [
        f'# These settings were converted and appended to {ctrlg.name}',
        f'# as an [esm] section. Edit them there from now on; this file is',
        f'# no longer read.'])
    _say(f'{esm} -> [esm] in {ctrlg.name}; original kept as {esm}.bk')
    return True


def absorb_pb(ctrlg: Path, sname: str) -> bool:
    """PB.<sname>.toml -> nlx / valence / core appended to [product_basis] of ctrlg."""
    pb = Path(f'PB.{sname}.toml')
    if not pb.exists():
        return False
    text = ctrlg.read_text()
    have = tomllib.loads(text).get('product_basis', {})
    if 'nlx' in have:
        _archive(pb, ctrlg, [
            f'# {ctrlg.name} already carries nlx/valence/core in [product_basis],',
            f'# so nothing was copied; the TOML ones are what runs.'])
        _say(f'{ctrlg.name} already has the tables; {pb} -> {pb}.bk (unused)')
        return True
    if not re.search(r'^\[product_basis\]', text, re.M):
        sys.exit(f'ctrlg_absorb: {ctrlg.name} has no [product_basis] section; '
                 f'run  ctrlgenToml.py {sname} --addgw  (or gwinit) first')
    if not text.rstrip('\n').split('\n[')[-1].startswith('product_basis]'):
        # ctrlg written before 2026-09 ([gw] [product_basis] [mlo] ... order): reorder first
        from pylib.toml_tidy import tidy_gw_sections
        new = tidy_gw_sections(text)
        if tomllib.loads(new) != tomllib.loads(text):
            sys.exit(f'ctrlg_absorb: internal error, tidy changed the TOML content of {ctrlg.name}')
        ctrlg.write_text(new)
        text = new
        _say(f'{ctrlg.name}: sections reordered so that [product_basis] comes last')
        if not text.rstrip('\n').split('\n[')[-1].startswith('product_basis]'):
            sys.exit(f'ctrlg_absorb: could not move [product_basis] to the end of {ctrlg.name}')
    pbtext = pb.read_text()
    blocks = []
    for key in ('nlx', 'valence', 'core'):
        m = re.search(rf'^{key}\s*=\s*\[\n.*?^\]', pbtext, re.S | re.M)
        if m:
            blocks.append(m.group(0).rstrip())
    if not blocks:
        _archive(pb, ctrlg, [f'# {pb.name} carried no nlx/valence/core tables; nothing to copy.'])
        _say(f'{pb} had no tables; -> {pb}.bk')
        return True
    ctrlg.write_text(text.rstrip('\n') + '\n' + '\n'.join(b + '\n' for b in blocks))
    _archive(pb, ctrlg, [
        f'# nlx / valence / core were appended to [product_basis] of {ctrlg.name}.',
        f'# Edit them there from now on; this file is no longer read.'])
    _say(f'{pb} -> [product_basis] of {ctrlg.name}; original kept as {pb}.bk')
    return True


def finish(ctrlg: Path):
    """Annotate + tidy, asserting the TOML content is unchanged."""
    from pylib.toml_comments import apply_toml_annotations
    from pylib.toml_tidy import tidy_gw_sections
    raw = ctrlg.read_text()
    new = tidy_gw_sections(apply_toml_annotations(raw))
    if tomllib.loads(raw) != tomllib.loads(new):
        sys.exit(f'ctrlg_absorb: internal error, tidy changed the TOML content of {ctrlg.name}')
    ctrlg.write_text(new)
