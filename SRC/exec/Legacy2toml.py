#!/usr/bin/env python3
"""
Legacy2toml.py — one-shot migration tool: legacy ecalj input -> TOML.

==============================================================================
  IMPORTANT for users of ecalj prior to 2026-05
==============================================================================
  As of 2026-05, the Fortran binaries (lmf / lmfa / lmchk / gwsc / hsfp0 ...)
  read ONLY structured TOML:
      ctrlG.<sname>.toml   (merged ctrl + GW driver sections + PB cut-offs)
      PB.toml              (per-atom product basis tables, sname-free)
  Legacy text inputs (ctrl.<sname>, GWinput) are NO LONGER read by Fortran.

  To migrate an existing working directory:
      $ cd <your-old-dir>          # contains ctrl.<sname> and GWinput
      $ Legacy2toml.py <sname>     # produces ctrlG.<sname>.toml + PB.toml
      # then resume your normal workflow (lmf, gwsc, ...) unchanged.

  Run-time tunables (-v) have moved from %const to TOML-path syntax:
      OLD:  lmf si -vnk=8 -vmetal=3
      NEW:  lmf si -v[bz.nkabc]=[8,8,8] -v[bz.metal]=3
  The -v[<path>]=val form is processed in-memory by m_toml_override.f90;
  it never rewrites the .toml file on disk.

==============================================================================
  Inputs (read in cwd)
==============================================================================
    ctrl.<sname>      legacy ctrl text  (required)
    GWinput           legacy GWinput    (optional; only for GW workflows)

==============================================================================
  Outputs (written / overwritten in cwd; .bakup of any prior file is kept)
==============================================================================
    ctrlG.<sname>.toml    [io]/[struc]/[[site]]/[[spec]]/[bz]/[iter]/[ham]/...
                          plus [gw]/[product_basis] (scalars only)/[blocks]
    PB.toml               [product_basis] per-atom: nlx, valence, core

  Both files include human-readable inline comments (units / role of each
  key); see toml_comments.py to update the wording in one place.

==============================================================================
  Usage
==============================================================================
    Legacy2toml.py <sname>                    # one-shot conversion
    Legacy2toml.py <sname> -vfoo=bar -vbaz=2  # bake %const overrides
    Legacy2toml.py -h | --help                # show this help

  -v handling (only useful when generating ctrlG variants):
    Each "-vNAME=VAL" overrides %const NAME=... before {NAME} substitution
    in the legacy ctrl, baking the resulting numbers into ctrlG.<sname>.toml.

    A 3-level diagnostic warns about each -v BEFORE conversion:
      [WARN]  NAME not in %const         -- the override is a no-op
      [INFO]  NAME maps to a TOML path   -- prefer runtime -v[<path>]=val
                                           (no reconversion needed)
      [ERROR] NAME affects topology      -- variant required:
              save the result as ctrlG.<sname>.<tag>.toml and switch via
              "cp" before each lmf invocation.  See Samples/TestInstall/te
              for a worked example.

==============================================================================
  Pipeline (internal)
==============================================================================
  1. ctrl2ctrltoml.py  : ctrl.<sname>  -> ctrlG.<sname>.toml (sections from
                                          ctrl_schema.py, typed)
  2. gwinput2toml.py   : GWinput       -> intermediate GWinput.toml.tmp
  3. split + append    : [gw] + [product_basis] (pb_tolerance / pb_lcutmx)
                         + [blocks] are appended to ctrlG.<sname>.toml;
                         per-atom nlx / valence / core go to PB.toml.
  4. apply annotations : toml_comments.py inserts SECTION_HEADER blocks and
                         unit-bearing inline comments (idempotent).

  Each step is echoed with a "=== Legacy2toml.py: ..." banner so failures
  point you straight at the offending file.
"""
import os, sys, re, shutil, subprocess
from pathlib import Path
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from toml_comments import apply_toml_annotations

def banner(msg):
    print(f'=== Legacy2toml.py: {msg}', flush=True)


# ----- -v override analyzer (3-level diagnostic) ----------------------------
# Tokens that inherently change schema topology / structure of TOML;
# overriding them via -v[<path>]=val cannot be expressed in default TOML mode.
TOPOLOGY_TOKENS = {
    ('STRUC', 'NBAS'),
    ('STRUC', 'NSPEC'),
    ('STRUC', 'NL'),
    ('STRUC', 'NCLASS'),
    ('STRUC', 'NLFIT'),
    ('STRUC', 'SHEAR'),
    ('STRUC', 'TET'),
    ('STRUC', 'SLAT'),
    ('STRUC', 'FILE'),
}

def _load_schema():
    """Lazy import ctrl_schema (used to look up legacy_cattok -> TOML path)."""
    here = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, here)
    try:
        from ctrl_schema import SCHEMA
        return SCHEMA
    except Exception:
        return {}

def _v_pairs_from_args(args):
    """Yield (name, val) for each '-vNAME=VAL' style argv token."""
    for a in args:
        m = re.match(r'^-v([A-Za-z_][A-Za-z0-9_-]*)=(.*)$', a)
        if m:
            yield m.group(1), m.group(2)

def _const_defs(ctrl_text):
    """Map %const VAR -> source line of definition. Picks last-wins."""
    out = {}
    for line in ctrl_text.splitlines():
        if not line.strip().startswith('%') or 'const' not in line:
            continue
        body = line.split('const', 1)[1].split('#', 1)[0]
        body = re.sub(r'=\s+', '=', body)
        for tok in body.split():
            if '=' in tok:
                name = tok.split('=', 1)[0].strip()
                if name:
                    out[name] = line.strip()
    return out

def _find_var_refs(ctrl_text, name):
    """Find each ctrl line where {name} is referenced, returning a list of
    (cat, tok, line_no, line) tuples. cat is the most recent category line
    starting at column 0; tok is the most recent KEY= token before {name}."""
    refs = []
    cat = None
    needle = '{' + name + '}'
    for i, line in enumerate(ctrl_text.splitlines(), start=1):
        # category lines start at column 0 with an upper-case word
        m = re.match(r'^([A-Z][A-Z_]+)\b', line)
        if m:
            cat = m.group(1)
        if needle in line:
            # find which TOK= (= sub-token name) precedes the reference
            before = line.split(needle, 1)[0]
            tk_matches = re.findall(r'([A-Z][A-Z0-9_]*)\s*=', before)
            tok = tk_matches[-1] if tk_matches else None
            # special-case: line begins with a sole {name} (e.g. "{dyn}     MODE=...")
            if line.lstrip().startswith(needle):
                tok = '<SECTION_HEADER>'
            refs.append((cat, tok, i, line.rstrip()))
    return refs

def _toml_path_for(cat, tok, schema):
    """Look up TOML path for (CAT, TOK). Returns string like 'ham.gmax' or
    None if not in schema."""
    if cat is None or tok is None:
        return None
    spec = schema.get(f'{cat}_{tok}')
    if spec is None:
        return None
    section, key, _typ, _legacy_sub = spec
    if section is None:
        return key
    if section in ('site', 'spec'):
        return f'{section}.<idx>.{key}'   # caller picks index
    return f'{section}.{key}'

def analyze_v_overrides(ctrl_path, v_args):
    """Walk every -v argument and emit one of three diagnostic levels:

      1. WARN if the var is not defined in any %const     (no-op override)
      2. ERROR if any reference site is structural/topology (cannot be
         expressed as -v[<path>]=val in default TOML mode -- requires a
         pre-generated ctrlG.<sname>.<tag>.toml variant).
      3. INFO otherwise: suggest the equivalent -v[<toml-path>]=val that
         can be used in default TOML mode without going through Legacy2toml.

    Returns the count of ERRORs (caller decides whether to abort).
    """
    ctrl_text = open(ctrl_path).read()
    consts = _const_defs(ctrl_text)
    schema = _load_schema()
    n_err = 0
    pairs = list(_v_pairs_from_args(v_args))
    if not pairs:
        return 0
    print('--- Legacy2toml.py: analyzing -v overrides ---')
    for name, val in pairs:
        if name not in consts:
            print(f'  [WARN] -v{name}={val}: '
                  f"'{name}' is not defined in any %const of {ctrl_path}; "
                  f"override has no effect.")
            continue
        refs = _find_var_refs(ctrl_text, name)
        if not refs:
            print(f'  [WARN] -v{name}={val}: defined in %const but never '
                  f"referenced as {{{name}}}; no-op.")
            continue
        any_topology = False
        info_paths = []
        for cat, tok, _ln, src in refs:
            if tok == '<SECTION_HEADER>':
                print(f'  [ERROR] -v{name}={val}: '
                      f"used as section header line ({src.lstrip()[:60]}...). "
                      f"Toggles a TOML section (e.g. enabling [DYN]); "
                      f"cannot be expressed as -v[<path>]=val. "
                      f"Use a ctrlG.<sname>.<tag>.toml variant.")
                n_err += 1
                any_topology = True
                continue
            if (cat, tok) in TOPOLOGY_TOKENS:
                print(f'  [ERROR] -v{name}={val}: '
                      f"changes topology ({cat}_{tok}={val}). "
                      f"TOML schema fixes [[site]]/[[spec]] count; "
                      f"cannot be expressed as -v[<path>]=val. "
                      f"Generate a ctrlG.<sname>.<tag>.toml variant via:\n"
                      f"      Legacy2toml.py <sname> -v{name}={val} ...   # save output\n"
                      f"      mv ctrlG.<sname>.toml ctrlG.<sname>.<tag>.toml")
                n_err += 1
                any_topology = True
                continue
            tp = _toml_path_for(cat, tok, schema)
            info_paths.append((cat, tok, tp))
        if not any_topology and info_paths:
            for cat, tok, tp in info_paths:
                if tp is None:
                    print(f'  [INFO] -v{name}={val}: '
                          f"used at {cat}_{tok}; legacy-only key "
                          f"(no TOML path in schema).")
                else:
                    print(f'  [INFO] -v{name}={val}: '
                          f"used at {cat}_{tok} -> TOML path '{tp}'. "
                          f"Equivalent: -v[{tp}]={val}")
    print('---')
    return n_err

def main():
    if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
        sys.exit(__doc__.strip())
    sname = sys.argv[1]
    extra_args = sys.argv[2:]

    here = Path(os.path.dirname(os.path.realpath(__file__)))

    ctrl_legacy = Path(f'ctrl.{sname}')
    if not ctrl_legacy.exists():
        sys.exit(f'Legacy2toml.py: {ctrl_legacy} not found')

    # Analyze -v overrides BEFORE invoking the converters, so the user can
    # see which overrides are well-defined, which need a variant, and which
    # already have an equivalent TOML-path form.
    n_err = analyze_v_overrides(str(ctrl_legacy), extra_args)
    if n_err > 0:
        print(f'Legacy2toml.py: {n_err} -v override(s) cannot be expressed '
              f'as -v[<path>]=val in default TOML mode.', flush=True)
        print('  Proceeding with conversion anyway -- the resulting '
              'ctrlG.{0}.toml reflects these overrides; save it with a '
              'descriptive suffix (e.g. ctrlG.{0}.<tag>.toml) and switch '
              'between variants via cp in your test.py.'.format(sname),
              flush=True)

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

    # Annotate output with help comments from toml_comments.py
    for fname in (str(out_path), 'PB.toml'):
        if Path(fname).exists():
            txt = Path(fname).read_text()
            Path(fname).write_text(apply_toml_annotations(txt))
    banner(f'wrote ctrlG.{sname}.toml + PB.toml (annotated)')

if __name__ == '__main__':
    main()
