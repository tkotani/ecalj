#!/usr/bin/env python3
"""
ctrlgenToml.py — generate ctrlG.<ext>.toml + PB.toml from ctrls.<ext>.

Companion to ctrlgenM1.py (which generates a legacy ctrl template).
Same ctrls.<ext> input and the same periodic-table atomlist defaults,
but emits the structured TOML pair that today's Fortran reads directly
(ctrlG.<ext>.toml via m_ctrl_toml_loader; PB.toml has the per-atom
product-basis tables — GW path only, normally not edited by hand).

Recommended usage (no flags):

    ctrlgenToml.py <ext>

Just run with the sname; everything else is sensible defaults baked
into ctrlG.<ext>.toml. The file is fully commented and TOML-typed —
**you are expected to edit ctrlG.<ext>.toml afterwards** (k-mesh,
nspin, so, xcfun, scaledsigma, [[spec]] U/J, ...). The CLI flags
below are bake-in equivalents of those edits, useful for scripting
or batch generation; for normal interactive use, leave them off.

Optional flags (all bake into ctrlG.<ext>.toml):
    --nspin=N            1 nonmag, 2 spin-polarized
    --so=N               0 none, 1 L.S, 2 Lz.Sz (forces nspin=2)
    --xcfun=vwn|bh|pbe   exchange-correlation
    --mmom='MMOM=...'    initial magnetic moments
    --nk1=N --nk2=N --nk3=N   LDA k-mesh (default 8 8 8)
    --insulator          metal=0 (default 3 with tetra)
    --systype=bulk|molecule    molecule disables tetra, biases FSMOM
    --fsmom=F            fixed-spin moment
    --ssig=F             ScaledSigma (1.0 = full QSGW, 0.8 = QSGW80)
    --tratio=F           sphere-touching ratio (default 0.97)
    --ehmol              EH=-1 group only (smaller basis for molecules)
    --skipgw             generate ctrlG.<ext>.toml only (no [gw] /
                         [product_basis] / [blocks]; no PB.toml).
                         Use for DFT-only / no-GW workflows.
    --addgw              append [gw] / [product_basis] / [blocks] (and
                         emit PB.toml) to an *existing* ctrlG.<ext>.toml
                         WITHOUT regenerating the ctrl-side keys —
                         preserves any hand-edits you made to [bz] /
                         [ham] / [[spec]] etc.  Errors out if [gw] is
                         already present (delete the GW tables and
                         PB.toml first, or omit --addgw to fully
                         regenerate from ctrls.<ext>).
    --showatomlist       print the periodic-table defaults table and exit
    -h, --help           show this docstring

Inputs:
    ctrls.<ext> — STRUC + SITE [+ optional SPEC] skeleton

Outputs:
    ctrlG.<ext>.toml — schema-typed TOML, ready for lmfa/lmf/gwsc.
    PB.toml          — per-atom product-basis tables (GW only;
                       normally not edited by hand).

Compatibility note: ctrlgenM1.py + ctrl2ctrltoml.py is preserved as
the legacy path. As of 2026-05 the Fortran reads ctrlG.<ext>.toml +
PB.toml only; ctrlgenM1.py exits immediately with a pointer to this
script.

2026-05-03 T.K. + Claude.
"""
import os, sys, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.toml_comments import fmt_section_header, fmt_key_inline, apply_toml_annotations

# ---------------------------------------------------------------------------
# Atomlist — extracted at runtime from ctrlgenM1.py so there's a single
# source of truth for the periodic-table defaults.
# ---------------------------------------------------------------------------

def _load_atomlist():
    here = os.path.dirname(os.path.realpath(__file__))
    m1path = os.path.join(here, 'ctrlgenM1.py')
    if not os.path.exists(m1path):
        sys.exit('ctrlgenToml: cannot find ctrlgenM1.py at ' + m1path)
    src = open(m1path).read()
    m = re.search(r'atomlist\s*=\s*"""(.*?)"""', src, re.DOTALL)
    if not m:
        sys.exit('ctrlgenToml: failed to extract atomlist from ctrlgenM1.py')
    return m.group(1)

atomlist = _load_atomlist()

# Build dicatom: 'Si' -> ' atomz=14@ eh=-1*4@ eh2=-2*3@ R=1.15@'
dicatom = {}
for line in atomlist.split('\n'):
    s = line.strip()
    if not s or s.startswith('#'):
        continue
    if '=' not in s:
        continue
    head, _, body = s.partition('=')
    sym = head.strip()
    body = body.strip()
    # Skip if not in the "Sym = "..."" form
    if not (body.startswith('"') and 'atomz=' in body):
        continue
    # Strip leading "  and trailing  " (or "@")
    body = body.lstrip('"').rstrip()
    if body.endswith('"'):
        body = body[:-1].rstrip()
    dicatom[sym] = body

z2sym = {}
for sym, body in dicatom.items():
    try:
        z = body.split('atomz=')[1].split('@')[0].strip()
        z2sym[z] = sym
    except IndexError:
        pass

# ---------------------------------------------------------------------------
# CLI parsing (mirrors ctrlgenM1.manip_argset)
# ---------------------------------------------------------------------------

def parse_args(argv):
    opts = dict(
        ext=None,
        nspin='1', so='0', xcfun_str='vwn',
        mmom='#MMOM=0 0 1 0',
        systype='bulk',
        nk1='8', nk2=None, nk3=None,
        metali=3,
        fsmom=0.0, ssig=1.0,
        touchingratio=0.97,
        eh1set=1,
        skipgw=False,
        addgw=False,
    )
    pos_args = []
    for a in argv[1:]:
        if a.startswith('--nspin='):       opts['nspin'] = a.split('=',1)[1]
        elif a.startswith('--so='):
            opts['so'] = a.split('=',1)[1]; opts['nspin'] = '2'
        elif a.startswith('--xcfun='):     opts['xcfun_str'] = a.split('=',1)[1]
        elif a.startswith('--mmom='):      opts['mmom'] = a.split('=',1)[1]
        elif a == '--insulator':           opts['metali'] = 0
        elif a.startswith('--nk1='):       opts['nk1'] = a.split('=',1)[1]
        elif a.startswith('--nk2='):       opts['nk2'] = a.split('=',1)[1]
        elif a.startswith('--nk3='):       opts['nk3'] = a.split('=',1)[1]
        elif a.startswith('--systype='):   opts['systype'] = a.split('=',1)[1]
        elif a.startswith('--fsmom='):     opts['fsmom'] = float(a.split('=',1)[1])
        elif a.startswith('--ssig='):      opts['ssig'] = float(a.split('=',1)[1])
        elif a.startswith('--tratio='):    opts['touchingratio'] = float(a.split('=',1)[1])
        elif a == '--ehmol':               opts['eh1set'] = 0
        elif a == '--skipgw':             opts['skipgw'] = True
        elif a == '--addgw':              opts['addgw'] = True
        elif a == '--showatomlist':
            # Print the atomlist (periodic-table defaults) and exit;
            # no <ext> required.  The atomlist is shared with the
            # legacy ctrlgenM1.py via runtime regex extraction.
            print(atomlist)
            sys.exit(0)
        elif a in ('-h', '--help'):
            sys.exit(__doc__.strip())
        elif a.startswith('-'):
            sys.exit('ctrlgenToml: unknown option ' + a)
        else:
            pos_args.append(a)
    if not pos_args:
        sys.exit(__doc__.strip())
    opts['ext'] = pos_args[0]
    if opts['nk2'] is None: opts['nk2'] = opts['nk1']
    if opts['nk3'] is None: opts['nk3'] = opts['nk1']
    if opts['xcfun_str'].upper() == 'PBE':  opts['xcfun'] = 103
    elif opts['xcfun_str'].upper() == 'VWN': opts['xcfun'] = 1
    elif opts['xcfun_str'].upper() == 'BH':  opts['xcfun'] = 2
    else: sys.exit('ctrlgenToml: unknown xcfun ' + opts['xcfun_str'])
    if opts['systype'].upper() == 'MOLECULE':
        opts['nk1'] = opts['nk2'] = opts['nk3'] = '1'
    elif opts['systype'].upper() != 'BULK':
        sys.exit('ctrlgenToml: unknown systype ' + opts['systype'])
    if opts['so'] != '0' and opts['nspin'] != '2':
        sys.exit('ctrlgenToml: so/=0 requires nspin=2')
    return opts

# ---------------------------------------------------------------------------
# ctrls.<ext> parser (lifted from ctrlgenM1.lineReadfile / GetCat)
# ---------------------------------------------------------------------------

def line_read(filename):
    out = []
    for line in open(filename):
        line = line.rstrip('\n')
        if line:
            out.append(line)
    return out

def get_category(lines, key):
    """Return list of lines forming the category starting with KEY."""
    res = []
    state = 0
    pat = re.compile(r'^(' + key.upper() + '|' + key.lower() + r')(\s|\Z)')
    for ln in lines:
        if state == 1 and re.match(r'^\w', ln):
            break
        if pat.match(ln):
            state = 1
        if state == 1:
            res.append(ln)
    return res

def site_atom_names(site_lines):
    names = []
    for ln in site_lines:
        for m in re.finditer(r'\bATOM=(\S+)', ln):
            names.append(m.group(1))
    return names

def uniq(lst):
    seen, out = set(), []
    for x in lst:
        if x not in seen:
            seen.add(x); out.append(x)
    return out

def get_field(body, key):
    """body = ' atomz=14@ eh=-1*3@ ...'; key='eh=' returns '-1*3 '."""
    return body.split(key, 1)[1].split('@', 1)[0].strip() + ' '

def get_quoted(body, key):
    """Pull pz='PZ=0,3.9' out of body."""
    return body.split(key + "'", 1)[1].split("'@", 1)[0].strip() + ' '

# ---------------------------------------------------------------------------
# %const evaluation (subset shared with ctrl2ctrltoml.py step 1)
# ---------------------------------------------------------------------------

from math import sin, cos, tan, log, exp, sqrt, pi
T_, t_, F_, f_ = 1, 1, 0, 0  # legacy bool aliases

def eval_consts_and_substitute(text, extra_vars=None):
    """Evaluate %const lines, substitute {var}, evaluate inline math tokens."""
    consts = {}
    if extra_vars:
        consts.update(extra_vars)
    out_lines = []
    for raw in text.split('\n'):
        if not raw.strip(): continue
        if raw.lstrip().startswith('#'): continue
        if raw.lstrip().startswith('%'):
            line = re.sub(r'=\s+', '=', raw)
            after = line.split('const', 1)[1].split('#', 1)[0]
            for tok in after.split():
                if not tok: continue
                tok = tok.replace('^', '**')
                try:
                    name, expr = tok.split('=', 1)
                    consts[name] = eval(expr, {"__builtins__": {}}, dict(consts, **_math_env))
                except Exception:
                    pass
            continue
        out_lines.append(raw)
    body = '\n'.join(out_lines)
    # substitute {var}
    for name, val in consts.items():
        body = body.replace('{' + name + '}', str(val))
    # evaluate numeric tokens
    eval_lines = []
    for ln in body.split('\n'):
        line = ln.split('#', 1)[0].split('!', 1)[0]
        if not line.strip():
            eval_lines.append('')
            continue
        toks = re.split(r'([=, ])', line)
        new = []
        skip_next = False
        for i, t in enumerate(toks):
            # Don't eval the token right after ATOM= (element name like F, N, etc.)
            if skip_next:
                new.append(t)
                skip_next = False
                continue
            if t == '=' and i > 0 and toks[i-1].rstrip().upper().endswith('ATOM'):
                new.append(t)
                skip_next = True
                continue
            v = t
            try:
                v = eval(t, {"__builtins__": {}}, _math_env)
            except Exception:
                pass
            new.append(str(v).replace(',', ' '))
        eval_lines.append(''.join(new))
    return '\n'.join(eval_lines), consts

_math_env = dict(sin=sin, cos=cos, tan=tan, log=log, exp=exp, sqrt=sqrt, pi=pi,
                 T=1, t=1, F=0, f=0)

# ---------------------------------------------------------------------------
# TOML emit helpers (mirror ctrl2ctrltoml v2)
# ---------------------------------------------------------------------------

def fmt_real(x):
    if x is None or x == 0: return '0.0'
    s = f'{float(x):.15g}'
    if '.' not in s and 'e' not in s and 'E' not in s: s += '.0'
    return s

def fmt_int(x): return str(int(x))

def fmt_bool(b): return 'true' if b else 'false'

def fmt_str(s):
    s = str(s).replace('\\', '\\\\').replace('"', '\\"')
    return '"' + s + '"'

def fmt_real_vec(vals):
    if not vals: return '[]'
    strs = [fmt_real(v) for v in vals]
    w = max(len(s) for s in strs)
    return '[' + ', '.join(s.rjust(w) for s in strs) + ']'

def fmt_int_vec(vals):
    if not vals: return '[]'
    strs = [fmt_int(v) for v in vals]
    w = max(len(s) for s in strs)
    return '[' + ', '.join(s.rjust(w) for s in strs) + ']'

def fmt_mat3x3(rows):
    cols = list(zip(*rows))
    widths = [max(len(fmt_real(v)) for v in col) for col in cols]
    lines = []
    for i, row in enumerate(rows):
        cells = [fmt_real(v).rjust(widths[j]) for j, v in enumerate(row)]
        ln = '[' + ', '.join(cells) + ']'
        if i == 0:           lines.append('[' + ln + ',')
        elif i == len(rows)-1: lines.append(' ' + ln + ']')
        else:                  lines.append(' ' + ln + ',')
    return '\n        '.join(lines)

# ---------------------------------------------------------------------------
# R upper limits (mirror ctrlgenM1)
# ---------------------------------------------------------------------------

def r_upper_limit(z, r):
    z = float(z)
    if r > 2.7 and 2.8 < z < 4.2:    r = 2.7
    if r > 2.4 and 10.8 < z < 12.2:  r = 2.4
    if r > 2.6 and 18.8 < z < 20.2:  r = 2.6
    if r > 2.8 and 36.8 < z < 38.2:  r = 2.8
    if r > 2.8 and 54.8 < z < 56.2:  r = 2.8
    if r > 3.0:                      r = 3.0
    return r

def lmxa_for_z(z):
    z = float(z)
    if z < 1.5:                  return 3
    if 56.1 < z < 71.9:          return 6
    return 4

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    opts = parse_args(sys.argv)
    ext = opts['ext']

    if opts.get('skipgw') and opts.get('addgw'):
        sys.exit('ctrlgenToml: --skipgw and --addgw are mutually exclusive')

    # --- --addgw: append [gw]/[product_basis]/[blocks] to an existing
    # ctrlG.<ext>.toml without regenerating the ctrl-side keys. ---
    if opts.get('addgw'):
        out_path = 'ctrlG.' + ext + '.toml'
        if not os.path.exists(out_path):
            sys.exit(f'ctrlgenToml --addgw: {out_path} does not exist; run '
                     f'`ctrlgenToml.py {ext}` first to generate it.')
        existing = open(out_path).read()
        # Match a top-level [gw] section header (ignore [[gw]] arrays;
        # ignore in-string accidental occurrences only loosely — TOML
        # does not allow [gw] inside a multiline string at column 0
        # in well-formed input).
        if re.search(r'(?m)^\s*\[gw\]\s*$', existing):
            sys.exit(f'ctrlgenToml --addgw: {out_path} already has a [gw] '
                     f'section.  Refusing to append (delete the existing '
                     f'[gw]/[product_basis]/[blocks] tables and PB.toml '
                     f'first, or just re-run `ctrlgenToml.py {ext}` '
                     f'without --addgw to regenerate from ctrls.{ext}).')
        _run_gw_append(ext, out_path)
        return

    ctrls_path = 'ctrls.' + ext
    if not os.path.exists(ctrls_path):
        sys.exit('ctrlgenToml: ' + ctrls_path + ' not found')

    # Pre-eval %const and {var} substitution + inline math so values are concrete
    raw_text = open(ctrls_path).read()
    evaled, consts = eval_consts_and_substitute(raw_text)
    eval_lines = [ln for ln in evaled.split('\n') if ln.strip()]

    # Re-flow continuation lines (indent = continuation)
    lll = ''
    init = False
    for line in eval_lines:
        if line[0] != ' ':
            if line.strip(): lll += '\n' + line
            init = True
        elif init:
            lll += line.rstrip(' ')
    if lll.startswith('\n'): lll = lll[1:]
    lll = re.sub(r'=\s+', '=', lll)
    flow_lines = [l for l in lll.split('\n') if l.strip()]

    listsite = get_category(flow_lines, 'SITE')
    listspec = get_category(flow_lines, 'SPEC')
    liststruc = get_category(flow_lines, 'STRUC')

    sitenames = site_atom_names(listsite)
    specnames = uniq(sitenames)

    # Build spec2z: explicit from listspec ATOM= ... Z=, else lookup atomlist
    spec2z = {}
    spec_extra = {}
    if listspec:
        for ln in listspec:
            for atomtok in ln.split('ATOM=')[1:]:
                name = atomtok.split()[0]
                m = re.search(r'\bZ\s*=\s*([0-9]+)', atomtok)
                if m:
                    spec2z[name] = m.group(1)
                # take rest after "Z=NN " for extra tokens user wrote (MMOM, IDU, etc.)
                rest = re.split(r'Z\s*=\s*[0-9]+\s', atomtok, 1)
                if len(rest) > 1:
                    spec_extra[name] = rest[1].strip()
    # Fill in from atomlist for missing
    for sym in specnames:
        if sym in spec2z: continue
        if sym in dicatom:
            spec2z[sym] = dicatom[sym].split('atomz=')[1].split('@')[0].strip()

    # ---- run lmchk --getwsr to get touching MT radii ----
    rdic = {}
    if opts['touchingratio'] > 0:
        # Build a tmp ctrl with placeholder R= for spec section
        head = ('IO     VERBOS=35 TIM=0,0\n'
                'SYMGRP find\n')
        specsec0 = 'SPEC\n'
        for sym in specnames:
            z = spec2z.get(sym, '?')
            specsec0 += f'    ATOM={sym} Z={z}\n'
        glist = lambda L: '\n'.join(L) + '\n'
        no_categories = [l for l in flow_lines
                         if not re.match(r'^(SITE|SPEC|STRUC)(\s|\Z)', l)]
        # Build a minimal ctrlG.tmp.toml for the new TOML-only lmchk
        _alat = None; _plat = None
        for ln in liststruc + flow_lines:
            if 'ALAT=' in ln and _alat is None:
                try: _alat = float(ln.split('ALAT=', 1)[1].split()[0])
                except (ValueError, IndexError): pass
            if 'PLAT=' in ln and _plat is None:
                tail = ln.split('PLAT=', 1)[1].split('#', 1)[0]
                nums = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', tail)
                if len(nums) >= 9:
                    vals = [float(x) for x in nums[:9]]
                    _plat = [vals[0:3], vals[3:6], vals[6:9]]
        tmp_toml = 'symgrp = "find"\n\n[io]\nverbose = 35\n\n'
        # nspec/nbas omitted: synthesized by the TOML loader from
        # [[spec]] / [[site]] array lengths (m_ctrl_toml_loader.f90).
        tmp_toml += f'[struc]\nalat = {fmt_real(_alat)}\n'
        tmp_toml += f'plat = {fmt_mat3x3(_plat)}\n\n'
        _si = 0
        for ln in listsite:
            for atomtok in ln.split('ATOM=')[1:]:
                _si += 1
                _name = atomtok.split()[0]
                _pos = None
                _mpos = re.search(r'POS\s*=\s*([^A-Z]*?)(?=[A-Z]|$)', atomtok)
                if _mpos:
                    _nums = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', _mpos.group(1))
                    if len(_nums) >= 3:
                        _pos = [float(x) for x in _nums[:3]]
                tmp_toml += f'[[site]]\natom = "{_name}"\n'
                if _pos: tmp_toml += f'pos = [{fmt_real(_pos[0])}, {fmt_real(_pos[1])}, {fmt_real(_pos[2])}]\n'
                tmp_toml += '\n'
        for sym in specnames:
            z = spec2z.get(sym, '?')
            tmp_toml += f'[[spec]]\natom = "{sym}"\nz = {z}\n\n'
        with open('ctrlG.tmp.toml', 'wt') as f:
            f.write(tmp_toml)
        os.system('mpirun -np 1 lmchk --getwsr tmp > llmchk_getwsr 2>&1')
        try:
            for ln in open('rmt.tmp'):
                parts = ln.split()
                if len(parts) >= 2:
                    rdic[parts[0]] = float(parts[1])
        except FileNotFoundError:
            sys.exit('ctrlgenToml: lmchk --getwsr failed; see llmchk_getwsr')

    # ---- Compose TOML output ----
    out = []
    out.append(f'# ctrl.{ext}.toml — auto-generated by ctrlgenToml.py from ctrls.{ext}')
    out.append('')
    out.append('symgrp = "find"')
    out.append('')
    out.append('[io]')
    out.append('verbose = 35')
    out.append('time    = [0, 0]')
    out.append('')

    # [struc] — ALAT, PLAT from listsite/liststruc evaled
    alat = None; plat = None
    for ln in liststruc + flow_lines:
        if 'ALAT=' in ln and alat is None:
            try: alat = float(ln.split('ALAT=', 1)[1].split()[0])
            except (ValueError, IndexError): pass
        if 'PLAT=' in ln and plat is None:
            tail = ln.split('PLAT=', 1)[1].split('#', 1)[0]
            nums = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', tail)
            if len(nums) >= 9:
                vals = [float(x) for x in nums[:9]]
                plat = [vals[0:3], vals[3:6], vals[6:9]]
    if alat is None: sys.exit('ctrlgenToml: STRUC ALAT not found')
    if plat is None: sys.exit('ctrlgenToml: STRUC PLAT not found')
    out.append('[struc]')
    # nspec/nbas omitted: synthesized by the TOML loader from
    # [[spec]] / [[site]] array lengths (m_ctrl_toml_loader.f90).
    out.append(f'alat  = {fmt_real(alat)}')
    out.append(f'plat  = {fmt_mat3x3(plat)}')
    out.append('')

    # [[site]] entries
    site_idx = 0
    for ln in listsite:
        for atomtok in ln.split('ATOM=')[1:]:
            site_idx += 1
            name = atomtok.split()[0]
            pos = None
            mpos = re.search(r'POS\s*=\s*([^A-Z]*?)(?=[A-Z]|$)', atomtok)
            if mpos:
                nums = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', mpos.group(1))
                if len(nums) >= 3:
                    pos = [float(x) for x in nums[:3]]
            af = None
            maf = re.search(r'\bAF\s*=\s*(-?\d+)', atomtok)
            if maf: af = int(maf.group(1))
            out.append(f'[[site]]   # @{site_idx}')
            out.append(f'atom = {fmt_str(name)}')
            if pos is not None:
                out.append(f'pos  = {fmt_real_vec(pos)}')
            if af is not None and af != 0:
                out.append(f'af   = {af}')
            out.append('')

    # [[spec]] entries
    rsma_factor = None  # not set; default (rval2 default 0.4*rmt)
    for spec_idx, sym in enumerate(specnames, 1):
        if sym in dicatom:
            body = dicatom[sym]
        elif sym in spec2z and spec2z[sym] in z2sym:
            # Custom name (e.g. Niup, Nidn) with explicit Z= in ctrls SPEC:
            # fall back to the periodic-table entry for that Z.
            body = dicatom[z2sym[spec2z[sym]]]
        else:
            sys.exit(f'ctrlgenToml: {sym} not in atomlist; add it or set ctrls SPEC explicitly')
        z = spec2z.get(sym, body.split('atomz=')[1].split('@')[0].strip())
        # R: from lmchk --getwsr × touchingratio, or atomlist R=
        if opts['touchingratio'] > 0 and sym in rdic:
            r_au = rdic[sym] * opts['touchingratio']  # already in a.u.
            r = r_au
        else:
            r_ang = float(get_field(body, 'R=').split('*')[0].strip().rstrip().replace('?', '0'))
            r = r_ang / 0.529177  # Angstrom -> a.u.
        r = r_upper_limit(z, r)
        r = round(r, 2)
        rh = max(0.5, r / 2.0)  # rsmh = (displayed r) / 2 exactly

        eh_data = get_field(body, 'eh=').strip()
        eh_val = float(eh_data.split('*')[0])
        eh_count = int(eh_data.split('*')[1])
        eh2_data = get_field(body, 'eh2=').strip()
        eh2_val = float(eh2_data.split('*')[0])
        eh2_count = int(eh2_data.split('*')[1])

        if opts['eh1set'] == 1: eh_val = -1.0

        lmxa = lmxa_for_z(z)

        out.append(f'[[spec]]   # @{spec_idx}')
        out.append(f'atom   = {fmt_str(sym)}')
        out.append(f'z      = {fmt_int(z)}')
        out.append(f'r      = {fmt_real(r)}')
        out.append(f'eh     = {fmt_real_vec([eh_val]*eh_count)}')
        out.append(f'rsmh   = {fmt_real_vec([rh]*eh_count)}')
        out.append(f'eh2    = {fmt_real_vec([eh2_val]*eh2_count)}')
        out.append(f'rsmh2  = {fmt_real_vec([rh]*eh2_count)}')
        out.append(f'lmxa   = {lmxa}')
        out.append(f'kmxa   = 5')
        out.append(f'nmcore = 1')
        # PZ / P from atomlist (semi-core orbitals). Skip if commented out (#-prefix).
        for src_tag, toml_key in [('pz=', 'pz'), ('p=', 'p')]:
            try:
                raw = get_quoted(body, src_tag).strip()
            except IndexError:
                continue
            # raw looks like "PZ=0,3.9" or "#PZ=0,3.9" or "P=0,0,4.5"
            if raw.startswith('#'):
                continue  # commented out in atomlist; no-op
            m = re.match(r'(PZ|P|Q)\s*=\s*([-0-9.,\s]+)', raw, re.IGNORECASE)
            if not m:
                continue
            vals = [float(v) for v in re.split(r'[, ]+', m.group(2).strip()) if v]
            out.append(f'{toml_key:<6} = {fmt_real_vec(vals)}')
        # ext='IDU=... UH=... JH=... MMOM=... Q=...' for rare-earths / LDA+U cases.
        try:
            ext_str = get_quoted(body, 'ext=').strip()
        except IndexError:
            ext_str = ''
        if ext_str and not ext_str.startswith('#'):
            for key in ('IDU', 'UH', 'JH', 'Q', 'MMOM'):
                m = re.search(rf'\b{key}\s*=\s*([-0-9.,\s]+?)(?=\b[A-Z]+\s*=|$)', ext_str)
                if not m:
                    continue
                vals_raw = m.group(1).strip()
                vals = [v for v in re.split(r'[, ]+', vals_raw) if v]
                if not vals:
                    continue
                if key == 'IDU':
                    out.append(f'idu    = {fmt_int_vec([int(float(v)) for v in vals])}')
                elif key == 'MMOM':
                    out.append(f'mmom   = {fmt_real_vec([float(v) for v in vals])}')
                elif key == 'Q':
                    out.append(f'q      = {fmt_real_vec([float(v) for v in vals])}')
                else:
                    out.append(f'{key.lower():<6} = {fmt_real_vec([float(v) for v in vals])}')
        # Optional initial MMOM from CLI (e.g. --mmom='MMOM=0 0 2 0')
        if not opts['mmom'].startswith('#'):
            mvals = re.findall(r'[-+]?\d+\.?\d*', opts['mmom'])
            if mvals:
                out.append(f'mmom   = {fmt_real_vec([float(x) for x in mvals[:lmxa+1]])}')
        out.append('')

    # [bz] — from CLI options, all concrete values
    out.append('[bz]')
    out.append(f'nkabc  = [{opts["nk1"]}, {opts["nk2"]}, {opts["nk3"]}]')
    out.append(f'metal  = {opts["metali"]}')
    out.append('tetra  = true')
    if opts['systype'].upper() == 'MOLECULE':
        out.append('n      = -1')
        out.append('w      = 0.01')
    out.append('savdos = true')
    out.append('# dosmax = 1.5     # uncomment for DOS plot range')
    out.append('# npts   = 2001    # uncomment for DOS resolution')
    if opts['fsmom']:
        out.append(f'fsmom  = {fmt_real(opts["fsmom"])}')
    out.append('')

    # [iter]
    out.append('[iter]')
    out.append('mix   = "B3"')
    out.append('b     = 0.2')
    out.append('umix  = 0.2')
    out.append('nit   = 80')
    out.append('conv  = 1.0e-5')
    out.append('convc = 1.0e-5')
    out.append('')

    # [ham]
    out.append('[ham]')
    out.append(f'nspin       = {opts["nspin"]}')
    out.append('forces      = 0')
    out.append('gmax        = 12.0')
    out.append('rel         = true')
    out.append(f'xcfun       = {opts["xcfun"]}')
    out.append('pwmode      = 11')
    out.append('pwemax      = 3.0')
    out.append('readp       = true')
    out.append('pnufix      = true')
    out.append('frzwf       = false')
    out.append('phispinsym  = false  # true = spin-averaged radial wfns (needed for SOC=1 perturbation)')
    out.append(f'scaledsigma = {fmt_real(opts["ssig"])}')
    out.append(f'so          = {opts["so"]}')
    out.append('oveps       = 1.0e-8')
    out.append('# rdsig = 12   # uncomment when sigm.<ext> exists (QSGW)')
    out.append('')

    out.append('[options]')
    out.append('hf = false')
    out.append('')

    out_path = 'ctrlG.' + ext + '.toml'

    # Backup any pre-existing ctrlG.<ext>.toml / PB.toml to *.bakup before
    # we write fresh ones. One-generation backup (we overwrite previous .bakup).
    import shutil
    for p in (out_path, 'PB.toml'):
        if os.path.exists(p):
            shutil.move(p, p + '.bakup')
            print(f'ctrlgenToml: backed up existing {p} -> {p}.bakup')

    text = apply_toml_annotations('\n'.join(out))
    with open(out_path, 'wt') as f:
        f.write(text)
    print(f'ctrlgenToml: wrote {out_path} ({len(specnames)} spec, {len(sitenames)} sites)')

    # cleanup tmp
    for p in ('ctrl.tmp', 'ctrlG.tmp.toml', 'rmt.tmp', 'llmchk_getwsr', 'exitcode'):
        try: os.unlink(p)
        except OSError: pass

    # ---------------------------------------------------------------------
    # Append [gw] / [product_basis] / [blocks] sections + emit PB.toml
    # by running the standard pipeline:   lmfa -> lmf --jobgw=0 -> gwinit
    # ---------------------------------------------------------------------
    if opts.get('skipgw'):
        print(f'ctrlgenToml: --skipgw, leaving ctrlG.{ext}.toml without [gw]/[product_basis]/[blocks].  PB.toml not generated.')
        return
    _run_gw_append(ext, out_path)


def _run_gw_append(ext, out_path):
    """Run lmfa -> lmf --jobgw=0 -> gwinit to append [gw] / [product_basis]
    / [blocks] to ctrlG.<ext>.toml (in place) and write PB.toml.  Used by
    the default flow and by --addgw."""
    import subprocess
    from pathlib import Path
    here = Path(__file__).parent  # do NOT resolve symlinks; ~/bin/lmfa etc. live here
    print(f'ctrlgenToml: running lmfa -> lmf --jobgw=0 -> gwinit  to fill GW sections')
    rc = subprocess.run(['mpirun', '-np', '1', str(here / 'lmfa'), ext],
                        stdout=open('llmfa', 'wt'), stderr=subprocess.STDOUT).returncode
    if rc != 0:
        sys.exit(f'ctrlgenToml: lmfa failed (rc={rc}); see llmfa')
    rc = subprocess.run(['mpirun', '-np', '1', str(here / 'lmf'),
                         '--jobgw=0', ext],
                        stdout=open('llmfgw00', 'wt'), stderr=subprocess.STDOUT).returncode
    if rc != 0:
        sys.exit(f'ctrlgenToml: lmf --jobgw=0 failed (rc={rc}); see llmfgw00')
    rc = subprocess.run(['mpirun', '-np', '1', str(here / 'gwinit'), ext]).returncode
    if rc != 0:
        sys.exit(f'ctrlgenToml: gwinit failed (rc={rc})')
    # Re-apply annotations: gwinit appended raw [gw]/[product_basis]/[blocks];
    # also wrote PB.toml. Both need section headers from toml_comments.
    for fname in (out_path, 'PB.toml'):
        if os.path.exists(fname):
            txt = open(fname).read()
            open(fname, 'w').write(apply_toml_annotations(txt))
    print(f'ctrlgenToml: done. {out_path} has [io]/[struc]/[[site]]/[[spec]]/...')
    print(f'             plus [gw]/[product_basis]/[blocks].  PB.toml has nlx/valence/core.')


if __name__ == '__main__':
    main()
