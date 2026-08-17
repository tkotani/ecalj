#!/usr/bin/env python3
"""
ctrl2ctrltoml v2 — schema-based typed TOML output.

Pipeline (steps 1-6 unchanged from v1):
  1. %const eval (math, ^ -> **)
  2. -v command-line overrides (-vlabel=val)
  3. {var} substitution
  4. math eval per token
  5. re-flow continuation lines
  6. tokenize "CAT TOK1=VAL1 TOK2=VAL2 ..." -> "CAT_TOK[@idx] vals" lines
  7. (NEW) lookup ctrl_schema, emit typed sectioned TOML

Unknown keys (CHARGE_*, DOS_*, GW_*, etc.) are dead emission from comments
or obsolete blocks. They are skipped with a warning to stderr.

2026-05-03 T.K. + Claude.
"""
import sys, re, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.toml_comments import fmt_section_header, fmt_key_inline
from math import sin, cos, tan, log, exp, sqrt, pi

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))  # follow symlink
from pylib.ctrl_schema import SCHEMA, SECTION_ORDER, ARRAY_SECTIONS

# ---------------- Steps 1-5 (verbatim from v1, comments trimmed) -------------

instr  = sys.stdin.read()
instrl = instr.split('\n')

labels = []
T, t, F, f = 1, 1, 0, 0
# Two distinct variable layers of the legacy format:
#   % const NAME=...  -> preprocessor vars, referenced ONLY as {NAME}  (pconstns)
#   CONST NAME=...    -> runtime vars, referenced bare in expressions
#                        (e.g. ALAT=a*fv)                              (constns)
# Keeping them separate matters: a % const named 'b' must NOT capture the
# 'b' token of "ITER MIX=B3,b=0.2".
pconstns = {}
constns = {}
pat = r'(\r?\n)|(\r\n?)'
aft = r'\n'

# -v overrides must be visible while scanning (they can flip %ifdef branches)
vconst = {}
for iarg in sys.argv[1:]:
    try:
        vin = iarg.split('-v')[1]
        label, val = vin.split('=')
    except Exception:
        continue
    vconst[label] = val

# Single sequential pass: %ifdef/%ifndef/%if/%else/%endif preprocessing
# (branch conditions evaluate % const values, undefined -> false) combined
# with % const collection. Inactive-branch lines are dropped entirely.
active_lines = []
ifstack = []  # per level: [this_branch_active, some_branch_taken]
for linei in instrl:
    line = re.sub(pat, aft, linei)
    stripped = line.strip()
    m = re.match(r'^%\s*(ifdef|ifndef|if)\b\s*(.*)$', stripped)
    if m:
        parent = all(lv[0] for lv in ifstack)
        cond = False
        if parent:
            try:
                cond = bool(eval(m.group(2).replace('^', '**'), globals(), pconstns))
            except Exception:
                cond = False
        if m.group(1) == 'ifndef':
            cond = not cond
        ifstack.append([cond, cond])
        continue
    if re.match(r'^%\s*elseif\b|^%\s*else\b', stripped):
        if ifstack:
            parent = all(lv[0] for lv in ifstack[:-1])
            ifstack[-1][0] = parent and not ifstack[-1][1]
            ifstack[-1][1] = ifstack[-1][1] or ifstack[-1][0]
        continue
    if re.match(r'^%\s*endif\b', stripped):
        if ifstack:
            ifstack.pop()
        continue
    if not all(lv[0] for lv in ifstack):
        continue  # inside a false branch
    active_lines.append(line)
    if len(line) == 0: continue
    if line[0] != '%': continue
    if 'const' not in line: continue  # skip %show etc.
    line = re.sub(r'=\s+', '=', line)
    constdata0 = (line.split('const')[1]).split('#')[0].split(' ')
    for ix in constdata0:
        if ix == '': continue
        label = ix.split('=')[0]
        if label in vconst:  # -v wins over the ctrl default
            ix = label + '=' + vconst[label]
        ix = ix.replace('^', '**')
        exec(ix, globals(), pconstns)
        labels.append(label)
constrep = {i: str(eval(i, globals(), pconstns)) for i in labels}

for label, val in vconst.items():
    constrep[label] = str(val)

midfile = '\n'.join(active_lines)

def _braceval(m):
    """{name} -> constrep (so -v overrides win); {expr} -> eval, e.g. ALAT={a1/au}."""
    expr = m.group(1)
    if expr in constrep:
        return constrep[expr]
    try:
        return str(eval(expr.replace('^', '**'), globals(), pconstns))
    except Exception:
        return m.group(0)  # leave unresolved braces as-is (caller warns downstream)

midfile = re.sub(r'\{([^{}]+)\}', _braceval, midfile)

# Legacy CONST category: bare-name variables referenced by later expressions
# (e.g. "CONST a=4.977/0.529177 fv=1" then "STRUC ALAT=a*fv"). Collect its
# assignments (indented continuation lines included) into constns so the
# per-token eval below resolves them; the category itself is dropped at the
# tokenize step.
incat = False
for cline in midfile.split('\n'):
    if not cline:
        incat = False
        continue
    if cline[0] not in ' \t':
        parts = cline.split(None, 1)
        incat = parts[0] == 'CONST'
        body = parts[1] if (incat and len(parts) > 1) else ''
    else:
        body = cline if incat else ''
    if not body:
        continue
    body = re.sub(r'=\s+', '=', body).split('#')[0]
    for ix in body.split():
        if '=' not in ix:
            continue
        ix = re.sub(r'(\d)[dD]([-+]?\d)', r'\1e\2', ix.replace('^', '**'))
        try:
            exec(ix, globals(), constns)
        except Exception:
            print(f'# WARNING: CONST: could not evaluate {ix!r}', file=sys.stderr)

outfile = ''
for ilinex in midfile.split('\n'):
    iline = re.sub(pat, aft, ilinex)
    if len(iline) == 0: continue
    if iline[0] == '%': continue
    iii = iline.split('#')[0].split('!')[0].split('%')[0]
    iii = re.sub(r'^\s+$', '', iii)
    if len(iii) == 0: continue
    mmm = re.split(r'([=, ])', iii)
    if 'SYMGR' in mmm[0]:
        iii = re.sub(r'\(', '( ', iii)
        iii = re.sub(r'\)', ' )', iii)
        mmm = re.split(r'([=, ])', iii)
    nnn = []
    for k, ix in enumerate(mmm):
        eout = ix
        if k + 1 < len(mmm) and mmm[k + 1] == '=':
            pass  # KEY position (TOK=...): a CONST var may share the name (e.g. R=R) — never evaluate
        else:
            try:
                eout = eval(ix, globals(), constns)
            except Exception:
                pass
        nnn.append(str(eout).replace(',', ' '))
    outfile += ''.join(nnn) + '\n'
outfile = re.sub(r'\t', ' ', outfile)

# Step 5: re-flow continuation lines
lll = ''
init = False
for line in outfile.split('\n'):
    if len(line) == 0: continue
    if line[0] != ' ':
        ladd = line
        if len(ladd) != 0:
            lll = lll + '\n' + ladd
        init = True
    elif init:
        lll = lll + line.rstrip(' ')
if lll.startswith('\n'):
    lll = lll[1:]
lll = re.sub(r'=\s+', '=', lll)

# Step 6: tokenize each "CAT TOK1=VAL1 TOK2=VAL2 ..." line
# Output: catok = list of "CAT_TOK[@idx] val val ..." entries
catok = ''
nbas = 0
nspec = 0
for iline in lll.split('\n'):
    if not iline: continue
    line = iline
    cat = iline.split(' ')[0]
    if cat == 'CONST':
        continue  # consumed into constns above; not an rval2 category
    if cat == 'ITER':
        line = line.replace(',', ' ')
    line = [i for i in line.split(' ') if i != '']
    if not line: continue
    cat = line[0].split(' ')[0]
    tok = ''
    tokk = []
    id_ = 0
    idx = ''
    line.append('EOL')
    for ix in line[1:]:
        i = ix
        if cat == 'SITE':
            if i[0:4] == 'ATOM':
                id_ += 1
                nbas += 1
        if cat == 'SPEC':
            if i[0:4] == 'ATOM':
                id_ += 1
                nspec += 1
        if i and re.match('[a-zA-Z]', i[0]):
            ic = '_'
            dat = cat + ic + tok
            dat = dat.replace('=', idx + ' ', 1)
            tokk.append(dat)
            idx = ''
            if id_ > 0:
                idx = '@' + str(id_)
            tok = i
        else:
            tok = tok + ' ' + i
    for itk in tokk:
        if itk[-1] != '_' and itk[0:6] != 'SYMGRP':
            catok += itk + '\n'

# ---------------- Step 7 (NEW): schema-driven typed emit ---------------------

def warn(msg):
    print(f'# WARNING: {msg}', file=sys.stderr)

def parse_value(val_str: str, typ: str):
    """Parse a whitespace-separated raw value string per type."""
    # Accept Fortran D-notation (1d-8 -> 1e-8) so legacy ctrl files round-trip.
    raw = re.sub(r'(\d)[dD]([-+]?\d)', r'\1e\2', val_str.strip())
    toks = raw.split()
    if typ == 'str':
        return raw
    if typ in ('int', 'bool'):
        if not toks:
            return None
        v = toks[0]
        # legacy T/F or 0/1
        if v.upper() in ('T','TRUE','.TRUE.'):
            return True if typ=='bool' else 1
        if v.upper() in ('F','FALSE','.FALSE.'):
            return False if typ=='bool' else 0
        try:
            return int(float(v))
        except ValueError:
            return None
    if typ == 'real':
        if not toks: return None
        try:
            return float(toks[0])
        except ValueError:
            return None
    try:
        if typ == 'int_vec':
            return [int(float(t)) for t in toks]
        if typ == 'int_vec3':
            vals = [int(float(t)) for t in toks]
            return (vals + [0,0,0])[:3]
        if typ == 'real_vec':
            return [float(t) for t in toks]
        if typ == 'real_vec3':
            vals = [float(t) for t in toks]
            return (vals + [0.0,0.0,0.0])[:3]
    except ValueError:
        return None  # unresolved symbol etc.; caller warns + rval2 default
    if typ == 'real_mat3x3':
        try:
            vals = [float(t) for t in toks]
        except ValueError:
            return None
        if len(vals) != 9:
            warn(f'real_mat3x3 expected 9 values, got {len(vals)}')
            vals = (vals + [0.0]*9)[:9]
        return [vals[0:3], vals[3:6], vals[6:9]]
    return raw

# ---- accumulate by section ----
top_scalars = {}  # toml_key -> (typ, value, legacy_subkey)
sections    = {}  # section -> dict of toml_key -> (typ, value, legacy_subkey)
arrays      = {}  # 'site' or 'spec' -> {idx_int -> dict of toml_key -> (typ, value, legacy_subkey)}

# Insert SYMGRP / SYMGRPAF top-level entries (these were never tokenized into catok)
for iline in lll.split('\n'):
    if 'SYMGRP' in iline:
        parts = iline.split(None, 1)
        cat = parts[0]
        val = parts[1] if len(parts) > 1 else ''
        spec = SCHEMA.get(cat)
        if spec is None:
            continue
        sec, key, typ, _ = spec
        # collapse whitespace
        val = re.sub(r'\s+', ' ', val).strip()
        if key in top_scalars:
            continue  # first wins (matches rval2 behaviour for duplicate keys)
        top_scalars[key] = (typ, val, spec[3])

unknown_keys = set()
seen_legacy = set()
for iline in catok.split('\n'):
    if not iline.strip():
        continue
    parts = iline.split(None, 1)
    head = parts[0]
    val_raw = parts[1] if len(parts) > 1 else ''
    # split off @N
    m = re.match(r'^([A-Za-z][A-Za-z0-9_/-]*)(?:@(\d+))?$', head)
    if not m:
        unknown_keys.add(head)
        continue
    legacy_full = m.group(1)
    idx_str = m.group(2)
    spec = SCHEMA.get(legacy_full)
    if spec is None:
        unknown_keys.add(legacy_full)
        continue
    if (legacy_full, idx_str) in seen_legacy:
        continue  # first wins
    seen_legacy.add((legacy_full, idx_str))
    sec, key, typ, legacy_sub = spec
    parsed = parse_value(val_raw, typ)
    if parsed is None or (isinstance(parsed, list) and len(parsed) == 0):
        # value couldn't be parsed (e.g. unsubstituted {var}); skip so rval2
        # falls back to its built-in default rather than reading garbage.
        warn(f'{legacy_full}: unparseable value {val_raw!r}; falling back to rval2 default')
        continue
    if sec is None:
        top_scalars[key] = (typ, parsed, legacy_sub)
    elif sec in ARRAY_SECTIONS:
        if idx_str is None:
            warn(f'{legacy_full} has no @N; expected for [[{sec}]]')
            continue
        idx = int(idx_str)
        arrays.setdefault(sec, {}).setdefault(idx, {})[key] = (typ, parsed, legacy_sub)
    else:
        sections.setdefault(sec, {})[key] = (typ, parsed, legacy_sub)

# Auto-fill STRUC_NSPEC / STRUC_NBAS if missing
struc = sections.setdefault('struc', {})
if 'nspec' not in struc:
    struc['nspec'] = ('int', nspec, 'NSPEC')
if 'nbas' not in struc:
    struc['nbas'] = ('int', nbas, 'NBAS')

if unknown_keys:
    warn('skipped dead/unknown keys: ' + ', '.join(sorted(unknown_keys)))

# ---------------- Emit TOML --------------------------------------------------

def fmt_real(x):
    """Format real preserving full double precision (~15 sig figs)."""
    if x is None:
        return '0.0'
    if x == 0:
        return '0.0'
    s = f'{x:.15g}'
    if '.' not in s and 'e' not in s and 'E' not in s:
        s += '.0'
    return s

def fmt_int(x):
    return str(int(x))

def fmt_bool(b):
    return 'true' if b else 'false'

def fmt_str(s):
    s = s.replace('\\', '\\\\').replace('"', '\\"')
    return '"' + s + '"'

def fmt_key(k):
    """TOML key: bare for [A-Za-z0-9_-]+, quoted otherwise (e.g. r/w)."""
    if re.match(r'^[A-Za-z0-9_-]+$', k):
        return k
    return '"' + k + '"'

def fmt_vec(vals, fmtfn, align=True):
    if not vals:
        return '[]'
    strs = [fmtfn(v) for v in vals]
    if align:
        w = max(len(s) for s in strs)
        strs = [s.rjust(w) for s in strs]
    return '[' + ', '.join(strs) + ']'

def fmt_mat3x3(rows, fmtfn):
    """Emit 3-row matrix with comma-aligned columns."""
    cols = list(zip(*rows))
    widths = [max(len(fmtfn(v)) for v in col) for col in cols]
    lines = []
    for i, row in enumerate(rows):
        cells = [fmtfn(v).rjust(widths[j]) for j, v in enumerate(row)]
        line = '[' + ', '.join(cells) + ']'
        if i == 0:
            lines.append('[' + line + ',')
        elif i == len(rows) - 1:
            lines.append(' ' + line + ']')
        else:
            lines.append(' ' + line + ',')
    return '\n        '.join(lines)

def fmt_typed(typ, val):
    if val is None:
        if typ == 'str':       return fmt_str('')
        if typ == 'int':       return '0'
        if typ == 'real':      return '0.0'
        if typ == 'bool':      return 'false'
        if typ.endswith('vec') or typ.endswith('vec3'): return '[]'
        if typ == 'real_mat3x3': return '[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]'
    if typ == 'str':
        return fmt_str(str(val))
    if typ == 'int':
        return fmt_int(val)
    if typ == 'real':
        return fmt_real(val)
    if typ == 'bool':
        if isinstance(val, bool):
            return fmt_bool(val)
        return fmt_bool(bool(int(val)))
    if typ == 'int_vec' or typ == 'int_vec3':
        return fmt_vec(val, fmt_int)
    if typ == 'real_vec' or typ == 'real_vec3':
        return fmt_vec(val, fmt_real)
    if typ == 'real_mat3x3':
        return fmt_mat3x3(val, fmt_real)
    return fmt_str(str(val))

# Order keys within a section by SCHEMA insertion order (stable)
def schema_key_order(section):
    return [k for legacy_full, (sec, k, typ, _) in SCHEMA.items() if sec == section]

print('# ctrlg.<sname>.toml -- auto-generated by ctrl2ctrltoml.py (v2)')
print('# Schema-typed TOML; m_lmfinit reads it via m_ctrl_toml_loader.')
print('# Inline comments come from toml_comments.py (edit there to update).')
print('')
print(fmt_section_header('top'), end='')

# Top-level scalars
top_order = schema_key_order(None)
for k in top_order:
    if k in top_scalars:
        typ, val, _ = top_scalars[k]
        ln = f'{fmt_key(k)} = {fmt_typed(typ, val)}'
        c = fmt_key_inline('top', k)
        if c and '\n' not in ln:
            ln = ln + '  ' + c
        print(ln)
print('')

# [io], [struc] (regular tables before arrays)
def emit_section(sec):
    if sec not in sections:
        return
    hdr = fmt_section_header(sec)
    if hdr:
        print(hdr, end='')
    print(f'[{sec}]')
    order = schema_key_order(sec)
    seen = set()
    keys_present = [k for k in order if k in sections[sec]]
    if not keys_present:
        print('')
        return
    w = max(len(fmt_key(k)) for k in keys_present)
    for k in order:
        if k in sections[sec]:
            typ, val, _ = sections[sec][k]
            ln = f'{fmt_key(k).ljust(w)} = {fmt_typed(typ, val)}'
            c = fmt_key_inline(sec, k)
            if c and '\n' not in ln:
                ln = ln + '  ' + c
            print(ln)
            seen.add(k)
    if sec == 'ham' and 'phispinsym' not in seen:
        # [ham] phispinsym is required by the current Fortran reader
        # (m_lmfinit rval2 nreq=1). Legacy ctrl has no PHISPINSYM token,
        # so emit the default explicitly (mirrors ctrlgenToml.py).
        print('phispinsym  = false  # true = spin-averaged radial wfns (needed for SOC=1 perturbation)')
    print('')

emit_section('io')
emit_section('struc')

# [[site]] / [[spec]]
def emit_array(sec):
    if sec not in arrays:
        return
    items = arrays[sec]
    indices = sorted(items.keys())
    if not indices:
        return
    order = schema_key_order(sec)
    hdr = fmt_section_header(sec)
    if hdr:
        print(hdr, end='')
    for i in indices:
        print(f'[[{sec}]]   # @{i}')
        present = [k for k in order if k in items[i]]
        if present:
            w = max(len(fmt_key(k)) for k in present)
            for k in order:
                if k in items[i]:
                    typ, val, _ = items[i][k]
                    ln = f'{fmt_key(k).ljust(w)} = {fmt_typed(typ, val)}'
                    c = fmt_key_inline(sec, k)
                    if c and '\n' not in ln:
                        ln = ln + '  ' + c
                    print(ln)
        print('')

emit_array('site')
emit_array('spec')

emit_section('bz')
emit_section('iter')
emit_section('ham')
emit_section('options')
emit_section('dyn')
emit_section('str')
emit_section('ewald')
