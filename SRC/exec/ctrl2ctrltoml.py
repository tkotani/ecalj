#!/usr/bin/env python3
# Convert ctrl file to ctrl.<sname>.toml (TOML format), stdin to stdout.
# Replaces the old ctrl2ctrlp.py which produced a flat "KEY value" file.
#
# Behavior preserved:
#   - %const evaluation (math, ^ -> **)
#   - {var} substitution from %const and -v args
#   - -v label=value command-line overrides
#   - Math expressions in numeric tokens (eval)
#   - Per-spec / per-site indexing (@N)
#   - SYMGRP lines printed verbatim (as "SYMGRP" key)
#   - _ATOM@ lines: " 0" -> " F" (relax flag)
#
# Output: TOML, one line per legacy "KEY VALUE" record:
#     "KEY" = "VAL"
# Each value is encoded as a TOML string. The Fortran loader (m_lmfinit)
# concatenates "KEY VAL" back into the recrd(:) array, so rval2's
# search/parse logic is untouched.
#
# 2026-05-03 T.K. + Claude.
import sys, re
from math import sin, cos, tan, log, exp, sqrt, pi

instr  = sys.stdin.read()
instrl = instr.split('\n')

# ---- %const evaluation ---------------------------------------------------
labels = []
T, t, F, f = 1, 1, 0, 0
pat = r'(\r?\n)|(\r\n?)'
aft = r'\n'
for linei in instrl:
    line = re.sub(pat, aft, linei)  # DOS compatible
    if len(line) == 0:    continue
    if line[0] != '%':    continue
    if line[0] == '#':    continue
    line = re.sub(r'=\s+', '=', line)
    constdata0 = (line.split('const')[1]).split('#')[0].split(' ')
    for ix in constdata0:
        if ix == '':
            continue
        ix = ix.replace('^', '**')
        exec(ix)
        labels.append(ix.split('=')[0])
constrep = {i: str(eval(i)) for i in labels}
for i in constrep.keys():
    exec('del ' + str(i))

# ---- -v command-line overrides -------------------------------------------
for iarg in sys.argv[1:]:
    try:
        vin = iarg.split('-v')[1]
        label, val = vin.split('=')
    except Exception:
        continue
    constrep[label] = str(val)

# ---- {var} substitution --------------------------------------------------
midfile = instr
for i, irep in constrep.items():
    midfile = midfile.replace('{' + i + '}', irep)

# ---- Pure-math evaluation of numeric tokens ------------------------------
outfile = ''
for ilinex in midfile.split('\n'):
    iline = re.sub(pat, aft, ilinex)
    if len(iline) == 0:    continue
    if iline[0] == '%':    continue
    iii = iline.split('#')[0].split('!')[0].split('%')[0]
    iii = re.sub(r'^\s+$', '', iii)
    if len(iii) == 0:    continue
    mmm = re.split(r'([=, ])', iii)
    if 'SYMGR' in mmm[0]:
        iii = re.sub(r'\(', '( ', iii)
        iii = re.sub(r'\)', ' )', iii)
        mmm = re.split(r'([=, ])', iii)
    nnn = []
    for ix in mmm:
        eout = ix
        try:
            eout = eval(ix)
        except Exception:
            pass
        nnn.append(str(eout).replace(',', ' '))
    outfile += ''.join(nnn) + '\n'
outfile = re.sub(r'\t', ' ', outfile)

# ---- Re-flow: glue continuation lines onto category line -----------------
lll = ''
init = False
for line in outfile.split('\n'):
    if len(line) == 0:    continue
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

# ---- Tokenize each category line into "CAT_TOK[@idx] value" entries ------
catok = ''
nbas = 0
nspec = 0
for iline in lll.split('\n'):
    if not iline:
        continue
    line = iline
    cat = iline.split(' ')[0]
    if cat == 'ITER':
        line = line.replace(',', ' ')
    line = [i for i in line.split(' ') if i != '']
    if not line:
        continue
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
# Auto-fill STRUC_NSPEC / STRUC_NBAS only if not already in ctrl (TOML must not
# have duplicate keys; legacy ctrlp tolerated dupes but rval2 picked the first).
if not re.search(r'(^|\n)STRUC_NSPEC\b', catok):
    catok += 'STRUC_NSPEC ' + str(nspec) + '\n'
if not re.search(r'(^|\n)STRUC_NBAS\b', catok):
    catok += 'STRUC_NBAS '  + str(nbas)  + '\n'

# ---- Emit TOML -----------------------------------------------------------
def toml_str(s: str) -> str:
    """Encode s as a TOML basic string (single-quoted)."""
    s = s.replace('\\', '\\\\').replace('"', '\\"')
    return '"' + s + '"'

def emit(key: str, val: str):
    """Emit one TOML line. Key is always quoted (it may contain @)."""
    print(f'{toml_str(key)} = {toml_str(val)}')

print('# ctrl.<sname>.toml -- auto-generated by ctrl2ctrltoml.py')
print('# Each entry mirrors a record of the legacy ctrlp file:')
print('#   "KEY" = "VAL1 VAL2 ..."')
print('# m_lmfinit reformats these back to "KEY VAL1 VAL2 ..." for rval2.')
print('')

# Keep insertion order; suppress duplicate keys (legacy rval2 picks first match).
seen = set()
def emit_once(key, val):
    if key in seen:
        return
    seen.add(key)
    # whitespace-normalize the value (rval2 retokenizes anyway, but TOML stays tidy)
    val = re.sub(r'\s+', ' ', val).strip()
    emit(key, val)

# (a) SYMGRP lines (verbatim from lll)
for iline in lll.split('\n'):
    if 'SYMGRP' in iline:
        parts = iline.split(None, 1)
        if len(parts) == 1:
            emit_once(parts[0], '')
        else:
            emit_once(parts[0], parts[1])

# (b) catok lines, with the legacy "_ATOM@ : 0 -> F" substitution
for iline in catok.split('\n'):
    if not iline.strip():
        continue
    ilinex = iline
    if re.search('_ATOM@', iline):
        ilinex = re.sub(r' 0', ' F', iline)
    parts = ilinex.split(None, 1)
    if not parts:
        continue
    if len(parts) == 1:
        emit_once(parts[0], '')
    else:
        emit_once(parts[0], parts[1])
