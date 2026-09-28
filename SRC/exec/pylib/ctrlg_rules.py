"""Bring a ctrlg.<sname>.toml up to the current rules (used by ctrlg_update.py and the legacy converters).

Rules (2026-09-28):
  1. [struc] carries no nbas / nspec: the numbers of sites and species are those of the [[site]] / [[spec]] tables
     (lmf stops on them). drop_counts() removes them; when nbas was smaller than the tables (the old way to use only
     the first sites) the [[site]] blocks beyond it are commented out, so the calculation stays the same.
  2. [gw] carries t_tetrakbt (K), always: 0 = no smearing of chi0, > 0 = finite-T tetrahedron at T, < 0 = Im chi0
     smeared by a Gaussian as wide as the Fermi-Dirac distribution at |T| (std = pi kB|T|/sqrt3), which replaces
     [gw] SmearX0 (no longer read). gw_tetrakbt() converts SmearX0 = s to t_tetrakbt = -T with the same width and
     adds t_tetrakbt = 0 (the former default) when it is missing, so the calculation stays the same.
"""
import math
import re

KB_EV = 8.6171e-5          # the constants of m_GWinput / m_tetrakbt
RY_EV = 13.605693

HEAD = re.compile(r'^\[')                    # any table header: [struc], [[site]], [bz], ...
COUNT = re.compile(r'^\s*(nbas|nspec)\s*=\s*(\d+)')
OLDNOTE = '# nbas/nspec are filled automatically from the [[site]]/[[spec]] arrays.'
NEWNOTE = ['# The numbers of sites and species are those of the [[site]]/[[spec]] tables: do not write',
           '# nbas/nspec (lmf stops on them). To leave a site out, comment out its [[site]] block.']


def smearx0_to_kelvin(sx):
    """SmearX0 (Ha, Gaussian std) -> T (K) with pi kB T / sqrt3 = sx; rounded to 0.1 K."""
    return round(sx * 2 * RY_EV / (math.pi / math.sqrt(3) * KB_EV), 1)


def _section(L, name):
    try:
        i0 = next(i for i, l in enumerate(L) if re.match(r'^\[' + name + r'\]\s*(#.*)?$', l))
    except StopIteration:
        return None, None
    i1 = next((i for i in range(i0 + 1, len(L)) if HEAD.match(L[i])), len(L))
    return i0, i1


def drop_counts(text, name=''):
    """Rule 1. Return (new text, message)."""
    L = text.split('\n')
    i0, i1 = _section(L, 'struc')
    if i0 is None:
        return text, ''
    given = {}
    for i in range(i0 + 1, i1):
        m = COUNT.match(L[i])
        if m:
            given[m.group(1)] = (int(m.group(2)), i)
    if not given:
        return text, ''
    heads = {arr: [i for i, l in enumerate(L) if re.match(r'^\[\[' + arr + r'\]\]', l)] for arr in ('site', 'spec')}
    comment, notes, msg = set(), {}, []
    for key, arr in (('nbas', 'site'), ('nspec', 'spec')):
        if key not in given:
            continue
        n, _ = given[key]
        nt = len(heads[arr])
        if n > nt:
            raise SystemExit(f'{name}: {key} = {n} but only {nt} [[{arr}]] tables; fix the file by hand')
        if n < nt:
            for h in heads[arr][n:]:
                e = next((j for j in range(h + 1, len(L)) if HEAD.match(L[j]) or L[j].startswith('# ===')), len(L))
                while e > h + 1 and not L[e - 1].strip():      # trailing blank lines stay blank
                    e -= 1
                comment.update(range(h, e))
            notes[heads[arr][n]] = (f'# [[{arr}]] {n + 1} and after: commented out instead of [struc] {key} = {n} '
                                    f'(ctrlg_update.py, 2026-09-28)')
            msg.append(f'{key} = {n} < {nt}: [[{arr}]] {n + 1}-{nt} commented out')
        else:
            msg.append(f'{key} = {n} removed')
    drop = {i for _, i in given.values()}
    out = []
    for i, l in enumerate(L):
        if i in drop:
            continue
        if l.strip() == OLDNOTE:                           # the [struc] comment of the generators until 2026-09-28
            out.extend(NEWNOTE)
            continue
        if i in notes:
            out.append(notes[i])
        out.append('# ' + l if i in comment and l.strip() else l)
    return '\n'.join(out), '; '.join(msg)


def gw_tetrakbt(text, name=''):
    """Rule 2. Return (new text, message)."""
    L = text.split('\n')
    i0, i1 = _section(L, 'gw')
    if i0 is None:
        return text, ''
    val = lambda i: float(L[i].split('=', 1)[1].split('#')[0].strip())
    isx = [i for i in range(i0 + 1, i1) if re.match(r'^\s*SmearX0q?0?\s*=', L[i])]
    itk = [i for i in range(i0 + 1, i1) if re.match(r'^\s*t_tetrakbt\s*=', L[i])]
    msg = []
    sx = 0.0
    for i in isx:
        if re.match(r'^\s*SmearX0\s*=', L[i]):
            sx = val(i)
    if sx > 0:
        tk = val(itk[0]) if itk else 0.0
        if tk > 0:
            raise SystemExit(f'{name}: [gw] SmearX0 > 0 and t_tetrakbt > 0 together; fix the file by hand')
        T = smearx0_to_kelvin(sx)
        new = (f't_tetrakbt = {-T}   # (K) < 0: Im chi0 smeared by a Gaussian as wide as Fermi-Dirac at |T| '
               f'(was SmearX0 = {sx} Ha; ctrlg_update.py, 2026-09-28)')
        L[isx[0]] = new
        drop = set(isx[1:]) | set(itk)
        msg.append(f'SmearX0 = {sx} -> t_tetrakbt = {-T}')
    else:
        drop = set(isx)
        if isx:
            msg.append('SmearX0 = 0 removed')
        if not itk:
            new = 't_tetrakbt = 0   # (K) 0: no smearing of chi0 (the former default; ctrlg_update.py, 2026-09-28)'
            its = [i for i in range(i0 + 1, i1) if re.match(r'^\s*t_sigmaw\s*=', L[i])]
            at = its[0] + 1 if its else i0 + 1
            L.insert(at, new)
            drop = {j + 1 if j >= at else j for j in drop}
            msg.append('t_tetrakbt = 0 added')
    out = [l for i, l in enumerate(L) if i not in drop]
    return '\n'.join(out), '; '.join(msg)


def update(text, name=''):
    """All rules. Return (new text, message)."""
    msgs = []
    for f in (drop_counts, gw_tetrakbt):
        text, m = f(text, name)
        if m:
            msgs.append(m)
    return text, '; '.join(msgs) if msgs else 'nothing to do'
