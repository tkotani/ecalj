#!/usr/bin/env python3
"""Remove nbas / nspec from [struc] of ctrlg.<sname>.toml (2026-09-28: the TOML loader of lmf stops on them).

  ctrlg_drop_counts.py [--no-backup] <ctrlg file> [...]      (in place; the original is kept as <file>.bak_counts)

The numbers of sites and species are those of the [[site]] / [[spec]] tables. A count written beside the tables used
to win over them, so sites added to the tables were silently left out (LiTi2O4 empty spheres, 2026-09-23).
  - a count equal to the tables is removed;
  - a smaller one (the old way to use only the first N sites or species) is replaced by commenting out the blocks
    beyond N, so the calculation stays the same;
  - a larger one is an error in the file: reported, the file is left as it is.
To leave sites out from now on, comment out their [[site]] blocks.
"""
import re
import shutil
import sys

HEAD = re.compile(r'^\[')                    # any table header: [struc], [[site]], [bz], ...
COUNT = re.compile(r'^\s*(nbas|nspec)\s*=\s*(\d+)')
OLDNOTE = '# nbas/nspec are filled automatically from the [[site]]/[[spec]] arrays.'
NEWNOTE = ['# The numbers of sites and species are those of the [[site]]/[[spec]] tables: do not write',
           '# nbas/nspec (lmf stops on them). To leave a site out, comment out its [[site]] block.']


def drop_counts(text, name=''):
    """Return (new text, message)."""
    L = text.split('\n')
    try:
        i0 = next(i for i, l in enumerate(L) if re.match(r'^\[struc\]\s*(#.*)?$', l))
    except StopIteration:
        return text, 'no [struc]'
    i1 = next((i for i in range(i0 + 1, len(L)) if HEAD.match(L[i])), len(L))
    given = {}
    for i in range(i0 + 1, i1):
        m = COUNT.match(L[i])
        if m:
            given[m.group(1)] = (int(m.group(2)), i)
    if not given:
        return text, 'nothing to do'
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
                                    f'(ctrlg_drop_counts.py, 2026-09-28)')
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


if __name__ == '__main__':
    args = sys.argv[1:]
    backup = '--no-backup' not in args
    files = [a for a in args if a != '--no-backup']
    if not files:
        sys.exit(__doc__)
    for f in files:
        text = open(f).read()
        new, msg = drop_counts(text, f)
        if new != text:
            if backup:
                shutil.copy(f, f + '.bak_counts')
            open(f, 'w').write(new)
        print(f'{f}: {msg}')
