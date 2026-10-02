#!/usr/bin/env python3
"""mdlinkify.py <md>...: turn a code-quoted repository path (`Samples/Magnon/Fe_mlo_magnon/README.md`) into a relative link
when the file or directory exists (2026-10-02, user: paths in the MD documents should be clickable). Paths are taken from the top
of the repository (also `ecalj/...` and `<ecalj>/...`), or relative to the document. Fenced code, link texts and paths with
wildcards or placeholders (<sname>, *, ...) are left alone. Prints what it changed; --dry to only print."""
import os, re, sys, subprocess
top = subprocess.run(['git', 'rev-parse', '--show-toplevel'], capture_output=True, text=True).stdout.strip()
dry = '--dry' in sys.argv
rx = re.compile(r'(?<!\[)`([A-Za-z0-9_.][^`\s<>*{}|]*?)`(?![^\[]*\]\()')
for f in [a for a in sys.argv[1:] if a != '--dry']:
    fa = os.path.abspath(f); d = os.path.dirname(fa)
    out, fence, n = [], False, [0]
    for l in open(fa, encoding='utf-8'):
        if re.match(r'\s*(```|~~~)', l): fence = not fence; out.append(l); continue
        if fence: out.append(l); continue
        def sub(m):
            p = m.group(1).rstrip('/')
            if '/' not in p and not p.endswith('.md'): return m.group(0)        # bare names (keys, programs) stay
            cands = [os.path.join(top, re.sub(r'^(ecalj/|<ecalj>/)', '', p)), os.path.join(d, p)]
            for c in cands:
                if os.path.exists(c):
                    rel = os.path.relpath(c, d)
                    if os.path.isdir(c) and os.path.exists(os.path.join(c, 'README.md')): rel = os.path.join(rel, 'README.md')
                    n[0] += 1
                    return f'[`{m.group(1)}`]({rel})'
            return m.group(0)
        # skip text already inside a markdown link: split on links and only touch the rest
        parts = re.split(r'(\[[^\]]*\]\([^)]*\))', l)
        out.append(''.join(pt if pt.startswith('[') and '](' in pt else rx.sub(sub, pt) for pt in parts))
    if n[0]:
        print(f'{f}: {n[0]} paths linked')
        if not dry: open(fa, 'w', encoding='utf-8').write(''.join(out))
