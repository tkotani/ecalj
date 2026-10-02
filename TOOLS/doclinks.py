#!/usr/bin/env python3
"""doclinks.py: check that the documents of ecalj can be followed (2026-10-02, user: 文書を矛盾なく辿れる構造か).

    python3 TOOLS/doclinks.py            (from the top of ecalj)

1. relative links [..](path) and ![..](path) in the tracked *.md (not *_work, not trash): the target must exist.
   In ecaljdoc/ (the site, not ecaljdoc/MD) links without an extension are VitePress pages: path.md must exist.
2. links in the site pages that leave ecaljdoc/ with ../ (they break on the site; use a GitHub URL)
3. Samples/*/README.md without a link to the main documents (ecaljdoc/ or the site URL)
4. files of ecaljdoc/MD not linked or named from CLAUDE.md, ecaljclaude.md, handover.md or README.md (MD)
5. links from the site pages (ecaljdoc/ but not ecaljdoc/MD) to ecaljdoc/MD: not allowed (MD is not on the site; a path in code
   quotes is fine). The other way, from MD to the site pages, is fine (user 2026-10-02)
"""
import os, re, subprocess, sys
from urllib.parse import unquote
top = subprocess.run(['git', 'rev-parse', '--show-toplevel'], capture_output=True, text=True).stdout.strip()
os.chdir(top)
files = [f for f in subprocess.run(['git', 'ls-files', '*.md'], capture_output=True, text=True).stdout.split('\n')
         if f and '_work/' not in f and not f.startswith('trash/')]
rx = re.compile(r'!?\[[^\]]*\]\(([^)\s]+)(?:\s+"[^"]*")?\)')
broken, leave, nolink, tomd = [], [], [], []
for f in files:
    d = os.path.dirname(f); site = f.startswith('ecaljdoc/') and not f.startswith('ecaljdoc/MD/')
    fence = False
    for n, l in enumerate(open(f, encoding='utf-8', errors='replace'), 1):
        if re.match(r'\s*(```|~~~)', l): fence = not fence; continue
        if fence: continue
        t = re.sub(r'`[^`]*`', '', l); t = re.sub(r'\$\$.*?\$\$', '', t); t = re.sub(r'\$[^$]*\$', '', t)   # code and math
        for m in rx.finditer(t):
            u = m.group(1)
            if site and re.search(r'tkotani/ecalj[^/]*/(tree|blob)/[^/]+/(ecaljdoc/)?MD/', u): tomd.append(f'{f}:{n}: {u}')
            if re.match(r'^[a-z]+:', u) or u.startswith('#') or u.startswith('mailto'): continue
            p = unquote(u.split('#')[0])
            if not p: continue
            if p.startswith('/'):                       # site-absolute (VitePress base) -> ecaljdoc/
                q = os.path.normpath(os.path.join('ecaljdoc', p.lstrip('/').replace('ecaljdoc/', '', 1)))
            else:
                q = os.path.normpath(os.path.join(d, p))
            if site and (q.startswith('ecaljdoc/MD') or '/ecaljdoc/MD/' in u or re.search(r'tkotani/ecalj[^/]*/(tree|blob)/[^/]+/(ecaljdoc/)?MD/', u)):
                tomd.append(f'{f}:{n}: {u}')
            if site and not os.path.normpath(q).startswith('ecaljdoc'):
                leave.append(f'{f}:{n}: {u}')
            cands = [q]
            if p.startswith('/'): cands.append(os.path.normpath(os.path.join('ecaljdoc', 'public', p.lstrip('/'))))   # assets of public/
            if site and not os.path.splitext(q)[1]: cands += [q + '.md', os.path.join(q, 'index.md')]
            if site and q.endswith('.html'): cands += [q[:-5] + '.md']
            if not any(os.path.exists(c) for c in cands): broken.append(f'{f}:{n}: {u}')
for r in sorted(set(os.path.dirname(f) for f in files if f.startswith('Samples/') and f.endswith('README.md'))):
    s = open(os.path.join(r, 'README.md'), encoding='utf-8', errors='replace').read()
    if 'ecaljdoc' not in s: nolink.append(os.path.join(r, 'README.md'))
entry = ''.join(open(p, encoding='utf-8').read() for p in ['CLAUDE.md', 'ecaljdoc/MD/ecaljclaude.md', 'ecaljdoc/MD/handover.md'] if os.path.exists(p))
orphan = [m for m in sorted(os.listdir('ecaljdoc/MD')) if m not in entry]
print(f'1. broken relative links: {len(broken)}'); print('\n'.join('   ' + b for b in broken[:60]))
print(f'2. site links leaving ecaljdoc/ with ../: {len(leave)}'); print('\n'.join('   ' + b for b in leave[:30]))
print(f'3. Samples README without a link to the main documents: {len(nolink)}'); print('\n'.join('   ' + b for b in nolink))
print(f'4. ecaljdoc/MD files not named from the entry documents: {orphan}')
print(f'5. links from the site pages to ecaljdoc/MD: {len(tomd)}'); print('\n'.join('   ' + b for b in tomd))
