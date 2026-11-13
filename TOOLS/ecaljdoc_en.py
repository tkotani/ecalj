#!/usr/bin/env python3
"""English mirror of ecaljdoc (2026-11-13).

ecaljdoc keeps every page twice: ecaljdoc/X.md (the reference, Japanese or English) and ecaljdoc/en/X.md (English).
Pages that are already in English are copied; Japanese pages are translated. Each en/X.md starts with a line that names
the source it was made from (the commit, and the blob hash of X.md), so a later change of X.md can be found.

    python3 TOOLS/ecaljdoc_en.py status          # list pages whose English version is missing or older than X.md
    python3 TOOLS/ecaljdoc_en.py copy [X.md ...] # write en/X.md for missing pages (or the named ones) as a copy of X.md

`copy` rewrites links so that the copy works from en/: links to other pages stay relative (they go to the English
pages, which have the same names), links to files that are not pages (figures, pdf, scripts) are pointed back to the
originals, and site-absolute page links /x get the /en prefix. Text inside ``` fences is left as it is.
The Japanese text of a copied page still has to be translated (by hand or by Claude); `copy` does not translate.
"""
import os
import re
import subprocess
import sys

DOC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ecaljdoc')
DOC = os.path.normpath(DOC)
SKIP_DIRS = {'node_modules', 'trash', 'en', '.vitepress', 'public'}
SKIP_FILES = {'index.md', 'README.md'}   # en/index.md is written by hand; README.md is the repository readme, not a page
HEADER = ('> English version of the page [{name}]({rel}). The pages at the top of ecaljdoc are the reference; '
          'this one was made from {name}.md of {commit} (source blob {blob}).')
PAGE_EXT = ('', '.md', '.html')


def pages():
    out = []
    for root, dirs, files in os.walk(DOC):
        rel = os.path.relpath(root, DOC)
        dirs[:] = [d for d in dirs if not (rel == '.' and d in SKIP_DIRS) and not d.startswith('.')]
        for f in files:
            if f.endswith('.md') and not (rel == '.' and f in SKIP_FILES):
                out.append(os.path.normpath(os.path.join(rel, f)))
    return sorted(out)


def git(*args):
    return subprocess.run(['git', '-C', DOC, *args], capture_output=True, text=True).stdout.strip()


def last_commit(page):
    # Some pages have no history under their present path (it was cut when ecaljdoc was moved into ecalj on
    # 2026-10-02), so the source is identified by its blob hash; the commit is given when git knows it.
    c = git('log', '-1', '--format=%h', '--', page)
    return (('commit ' + c) if c else 'the tree of ' + git('rev-parse', '--short', 'HEAD')), blob(page)


def blob(page):
    return git('hash-object', page)[:12]


def is_page_link(path):
    return os.path.splitext(path)[1] in PAGE_EXT


def rewrite_target(t, src_dir, en_dir):
    if re.match(r'^(https?:|mailto:|#|data:)', t) or not t:
        return t
    path, sep, frag = t.partition('#')
    if path.startswith('/'):
        return ('/en' + path if is_page_link(path) and not path.startswith('/en/') else path) + sep + frag
    if is_page_link(path):
        return t                                   # the English page with the same name
    target = os.path.normpath(os.path.join(src_dir, path))
    return os.path.relpath(target, en_dir) + sep + frag


def convert(page, text, commit, blobh):
    src_dir = os.path.join(DOC, os.path.dirname(page))
    en_dir = os.path.join(DOC, 'en', os.path.dirname(page))
    fix = lambda m: m.group(1) + rewrite_target(m.group(2), src_dir, en_dir) + m.group(3)
    lines, fence = [], False
    for ln in text.split('\n'):
        if ln.lstrip().startswith('```'):
            fence = not fence
        if not fence:
            ln = re.sub(r'(\]\()([^)\s]*)(\))', fix, ln)
            ln = re.sub(r'((?:src|href)=")([^"]*)(")', fix, ln)
        lines.append(ln)
    body = '\n'.join(lines)
    rel_src = os.path.relpath(os.path.join(DOC, page), en_dir)
    head = HEADER.format(name=os.path.splitext(page)[0], rel=rel_src, commit=commit, blob=blobh)
    m = re.match(r'(---\n.*?\n---\n)', body, re.S)           # keep a YAML front matter first
    front, rest = (m.group(1), body[m.end():]) if m else ('', body)
    m = re.match(r'(\s*# [^\n]*\n)', rest)                  # header line just under the title
    if m:
        return front + m.group(1) + '\n' + head + '\n' + rest[m.end():]
    return front + head + '\n\n' + rest


def status():
    for p in pages():
        en = os.path.join(DOC, 'en', p)
        if not os.path.exists(en):
            print(f'missing  {p}')
            continue
        m = re.search(r'source blob ([0-9a-f]{7,})', open(en).read(3000))
        if not m:
            print(f'no-blob  {p}')
            continue
        if not blob(p).startswith(m.group(1)):
            print(f'stale    {p}  ({p} has changed since the English version was made)')


def copy(names):
    todo = names or [p for p in pages() if not os.path.exists(os.path.join(DOC, 'en', p))]
    for p in todo:
        commit, blobh = last_commit(p)
        out = os.path.join(DOC, 'en', p)
        os.makedirs(os.path.dirname(out), exist_ok=True)
        open(out, 'w').write(convert(p, open(os.path.join(DOC, p)).read(), commit, blobh))
        print(f'wrote en/{p}  (from {commit}, blob {blobh})')


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'status'
    if cmd == 'status':
        status()
    elif cmd == 'copy':
        copy(sys.argv[2:])
    else:
        sys.exit(__doc__)
