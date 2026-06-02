#!/usr/bin/env python3
"""Migrate existing ctrlG.<sname>.toml to the new schema:

  - drop the [io] section header; its two keys (verbose, time) move
    to the top of the file (just after symgrp / symgrpaf)
  - rename truncated spellings `verbos` -> `verbose`, `tim` -> `time`
    inside what used to be the [io] section
  - inject `phispinsym = false` into [ham] if missing
    (mandatory in the new schema; old files may not declare it)

Idempotent: a file that already follows the new schema is left alone.

Usage:
  ./migrate_ctrlG_v2.py [ctrlG.*.toml ...]
"""
import sys, re, pathlib

def split_sections(text):
    """Return a list of (header_or_None, body_text) tuples in file order.
    A `header_or_None` is the literal `[section]` / `[[section]]` line; None
    for the pre-section preamble at file start.
    """
    sections = []
    cur_head = None
    cur_lines = []
    for line in text.splitlines(keepends=True):
        m = re.match(r'\s*(\[\[[^\]]+\]\]|\[[^\]]+\])\s*$', line.rstrip('\n'))
        if m:
            sections.append((cur_head, ''.join(cur_lines)))
            cur_head = line
            cur_lines = []
        else:
            cur_lines.append(line)
    sections.append((cur_head, ''.join(cur_lines)))
    return sections

def migrate(text):
    sections = split_sections(text)
    # 1) Find [io] section, extract verbose/time lines, drop the header.
    io_lines = []
    new_sections = []
    for head, body in sections:
        if head and head.strip() == '[io]':
            # Rename legacy spellings while we are moving the lines up.
            body = re.sub(r'^(\s*)verbos(\s*=)', r'\1verbose\2', body, flags=re.M)
            body = re.sub(r'^(\s*)tim(\s*=)',    r'\1time\2',    body, flags=re.M)
            io_lines.append(body)
            continue
        new_sections.append((head, body))

    # 2) Splice the io lines into the top-of-file preamble (first None section).
    if io_lines:
        merged_io = ''.join(io_lines).rstrip() + '\n'
        # Find the preamble (first entry, head=None) and append io_lines to it.
        head0, body0 = new_sections[0]
        if not body0.endswith('\n'):
            body0 += '\n'
        new_sections[0] = (head0, body0 + merged_io)

    # 3) Inject phispinsym = false into [ham] if not present.
    out = []
    for head, body in new_sections:
        if head and head.strip() == '[ham]':
            if 'phispinsym' not in body:
                # Insert just after the header (preserving any leading
                # comments) on a fresh line.
                body = 'phispinsym = false  # spin-averaged radial wfns (needed for SO=1 perturbation)\n' + body
        out.append(head or '')
        out.append(body)
    return ''.join(out)

def main():
    paths = [pathlib.Path(p) for p in sys.argv[1:]]
    if not paths:
        print(__doc__)
        sys.exit(1)
    changed = 0
    for p in paths:
        if not p.exists():
            print(f'skip (missing): {p}', file=sys.stderr)
            continue
        old = p.read_text()
        new = migrate(old)
        if new != old:
            p.with_suffix(p.suffix + '.v1bak').write_text(old)
            p.write_text(new)
            print(f'migrated: {p}')
            changed += 1
        else:
            print(f'unchanged: {p}')
    print(f'{changed}/{len(paths)} files updated')

if __name__ == '__main__':
    main()
