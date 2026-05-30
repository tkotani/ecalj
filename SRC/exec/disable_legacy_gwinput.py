#!/usr/bin/env python3
"""Comment out legacy GWinput else-branches in migrated callers.

Pattern:
    if (gwinput_loaded) then
       <TOML body>
    else
       <legacy getkeyvalue body>      <-- comment out
    endif

Strategy:
    - Walk lines, track if-stack scoped to "gwinput_loaded" branches.
    - Inside the matching else block, prefix every non-blank line with '!'.
    - Insert a `call rx(...)` as the first effective else line so it aborts
      if somehow reached.

Safe properties:
    - Original lines are preserved as-is, just prefixed with '! ' so
      reverting is one substitution away.
    - Only the immediate `else` of `if (gwinput_loaded) then` is touched;
      nested ifs inside the else are also commented out as part of the
      block, which is what we want.
"""
import re
import sys
from pathlib import Path

IF_GWLOADED   = re.compile(r'^\s*if\s*\(\s*gwinput_loaded\s*\)\s*then\b', re.IGNORECASE)
IF_NOTGWLOADED = re.compile(r'^\s*if\s*\(\s*\.not\.\s*gwinput_loaded\s*\)\s*', re.IGNORECASE)
ELSE_LINE     = re.compile(r'^(\s*)else\s*(?:!.*)?$', re.IGNORECASE)
ELSEIF_LINE   = re.compile(r'^\s*else\s*if\b', re.IGNORECASE)
IF_ANY_THEN   = re.compile(r'^\s*if\s*\(.*\)\s*then\b', re.IGNORECASE)
ENDIF_LINE    = re.compile(r'^\s*end\s*if\b', re.IGNORECASE)

def process(text: str, file_label: str = '') -> str:
    lines = text.split('\n')
    out = []
    i = 0
    n = len(lines)
    nested_depth = 0  # depth of any nested `if-then`s in the else branch
    in_else = False    # are we inside the else of a top-level gwinput_loaded if?
    rx_emitted = False  # have we written the rx() guard at the start of the else?
    indent_for_rx = ''
    file_done = False
    edits = 0

    while i < n:
        line = lines[i]
        if not in_else:
            m = IF_GWLOADED.match(line)
            if m:
                # Look ahead to find the matching `else` and `endif`.
                # We scan forward, tracking nested if-then.
                depth = 1
                j = i + 1
                else_idx = -1
                endif_idx = -1
                while j < n and depth > 0:
                    lj = lines[j]
                    if depth == 1 and ELSE_LINE.match(lj):
                        else_idx = j
                    if IF_ANY_THEN.match(lj) and not IF_GWLOADED.match(lj):
                        depth += 1
                    elif ENDIF_LINE.match(lj):
                        depth -= 1
                        if depth == 0:
                            endif_idx = j
                            break
                    j += 1
                if else_idx == -1 or endif_idx == -1:
                    # No else, nothing to do
                    out.append(line)
                    i += 1
                    continue
                # Emit the if-then header and the THEN body verbatim
                out.append(line)
                # then body: lines[i+1 .. else_idx-1] verbatim
                for k in range(i+1, else_idx):
                    out.append(lines[k])
                # else line: keep as-is so syntax is preserved
                out.append(lines[else_idx])
                # else body: lines[else_idx+1 .. endif_idx-1]
                # Insert a rx() guard, then comment out the rest.
                indent = re.match(r'^(\s*)', lines[else_idx]).group(1) + '   '
                out.append(f"{indent}call rx('m_GWinput: legacy GWinput reader is disabled. GWinput.toml is required.')")
                for k in range(else_idx+1, endif_idx):
                    body_line = lines[k]
                    if body_line.strip() == '':
                        out.append(body_line)
                    else:
                        out.append('!' + body_line)
                # endif line
                out.append(lines[endif_idx])
                edits += 1
                i = endif_idx + 1
                continue
        out.append(line)
        i += 1

    new_text = '\n'.join(out)
    return new_text, edits


def main():
    if len(sys.argv) < 2:
        print('usage: disable_legacy_gwinput.py FILE [FILE ...]')
        sys.exit(1)
    total = 0
    for fname in sys.argv[1:]:
        p = Path(fname)
        text = p.read_text()
        new_text, edits = process(text, file_label=str(p))
        if edits == 0:
            print(f'{fname}: no `if (gwinput_loaded) then ... else ... endif` patterns found')
            continue
        p.write_text(new_text)
        total += edits
        print(f'{fname}: {edits} block(s) commented out')
    print(f'Total: {total} block(s)')


if __name__ == '__main__':
    main()
