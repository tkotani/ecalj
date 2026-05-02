#!/usr/bin/env python3
"""Systematically check m_GWinput compile-time defaults vs legacy getkeyvalue
default= values. Flags any mismatch.

Usage: python3 check_defaults.py <ecalj-src-dir>
  Default ecalj-src-dir: ~/ecaljdeveloper/SRC
"""
import re
import sys
import os
from pathlib import Path
from collections import defaultdict


# Pattern for m_GWinput declarations:
#   real(8), protected, public :: KEY = VALUE
#   integer, protected, public :: KEY = VALUE
#   logical, protected, public :: KEY = VALUE
DECL = re.compile(
    r'^\s*(real\(8\)|integer|logical)\s*,\s*protected\s*,\s*public\s*::\s*(\w+)(?:\((\d+)\))?\s*=\s*([^!\n]+?)\s*(?:!.*)?$',
    re.MULTILINE,
)

# Pattern for legacy getkeyvalue calls (live or commented):
#   call getkeyvalue("GWinput","KEY",VAR, default=VALUE )
GET = re.compile(
    r'^[^"\n]*?call\s+getkeyvalue\s*\(\s*["\']GWinput["\']\s*,\s*["\']([^"\']+)["\']\s*,[^,)]+(?:\s*,\s*\d+)?\s*(?:,\s*default\s*=\s*([^,)\n]+))?',
    re.IGNORECASE | re.MULTILINE,
)


def parse_m_gwinput(src_dir: Path):
    """Return dict: key -> (type, default_str)."""
    p = src_dir / 'subroutines' / 'm_GWinput.f90'
    text = p.read_text()
    out = {}
    for m in DECL.finditer(text):
        typ, key, _arr, default = m.groups()
        out[key] = (typ.replace(' ', ''), default.strip())
    return out


def parse_legacy_defaults(src_dir: Path):
    """Scan all .f90 files for legacy getkeyvalue() calls with default=.
    Return dict: key -> set of default strings (across all callsites).
    Both live calls and commented-out ones are included.
    """
    out = defaultdict(set)
    for f in (src_dir / 'subroutines').glob('*.f90'):
        text = f.read_text()
        for line in text.split('\n'):
            stripped = line.strip()
            # Strip leading `!` so commented calls are also picked up
            probe = re.sub(r'^!+\s*', '', stripped)
            m = re.search(
                r'call\s+getkeyvalue\s*\(\s*["\']GWinput["\']\s*,\s*["\']([^"\']+)["\']\s*,(.*)$',
                probe,
                re.IGNORECASE,
            )
            if not m:
                continue
            key, rest = m.groups()
            d = re.search(r'default\s*=\s*([^,)\n]+)', rest)
            if not d:
                continue
            out[key].add(d.group(1).strip())
    return out


def normalize(s: str, typ: str) -> str:
    """Normalize a Fortran literal for comparison."""
    s = s.strip().lower().rstrip(' ')
    s = s.rstrip(')')
    s = re.sub(r'\s+', '', s)
    s = s.replace('.true.', 'T').replace('.false.', 'F')
    # numeric forms: 1.0d-5, 1d-5, 1e-5, 1.0e-5 -> 1.0e-05
    if typ in ('real(8)', 'integer'):
        try:
            v = float(s.replace('d', 'e').replace('+', ''))
            return f'{v:.6e}'
        except ValueError:
            pass
    # array literal like (/...,...,.../)
    s = s.replace('(/', '[').replace('/)', ']')
    return s


def main():
    src_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home() / 'ecaljdeveloper' / 'SRC'
    if not src_dir.is_dir():
        print(f'ERROR: {src_dir} is not a directory', file=sys.stderr)
        sys.exit(1)

    defaults_m = parse_m_gwinput(src_dir)
    defaults_legacy = parse_legacy_defaults(src_dir)

    keys = set(defaults_m) | set(defaults_legacy)
    mismatches = []
    only_legacy = []
    only_m = []
    for key in sorted(keys):
        if key not in defaults_m:
            only_legacy.append((key, defaults_legacy[key]))
            continue
        if key not in defaults_legacy:
            # m_GWinput has it but no getkeyvalue caller -- not a regression
            only_m.append((key, defaults_m[key]))
            continue
        typ, m_def = defaults_m[key]
        m_norm = normalize(m_def, typ)
        legacy_set = {normalize(d, typ) for d in defaults_legacy[key]}
        if m_norm not in legacy_set:
            mismatches.append((key, typ, m_def, defaults_legacy[key], m_norm, legacy_set))

    print(f'==== Default mismatches between m_GWinput and legacy getkeyvalue ====')
    if not mismatches:
        print('No mismatches found.')
    else:
        for key, typ, m_def, leg_set, m_norm, leg_norm in mismatches:
            print(f'\n  {key} ({typ})')
            print(f'    m_GWinput default : {m_def!r:20}  -> {m_norm}')
            print(f'    legacy default(s) : {sorted(leg_set)}')
            print(f'                        raw: {sorted(legacy_set:=leg_set)}')

    print(f'\n==== Keys present in legacy callers but NOT in m_GWinput ====')
    for key, defs in only_legacy:
        print(f'  {key}: {sorted(defs)}')

    print(f'\n==== Keys in m_GWinput but no getkeyvalue caller ====')
    for key, (typ, val) in only_m:
        print(f'  {key} ({typ}) = {val}')


if __name__ == '__main__':
    main()
