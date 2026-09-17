"""Tidy the GW-side sections of ctrlg.<sname>.toml.

Only [gw], [mlo], [product_basis] and [blocks] are touched; the ctrl sections
written by ctrlgenToml.py ([struc], [[site]], [[spec]], [bz], [iter], [ham],
[esm]) are left byte-for-byte as they are.

What it does, per section:
  - aligns   key = value   # comment   into columns (like [esm])
  - one blank line before every section header; order [gw] [mlo] [blocks]
    [product_basis] ([product_basis] last, whatever order the file had)
  - one blank line before and after every multi-line \"\"\" block
  - drops header comments that got orphaned by earlier key moves
    ("# === PRODUCT_BASIS ===" outside [product_basis], "# === BLOCKS: ..."
    with nothing under it)
The parsed TOML is unchanged; callers should verify that (Legacy2toml.py does).
"""
from __future__ import annotations
import re

GW_SECTIONS = ("gw", "mlo", "product_basis", "blocks")
_HDR = re.compile(r"^\[(\w+)\]\s*$")
_KEY = re.compile(r"^(\s*)([A-Za-z_][A-Za-z0-9_]*)(\s*=\s*)(.*?)(\s*#.*)?$")
# header comments that belong to [product_basis]; dropped when found elsewhere
_ORPHAN_PREFIX = (
    "# === PRODUCT_BASIS",
    "# pb_tolerance / pb_lcutmx",
    "# nlx / valence / core",
    "# Per-atom product-basis tables",
)
_BLOCKS_HDR = re.compile(r"^# === BLOCKS:")


def _split_sections(text: str):
    """[(name_or_None, [lines])] in file order; name None = preamble."""
    out = []; cur = (None, [])
    for line in text.split("\n"):
        m = _HDR.match(line)
        if m:
            out.append(cur); cur = (m.group(1), [line])
        else:
            cur[1].append(line)
    out.append(cur)
    return out


def _items(lines):
    """Group a section body into items: ('hdr',[line]) ('blank',[]) ('comment',[..])
    ('key',[line]) ('block',[lines])."""
    items = []; i = 0; n = len(lines)
    while i < n:
        L = lines[i]
        if i == 0 and _HDR.match(L):
            items.append(("hdr", [L])); i += 1; continue
        if L.strip() == "":
            items.append(("blank", [])); i += 1; continue
        if L.lstrip().startswith("#"):
            j = i
            while j < n and lines[j].lstrip().startswith("#"):
                j += 1
            items.append(("comment", lines[i:j])); i = j; continue
        m = _KEY.match(L)
        if m and m.group(4).startswith('"""'):
            j = i + 1
            while j < n and '"""' not in lines[j]:
                j += 1
            items.append(("block", lines[i:j + 1])); i = j + 1; continue
        items.append(("key", [L])); i += 1
    return items


def _drop_orphans(items, secname):
    """Drop [product_basis]/[blocks] header comments that sit in the MIDDLE of
    another section (left behind by earlier key moves). A comment run at the
    very end of a section that starts with '# ===' is the next section's
    header and is kept."""
    body = [it for it in items if it[0] != "blank"]
    last_comment = body[-1] if body and body[-1][0] == "comment" and body[-1][1][0].lstrip().startswith("# ===") else None
    out = []
    for it in items:
        kind, ls = it
        if kind == "comment" and it is not last_comment:
            keep = [l for l in ls if not (secname != "product_basis" and l.strip().startswith(_ORPHAN_PREFIX))]
            keep = [l for l in keep if not (secname != "blocks" and _BLOCKS_HDR.match(l.strip()))]
            if secname != "blocks":
                keep = [l for l in keep if not l.strip().startswith("# QforEPS / QforEPSL : q-point lists")]
            if not keep:
                continue
            ls = keep
        out.append((kind, ls))
    return out


def _align_keys(items):
    keys = [ls[0] for k, ls in items if k == "key"]
    parsed = {}
    for L in keys:
        m = _KEY.match(L)
        if m:
            parsed[L] = (m.group(2), m.group(4).rstrip(), (m.group(5) or "").strip())
    if not parsed:
        return items
    kw = min(max(len(p[0]) for p in parsed.values()), 14)
    vw = min(max(len(p[1]) for p in parsed.values()), 26)
    out = []
    for kind, ls in items:
        if kind == "key" and ls[0] in parsed:
            key, val, cmt = parsed[ls[0]]
            line = f"{key:<{kw}} = {val}"
            if cmt:
                line = f"{line:<{kw + 3 + vw}}  {cmt}"
            out.append(("key", [line.rstrip()]))
        else:
            out.append((kind, ls))
    return out


def _emit(items):
    """Re-join with the spacing rule: blank around blocks, no doubled blanks.
    A trailing '# === X ===' comment run is the header of the NEXT section and
    is returned separately so the caller can glue it to that header."""
    lines = []
    prev = None
    body = [it for it in items if it[0] != "blank"]
    trailer = []
    if body and body[-1][0] == "comment" and body[-1][1][0].lstrip().startswith("# ==="):
        trailer = body[-1][1]; body = body[:-1]
    for kind, ls in body:
        if kind == "block":
            key = _KEY.match(ls[0]).group(2)
            attached = prev == "comment" and key in " ".join(lines[-3:])
            if lines and lines[-1] != "" and not attached and prev != "hdr":
                lines.append("")
        if kind == "comment" and prev in ("key", "block") and lines and lines[-1] != "":
            lines.append("")
        lines.extend(ls)
        if kind == "block":
            lines.append("")
        prev = kind
    while lines and lines[-1] == "":
        lines.pop()
    return lines, trailer


ORDER = ("gw", "mlo", "blocks", "product_basis")   # [product_basis] goes last


def tidy_gw_sections(text: str) -> str:
    secs = _split_sections(text)
    # A '# === X ===' comment run at the END of one section is the header of
    # the NEXT section: move it there so sections can be reordered freely.
    moved = []
    for name, lines in secs:
        if moved and moved[-1][0] in GW_SECTIONS + (None,):
            pass
        moved.append([name, lines])
    for i in range(len(moved) - 1):
        _, lines = moved[i]
        k = len(lines)
        while k > 0 and lines[k-1].strip() == "":
            k -= 1
        j = k
        while j > 0 and lines[j-1].lstrip().startswith("#"):
            j -= 1
        run = lines[j:k]
        if run and any(l.lstrip().startswith("# ===") for l in run):
            moved[i][1] = lines[:j]
            moved[i+1][1] = run + moved[i+1][1]        # header run now precedes '[section]'
    # split into the leading non-GW part and the GW-side sections
    head = []; gw = {}
    for name, lines in moved:
        if name in GW_SECTIONS:
            gw[name] = lines
        else:
            head.extend(lines)
    while head and head[-1].strip() == "":
        head.pop()
    out = list(head)
    for name in ORDER:
        if name not in gw:
            continue
        lines = gw[name]
        # header comment run (moved above) comes first, then the items
        h = 0
        while h < len(lines) and lines[h].lstrip().startswith("#"):
            h += 1
        hdr_run, body_lines = lines[:h], lines[h:]
        items = _items(body_lines)
        items = _drop_orphans(items, name)
        items = _align_keys(items)
        body, trailer = _emit(items)       # trailer: stray header run, discard (re-added by rule above)
        out.append("")
        out.extend(hdr_run)
        out.extend(body)
    res = "\n".join(out)
    res = re.sub(r"\n{3,}", "\n\n", res)
    return res.rstrip("\n") + "\n"


if __name__ == "__main__":
    import sys, tomllib
    for p in sys.argv[1:]:
        s = open(p, encoding="utf-8").read()
        t = tidy_gw_sections(s)
        assert tomllib.loads(s) == tomllib.loads(t), f"content changed: {p}"
        open(p, "w", encoding="utf-8").write(t)
        print("tidied", p)
