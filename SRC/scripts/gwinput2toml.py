#!/usr/bin/env python3
"""
gwinput2toml.py - Convert legacy GWinput to GWinput.toml

Usage: gwinput2toml.py [GWinput] [-o GWinput.toml]

Reads ecalj's tag-based GWinput format and emits a TOML file with the
same content, structured under [gw] (scalar/vector keys) and
[product_basis] (PRODUCT_BASIS block). Other blocks (<QPNT>, <QforGW>,
<QforEPS>, <QforEPSL>, <Worb>) are passed through as raw multi-line
strings under [blocks] for future structured handling.
"""

from __future__ import annotations
import re
import sys
import argparse
from pathlib import Path

# Keys that take on/off (boolean) or .true./.false.
BOOL_KEYS = {
    "GaussSmear", "KeepEigen", "KeepPPOVL", "NormChk",
    "unit_2pioa", "CoreOrth", "LFC@Gamma",
    "AnyQ", "QforEPSau", "QforEPSunita", "QforEPSLIncLeft",
    "tetrakbt", "wan_in_ewin",
}

# Keys that are vectors of integers (BZ mesh etc.)
INT_VEC_KEYS = {"n1n2n3", "BZmesh", "multitet", "n1n2n3eps"}


def fortran_num_to_python(s: str) -> str:
    """Convert Fortran-style scientific notation to Python: 1d-3 -> 1e-3, 0.10D-05 -> 0.10E-05."""
    return re.sub(r"([0-9.])(d|D)([+-]?[0-9])", r"\1e\3", s)


def strip_comment(line: str) -> str:
    """Remove trailing comment after '!' or '#' (whichever comes first).

    GWinput uses '!' as the documented comment marker but some files also
    use '#' inline (e.g. '2 2 2 2  #4 4 3 3' on the lcutmx line). Strip both.
    """
    cuts = [i for i in (line.find("!"), line.find("#")) if i >= 0]
    if not cuts:
        return line
    return line[:min(cuts)]


def _try_int(t: str):
    try: return int(t)
    except ValueError: return None

def _try_float(t: str):
    try: return float(fortran_num_to_python(t))
    except ValueError: return None

def parse_value(tokens: list[str], key: str) -> object:
    """Parse value(s) from whitespace tokens. Strips trailing non-numeric
    annotations (e.g., 'HistBin_dw 0.01 (a.u.)' -> 0.01). Returns int,
    float, bool, str, or homogeneous numeric list."""
    if key in BOOL_KEYS:
        v = tokens[0].lower()
        return v in ("on", "true", "yes", "1", ".true.")

    if not tokens:
        return None

    # Take only the leading run of numeric tokens. Drop trailing annotations.
    numeric_run = []
    for t in tokens:
        if _try_float(t) is not None:
            numeric_run.append(t)
        else:
            break

    if not numeric_run:
        # Pure string scalar
        return tokens[0]

    if key in INT_VEC_KEYS:
        ints = [_try_int(t) for t in numeric_run]
        if all(v is not None for v in ints):
            return ints

    if len(numeric_run) == 1:
        t = numeric_run[0]
        iv = _try_int(t)
        if iv is not None and "." not in t and "d" not in t.lower() and "e" not in t.lower():
            return iv
        return _try_float(t)

    # multi-numeric: prefer ints if all integral
    ints = [_try_int(t) for t in numeric_run]
    if all(v is not None for v in ints):
        return ints
    return [_try_float(t) for t in numeric_run]


def parse_product_basis(block_lines: list[str]) -> dict:
    """Parse <PRODUCT_BASIS> ... </PRODUCT_BASIS> content into structured dict."""
    # Strip comments and empty lines, preserve order
    clean = []
    for line in block_lines:
        s = strip_comment(line).strip()
        if s:
            clean.append(s)

    pb = {}
    idx = 0

    # First descriptive line (skip text), then tolerance value(s)
    # Heuristic: skip lines until we find one that's purely numeric
    while idx < len(clean) and not re.match(r"^[\d\s.deDE+-]+$", clean[idx]):
        idx += 1
    if idx < len(clean):
        toks = clean[idx].split()
        pb["tolerance"] = [float(fortran_num_to_python(t)) for t in toks]
        idx += 1

    # Skip until next numeric line: lcutmx values
    while idx < len(clean) and not re.match(r"^[\d\s+-]+$", clean[idx]):
        idx += 1
    if idx < len(clean):
        pb["lcutmx"] = [int(t) for t in clean[idx].split()]
        idx += 1

    def line_starts_with_n_ints(line: str, n: int) -> bool:
        """Return True iff line begins with at least n whitespace-separated ints."""
        toks = line.split()
        if len(toks) < n:
            return False
        for t in toks[:n]:
            try:
                int(t)
            except ValueError:
                return False
        return True

    # nlx table: "atom l nnvv nnc" — at least 4 ints (some files have trailing fields)
    while idx < len(clean) and not line_starts_with_n_ints(clean[idx], 4):
        idx += 1
    nlx = []
    while idx < len(clean) and line_starts_with_n_ints(clean[idx], 4):
        nlx.append([int(t) for t in clean[idx].split()[:4]])
        idx += 1
    pb["nlx"] = nlx

    # Valence table: "atom l n occ unocc" — at least 5 ints
    while idx < len(clean) and not line_starts_with_n_ints(clean[idx], 5):
        idx += 1
    valence = []
    while idx < len(clean) and line_starts_with_n_ints(clean[idx], 5):
        valence.append([int(t) for t in clean[idx].split()[:5]])
        idx += 1
    pb["valence"] = valence

    # Core table: "atom l n occ unocc forX0 forSxc" — at least 7 ints
    # (some files have extra fields like "1 0 3  0 0 0 0   1 0  1 1")
    while idx < len(clean) and not line_starts_with_n_ints(clean[idx], 7):
        idx += 1
    core = []
    while idx < len(clean) and line_starts_with_n_ints(clean[idx], 7):
        core.append([int(t) for t in clean[idx].split()[:7]])
        idx += 1
    pb["core"] = core

    return pb


def parse_gwinput(text: str) -> dict:
    """Parse a GWinput file content into a structured dict."""
    out = {"gw": {}, "product_basis": {}, "blocks": {}}
    lines = text.splitlines()
    i = 0

    while i < len(lines):
        raw = lines[i]
        s = strip_comment(raw).strip()

        if not s:
            i += 1
            continue

        # Block start: <TAG>
        m = re.match(r"^<([A-Za-z0-9_@]+)\s*>", s)
        if m and not s.startswith("</"):
            tag = m.group(1)
            block_lines = []
            i += 1
            # Collect until </TAG>
            while i < len(lines):
                end = re.match(r"^\s*</" + tag + r"\s*>", lines[i])
                if end:
                    i += 1
                    break
                block_lines.append(lines[i])
                i += 1

            if tag == "PRODUCT_BASIS":
                out["product_basis"] = parse_product_basis(block_lines)
                # DEBUG: also keep raw for legacy parser fallback
                if tag not in out["blocks"]:
                    out["blocks"][tag] = "\n".join(block_lines)
            else:
                # Keep the FIRST occurrence of duplicate tags (e.g. <Worb>):
                # legacy GWinput files often have a real first block plus a
                # commented-out example block at the end as documentation.
                # Last-write-wins would silently drop the real data.
                if tag not in out["blocks"]:
                    out["blocks"][tag] = "\n".join(block_lines)
            continue

        # Simple key value. Some files use 'key=value' form with no space
        # (e.g. wan_conv_1st=1d-7). Normalize: replace first '=' with space.
        if "=" in s and not s.split()[0].replace("=", "").replace(".", "").replace("-", "").isdigit():
            # Heuristic: if first token has '=', split it
            first = s.split()[0]
            if "=" in first:
                eq = first.index("=")
                rest = s[s.index(first) + len(first):]
                s = first[:eq] + " " + first[eq+1:] + rest
        toks = s.split()
        if len(toks) >= 2:
            key = toks[0]
            val = parse_value(toks[1:], key)
            out["gw"][key] = val
        elif len(toks) == 1:
            out["gw"][toks[0]] = True

        i += 1

    return out


def emit_value(v: object) -> str:
    """Emit a single TOML value (no newline)."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        # Python's repr usually gives a clean form
        return repr(v)
    if isinstance(v, list):
        if not v:
            return "[]"
        return "[" + ", ".join(emit_value(x) for x in v) + "]"
    # String: TOML basic string with escape
    s = str(v).replace("\\", "\\\\").replace('"', '\\"')
    return '"' + s + '"'


def emit_toml(parsed: dict) -> str:
    """Render the parsed dict as TOML text."""
    out = []

    out.append("# GWinput.toml — auto-generated from legacy GWinput by gwinput2toml.py")
    out.append("")

    if parsed["gw"]:
        out.append("[gw]")
        for k, v in parsed["gw"].items():
            out.append(f"{k} = {emit_value(v)}")
        out.append("")

    pb = parsed["product_basis"]
    if pb:
        out.append("[product_basis]")
        if "tolerance" in pb:
            out.append(f"tolerance = {emit_value(pb['tolerance'])}")
        if "lcutmx" in pb:
            out.append(f"lcutmx = {emit_value(pb['lcutmx'])}")
        out.append("")
        if "nlx" in pb:
            out.append("# [iatom, l, nnvv, nnc]")
            out.append("nlx = [")
            for row in pb["nlx"]:
                out.append(f"  {emit_value(row)},")
            out.append("]")
            out.append("")
        if "valence" in pb:
            out.append("# [iatom, l, n, occ, unocc]")
            out.append("valence = [")
            for row in pb["valence"]:
                out.append(f"  {emit_value(row)},")
            out.append("]")
            out.append("")
        if "core" in pb:
            out.append("# [iatom, l, n, occ, unocc, forX0, forSxc]")
            out.append("core = [")
            for row in pb["core"]:
                out.append(f"  {emit_value(row)},")
            out.append("]")
            out.append("")

    if parsed["blocks"]:
        out.append("# Blocks not yet structured — kept as raw text for round-trip")
        out.append("[blocks]")
        for tag, body in parsed["blocks"].items():
            out.append(f'{tag} = """')
            out.append(body)
            out.append('"""')
            out.append("")

    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description="Convert GWinput to GWinput.toml")
    ap.add_argument("input", nargs="?", default="GWinput", help="Input GWinput file (default: GWinput)")
    ap.add_argument("-o", "--output", default=None, help="Output TOML (default: <input>.toml)")
    args = ap.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output) if args.output else in_path.with_suffix(in_path.suffix + ".toml")

    if not in_path.exists():
        print(f"Error: input file '{in_path}' not found", file=sys.stderr)
        sys.exit(1)

    text = in_path.read_text()
    parsed = parse_gwinput(text)
    toml_text = emit_toml(parsed)
    out_path.write_text(toml_text)

    print(f"Wrote {out_path} ({len(toml_text)} bytes)", file=sys.stderr)


if __name__ == "__main__":
    main()
