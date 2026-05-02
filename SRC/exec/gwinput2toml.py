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

# Keys that take on/off (boolean)
BOOL_KEYS = {
    "GaussSmear", "KeepEigen", "KeepPPOVL", "NormChk",
    "unit_2pioa", "CoreOrth", "LFC@Gamma",
}

# Keys that are vectors of integers (BZ mesh etc.)
INT_VEC_KEYS = {"n1n2n3", "BZmesh", "multitet"}


def fortran_num_to_python(s: str) -> str:
    """Convert Fortran-style scientific notation to Python: 1d-3 -> 1e-3, 0.10D-05 -> 0.10E-05."""
    return re.sub(r"([0-9.])(d|D)([+-]?[0-9])", r"\1e\3", s)


def strip_comment(line: str) -> str:
    """Remove trailing comment after '!' (preserve leading whitespace before any '!')."""
    idx = line.find("!")
    return line[:idx] if idx >= 0 else line


def parse_value(tokens: list[str], key: str) -> object:
    """Parse a value (or vector) from whitespace-separated tokens. Returns int, float, bool, str, or list."""
    if key in BOOL_KEYS:
        v = tokens[0].lower()
        return v in ("on", "true", "yes", "1", ".true.")

    if not tokens:
        return None

    # Try multi-value: integer vector
    if key in INT_VEC_KEYS:
        try:
            return [int(t) for t in tokens]
        except ValueError:
            pass

    # Single value: try int, then float, then string
    if len(tokens) == 1:
        t = tokens[0]
        try:
            return int(t)
        except ValueError:
            pass
        try:
            return float(fortran_num_to_python(t))
        except ValueError:
            pass
        return t

    # Multi-token: try float vector
    try:
        return [float(fortran_num_to_python(t)) for t in tokens]
    except ValueError:
        pass

    # Fallback: list of strings
    return tokens


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

    # nlx table: "atom l nnvv nnc"
    # Skip header, collect rows of 4 ints
    while idx < len(clean) and not re.match(r"^\d+\s+\d+\s+\d+\s+\d+\s*$", clean[idx]):
        idx += 1
    nlx = []
    while idx < len(clean) and re.match(r"^\d+\s+\d+\s+\d+\s+\d+\s*$", clean[idx]):
        nlx.append([int(t) for t in clean[idx].split()])
        idx += 1
    pb["nlx"] = nlx

    # Valence table: "atom l n occ unocc" (5 ints)
    while idx < len(clean) and not re.match(r"^\d+\s+\d+\s+\d+\s+\d+\s+\d+\s*$", clean[idx]):
        idx += 1
    valence = []
    while idx < len(clean) and re.match(r"^\d+\s+\d+\s+\d+\s+\d+\s+\d+\s*$", clean[idx]):
        valence.append([int(t) for t in clean[idx].split()])
        idx += 1
    pb["valence"] = valence

    # Core table: "atom l n occ unocc forX0 forSxc" (7 ints)
    while idx < len(clean) and not re.match(
        r"^\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s*$", clean[idx]
    ):
        idx += 1
    core = []
    while idx < len(clean) and re.match(
        r"^\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s*$", clean[idx]
    ):
        core.append([int(t) for t in clean[idx].split()])
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
            else:
                # Pass through as raw
                out["blocks"][tag] = "\n".join(block_lines)
            continue

        # Simple key value
        toks = s.split()
        if len(toks) >= 2:
            key = toks[0]
            val = parse_value(toks[1:], key)
            out["gw"][key] = val
        elif len(toks) == 1:
            # Bare keyword (rare, e.g. flag-only)
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
