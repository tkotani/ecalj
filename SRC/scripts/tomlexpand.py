#!/usr/bin/env python3
"""
tomlexpand.py - Convert GWinput.toml back to legacy GWinput format

Usage: tomlexpand.py [GWinput.toml] [-o GWinput]

Reads the TOML produced by gwinput2toml.py and emits a legacy-format
GWinput. Round-trip with gwinput2toml.py should yield an equivalent
TOML structure (idempotent at the TOML level).
"""

from __future__ import annotations
import sys
import argparse
from pathlib import Path

try:
    import tomllib  # Python 3.11+
except ImportError:
    import tomli as tomllib  # type: ignore

# Keys that are booleans (on/off in legacy)
BOOL_KEYS = {
    "GaussSmear", "KeepEigen", "KeepPPOVL", "NormChk",
    "unit_2pioa", "CoreOrth", "LFC@Gamma",
}


def fmt_value(v: object, key: str = "") -> str:
    """Format a TOML value as legacy GWinput would have it."""
    if isinstance(v, bool) or key in BOOL_KEYS:
        return "on" if v else "off"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        # Match Fortran-friendly form. Avoid 'e-05' for small numbers? Use general %g
        if abs(v) < 1e-3 and v != 0:
            return f"{v:.5e}".replace("e", "d")
        return f"{v:.6f}"
    if isinstance(v, list):
        return " ".join(fmt_value(x, key) for x in v)
    return str(v)


def emit_gw_section(gw: dict) -> list[str]:
    """Emit [gw] section as 'key value' lines."""
    out = []
    for k, v in gw.items():
        out.append(f"{k}    {fmt_value(v, k)}")
    return out


def emit_product_basis(pb: dict) -> list[str]:
    """Emit <PRODUCT_BASIS>...</PRODUCT_BASIS>."""
    out = ["<PRODUCT_BASIS>"]
    out.append("  tolerance to remove products due to poor linear-independency")
    if "tolerance" in pb:
        out.append("  " + " ".join(f"{x:.5e}".replace("e", "d") if abs(x) < 1e-3 else f"{x:g}" for x in pb["tolerance"]))
    out.append("  lcutmx(atom) = maximum l-cutoff for the product basis.")
    if "lcutmx" in pb:
        out.append("  " + " ".join(str(x) for x in pb["lcutmx"]))
    out.append("  atom   l  nnvv  nnc")
    if "nlx" in pb:
        for row in pb["nlx"]:
            out.append("  " + " ".join(f"{x:4d}" for x in row))
    out.append("  atom   l    n  occ unocc")
    if "valence" in pb:
        for row in pb["valence"]:
            out.append("  " + " ".join(f"{x:4d}" for x in row))
    out.append("  atom   l    n  occ unocc  ForX0 ForSxc")
    if "core" in pb:
        for row in pb["core"]:
            out.append("  " + " ".join(f"{x:4d}" for x in row))
    out.append("</PRODUCT_BASIS>")
    return out


def emit_block(tag: str, body: str) -> list[str]:
    """Emit <TAG>...</TAG> with raw body."""
    return [f"<{tag}>", body.rstrip("\n"), f"</{tag}>"]


def expand_toml_to_gwinput(parsed: dict) -> str:
    """Render the parsed dict as a legacy-format GWinput."""
    lines = []
    lines.append("!!! Auto-generated from GWinput.toml by tomlexpand.py")
    lines.append("")
    if "gw" in parsed:
        lines.extend(emit_gw_section(parsed["gw"]))
        lines.append("")
    if "product_basis" in parsed and parsed["product_basis"]:
        lines.extend(emit_product_basis(parsed["product_basis"]))
        lines.append("")
    if "blocks" in parsed:
        for tag, body in parsed["blocks"].items():
            lines.extend(emit_block(tag, body))
            lines.append("")
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description="Convert GWinput.toml to legacy GWinput")
    ap.add_argument("input", nargs="?", default="GWinput.toml", help="Input TOML")
    ap.add_argument("-o", "--output", default=None, help="Output GWinput (default: <input>.expanded)")
    args = ap.parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        print(f"Error: {in_path} not found", file=sys.stderr)
        sys.exit(1)

    out_path = Path(args.output) if args.output else in_path.with_suffix(".expanded")

    parsed = tomllib.loads(in_path.read_text())
    text = expand_toml_to_gwinput(parsed)
    out_path.write_text(text)
    print(f"Wrote {out_path} ({len(text)} bytes)", file=sys.stderr)


if __name__ == "__main__":
    main()
