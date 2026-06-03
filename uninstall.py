#!/usr/bin/env python3
"""Remove everything InstallAll.py placed outside the ecalj/ checkout.

Reads BIN_DIR/ecalj_install_manifest.txt (written by InstallAll.py) and
deletes every listed path: symlinks dropped by InstallAll.py + the
cmake `deliver` target (lmf, lmfa, lmchk, gwsc, hsfp0, libecaljF*.so,
viewvesta, getsyml, ...) plus the real files ecalj_cmdopts.list and
libgemmul8.* (if --gemmul8 was used).

Also strips the auto-installed bash-completion block from ~/.bashrc.

Does NOT delete the ecalj/ source checkout itself nor SRC/build_*; those
are inside the clone and disappear when you `rm -rf` the checkout.
"""

import argparse
import os
import sys
from pathlib import Path

DEFAULT_BIN_DIR = Path.home() / "bin"
MANIFEST_NAME   = "ecalj_install_manifest.txt"
BASHRC_MARKER   = "# >>> ecalj bash completion (auto-installed by InstallAll.py) >>>"
BASHRC_END      = "# <<< ecalj bash completion <<<"


def parse_manifest(manifest_path: Path):
    """Return (paths, header) from the manifest file."""
    if not manifest_path.is_file():
        sys.exit(
            f"ERROR: manifest not found at {manifest_path}\n"
            "       Re-run InstallAll.py once to regenerate it, or pass\n"
            "       --bindir <path> if you installed to a non-default location."
        )
    paths, header = [], []
    for line in manifest_path.read_text().splitlines():
        if line.startswith("#"):
            header.append(line)
        elif line.strip():
            paths.append(Path(line))
    return paths, header


def remove_path(p: Path, dry_run: bool) -> bool:
    """Best-effort remove of a single path. Returns True if (would be) removed."""
    if not (p.is_symlink() or p.exists()):
        return False
    if dry_run:
        print(f"  would rm {p}")
        return True
    try:
        p.unlink()
        print(f"  rm {p}")
        return True
    except IsADirectoryError:
        import shutil
        shutil.rmtree(p, ignore_errors=True)
        print(f"  rmtree {p}")
        return True
    except OSError as e:
        print(f"  WARN: failed to remove {p}: {e}", file=sys.stderr)
        return False


def strip_bashrc_block(bashrc: Path, dry_run: bool) -> bool:
    """Remove the [MARKER ... END] block (inclusive) from bashrc. Returns True if changed."""
    if not bashrc.is_file():
        return False
    text = bashrc.read_text()
    if BASHRC_MARKER not in text:
        return False
    out, skipping = [], False
    for line in text.splitlines(keepends=True):
        if not skipping and BASHRC_MARKER in line:
            skipping = True
            continue
        if skipping:
            if BASHRC_END in line:
                skipping = False
            continue
        out.append(line)
    new_text = "".join(out).rstrip() + "\n"
    if dry_run:
        print(f"  would strip ecalj bash-completion block from {bashrc}")
        return True
    bashrc.write_text(new_text)
    print(f"  stripped ecalj bash-completion block from {bashrc}")
    return True


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Undo InstallAll.py: remove every BIN_DIR entry listed in "
            "the install manifest, and strip the bash-completion block "
            "from ~/.bashrc.  Leaves the ecalj source checkout intact."
        )
    )
    ap.add_argument(
        "--bindir",
        type=Path,
        default=DEFAULT_BIN_DIR,
        help=f"directory InstallAll.py was told to populate (default: {DEFAULT_BIN_DIR})",
    )
    ap.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help=f"path to the install manifest (default: <bindir>/{MANIFEST_NAME})",
    )
    ap.add_argument(
        "--no-bashrc",
        action="store_true",
        help="leave ~/.bashrc untouched",
    )
    ap.add_argument(
        "-y", "--yes",
        action="store_true",
        help="skip the y/N confirmation",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="print what would be removed without touching anything",
    )
    args = ap.parse_args()

    bin_dir = args.bindir.expanduser().resolve()
    manifest = (args.manifest or (bin_dir / MANIFEST_NAME)).expanduser().resolve()
    paths, header = parse_manifest(manifest)

    print(f"Install manifest : {manifest}")
    for line in header:
        print(f"  {line}")
    print(f"Entries to remove : {len(paths)}")
    bashrc = Path.home() / ".bashrc"
    will_touch_bashrc = (
        not args.no_bashrc
        and bashrc.is_file()
        and BASHRC_MARKER in bashrc.read_text()
    )
    print(f"Bashrc cleanup    : {'yes (' + str(bashrc) + ')' if will_touch_bashrc else 'no'}")
    print()

    if not args.yes and not args.dry_run:
        ans = input("Proceed? [y/N] ").strip().lower()
        if ans not in ("y", "yes"):
            print("aborted.")
            sys.exit(0)

    print("Removing manifest entries:")
    removed = 0
    for p in paths:
        if remove_path(p, args.dry_run):
            removed += 1

    if will_touch_bashrc:
        print("Cleaning bashrc:")
        strip_bashrc_block(bashrc, args.dry_run)

    if not args.dry_run and manifest.exists():
        try:
            manifest.unlink()
            print(f"Removed manifest: {manifest}")
        except OSError as e:
            print(f"WARN: failed to remove manifest {manifest}: {e}", file=sys.stderr)

    print()
    print(f"Done. {removed}/{len(paths)} entries removed.")
    print("Source checkout (ecalj/) is untouched. To finish removing ecalj entirely:")
    print("  rm -rf <path-to-ecalj-checkout>")
    if args.dry_run:
        print("(dry run — nothing was actually changed.)")


if __name__ == "__main__":
    main()
