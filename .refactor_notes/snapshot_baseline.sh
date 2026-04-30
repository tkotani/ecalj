#!/bin/bash
# Snapshot of internals for sxcf refactor regression check.
# Usage:
#   ./snapshot_baseline.sh capture <label>   — record snapshot of last test run
#   ./snapshot_baseline.sh diff <label1> <label2>  — compare two snapshots
#
# The snapshot captures, for each gwsc test workdir:
#   - md5 of SECU, SEC2U, SEXU, SEX2U (binary self-energy outputs)
#   - selected diagnostic lines from lcombined (ntq, nbandmx, ncount, etc.)
#
# By default scans all *_gwsc_work directories under TestInstall.
set -uo pipefail
# Note: -e is intentionally OMITTED; grep returning non-zero on no-match would
# kill this loop. Each command tolerates missing files individually.
TESTROOT=${TESTROOT:-$HOME/ecaljdeveloper/Samples/TestInstall}
SNAPDIR=${SNAPDIR:-$HOME/ecaljdeveloper/.refactor_notes/snapshots}
mkdir -p "$SNAPDIR"

capture_one() {
    local workdir="$1" outfile="$2"
    cd "$workdir"
    {
        echo "=== md5 of self-energy outputs (in SEBK after gwsc) ==="
        for f in SEBK/SECU SEBK/SEC2U SEBK/SEXU SEBK/SEX2U; do
            if [ -f "$f" ]; then
                md5sum "$f"
            else
                echo "(missing) $f"
            fi
        done
        echo
        echo "=== diagnostic lines from lcombined ==="
        if [ -f lcombined ]; then
            grep -E "ntq=|nspin nq ntq|=nbandmx|ncount=|MPI: worker_intask|Imag omega mesh|Real omega mesh|nblochpmx|niw" lcombined 2>/dev/null | head -50 || true
        fi
        echo
        echo "=== first 100 SECU rows (text, fp-precise) ==="
        if [ -f SEBK/SECU ]; then
            head -110 SEBK/SECU | tail -100
        fi
        echo
        echo "=== last 50 lines of lcombined ==="
        if [ -f lcombined ]; then
            tail -50 lcombined
        fi
    } > "$outfile"
}

case "${1:-}" in
    capture)
        label="${2:?usage: capture <label>}"
        outroot="$SNAPDIR/$label"
        rm -rf "$outroot"; mkdir -p "$outroot"
        for w in "$TESTROOT"/*_gwsc_work; do
            [ -d "$w" ] || continue
            name=$(basename "$w")
            capture_one "$w" "$outroot/$name.txt"
            echo "captured $name -> $outroot/$name.txt"
        done
        ;;
    diff)
        a="${2:?usage: diff <a> <b>}"
        b="${3:?usage: diff <a> <b>}"
        for f in "$SNAPDIR/$a"/*.txt; do
            name=$(basename "$f")
            other="$SNAPDIR/$b/$name"
            [ -f "$other" ] || { echo "missing $other"; continue; }
            if diff -q "$f" "$other" >/dev/null; then
                echo "OK  $name"
            else
                echo "DIFF $name:"
                diff "$f" "$other" | head -20
            fi
        done
        ;;
    *)
        echo "usage: $0 capture <label> | diff <label1> <label2>" >&2
        exit 1
        ;;
esac
