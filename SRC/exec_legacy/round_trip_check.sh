#!/bin/bash
# round_trip_check.sh - Validate gwinput2toml + tomlexpand idempotence at TOML level
#
# For each GWinput in $1 (find search root), do:
#   GWinput  -[gwinput2toml.py]->  toml1
#   toml1    -[tomlexpand.py]  ->  GWinput'
#   GWinput' -[gwinput2toml.py]->  toml2
#   diff toml1 toml2
#
# Prints PASS/FAIL per file and a summary.

set -u
ROOT="${1:-/home/takao/ecaljdeveloper/Samples}"
TOMLGEN="${TOMLGEN:-$HOME/bin/gwinput2toml.py}"
TOMLEXPAND="${TOMLEXPAND:-$HOME/bin/tomlexpand.py}"

pass=0
fail=0
failed_files=()

while IFS= read -r gwfile; do
  dir=$(dirname "$gwfile")
  workdir=$(mktemp -d)

  cp "$gwfile" "$workdir/GWinput"
  ( cd "$workdir" && python3 "$TOMLGEN" GWinput >/dev/null 2>&1 ) || {
    echo "FAIL (tomlgen step 1): $gwfile"
    fail=$((fail+1)); failed_files+=("$gwfile (tomlgen1)"); rm -rf "$workdir"; continue
  }
  cp "$workdir/GWinput.toml" "$workdir/toml1.toml"

  ( cd "$workdir" && python3 "$TOMLEXPAND" GWinput.toml -o GWinput.round >/dev/null 2>&1 ) || {
    echo "FAIL (tomlexpand): $gwfile"
    fail=$((fail+1)); failed_files+=("$gwfile (tomlexpand)"); rm -rf "$workdir"; continue
  }

  ( cd "$workdir" && python3 "$TOMLGEN" GWinput.round -o toml2.toml >/dev/null 2>&1 ) || {
    echo "FAIL (tomlgen step 2): $gwfile"
    fail=$((fail+1)); failed_files+=("$gwfile (tomlgen2)"); rm -rf "$workdir"; continue
  }

  # Compare on parsed-data level (ignores formatting + block ordering)
  python3 - "$workdir/toml1.toml" "$workdir/toml2.toml" <<'EOF' >/dev/null 2>&1
import sys, tomllib
a = tomllib.loads(open(sys.argv[1]).read())
b = tomllib.loads(open(sys.argv[2]).read())
sys.exit(0 if a == b else 1)
EOF
  if [ $? -eq 0 ]; then
    pass=$((pass+1))
  else
    fail=$((fail+1))
    failed_files+=("$gwfile")
    echo "FAIL (diff): $gwfile"
    diff "$workdir/toml1.toml" "$workdir/toml2.toml" | head -15 | sed 's/^/    /'
  fi
  rm -rf "$workdir"
done < <(find "$ROOT" -name "GWinput" -type f 2>/dev/null)

echo
echo "===== SUMMARY ====="
echo "PASS: $pass"
echo "FAIL: $fail"
if [ $fail -gt 0 ]; then
  echo "Failed files:"
  for f in "${failed_files[@]}"; do echo "  $f"; done
fi
