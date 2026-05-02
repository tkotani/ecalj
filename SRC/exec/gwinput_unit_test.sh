#!/bin/bash
# Unit-test harness for m_GWinput TOML reader.
# For each Samples/*/GWinput, convert to GWinput.toml via gwinput2toml.py,
# then load via gwinput_dump and assert no parse errors.
#
# Usage: ./gwinput_unit_test.sh [path_to_gwinput_dump]
set -e

ECALJ_DIR="${ECALJ_DIR:-$HOME/ecaljdeveloper}"
DUMP="${1:-$ECALJ_DIR/SRC/exec/build/gwinput_dump}"
CONV="$ECALJ_DIR/SRC/exec/gwinput2toml.py"

if [ ! -x "$DUMP" ]; then
  echo "ERROR: gwinput_dump not found at $DUMP" >&2
  exit 1
fi

WORK=$(mktemp -d /tmp/gwinput_unit_test.XXXXXX)
trap "rm -rf $WORK" EXIT

PASS=0
FAIL=0
FAIL_NAMES=()

# Iterate over every GWinput in Samples (skip *_work directories used during tests)
while IFS= read -r gw; do
  name=$(echo "$gw" | sed "s|$ECALJ_DIR/Samples/||" | sed 's|/GWinput||' | tr '/' '_')
  cd "$WORK"
  rm -f GWinput GWinput.toml
  cp "$gw" ./GWinput

  if ! python3 "$CONV" GWinput >/dev/null 2>&1; then
    echo "FAIL [convert]: $name"
    FAIL=$((FAIL+1)); FAIL_NAMES+=("$name [convert]")
    continue
  fi

  out=$("$DUMP" GWinput.toml 2>&1) || true
  if echo "$out" | grep -q "^OK!" && ! echo "$out" | grep -q "Error:"; then
    PASS=$((PASS+1))
  else
    FAIL=$((FAIL+1)); FAIL_NAMES+=("$name [load]")
    echo "FAIL [load]: $name"
    echo "$out" | head -5
  fi
done < <(find "$ECALJ_DIR/Samples" -name GWinput | grep -v "_work" | sort)

echo ""
echo "===== Result ====="
echo "PASS: $PASS"
echo "FAIL: $FAIL"
if [ $FAIL -gt 0 ]; then
  echo "Failed cases:"
  for n in "${FAIL_NAMES[@]}"; do echo "  - $n"; done
  exit 1
fi
echo "ALL OK"
