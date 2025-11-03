#!/bin/bash
cd "$(dirname "$0")/build"
for f in *; do
  if [ -f "$f" ] && [ -x "$f" ]; then
    target="../$f"
    [ -e "$target" ] && rm -f "$target"
    cp "$f" "$target"
  fi
done
