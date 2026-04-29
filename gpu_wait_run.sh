#!/bin/bash
# gpu_wait_run.sh — Wait for GPU lock, then run command
# Uses the same /tmp/gpu.lock as the job script
LOCKFILE=/tmp/gpu.lock
(
  flock 9
  echo "GPU lock acquired. Running: $*"
  "$@"
) 9>"$LOCKFILE"
