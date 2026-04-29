#!/bin/bash
# Slot-based scheduler: 2 CPU slots + 2 GPU slots, dynamic CUDA_VISIBLE_DEVICES
# N workers compete for slots via flock; no fixed GPU/worker binding
WORKDIR=~/DATA/gw1500
N_WORKERS=6

echo "=== GW1500 QSGW80 slot-scheduler started at $(date) ===" | tee $WORKDIR/run.log

# Cleanup stale MPI semaphores and slot locks
rm -f /dev/shm/sem.OMPIO* /tmp/cpu_slot_*.lock /tmp/gpu_slot_*.lock 2>/dev/null

PIDS=()
for i in $(seq 1 $N_WORKERS); do
    WORKER_ID=W$i bash $WORKDIR/worker.sh - W$i >> $WORKDIR/run.log 2>&1 &
    PIDS+=($!)
    echo "Worker W$i: PID=${PIDS[-1]}" | tee -a $WORKDIR/run.log
    echo "${PIDS[-1]}" > $WORKDIR/pid_W$i
done

wait "${PIDS[@]}"
echo "=== GW1500 QSGW80 finished at $(date) ===" | tee -a $WORKDIR/run.log
