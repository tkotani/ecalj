#!/bin/bash
# Slot-based scheduler: 2 CPU slots + 2 GPU slots, dynamic CUDA_VISIBLE_DEVICES
# N workers compete for slots via flock; no fixed GPU/worker binding
WORKDIR=~/DATA/gw1500
N_WORKERS=6

echo "=== GW1500 QSGW80 slot-scheduler started at $(date) ===" | tee $WORKDIR/run.log

# Cleanup stale MPI semaphores, slot locks, and stale scheduler socket
rm -f /dev/shm/sem.OMPIO* /tmp/cpu_slot_*.lock /tmp/gpu_slot_*.lock /tmp/slot_scheduler.sock /tmp/worker_*.state 2>/dev/null

# Start slot scheduler daemon (FIFO ordered slot allocation)
nohup python3 ~/bin2/slot_scheduler_daemon.py > /tmp/slot_scheduler.log 2>&1 &
SCHED_PID=$!
echo "slot_scheduler_daemon PID=$SCHED_PID" | tee -a $WORKDIR/run.log
# wait for socket to appear
for i in $(seq 1 30); do
    [ -S /tmp/slot_scheduler.sock ] && break
    sleep 0.2
done
[ -S /tmp/slot_scheduler.sock ] || { echo "FATAL: scheduler socket not created"; exit 1; }
trap "kill $SCHED_PID 2>/dev/null; rm -f /tmp/slot_scheduler.sock" EXIT

PIDS=()
for i in $(seq 1 $N_WORKERS); do
    WORKER_ID=W$i bash $WORKDIR/worker.sh - W$i >> $WORKDIR/run.log 2>&1 &
    PIDS+=($!)
    echo "Worker W$i: PID=${PIDS[-1]}" | tee -a $WORKDIR/run.log
    echo "${PIDS[-1]}" > $WORKDIR/pid_W$i
done

wait "${PIDS[@]}"
echo "=== GW1500 QSGW80 finished at $(date) ===" | tee -a $WORKDIR/run.log
