#!/bin/bash
# GW1500 オプションラン: 既存の done 物質に対して QSGW iter を NADD 回追加実行する。
#
# 使い方:
#   bash run_gw1500_addrun.sh <list_file> [NADD] [N_WORKERS]
#       list_file : 1行 1 mpid のテキスト (空行・# 行は無視)
#       NADD      : 追加 QSGW iter 数 (default 2)
#       N_WORKERS : 並列ワーカー数 (default 4)
#
# 前提:
#   - production の worker.sh が cleanup で QPU.{n}run を消すため、
#     gwsc の Iter0 検出が動かない。本スクリプトは QSGW.{n}run/ から
#     最大 n を割り出し、QPU.{n}run sentinel を作って gwsc に
#     "Iter0=n" を伝える。
#   - rst.<mid>, sigm (or sigm.<mid>), ctrlG.<mid>.toml, PB.toml が
#     mp-XXXX/ 直下に必要。 production 終了後の状態ならそのまま OK。
#
# 出力:
#   - 既存 QSGW.{n}run/ は保持、新たに QSGW.{n+1}run, QSGW.{n+2}run 追加
#   - rst.<mid>, sigm が in-place で更新される (オリジナルは QSGW.{n}run/ に
#     スナップショットとして残ってるので復元可)
#   - ログ: addrun.log と各 mp-XXXX/ 内 osgw.addrun.out

set -u

LIST=${1:?"usage: $0 <list_file> [NADD=2] [N_WORKERS=4]"}
NADD=${2:-2}
N_WORKERS=${3:-4}

EPATH=~/bin
WORKDIR=~/DATA/gw1500
NCORE=30
NP2=1
LOG=$WORKDIR/addrun.log

[ ! -f "$LIST" ] && { echo "ERROR: list file not found: $LIST"; exit 1; }

# Read mpids (skip blank / comment lines)
mapfile -t MIDS < <(grep -vE '^\s*(#|$)' "$LIST" | awk '{print $1}')
echo "$(date '+%F %T') === addrun start: ${#MIDS[@]} mpids, NADD=$NADD, N_WORKERS=$N_WORKERS ===" | tee -a $LOG

# Spawn slot_scheduler if not running (workers respect CPU/GPU slots)
if [ ! -S /tmp/slot_scheduler.sock ]; then
    rm -f /dev/shm/sem.OMPIO* /tmp/cpu_slot_*.lock /tmp/gpu_slot_*.lock /tmp/worker_*.state 2>/dev/null
    nohup python3 $EPATH/slot_scheduler_daemon.py > /tmp/slot_scheduler.log 2>&1 &
    SCHED_PID=$!
    echo "spawned slot_scheduler_daemon PID=$SCHED_PID" | tee -a $LOG
    for i in $(seq 1 30); do
        [ -S /tmp/slot_scheduler.sock ] && break
        sleep 0.2
    done
    [ -S /tmp/slot_scheduler.sock ] || { echo "FATAL: scheduler socket not created"; exit 1; }
    trap "kill $SCHED_PID 2>/dev/null; rm -f /tmp/slot_scheduler.sock" EXIT
else
    echo "slot_scheduler already running, reusing" | tee -a $LOG
fi

# Per-material work queue — use a tmp file as shared queue with flock
QUEUE_TMP=$(mktemp)
printf '%s\n' "${MIDS[@]}" > "$QUEUE_TMP"
LOCK_TMP=${QUEUE_TMP}.lock

run_one() {
    local mpid=$1
    local dir=$WORKDIR/$mpid
    local t0=$(date +%s)

    if [ ! -d "$dir" ]; then
        echo "$(date '+%F %T') $mpid SKIP: dir missing" | tee -a $LOG
        return 1
    fi
    cd "$dir"

    # Sanity check files
    for f in ctrlG.$mpid.toml PB.toml rst.$mpid; do
        if [ ! -e "$f" ]; then
            echo "$(date '+%F %T') $mpid SKIP: $f missing" | tee -a $LOG
            return 1
        fi
    done
    # sigm.$mpid OR sigm
    if [ ! -e "sigm.$mpid" ] && [ ! -e "sigm" ]; then
        echo "$(date '+%F %T') $mpid SKIP: sigm missing" | tee -a $LOG
        return 1
    fi

    # Restore QPU.{n}run sentinels from QSGW.{n}run dirs
    local maxn=0
    for d in QSGW.*run; do
        [ -d "$d" ] || continue
        local n=${d#QSGW.}; n=${n%run}
        if [[ "$n" =~ ^[0-9]+$ ]]; then
            touch "QPU.${n}run"
            (( n > maxn )) && maxn=$n
        fi
    done
    if [ "$maxn" -eq 0 ]; then
        echo "$(date '+%F %T') $mpid SKIP: no QSGW.*run/ found" | tee -a $LOG
        return 1
    fi

    echo "$(date '+%F %T') $mpid START addrun (Iter0=$maxn, NADD=$NADD)" | tee -a $LOG

    # Run gwsc with watchdog (NaN check + 8h timeout, mirroring worker.sh)
    $EPATH/gwsc $NADD -np $NCORE -np2 $NP2 --gpu --mp $mpid '-v[ham.scaledsigma]=0.8' > osgw.addrun.out 2>&1 &
    local gwsc_pid=$!
    local kill_reason=""
    while kill -0 $gwsc_pid 2>/dev/null; do
        sleep 300
        if grep -l "NaN" lsc lsx lqpe llmf 2>/dev/null | head -1 | grep -q .; then
            kill_reason="NaN detected"
            kill $gwsc_pid 2>/dev/null
            pkill -P $$ -f "lmf $mpid" 2>/dev/null
            pkill -P $$ -f "hsfp0|hrcxq|hqpe|hvccfp0" 2>/dev/null
            break
        fi
        if [ $(( $(date +%s) - t0 )) -gt 28800 ]; then
            kill_reason="TIMEOUT 8h"
            kill $gwsc_pid 2>/dev/null
            pkill -P $$ -f "lmf $mpid" 2>/dev/null
            pkill -P $$ -f "hsfp0|hrcxq|hqpe|hvccfp0" 2>/dev/null
            break
        fi
    done
    wait $gwsc_pid 2>/dev/null
    local last=$(tail -1 osgw.addrun.out 2>/dev/null)

    if [ -n "$kill_reason" ]; then
        echo "$(date '+%F %T') $mpid FAIL: $kill_reason" | tee -a $LOG
        return 1
    fi
    if echo "$last" | grep -q "All calclation finished"; then
        local elapsed=$(( $(date +%s) - t0 ))
        echo "$(date '+%F %T') $mpid OK ${elapsed}s (now iter=$((maxn+NADD)))" | tee -a $LOG
        return 0
    fi
    echo "$(date '+%F %T') $mpid FAIL: $last" | tee -a $LOG
    return 1
}

worker_loop() {
    local wid=$1
    while true; do
        # Atomic dequeue
        exec 9>$LOCK_TMP
        flock 9
        local mpid=$(head -1 $QUEUE_TMP)
        if [ -n "$mpid" ]; then
            sed -i '1d' $QUEUE_TMP
        fi
        exec 9>&-
        [ -z "$mpid" ] && break
        WORKER_ID=$wid run_one "$mpid"
    done
}

# Spawn N_WORKERS in parallel
PIDS=()
for i in $(seq 1 $N_WORKERS); do
    worker_loop "AW$i" &
    PIDS+=($!)
done
wait "${PIDS[@]}"

rm -f "$QUEUE_TMP" "$LOCK_TMP"
echo "$(date '+%F %T') === addrun done ===" | tee -a $LOG
