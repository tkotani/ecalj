#!/bin/bash
# GW1500 オプションラン (収束版): 既存 done 物質に対して
# (A) gwscconv で QSGW iter を収束 or max-iter まで継続
# (B) LDA band を新規計算 (元 worker.sh ではやってなかった)
# (C) QSGW80 band を新しい sigm で再計算
#
# 使い方:
#   bash run_gw1500_addrun_conv.sh <list_file> [N_WORKERS=6] [CONV_TOL=0.1] [MAX_ITER=10]
#
# slot 構成 (元 production と同じ):
#   - 2 CPU slot (lmf -np 30) × 2
#   - 2 GPU slot (hgw_combined_mp_gpu, hvccfp0_mp_gpu, hsfp0_sc_mp_gpu)
#   slot_scheduler_daemon が FIFO で割当。 worker 6 個が競合取得。

set -u

LIST=${1:?"usage: $0 <list_file> [N_WORKERS=6] [CONV_TOL=0.1] [MAX_ITER=10]"}
N_WORKERS=${2:-6}
CONV_TOL=${3:-0.1}
MAX_ITER=${4:-10}

EPATH=~/bin
WORKDIR=~/DATA/gw1500
POSCAR_DIR=~/ecaljdeveloper/ecalj_auto/INPUT/gw1500/POSCARALL
NCORE=30
NP2=1
LOG=$WORKDIR/addrun_conv.log

[ ! -f "$LIST" ] && { echo "ERROR: list not found: $LIST"; exit 1; }
mapfile -t MIDS < <(grep -vE '^\s*(#|$)' "$LIST" | awk '{print $1}')
echo "$(date '+%F %T') === addrun_conv start: ${#MIDS[@]} mpids, N_WORKERS=$N_WORKERS, CONV=$CONV_TOL, MAX=$MAX_ITER ===" | tee -a $LOG

# slot scheduler
if [ ! -S /tmp/slot_scheduler.sock ]; then
    rm -f /dev/shm/sem.OMPIO* /tmp/cpu_slot_*.lock /tmp/gpu_slot_*.lock /tmp/worker_*.state 2>/dev/null
    nohup python3 $EPATH/slot_scheduler_daemon.py > /tmp/slot_scheduler.log 2>&1 &
    SCHED_PID=$!
    echo "spawned slot_scheduler_daemon PID=$SCHED_PID" | tee -a $LOG
    for i in $(seq 1 30); do
        [ -S /tmp/slot_scheduler.sock ] && break
        sleep 0.2
    done
    [ -S /tmp/slot_scheduler.sock ] || { echo "FATAL: scheduler socket"; exit 1; }
    trap "kill $SCHED_PID 2>/dev/null; rm -f /tmp/slot_scheduler.sock" EXIT
else
    echo "slot_scheduler reused" | tee -a $LOG
fi

QUEUE=$(mktemp /tmp/addrun_conv_queue.XXXXXX)
printf '%s\n' "${MIDS[@]}" > "$QUEUE"
LOCK=${QUEUE}.lock

run_one() {
    local mpid=$1 wid=$2
    local dir=$WORKDIR/$mpid
    local t0=$(date +%s)

    # Create dir from POSCARALL if not present (scratch entry, e.g. REDO_FAIL
    # mids whose dir was wiped while FAILBACKUP preserves the original).
    if [ ! -d "$dir" ]; then
        local poscar=$POSCAR_DIR/POSCAR.$mpid
        if [ ! -e "$poscar" ]; then
            echo "$(date '+%F %T') $wid $mpid SKIP: no dir + no POSCAR.$mpid" | tee -a $LOG
            return 1
        fi
        mkdir -p "$dir"
        cp "$poscar" "$dir/POSCAR"
    fi
    cd "$dir"

    # If TOML inputs missing (legacy-era materials), regenerate from POSCAR via
    # vasp2ctrl + ctrlgenToml.py — same flow as production worker.sh used post
    # TOML migration. ctrlgenToml.py's default 8x8x8 k-mesh matches the
    # GW1500 production convention.
    if [ ! -e "ctrlG.$mpid.toml" ] || [ ! -e "PB.toml" ]; then
        if [ ! -e "POSCAR" ]; then
            echo "$(date '+%F %T') $wid $mpid SKIP: no POSCAR" | tee -a $LOG
            return 1
        fi
        $EPATH/vasp2ctrl POSCAR > l_vasp2ctrl 2>&1 || {
            echo "$(date '+%F %T') $wid $mpid SKIP: vasp2ctrl failed" | tee -a $LOG
            return 1
        }
        cp ctrls.POSCAR.vasp2ctrl "ctrls.$mpid"
        $EPATH/ctrlgenToml.py $mpid > l_ctrlgen 2>&1 || {
            echo "$(date '+%F %T') $wid $mpid SKIP: ctrlgenToml failed" | tee -a $LOG
            return 1
        }
        if [ ! -e "ctrlG.$mpid.toml" ] || [ ! -e "PB.toml" ]; then
            echo "$(date '+%F %T') $wid $mpid SKIP: ctrlgenToml output missing" | tee -a $LOG
            return 1
        fi
    fi

    # sanity
    if [ ! -e "PB.toml" ]; then
        echo "$(date '+%F %T') $wid $mpid SKIP: PB.toml missing" | tee -a $LOG
        return 1
    fi

    # If rst.<mid> missing → scratch material: run LDA SCF first.
    # (This covers REDO_FAIL mids whose mp-XXXX/ was reset to bare POSCAR+TOML,
    # plus any future scratch addition.)
    if [ ! -e "rst.$mpid" ]; then
        mpirun -np 1 $EPATH/lmfa $mpid > llmfa.scratch 2>&1
        if [ $? -ne 0 ]; then
            echo "$(date '+%F %T') $wid $mpid SKIP: lmfa scratch failed" | tee -a $LOG
            return 1
        fi
        $EPATH/slot_run.py cpu "$wid" lmf_lda_init "$mpid" -- \
            mpirun -np $NCORE $EPATH/lmf $mpid '-v[iter.nit]=80' > llmf_lda.scratch 2>&1
        local ldac=$(tail -1 save.$mpid 2>/dev/null | awk '{print $1}')
        if [ "$ldac" != "c" ] && [ "$ldac" != "x" ]; then
            echo "$(date '+%F %T') $wid $mpid SKIP: LDA-SCF not converged (status=$ldac)" | tee -a $LOG
            return 1
        fi
    fi

    # Restore QPU sentinels for gwscconv Iter0 detection (only if QSGW.*run/ exists)
    local maxn=0 d n
    for d in QSGW.*run; do
        [ -d "$d" ] || continue
        n=${d#QSGW.}; n=${n%run}
        if [[ "$n" =~ ^[0-9]+$ ]]; then
            touch "QPU.${n}run"
            (( n > maxn )) && maxn=$n
        fi
    done
    # maxn=0 (Iter0=0) is OK for scratch processing — gwscconv handles it.

    echo "$(date '+%F %T') $wid $mpid START Iter0=$maxn" | tee -a $LOG

    # ===== Phase A: gwscconv (continue until conv or max-iter) =====
    $EPATH/gwscconv -np $NCORE -np2 $NP2 --gpu --mp --fp32 $mpid '-v[ham.scaledsigma]=0.8' \
        --conv-tol $CONV_TOL --max-iter $MAX_ITER > osgw.conv.out 2>&1 &
    local pid=$!
    local kill_reason=""
    while kill -0 $pid 2>/dev/null; do
        sleep 300
        if grep -l "NaN" lsc lsx lqpe llmf 2>/dev/null | head -1 | grep -q .; then
            kill_reason="NaN"
            kill $pid 2>/dev/null
            pkill -P $$ -f "lmf $mpid" 2>/dev/null
            pkill -P $$ -f "hsfp0|hrcxq|hqpe|hvccfp0" 2>/dev/null
            break
        fi
        if [ $(( $(date +%s) - t0 )) -gt 57600 ]; then  # 16h hard cap
            kill_reason="TIMEOUT 16h"
            kill $pid 2>/dev/null
            pkill -P $$ -f "lmf $mpid" 2>/dev/null
            pkill -P $$ -f "hsfp0|hrcxq|hqpe|hvccfp0" 2>/dev/null
            break
        fi
    done
    wait $pid 2>/dev/null
    local phaseA_rc=$?
    if [ -n "$kill_reason" ]; then
        echo "$(date '+%F %T') $wid $mpid PHASE-A FAIL $kill_reason" | tee -a $LOG
        return 1
    fi
    if [ "$phaseA_rc" -ne 0 ]; then
        local errtail=$(tail -3 osgw.conv.out 2>/dev/null | tr '\n' ' ')
        echo "$(date '+%F %T') $wid $mpid PHASE-A FAIL rc=$phaseA_rc: $errtail" | tee -a $LOG
        return 1
    fi
    local conv_line
    conv_line=$(grep -E "CONVERGED|max_iter|max-iter|reached" osgw.conv.out | head -1)
    [ -n "$conv_line" ] && echo "$(date '+%F %T') $wid $mpid PHASE-A: $conv_line" | tee -a $LOG

    # ===== Phase B: LDA band (new dir, fresh LDA SCF) =====
    mkdir -p PlotBand_LDA
    cp "ctrlG.$mpid.toml" PB.toml POSCAR PlotBand_LDA/ 2>/dev/null
    if [ -f "PlotBand/syml.$mpid" ]; then
        cp "PlotBand/syml.$mpid" PlotBand_LDA/
    fi
    cd PlotBand_LDA
    mpirun -np 1 $EPATH/lmfa $mpid '-v[ham.scaledsigma]=0' > llmfa 2>&1
    if [ $? -ne 0 ]; then
        echo "$(date '+%F %T') $wid $mpid PHASE-B lmfa FAIL" | tee -a $LOG
        cd "$dir"; return 1
    fi
    $EPATH/slot_run.py cpu "$wid" lmf_lda_scf "$mpid" -- \
        mpirun -np $NCORE $EPATH/lmf $mpid '-v[ham.scaledsigma]=0' > llmf 2>&1
    if [ $? -ne 0 ]; then
        echo "$(date '+%F %T') $wid $mpid PHASE-B lmf SCF FAIL" | tee -a $LOG
        cd "$dir"; return 1
    fi
    if [ ! -f "syml.$mpid" ]; then
        $EPATH/getsyml $mpid --nobzview > lgetsyml 2>&1
    fi
    $EPATH/slot_run.py cpu "$wid" lmf_lda_band "$mpid" -- \
        mpirun -np $NCORE $EPATH/lmf $mpid '-v[ham.scaledsigma]=0' --band > llmf_band 2>&1
    cd "$dir"

    # ===== Phase C: QSGW80 band re-run with new sigm =====
    mkdir -p PlotBand
    cp rst.$mpid "ctrlG.$mpid.toml" PB.toml atmpnu.*.$mpid PlotBand/ 2>/dev/null || true
    cp sigm "sigm.$mpid" PlotBand/ 2>/dev/null || true
    cd PlotBand
    [ -f "syml.$mpid" ] || $EPATH/getsyml $mpid --nobzview > lgetsyml 2>&1
    # job_band 内部で run_cmd.py が slot 取得するので外側 slot_run 不要
    $EPATH/job_band $mpid -np $NCORE --NoGnuplot '-v[ham.scaledsigma]=0.8' > ljobband 2>&1
    cd "$dir"

    # cleanup intermediate
    $EPATH/cleargw . > /dev/null 2>&1
    rm -rf SEBK LDA STDOUT __*
    rm -f QPU.*run llmf.*run lsx lsc lsxC lrcxq

    local elapsed=$(( $(date +%s) - t0 ))
    echo "$(date '+%F %T') $wid $mpid OK ${elapsed}s" | tee -a $LOG
    return 0
}

worker_loop() {
    local wid=$1
    while true; do
        exec 9>$LOCK
        flock 9
        local mpid=$(head -1 $QUEUE)
        if [ -n "$mpid" ]; then
            sed -i '1d' $QUEUE
        fi
        exec 9>&-
        [ -z "$mpid" ] && break
        run_one "$mpid" "$wid"
    done
}

PIDS=()
for i in $(seq 1 $N_WORKERS); do
    worker_loop "AC$i" &
    PIDS+=($!)
done
wait "${PIDS[@]}"

rm -f "$QUEUE" "$LOCK"
echo "$(date '+%F %T') === addrun_conv done ===" | tee -a $LOG
