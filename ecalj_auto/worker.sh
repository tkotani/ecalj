#!/bin/bash
# Worker script: picks materials from shared queue, runs LDA+QSGW
# Usage: worker.sh <gpu_id> <worker_id>

GPU_ID=$1
WORKER_ID=$2
EPATH=~/bin2
WORKDIR=~/DATA/gw1500
POSCAR_DIR=~/ecaljdeveloper/ecalj_auto/INPUT/gw1500/POSCARALL
QUEUE=$WORKDIR/queue.txt
LOCKFILE=$WORKDIR/queue.lock
NCORE=30
NP2=1
NITER=5
LOG=$WORKDIR/worker${WORKER_ID}.log

# CUDA_VISIBLE_DEVICES is set per-subprocess inside run_cmd.py (dynamic GPU slot binding)

echo "$(date '+%Y-%m-%d %H:%M:%S') Worker$WORKER_ID started GPU=$GPU_ID" | tee -a $LOG

pick_next() {
    # Atomically pick next material from queue using flock
    local mpid=""
    exec 9>$LOCKFILE
    flock 9
    mpid=$(head -1 $QUEUE)
    if [ -n "$mpid" ]; then
        sed -i '1d' $QUEUE
    fi
    exec 9>&-
    echo "$mpid"
}

run_material() {
    local mpid=$1
    local dir=$WORKDIR/$mpid
    local t0=$(date +%s)
    # Disk safety: exit if <20GB free
    avail_kb=$(df --output=avail / | tail -1)
    if [ "$avail_kb" -lt 20971520 ]; then
        exec 9>$LOCKFILE; flock 9
        sed -i "1i $mpid" $QUEUE
        exec 9>&-
        echo "$(date +%F' '%T) Worker$WORKER_ID: disk low (${avail_kb}KB), put $mpid back and exiting" | tee -a $LOG
        exit 99
    fi
    mkdir -p $dir
    cd $dir

    echo "$(date '+%Y-%m-%d %H:%M:%S') START $mpid" | tee -a $LOG

    # --- LDA ---
    if [ ! -f rst.$mpid.lda ]; then
        # Setup
        if [ ! -f atm.$mpid ]; then
            cp $POSCAR_DIR/POSCAR.$mpid POSCAR || return 1
            $EPATH/vasp2ctrl POSCAR > lvasp2ctrl 2>&1 || { echo "ERROR vasp2ctrl"; return 1; }
            cp ctrls.POSCAR.vasp2ctrl ctrls.$mpid
            $EPATH/ctrlgenM1.py $mpid > lctrlgen 2>&1 || { echo "ERROR ctrlgen"; return 1; }
            if [ ! -f ctrlgenM1.ctrl.$mpid ]; then echo "ERROR no ctrl"; return 1; fi
            cp ctrlgenM1.ctrl.$mpid ctrl.$mpid
            mpirun -np 1 $EPATH/lmchk $mpid > llmchk 2>&1
            mpirun -np 1 $EPATH/lmfa $mpid > llmfa 2>&1 || { echo "ERROR lmfa"; return 1; }
        fi
        # LDA SCF: acquire one of 2 CPU slots
        (
            t0=$(date +%s)
            while true; do
                for i in 0 1; do
                    exec 7>/tmp/cpu_slot_$i.lock
                    if flock -nx 7; then
                        wait_s=$(( $(date +%s) - t0 ))
                        echo "$(date '+%Y-%m-%d %H:%M:%S') $WORKER_ID $mpid cpu acquire slot=$i bin=lmf_lda wait=${wait_s}s" >> ~/DATA/gw1500/slot_history.log
                        mpirun -np $NCORE $EPATH/lmf $mpid -vnit=80 > llmf_lda 2>&1
                        rc=$?
                        echo "$(date '+%Y-%m-%d %H:%M:%S') $WORKER_ID $mpid cpu release slot=$i bin=lmf_lda" >> ~/DATA/gw1500/slot_history.log
                        exit $rc
                    fi
                    exec 7>&-
                done
                sleep 1
            done
        )
        local lda_status=$(tail -1 save.$mpid 2>/dev/null | awk '{print $1}')
        if [ "$lda_status" != "c" ] && [ "$lda_status" != "x" ]; then
            echo "ERROR lda_conv=$lda_status"
            return 1
        fi
        cp rst.$mpid rst.$mpid.lda
        echo "$(date '+%Y-%m-%d %H:%M:%S') LDA_DONE $mpid status=$lda_status" | tee -a $LOG
    else
        echo "$(date '+%Y-%m-%d %H:%M:%S') LDA_SKIP $mpid (rst exists)" | tee -a $LOG
    fi

    # --- QSGW ---
    $EPATH/mkGWinput $mpid > lgwin 2>&1 || { echo "ERROR mkGWinput"; return 1; }
    cp GWinput.tmp GWinput

    # In-flight watchdog: NaN check every 5min, 8h hard timeout
    $EPATH/gwsc $NITER -np $NCORE -np2 $NP2 --gpu --mp $mpid -vssig=0.8 > osgw.out 2>&1 &
    local gwsc_pid=$!
    local kill_reason=""
    while kill -0 $gwsc_pid 2>/dev/null; do
        sleep 300
        # NaN check during execution
        if grep -l "NaN" lsc lsx lqpe llmf 2>/dev/null | head -1 | grep -q .; then
            kill_reason="NaN detected in $(grep -l "NaN" lsc lsx lqpe llmf 2>/dev/null | head -1) (in-flight)"
            kill $gwsc_pid 2>/dev/null
            pkill -P $$ -f "lmf $mpid" 2>/dev/null
            pkill -P $$ -f "hsfp0|hrcxq|hqpe|hvccfp0" 2>/dev/null
            break
        fi
        # 8h hard timeout
        if [ $(( $(date +%s) - t0 )) -gt 28800 ]; then
            kill_reason="TIMEOUT after 8h"
            kill $gwsc_pid 2>/dev/null
            pkill -P $$ -f "lmf $mpid" 2>/dev/null
            pkill -P $$ -f "hsfp0|hrcxq|hqpe|hvccfp0" 2>/dev/null
            break
        fi
    done
    wait $gwsc_pid 2>/dev/null
    local gwsc_rc=$?
    local gwsc_last=$(tail -1 osgw.out 2>/dev/null)
    if [ -n "$kill_reason" ]; then
        echo "ERROR gwsc: $kill_reason"
        return 1
    fi
    # Post-run NaN check (in case gwsc finished but with NaN)
    if grep -l "NaN" lsc lsx lqpe llmf QSGW.*run/l* 2>/dev/null | head -1 | grep -q .; then
        nanfile=$(grep -l "NaN" lsc lsx lqpe llmf QSGW.*run/l* 2>/dev/null | head -1)
        echo "ERROR gwsc: NaN detected in $nanfile"
        return 1
    fi
    if echo "$gwsc_last" | grep -q "All calclation finished"; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') QSGW_DONE $mpid" | tee -a $LOG
    else
        echo "ERROR gwsc: $gwsc_last"
        return 1
    fi

    # --- Band plot ---
    mkdir -p PlotBand
    cp rst.$mpid ctrl.$mpid atmpnu.*.$mpid PlotBand/ 2>/dev/null || true
    cp sigm sigm.$mpid PlotBand/ 2>/dev/null || true
    cd PlotBand
    python3 $EPATH/getsyml $mpid --nobzview > lgetsyml 2>&1
    (
        t0=$(date +%s)
        while true; do
            for i in 0 1; do
                exec 7>/tmp/cpu_slot_$i.lock
                if flock -nx 7; then
                    wait_s=$(( $(date +%s) - t0 ))
                    echo "$(date '+%Y-%m-%d %H:%M:%S') $WORKER_ID $mpid cpu acquire slot=$i bin=job_band wait=${wait_s}s" >> ~/DATA/gw1500/slot_history.log
                    $EPATH/job_band $mpid -np $NCORE --NoGnuplot -vssig=0.8 > ljobband 2>&1
                    rc=$?
                    echo "$(date '+%Y-%m-%d %H:%M:%S') $WORKER_ID $mpid cpu release slot=$i bin=job_band" >> ~/DATA/gw1500/slot_history.log
                    exit $rc
                fi
                exec 7>&-
            done
            sleep 1
        done
    )
    cd $dir

    # --- Cleanup large GW temp files ---
    $EPATH/cleargw . > /dev/null 2>&1

    # --- Cleanup intermediate files (keep rst/sigm in QSGW.Xrun) ---
    rm -rf SEBK LDA STDOUT __*
    rm -f QPU.*run llmf.*run lsx lsc lsxC lrcxq
    rm -f GWinput GWinput.tmp *.chk ctrlgenM1.* ctrls.* ctrlp.*
    rm -f lbas lbasC leftet llmfgw00 llmfgw01 lvcc lvccC EFERMI
    rm -f SiteInfo.* @MNLA_* NLAindx.* PlatQlat.* QPLIST.* QBZ.* hbe.* freq_r estaticpot.dat efermi.lmf
    for run in QSGW.*run; do
      [ -d "$run" ] || continue
      find "$run" -type f ! -name "rst.*" ! -name "sigm.*" -delete
    done

    local t1=$(date +%s)
    local elapsed=$(( t1 - t0 ))
    echo "$(date '+%Y-%m-%d %H:%M:%S') ALL_DONE $mpid ${elapsed}s" | tee -a $LOG
    return 0
}

# Main loop
while true; do
    mpid=$(pick_next)
    if [ -z "$mpid" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') Worker$WORKER_ID: queue empty, exiting" | tee -a $LOG
        break
    fi

    errmsg=$(run_material "$mpid" 2>&1)
    rc=$?

    if [ $rc -eq 0 ]; then
        # Success
        exec 9>$LOCKFILE; flock 9
        echo "$mpid $(date '+%Y-%m-%d %H:%M:%S')" >> $WORKDIR/done.log
        exec 9>&-
    else
        # Failure: cleanup intermediate files + log single-line summary
        errsummary=$(echo "$errmsg" | grep -E "ERROR|TIMEOUT|NaN|disk low" | head -1)
        [ -z "$errsummary" ] && errsummary=$(echo "$errmsg" | tail -1)
        exec 9>$LOCKFILE; flock 9
        echo "$mpid $(date '+%Y-%m-%d %H:%M:%S') $errsummary" >> $WORKDIR/failed.log
        exec 9>&-
        echo "$(date '+%Y-%m-%d %H:%M:%S') FAILED $mpid: $errsummary" | tee -a $LOG
        # Cleanup huge intermediate files in failed dir
        if [ -d "$WORKDIR/$mpid" ]; then
            cd "$WORKDIR/$mpid"
            rm -f __WV* __PP* __BASFP* __atm.* __mixm* __mixsig __vxcevec* __GEIG __VXCFP __BZDATA __CPHI __EValue __HAMindex* __MTOindex __PHIVC __QGcou __QGpsi __Vcoud.* __PPOVLG* __PPBRD* __PPOVLGG SEX2U SEXcore2U SEC2U 2>/dev/null
            cd "$WORKDIR"
        fi
    fi
done

echo "$(date '+%Y-%m-%d %H:%M:%S') Worker$WORKER_ID finished" | tee -a $LOG
