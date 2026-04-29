# GW1500 Slot Scheduler

## 構成

- 2 CPU slots + 2 GPU slots (グローバル)
- 6 worker.sh 並列、GPU 固定割当なし
- 動的 CUDA_VISIBLE_DEVICES 設定

## ロック

| Lock | 対象 | 数 |
|---|---|---|
| /tmp/cpu_slot_{0,1}.lock | lmf | 2 |
| /tmp/gpu_slot_{0,1}.lock | hsfp0_sc_mp_gpu, hrcxq_mp_gpu, hvccfp0_mp_gpu | 2 |

## ファイル

- `~/bin2/run_cmd.py` — _run_mpi で flock + slot 取得
- `~/DATA/gw1500/worker.sh` — LDA SCF と job_band も flock CPU slot
- `~/DATA/gw1500/run_gw1500.sh` — N_WORKERS=6 起動
- `~/DATA/gw1500/slot_history.log` — acquire/release 全記録

## 起動

```bash
cd ~/DATA/gw1500
nohup bash run_gw1500.sh > run_$(date +%Y%m%d-%H%M).log 2>&1 &
```

## 停止 + 再起動 (clean)

```bash
# Kill
ps -ef | grep -E "worker.sh|run_gw1500" | grep -v grep | awk "{print \$2}" | xargs -r kill
sleep 3
pkill -9 -f "lmf mp-|mpirun.*mp-|gwsc.*mp-"
pkill -9 -f "hsfp0_sc|hrcxq|hvccfp0|hbasfp0|hqpe_sc|qg4gw|heftet"
sleep 3

# Cleanup locks/sems
rm -f /dev/shm/sem.OMPIO* /tmp/cpu_slot_*.lock /tmp/gpu_slot_*.lock /tmp/worker_*.state

# Cleanup half-done dirs (重要: QPU.*run も消すこと)
for mpid in <list>; do
  cd ~/DATA/gw1500/$mpid 2>/dev/null && {
    rm -f __WV* __PP* __BASFP* __atm.* __mixm* __mixsig __vxcevec* __GEIG __VXCFP \
          __BZDATA __CPHI __EValue __HAMindex* __MTOindex __PHIVC __QGcou __QGpsi \
          __Vcoud.* __PPOVLG* __PPBRD* __PPOVLGG SEX2U SEXcore2U SEC2U sigm \
          QPU.*run llmf.*run
    rm -rf QSGW.*run
    cd ~/DATA/gw1500
  }
done

# Restart
nohup bash run_gw1500.sh > run_$(date +%Y%m%d-%H%M).log 2>&1 &
```

## 監視

### スロット使用状況
```bash
for f in cpu_slot_0 cpu_slot_1 gpu_slot_0 gpu_slot_1; do
  pid=$(lsof "/tmp/${f}.lock" 2>/dev/null | tail -n +2 | awk "{print \$2}" | head -1)
  if [ -n "$pid" ]; then
    cwd=$(readlink /proc/$pid/cwd 2>/dev/null | xargs -I{} basename {})
    wid=$(cat /proc/$pid/environ 2>/dev/null | tr "\0" "\n" | grep ^WORKER_ID= | cut -d= -f2)
    echo "$f: $wid mpid=$cwd"
  else
    echo "$f: (free)"
  fi
done
```

### スロットヒストリの overlap audit
```bash
python3 -c "
import sys
slots={}; ov=[]
for line in open(\"slot_history.log\"):
    p = line.split()
    if len(p)<7 or p[4] not in (\"cpu\",\"gpu\"): continue
    k=(p[4], p[6].split(\"=\")[1])
    if p[5]==\"acquire\":
        if k in slots: ov.append((line.strip(), slots[k]))
        slots[k]=line.strip()
    elif p[5]==\"release\":
        slots.pop(k, None)
print(f\"overlaps={len(ov)}\")
"
```

## 性能上限

- CPU: 2 lmf × 30 cores = 60 cores 使用 (64 中)
- GPU: 2 GPUs フル稼働
- 期待スループット: 1物質あたり約 0.5-1h (8原子)、6ワーカーで 6-12 物質/h 理論値、実効 4物質/h 程度

## 既知の罠

- **QPU.\*run 残骸**: 再起動時に消さないと gwsc が iter 番号オフセット (label のみ、計算は valid)
- **NCORE x 同時 lmf 数 ≤ 64**: スロット 2 で守られる
- **Python プロセス kill**: fd 自動 close で flock 自動解放
