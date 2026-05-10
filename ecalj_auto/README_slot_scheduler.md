# GW1500 QSGW80 Production

## 概要

1546物質の QSGW80 (scaledsigma=0.8) を kt1 (Dell G15, 64 cores, 2×RTX5090) で量産。
TOML 入力 (ctrlG.<sname>.toml + PB.toml) + hgw_combined (in-memory W) で実行。

## ディレクトリ構成

```
~/DATA/gw1500/              作業ディレクトリ
  queue_all.txt             全1546物質リスト
  queue.txt                 残りキュー (mpid niter)
  done.log                  完了ログ
  failed.log                失敗ログ
  worker.sh                 ワーカースクリプト
  run_gw1500.sh             起動スクリプト (daemon + 6 workers)
  slot_history.log          スロット acquire/release 全記録
  mp-XXXX/                  各物質の作業ディレクトリ

~/ecaljdeveloper/ecalj_auto/INPUT/gw1500/POSCARALL/
  POSCAR.mp-XXXX            入力 POSCAR (1546ファイル)

~/bin/                      ecalj 本体は build/ への symlink (InstallAll.py 経由)
  slot_scheduler_daemon.py  スロットスケジューラ (Unix socket、 ecalj_auto/ への symlink)
  slot_run.py               スロット取得 → 子プロセス exec (ecalj_auto/ への symlink)
  run_cmd.py                MPI 実行 (スロット自動取得)
  gwsc                      QSGW ドライバ (hgw_combined 使用)
  ctrlgenToml.py            POSCAR → ctrlG.toml + PB.toml 生成
  clusters.toml             MPI launcher 設定

~/ecaljdeveloper/ecalj_auto/  GW1500 専用ヘルパの実体 (ecalj 本体ではない)
  slot_scheduler_daemon.py    (~/bin から symlink される)
  slot_run.py                 (同上)
  worker.sh / run_gw1500.sh   起動・ワーカースクリプトの canonical コピー
```

## スロット構成

| Slot | 対象バイナリ | 数 | 備考 |
|------|-------------|---|------|
| CPU  | lmf (np=30) | 2 | 30×2=60 cores ≤ 64 |
| GPU  | hgw_combined_mp_gpu, hvccfp0_mp_gpu, hsfp0_sc_mp_gpu | 2 | GPU0, GPU1 |

6 workers が daemon 経由で FIFO 順にスロットを取得。
その他のバイナリ (qg4gw, hbasfp0, heftet, hqpe_sc, lmfa) はスロット不要。

## ワークフロー (worker.sh, 1物質あたり)

```
1. POSCAR → vasp2ctrl → ctrls.<mpid>
2. ctrlgenToml.py <mpid>  → ctrlG.<mpid>.toml + PB.toml
   (内部で lmchk + lmfa + lmf --jobgw=0 + gwinit を実行)
3. lmf <mpid> -v[iter.nit]=80  (LDA SCF, CPU slot)
4. gwsc 5 -np 30 -np2 1 --gpu --mp <mpid> -v[ham.scaledsigma]=0.8
   gwsc 内部:
     lmf --jobgw=0 → qg4gw → lmf --jobgw=1
     heftet → hbasfp0 --job=3 → hvccfp0 --job=3 → hsfp0_sc --job=3 (core Sx)
     hbasfp0 --job=0 → hvccfp0 --job=0 (valence basis)
     hgw_combined --jobgw=1  (Sx + W + Sc, in-memory, GPU)
     hqpe_sc → lmf (QSGW SCF, CPU slot)
   × 5 iterations
5. job_band (バンドプロット)
6. cleargw (中間ファイル削除)
```

## 起動

```bash
cd ~/DATA/gw1500
nohup bash run_gw1500.sh > run_$(date +%Y%m%d-%H%M).log 2>&1 &
```

run_gw1500.sh が自動で:
- stale semaphore/lock/socket を掃除
- slot_scheduler_daemon を起動
- 6 workers (W1-W6) を起動

## 停止 (graceful)

```bash
# 全プロセス kill
pkill -9 -f run_gw1500.sh
pkill -9 -f 'worker\.sh'
pkill -9 -f slot_scheduler_daemon
pkill -9 -f 'gwsc.*mp-'
pkill -9 -f 'mpirun.*/bin/'
sleep 2

# Stale lock/socket 掃除
rm -f /tmp/slot_scheduler.sock /dev/shm/sem.OMPIO* \
      /tmp/cpu_slot_*.lock /tmp/gpu_slot_*.lock /tmp/worker_*.state
```

## キュー管理

```bash
# 状態確認
cd ~/DATA/gw1500
wc -l done.log queue.txt failed.log

# キュー再構築 (done/failed 以外を全てキューに)
awk '{print $1}' queue_all.txt | sort -u > /tmp/all.txt
awk '{print $1}' done.log | sort -u > /tmp/done.txt
awk '{print $1}' failed.log | sort -u > /tmp/fail.txt
comm -23 /tmp/all.txt /tmp/done.txt | comm -23 - /tmp/fail.txt | \
  awk '{print $1, 5}' > queue.txt

# 中断物質のディレクトリ削除 (再起動前)
while read mpid niter; do rm -rf "$mpid"; done < queue.txt
```

## 監視

ワーカー状態表:
```bash
cd ~/DATA/gw1500
for w in W1 W2 W3 W4 W5 W6; do
  mpid=$(tail -1 worker${w}.log | grep -oP 'mp-\d+')
  if [ -f "$mpid/osgw.out" ]; then
    iter=$(grep -c 'iteration end' "$mpid/osgw.out")
    phase=$(tail -1 "$mpid/osgw.out" | grep -oP "bin/\K[^ ]+" | head -1)
  fi
  echo "$w $mpid iter=$iter $phase"
done
```

## 既知の罠

- **QPU.\*run 残骸**: 再起動前に消す。残すと gwsc の iter 番号がオフセット
- **CPU oversubscription**: NCORE×同時lmf数 ≤ 64。CPU slot 2 で保証
- **GPU async 禁止**: WB.4 async(1) が sigm 破壊 (005ba221 で revert 済み)
- **nvfortran signal 11**: BUILD_MP=OFF にする。ビルド時 -j4 リトライ必要
- **/dev/shm/sem.OMPIO\* 残骸**: kill 後に必ず掃除。残ると hrcxq がハング
- **disk 監視**: worker.sh が 20GB 未満で自動停止。__WV* 残骸は cleargw で削除
