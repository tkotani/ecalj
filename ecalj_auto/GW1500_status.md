# GW1500 の状況（2026-09-30）

GW1500 は 1546 物質の QSGW80（`scaledsigma = 0.8`）の計算で、2026-04〜05 に kt1 で回した。ここには物質ごとの最終の状態、
落ちたものの分類と原因、2026-09-30 に始めた回し直しを書く。物質ごとの表は [`gw1500_status_20260930.tsv`](gw1500_status_20260930.tsv)
（1546 行）。計算の仕組み（スロットスケジューラなど）は [`README_slot_scheduler.md`](README_slot_scheduler.md)。

## 1. 計算の設定と経過

- 入力: Materials Project の構造（`INPUT/gw1500/POSCARALL/POSCAR.<mpid>`、1〜8 原子）から `vasp2ctrl` と `ctrlgenToml.py`（4 月は `ctrlgenM1.py`）で作る。
  `nkabc = [8, 8, 8]`、`n1n2n3 = [4, 4, 4]`、`QpGcut_psi = 4.0`、`QpGcut_cou = 3.0`、`pwmode = 11`、`pwemax = 3.0`、`xcfun = 1`、スピン分極なし
- データ: kt1 の `~/DATA/gw1500/<mpid>/`（288 GB）。GPU の精度はすべて TF32（`--gpu --mp`）

**表 1**. 経過

| 期間（2026） | 中身 | 結果 |
| --- | --- | --- |
| 04-03〜05-10 | 本計算（production）。LDA から `gwsc 5` | 完了 1484、失敗 62（`done.log`、`failed.log`） |
| 05-11〜05-27 | 追加の計算（AC）。本計算の続きを `gwscconv` で、ギャップが 0.1 eV で収束するまで（最大 10 反復）。失敗した 122 物質は最初からやり直し（REDO） | 表 2 |
| 06-04 | 精度が原因で落ちた 36 物質を fp32 で 1 反復 | 36 物質とも NaN なし |
| 09-30〜 | 落ちたものを今のコードで最初から回し直し（§4） | §4 の表 |

## 2. 最終の状態

**表 2**. 1546 物質の最終の状態（2026-09-30 に、各物質のディレクトリの `osgw.conv.out` とログから集計）

| 状態 | 数 | 意味 |
| --- | --- | --- |
| GOOD | 1210 | 収束した（最後の 3 反復のギャップの幅が 0.1 eV 未満） |
| NOTCONV | 179 | 10 反復まで正常に回ったが、収束の条件を満たさない。落ちてはいない |
| FAILED | 147 | 使える結果が無い。回し直しが要る |
| UNKNOWN | 10 | `lmf` がギャップを出さない（LDA でもギャップが無いか 0.3 eV 未満）。金属か半金属で、ギャップの条件では判定できない |

- 追加の計算のログ（`addrun_conv.log`）の「OK」は収束を意味しない。OK の 1388 物質のうち `CONVERGED` は 1207、179 は最大反復に達しただけ。
  判定は各物質の `osgw.conv.out` にある
- 最初の追加の計算のログ `addrun_conv.log.broken`（05-10〜11）は無効（`gwscconv` が `-np 30` の 30 を反復の数と取り、何も回さずに「OK」と書いた）
- GOOD は「ギャップが収束した」という意味で、値が正しいことは確かめていない。すべて TF32 の計算で、fp32 や CPU との比較は無い。
  **mp-546711（CsClO₄）は LDA のギャップ 5.41 eV に対して QSGW80 が 0.26 eV で、誤り**（回し直しに入れた）
- 4 月の本計算は旧形式の入力（`ctrl` + `GWinput`）、追加の計算は 05-11 に作った TOML の入力で、Rb・Cs・K・Ba・Sr などの MT 半径が違う（3.0 と 2.8）。
  GOOD の収束値は TOML の入力のもの

## 3. 落ちたものの分類と原因

**表 3**. FAILED の 147 物質と UNKNOWN の 10 物質

| 分類 | 数 | 起きたこと | 原因 |
| --- | --- | --- | --- |
| 本計算: `lsc`・`lqpe` の NaN | 30 | Σc が NaN、または 1 反復前から Σ の大きさ（`sigma_m`）が 1〜2 桁大きい | GPU の行列積の精度（TF32）。fp32 の 1 反復では 30 物質とも正常 |
| 本計算: 「lmf failed even after reducing bmix to minimum」 | 6 | Σ が壊れていて、`lmf` が最初のバンド計算で止まる（混合の問題ではない） | 上と同じ（4 物質）。2 物質は前の実行の壊れた `sigm` を読んで止まった |
| 本計算: 8 時間の制限 | 22 | すべて 8 原子。2〜4 反復で打ち切り。エラーは無い | 遅いだけ（1 反復 1.5〜2.5 時間。GPU の順番待ちを含む） |
| 本計算: 停止時に実行中 | 3 | 05-10 に全体を止めたときに回っていた | — |
| 本計算: `lmf --jobgw=1` が落ちた | 1 | mp-1179832（Rb 8 原子）。プロセスが signal 9 で終了 | メモリ不足と見ている（確かめていない） |
| 追加: `rst` が入力に合わない | 44 | `llmf_start` に `File mismatch: rmt`。ログには「bmix to minimum」と出る | 04-03・04 に作った `rst` の MT 半径が 3.0、05-11 の TOML の入力は 2.8 |
| 追加: ギャップが消える（rc=3） | 38 | LDA ではギャップがある（0.17〜7.9 eV）のに、QSGW の途中で `lmf` がギャップを出さなくなる | 12 物質ほどは `sigma_m` が異常に大きい（精度）。残りは分かっていない |
| 追加: NaN | 2 | 本計算の NaN と同じ | 同じ |
| 追加: SCF が発散 | 1 | mp-18921。反復 7 の `lmf` が 80 回で収束せず、次の `lmf` が止まった | 分かっていない |
| UNKNOWN | 10 | LDA でもギャップが無い（8）か 0.3 eV 未満（2） | 金属・半金属 |

- 追加の計算の 2 回目（05-11〜12）では 83 物質が GPU のメモリ不足などで落ちた。6 本のワーカーの GPU のプログラムがすべて GPU 0 に載っていたためで
  （GPU 1 は空いていた）、スクリプトを直した 3 回目で 78 物質が収束した。表 3 には入っていない
- 物質の名前の一覧は `gw1500_status_20260930.tsv` の `final_state`・`final_detail` の列

## 4. 回し直し（2026-09-30〜）

FAILED の 147 物質と mp-546711 を、今のコードで最初から回し直している。

- スクリプト: [`gw1500_rerun.sh`](gw1500_rerun.sh)（キューから 1 物質ずつ取り、POSCAR → `vasp2ctrl` → `ctrlgenToml.py --ssig=0.8` → `gwscconv`）、
  バンドの図は [`gw1500_bandplot.sh`](gw1500_bandplot.sh)
- 設定: `--prec=fp32`、収束の条件 0.1 eV（2 回続けて）、最大 10 反復、1 物質 8 時間まで。`[gw]` はテンプレートのまま（`t_tetrakbt = 300`、`t_sigmaw = 300`）。
  4〜5 月の計算は χ0 が T = 0、`esmr = 0.003` Ry（262 K に当たる）で、ギャップのある物質ではこの違いは 1 meV 程度
- 場所: kt1 の `/mnt/data1/gw1500_rerun/run1/<mpid>/`、ログは `run1/rerun.log`。バイナリは `~/bin_frozen_b81da2342`（実体のコピー）
- 順番: 原子数の少ないものから

```bash
# kt1 での起動（4 本のワーカー、1 本に 12 コア）
export RUN_DIR=/mnt/data1/gw1500_rerun/run1 BIN=$HOME/bin_frozen_b81da2342 NP=12 GPUS=0,1 PREC=fp32 LIMIT=28800
export ENVSH=<PATH と LD_LIBRARY_PATH を設定するファイル>
for i in 0 1 2 3; do
  a=$((16+12*i)); nohup setsid taskset -c $a-$((a+11)) bash gw1500_rerun.sh queue.txt W$i > worker_W$i.out 2>&1 < /dev/null &
done
```

**複数のワーカーを回すときは、ワーカーごとに別のコアを割り当て（`taskset`）、ランクを固定しない（`OMPI_MCA_hwloc_base_binding_policy=none`。
`gw1500_rerun.sh` が設定する）。** OpenMPI は、起動した側の `taskset` にかかわらずランク i を機械のコア i に固定する。固定したままだと、
どのワーカーのランクもコア 0〜11 に、1 ランクの GPU のプログラムはコア 0 に重なる（2026-09-30: ほとんど進まなくなり、一度は全体が 17 分止まった）。
固定しなければランクはワーカーの `taskset` の中にとどまる。GPU のプログラムは `pylib/run_cmd.py` の GPU のロックで順番に GPU を使う。

**表 4**. 回し直しの結果（`rerun.log` から。回し終わったものから足す）

| 状態 | 数 |
| --- | --- |
| （実行中） | |

## 5. 残っていること

- NOTCONV の 179 物質: 続きを回す（入力を今の形に直し、5 月の `rst`・`sigm` を今の `lmf` が読めるかを確かめる）か、最初から 20 反復まで回す。どちらも試していない
- UNKNOWN の 10 物質: ギャップでは判定できないので、QP エネルギーの変化で収束を見る
- GOOD の値の検証（fp32 との比較）
- kt1 のディスク（`/` の空きは 50 GB）: `~/DATA/gw1500` の `*_FAILBACKUP`（121 個、69 GB）と、落ちた物質のディレクトリに残った作業ファイル `__*`（59 GB）、
  `~/DATA/gw1500fp32`・`gw1500_test`・`gw1500_bench`（合わせて 130 GB、ほとんどが `__*`）。消してよいかは未定（何も消していない）

## 6. `gw1500_status_20260930.tsv` の列

| 列 | 中身 |
| --- | --- |
| `mpid`、`formula`、`natom` | POSCAR から |
| `prod_status`、`prod_category`、`prod_date`、`prod_input`、`prod_bindir`、`prod_iter` | 本計算の結果、入力の形式、バイナリ、終わった反復 |
| `ac_status`、`ac_when`、`redo_status`、`redo_when`、`n_ac_runs` | 追加の計算（続き）と、最初からのやり直しの最新の結果 |
| `qsgw_iter_final`、`gap_history_eV`、`last3_swing_eV` | 最後の反復、最新の実行のギャップの履歴、最後の 3 反復の幅 |
| `gap_LDA_eV`、`gap_QSGW80_last_llmf_eV` | LDA のギャップと、最後の `lmf` のギャップ（**収束した値は GOOD の行だけ**） |
| `sigma_m` | `lqpe` の `sum check of sigma_m`（Σ の大きさ。収束した 1207 物質では 4〜3307） |
| `final_state`、`final_detail`、`rerun_action` | 表 2・表 3 の分類 |
