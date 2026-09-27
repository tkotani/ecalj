# bench_hgw — `hgw` だけを 1 回回すベンチ（kt1、2026-09-27〜28）

GW の GPU 高速化（報告 `../gpu_fp32_report.md`、ecaljdoc の ForDevelopers §11）で使った道具。kt1 に置いてあったものを 2026-09-28 にここへ写した。

| ファイル | 中身 |
| --- | --- |
| `bench3.sh <tag> <np>` | LiTi₂O₄ の 10 反復後の状態で `hgw` を 1 回だけ回し、`res_<tag>/` に `lgw` と Σ のファイル、`bench.log` に 1 行（時間）。`np` は GPU の数（1 ランク 1 枚）。`GPUS`（既定 0）、`BENCH_BIN`、`BENCH_PREC`、`BENCH_DIR` |
| `cmpef.py <ref> <res>...` | Re Σc（`SECU`）と Σx（`SEXU`）の基準との差を、全状態・$\vert E-E_F\vert$ < 3 eV・< 1 eV で最大と rms（meV） |
| `cmpse.py <res_A> <res_B>` | `SEXU`・`SECU` の列ごとの最大の差 |

kt1 の場所: `/mnt/data1/LiTi2O4_kbt_runs/bench_hgw666`（6³、`liti_mlo_v9` の反復 10 のコピー）、`bench_hgw999`（9³）。倍精度の基準は `res_fp64_oz`。
その前に `hgw --tetwt_write` で四面体の重みを書いておくと本番と同じ条件（無ければ `hgw` が自分で計算する）。

```bash
cd /mnt/data1/LiTi2O4_kbt_runs/bench_hgw666
BENCH_BIN=~/bin_dev BENCH_PREC="--use_fp32 --sigma_tf32" ./bench3.sh tf32test 1      # GPU 0 の 1 枚
python3 cmpef.py res_fp64_oz res_tf32test
```
