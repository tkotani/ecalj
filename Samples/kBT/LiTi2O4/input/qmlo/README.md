# LiTi₂O₄ の MLO-QSGW の入力（2026-09-25〜）

kt1 の `qmlo_k6_*`（6³）と `qmlo_k9_*`（9³）の計算は全部この `ctrlg.liti2o4.toml` から始めている
（kt1 では `/mnt/data1/LiTi2O4_kbt_runs/liti_src_full9`、9³ は `liti_src_full9_k9`）。

- 6³: このファイルのまま。**9³ は `nkabc`、`n1n2n3`、`mlo_nkabc` の 3 行を `[9, 9, 9]` にする**（違いはそれだけ）
- 主な設定: `pwmode = 11`、`mixbeta = 0.5`、`t_tetrakbt = 0`（χ0 は T=0）、`SmearX0 = 0.0057` Ha（Im χ0 の Gaussian、FD 1000 K の幅）、
  `t_sigmaw = 1000`、`wcsmear`（既定）、`deltaq_scale = 0.1`、`mlo_method = 4`、
  MLO は 126 軌道（14 原子 × (s+p+d)）

## 回し方（LDA から `gwsc 10`）

```bash
mkdir my_input && cp ctrlg.liti2o4.toml my_input/
# env.sh: ecalj が動くシェルの PATH と LD_LIBRARY_PATH（ジョブの中で source される）
{ echo "export PATH='$PATH'"; echo "export LD_LIBRARY_PATH='$LD_LIBRARY_PATH'"; } > my_input/env.sh
PREC=tf32 RUNS_DIR=/path/to/runs bash ../../run_gwsc10.sh <tag> $PWD/my_input <bindir>     # GPUS=0,1 で 2 枚
```

- GPU は既定で GPU 0 の 1 枚（`-np 60 -np2 1`）。2 枚は `GPUS=0,1`（kt1 の GPU 1 は user のジョブ用なので、使ってよいと言われたときだけ。
  ecaljdoc の ForDevelopers §7）
- かかる時間（kt1、2026-09-27 夜のコード、GPU 2 枚）: 6³ tf32 は 39 分、9³ tf32 は 3.5 時間（1 反復 21 分、うち `hgw` 19.4 分。研究ログ 2026-09-28 03:05）。
  1 枚ならおよそ 2 倍
- 結果: `<tag>/steps.log`（反復ごとの ehf、最後にバンドの描画の成否）、`llmf.<N>run`、`QPU.<N>run`、`bnd_mlo_final.dat`（MLO バンド）、
  `bndPMT_final.dat`（従来の sigm バンド）
- 比べる: `../../cmp_gwsc10.py <new> <old>`（反復ごとの ehf と QP、最後のバンド）、図は `../../plot_band_pair.py`

経過と結果は `../../../kBT_research.md`（2026-09-27 22:50、23:35 など）。

## 反復ごとの確認（`check_iter.sh`）と 9³ の正常値

研究ログ 2026-09-26 の表から写した（9³ は 2026-09-26 の `gwsc 1` を 10 回の鎖で確かめた値。`run_gwsc10.sh` の `steps.log` の `mloON=` も同じ数になる）。

**9³ での正常値**（6³ と違うので注意。`check_iter.sh` は `CHECK_MLOON=35 CHECK_GWDRV=50 CHECK_MAXSECS=12000` で使う）

どちらの印も `sigmlo_init` が**プロセスごとに最初の 1 回**出す行で、`getsenex` を呼ぶのは担当する k 点を持つプロセスだけ。
`lmf` は k 点を大きさ ⌈nk/np⌉ のかたまりで配るので、印の数 = ⌈nk / ⌈nk/np⌉⌉ になる。

| 印 | `lmf` の np | k 点数 6³ → 9³ | かたまり | 印の数 6³ → 9³ |
|---|---|---|---|---|
| `mloON`（SCF）| 60 | 既約 16 → 35（`BZMESH: ngrp nq 48 ...`）| 1 → 1 | **16 → 35** |
| `gwdrv`（GW ドライバ）| 60 | `QPLIST.jobgw1` 72 → 200（`ZmloSig` のレコード数と同じ）| 2 → 4 | **36 → 50** |
| MLO バンドの `writeham` | 8 | 16 → 35 | 2 → 5 | **8 → 7**（9³ は 5×7=35 で 8 番目が空）|

意味のある異常はこの数より**少ない**場合（一部のプロセスが MLO 経路に入らず従来の `sigm` に落ちた）。
ほか: 1 反復 ≈ 5300 秒（6³ は ≈ 950。いずれも 2026-09-26 の fp32 のコード。09-27 の tf32 では 6³ が約 220 秒）、スナップショット 2.4 GB（6³ は 0.9）。
