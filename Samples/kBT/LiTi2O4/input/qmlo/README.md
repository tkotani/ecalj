# LiTi₂O₄ の MLO-QSGW の入力（2026-09-25〜）

kt1 の `qmlo_k6_*`（6³）と `qmlo_k9_*`（9³）の計算は全部この `ctrlg.liti2o4.toml` から始めている
（kt1 では `/mnt/data1/LiTi2O4_kbt_runs/liti_src_full9`、9³ は `liti_src_full9_k9`）。

- 6³: このファイルのまま。**9³ は `nkabc`、`n1n2n3`、`mlo_nkabc` の 3 行を `[9, 9, 9]` にする**（違いはそれだけ）
- 主な設定: `pwmode = 11`、`mixbeta = 0.5`、`t_tetrakbt = 0`（χ0 は T=0）、`t_sigmaw = 1000`、`mlo_method = 4`、
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
- かかる時間（kt1、2026-09-27 のコード、GPU 2 枚）: 6³ tf32 は 39 分、9³ tf32 は約 6〜7 時間。1 枚ならおよそ 2 倍
- 結果: `<tag>/steps.log`（反復ごとの ehf、最後にバンドの描画の成否）、`llmf.<N>run`、`QPU.<N>run`、`bnd_mlo_final.dat`（MLO バンド）、
  `bndPMT_final.dat`（従来の sigm バンド）
- 比べる: `../../cmp_gwsc10.py <new> <old>`（反復ごとの ehf と QP、最後のバンド）、図は `../../plot_tf32_fp32.py`

経過と結果は `../../../kBT_research.md`（2026-09-27 22:50、23:35 など）。
