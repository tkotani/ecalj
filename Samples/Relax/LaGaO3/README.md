# Relax/LaGaO3: 原子位置の緩和（LaGaO3、20 原子）

斜方晶ペロブスカイト LaGaO3（単位胞に 20 原子）の原子位置を、力を使って緩和する。`lmf` は二重のループで動く:
内側が電子状態の自己無撞着計算、外側が原子を動かすループである。`[dyn]` を書くと力の計算は自動で入る。
格子定数と単位胞の形は動かさない。

## 入力の要点（`ctrlg.lagao3.toml`）

- `[dyn] mode = 5`（Fletcher-Powell）、`nit = 20`（原子を動かす回数の上限）、`xtol = 0.0001`、`step = 0.015`。
- `[iter] nit = 20`, `conv = 1e-5`: 一つの配置での自己無撞着計算。
- `[bz] nkabc = [3, 2, 3]`, `pwmode = 11`, `xcfun = 1`。この系では 1 反復の時間は k 点の数にほとんどよらない（表 3）。

## 実行

```bash
lmfa lagao3 > llmfa
mpirun -np 4 lmf lagao3 > llmf                    # 緩和を最後まで行う（長い）
python3 relax_result.py lagao3                     # relax_energy.txt, relax_force.txt
```

途中で止めるときや試すときは、コマンド行で回数を絞る。表 1 は

```bash
mpirun -np 4 lmf lagao3 --ctrlg:iter.conv=1e-4 --ctrlg:iter.convc=1e-4 --ctrlg:dyn.nit=2 > llmf
```

の結果である。試験は `Samples/Relax` で `testecalj LaGaO3 -np 4`（表 2 の短い計算）。
手で実行するときは、このディレクトリを別の場所に写してその中で行う（`testecalj` はこのディレクトリをまるごと写して使う）。

## 結果の見方

- `save.lagao3`: 反復ごとに 1 行。`c`（収束）か `x`（`iter.nit` に達した）の行が一つの配置の終わり。
  緩和が収束すると最後に `C` の行が付く。
- `llmf` の `Forces, with eigenvalue correction` の表と `Maximum Harris force`: 力（mRy/bohr）。
- `AtomPos.lagao3`: 配置ごとの原子位置（alat 単位、`itrlx` が配置の番号）。**このファイルがあると、`lmf` は入力の
  位置の代わりにここから位置を読む。** 初めからやり直すときは `AtomPos.lagao3` と `rst.lagao3` を消す。
  続きを計算するときはそのまま `lmf` を実行する。
- `relax_result.py` は配置ごとの終わりのエネルギーと力を `relax_energy.txt`, `relax_force.txt` にまとめる。

表 1: `conv = 1e-4` で収束させた二つの配置（gfortran、MPI 4 並列、2026-09-30）。配置 1 は入力の位置。

| 配置 | 反復 | ehf (eV) | 最大の力 (mRy/bohr) | 力が最大の原子 |
|---|---|---|---|---|
| 1 | 10 | -1159548.603892 | 17.76 | 11 (O) |
| 2 | 10 | -1159548.668130 | 12.34 | 16 (O) |

- 配置 1 から 2 へエネルギーは 64.2 meV（4.72 mRy）下がる。二つの配置の力の平均と変位から見積もった変化
  $-\sum_i \frac{1}{2}(\mathbf F_i^{(1)}+\mathbf F_i^{(2)})\cdot\Delta\mathbf R_i$ は −4.90 mRy で、力とエネルギーは整合している。
- Ga（原子 1〜4）は対称性で力が 0。変位の最大は 0.155 bohr。

表 2: 試験の短い計算 `--ctrlg:iter.nit=5 --ctrlg:dyn.nit=2`（各配置 5 反復で打ち切り、収束していない）。

| 配置 | 反復 | ehf (eV) | 最大の力 (mRy/bohr) |
|---|---|---|---|
| 1 | 5 | -1159548.659520 | 15.23 |
| 2 | 5 | -1159548.673665 | 15.18 |

表 3: 1 反復の時間（MPI 4 並列、他のジョブと同時）。

| nkabc | 既約 k 点 | 3 反復の時間 (秒) |
|---|---|---|
| 2×1×2 | 4 | 88 |
| 1×1×1 | 1 | 97 |

## 計算時間

MPI 4 並列（16 コアの計算機）で 1 反復が 25〜50 秒（他のジョブの量による）。表 1 の計算は 18 分、
試験（表 2、10 反復）は 6.3 分。緩和を最後まで行うと、旧サンプルの記録（`Samples/Legacy/LaGaO3_relax/xxx/save.lagao3~`）では配置が 12 個、反復が 181 回である。

## 試験の内容（`test.py`）

- `relax_energy.txt`: 二つの配置の終わりの ehf, ehk（許容 2e-3 eV）
- `relax_force.txt`: 二つの配置の最大の力と各原子の力（許容 0.05 mRy/bohr）
- `AtomPos.lagao3` を参照 `AtomPos.lagao3.ref` と比べる（許容 2e-5）

## 由来

`Samples/Legacy/LaGaO3_relax` を今の入力形式と命令に直したもの。入力の数値は変えていない。
