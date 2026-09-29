# AHC/Fe: bcc Fe の異常ホール伝導度

磁化を z 方向に向けた bcc Fe の異常ホール伝導度（AHC）$\sigma_{xy}$ を、フェルミエネルギーのずらし量 $\Delta E_F$ の関数として求める。
k 点メッシュを粗くした動作確認用のサンプルで、**値は収束していない**（§4）。

## 1. 計算する量

占有バンドのベリー曲率を k 点で和をとる。

$$ \sigma_{xy} = -\frac{e^2}{\hbar}\,\frac{1}{N_k V}\sum_{\mathbf{k}}\sum_{n}^{\rm occ}\Omega^{xy}_{n}(\mathbf{k}),\qquad
\Omega^{xy}_{n}(\mathbf{k}) = -2\,{\rm Im}\sum_{m}^{\rm unocc}
\langle \partial_{k_x}u_{n\mathbf{k}}|u_{m\mathbf{k}}\rangle\langle u_{m\mathbf{k}}|\partial_{k_y}u_{n\mathbf{k}}\rangle \tag{1}$$

$V$ は単位胞の体積、$N_k$ は k 点の数。k 微分は、隣の k 点との重なり $\langle u_{n\mathbf{k}}|u_{m\mathbf{k}+\mathbf{b}}\rangle$
（`huumat_MPI --ahc` が書く `UUU`, `UUD`）の差分でつくる。占有と非占有の組の重みは、四面体法（`ahc_tet.*`）と
k 点ごとの単純な和（`ahc_sp.*`）の二通りで出る。符号と規格化は `SRC/subroutines/x0kf_ahc.f90` の `sigAHC` ブロックにある。

## 2. 入力の要点（`ctrlg.fe.toml`）

| 表 1 | キー | 値 | 意味 |
|---|---|---|---|
| | `symgrp` | `"r4z"` | 磁化が z を向くので、対称操作を z 軸まわりの 4 回回転だけにする |
| | `[ham] so` | 2 | スピン軌道相互作用のうち $L_zS_z$ だけ（スピンは対角のまま、上向きと下向きの二つのチャンネル） |
| | `[ham] phispinsym` | true | 動径関数を上向きと下向きで共通にする |
| | `[bz] nkabc` | [8, 8, 8] | SCF の k 点メッシュ |
| | `[gw] n1n2n3` | [4, 4, 4] | AHC の式 (1) の k 点メッシュ |
| | `[gw] MagAtom` | 1 | 磁性サイト。`hahc` は一つ以上を要する |
| | `[gw] wan_out_emin`, `wan_out_emax` | -60, 180 (eV) | `hmaxloc` が b ベクトル（`BBVEC`）を作るときに読む窓。無いと止まる |

## 3. 実行

```bash
lmfa fe > llmfa
mpirun -np 4 lmf fe > llmf
job_band fe -np 4 --NoGnuplot > ljob_band     # job_AHC は bnd001.spin1 があることを確かめる
job_AHC fe -np 4 > ljob_AHC                   # ΔE_F = 0 の AHC まで
python3 `which hx0ahc.py` -1. 1. 11 > lhx0ahc # ΔE_F = -1 ... 1 eV の 11 点
python3 ahc_table.py fe                       # ahc.txt, scf.txt
```

- `job_AHC` は `lmf --jobgw=1 --wrhomt`, `hbasfp0 --job=8`, `hmaxloc`, `huumat_MPI --ahc`, `hahc --job=202 --ahc --interbandonly` を順に走らせる。
- `hx0ahc.py <始点> <終点> <点数> [hahc へのオプション]` は、`hahc` に `--EfermiShifteV=<値>` を付けて点の数だけ走らせ、
  `ahc_tet.isp11.dat` などを $\Delta E_F$ の順に並べ直す。`mpirun -np 4 python3 ...` とすれば点をランクに分けるが、
  これは mpi4py が `mpirun` と同じ MPI で作られているときに限る。違う MPI だと全ランクが全部の点を計算して出力が壊れるので、
  そのときは上のように `mpirun` なしで走らせる（`hahc` は `mpirun -np 1` で順に起動される。mpi4py は要らない）。
- AHC の k 点メッシュを変えるには、`job_AHC` と `hx0ahc.py` の両方に同じオプションを付ける。

```bash
job_AHC fe -np 4 '--ctrlg:gw.n1n2n3=[6,6,6]' > ljob_AHC
python3 `which hx0ahc.py` -1. 1. 11 '--ctrlg:gw.n1n2n3=[6,6,6]' > lhx0ahc
```

## 4. 結果

出力 `ahc_tet.isp11.dat`（スピン 1）, `ahc_tet.isp22.dat`（スピン 2）の列は、1: $\Delta E_F$ (eV)、2〜10: $\sigma_{xx},\sigma_{xy},\sigma_{xz},
\sigma_{yx},\dots,\sigma_{zz}$（$\Omega^{-1}{\rm cm}^{-1}$）。`ahc_table.py` が $\sigma_{xy}$（第 3 列）を取り出して `ahc.txt` に並べる（表 2）。

表 2: $\sigma_{xy}$（$\Omega^{-1}{\rm cm}^{-1}$）。AHC のメッシュ 4×4×4、gfortran、2026-09-30。

| $\Delta E_F$ (eV) | 四面体 スピン 1 | 四面体 スピン 2 | 四面体 和 | 単純和 |
|---:|---:|---:|---:|---:|
| -1.0 |  52.29 | 266.43 | 318.72 |   35.92 |
| -0.8 |  32.39 | 259.30 | 291.70 |   40.44 |
| -0.6 |  31.93 | 260.15 | 292.08 |  371.39 |
| -0.4 |  31.93 | 260.15 | 292.08 |  371.39 |
| -0.2 |  31.93 | 253.88 | 285.81 |  374.57 |
|  0.0 |  22.81 | 243.35 | 266.15 |  294.10 |
|  0.2 |  21.96 | 242.50 | 264.46 |    7.13 |
|  0.4 | 156.67 | 377.21 | 533.87 |    0.09 |
|  0.6 | 156.67 | 300.99 | 457.66 |   -7.84 |
|  0.8 | 156.67 | 300.99 | 457.66 |   -3.12 |
|  1.0 | 156.67 | 300.89 | 457.56 | -165.00 |

- SCF の磁気モーメントは 2.2755 $\mu_B$（`scf.txt`）。
- メッシュ 4×4×4（64 点）では $\sigma_{xy}$ は $\Delta E_F$ の階段関数になり、値も収束していない。同じポテンシャルで
  メッシュを 6×6×6 にすると $\Delta E_F=0$ の四面体の和は 362.3 になる。ベリー曲率はフェルミ面の近くの狭い領域で大きいので、
  文献の計算値（bcc Fe で 750 $\Omega^{-1}{\rm cm}^{-1}$ 程度）と比べるには、ずっと細かいメッシュが要る。
- 単純和は k 点ごとに占有を 0 か 1 で数えるので、粗いメッシュでは四面体法よりも大きく揺れる。
- 粗いメッシュでは、値がビルド（コンパイラと LAPACK）によっても変わる（表 4）。隣の k 点とのバンドの対応付けが、縮退した準位の固有ベクトルの
  取り方で変わるためで、メッシュを細かくすると差は縮む。

表 4: $\Delta E_F=0$ の $\sigma_{xy}$（四面体法、$\Omega^{-1}{\rm cm}^{-1}$）とビルド。同じ入力、2026-09-30。ifx は 4×4×4 で gfortran と同じ値。

| AHC のメッシュ | スピン 1: gfortran | スピン 1: nvfortran | スピン 2: gfortran | スピン 2: nvfortran |
|---|---:|---:|---:|---:|
| 4×4×4 | 22.79 | 19.63 | 243.32 | 287.87 |
| 6×6×6 | 125.69 | 137.15 | 236.59 | 246.62 |
| 8×8×8 | 885.36 | 884.84 | 1131.06 | 1148.30 |

## 5. スピン軌道相互作用の全部を入れる（`so = 1`）

$\mathbf{L}\cdot\mathbf{S}$ の全部を入れるには、§3 の計算のあと、同じディレクトリで `--ctrlg:ham.so=1` を付けて続ける。
`job_AHC` の中の lmf が `so = 1` で自己無撞着計算をやり直す（§3 の `rst.fe` から始まるので数反復）。

```bash
job_band fe -np 4 --NoGnuplot --ctrlg:ham.so=1 > ljob_band
job_AHC fe -np 4 --ctrlg:ham.so=1 > ljob_AHC
python3 `which hx0ahc.py` -1. 1. 11 --ctrlg:ham.so=1 > lhx0ahc
python3 ahc_table.py fe
```

`so = 1` ではスピンが混じるのでチャンネルは一つになり、出力は `ahc_tet.isp11.dat`, `ahc_sp.isp11.dat` だけになる。
4 コアで `job_AHC` が約 50 秒、`hx0ahc.py` が約 25 秒。この部分はテストに入っていない。

表 3: `so = 1` の $\sigma_{xy}$（$\Omega^{-1}{\rm cm}^{-1}$、四面体法）。AHC のメッシュ 4×4×4、gfortran、2026-09-30。磁気モーメントは 2.2820 $\mu_B$。

| $\Delta E_F$ (eV) | -1.0 | -0.8 | -0.6 | -0.4 | -0.2 | 0.0 | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| $\sigma_{xy}$ | 271.3 | 309.5 | 327.5 | 327.5 | 327.8 | 283.4 | 283.4 | 312.5 | 443.2 | 443.2 | 443.2 |

## 6. 実行時間とテスト

4 コアで 2〜3 分（SCF 1 分、`job_band` 20 秒、`job_AHC` 50 秒、`hx0ahc.py` 35 秒）。

```bash
cd Samples/AHC
testecalj Fe -np 4
```

テストは `scf.txt`（磁気モーメントと全エネルギー、許容差 1e-3）と `ahc.txt`（表 2 の全部の数、許容差 60 $\Omega^{-1}{\rm cm}^{-1}$。表 4 のビルドによる差を含む幅）を比べる。

## 7. 元のサンプル

`Samples/Legacy/AHC/Fe`（`ctrl.Fe` + `GWinput`、2025-09）を `ctrlg.fe.toml` に直し、SCF のメッシュを 12 から 8 に減らしたもの。
