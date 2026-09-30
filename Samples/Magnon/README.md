# Samples/Magnon — マグノンのスペクトル

強磁性体の横スピン帯磁率 χ^{+-}(q, ω) を、局在した基底の上で作った相互作用 W から計算し、マグノンの分散を出す。
基底は最局在 Wannier 関数（`job_magnon`）か MLO（`job_mlo_magnon`。設定が少なく、Wannier の窓や反復が要らない）。

**表 1**. サンプル（どれも `testecalj <名前> -np 8`）

| ディレクトリ | 中身 |
| --- | --- |
| `Fe_magnon` | bcc Fe |
| `Ni_magnon` | fcc Ni（経路 2 本） |
| `FeCo_magnon` | B2 FeCo（2 原子） |
| `Fe_bcc_in_sc_magnon` | bcc Fe を 2 原子の単純立方のセルで計算したもの |
| [`Fe_mlo_magnon`](Fe_mlo_magnon/README.md) | bcc Fe を MLO の模型で（`job_mlo_magnon`）。Wannier 版との比較の表と図がある |

Wannier の 4 つ合わせて 8 分（kr7、`-np 8`）、MLO 版は 2 分。作業ディレクトリは Wannier 版で 1 つ 1〜5 GB になる。

## 流れ

```bash
lmfa fe
mpirun -np 8 lmf fe                    # LDA
job_band fe -np 8 --NoGnuplot          # バンド
job_magnon -np 8 fe                    # Wannier 関数、W、χ^{+-}。TrKpm.syml001 と TrRpm.syml001 を書く
gnuplot mag3d.glt                      # 図（gnuplot が無ければ飛ばしてよい）
```

- `TrKpm.syml<n>`: 経路 n の上の K(q, ω) の対角和。`TrRpm.syml<n>`: R(q, ω) の対角和で、Im の山がマグノン
- 入力は `ctrlg.<sname>.toml` と経路 `syml.<sname>`
- MLO 版は `job_magnon` の代わりに `job_mlo_magnon fe -np 8`（入力に `[mlo] mlo_nkabc` が要る）。K と R は `MagSuscep.syml<n>` に出る

## 試験が比べるもの

- `Fe_magnon`、`Ni_magnon`、`FeCo_magnon`: `TrKpm` と `TrRpm` を点ごとに比べる（絶対・相対とも 1e-3、ω = 0 の行は除く）
- `Fe_bcc_in_sc_magnon`: K は 10 % の幅で、R はマグノンの山の位置（各 q で |Im Tr R| が最大になる ω）を 10 % の幅で比べる。
  2 原子のセルではバンドが折り返されて縮退し、K と R がコンパイラとランクの数で変わるため（K で最大 7 %、山の位置で最大 6 %。
  kt1 と kr7 の nvfortran の `-np 8` どうしは一致する）

- `Fe_mlo_magnon`: `MagSuscep.syml001` の値を 10 % の幅で、Im R の山の位置を 5 % の幅で比べる

参照は Wannier の 4 つが 2026-08-19 に kt1（nvfortran、`-np 60`）、`Fe_mlo_magnon` が 2026-09-30 に手元（gfortran、`-np 8`）で作ったもの。

図の見本は `eps/`、`Ni_magnon/ref/`。
