# Samples/Magnon — マグノンのスペクトル

強磁性体の横スピン帯磁率 χ^{+-}(q, ω) を、MLO の模型の上で作った相互作用 W から計算し、マグノンの分散を出す（`job_mlo_magnon`）。
2026-10-02 までは最局在 Wannier 関数による版（`job_magnon`、試料 `Fe_magnon`・`Ni_magnon`・`FeCo_magnon`・`Fe_bcc_in_sc_magnon`）もあった。
それは git のタグ `last-wannier` にある（`MD/past_log.md` 表 1）。両者の比較は `MD/wannier_vs_mlo.md` と ecaljdoc mlo §6。

**表 1**. サンプル（`testecalj <名前> -np 8`）

| ディレクトリ | 中身 |
| --- | --- |
| [`Fe_mlo_magnon`](Fe_mlo_magnon/README.md) | bcc Fe。Wannier 版との比較の表と図（Wannier 版の結果 `wannier_TrRpm.syml001` を添えてある） |

## 流れ

```bash
lmfa fe
mpirun -np 8 lmf fe                    # LDA
job_band fe -np 8 --NoGnuplot          # バンド（MLO の模型を合わせる相手）
job_mlo_magnon fe -np 8                # MLO の模型、MLO 基底の W、k と k+q の重なり、χ^{+-}。MagSuscep.syml<n> を書く
```

- 入力は `ctrlg.<sname>.toml`（`[mlo] mlo_lm` と `mlo_nkabc` が要る）と経路 `syml.<sname>`
- `MagSuscep.syml<n>`: 経路 n の上の K(q, ω) と R(q, ω)。Im R の山がマグノン

## 試験が比べるもの

- `Fe_mlo_magnon`: `MagSuscep.syml001` の値を 10 % の幅で、Im R の山の位置を 5 % の幅で比べる。参照は 2026-09-30 に手元（gfortran、`-np 8`）で作ったもの
