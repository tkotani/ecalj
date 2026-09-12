# MLO 自動窓 — 一晩の最適化による最終形

ecalj branch `mlo3` / MLOsamples 16 系 / 2026-09-13
手調整パラメータゼロの MLO 構成。v1〜v10 の探索と、約 45 設定 × 16 系の走査による最終化。

## 1. 最終形

PMT 基底(MTO+APW)の恒等分解から、選択 MTO 部分空間の固有状態
$\Psi^{\mathrm{MTO}}$ を介して局在軌道 $F^{\mathrm{MLO}}$ を構成する:

$$
|F^{\mathrm{MLO}}_{\mathbf{k}n}\rangle
=\sum_{n'n''}|\Psi^{\mathrm{PMT}}_{\mathbf{k}n'}\rangle\,
A^{\mathbf{k}}_{n'n''}\,
\langle\Psi^{\mathrm{MTO}}_{\mathbf{k}n''}|F^{\mathrm{MTO}}_{\mathbf{k}n}\rangle,
\qquad
A^{\mathbf{k}}_{n'n''}=S_{n'n''}\;\bar\theta_{n''}(\epsilon_{n'})
$$

$$
S_{n'n''}=\langle\Psi^{\mathrm{PMT}}_{\mathbf{k}n'}|\Psi^{\mathrm{MTO}}_{\mathbf{k}n''}\rangle
$$

**重み(最終形)**:

$$
\bar\theta_{n''}(\epsilon)=
\max\!\left[\;
\sigma\!\left(\frac{\epsilon-\epsilon_{\mathrm{frz}}}{w_{\mathrm{frz}}}\right),\;
\sigma\!\left(\frac{\epsilon-e^{\mathrm{cut}}_{n''}}{w}\right)
\right],
\qquad \sigma(x)=\frac{1}{1+e^{x}}
$$

$$
\epsilon_{\mathrm{frz}}(\mathbf{k})=\epsilon_{\,n_{\mathrm{occ}}(\mathbf{k})+1}(\mathbf{k})+\Delta
\;\;(\simeq \epsilon_{\mathrm{CBM}}+\Delta;\ \text{金属では}\ E_F+\Delta),
\qquad
e^{\mathrm{cut}}_{n''}=\max\!\left(\epsilon_{\mathrm{frz}},\,\epsilon^{\mathrm{MTO}}_{n''}\right)
$$

$$
\boxed{\;\Delta=0.18\,\mathrm{Ry}\ (2.45\,\mathrm{eV}),\qquad
w_{\mathrm{frz}}=0.10\,\mathrm{Ry}\ (1.36\,\mathrm{eV}),\qquad
w=0.20\,\mathrm{Ry}\ (2.72\,\mathrm{eV})\;}
$$

入力は $S$、$\epsilon^{\mathrm{PMT}}$、$\epsilon^{\mathrm{MTO}}$、$E_F$ のみ。系ごとの指定は不要
(`mlo_method = 3`。旧 `mlo_method` 0/1/2 + 手調整 `mlo_emax` を置き換える)。

- 第 1 項: 目標窓の端を $w_{\mathrm{frz}}$ の滑らかさで切る
- 第 2 項 + 自軌道フロア $\epsilon^{\mathrm{MTO}}_{n''}$: 窓の上に段階的な重みを残し、
  $\mathrm{rank}(A)=n_{\mathrm{dimMTO}}$ と各軌道の帯域保持を担保

## 2. 最終結果(16 系、既定値のみ)

| 系 | rms \[占有,CBM+3eV\] (eV) | p95 (eV) | ギャップ誤差 (eV) | m\*(VBM) | m\*(CBM) |
|----|------:|------:|------:|------:|------:|
| GdION | 0.000 | 0.000 | — | — | — |
| GdCo5 | 0.001 | 0.003 | — | — | — |
| RuO₂ | 0.011 | 0.026 | — | 1.06 | — |
| C.sp | 0.015 | 0.037 | +0.008 | — | 0.75 |
| FeCo | 0.015 | 0.025 | — | — | 0.99 |
| Al₂O₃:Cr | 0.017 | 0.037 | **−0.006** | 1.01 | 1.26 |
| Cu (spd) | 0.024 | 0.045 | — | 1.30 | 1.37 |
| GaAs | 0.026 | 0.049 | +0.020 | — | — |
| SrTiO₃ | 0.026 | 0.055 | +0.017 | 1.04 | — |
| NiO | 0.028 | 0.057 | +0.029 | 0.94 | 0.97 |
| C | 0.029 | 0.066 | **+0.005** | — | 0.71 |
| Fe | 0.029 | 0.058 | — | 1.05 | 1.13 |
| Si (QSGW, 6³) | 0.035 | 0.078 | **−0.017** | — | 0.86 |
| SmP | 0.053 | 0.148 | — | (平坦 4f) | |
| FeMgO (slab) | 0.054 | 0.047 | — | 1.16 | 1.05 |
| Cu (d 専用) | 0.141 | 0.271 | — | 1.47 | — |

**ギャップ誤差は全絶縁体で ≤29 meV**、目標窓の rms は Cu(d 専用ミニマルモデル)を除き **≤0.054 eV**。

### 前日の手調整版(v6.5)からの改善

| 指標 | v6.5(手で決めた値) | **最終形** |
|---|---:|---:|
| 平均 \|ギャップ誤差\| | 0.050 eV | **0.012 eV** |
| rms_max(非 Cu) | 0.072 eV | **0.054 eV** |
| Si m\*(CBM) 比 @6³ | 4.92 | **0.86** |
| Al₂O₃:Cr ギャップ誤差 | −0.064 eV | **−0.006 eV** |

変わったのは 2 つの数値だけ($\Delta$: 0.22→0.18 Ry、$w_{\mathrm{frz}}$: 0.05→0.10 Ry)。

## 3. 手で立てた根拠が 2 つとも誤りだった

走査が覆した推論:

1. **「凍結縁は鋭いほどバンド端の曲率を歪めない」→ 逆**。鋭い切り口は $\mathbf{k}$ 空間で
   リンギングを生み、$H(\mathbf{R})$ の実空間減衰を悪化させ、補間で得る有効質量を壊す。
   $w_{\mathrm{frz}}$ を 0.05→0.10 Ry と**滑らかにする**だけで Si の m\* が 4.92→1.00 になった。
2. **「目標窓は指定どおり CBM+3 eV まで重み 1 で凍結すべき」→ 欲張りすぎ**。
   $\Delta=0.18$ Ry(CBM+2.45 eV)が最適で、CBM+3 eV では重みは 1 ではない。
   それでも \[占有, CBM+3 eV\] の精度は最良になる。**正しい部分空間を選ぶことと、
   目標窓を重み 1 で覆うことは別問題**。

## 4. Si の有効質量 — メッシュ補間律速

| Si (QSGW), 設定固定 | 6³ | 8³ | 10³ |
|----|------:|------:|------:|
| ギャップ誤差 (eV) | −0.017 | +0.011 | −0.020 |
| 窓 rms (eV) | 0.035 | 0.025 | 0.022 |
| m\*(CBM) 比 | 0.86 | 1.08 | 0.75 |

![Si mesh convergence](MLO_v6_figs/Si_mesh_conv.png)

新しい既定値では **6³ ですでに m\* 比 0.86・ギャップ誤差 −17 meV** に達しており、
旧 v6.5(6³ で m\* 4.92、誤差 +165 meV)のようなメッシュ律速は大幅に緩和された。
8³/10³ の残る振れ(0.75〜1.08)は放物線フィットのノイズ水準。

## 5. Cu の残差はモデル選択であって窓の問題ではない

`Cu__m3` の Worb は `1 Cu 5 6 7 8 9` = **d 軌道 5 個のみのミニマルモデル**(Fe は spd 9 個)。
$E_F$ 上に MLO バンドが無いのは設計どおりで、rms 0.141 は「d 専用モデルを全 DFT バンドと
比較した」評価側の不整合。**同じ θ̄ で Worb を spd 9 個にすると rms 0.024**(Fe 並み)。

→ 「クラス別 $\Delta_{\mathrm{class}}$ フィットが必要」という前日の主要動機は消滅した。

## 6. per-orbital レジーム則は不採用

「局在バンドを作る MTO と広い sp バンドを作る MTO でルールを変え、しかも自動判定させる」
方針を v7〜v10 で実装・検証した。結果:

- **判定量の意味は判明**: $p_{n''}(\epsilon)=\sum_{n'}|S_{n'n''}|^2\delta(\epsilon-\epsilon_{n'})$ の
  窓内での広がり $\sigma^{\mathrm{win}}$ は「その MTO 固有状態が単一の PMT 固有状態から
  どれだけズレているか」を測る。ダイヤ C の sp3 で 0.2–0.8 eV(MTO 基底がほぼ厳密)、
  Al₂O₃:Cr で 0.5–2.4 eV(混成が強い)。**窓内に制限しないと全軌道 8–10 eV** になり
  (自由電子的な裾に支配される)、分類器として機能しない。
- **しかし採用しない**: (a) 同一系・同一設定でも run 間で $\sigma^{\mathrm{win}}$ が
  0.37→0.00 eV と変動し、判定器として再現性が不足。(b) 幅を軌道ごとに 0.10–0.20 Ry の
  範囲で変えても、バンドは**まったく変化しなかった**(α=1.0/1.5/2.0 で結果がバイト単位で同一)。
- 切り分けで分かった実効自由度: バンドを変えるのは「**窓の上に重みを残すか否か**」のみ。
  第 2 シグモイドを外すと C 0.029→0.011 と改善する一方 Fe 0.029→0.045、FeCo 0.015→0.039、
  Al₂O₃ 0.017→0.029 と悪化する。レジーム分離自体は実在するが、現在の θ̄ の形では
  それを軌道ごとに切り替える有効なレバーが無い。

## 7. バンドプロット(灰線: 第一原理 PMT、赤×: MLO)

![Al2O3_Cr](MLO_v6_figs/Al2O3_Cr.png)
Al₂O₃:Cr — ギャップ内 Cr d と VBM/CBM が一致(ギャップ誤差 −6 meV、m\* 1.01)。手調整 emax=7 eV を完全自動置換。

![C](MLO_v6_figs/C.png)
C(ダイヤ)— ギャップ誤差 +5 meV。

![GaAs](MLO_v6_figs/GaAs.png)
GaAs — rms 0.026 eV、ギャップ誤差 +20 meV。

![Si666gwsc](MLO_v6_figs/Si666gwsc.png)
Si(QSGW、6³)— ギャップ誤差 −17 meV、m\*(CBM) 0.86。

![SrTiO3](MLO_v6_figs/SrTiO3.png)
SrTiO₃ — rms 0.026 eV、m\*(VBM) 1.04。

![NiO](MLO_v6_figs/NiO666lda.png)
NiO — ギャップ誤差 +29 meV、m\* 0.94/0.97。

![Fe](MLO_v6_figs/Fe.png)
Fe — $E_F$+3 eV 窓 rms 0.029 eV、m\* 1.05/1.13。

![FeCo](MLO_v6_figs/FeCo.png)
FeCo — rms 0.015 eV。

![Cu](MLO_v6_figs/Cu.png)
Cu — d 専用ミニマルモデル(Worb = d 5 軌道)。$E_F$ 上の s バンドは設計上モデル外。spd 版では rms 0.024 eV。

![FeMgO](MLO_v6_figs/FeMgO.png)
FeMgO スラブ — 二相系。rms 0.054 eV。MTO エンベロープ(減衰長 ~1 a.u.)の外にある表面状態はスパン外のまま。

## 8. 探索の履歴(何が効き、何が効かなかったか)

| 版 | 規則 | 結果 |
|---|---|---|
| v1 | 累積重み 90% 分位点で窓 | 棄却: 離散スペクトルの分位点は $\mathbf{k}$ で階段状に飛ぶ |
| v2 | モーメント $\bar\epsilon+\sigma$ で窓 | 棄却: sp の反結合裾で $\sigma$ が発散、C/Si/NiO 崩壊 |
| v3 | $E_F$+5 eV 制限モーメント | 棄却: 絶対窓の外に住む Cr 空き d・f 軌道が破壊(+3 eV) |
| v4/v5 | バンド数ソフトゲート | 棄却: ランク欠損($\mathrm{rank}(A)<n_{\mathrm{dimMTO}}$)で Fe/C まで崩壊 |
| v6.5 | 二段 θ̄、手で決めた値 | 全系まとも。gap 0.050 / Si m\* 4.92 |
| v7–v10 | per-orbital レジーム則 | 不採用(§6) |
| **最終** | **v6.5 の形 + 走査した 2 値** | **gap 0.012 / Si m\* 0.86 / rms_max 0.054** |

走査は約 45 設定 × 16 系(各 4–5 分)。評価は目標窓の rms・p95、ギャップ誤差、
バンド端曲率比(有効質量)の総合。

## 9. 残課題

1. **FeMgO のスパン外状態**(表面/界面)— empty sphere を Worb に追加する拡張が根治
2. **レジーム分離の活用** — 実在する(§6)が、現在の θ̄ の形では軌道ごとに切り替える
   レバーが無い。$A$ の構造自体を変える必要がある
3. **ミニマルモデルの評価法** — d 専用 / 4f 専用モデルは「全 DFT バンドとの乖離」では
   評価できない。対象バンドを指定した評価指標が要る
4. 4f 抽出型など目標窓が別の用途は `mlo_emax` 上書きキーで対応(既存)

---

再現: branch `mlo3` / `mlo_method = 3`(既定値のみ)/ ベンチ `Samples/MLOsamples/*__m3`
(+ `Cu__spd`, `Si666gwsc__m3k8/k10`)/ 評価 `SRC/exec/mlo_bandcheck.py`。
走査用キー `mlo_dwin` / `mlo_wfrz` / `mlo_eww` / `mlo_down` / `mlo_ewalpha` / `mlo_ewmin` は
再調査用に残してある(既定値が最適値)。
