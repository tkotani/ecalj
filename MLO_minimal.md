# 課題 (1) — $E_F$ を含む領域を再現するミニマル模型

ecalj branch `mlo3` / 2026-09-13

MLO には性格の違う 2 つの課題がある。

1. **$E_F$ を含むエネルギー領域を再現するミニマル模型。** 半導体・金属・
   Al₂O₃:Cr(ギャップ中の Cr d + バンド端)はすべてこれ。**本書の対象。**
2. **4f(や 3d)だけを取り出す。** これは別課題で、比較的やさしい。

FeMgO は空格子球が要るので保留(§5)。

---

## 1. 重み行列 — 4 つの method は同じ一本の式

$$
|F^{\mathrm{MLO}}_{\mathbf{k}n}\rangle
=\sum_{n'n''}|\Psi^{\mathrm{PMT}}_{\mathbf{k}n'}\rangle\;
A^{\mathbf{k}}_{n'n''}\;
\langle\Psi^{\mathrm{MTO}}_{\mathbf{k}n''}|F^{\mathrm{MTO}}_{\mathbf{k}n}\rangle,
\qquad
A^{\mathbf{k}}_{n'n''}=S^{\mathbf{k}}_{n'n''}\,\theta^{\mathbf{k}}_{n'n''}
$$

$$
S^{\mathbf{k}}_{n'n''}=\langle\Psi^{\mathrm{PMT}}_{\mathbf{k}n'}|\Psi^{\mathrm{MTO}}_{\mathbf{k}n''}\rangle
$$

実装 (`m_hreduction.f90`) はどの method でも

$$
\theta^{\mathbf{k}}_{n'n''}=
\max\Bigl[\;
\sigma\Bigl(\tfrac{\varepsilon_{n'}-\varepsilon_{\mathrm{frz}}}{w_{\mathrm{frz}}}\Bigr),\;
\sigma\Bigl(\tfrac{\varepsilon_{n'}-e^{\mathrm{cut}}_{n''}}{w}\Bigr)
\Bigr],
\qquad \sigma(x)=\frac{1}{1+e^{x}}
$$

を組み立てる。method 0/1/2/4 では $\varepsilon_{\mathrm{frz}}=-\infty$ なので
第 1 項は恒等的に 0 で、**第 2 項だけ**が働く。違いは $e^{\mathrm{cut}}_{n''}$ だけ:

| method | $e^{\mathrm{cut}}_{n''}$ | `mlo_emax` | 性格 |
|---|---|---|---|
| **2** | $\varepsilon^{\mathrm{MTO}}_{n''}$ | 無視 | 自己参照のみ。調整パラメータなし |
| **1** | $e_{\max}$ | 唯一のカット | 大域カットのみ |
| **0** | $\max(e_{\max},\varepsilon^{\mathrm{MTO}}_{n''})$ | **床** | 大域の床 + 自己参照 |
| **4** | $\max(E_F+\Delta,\varepsilon^{\mathrm{MTO}}_{n''})$ | 無視 | 床を自動化した method 0 |
| **3** | $\max(\varepsilon_{\mathrm{frz}},\varepsilon^{\mathrm{MTO}}_{n''})$、第 1 項も有効 | 無視 | 二段シグモイド |

$$
\varepsilon_{\mathrm{frz}}=\varepsilon_{n_{\mathrm{occ}}+1}(\mathbf{k})+\Delta
\quad(\text{method 3 のみ、}\mathbf{k}\text{ ごとの最低非占有準位})
$$

### 二つの項の役割は違う

- **自己参照項** $\varepsilon^{\mathrm{MTO}}_{n''}$ は $\mathrm{rank}(A)=n_{\mathrm{dimMTO}}$ と
  各軌道の帯域保持を担う。局在バンドはこれだけで済む(課題 2)。
- **床**(大域カット)は $n''$ に依らず、**窓の中の全状態の重みを 1 にする**。
  「$E_F$ の少し上まではきっちり合わせたい」という系全体の要求で、
  これが課題 1 の本体である。

---

## 2. 原本サンプルは何をしていたか

MLOsamples の 19 サンプルは**すべて method 0**、コマンドラインの上書きも無し。
物質ごとの調整は `mlo_emax` 一本で、3 つのレジームに分かれていた。

| `mlo_emax` | 物質 | 床の位置 |
|---|---|---|
| **0** | Si, GaAs, GaAsSoc, C, C.sp | $E_F$ |
| **5〜7 eV** | Al₂O₃:Cr (7), FeCo (7), FeMgO (5), FeMgOSoc (5) | $E_F+5\ldots7$ eV |
| **未設定** | Cu, Fe, FeSoc, GdCo5, GdION, NiO, RuO₂, SmP, SrTiO₃ | auto $=\varepsilon_{n_{\mathrm{dimMTO}}}$ |

### auto は原理のある窓ではなかった

床 $\varepsilon_{n_{\mathrm{dimMTO}}}-E_F$ を実測すると:

| 模型 | 系 | auto 床 |
|---|---|---:|
| ミニマル | Cu (d のみ, 5) | **−2.15 eV** |
| | SmP (4f, 7) | **−6.56 eV** |
| | SrTiO₃ (15) | **−3.35 eV** |
| | GdION (4f, 7) | −0.03 eV |
| ちょうど良い | Fe (spd, 9) | +2.74 eV |
| | RuO₂ (26) | +3.58 eV |
| 大きい | Si (19) | **+22.8 eV** |
| | FeMgO (85) | +17.6 eV |
| | C (19) | **+63.2 eV** |

小さい模型では床が $E_F$ より下に落ちて**何もしていない**(だから Cu・SmP・GdION では
method 0 と method 2 の数値が一致する)。大きい模型では自由電子連続帯の中に飛び、
やはり窓にならない — 手置き `emax` が必要だった列がこれである。
Fe と RuO₂ でたまたま有用な位置に来たのが「auto で動いていた」正体で、
**原理のある窓ではなく本数勘定の偶然**だった。

---

## 3. 固定設定 — method 0, `mlo_emax = 2.449` eV

床を $E_F+\Delta$($\Delta=0.18$ Ry $=2.449$ eV)と置く。物質ごとの入力はゼロ。
評価窓は $E_F\pm1.5$ eV(課題 1 の目的そのもの)と VBM−2 eV 〜 CBM+2 eV。

| 系 (MTO 数) | $E_F\pm1.5$ | 窓全体 | $\Delta v/v$ | ギャップ誤差 | $m^*$ 比 |
|---|---:|---:|---:|---:|---:|
| C (sp のみ, 8) | 0.0072 | 0.0189 | 0.011 | −29 | 0.69 |
| Cu (spd, 9) | 0.0498 | 0.0274 | 0.054 | −77 | 0.98 |
| Fe (spd, 9) | 0.0163 | 0.0166 | 0.046 | −7 | — |
| SrTiO₃ (15) | 0.0054 | 0.0124 | 0.009 | +17 | 1.07 |
| NiO (17) | 0.0211 | 0.0207 | 0.034 | +14 | 1.07 |
| Si (19) | 0.0180 | 0.0211 | 0.019 | +23 | 1.16 |
| C (spd, 19) | 0.0101 | 0.0179 | 0.009 | −6 | 1.00 |
| GaAs (19) | 0.0123 | 0.0173 | 0.011 | +26 | 1.00 |
| FeCo (20) | 0.0064 | 0.0059 | 0.013 | −3 | — |
| RuO₂ (26) | 0.0091 | 0.0089 | 0.030 | +0 | — |
| Al₂O₃:Cr (53) | 0.0052 | 0.0350 | 0.017 | **−128** | **0.58** |
| **平均** | **0.0146** | **0.0184** | **0.023** | | |

ΔE rms は eV。$m^*$ 比は $m^*_{\rm MLO}/m^*_{\rm DFT}$。金属に $m^*$ は無い。

### 手調整版・method 3 との比較($E_F\pm1.5$ eV の ΔE rms)

| 床の取り方 | 平均 | 最大 |
|---|---:|---:|
| **$E_F+\Delta$(定数、全物質共通)** | **0.0111** | **0.0211** |
| CBM$+\Delta$(定数、全物質共通) | 0.0110 | 0.0222 |
| $\varepsilon_{n_{\mathrm{occ}}+1}(\mathbf{k})+\Delta$($\mathbf{k}$ 依存) | 0.0240 | 0.0749 |
| method 3(二段シグモイド) | 0.0184 | 0.0325 |
| 原本(物質ごと手調整) | 0.0206 | 0.0556 |

(Cu を除く 10 系。Cu spd は別枠)

**固定設定 1 つが、物質ごとの手調整より良い。**

---

## 4. 設計上の 3 つの判断

### (a) シグモイドは 1 つ。$\varepsilon^{\mathrm{MTO}}_{n''}$ はカット位置に入れる

method 3 は 2 つのシグモイドを $\max$ で結ぶが、
$\varepsilon_{\mathrm{frz}}\le e^{\mathrm{cut}}_{n''}$ が常に成り立つので、
$w_{\mathrm{frz}}=w$ とすると**第 1 項が常に勝ち、$\varepsilon^{\mathrm{MTO}}_{n''}$ が
完全に消える**(Fe: 0.031 → 0.075 eV)。
$\varepsilon^{\mathrm{MTO}}_{n''}$ は裾を通してしか入らず、しかも $w>w_{\mathrm{frz}}$ のときだけ。
method 0/4 のように**カット位置そのもの**に入れる方が強い
(同じ床で Fe: 0.017 対 0.031 eV)。

### (b) 床は $\mathbf{k}$ に依らせない

$\varepsilon_{n_{\mathrm{occ}}+1}(\mathbf{k})$ を使うと窓の意味が $\mathbf{k}$ ごとに変わる。
金属では最低非占有準位の振れが大きく、Fe で 0.075 eV と大きく悪化する。

### (c) ギャップの広い系では床を CBM 基準に

$E_F+\Delta$ と CBM$+\Delta$ はほとんどの系で同じだが、**Al₂O₃:Cr だけ決定的に違う**:

| 系 | 床 $=E_F+\Delta$ | 床 $=$ CBM$+\Delta$ |
|---|---|---|
| Al₂O₃:Cr ギャップ誤差 | **−128 meV** | **+7 meV** |
| Al₂O₃:Cr $m^*$ 比 | 0.58 | 1.08 |
| 他 10 系 $E_F\pm1.5$ 平均 | 0.0111 | 0.0110 |

Al₂O₃:Cr は CBM が $E_F+6.23$ eV にあり、$E_F+2.449$ eV の床では CBM が拘束されない。
原本が `emax = 7 eV` としていた理由そのものである。

**したがって課題 1 の最終形は**

$$
\boxed{\;
\theta^{\mathbf{k}}_{n'n''}=\sigma\Bigl(\frac{\varepsilon_{n'}-e^{\mathrm{cut}}_{n''}}{w}\Bigr),
\qquad
e^{\mathrm{cut}}_{n''}=\max\bigl(\varepsilon_{\mathrm{CBM}}+\Delta,\;\varepsilon^{\mathrm{MTO}}_{n''}\bigr)
\;}
$$

$\varepsilon_{\mathrm{CBM}}$ は**大域の**伝導帯下端(金属では $E_F$)。
パラメータは $\Delta=0.18$ Ry、$w=0.20$ Ry の 2 つだけで、物質ごとの入力は無い。

> 現在の `mlo_method = 4` は床が $E_F+\Delta$ で、CBM 基準になっていない。
> 大域 CBM は $\mathbf{k}$ 全点を見ないと決まらないので、`m_HamPMT` 側で
> 先に求めて渡す実装が要る。**未実装。**

---

## 5. 残る問題

- **Al₂O₃:Cr のギャップ** — 上記 (c) の実装待ち。
- **C.sp(s,p のみ 8 軌道)の $m^*$ 比 0.69、ギャップ −29 meV。**
  ミニマル模型ほど自己参照項の裾($w$)に依存し、$w$ を狭めると痩せる
  (C.sp は $w=0.12$ で損失が 6 倍悪化、ふつうの半導体は $w=0.12$ が最良)。
  $w$ の最適値が模型の大きさで割れており、全系共通なら $w\simeq0.15$ が折衷点。
- **Cu (spd) のギャップ誤差 −77 meV**、$E_F\pm1.5$ で 0.050 eV と本表の最大。未精査。
- **FeMgO。** 二乗誤差の 99% を準位の 10% が担い、その全部が $E_F+0.7$〜$1.0$ eV に集中。
  Fe/MgO 積層で Worb は各原子 15 チャネル(ほぼフル MTO)だが**空格子球が無い**。
  障壁領域に基底関数が一つも無いので界面状態は原理的に表現できない。空格子球が要る。
- **Si の $6^3$ の $m^*$ は信用できない。** $8^3$ で取り直すと損失最小の設定が
  ギャップも $m^*$ も最良になる(ギャップ −3 meV、$m^*$ 0.98)。

---

## 6. コードの訂正 — レジーム則は一度も走っていなかった

`m_hreduction.f90` の method 3 の枝は、`regime7` ブロックで軌道ごとの幅
`ewuse` を計算した直後に

```fortran
              endblock regime7
            endif
            ewuse = eww          ! ← 計算結果を無条件に破棄
```

と上書きしていた。**`mlo_ewalpha` は一度も効いていない。**
以前「MTO ごとに幅を変えてもバンドがバイト単位で同一だった」ので
レジーム則は無効と報告したが、それはこの死んだコードのせいで、
**レジーム則は検証されたのではなく一度も試されていない**。

この死んだコードと、互いに矛盾する v5/v7/v8/v9 のコメント地層、
未使用の `itgt` / `regime_report` / `mlo_tau` / `mlo_ewalpha` / `mlo_ewmin` を削除した。
削除後に全 11 系のバンドがバイト単位で一致することを確認済み。
