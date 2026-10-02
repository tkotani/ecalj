# contour_test — Σc の contour 分解（虚軸積分 + 実軸極項）と準位 smearing の単体試験

> 説明の本体: ecaljdoc の [manual/kBT.md](../../../ecaljdoc/manual/kBT.md)「3.6 Σc の contour 分解と準位 smearing の整合」（サイト https://ecalj.github.io/ecaljdoc/manual/kBT）。サンプルの一覧は [manual/samples.md](../../../ecaljdoc/manual/samples.md)。

2026-09-22 01:30。1 中間準位 ε′、モデル $W_c(z) = -v\,\omega_p^2/(\omega_p^2 - z^2 - i\gamma z)$（Drude–Lorentz、上半面で解析的、
$W_c(i\omega') = -v\omega_p^2/(\omega_p^2+\omega'^2+\gamma\omega')$ は $\omega'$ の 1 次項を持つ = 金属型）、行列要素 = 1、
$\omega_p$ = 1.8 eV。実コードの経路そのもの — `m_wfac::pole_weights`（窓重み Φ、wcsmear の bin 重み）と
m_sxcf_sc の虚軸重み（FD 40 点平均、niw GL 点、$u_a$ 引き算、段差の解析処理。ブロックをそのまま複写）— で
$\Sigma_c(\omega) = I + P$ を作り、厳密求積（$I_{ex}$ は $\omega'=|w_e|\tan t$ の置換で 20000 点、$P_{ex}$ は FD 4000 点）と比較。
$\omega$ を $\varepsilon'\pm2$ eV で掃引（34 meV 刻み）。$\varepsilon'-E_F$ = +0.27 eV（非占有、t2g 想定）。単位は $W_c(0)=-1$。

```
gfortran -O2 -ffree-line-length-none -o contour_test stub_cmdopt.f90 ../../../SRC/subroutines/wfacx.f90 ../../../SRC/subroutines/mate.f90 contour_test.f90
./contour_test [niw] [T_K] [gamma_eV] [wcsmear 0/1]   > out_niw10_T1000_g0.1_wcs1.txt
```

## 結果

| case | max\|dSum\| for \|ω−ε′\|<0.5 eV | max\|dSum\| at the pole (ω−ε′ = 1.3–2.1 eV) | 2nd diff of Sum_code near ε′ | same, reference |
|---|---|---|---|---|
| niw10_T1000_g0.02_wcs1 | 0.00275 | 1.170 | 0.0009 | 0.0037 |
| niw10_T1000_g0.1_wcs0 | 0.01021 | 6.812 | 0.0008 | 0.0036 |
| **niw10_T1000_g0.1_wcs1**（標準） | **0.00236** | 0.051 | 0.0008 | 0.0036 |
| niw10_T1000_g0.5_wcs1 | 0.00350 | 0.006 | 0.0006 | 0.0034 |
| niw10_T3000_g0.1_wcs1 | 0.00184 | 0.075 | 0.0005 | 0.0040 |
| niw10_T300_g0.1_wcs1 | 0.00182 | 0.493 | 0.0009 | 0.0029 |

（dSum = コード − 厳密。2nd diff = 34 meV 刻みの 2 階差分の最大、滑らかさの指標。）

1. **ω = ε′ の跨ぎ（中間準位が最終エネルギーと縮退する場所）**: I と P はそれぞれ $\pm\tfrac12 W_c(0)$ を kBT 幅で入れ替え、
   和は連続（2 階差分 ≤ 0.001、参照と同じ）。コードと厳密求積の差は ≤ 0.0035（$W_c(0)$ の 0.35 %）。ω = E_F の跨ぎ
   （窓の反転）も滑らか。T = 300 / 1000 / 3000 K、γ = 0.02 / 0.1 / 0.5 eV で同じ。**contour 分解と FD smearing の実装は正しく、
   niw = 10 の虚軸求積で十分**（niw を増やす試験は不要）。
2. **W_c の極を踏む場所（ω − ε′ ≈ ω_p）**: wcsmear なし（3 点内挿）はコードが厳密から最大 6.8（$W_c(0)$ の 7 倍）ずれる —
   針の機構そのもの。wcsmear ありで 0.05（γ = 0.1）、極が bin より鋭い γ = 0.02 で 1.2、核が狭い 300 K で 0.5。
   つまり wcsmear の残差は「極の幅 vs bin 幅・核幅」で決まり、極そのものを抜くか（chi0_filterw）、極を鈍らせる（温度、SmearX0）以外に
   下げる手は無い。

結論: 反復で残る E_F 直上の凸凹は、1 準位・1 行の $\Sigma(\omega)$ の contour 実装（縮退跨ぎ）では説明できない。
