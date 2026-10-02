# kBT — 有限温度と均し方 (`t_tetrakbt` / `t_sigmaw` / `wcsmear`)

GW と QSGW で、電子温度と数値的な均し方を決める `[gw]` のキーを説明する。$W$ 側($\chi_0$)と $\Sigma$ 側($G$)に別のキーがある(表 4)。
$\Sigma$ 側の `t_sigmaw` は常に有効な数値的な幅、$\chi_0$ 側の `t_tetrakbt` は必須で、正なら物理的な電子温度、
負なら $\mathrm{Im}\chi_0$ を均す Gaussian の幅(温度で書く)、0 なら均さない。

この頁にはいまの仕様を書く。古い入力や記録を読むときのキーの変遷は §2.5 の表 2。2026-06〜09 の経緯と当時の計算の説明は
[kBT の記録](./kBT_history)。

**表 4**. `[gw]` のキー(温度の単位は K)

| キー | どこに効くか | 既定 |
|---|---|---|
| `t_tetrakbt`(必須) | $\chi_0$ の均し方を符号で選ぶ(0: 均さない、$T>0$: 有限温度のテトラヘドロン、$-T$: Gaussian。§2 の表 1) | 無し(書かないと止まる) |
| `t_sigmaw` | $\Sigma_x=Gv$、$\Sigma_c=G(W-v)$ の中間準位を均す Fermi-Dirac 核の幅(§3) | 1000 |
| `wcsmear` | $\Sigma_c$ の実軸極項で、極の位置にも同じ核を使って $W_c$ を積分する(§3.5) | true |
| `chi0_skip_window` = [emin, emax] (eV) | 両端がこの窓に入るバンド対を $\chi_0$ から除く(cRPA 風。試験中) | off |

gwinit のテンプレート(`ctrlgenToml.py` もこれを使う)は `t_tetrakbt = 300`、`t_sigmaw = 300` を書く。

**使い分け**

- 絶縁体・半導体: テンプレートのままでよい(§4)。
- 電子温度の効果を見る: `t_tetrakbt = t_sigmaw = T`(§5.1 の Si、§5.2 の Fe)。
- 金属で反復や k メッシュによってバンドが荒れる($W$ のプラズモン極、§3.5): `wcsmear`(既定)と `t_sigmaw = 1000`。
  $\chi_0$ を T=0 のまま均すなら `t_tetrakbt` を負にする。

**LiTi₂O₄ の設定**: $\chi_0$ は T=0 のテトラヘドロンで、$\mathrm{Im}\chi_0$ を Gaussian で均す `t_tetrakbt = -992.4`
(Fermi-Dirac 992.4 K と同じ幅、標準偏差 0.155 eV)、`t_sigmaw = 1000`、`wcsmear = true`、`mixbeta = 0.5`。
$6^3$ はこれで 8 反復で収束し(反復 7→8 の QP の変化は最大 10、平均 2 meV)、$6^3$ と $9^3$ の一発の差は `dSEnoZ` 最大 13 meV。
MLO-QSGW の計算(§5.3、ecalj `Samples/kBT/LiTi2O4`)も同じ設定である。

---

## 0. T=0 と有限温度で何が変わるか(一覧)

温度は $\chi_0$ 側と $\Sigma$ 側で独立のキー。`t_tetrakbt > 0` を使うなら $\Sigma$ 側の `t_sigmaw` と同じ T にするのが自然
(必須ではない)。LiTi₂O₄ の設定は $\chi_0$ を T=0 のまま `t_tetrakbt < 0` の Gaussian で均す(冒頭の「LiTi₂O₄ の設定」)。

### $\chi_0$ 側(hx0fp0 / hgw の W-build)

`t_tetrakbt` の 3 通り(0、$T>0$、$-T$)の違いは §2 の表 1 と式 (3)〜(8)。どの場合も `chi0_skip_window = [emin, emax]`
(eV、両端が窓に入るバンド対を $\chi_0$ から除く。試験中)が使える。

### $\Sigma$ 側(hsfp0_sc / hgw の $\Sigma_c$、hsfp0 の $\Sigma_x$)

$\Sigma$ 側に「T=0 と有限温度」の区別は無い。核は常に Fermi-Dirac、幅は `t_sigmaw`。
唯一 $\chi_0$ 側の温度に連動するのは Fermi 準位である。

**表 5**. $\Sigma$ 側で $\chi_0$ の温度に連動するもの、しないもの

| | `t_tetrakbt` ≤ 0 | `t_tetrakbt > 0` |
|---|---|---|
| Fermi 準位 | k メッシュの上の準位を、幅 $k_BT$(`t_sigmaw`)の Gaussian で均して数えて決める(`efsimplef2ax`。$\chi_0$ が使う `EFERMI` とは別に決まる) | `EFERMI_kbt`(無ければ止まる) |
| 準位の占有核($\Sigma_x$ の重み、極項の窓重み) | Fermi-Dirac、幅 $k_BT$ = `t_sigmaw` | 同左 |
| 極項で $W_c$ を拾う位置 | `wcsmear = true`(既定): FD 核で $W_c$ を積分。`false`: 平均 $\bar\omega$(3 点内挿) | 同左 |
| 中間状態の候補窓・$\omega$ メッシュの余裕 | $E_F \pm 15\,k_BT$(`sig_window`) | 同左 |
| 虚軸積分、$W_c$ 自体、sigm の混合(`mixbeta`) | 共通 | 共通 |

- `t_sigmaw` は数値的な幅であって、完全な有限温度自己エネルギー(Bose 項、§8.7)ではない。
  核 $-\partial f/\partial\varepsilon$ の標準偏差は $\pi k_BT/\sqrt3 \simeq 1.81\,k_BT$(1000 K で 0.156 eV)。
- `wcsmear` は温度スイッチではなく「極の位置にも占有核と同じ分布を使う」スイッチ(§3.5)。
- $\chi_0$ 側は T=0 で本当に鋭い $\theta$ なので、`t_tetrakbt` を 0 でなくして初めて幅が付く。

---

## 1. なぜ要るか

用途は二つある。

1. **電子温度の物理**。$T>0$ では Fermi 準位の近くの占有が Fermi-Dirac 分布になり、半導体ではギャップを越えて熱励起された
   キャリアが遮蔽に加わる。$\chi_0$ 側の `t_tetrakbt = T`(§2.3)と $\Sigma$ 側の `t_sigmaw = T`(§3)で扱う。
   Si では電子温度 3000 K でギャップが 0.16 eV 縮む(§5.1)。
2. **金属の QSGW の数値的な安定化**。金属では小さな $q$ の $W_c(\omega)$ に鋭いプラズモン極があり、$\Sigma_c$ の実軸極項が
   それを 1 点で拾うと、k メッシュと反復のたびに値が大きく変わる(§3.5)。これを均すのが $\Sigma$ 側の `wcsmear` と `t_sigmaw`、
   $\chi_0$ 側の `t_tetrakbt < 0`(Gaussian、§2.4)である。このときの温度は数値的な幅で、物理の温度ではない。

---

## 2. $\chi_0$ 側 — `t_tetrakbt`

`[gw] t_tetrakbt` (K) は**必須**である(書かないと止まる)。符号で $\chi_0$ の均し方を選ぶ(表 1)。
gwinit のテンプレートは `t_tetrakbt = 300` を書く。キーの名前と意味は何度か変わったので、古い入力や記録を読むときは §2.5 の表 2 を見よ。

**表 1**. `t_tetrakbt` の 3 通り

| `t_tetrakbt` | 占有の差 $f_a-f_b$ | $\mathrm{Im}\chi_0$ の $\omega$ 方向 | $\chi_0$ の Fermi 準位($\Sigma$ は §0 の表 5) | 使いどころ |
|---|---|---|---|---|
| 0 | 鋭い階段、式 (3) | そのまま | `EFERMI` | 絶縁体・半導体 |
| $T>0$ | 温度 $T$ の Fermi-Dirac、式 (4)(5) | そのまま | `EFERMI_kbt`、式 (6) | 電子温度の物理(§1、§5.1) |
| $-T$ | 鋭い階段、式 (3) | Gaussian で均す、式 (7)(8) | `EFERMI` | 金属で $W$ のプラズモン極が荒れるとき(§3.5) |

```toml
[gw]
t_tetrakbt = 0        # (K) 必須。0: 均さない。T > 0: 温度 T の有限温度テトラヘドロン。-T: Im chi0 を Gaussian で均す
```

### 2.1 共通 — $\mathrm{Im}\chi_0$ のヒストグラムと Hilbert 変換

どの場合も、$\chi_0$ はまず $\mathrm{Im}\chi_0$ を $\omega>0$ の**ヒストグラム**として作り、そこから Hilbert 変換で実軸の実部と
虚軸の値を作る。ビン $j$ は $[\omega_j,\omega_{j+1}]$ で、中心を $\bar\omega_j$、幅を $\Delta_j=\omega_{j+1}-\omega_j$ と書く
(`frhis`。`HistBin_dw` と `HistBin_ratio` で決まる指数メッシュで、$\omega$ が大きいほど広い)。$\mathbf k$ の準位 $e_a$ から
$\mathbf k+\mathbf q$ の準位 $e_b$ への遷移をテトラヘドロン法で積んだビンの重みは

$$
h_j \propto \sum_{\mathbf k}\sum_{a,b}\int_{\omega_j}^{\omega_{j+1}}\!d\omega\;
\bigl|M_{ab}\bigr|^2\,\bigl(f_a-f_b\bigr)\,\delta\bigl(\omega-(e_b-e_a)\bigr),
\qquad f_a=f(e_a)
\tag{1}
$$

である($M_{ab}$ は積基底との行列要素。係数とスピンの和は省いた)。テトラヘドロン法は四面体の中で準位と行列要素を
線形に補間し、その範囲で $\delta$ 関数のビンへの積分を厳密に行う。$h_j/\Delta_j$ がビン $j$ での $\mathrm{Im}\chi_0$ の平均になる。
$\mathrm{Im}\chi_0$ は $\omega$ の奇関数($\mathrm{Im}\chi_0(-\omega)=-\mathrm{Im}\chi_0(\omega)$)なので $\omega>0$ だけを持ち、
ビンの平均を中心の間で区分線形につないで

$$
\chi_0(z)=\frac1\pi\int_{-\infty}^{\infty}\!d\omega'\,\frac{\mathrm{Im}\chi_0(\omega')}{\omega'-z}
\tag{2}
$$

をビンごとに解析的に積分する(`hilbertmat`、`dpsion5.f90` の `dpsion_init` と `dpsion_chiq_*`)。実軸では $z=\omega+i0$ の
実部が $\mathrm{Re}\chi_0(\omega)$ で、虚部はビンの平均をそのまま使う。虚軸では $z=i\nu$ で、奇関数なので値は実数になる。
3 通りの違いは、式 (1) の占有の差 $f_a-f_b$ と、式 (2) に渡す前の $h_j$ の扱いだけである。

### 2.2 `t_tetrakbt = 0` — T=0 のテトラヘドロン

占有は Fermi 準位 $E_F$(`EFERMI`)での階段で、式 (1) の占有の差は

$$
f_a-f_b\;\to\;\theta(E_F-e_a)\,\theta(e_b-E_F)
\tag{3}
$$

である(`lindtet6`。占有された $a$ から空の $b$ への遷移で、$\omega=e_b-e_a>0$ に入る)。

### 2.3 `t_tetrakbt = T > 0` — 有限温度のテトラヘドロン(method B′)

有限温度の占有の差は、式 (3) の階段を Fermi 準位 $E$ について熱核 $-f'$ で平均したものに等しい。$f(x)=1/(e^{x/k_BT}+1)$ として

$$
\int dE\,\bigl(-f'(E-\mu)\bigr)\,\theta(E-e_a)\,\theta(e_b-E)=f(e_a-\mu)-f(e_b-\mu)\qquad(e_a<e_b)
\tag{4}
$$

式 (1) は $f_a-f_b$ について線形なので、有限温度のビンの重みは、Fermi 準位を $E$ に置いた T=0 の重み $h^{(0)}_j(E)$ を
$E$ について平均したものになる。$E=\mu+2k_BT\,t$ と置くと $(-f')\,dE=\tfrac12\,\mathrm{sech}^2t\,dt$ で、これを 20 点で積分する:

$$
h_j(T)=\sum_{i=1}^{20}c_i\,h^{(0)}_j\bigl(\mu+2k_BT\,t_i\bigr),\qquad
c_i=\frac{w_i\,\tfrac12\mathrm{sech}^2t_i}{\sum_{i'}w_{i'}\,\tfrac12\mathrm{sech}^2t_{i'}}
\tag{5}
$$

節点 $t_i$ と重み $w_i$ は、$t\in[-6,6]$(すなわち $|E-\mu|\le12\,k_BT$)を $[-6,-1],[-1,0],[0,1],[1,6]$ の 4 区間に分けた
5 点 Gauss-Legendre(`gausq_fd`)。$c_i$ の和を 1 にしてあるので、区間の外を切っても占有の和は保たれる。
四面体の 4 頂点で $a$ も $b$ も窓 $|e-\mu|\le12\,k_BT$ の片側の外にそろっている(窓の中で占有が平ら)ときは、T=0 の値
そのものなので 20 回の計算を省く(`lindtet6_kbt`、`tetwt5.f90`)。式 (5) は $\mathbf k$ 積分と $E$ 積分の順を入れ替えただけで、
誤差は $E$ の数値積分だけである(§7.1、§7.2)。得られる分子は Adler–Wiser の $f_a-f_b$ であって $f_a(1-f_b)$ ではない
(応答関数として正しいのはこちら。$\Sigma$ 側で落ちるものは §7.4)。

$\mu$ = `EFERMI_kbt` は `heftet` が、テトラヘドロンの状態数 $N(E)$ を同じ熱核で均して

$$
\sum_{i=1}^{20}c_i\,N\bigl(\mu+2k_BT\,t_i\bigr)=N_{\rm val}
\tag{6}
$$

を二分法で解き、ファイル `EFERMI_kbt` に書く(`fermi_kbt_tetra`、`main_heftet.f90`)。$T\to0$ でテトラヘドロンの $E_F$ に戻る。
ギャップの中で解が区間になる(左辺がギャップにわたって $N_{\rm val}$ で平ら)ときは区間の中央を取る(§9 の 2)。
`t_tetrakbt > 0` のときは $\chi_0$ も $\Sigma$ もこれを使い、ファイルが無ければ止まる(§0、§3)。`t_tetrakbt` ≤ 0 のときは、$\chi_0$ は `heftet` のテトラヘドロン法の `EFERMI` を、$\Sigma$ は k メッシュの上で数えた Fermi 準位(§0 の表 5)を使う。

出力での確認(`heftet` の標準出力):

```
tetrakbt_init: T[K], kbt[Ry], kbt[eV]  2000  0.12667E-01 0.17234E+00
```

### 2.4 `t_tetrakbt = -T` — $\mathrm{Im}\chi_0$ を Gaussian で均す

占有と Fermi 準位は T=0 のまま(式 (3)、`EFERMI`)で、式 (2) に渡す前にビンの重み $h_j$ を $\omega$ 方向に Gaussian で均す。
Gaussian の標準偏差は、温度 $T$ の Fermi-Dirac の核 $-f'$ の標準偏差に合わせる:

$$
\sigma=\Bigl[\int dx\,x^2\bigl(-f'(x)\bigr)\Bigr]^{1/2}=\frac{\pi k_BT}{\sqrt3}\simeq1.81\,k_BT
\tag{7}
$$

($T=1000$ K で $\sigma=0.1563$ eV $=0.005744$ Ha)。均し方は

$$
\tilde h_i=\sum_j G_{ij}\,h_j,\qquad
G_{ij}=\frac{\bigl[K(\bar\omega_i-\bar\omega_j)-K(\bar\omega_i+\bar\omega_j)\bigr]\,\Delta_i}{\sum_{i'}K(\bar\omega_{i'}-\bar\omega_j)\,\Delta_{i'}},\qquad
K(x)=e^{-x^2/2\sigma^2}
\tag{8}
$$

で、和はどれも $\omega>0$ のビンを走る(`gaussianfilterhis` と `smearx0_apply`、`dpsion5.f90`)。列 $j$ はビン $j$ の重みを
周りのビンへ Gaussian で配り直す:

- 分子の $\Delta_i$(行き先のビンの幅): 幅の違う指数メッシュの上でも、平均 $h/\Delta$ を Gaussian で畳み込んだことになる
- 分子の第 2 項 $-K(\bar\omega_i+\bar\omega_j)$: 奇関数への拡張。負の $\omega$ にこぼれる分を符号を変えて折り返す。
  $\tilde h$ は $\omega\to0$ で 0 に行く
- 分母: 対称な核の $\omega>0$ での和。$\bar\omega_j\gg\sigma$ では $\sqrt{2\pi}\,\sigma$ になり、式 (8) は連続の式

$$
\widetilde{\mathrm{Im}\chi_0}(\omega)=\int_0^\infty\!d\omega'\,\bigl[g_\sigma(\omega-\omega')-g_\sigma(\omega+\omega')\bigr]\,\mathrm{Im}\chi_0(\omega'),
\qquad g_\sigma(x)=\frac{e^{-x^2/2\sigma^2}}{\sqrt{2\pi}\,\sigma}
\tag{9}
$$

  すなわち奇関数に延ばした $\mathrm{Im}\chi_0$ と $g_\sigma$ の畳み込みになる。式 (9) は f-sum $\int_0^\infty\omega\,\mathrm{Im}\chi_0\,d\omega$ を保つ
  ($\omega>0$ での重みの和は、0 の近くで折り返す分だけ減る)。$\bar\omega_j\lesssim\sigma$ では Gaussian の一部しか $\omega>0$ に
  無いので分母が $\sqrt{2\pi}\,\sigma$ より小さく、式 (8) は式 (9) からずれる
- $\pm\omega$ を両方持つ場合(`npm = 2`、マグノンの $\chi^{+-}$ など): 折り返しの無い対称な核を $\omega>0$ と $\omega<0$ に別々にかける

均した $\tilde h_j$ から式 (2) で実部と虚軸の値を作るので、実軸の $W(\omega)$、虚軸の $W(i\nu)$、静的な $W(0)$ はどれも同じ
均した $\mathrm{Im}\chi_0$ から来る。

$-T$ は Gaussian の幅を温度で書く約束であって、有限温度の $\chi_0$(式 (5))とは別物である。式 (5) は準位の占有を均し、
式 (8) は遷移エネルギー $e_b-e_a$ の分布を均す。

入力の値は `m_GWinput` で $\sigma$(Ha。内部の変数名は `SmearX0`)に直され、出力に次の行が出る:

```
 SmearX0 (Ha) =   5.7438E-03  : Gaussian smearing of Im chi0 along omega ([gw] t_tetrakbt < 0)
```

### 2.5 キーの変遷

**表 2**. $\chi_0$ と $\Sigma$ の均し方の入力の変遷

| 期間 | $\chi_0$ の均し方 | $\Sigma$ の準位の幅(§3) |
|---|---|---|
| 〜2026-06 | `GaussianFilterX0`(GWinput、$\mathrm{Im}\chi_0$ の Gaussian) | `esmr`(Ry、Gaussian) |
| 2026-06 〜 09-18 | 有限温度は論理キー `tetrakbt = true` と `t_tetrakbt` (K)。Gaussian は `SmearX0`(Ha、$\sigma$。`GaussianFilterX0` の改名) | `esmr`、または `t_sigmakbt`(K、Fermi-Dirac) |
| 09-19 〜 09-27 | `t_tetrakbt > 0` だけで有限温度(既定 0 = 均さない)。`SmearX0` は 09-20 から `t_tetrakbt > 0` と排他 | 09-20 から `t_sigmaw`(K、Fermi-Dirac)。`esmr` と `t_sigmakbt` は換算して読む |
| 09-28 〜 | `t_tetrakbt` は必須で、符号で 3 通り(表 1)。`SmearX0` は読まない(書いてあれば止まる) | `t_sigmaw` |

- `SmearX0 = s` と `t_tetrakbt = -T` は、式 (7) で $s\,[\mathrm{Ha}]=\pi k_BT/\sqrt3$ のとき同じ計算になる(0.0057 Ha と $-992.4$ K)。
  `ctrlg_update.py ctrlg.<sname>.toml` がこの換算で書き換え、`t_tetrakbt` の無いファイルには 0(以前の既定)を足す
- `esmr` = $e$ (Ry、Gaussian の標準偏差)と `t_sigmaw` = $T$ は、$e=\pi k_BT/\sqrt3$ のとき同じ幅になる(0.003 Ry と 262 K、0.01 Ry と 872 K)。
  `t_sigmaw` が無い入力に `esmr` や `t_sigmakbt` があれば、`esmr` はこの換算で、`t_sigmakbt` は値のまま `t_sigmaw` として読む
- gwinit のテンプレートは 2026-06 から 09-16 まで `tetrakbt = true`、`t_tetrakbt = 300` を書いていた(09-17 に外した)。
  09-28 からは `t_tetrakbt = 300`、`t_sigmaw = 300` を書く

---

## 3. $\Sigma$ 側 — `t_sigmaw`

$\Sigma$ の中間準位($G$ の極)は Fermi-Dirac 核で均される。その幅が `t_sigmaw` (K) で、
$\Sigma_x = Gv$ の占有重み、$\Sigma_c = G(W-v)$ の実軸極項の窓重み、`wcsmear` の極位置の分布の
すべてに同じ核が使われる。

```toml
[gw]
t_sigmaw   = 1000.0   # (K) Sigma の準位を均す FD 核の幅。既定 1000
t_tetrakbt = 1000.0   # chi0 も有限温度にするなら。Sigma の Fermi 準位が EFERMI_kbt になる
```

窓の重みは Fermi–Dirac の累積分布

$$
\Phi(x)=\frac{1}{1+e^{-x/k_BT}},
\qquad
w(e_l,e_h;\varepsilon_k)=\Phi(e_h-\varepsilon_k)-\Phi(e_l-\varepsilon_k)
\tag{10}
$$

である。**同じ核が $\Sigma_c$ の虚軸積分にも入る**
— 極項と虚軸積分は一つの contour の二つの半分で、準位の smearing は両方に同じ核で入れないと
$\Sigma_c$ が壊れる。これは重要な点なので [§3.6](#_3-6-sigma-c-の-contour-分解と準位-smearing-の整合) に独立に書く。
中間状態の候補窓と実軸 $\omega$ メッシュの
余裕は核の裾 $E_F \pm 15\,k_BT$(`sig_window`)で、`t_sigmaw` を上げれば自動的に広がる(§7.3)。

### 温度の選び方

- 絶縁体: ギャップより十分小さければ効かない(Si で 262 K と 1000 K の差は数 meV)。
- 金属: 原理的には $n_1n_2n_3 \to \infty$、`t_sigmaw` → 0 の極限だが、k メッシュの
  Fermi 面での準位間隔と同程度の幅が要る。300〜1000 K を物質ごとに。LiTi₂O₄ $9^3$ は 1000 K。
- `t_tetrakbt` を使うときは同じ温度にするのが自然(必須ではない)。

### 出力での確認

```
    t_sigmaw =   1000.0 K  (Sigma kernel width kBT =    0.006333 Ry)
 sigmaw_setup: Sigma kernel Fermi-Dirac, t_sigmaw[K]=  1000.0  ef: counted on the k mesh (efsimplef2ax); chi0 at T=0
```

`t_tetrakbt > 0` なら `ef<-EFERMI_kbt=...`。このとき `EFERMI_kbt`(`heftet` が書く)が無ければ abort する。

---

## 3.5 金属の荒れの正体と `wcsmear`

金属の QSGW では、バンドが k メッシュと反復によって数十 meV 荒れたり、LDA からの一発 GW で $E_F$ の近くに eV の大きさの針が
立ったりすることがある(LiTi₂O₄ の 1000 K 以下、$9^3$ で O 2p 帯の底が 30〜100 meV 荒れ、一発 GW で −2〜−3 eV の針)。原因は次のとおり
(調べた記録は ecalj の `MD/research_log.md` の 2026-09-18〜19)。

- **犯人は $\Sigma_c$ の実軸極項**(切ると異常な非対角 ⟨48|Σc|33⟩ が −13 → −0.35 eV)。
  虚軸積分、sigm の混合、offset-Γ の頭、$\omega$ 内挿法は無関係。
- 極項が踏んでいるのは、第一殻 $q$ の $W_c(q,\omega)$ にある **t2g プラズモン極**
  ($\omega_p \simeq 1.78$ eV、幅 ≈ 0.1 eV。`--dumpW` で直接確認)。極項は各中間状態
  $n'$ について $W_c(\mathbf k,\, \varepsilon_{\mathbf q n} - \varepsilon_{\mathbf q-\mathbf k, n'})$ を
  1 点で拾うので、その $\varepsilon$ が極の ±0.1 eV に落ちる項が ±50 a.u. を受け取り、
  反復で固有値が動くたびに極の反対側へ渡る(シーソー)。
- 温度で単調に消える: ⟨48|Σc|33⟩ = −27 / −13 / −3.2 / −1.9 eV (500 / 1000 / 1500 / 2000 K)。
  2000 K 以上の run に事故が無かったのはこのため。

**対策 `wcsmear = true`**(既定): 中間準位の smearing 核 $g$(Fermi-Dirac の
$-\partial f/\partial\varepsilon$、幅 `t_sigmaw`)を、窓の重み $w_{n'} = \int g$ にだけ使い、$W_c$ は平均エネルギー $\bar\omega$ の
1 点で評価するのが PRB 76, 165106 Eq. 58 のやり方(`wcsmear = false`)である。$k$ 積分の離散化として重みと極の位置に同じ核を使うのが一貫しており、

$$
\Sigma^{\rm pole}_{nn''} = \sum_{\mathbf k n'} s\,M^{*}_{n'n}
\left[\int_{e_l}^{e_h} de\; g(e-\varepsilon')\,W_c(\mathbf k,\,\omega-e)\right] M_{n'n''}
\tag{11}
$$

と $W_c$ 側を積分する(式 (11))。$W_c$ が核の幅の中で滑らかなら 1 点評価と一致し(Si で最大 77 meV、
低い状態は数 meV)、極は核の幅で均される(表 6)。

**表 6**. `wcsmear` の効果

| LiTi₂O₄ 1000 K | `wcsmear = false` | `wcsmear = true` |
|---|---|---|
| $6^3$ 一発、Γ–L の針 | −3.4 eV | 消える |
| $6^3$ 一発、O 2p 荒れ | 12.2 meV | 7.2 meV |
| $9^3$ iter 2→3、O 2p 荒れ(混合なし) | 47.8 meV | 12.7 meV |
| 実軸極項のコスト | — | ×6.7(全体で $6^3$ +7 %、$9^3$ +26 %) |

```toml
[gw]
t_sigmaw   = 1000.0   # 既定値。wcsmear = true も既定
t_tetrakbt = 1000.0   # chi0 側も有限温度にするなら
```

注意: (1) 概念的には「G の極を幅 $k_BT$ のスペクトル関数として扱う」近似を極項全体に通したもので、
contour の縁(E_F と ω)は鋭いまま。(2) 有限温度で Fermi-Dirac の $-\partial f/\partial\varepsilon$ を
極の位置の分布に使うのは発見的な正則化(厳密な有限温度の式は準位を鋭いまま θ → f にするだけ)。
本来の幅は $k$ 積分(テトラヘドロン)から出るべきもの。(3) 実装は QSGW 経路(`hsfp0_sc` / `hgw`)
と、一発 GW の `hsfp0`(`sxcf_fal2`、`gw_lmfh`)の極項。

**限界**: `wcsmear` は極を「跨いで拾う」状態(O 2p、一発目の針)の荒れを均すが、
$E_F$ 直下の準位と +3 eV の準位のように**2 準位の間にプラズモン極が挟まる対**の非対角
$\tfrac12[\Sigma(\varepsilon_i)+\Sigma(\varepsilon_j)]_{ij}$ は消せない($\Sigma(\omega)$ が急変する区間を
跨ぐ対では静的化の前提が破れる)。LiTi₂O₄ の $9^3$、1000 K では反復 2 以降に $E_F$ 近傍に
この型の針が残る。RPA の $W$ に減衰の無い鋭い極がある以上、静的 QSGW の適用限界であり、
1500 K 以上(極が Landau 減衰)では出ない。

`wcsmear = false` で 1 点評価に戻せる(絶縁体で数 meV、金属で数十 meV の違い)。極項の重みを $\chi_0$ と同じテトラヘドロン
積分で作る方法は実装していない。

---

## 3.6 $\Sigma_c$ の contour 分解と準位 smearing の整合

`t_sigmaw` の核は「占有数」だけの話ではない。$\Sigma_c = G(W-v)$ を虚軸積分と実軸極項に分ける
contour 分解では、**中間準位を均す核が両方の項に同じ形で入っていないと、$\omega$ の近くにある準位で
段差が打ち消さず、$|M|^2 W_c(0)$ 級(eV 級)の誤差が出る**。核を変える・別の核を試すときは必ずここに戻ること。

### 鋭い準位の式

中間準位 $\varepsilon'$ 一つ、$w_e \equiv (\omega-\varepsilon')/2$(Hartree)、$W_c$ の行列要素 $|M|^2$ を
省いて書く(PRB 76, 165106 Eq. 57–58 の骨格):

$$
\Sigma_c^{(\varepsilon')}(\omega)
= \underbrace{-\frac{1}{\pi}\int_0^\infty d\omega'\,\frac{w_e\,W_c(i\omega')}{w_e^2+\omega'^2}}_{I(w_e)\ \text{虚軸積分}}
\;+\;\underbrace{\theta\bigl[(\varepsilon'-\omega)(\varepsilon'-E_F)<0\bigr]\; s\,W_c(|w_e|)}_{P(w_e)\ \text{実軸極項}}
\tag{12}
$$

($s=\pm1$ は $\omega \gtrless E_F$、$P$ の $W_c$ は実軸の値、$\theta$ は「準位が $\omega$ と $E_F$ の間にある」)。
$w_e \to 0$ で $I \to -\tfrac12\,\mathrm{sign}(w_e)\,W_c(0)$ という**段差**があり、極項の窓の縁
($\varepsilon'=\omega$)で $P$ が $0 \to W_c(0)$ と跳ぶ段差とちょうど打ち消して、和は $\omega$ の連続関数になる。
コードは $W_c(i\omega')$ を $W_c(0)\,e^{-(u_a\omega')^2}$ とその残りに分け、段差を含む前者を解析的に

$$
I(w_e) = -\frac1\pi\int_0^\infty d\omega'\,\frac{w_e\,[W_c(i\omega')-W_c(0)e^{-(u_a\omega')^2}]}{w_e^2+\omega'^2}
\;-\;\frac12\,\mathrm{sign}(w_e)\,e^{a^2}\mathrm{erfc}(a)\,W_c(0),\qquad a=u_a|w_e|
\tag{13}
$$

と書く(`wintz_npm`、m_sxcf_sc の core 準位用の分岐)。残りは滑らかなので `niw` 点の Gauss–Legendre で足りる。

### Gaussian の核の場合(2026-09-19 までのコード、`esmr`)

価電子準位には (13) ではなく、被積分関数から 2 次元 Gaussian の一片を落とした

$$
\frac{w_e}{w_e^2+\omega'^2}\;\to\;\frac{w_e\bigl(1-e^{-(w_e^2+\omega'^2)/2\sigma^2}\bigr)}{w_e^2+\omega'^2},
\qquad \sigma = \texttt{esmr}/2\ [\mathrm{Ha}]
\tag{14}
$$

が使われていた。落とした一片の解析積分は $+\tfrac12\mathrm{sign}(w_e)\,e^{a^2}\mathrm{erfc}\bigl(\sqrt{a^2+w_e^2/2\sigma^2}\bigr)W_c(0)$ で、
$u_a \to 0$ の極限で (13) の段差項と合わせると

$$
-\tfrac12\,\mathrm{sign}(w_e) + \tfrac12\,\mathrm{sign}(w_e)\,\mathrm{erfc}\!\Bigl(\frac{|w_e|}{\sqrt2\sigma}\Bigr)
= -\tfrac12\,\mathrm{erf}\!\Bigl(\frac{w_e}{\sqrt2\sigma}\Bigr)
\tag{15}
$$

— **虚軸側の段差を Gaussian(標準偏差 $\sigma$)で均した**ものである。極項側の窓重み (10) も同じ $\sigma$ の
Gaussian 累積 $\Phi_G$ だったから、両側が整合して「Gaussian で均した準位の $\Sigma_c$」になっていた。
つまり Eq. 57 の `sig` は数値正則化ではなく、**準位 smearing そのもの**だった。

### 極項だけ核を変えるとどうなるか

極項の重みを Fermi–Dirac の $\Phi_{FD}$ にし、虚軸側を (14) のままにすると、準位ごとに
$[\Phi_{FD}-\Phi_G](\varepsilon'-\omega)\times|M|^2 W_c(0)$ が打ち消されずに残る。$|\varepsilon'-\omega|$ が
数 $k_BT$ の準位でだけ効くが、$q\to0$ の頭 $W_c(0)$ が乗るので大きさは eV 級になりうる。
実測(NiO 2³、L 点の O 2s 対、LDA で 70 meV 差、`t_sigmaw = 262 K`):

**表 7**. 核を片側だけ変えたときの $\Sigma_c$(NiO 2³、L 点の O 2s の対)

| SEc [eV] (state 1 / 2) | 両側 Gaussian | 極項だけ FD | 両側 FD(現在) wcsmear=false | 両側 FD、wcsmear=true 30 / 262 / 473 K |
|---|---|---|---|---|
| | 10.523 / 10.541 | 9.909 / 10.992 | **10.523 / 10.541**(QSGW 固有値も一致) | 10.495 / 10.503 / 10.510 ; 10.524 / 10.528 / 10.533 |

極項だけ FD では ±0.6 eV(反相関)、QSGW 固有値で 0.47 eV、しかも幅に強く依存した
(30/100/150/262/473 K で 10.54/10.48/10.31/9.91/9.72 eV)。両側を揃えると幅依存はほぼ消える。

### 現在の式 — 両側とも FD 核で平均する

核 $g(x) = -\partial f/\partial x$(幅 $k_BT$ = `t_sigmaw`)で均した準位の $\Sigma_c$ を、

$$
\bar\Sigma_c(\omega) = \int dx\,g(x)\,\Sigma_c^{(\varepsilon'+x)}(\omega)
= \int dx\,g(x)\,I\!\Bigl(\tfrac{\omega-\varepsilon'-x}{2}\Bigr)
+ \int dx\,g(x)\,P\!\Bigl(\tfrac{\omega-\varepsilon'-x}{2}\Bigr)
\tag{16}
$$

と**両方の項で**平均する。第 2 項が `wcsmear`(窓重み $\Phi_{FD}$、$W_c$ を核で積分、§3.5)。
第 1 項は鋭い式 (13) を使い、不連続な $-\tfrac12\mathrm{sign}(w_e)$ だけを解析的に

$$
\int dx\,g(x)\Bigl[-\tfrac12\,\mathrm{sign}(\omega-\varepsilon'-x)\Bigr]
= -\tfrac12\bigl[2\,\Phi_{FD}(\omega-\varepsilon')-1\bigr],\qquad
\Phi_{FD}(y) = \frac{1}{1+e^{-y/k_BT}}
\tag{17}
$$

で置き、連続な残り($-\tfrac12\mathrm{sign}(w_e)\,(e^{a^2}\mathrm{erfc}(a)-1)$ と数値部)は
$x_j = k_BT\ln\frac{u_j}{1-u_j}$、$u_j = (j-\tfrac12)/40$ の中点則で平均する
(FD の累積分布が測度なので $u$ が一様)。$|\omega-\varepsilon'| > 30\,k_BT$ の準位は被積分関数が滑らかなので
鋭い式を 1 点で使う。core 準位($k_BT = 0$)は (13) のまま。`sig` は消えた。
実装: `m_sxcf_sc.f90` の `CorrelationSelfEnergyImagAxis` ブロック(QSGW)、`wintzsg.f90` の `wintzsg_npm`(一発 GW)。

注意:
- `wcsmear = false` では第 2 項が「重みだけ核で平均、$W_c$ は平均エネルギー 1 点」になるので (16) の厳密な
  平均ではない(差は $O(k_BT^2\,W_c'')$)。段差の打ち消しは重みだけで決まるので、こちらは崩れない。
- $W_c(i\omega') \simeq W_c(0)e^{-(u_a\omega')^2}$ の分離と `niw` 点の精度は従来どおり。
- 求積のコストは (itp, it) 対あたり $40\times$`niw` の scalar 演算(CPU 側、GPU 転送前)。9³ LiTi₂O₄ でも rank あたり数秒。
- **核を変えるときの規則**: 占有重み(`wfacx`)、極項の窓重み(`wfacx2` / `pole_weights`)、虚軸の段差 (17) の
  3 か所が同じ $\Phi$ を指していること。Gaussian に戻すなら (14)–(15) のように 3 か所とも戻す。

経緯と数値は `MD/research_log.md` 2026-09-20 12:17 / 13:10。

---

## 4. 幅と温度の選び方

- **絶縁体・半導体**: ギャップが $k_BT$ より十分大きければ、幅は結果に効かない。テンプレートの 300 K のままでよい
  (Si: `t_tetrakbt` = 0 と 300 K のギャップの差は 0.001 meV 以下、1000 K で −1 meV。§5.1 の表 9)。
- **金属の `t_sigmaw`**: 原理的には $n_1n_2n_3\to\infty$、`t_sigmaw` → 0 の極限が求めるものだが、有限の k メッシュでは、
  メッシュが Fermi 面の近くで持つ準位の間隔と同じ程度以上の幅が要る。目安は

$$
k_BT \gtrsim \frac{v_F \cdot (\text{BZ の大きさ})}{N}
\tag{18}
$$

  ($N$ は一方向のメッシュの数)。LiTi₂O₄ の $6^3$・$9^3$ では 1000 K($k_BT = 0.086$ eV)。
- **`t_tetrakbt`**: 電子温度の物理を見るなら正の値にして `t_sigmaw` も同じ温度にする。$W$ の極を数値的に均すだけなら負の値
  (Gaussian)にする。このとき Fermi 準位と占有は T=0 のままである。
- 求めたい量が幅について収束していること、できれば同じ幅で 2 つの k メッシュが一致することを確かめる。
- 電子温度が入るのは $\chi_0$ と $\Sigma$ だけである。`lmf` が作る密度は T=0 の占有(テトラヘドロン法)のままで、格子の温度(フォノン)も入らない。

---

## 5. サンプル — `Samples/kBT/`

ecalj の `Samples/kBT/` にある（表 8）。数値と図の全部、回し方は各ディレクトリの README。

**表 8**. サンプル

| ディレクトリ | 中身 | 時間 |
|---|---|---|
| `scanT/` | 温度のスキャン。`scan_T.sh` が同じ入力を `t_tetrakbt` と `t_sigmaw` だけ変えて QSGW で回し、`collect_T.py` が表と図にする。Si・GaAs(§5.1)と bcc Fe・fcc Cu(§5.2) | 1 つの設定が 2〜10 分 |
| `LiTi2O4/` | 金属スピネル LiTi₂O₄ の MLO-QSGW(§5.3) | $9^3$ の 1 反復が 21 分(RTX 5090 × 2) |
| `scanT/Si`、`Samples/TestInstall/fe_kbt` | `testecalj` の試験(Si と Fe の 3000 K、1 反復) | 20〜40 秒 |

どのサンプルでも、電子温度が入るのは $\chi_0$ と $\Sigma$ だけである(§4)。`scanT` は k メッシュが粗く反復も少ないので、
数値は温度による変化を見るためのもので、収束した値ではない。

### 5.1 Si と GaAs — ギャップと電子温度

`nkabc` = `n1n2n3` = $4^3$、LDA から QSGW 5 反復。`t_tetrakbt = t_sigmaw = T`。

**表 9**. Si のギャップ(SCF の k メッシュの上での最小。QSGW 5 反復後)

| `t_tetrakbt` (K) | `t_sigmaw` (K) | ギャップ (eV) | $T=0$ からの変化 (meV) |
|---|---|---|---|
| 0 | 300 | 1.223 | 0 |
| 300 | 300 | 1.223 | 0 |
| 1000 | 1000 | 1.222 | −1 |
| 2000 | 2000 | 1.188 | −35 |
| 3000 | 3000 | 1.066 | −157 |
| 5000 | 5000 | 0.913 | −310 |
| 3000 | 300 | 1.090 | −133 |
| 0 | 3000 | 1.176 | −46 |
| −3000 | 3000 | 1.173 | −50 |

![Si のギャップと電子温度](kBT/si_T_gap.png)

**図 1**. 左: Si のギャップと電子温度(`t_tetrakbt = t_sigmaw = T`)。右: 反復ごとのギャップ

- 300 K では $T=0$ と変わらない(差は 0.001 meV 以下)。2000 K から縮みが目立ち、3000 K で 0.16 eV 縮む。
- 価電子帯の形はほとんど変わらず、伝導帯が価電子帯に対して全体に下がる。
- 縮みの主な原因は $\chi_0$ 側である(3000 K の −157 meV のうち、$\chi_0$ だけを 3000 K にすると −133 meV、$\Sigma$ だけだと −46 meV)。
  ギャップを越えて熱励起されたキャリアが遮蔽に加わり、$W$ が小さくなる。
- `t_tetrakbt = -3000`(Gaussian)では、$\Sigma$ だけの場合とほぼ同じ(−50 meV)。Gaussian は占有を変えないので、キャリアの遮蔽は入らない(§2.4)。
- GaAs($4^3$、QSGW 5 反復、直接ギャップ 1.415 eV)も同じ傾向である: 300 K で変化なし、1000 K で −4 meV、2000 K で −20 meV、3000 K で −95 meV、
  5000 K で −267 meV。3000 K では $\chi_0$ だけで −83 meV、$\Sigma$ だけで −39 meV(数値と図は `Samples/kBT/scanT/README.md` の表 3)。

### 5.2 bcc Fe と fcc Cu — 金属のバンドと磁気モーメント

`nkabc` = $7^3$、`n1n2n3` = $5^3$、スピン分極、LDA から QSGW 3 反復。

**表 10**. bcc Fe(QSGW 3 反復後)。交換分裂は Γ の d バンド(バンド 5・6)のスピン 2 とスピン 1 の差

| `t_tetrakbt` (K) | `t_sigmaw` (K) | 磁気モーメント ($\mu_B$) | `EFERMI_kbt` − `EFERMI` (eV) | 交換分裂 (eV) |
|---|---|---|---|---|
| 0 | 300 | 2.429 | — | 2.619 |
| 300 | 300 | 2.426 | +0.006 | 2.604 |
| 1000 | 1000 | 2.443 | +0.045 | 2.688 |
| 3000 | 3000 | 2.420 | +0.272 | 2.492 |
| 5000 | 5000 | 2.252 | +0.279 | 1.842 |
| 3000 | 300 | 2.441 | +0.258 | 2.676 |
| 0 | 3000 | 2.409 | — | 2.450 |
| −1000 | 1000 | 2.443 | — | 2.685 |

- 300 K では $T=0$ とほとんど変わらない(対称点でのバンドの差は最大 14 meV)。
- 1000 K までのバンドの動きは最大 63 meV で、温度について単調でない。$5^3$ のメッシュでは、$E_F$ の近くのどの状態が分数占有になるかで
  この程度は動く(§9 の 1)。
- 3000 K から交換分裂とモーメントが減る。金属では $\Sigma$ 側が効く: 3000 K の交換分裂の変化 −0.13 eV は、$\Sigma$ だけを 3000 K にした場合(−0.17 eV)に近く、
  $\chi_0$ だけでは +0.06 eV。Si と逆である。
- `t_tetrakbt = -1000`(Gaussian)と `t_tetrakbt = 1000`(有限温度)は、この温度ではほとんど同じ結果になる(バンドの差は最大 7 meV)。
- 非磁性の fcc Cu(`nkabc` = $8^3$、`n1n2n3` = $6^3$、QSGW 3 反復)では、1000 K までのバンドの動きは 11 meV 以下。3000 K から d バンドが $E_F$ に近づく
  (Γ の d バンドが 3000 K で +61 meV、5000 K で +130 meV 前後)。$\chi_0$ だけで +46 meV、$\Sigma$ だけで +24 meV で、両方が効く
  (`Samples/kBT/scanT/README.md` の表 4)。

### 5.3 LiTi₂O₄ — MLO-QSGW

金属スピネル LiTi₂O₄(14 原子)を、冒頭の「LiTi₂O₄ の設定」で、$\Sigma$ を MLO で内挿する QSGW([mlo_gwsc](./mlo_gwsc))で
LDA から 40 反復まで回した。まとめは ecalj の `Samples/kBT/LiTi2O4/README.md`。

- GPU の行列積の精度(tf32、fp32、fp64)による MLO バンドの差は数 meV 以下。
- $6^3$ は 10 反復まではバンドがなめらかで、20 反復から占有側のバンドに肩ができ、40 反復でも残る。
- $9^3$ はメッシュ点の間の波打ちがなかなか減らないが、40 反復でかなりなめらかになる。メッシュ点を結ぶ内挿の曲線は、わずかに波打って見える。

---

## 6. 実装の場所

**表 3**. 有限温度と均し方の実装(`SRC/subroutines/`)

| ファイル | 役割 |
|---|---|
| `m_GWinput.f90` | `[gw]` の `t_tetrakbt`(必須、表 1)、`t_sigmaw`、`wcsmear` を読む。`t_tetrakbt < 0` を $\sigma$(内部の `SmearX0`、Ha)に直す(式 (7)) |
| `tetwt5.f90` | $\chi_0$ のテトラヘドロンの重み。`lindtet6`(式 (3))と `lindtet6_kbt`(式 (5)) |
| `fpiint.f90` | `gausq_fd`: 熱核の 20 点(式 (5)(6)) |
| `main_heftet.f90` | `fermi_kbt_tetra`: `EFERMI_kbt`(式 (6)) |
| `m_tetwt.f90` | `t_tetrakbt > 0` なら `lindtet6_kbt` と `EFERMI_kbt` を使う |
| `dpsion5.f90`、`hilbertmat.f90` | Hilbert 変換(式 (2))と Gaussian(`gaussianfilterhis`、`smearx0_apply`、式 (8)) |
| `genallcf_mod.f90` | `sigmakbt_setup`: $\Sigma$ の Fermi 準位(`t_tetrakbt > 0` なら `EFERMI_kbt`、無ければ止まる) |
| `wfacx.f90` | $\Sigma$ 側の Fermi–Dirac 核(`fd_cdf`、`fd_iav`。幅 `t_sigmaw`、式 (10)) |
| `m_sxcf_sc.f90`、`sxcf_fal2.f90` | $\Sigma_c$ の実軸極項と `wcsmear`(§3.5、§3.6) |

---

## 7. 実装の検証と限界

実装を読み直して確かめたことと、残っている限界。$\Sigma$ 側の contour の恒等式(虚軸 + 極項)は
[§3.6](#_3-6-sigma-c-の-contour-分解と準位-smearing-の整合)。

### 7.1 恒等式は厳密 (確認済み)

`lindtet6_kbt`(`tetwt5.f90`)の method B′ は式 (4)(§2.3)を **k 積分の内側**で使っている。$E$ 積分と $k$ 積分が
交換するので厳密である。得られる分子は Adler–Wiser の $f_a-f_b$ であって $f_a(1-f_b)$ ではない
— 応答関数として正しいのはこちら。

個別に確認したもの:

- **ゲート**(`lindtet6_kbt` の冒頭) — a/b が窓の外で平坦になる場合分けを全部たどったが、
  どれも $T=0$ の値と一致する
- **node zero-skip**(`lindtet6_kbt` の $E$ のループ) — `lindtet6` が恒等的に 0 になる条件そのもの。
  `knorm` は全節点で正規化されるので、飛ばしても整合する
- `efermia` ≠ `efermib` でも (4) は成立する(それぞれの Fermi 準位での $f$ になる)
- **`EFERMI_kbt`**(`fermi_kbt_tetra`、`main_heftet.f90`、式 (6))
  — NOS を熱核で畳み込んで二分法。$T\to0$ でテトラヘドロンの $E_F$ に厳密に戻る。
  χ0 側もこれを読んでいる(`m_tetwt`)

### 7.2 限界 — $E$ 積分の求積

式 (5)(6) の $E$ 積分は 20 点で、$t\in[-6,6]$ を $[-6,-1],[-1,0],[0,1],[1,6]$ の 4 区間に分けた 5 点 Gauss-Legendre である
(`gausq_fd`)。$E_F$ に最も近い節点は $|t|=0.047$($|E-E_F|=0.094\,k_BT$)。テトラヘドロンの中のバンドの幅を $\Delta$ として、
占有の因子の最悪の誤差を測ったものが表 11 である。

**表 11**. $E$ 積分の求積と占有の因子の最悪の誤差

| 求積 | 節点 | $\Delta=0.4k_BT$ | $1k_BT$ | $2k_BT$ | $4k_BT$ | $8k_BT$ |
|---|---|---|---|---|---|---|
| **いまの実装**: 4 区間 $[-6,-1,0,1,6]$ × GL5 | 20 | 2.9e-2 | 1.6e-2 | 4.8e-3 | 2.6e-3 | 1.3e-3 |
| 6 区間 × GL4 | 24 | 2.9e-2 | 1.1e-2 | 2.8e-3 | 1.7e-3 | 9.6e-4 |
| 2 区間 $[-6,0,6]$ × GL10 | 20 | 5.5e-2 | 1.9e-2 | 1.1e-2 | 3.8e-3 | 1.9e-3 |
| 1 区間 $[-6,6]$ × GL20(2026-09-26 までのコード) | 20 | 1.7e-1 | 9.7e-2 | 2.3e-2 | 1.3e-2 | 6.8e-3 |

分散の大きいバンドでは誤差は小さい。平坦なバンド、バンドの端、van Hove 特異点が $E_F$ にあると、占有の因子が 0.01〜0.03 ずれる。
節点を増やしても誤差は $1/N$ でしか減らず、区間の切り方のほうが効く。被積分関数の折れ点(四面体の 8 個の角のエネルギー)で
区間を切るのが本筋だが、実装していない。Fe 3000 K $5^3$ の `EFERMI_kbt` の誤差(600 点の求積に対して)は 1.2 meV。

### 7.3 $\Sigma$ 側の状態の窓

中間状態の候補の窓と実軸の $\omega$ メッシュの余裕は、核の裾 $E_F\pm15\,k_BT$($k_BT$ は `t_sigmaw`、`sig_window`)で決まる。
切り捨てられる占有は $f(15)=3\times10^{-7}$ で、`t_sigmaw` を上げれば窓も広がる。窓の中心と重みの中心は同じ Fermi 準位
(`t_tetrakbt > 0` なら `EFERMI_kbt`)である。2026-09-19 までのコードでは窓が `esmr` で決まっていて、高温で占有が切り捨てられた
([記録](./kBT_history) §4)。

### 7.4 範囲 — $\Sigma$ が要るのは「差」ではなく「和」

$\chi_0$ を有限温度化しただけでは $\Sigma$ は揃わない。理由は一行で書ける。

中間状態の粒子‑正孔対を $(2,3)$、$\omega'=\varepsilon_2-\varepsilon_3>0$ として、

**表 12**. 粒子‑正孔対の重み

| | 重み |
|---|---|
| 対を**作る** | $f_3(1-f_2)$ |
| すでにある対を**壊す** | $f_2(1-f_3)$ |

$\mathrm{Im}\,\chi_0$ が持っているのはこの**差** $f_3(1-f_2)-f_2(1-f_3)=f_3-f_2$ だが、
$\Sigma$ の中間状態の重みは**和**である(どちらも同じ分母を持つ)。
いまの実装は $\mathrm{Im}W$ から差を取ってきて外側の $f_j$ を掛けるので、
**$f_2(1-f_3)$ の分が落ちる。** 詳細釣り合い $f_2(1-f_3)/f_3(1-f_2)=e^{-\beta\omega'}$ から

$$
f_2(1-f_3)=\bigl(f_3-f_2\bigr)\,n_B(\omega')
\tag{19}
$$

なので、これは普通 Bose 因子と呼ばれるものである。**ただし出どころは Fermi の
占有数**であって、ボソンの熱浴を外から入れたわけではない。導出は
[§8](#_8-付録-—-虚時間を使わない定式化) を見よ。物理は単純で、$T>0$ では
粒子‑正孔対がすでに熱励起されていて電子はそれを吸収できる、というだけのこと。

**大きさ。** $f_2(1-f_3)$ は $\omega'\gg k_BT$ で指数的に小さいので、効くのは
$\omega'\lesssim k_BT$ の粒子‑正孔連続体だけである。金属で
$\mathrm{Im}W_c\simeq c\,\omega'$ とすると

$$
\int_0^\infty\!\! d\omega'\,\mathrm{Im}W_c(\omega')\,n_B(\omega')
= c\,(k_BT)^2\frac{\pi^2}{6}
\tag{20}
$$

2000 K で $(k_BT)^2=0.029$ eV²、$\omega-\varepsilon_j\sim1$ eV なら数 meV–数十 meV。

**これは欠陥ではなく定義である。** ecalj のスキームは「温度 $T$ の占有数で作る
有効一体ハミルトニアン」であって、Mermin の有限温度 DFT が Fermi 占有数だけで
閉じているのと同じ立場である。QSGW は動的な $\Sigma$ ではなく**静的エルミートな**
一体ハミルトニアンを作るものなので、落ちているのが主に詳細釣り合い(= 寿命)側で
あることもあって、影響は一発 GW のスペクトル関数を出す場合よりずっと軽い。
平衡の多体摂動論の $\Sigma$ が欲しいなら (19) が要る、というだけのことである。

### 7.5 core の交換 (CoreEx)

core との交換(`hsfp0_sc --job=3`)は、有限温度でも Fermi 準位を価電子の底の下に置いて、和に入る状態を core だけにする。
core の準位は $k_BT$ よりずっと深いので、温度によらず T=0 の値になる。

**2026-06-26 から 09-27 までのコードの有限温度の結果についての注意**: この間のコードには CoreEx に誤りが 2 つあり、
core の状態の数 $n_\mathrm{ctot}$ が占有バンドの数より少ない系(1s だけが core の軽元素の化合物など)では、core との交換に
価電子の状態の一部が入っていた。`lsxC` の `CoreEx mode: ef nspin nctot=` と `valn=`(価電子数)を比べれば判定できる
(semicore を価電子に入れた系は core が減るので注意)。Fe($n_\mathrm{ctot}=9$)と LiTi₂O₄($n_\mathrm{ctot}=46$)には影響が無い。
誤りの中身は [記録](./kBT_history) §5。

---

## 8. 付録 — 虚時間を使わない定式化

$n_B$ がどこから来るのかは、松原形式を経由しなくてもはっきりする。むしろ
そちらのほうが**すべてが Fermi の占有数から出る**ことが見えてよい。
教科書では閉時間径路と平衡の話が分かれて書かれていることが多いので、
ここで通して書いておく。

以下 $\hbar=1$、エネルギーは $\mu$ から測る。
$f(\varepsilon)=1/(e^{\beta\varepsilon}+1)$、$n_B(\omega)=1/(e^{\beta\omega}-1)$。

### 8.1 $T=0$ の議論がなぜそのまま使えないか

$T=0$ の実時間摂動論は Gell-Mann–Low に依っている。断熱的に相互作用を入れると
$U(\infty,-\infty)|\Phi_0\rangle = e^{i\theta}|\Phi_0\rangle$、つまり
**基底状態は位相を除いて自分自身に戻る**ので、$+\infty$ 側を $-\infty$ 側と
取り替えられて

$$
\langle\Psi_0|T\{\cdots\}|\Psi_0\rangle
=\frac{\langle\Phi_0|T\{S\cdots\}|\Phi_0\rangle}{\langle\Phi_0|S|\Phi_0\rangle}
\tag{21}
$$

と片道の時間順序積で書ける。

$T>0$ ではこれが使えない。$\rho_0=e^{-\beta(H_0-\mu N)}/Z_0$ は固有状態ではなく、
断熱的に発展させても**重みが非相互作用系のまま**($e^{-\beta E_n^{(0)}}$ であって
$e^{-\beta E_n}$ ではない)なので、$+\infty$ 側を $-\infty$ 側と同一視できない。
(21) の分母に相当するものが書けない、というのが問題の本質である。

### 8.2 閉時間径路 — 行って戻る

そこで **$+\infty$ の状態を一切使わない**。$-\infty\to+\infty\to-\infty$ と
往復する径路 $C$ を取れば、演算子を挟まない限り

$$
\mathrm{Tr}\bigl[\rho_0\,U(-\infty,+\infty)\,U(+\infty,-\infty)\bigr]
=\mathrm{Tr}\,\rho_0 = 1
\tag{22}
$$

が**恒等的に**成り立つ。分母が要らない。これが閉時間径路(Keldysh)の全部である。

径路順序積 $T_C$(往路の後に復路が来る順序)を使って

$$
G(1,2)=-i\,\mathrm{Tr}\Bigl[\rho_0\,T_C\Bigl\{
e^{-i\int_C dt\,H_1(t)}\;\psi(1)\psi^\dagger(2)\Bigr\}\Bigr]
\tag{23}
$$

場は $H_0$ の相互作用表示。$\tau$ ではなく実時間で、$T$ 積が $T_C$ 積になっただけ
である。$1,2$ をどちらの枝に置くかで 4 つの成分が出るが、独立なのは 2 つで、
以下では

$$
G^<(1,2)=+i\langle\psi^\dagger(2)\psi(1)\rangle,\qquad
G^>(1,2)=-i\langle\psi(1)\psi^\dagger(2)\rangle
\tag{24}
$$

を使う。

### 8.3 Wick の定理 — 入力は $f$ だけ

$\rho_0$ は $H_0$ について Gauss 的なので Wick の定理がそのまま成立し、
縮約は自由な径路伝播関数になる。エネルギー $\varepsilon_j$ の準位について

$$
g_j^<(\omega)=2\pi i\,f_j\,\delta(\omega-\varepsilon_j),\qquad
g_j^>(\omega)=-2\pi i\,(1-f_j)\,\delta(\omega-\varepsilon_j)
\tag{25}
$$

**温度が入るのはここだけで、入るのは Fermi 分布だけである。**
Bose 分布はこの段階でどこにも無い。

### 8.4 GW を径路上で書く

径路引数のまま、形は $T=0$ と同じ:

$$
\chi_0(1,2)=-i\,G(1,2)G(2,1),\qquad
W=v+v\chi_0 W,\qquad
\Sigma(1,2)=i\,G(1,2)W(2,1)
\tag{26}
$$

実時間成分に落とすには Langreth 則を使う。同じ引数の積
$C(1,2)=A(1,2)B(2,1)$ に対して

$$
C^{\gtrless}(1,2)=A^{\gtrless}(1,2)\,B^{\lessgtr}(2,1)
\tag{27}
$$

**$\gtrless$ がひっくり返る**のが要点である。

### 8.5 $\chi_0$ の $\gtrless$ 成分 — ここで $n_B$ が出る

(26)(27) と (25) から、時間並進対称性を使って

$$
\chi_0^>(\omega)=-2\pi i\sum_{23} f_3(1-f_2)\,\delta(\omega-\varepsilon_2+\varepsilon_3)
$$
$$
\chi_0^<(\omega)=-2\pi i\sum_{23} f_2(1-f_3)\,\delta(\omega-\varepsilon_2+\varepsilon_3)
\tag{28}
$$

$\chi_0^>$ が**対を作る**過程、$\chi_0^<$ が**すでにある対を壊す**過程である。
$\omega>0$ では $\varepsilon_2>\varepsilon_3$。

スペクトル関数は差のほうで、

$$
\chi_0^>-\chi_0^< \;\propto\; f_3(1-f_2)-f_2(1-f_3)=f_3-f_2
\;=\;2i\,\mathrm{Im}\chi_0^R
\tag{29}
$$

これが $T=0$ のコードが計算している量である。一方 (28) の比は
$f_j=1/(e^{\beta\varepsilon_j}+1)$ を代入するだけで

$$
\frac{\chi_0^<(\omega)}{\chi_0^>(\omega)}
=\frac{f_2(1-f_3)}{f_3(1-f_2)}=e^{-\beta\omega}
\tag{30}
$$

となり、$1/(1-e^{-\beta\omega})=1+n_B(\omega)$ を使えば

$$
\chi_0^>=\bigl(1+n_B\bigr)\bigl(\chi_0^>-\chi_0^<\bigr),\qquad
\chi_0^<=n_B\,\bigl(\chi_0^>-\chi_0^<\bigr)
\tag{31}
$$

**$n_B$ はここで初めて現れる。仮定ではなく (25) の Fermi 分布からの帰結**である。
ボソンの熱浴を外から入れた覚えはないのに Bose 分布が出るのは、
$\chi_0$ がボソン的な相関関数だから((30) がボソンの KMS 条件そのもの)。

RPA の衣を着せても比は変わらない。$W^{\gtrless}=\epsilon^{-1,R}\,v\chi_0^{\gtrless}v\,
\epsilon^{-1,A}$ で $\epsilon^{-1,A}=(\epsilon^{-1,R})^\dagger$ だから、
(30) はそのまま $W^{\gtrless}$ に受け継がれる。

### 8.6 $\Sigma$ の $\gtrless$ 成分

(26)(27) より $\Sigma^{\gtrless}(t)=i\,G^{\gtrless}(t)\,W^{\lessgtr}(-t)$、
振動数では

$$
\Sigma^{\gtrless}(\omega)=i\sum_j\int\!\frac{d\nu}{2\pi}\,
g_j^{\gtrless}(\nu)\,W^{\lessgtr}(\nu-\omega)
\tag{32}
$$

(25) を入れ、$B(\Omega)\equiv W^>(\Omega)-W^<(\Omega)=2i\,\mathrm{Im}W^R(\Omega)$、
$W^>=(1+n_B)B$、$W^<=n_B B$ とすると、
$2i\,\mathrm{Im}\Sigma^R=\Sigma^>-\Sigma^<$ から

$$
\mathrm{Im}\,\Sigma^R(\omega)=\sum_j
\mathrm{Im}W^R(\varepsilon_j-\omega)\;
\bigl[\,n_B(\varepsilon_j-\omega)+f_j\,\bigr]
\tag{33}
$$

$\varepsilon_j-\omega$ の符号で分けて $n_B(-\Omega)=-(1+n_B(\Omega))$ を使えば、
見慣れた放出因子 $1-f_j+n_B$ と吸収因子 $f_j+n_B$ になる。実部は
Kramers–Kronig で決まる。

$T\to0$ の確認: $n_B(\Omega>0)\to0$、$n_B(\Omega<0)\to-1$、$f_j\to\theta(-\varepsilon_j)$。
$\omega>0$ の準粒子について、$j$ が占有だと $[-1+1]=0$(Pauli 阻止)、
$j$ が空だと $[-1+0]=-1$ で寄与する。従来の $T=0$ GW に戻る。

### 8.7 いまの実装が保っているもの・落としているもの

ecalj が計算しているのは $\mathrm{Im}W^R$、すなわち (29) の**差**である。
$\Sigma$ を組むときに外側の $f_j$ だけを有限温度にし、(33) の $n_B$ は入れていない。
つまり

$$
\Sigma^{\text{ecalj}}:\quad n_B\to 0,\qquad f_j\to f_j(T)
$$

これは $\Sigma$ の詳細釣り合い

$$
\Sigma^<(\omega)=-e^{-\beta\omega}\,\Sigma^>(\omega)
\tag{34}
$$

を破る((32) と (30) から (34) は**両方**を残したときにのみ成り立つ)。
したがって得られる $\Sigma$ は、厳密にはどんな温度 $T$ の平衡状態のものでもない。

一方でこれは QSGW では軽い。QSGW が作るのは**静的エルミートな**一体
ハミルトニアンで、$\mathrm{Im}\Sigma$ は最初から捨てているからである。
(34) が壊れているというのは主に寿命側の話で、QSGW が使う
$\mathrm{Re}\,\Sigma(\varepsilon_j)$ への影響は (20) の $O((k_BT)^2)$ にとどまる。

### 8.8 初期相関についての注意

8.2 では $\rho_0$(非相互作用の熱平衡)から断熱的に相互作用を入れると書いたが、
これは**初期相関を落としている**。厳密には径路に虚時間の縦枝
$[t_0,\,t_0-i\beta]$ を足した Kadanoff–Baym 径路を使い、そこに初期相関を
持たせる(Danielewicz)。平衡かつ断熱的な場合には縦枝が実時間部分から
分離し、結果は松原形式と一致する。

ここで $\tau$ が顔を出すのは**初期条件を指定するため**であって、
摂動展開そのものは実時間のままである。$n_B$ が出るかどうかとは関係がない —
それは 8.5 で見たとおり (25) の Fermi 分布だけから出る。

---

## 9. 残っている課題

1. **金属での `t_sigmaw` の選び方の検証**。$E_F$ から $\pm k_BT$ 以内の少数の状態が分数占有になり、その状態の SEx と SEc が eV の大きさで動いて
   打ち消し合う。粗いメッシュでは「どの状態が $E_F$ を跨ぐか」に敏感なので、メッシュについての収束($5^3\to9^3$)を数値で押さえる必要がある
   (Fe の温度依存は §5.2)。
2. **絶縁体の `EFERMI_kbt`**。ギャップの中で解が区間になるときは区間の中央を返す。物理的には状態密度の非対称で決まる真性の $E_F$ を返すべきで、
   `fermi_kbt_tetra` の窓($E_F\pm60\,k_BT$)と熱核の届く幅($\pm12\,k_BT$)の関係で結果が変わる(GaAs 300 K で 0.15345 Ry、T=0 の中央は 0.15221 Ry)。
3. **$O((k_BT)^2)$ の Bose 項**(§7.4、§8.7)。$f_2(1-f_3)=(f_3-f_2)\,n_B$ を落としている。3000 K で効くかは見積もっていない。
4. **$E$ 積分を被積分関数の折れ点で切る**(§7.2)。
5. **2 準位の間にプラズモン極が挟まる対の非対角要素**(§3.5 の限界)。`wcsmear` では消えない。
6. **`lmf` の密度は T=0 の占有のまま**(§4)。高温で熱励起されたキャリアの密度への寄与は入っていない。
7. **LiTi₂O₄ の $9^3$、3000 K** は反復を追える計算が無い。
8. **`t_tetrakbt` ≤ 0 のときの $\Sigma$ の Fermi 準位**は Gaussian の数え上げで決め(§0 の表 5)、準位を均す核は Fermi-Dirac である。
   同じ核に揃えていない。

---

## 10. gwkbt-dev について

この頁の内容はすべて `main` に入っている。ブランチを使い分ける必要はない。
別系統の `gwkbt-dev` は、有限温度の $\Sigma$ を別の設計(gwkbt Stage A/B)で作りかけたもので、使える状態ではない。
この頁の `t_sigmaw` とは別物である。

---

## 関連

- [QSGW の計算 (gwsc)](./gwsc)
- サンプル `ecalj/Samples/kBT/`(§5)
- [kBT の記録](./kBT_history): 2026-06〜09 の経緯、当時の計算の説明、直した誤り
- [`ctrlg.<sname>.toml` の `[gw]` セクション](./lmf#file-structure-sections)
- 一発 GW の QP エネルギー、任意 k 線上の評価、$W$ と実軸積分の診断は
  ecalj の `MD/past_log.md` §5 を見よ(2026-06 の手引き `FiniteT_and_QPE_HOWTO.md` の要点。本頁はその §1 を発展させたもの)
