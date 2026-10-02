# GWinput (legacy)

> ⚠️ The Fortran programs do not read the `GWinput` text file. The settings of this page are keys of `ctrlg.<sname>.toml` (sections `[gw]` / `[mlo]` / `[blocks]` / `[product_basis]`; the last one also holds the per-atom tables), written in TOML syntax (`key = value`). For a directory with the legacy `ctrl.<sname>` + `GWinput`, run `Legacy2toml.py <sname>` to convert. See [TOML migration](./toml_migration) for the key mapping.

The `[gw]` section of `ctrlg.<sname>.toml` sets the computatioanl conditions for  GW/QSGW calculations.

## mkGWinputによる自動生成
- GW の設定（`ctrlg.<sname>.toml` の `[gw]`・`[mlo]`・`[blocks]`・`[product_basis]`）は自動生成することができ、殆どのパラメータはデフォルト設定で使用できる。

`ctrls.<sname>` から始めるとき（`<sname>` は物質名）は
```bash
ctrlgenToml.py <sname>
```
が GW の設定も含めて `ctrlg.<sname>.toml` を書く。すでにある `ctrlg.<sname>.toml` に GW の設定を付けるには
```bash
mkGWinput <sname>                # GW の設定を作り直す。元のファイルは ctrlg.<sname>.toml.bakup に残る
ctrlgenToml.py <sname> --addgw   # [gw] の無いファイルに足す（[gw] があれば止まる）
```
> [!INFO]
> どちらも内部で ecalj package の実行ファイル（`lmfa`、`lmf --jobgw=0`、`gwinit`）を `mpirun` で実行する。
> 実行ファイルはスクリプトと同じディレクトリ(ecalj package のバイナリinstall ディレクトリ)のものが使われる。

生成される `[gw]` は `n1n2n3 = [4, 4, 4]`、`t_sigmaw = 300`、`t_tetrakbt = 300` である。計算に合わせて書き換える
（絶縁体・半導体で $\chi_0$ を均さないなら `t_tetrakbt = 0`）。


## Format
- TOML の文法で書く。`#` から行末まではコメント
- `key = value` の形式であり記載順序は問わないが, 同じkeyを複数記述しないようにする。論理値は `true` / `false`、
  リストは `[4, 4, 4]`、実数は `1e-5`（`1d-5` とは書けない）、複数行の表（`QforGW` など）は `"""` で囲む

## Parameters in GWinput
生成された `[gw]` に, パラメータの意味は記載されていますが, 以下に重要な変数と記載がない変数について記述します。

### `n1n2n3`
$GW$自己エネルギー計算におけるk点の数。`[bz] nkabc`（`lmf` のk点）とは異なる点を指定できる（`gwsc --mlo` では `nkabc` と同じにする。違うと止まる）。 計算コストは, このk点数の2乗に比例することに注意。
- **type** : integer list

We set number of k points for self-energy in `[gw]`. These can be 1/2 or 2/3 of k points specified in `[bz] nkabc`, in order to reduce computational time.
The 6x6x6 k points is feasible setting for ZB structure (2 atoms per cell). That is, 6x6x6x2 =432 \sim `k points \times atom number`.
For example, when we try 8 atoms per cell, 4x4x4 or 3x3x3
is fine because 4x4x4x8 \sim 400.
For metallic systems, larger is fine, but limited by computational
time. See takao kotani's papers, for example, https://doi.org/10.1103/PhysRevB.93.075125
Good news is that we don't need to use so many k points as in `[bz] nkabc`.
In cases, 4x4x4 for ZB is not so bad --- roughly speaking, this is a low limit for
publication. This means 3x3x3x8 (for 8 atom case) is not so bad. (3x3x3x8 > 4x4x4x2).
See our examinations shown in  https://doi.org/10.7566/JPSJ.83.094711 and
http://doi.org/10.7567/JJAP.55.051201

### `KeepEigen`  
波動関数をメモリに保持するかどうか。波動関数は`CPHI` および`GEIG`というファイルで出力される。それらが1MPIが使用できるメモリに対して, 同等もしくはそれ以上である場合は必ず`false`にする。 $k$点や軌道数が多い場合は`false`にする。
- **type**: boolean
- **default**: true

### `zmel_batch_gb`
行列要素(zmel)を一度に作るバッチサイズ:単位GB（1MPIランクあたり）。1MPIが使用できるメモリサイズの3分の1から4分の1程度を記載すると良い。自己エネルギーの計算と $W$（分極関数）の計算の両方で使用される。大規模系では大きい方が計算は早くなるが，メモリ不足で計算が落ちる場合は小さくする。
- **type** : float
- **default** : 0.4（CPU では 0.4 未満、GPU では 2.0 未満を書いても、それぞれ 0.4、2.0 に上げられる）

系に対して値が小さいと`sxcf_fal2_count.sc. Too small memory for nmbatch. Enlarge zmel_batch_gb in ctrlg.<sname>.toml [gw]` と出て計算が終了する場合がある。その場合は大きくする。

### `zmel_max_size`
`zmel_batch_gb` の古い名前。今も `zmel_batch_gb` として読まれる。古い入力にある `MEMnmbatch` は読まれない（`zmel_batch_gb` を書く）。

### `KeepPpb`
MT内の積基底とMT内波動関数基底の行列要素(`ppb`変数)を、対称操作のすべてについてメモリに保持するかどうか。通常はデフォルト値`false`（そのとき使う対称操作の分だけを持つ）で良い。`true` にすると使用メモリ量が対称操作の数だけ増える。
- **type**: boolean
- **default**: false

### `t_tetrakbt`  *(必須、2026-09-28)*
$\chi_0$ をどう均すかを温度 (K) の符号で選ぶ。書かないと止まる。
- `0`: 均さない(T=0 のテトラヘドロン)。絶縁体・半導体。
- `T > 0`: 有限温度のテトラヘドロン。$\chi_0$ の占有数が温度 $T$ の Fermi-Dirac になり、Fermi 準位は `EFERMI_kbt`($\Sigma$ 側も)。[kBT](kBT) §2。
- `-T`: T=0 のテトラヘドロンのまま、$\mathrm{Im}\chi_0$ を周波数方向に Gaussian で均す。幅は温度 $T$ の Fermi-Dirac と同じ(標準偏差 $\pi k_BT/\sqrt3$、$-1000$ K = 0.00574 Ha = 0.156 eV)。金属で $W$ のプラズモン極が荒れるときに使う(LiTi₂O₄ の実用の設定は $-992.4$。[kBT](kBT) の「実用の設定」)。$W(i\omega)$ や静的 $W(0)$ にも効く。
gwinit のテンプレートは 300 を書く。
- **type**: float (K)
- **default**: 無し(必須)

古い入力に出てくるキーの扱いは表 1 のとおり。`ctrlg_update.py ctrlg.<sname>.toml` は `SmearX0`(Ha、$\mathrm{Im}\chi_0$ の
Gaussian の標準偏差)を同じ幅の `t_tetrakbt = -T` に書き換え(0.0057 Ha → $-992.4$ K)、`t_tetrakbt` の無いファイルには
`t_tetrakbt = 0` を足す。計算は変わらない。`Legacy2toml.py` も変換のときに同じ書き換えをする。キーの変遷は [kBT](kBT) §2.5 の表 2。

**表 1**. `[gw]` の古いキーと今の扱い

| 古いキー | 今の扱い |
|---|---|
| `SmearX0`、`SmearX0q0` | 書いてあれば止まる。`t_tetrakbt = -T` を書く |
| `tetrakbt`(論理値) | 読まれない。有限温度は `t_tetrakbt > 0` |
| `esmr`(Ry、Gaussian) | `t_sigmaw` が無ければ、同じ標準偏差の `t_sigmaw` に換算して読まれる |
| `t_sigmakbt`(K) | `t_sigmaw` が無ければ、値がそのまま `t_sigmaw` になる |
| `GaussianFilterX0`、`GaussSmear`、`delta`、`dw`、`omg_c`、`WgtQ0P`、`MEMnmbatch` | 読まれない(書いてあっても効かない) |
| `zmel_max_size` | `zmel_batch_gb` として読まれる |

## `t_sigmaw`  *(2026-09-20、旧 `esmr` / `t_sigmakbt`)*
自己エネルギー $\Sigma = G(W)$ の中間準位(G の極)を均す Fermi-Dirac 核の幅。単位 K。
$\Sigma_x = Gv$ の占有重み、$\Sigma_c = G(W-v)$ の実軸極項の窓重みと極の位置(`wcsmear`)の
すべてに同じ核 $-\partial f/\partial\varepsilon$($k_BT$ = `t_sigmaw`)が使われる。
**これは数値的な幅であって完全な有限温度自己エネルギーではない**(Bose 因子などは入らない)。
Fermi 準位は `t_tetrakbt > 0`($\chi_0$ も有限温度)なら `EFERMI_kbt`、そうでなければ k メッシュの上の準位を数えて決める([kBT](kBT) §0 の表 5)。
中間状態の候補窓と実軸 $\omega$ メッシュの余裕は $E_F \pm 15\,k_BT$。同じ核が $\Sigma_c$ の虚軸積分の
半留数の段差にも入る(契約: 占有重み・極項の窓重み・虚軸の段差の 3 か所が同じ核。[manual/kBT](kBT) §3.6)。

- 既定 1000 K は旧 `esmr = 0.01` Ry(Gaussian $\sigma$)に近い幅(同じ標準偏差になるのは 872 K。FD の std は
  $\pi k_BT/\sqrt3 = 1.81\,k_BT$)。旧 `esmr = 0.003` Ry(Si など TestInstall の値)は 262 K に相当。
- 絶縁体: 幅はバンドギャップより十分小さければ結果に効かない(Si で 262 K と 1000 K の差は数 meV)。
- 金属: 原理的には $n_1n_2n_3 \to \infty$、`t_sigmaw` $\to 0$ の極限だが、有限の幅が収束に要る。
  細かい k メッシュと合わせて 300〜1000 K。物質ごとに選ぶ(LiTi₂O₄ は 1000 K)。
- `esmr` が ctrlg に残っていれば `t_sigmaw = esmr[Ry]·13.6057/(1.81·8.617e-5)` に換算して読む
  (`t_sigmaw` があればそちらが優先)。`t_sigmakbt` は値をそのまま `t_sigmaw` に写す。
- **type**: float (K)
- **default**: 1000.0(gwinit のテンプレートは 300 を書く)

### `wcsmear`  *(2026-09-19、2026-09-20 から既定 true)*
$\Sigma_c$ の実軸極項で、中間準位の Fermi-Dirac 核(`t_sigmaw`)を極の位置 $\omega_\varepsilon$ にも
適用して $W_c(\omega)$ を核で積分する。`false` で従来の平均エネルギー 1 点評価(3 点 Lagrange 内挿、
PRB76,165106 Eq.58)。$W_c$ が滑らかな系では両者はほぼ同じ(Si で数 meV〜最大 77 meV)。
鋭いプラズモン極を持つ金属(LiTi₂O₄ ≲ 1500 K)では必須。実軸極項のコストが数倍(全体で +7〜26 %)。
- **type**: boolean
- **default**: true

> テトラヘドロン重み(`tetwt5`)の並列: 2026-09-22 から **MPI の k 並列**(W-build の各ランクが自分の k 範囲だけ計算)。
> hgw を `-np2` で GPU 数より多いランクで起こす(GPU は round-robin で共有)とその分だけ速くなる。χ₀ を有限温度にする
> (`t_tetrakbt > 0`、9³ で 1 q 150 s)ときに効く。09-19〜22 に試した OpenMP(`omp_tetwt`)は削除。

### `unit_2pioa`
  `unit_2pioa = false`   # false --> a.u.; true --> unit of QpGcut_* are in 2*pi/alat

  普通 `false` でよい.

###  `alpha_OffG`
  `alpha_OffG = 1.0`   # (a.u.) Used in auxially function in the offset-Gamma method.

  Not need to change.

### `QpGcut_psi`

波動関数のIPW基底のカットオフ~~エネルギー: 単位 Ryd~~。 `lmf`計算で用いるPMT手法とは波動関数の基底関数が異なることに注意。原子間領域の波動関数はIPWで展開し直して記述される。IPWとはinterstitial plane wave.
- **type**: float
- **default** : 4.0
`|q+G| < QpGcut_psi` 単位 1/a.u.
個数はqg4gwの結果（lqg4gw)にかかれている。
エネルギー(Ry)にするには**2する。

波動関数のIPWのカットオフは `[ham] pwemax` とも関係する。pwemaxとしては2,3を試すことが多い
(こちらの単位は**2のエネルギーカットオフになっていて単位Ryd)。

### `QpGcut_cou`
積基底のIPW部分のカットオフ~~エネルギー: 単位 Ryd~~ 
`|q+G| < QpGcut_cou` 単位 1/a.u.
エネルギー(Ry)にするには**2する。
- **type**: float
- **default** : 3.0

個数はqg4gwの結果（lqg4gw)にかかれている。

### `emax_sigm`
フェルミ準位から測った, 計算する準粒子状態の自己エネルギーの最大値: 単位 Ryd。QSGW計算では自己エネルギーはポテンシャル構築に使用されるため小さくするとポテンシャルの精度も下がる。
- **type**: float
- **default** : 3.0


## Change `HistBin` for high resolution energy mesh near Ef for metal 
This defines histogram bins to accumulate weight (imaginary part) of the polarization functions.

For metal, it might be better to test
```toml
HistBin_dw    = 1e-5   # 1e-5 is fine mesh (good for metal?)  (a.u.) BinWidth along real axis at omega=0.
HistBin_ratio = 1.03   # 1.03 maybe safer. frhis(iw)= b*(exp(a*(iw-1))-1), where a=ratio-1.0 and dw=b*a
```
m_freq.f90のgetfreqでこのメッシュを作っている。

* HistBin_dw and HistBin_ratio specify real space bins which we accumulate imaginary part weight of polarization functions. The bins are (see frhis in m_freq.f90) 
$\omega_i = {{dw}} ∗ (\exp({\rm ratio} ∗ (i − 1)) − 1) $
where $[\omega_i, \omega_{i+1}]$ $(i = 1, 2, ...,\verb#nwhis#+1)$ is the i-th bin. HistBin_dw is bin width at ω = 0.
The ratio ω(i + 1)/ω(i) for large ω is exp(a)=HistBin_ratio. This choice of getting coarser at high energy is because we think W (ω) around ω ∼ 0 gives most important contribution to the GW approximation. If histogram bins are too wide, dielectric function can be less accurate, but results may be not so much affected.
In GW calculation, the Plasmon pole is important. It is determined not by the Drude weight; it usually gives small contribution to the Plasmon pole (for example, Si is described well by a Plasmon pole model, but Si has no Fermi surface). We expect that the GW results are not so sensitive to the choice of HistBin_dw, HistBin_ratio usually. We may use fine mesh when we plot quantities such as W (ω) near ω = 0.
<!-- The ecalj gives ¯W (ω = 0) ∼ 0 for metal; where ¯W (ω) is the effective interaction averaged in
the Γ-cell [?]. And ¯W (ω get closer to v for larger ω. -->

## `niw`
 Number of frequencies along Im axis. Used for integration along imaginary axis. Probably 10 is goog enough. 
See routines `wintzsg.f90`. 

The integration points are $i \omega'(n)= i( 1/x(n) -1)$, where $x(n)$ is the
usual Gaussian-integration points for the interval $[0,1]$. In addition, we give the special analytical treatment for the peaky part at $\omega'= 0$. Out tests shows niw=6 for Si is good enough for 0.01 eV accuracy. The convergence as for niw is quite good. This integration scheme has been developed by Ferdi Aryasetiawan. The number of points should be the one of 6,10,12,16,20,24,32,40,or 48, because we use a subroutine gauss in mate.f90 prepared by Ferdi.

## Q0PChoice
`Q0Pchoice` sets the size of the offset-Gamma points Q0P (the factor `deltaq_scale` on half of the mesh spacing).
`Q0Pchoice = 1` (default) gives `deltaq_scale` = 0.1 (close to the $q\to0$ limit), `Q0Pchoice = 2` gives $1/\sqrt3$ (mean value in a sphere).
`deltaq_scale = x` ($x>0$) in `[gw]` sets the factor directly and overrides `Q0Pchoice`. See lqg4gw for the Q0P in use.

### `deltaw`
real (a.u.) only for one-shot case.
deltaw is the interval for the numerical derivative ∂Σ(ω)/∂ω.  We calculate $\Sigma(\omega),\Sigma(\omega\pm {\rm deltaw})$.
From these values, we can calculate two Z (or second-derivative of Σ(ω)). 


## `<QforGW>` section
This is only for one-shot/spectrum function mode. 
In `[gw]`, `QforGW` gives the q vectors, one per line between `"""` (in the unit of `2pi/alat`. See BZ.html and syml file to check):
```toml
QforGW = """
 0.0 0.0 0.0
 0.5 0.0 0.0
"""
```
If no `QforGW` (or `QforGWIBZ = true`), all irreducible q points are used.

## EMAXforGW
one-shot GW. Set this (eV,above the Fermi energy) to specify up to which bands you calculate. 
(we have EMINforGW).


## `--readQforGW` option
We can use any q in the QforGW section
Because of the shifted-mesh method, we can use any q point which is not on the mesh points.
If on the mesh point, we need to prepare less number of eigenfunctions for GW calculations.

## `<QforEPS>`
This is for dielectric function mode. Set in `[gw]` as
```toml
QforEPS = """
 0 0 0.00050
 0 0 0.00100
 0 0 0.00200
"""
```
Because of the shifted-mesh method, we can use any q point which is not on the mesh points.
### QforEPSau 
Whether q is in the unit of `2pi/alat` or in the unit of a.u.(=`bohr^{-1}`). I think `QforEPSau = true` might be convenient.



## `Smaller lcutmx to reduce computatinal time`
  To reduce the computational time, we reduce number of MPB (mixed product basis).
  One is lcutoff of Product Basis (PB) within MT. Use lcutmx=2 for oxygen or something (s,p block atoms). 
  Thus it  is like
```toml
[product_basis]
pb_lcutmx = [4, 4, 4, 2, 2, 2]   # maximum l-cutoff for the product basis, one per atom. 4 is required for atoms with valence d
```
in `ctrlg.<sname>.toml`. The atoms are in the order of the `[[site]]` tables (`lmchk` shows it, also in `SiteInfo.lmchk`).

* NOTE: we need lcutmx =6 is requied for 4f/5f atoms, while 2 is enough for Oxygen (probably other species on the same row).


## `Smaller IPW related part to reduce computational time`
To reduce compuational time, we may reduce MPB coming from IPS, reduce the size of IPW for psi, reduce emax_sigm and pwemax as follows.
   + QpGcut_cou is fort the Interstitial plane wave (IPW) for MPB. 
   + QpGcut_psi is for expantion of eigenfunctions.
   + emax_sigm is the upper cutoff (relative to the Fermi energy) to calculate self energy.
   + `[ham] pwemax` is the APW basis cutoff for the eigenfunciton.

To reduce computational time, we may use
```toml
[gw]
QpGcut_psi = 3.0
QpGcut_cou = 2.5
emax_sigm  = 2.0
```
In addition, we may use
```toml
[ham]
pwemax = 2
```
This is for the APWs in the band calculation.

We sometimes use this setting. It is better to check how the numerical results are affected. 
