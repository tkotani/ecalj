# 実装メモ（GW の実装の詳細） {#dev-impl}

::: info このページについて
このページは、開発者向けの古いメモ `ecaljdetails.tex`（T. Kotani、最終版 2024-11。原本は ecalj のコミット `3284d0b` の `ecaljdetails/`、`git show 3284d0b:ecaljdetails/ecaljdetails.tex`）から、
いまも役に立つ節を写し、いまのコードと M. Obata et al. (2025) に合わせて見直したものです。
式番号 (I1)〜 は写したときに付け直したもので、消した式の番号は欠番です。`[key]` は文献のキー（同じコミットの `ecaljdetails/ecaljrefs_detail.bib`）です。
理論の説明は [ecalj のバンド計算法](../theory/pmt)、[QSGW](../theory/qsgw) と、そこにある文献を見てください。数学の補足は [付録](../theory/appendix) に移しました。
:::


GW の部分の式の主な文献は [Obata2025] である: M. Obata, T. Kotani, T. Oda, "Efficient implementation of the quasiparticle self-consistent $GW$ method on GPU", Comput. Mater. Sci. **260**, 114190 (2025), [doi:10.1016/j.commatsci.2025.114190](https://doi.org/10.1016/j.commatsci.2025.114190), 第 2 節。以下で「[Obata2025] Eq. (n)」と付けた式はそこから取ったものである。

## 固有関数の表し方

PMT 法 [kotani_fusion_2010] では、
与えられた $H^0$ の価電子の固有関数を、
Bloch 和をとった MTO
$\chi^{ {\bf k} }_{ {\bf R}L{j} }({\bf r})$ と APW $\chi^{ {\bf k} }_{\bf G}({\bf r})$ の線形結合で表す。
$$\begin{aligned}
{\Psi_{ {\bf k}n} }({\bf r}) = \sum_{ {\bf R}L{j} } z^{ {{\bf k} }n}_{ {\bf R}L{j} }
\chi^{ {\bf k} }_{ {\bf R}L{j} }({ {\bf r} })+ \sum_{ {\bf G} } z^{ {\bf k}n}_{ {\bf G} }
\chi^{ {\bf k} }_{\bf G}({\bf r}),
\end{aligned}
\tag{I1}
$$
ここで、添字は波数ベクトル ${\bf k}$、バンドの添字 $n$、
逆格子ベクトル ${\bf G}$ である。基本セルの中の MTO は、
MT サイトの添字 ${\bf R}$、
角運動量 $L=(l,m)$、動径関数の添字 $j$ で
指定する。
コアの固有関数は、MT の中に閉じこめるという条件の下で
計算する。
そのうえで、コアの寄与を、後の式 (I15) で定める交換の部分には入れる。相関の部分には入れない。
コアはすべて Ref.[kotani07a] の「core1」として扱う（`[product_basis]` の `core` の行はすべて 0）。

###  MTO の部分

MT の中では、MTO の Bloch 和
$\chi^{ {\bf k} }_{ {\bf R}L{j} }({\bf r})$ は、
原子に似た軌道
$A_{ {\bf R}{Lj} }({\bf r})
\equiv \{ 
\phi_{ {\bf R}{Lj} }(r),\dot{\phi}_{ {\bf R}{Lj} }(r),\phi^z_{ {\bf R}{Lj} }(r) \}
\times Y_L$ の線形結合で表す（$\phi^z$ は局在軌道（local orbital）を表す）。
これらの動径関数は、${\bf R}$ の中での
動径 Schrödinger 方程式の解（またはそのエネルギー微分）である。
MTO 基底は smooth Hankel 関数で指定し、
smooth Hankel 関数は 2 つのパラメータ（$E=-|\kappa|^2,R_{\rm sm}$）をもつ。

$A_{ {\bf R}{Lj} }({\bf r})$ は各 MT ${\bf R}$ で正規直交基底をなす。
すると、すその部分を含めた MTO は次のように書ける。
$$\begin{aligned}
\chi^{ {\bf k} }_{ {\bf R}{Lj} }({\bf r}) & =  \sum_{ {\bf R}{Lj} } C^{ {\bf k} }_{ {\bf R}{Lj} } A^{\bf k}_{ {\bf R}{Lj} }({\bf r})   \ \ \ {\rm \ if \ {\bf r} \in any \ MT}  \\
        & =    H^{\kappa,R_s,{\bf k} }_{ {\bf R}{L} }({\bf r}) \ \ \ {\rm otherwise},
\end{aligned}
\tag{I2}
$$
ここで Bloch 和は次のとおり。
$$\begin{aligned}
A^{\bf k}_{ {\bf R}{Lj} }({\bf r}) & \equiv  \sum_{\bf T} A_{ {\bf R}{Lj} }({\bf r-{\bf R}-{\bf T} }) \exp(i {\bf kT}), \\
H^{ {\bf k}s}({\bf r})   & \equiv  \sum_{\bf T} H_{s}({\bf r-{\bf R}-{\bf T} })
 \exp(i {\bf kT}).
\end{aligned}
\tag{I3}
$$
ここで smooth Hankel 関数 $H^{ {\bf k}s}({\bf r})$ は
MTO の包絡関数（envelope function）である。

MTO は角運動量 $L$ の対称性をもつ。
原子に似た局在関数はどれも同じである。
そうした関数はすべて、空間群の操作で同じ写像の式に従う。
節（[§](#spacegroup-rotation)）を見よ。

###  APW の部分

APW $\chi^{ {\bf k} }_{\bf G}({\bf r})$ は、
MT の中では原子に似た軌道
$A_{ {\bf R}{Lu} }({\bf r})\equiv 
\{ \phi_{ {\bf R}{l}u}(r) Y_L(\hat{\bf r}),
    \dot{\phi}_{ {\bf R}{l}u}(r) Y_L(\hat{\bf r}) \}$
の線形結合で、
格子間領域ではただの平面波で与える。ここで
$\phi_{ {\bf R}{l}u}(r)$ と $\dot{\phi}_{ {\bf R}{l}u}(r)$
は、MT 半径での対数微分を与えたときの、各 $l$ の
動径 Schrödinger 方程式の 2 つの解を表す。
`[ham] readp = true` と `pnufix = true`（既定）では、対数微分は `lmfa` の自由原子のもの（`atmpnu.*`）で、固定したままにする
（`pnufix = false` では、`pnunew.f90` が占有された PDOS の重心に置く。[lmf § readp](./lmf#readp-option)）。$\dot{\phi}$ はエネルギー微分
（またはそれに似たもの）を表す。
$u$ は $\phi$ と $\dot{\phi}$ を区別する
複合添字である。
${\bf R}$ は基本セルの中の MT を指定する添字である。
APW 基底は $s \equiv { {\bf R}{j}L}$ で指定する。ここで
$L\equiv(l,m)$ は角運動量の添字、
$j$ は追加の添字（主量子数など）である。
$A_{ {\bf R}{L}u}({\bf r})$ は各 MT ${\bf R}$ で
正規直交基底をなす。
APW は次のように書ける。
$$\begin{aligned}
\chi^{\bf k+G}({\bf r}) & =  \sum_{a u} C^{\bf k+G}_{a u} A^{\bf k}_{a u}({\bf r})   \ \ \ {\rm \ if \ {\bf r} \in any \ MT}  \\
        & =    \exp(i ({\bf k+G}){\bf r}) \ \ \ {\rm otherwise},
\end{aligned}
\tag{I4}
$$
ここで Bloch 和は次のとおり。
$$\begin{aligned}
A^{\bf k}_{ {\bf R}u}({\bf r}) & \equiv  \sum_{\bf T} 
A_{ {\bf R}u}({\bf r-R-T}) \exp(i {\bf kT}),
\end{aligned}
\tag{I5}
$$

${\bf G}$ の数は条件
$|{\bf k+G}|\lt  {\tt QpGcut\_psi}$（IPWpsi）で制限する。
係数 $\alpha^{ {\bf k}n}_{au}$ は次のように計算できる。
$$\begin{aligned}
\alpha^{ {\bf k}n}_{au} = \sum_{\bf G} C^{ {\bf k+G} }_{a u} z^{\bf k+G}_n.
\end{aligned}
\tag{I6}
$$

## 固有関数の展開し直し: CPHI と GEIG

$GW$ 計算をするには、まず
与えた BZ のメッシュについて、すべての固有関数（と固有値）を
用意しなければならない。そのうえで固有関数を
次のように表す。
式 (I1) の ${\Psi_{ {\bf k}n} }({\bf r})$ を、
MT の中の augmentation の部分と、格子間領域の
PW の部分との和として展開し直す。
$$\begin{aligned}
{\Psi_{ {\bf k}n} }({\bf r})
= \sum_{ {\bf R}u}  \alpha^{ {{\bf k} }n}_{ {\bf R}u} \varphi^{\bf k}_{ {\bf R}u}({\bf r})
 + \sum_{\bf G}  \beta^{ {{\bf k} }n}_{\bf G} P^{\bf k}_{\bf G}({\bf r}),
\end{aligned}
\tag{I7}
$$
ここで格子間平面波（IPW）を次で定める。
$$\begin{aligned}
P^{\bf k}_{\bf G}({\bf r}) =
\begin{cases}
0    & \text{if {\bf r} } \in \text{any MT} \\
\exp(i ({\bf k+G})\cdot{\bf r})& \text{otherwise}
\end{cases}
\end{aligned}
\tag{I8}
$$
また $\varphi^{\bf k}_{R u}({\bf r})$ は、$R$ の MT の中で定めた原子関数
$\varphi_{R u}({\bf r})$ の Bloch 和である。
$$\begin{aligned}
\varphi^{\bf k}_{R u}({\bf r}) & \equiv  \sum_{\bf T} \varphi_{R u}({\bf r-R-T}) \exp(i {\bf k\cdot{}T}).
\end{aligned}
\tag{I9}
$$
**T** と **G** は、それぞれ実空間と逆空間の格子並進ベクトルである。
式 (I7) は Ref.[kotani_quasiparticle_2014] の式 (17) である。
`lmf --jobgw=1`（`sugw.f90`）は、$\alpha_{ {\bf R}u}^{ {\bf k}n}$ を `__CPHI` に、$\beta_{\bf G}^{ {\bf k}n}$ を `__GEIG` に、固有値を `__EValue` に書く。これらは `readeigen.f90` の `readcphi` と `readgeig` が読む（たとえば `m_zmel.f90` で使う）。

$\Psi_{ {\bf k}n}$ は、**qg4gw** が作る ${\bf k}$ 点（正則メッシュと offset-Γ 点）で要る。
すると $P_{\bf G}^{\bf k}({\bf r})$ は
${\bf G}$ だけで指定される。これも **qg4gw** が作る。
$\varphi_{ {\bf R}u}^{\bf k}({\bf r})$ は
動径関数で指定する。動径関数は `__PHIVC` にあり、`sugw.f90` が書き、`main_hbasfp0.f90` で読む。動径関数の数は
`ncore(ic)+nrad(ic)` である
（$l,n,\sigma$ には依存するが、$m$ には依存しない）。
`lmxa` は種ごとに与え（`[[spec]] lmxa`）、異なってよい。配列は最大値 `lmxax` で大きさを決める（`sugw.f90`）。

ecalj の GW 計算では、固有関数は行列要素（[Obata2025] Eq. (8)）
$$
Z^{\mathbf q-\mathbf k n',\mathbf q n}_{\mathbf k\mu} = \langle E^{\mathbf k}_\mu \Psi_{\mathbf q-\mathbf k n'} | \Psi_{\mathbf q n}\rangle ,
\tag{I11}
$$
を通してだけ入る。
ここで $E^{\mathbf k}_\mu$ は、正規直交化し Coulomb 相互作用を対角にした MPB である（後の式 (I14)）。
分極関数では、ラベルを付け替えた $Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\mu}=\langle E^{\mathbf q}_\mu\Psi_{\mathbf k n}|\Psi_{\mathbf q+\mathbf k n'}\rangle$ として使う。
これらはその場で `build_zmel` が作り、`m_zmel.f90` の `zmel` に置く（[GW で使う行列要素](../theory/gw_zmel)）。

式 (I7) の係数（ここでは MTO の部分だけ）は
次のように計算する。
$$\begin{aligned}
 &\alpha^{ {\bf k} n}_{au} = \sum_s C^{ {\bf k}s}_{a u} z^{ {\bf k}n}_s \\
 &\beta^{ {\bf k}n}_{\bf G} = \sum_{ {\bf G}'s} 
\langle P^{\bf k}_{\bf G}|P^{\bf k}_{ {\bf G}'}\rangle^{-1}
\langle P^{\bf k}_{ {\bf G}'}|H^{ {\bf k}s}\rangle z^{ {\bf k}n}_s,
\end{aligned}
\tag{I12}
$$
ここで ${\bf G}$ の数は条件
$|{\bf k+G}|\lt  {\tt QpGcut\_psi}$ で制限する。
${\bf G}'$ は `lmf` の `gmax` で制限する（重なりは `sugw.f90` から呼ぶ `pwmat` が作る）。
基底 $s$ には MTO だけでなく APW も入る。そのうえで係数を Gram-Schmidt で正規直交化する（`--skipGS` で省く）。

## GW 計算の概要

$GW$ 計算では、固有関数の基底だけでなく、
固有関数の積を展開する基底も要る。
この基底を混合積基底（mixed product basis、MPB）$\{M^{\bf k}_I({\bf r}) \}$ と呼ぶ。Kotani が Ref.[kotani_all-electron_2002] で初めて導入した。MPB は、MT の中の積基底（PB）[aryasetiawan_product-basis_1994]
と格子間領域の IPW からなる。
$\{M^{\bf k}_I({\bf r}) \}$ は直交しない IPW を含むので、重なり行列と双対基底をもつ。
$$\begin{aligned}
 & O^{\bf k}_{IJ} = \langle M^{\bf k}_{I} |  M^{\bf k}_J \rangle, \\
 & |\tilde{M}^{\bf k}_{I} \rangle \equiv \sum_{J}
   |M^{\bf k}_{J} \rangle (O^{\bf k})^{-1}_{JI} .
\end{aligned}
\tag{I13}
$$

$v_{IJ}^{\bf k}= \langle M^{\bf k}_{I} |v|  M^{\bf k}_J \rangle$ として、一般化固有値問題
$\sum_J (v_{IJ}^{\bf k}- v_\mu({\bf k}) O^{\bf k}_{IJ} ) z^{\bf k}_{\mu J}= 0$
は、固有値 $v_\mu({\bf k})$ と、正規直交な MPB
$|E^{\bf k}_\mu\rangle=\sum_J z^{\bf k}_{\mu J}|M^{\bf k}_J\rangle$ の係数 $z^{\bf k}_{\mu J}$（コードでは `zcousq`）を与える。この基底は Coulomb 相互作用を対角にする（[Obata2025] Eq. (7)）。
$$\begin{aligned}
v({\bf k})=\sum_{\mu} | E^{ {\bf k} }_\mu \rangle {v_\mu({\bf k})}
\langle E^{ {\bf k} }_\mu |.
\end{aligned}
\tag{I14}
$$

**コードの中の 2 つの基底。** コードは MPB $M^{\bf k}_I$ と正規直交な $E^{\bf k}_\mu$ の両方をもち、$z^{\bf k}$ で互いに移す（`m_zmel.f90` の `m2e_prod_basis`。`set_m2e_prod_basis` が `zcousq` から作る）。
`build_zmel` はまず $\langle M^{\bf k}_I\Psi_{\mathbf q-\mathbf k n'}|\Psi_{\mathbf q n}\rangle$（`zmelt`）を作り、`is_m_basis=.true.` ではそのまま返し、
`is_m_basis=.false.` では式 (I11) の $Z^{\mathbf q-\mathbf k n',\mathbf q n}_{\mathbf k\mu}=\sum_I (z^{\bf k}_{\mu I})^*\langle M^{\bf k}_I\Psi_{\mathbf q-\mathbf k n'}|\Psi_{\mathbf q n}\rangle$ として返す。
`hgw`（`main_hgw.f90`）では、$\mathbf k$ とバンドについての大きな和を M 基底で行い、各 $\mathbf q$ と $\omega$ で $\mathrm{ngb}\times\mathrm{ngb}$ の行列だけを変換する。
$$
\begin{aligned}
\Pi^{E}_{\mu\nu} &= \sum_{IJ} (z_{\mu I})^*\,\Pi^{M}_{IJ}\, z_{\nu J}, \\
W^{{\rm c},M}_{IJ} &= \sum_{\mu\nu} z_{\mu I}\, W^{{\rm c},E}_{\mu\nu}\, (z_{\nu J})^*, \\
W^{\rm c} &= \sum_{IJ}|M_I\rangle W^{{\rm c},M}_{IJ}\langle M_J| .
\end{aligned}
\tag{I14a}
$$
- $\Pi^M$ は `zmelt` から積み上げ（`m_screened_coulomb.f90` から呼ぶ `x0kf_zxq`、`is_m_basis=.true.`）、$\Pi^E$ に直し、式 (I16) と (I19) を E 基底で計算する（`m_llw.f90` の `WVRllwR`、`WVIllwI`。`is_x0_m_basis`、`is_wc_m_basis` は true）。
  そのうえで $W^{\rm c}$ を M 基底で置く。$\mathbf q=0$ では、M 基底への変換を offset-$\Gamma$ の補正の後に行う（`w0w0i`、`is_wc_m_basis=.true.`）。
- $\Sigma^{\rm c}$（式 (I20)）は、M 基底の `zmelt` と $W^{{\rm c},M}$ を使う（`m_sxcf_sc.f90`、`build_zmel(..., is_m_basis=.true.)`）。
- $\Sigma^{\rm x}$（式 (I15)）は、$v$ が対角になる E 基底を使う（`build_zmel(..., is_m_basis=.false.)`）。コアの状態（コアの交換）では、MPB の各原子のブロックを別々に変換する。コアの状態の `zmelt` はその MT の外では 0 だからである。
- `hx0fp0`（`main_hx0fp0.f90`）は終始 E 基底で計算する（`is_x0_m_basis=.false.`）。

論文 [Obata2025] はすべてを E 基底で書いている。式 (I14a) は同じ $\Pi^E$、$W^{\rm c}$、$\Sigma$ を与える。

全電子フルポテンシャルの $GW$ 近似では、
式 (I14) は Ref.[friedrich_efficient_2010] で導入された。
これは平面波展開での表現
$v({\bf k}+{\bf G},{\bf k}+{\bf G}')=\frac{4 \pi \delta_{ {\bf G}{\bf G}'} }{|{\bf k}+{\bf G}|^2}$ に対応する。
$\mu=1$ は $v_{\mu}$ の最大の固有値に対応し、
$v_{\mu=1}({\bf k})$ は $\sim \frac{4 \pi}{|{\bf k}|^2}$ である。これは
節（[§](#offset-gamma)）で述べる発散項と関係する。

$\langle A| B\rangle =\int d^3r A^*({\bf r}) B({\bf r})$ として、
交換の自己エネルギーは次である（[Obata2025] Eq. (9)）。
$$\begin{aligned}
\Sigma^{\rm x}_{nm}({\bf q})&=
\langle \Psi_{ {\bf q}n} |\Sigma^{\rm x} |\Psi_{ {\bf q}m} \rangle \\
&=-\sum^{\rm BZ}_{ {\bf k} }  \sum^{\rm  occ}_{n'} \sum_\mu
[Z^{\mathbf q-\mathbf k n',\mathbf q n}_{\mathbf k\mu}]^*\, v_{\mu}({\bf k})\,
Z^{\mathbf q-\mathbf k n',\mathbf q m}_{\mathbf k\mu}.
\end{aligned}
\tag{I15}
$$
占有された価電子の状態 $\Psi_{\mathbf q-\mathbf k n'}$ をコアの状態に置きかえると、コアの交換 $\Sigma^{\rm Cx}$ になる。
その MPB（MPB$_{\rm C}$。ここでは $M_{{\rm C}I}$、$E_{{\rm C}\mu}$ と書く）は価電子のものとは異なる。コアと価電子の固有関数の積を表すからである。
`gwsc` はこれを価電子の部分の前の別の段で作る。`hbasfp0 --job=3` が積基底 $M_{\rm C}$ を作り（`genallcf_v3`、`incwfin=-2`: コアの状態と価電子の動径関数の積）、
`hvccfp0 --job=3` がその基底での Coulomb 行列とその固有ベクトル（$E_{\rm C}$、`__Vcoud.*` に置く。後で価電子の `hvccfp0 --job=0` が上書きする）を作り、
`hsfp0_sc --job=3` がコアの状態で式 (I15) を $E_{\rm C}$ 基底で計算する。
このように `gwsc` は各反復で交換を 2 回計算する。ここでのコアの交換 $\Sigma^{\rm Cx}$（$E_{\rm C}$ で）と、`hbasfp0 --job=0` と `hvccfp0 --job=0` の後に `hgw` の中で計算する価電子の交換 $\Sigma^{\rm x}$（$E$ で）である。
コアの状態は MT の中にあるので、コアの交換は積基底を使わずに計算できる（今後の課題）。

遮蔽された Coulomb 相互作用は対称な形で与える（[Obata2025] Eqs. (10), (11)）。
$$
\begin{aligned}
\epsilon_{\mu\nu}(\mathbf q,\omega) &= \delta_{\mu\nu} - \sqrt{v_\mu(\mathbf q)}\,\Pi_{\mu\nu}(\mathbf q,\omega)\sqrt{v_\nu(\mathbf q)}, \\
W_{\mu\nu}(\mathbf q,\omega) &= \sqrt{v_\mu(\mathbf q)}\,\epsilon^{-1}_{\mu\nu}(\mathbf q,\omega)\sqrt{v_\nu(\mathbf q)},
\end{aligned}
\tag{I16}
$$
ここで Lindhard の分極関数は次である（[Obata2025] Eq. (12)）。
$$\begin{aligned}
\Pi_{\mu \nu}({\bf q},\omega)
 &=
\sum^{\rm BZ}_{ {\bf k} } \sum^{\rm occ}_{n } \sum^{\rm unocc}_{n'}
\frac{Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\mu}[Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\nu}]^*}
{\omega-(\varepsilon_{ {\bf q+k} n'}-\varepsilon_{ {\bf k}n})+i \delta} \\
 &+ \sum^{\rm BZ}_{ {\bf k} } \sum^{\rm  unocc}_{n } \sum^{\rm occ}_{n'}
\frac{Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\mu}[Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\nu}]^*}
{-\omega-(\varepsilon_{ {\bf k}n}-\varepsilon_{ {\bf q+k} n'})+i \delta}.
\end{aligned}
\tag{I17}
$$
時間反転対称性がある（$\Psi^*_{\mathbf k n}$ が $\Psi_{-\mathbf k n}$ と縮退した固有関数である）と、第 2 項は第 1 項で $n\leftrightarrow n'$ としたものに等しく、次になる（[Obata2025] Eq. (14)）。
$$\begin{aligned}
\Pi_{\mu \nu}({\bf q},\omega)
 &=\sum^{\rm BZ}_{ {\bf k} }  \sum^{\rm  occ}_{n} \sum^{\rm  unocc}_{n'}
Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\mu}[Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\nu}]^* \\
 &\quad\times\Big(\frac{1}{\omega-\varepsilon_{ {\bf q+k}n'}+\varepsilon_{ {\bf k}n}+i \delta} \\
 &\qquad -\frac{1}{\omega+\varepsilon_{ {\bf q+k}n'}-\varepsilon_{ {\bf k}n}-i \delta}\Big).
\end{aligned}
\tag{I18}
$$
この定義は Fetter-Walecka の教科書のものと同じである。
$\omega=0$ では負定値の行列で、$\Pi_{\mu\nu}(\mathbf q,\omega)=\Pi_{\mu\nu}(\mathbf q,-\omega)$ である。
時間反転対称性がないときは、式 (I17) の負エネルギーの部分も積み上げる（近似として式 (I18) を使ってもよい）。

まず $\Pi$ の反 Hermite の部分を、実軸のヒストグラムのビン $\omega_j=[\omega_{{\rm L}j},\omega_{{\rm R}j}]$ の上に
tetrahedron method [rath_generalized_1975] で積み上げ、次に Hilbert 変換で Hermite の部分を得る（[Obata2025] Eqs. (15)–(18)）。
$$
\begin{aligned}
&\mathrm{Im}[\Pi_{\mu\nu}(\mathbf q,\omega_j)] = \frac{\pi}{i}\sum^{\rm BZ}_{\mathbf k}\sum^{\rm occ}_{n}\sum^{\rm unocc}_{n'} \\
&\qquad\qquad Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\mu}\,[Z^{\mathbf k n,\mathbf q+\mathbf k n'}_{\mathbf q\nu}]^*\, w^{\mathbf q,\mathbf k}_{n'n}(\omega_j), \\
&w^{\mathbf q,\mathbf k}_{n'n}(\omega_j) = \sum_{\mathbf k\in{\rm BZtet}(i_{\rm tet})}\frac14\, w(n',n,\mathbf q,\mathbf k,i_{\rm tet},\omega_j), \\
&w(n',n,\mathbf q,\mathbf k,i_{\rm tet},\omega_j) \\
&\qquad\qquad = \int_{{\rm BZtet}(i_{\rm tet})} d^3k \int_{\omega_{{\rm L}j}}^{\omega_{{\rm R}j}} \delta(\omega-\varepsilon_{\mathbf q+\mathbf k n'}+\varepsilon_{\mathbf k n})\,d\omega, \\
&\mathrm{Re}[\Pi_{\mu\nu}(\mathbf q,\omega_i)] = \sum_j K_{ij}\,\mathrm{Im}[\Pi_{\mu\nu}(\mathbf q,\omega_j)].
\end{aligned}
\tag{I18a}
$$
ビンは対数メッシュで、$\omega=0$ の近くで密である（`[gw] HistBin_dw`、`HistBin_ratio`。[kotani_quasiparticle_2007]）。
$K_{ij}$ は Hilbert 変換の行列で、ビンの大きさを考慮して $1/(\omega_i-\omega_j)$ を離散化したものである。
tetrahedron の重みは `gettetwt`（`m_tetwt.f90`）が作る。
これは実部を直接計算するより効率がよく、安全である。

Ref.[friedrich_efficient_2010] の拡張既約領域（EIBZ）での対称化は、いまは使っていない。`eibzgen.f90` はあるが呼ばれず、
`[gw] EIBZmode` は読むが使わず、`main_hx0fp0.f90` は `ngrpx=1`（対称操作なし）で走る。将来また使うかもしれない（[IBZ と EIBZ](#ibz-and-eibz-scheme)）。

相関の部分 $W^{\rm c}=W-v$ は次である（[Obata2025] Eq. (19)）。
$$
W^{\rm c}_{\mu\nu}(\mathbf k,\omega) = \sqrt{v_\mu(\mathbf k)}\,(\epsilon^{-1}_{\mu\nu}(\mathbf k,\omega)-\delta_{\mu\nu})\sqrt{v_\nu(\mathbf k)},
\tag{I19}
$$
相関の自己エネルギーは次である（[Obata2025] Eq. (20)）。
$$
\begin{aligned}
\Sigma^{\rm c}_{nm}({\bf q},\omega) &= \sum_{\mathbf k,n'}\sum_{\mu\nu} \int_{-\infty}^\infty \frac{i\,d\omega'}{2\pi}\,
\frac{e^{-i\delta\omega'}}{\omega-\omega'-\varepsilon_{\mathbf q-\mathbf k n'}\pm i \delta} \\
&\qquad\times [Z^{\mathbf q-\mathbf k n',\mathbf q n}_{\mathbf k\mu}]^*\,W^{\rm c}_{\mu\nu}(\mathbf k,\omega')\,Z^{\mathbf q-\mathbf k n',\mathbf q m}_{\mathbf k\nu},
\end{aligned}
\tag{I20}
$$
占有状態 $\mathbf q-\mathbf k n'$ には $-i\delta$、非占有状態には $+i\delta$ をとる。
QSGW では、[Obata2025] Eq. (3)（[QSGW](../theory/qsgw) の式 (Q3)）で $V^{\rm xc}$ を得るために、$\Sigma^{\rm c}_{nm}(\mathbf q,\varepsilon_{\mathbf q n})$ の Hermite の部分が要る。
$\omega'$ の積分路を変形すると、$\Sigma^{\rm c}$ は虚軸の部分 $I^{\rm c}$ と極の部分 $R^{\rm c}$ の和になる（[Obata2025] Eqs. (21)–(24)）。
$$
\begin{aligned}
\Sigma^{\rm c}_{nm}(\mathbf q,\varepsilon_{\mathbf q n}) &= I^{\rm c}_{nm}(\mathbf q,\varepsilon_{\mathbf q n}) + R^{\rm c}_{nm}(\mathbf q,\varepsilon_{\mathbf q n}), \\
I^{\rm c}_{nm}(\mathbf q,\varepsilon_{\mathbf q n}) &= -\sum_{\mathbf k,n'}\sum_{\mu\nu}\int_{-\infty}^{\infty}\frac{d\omega}{2\pi}
 \frac{1}{\varepsilon_{\mathbf q n}-i\omega-\varepsilon_{\mathbf q-\mathbf k n'}} \\
 &\qquad\times [Z^{\mathbf q-\mathbf k n',\mathbf q n}_{\mathbf k\mu}]^*\,W^{\rm c}_{\mu\nu}(\mathbf k,i\omega)\,Z^{\mathbf q-\mathbf k n',\mathbf q m}_{\mathbf k\nu} \\
 &\simeq -\sum_{\mathbf k,n'}\sum_{i}\sum_{\mu\nu} w^{\rm I}_i(\varepsilon_{\mathbf q n}-\varepsilon_{\mathbf q-\mathbf k n'}) \\
 &\qquad\times [Z^{\mathbf q-\mathbf k n',\mathbf q n}_{\mathbf k\mu}]^*\,W^{\rm c}_{\mu\nu}(\mathbf k,i\omega^{\rm I}_i)\,Z^{\mathbf q-\mathbf k n',\mathbf q m}_{\mathbf k\nu}, \\
R^{\rm c}_{nm}(\mathbf q,\varepsilon_{\mathbf q n}) &= \sum_{\mathbf k,n'}\sum_{\mu\nu}[Z^{\mathbf q-\mathbf k n',\mathbf q n}_{\mathbf k\mu}]^* \\
 &\qquad\times W^{\rm c}_{\mu\nu}(\mathbf k,\varepsilon_{\mathbf q n}-\varepsilon_{\mathbf q-\mathbf k n'})\,Z^{\mathbf q-\mathbf k n',\mathbf q m}_{\mathbf k\nu} \\
 &\quad\times\big[\theta(\varepsilon_{\mathbf q n}-E_{\rm F})\,\theta(\varepsilon_{\mathbf q n}-\varepsilon_{\mathbf q-\mathbf k n'}) \\
 &\qquad\qquad\times\theta(\varepsilon_{\mathbf q-\mathbf k n'}-E_{\rm F}) \\
 &\qquad -\theta(E_{\rm F}-\varepsilon_{\mathbf q n})\,\theta(\varepsilon_{\mathbf q-\mathbf k n'}-\varepsilon_{\mathbf q n}) \\
 &\qquad\qquad\times\theta(E_{\rm F}-\varepsilon_{\mathbf q-\mathbf k n'})\big] \\
 &\simeq \sum_{\mathbf k,n'}\sum_{\omega^{\rm R}_i}\sum_{\mu\nu} w^{\rm R}_i(\varepsilon_{\mathbf q n}-E_{\rm F},\varepsilon_{\mathbf q-\mathbf k n'}-E_{\rm F}) \\
 &\qquad\times[Z^{\mathbf q-\mathbf k n',\mathbf q n}_{\mathbf k\mu}]^* \\
 &\qquad\times W^{\rm c}_{\mu\nu}(\mathbf k,\omega^{\rm R}_i)\,Z^{\mathbf q-\mathbf k n',\mathbf q m}_{\mathbf k\nu}.
\end{aligned}
\tag{I20a}
$$
虚軸の点 $i\omega^{\rm I}_i$ での重み $w^{\rm I}_i(\Delta E)$ は Aryasetiawan と Gunnarsson に従う。
$W^{\rm c}(\omega)=W^{\rm c}(-\omega)$ なので $\omega^{\rm I}_i\ge0$ だけを使い、ふつう 10 点ほどである。
$\omega^{\rm R}_i$ は式 (I18a) のヒストグラムのビンの境界で、$w^{\rm R}_i$ は $\varepsilon_{\mathbf q n}-\varepsilon_{\mathbf q-\mathbf k n'}$ での $W^{\rm c}$ の 3 点補間の係数（階段関数を含む）である。
寿命（衝突電離率、impact ionization rate）を与えるのは $R^{\rm c}$ だけである。QSGW の静的な $\Sigma^{\rm c}$ では、$R^{\rm c}$ の中の $W^{\rm c}$ を Hermite 化する。
いまのコードでは、階段関数を無限小の $\delta$ ではなく `[gw] t_sigmaw`（K）の Fermi-Dirac 核でなまらせ、
既定（`[gw] wcsmear = true`）では $W^{\rm c}(\omega)$ そのものをその核で積分する。[Obata2025] の 3 点補間は `wcsmear = false` である（[kBT](./kBT)）。核での積分（`wcsmear = true`）は [Obata2025] の後に入れた改良で、$W^{\rm c}(\omega)$ の鋭いプラズモンの極への敏感さを取り除く。

第 1 の要点は、節（[§](#offset-gamma)）で述べる
改良した offset-$\Gamma$ 法である。これは
式 (I20) の $W^{\rm c}({\bf k}\to 0,\omega)$ の発散を扱う。
そのために、
${W^{\rm c} }({\bf k}=0,\omega)$ の代わりに、発散しない
有効相互作用 $\overline{W^{\rm c} }({\bf k}=0,\omega)$ を定める。
すると、式 (I15) と (I20) のどちらも
単純な離散和でとれる。

第 2 の要点は節（[§](#sigma-interpolation)）で述べる。
限られた数の ${\bf q}$ 点でだけ計算した $V^{\rm xc}_{\bf q}$ から、
BZ 全体の任意の ${\bf q}$ での $V^{\rm xc}_{\bf q}$ を
補間で与える方法である。
これは節（[§](#offset-gamma)）の offset-$\Gamma$ 法で要る。
つまり、${\bf q}=0$ の近くのいくつかの ${\bf q}$ 点で固有関数を計算しなければならない。
補間のために、[Obata2025] Eq. (3) の静的な
非局所ポテンシャル $V^{\rm xc}$ を、実空間でよく局在した
MTO で展開する。このように、そうした MTO は 2 つの目的に使う。
1 つは固有関数の基底として、もう 1 つは
$V^{\rm xc}$ を展開する基底としてである。$V^{\rm xc}_{\bf k}({\bf r},{\bf r}')$ の補間の手順は、
Ref.[kotani_quasiparticle_2007] の
複雑な補間の手順より
安定で単純になる。
これは、いまはよく局在した MTO を使うからである。`sigm.<sname>` は $\Sigma-V_{\rm xc}$ を、既約な $\mathbf q$ での MTO 基底の行列としてもつ。
`gwsc --mlo`（開発中）は、代わりに MLO の表現でもつ（[MLO-gwsc](./mlo_gwsc)）。
Hamann と Vanderbilt の平面波に基づく QSGW 法 [hamann09] では、
$V^{\rm xc}$ を MTO の代わりに最大局在 Wannier 関数で展開する。

実際の実装では、効率よく数値計算するために、LDA または GGA の交換相関
ポテンシャル $V^{\rm xc}_{\rm LDA}$ を使う。
すなわち、コアの固有関数と MT の中の動径関数を作るのに使う
（ここでは GGA を使うときも添字 LDA を使う。「LDA/GGA」は LDA または GGA を意味する）。
差 $V^{\rm xc}-V^{\rm xc}_{\rm LDA}$ を
BZ での補間（節（[§](#sigma-interpolation)）で説明する）に使う。
$V^{\rm xc}_{\rm LDA}$ が $V^{\rm xc}$ をおおまかに近似している限り、
この差は数値的に小さいからである。
効率よく数値計算するために $V^{\rm xc}_{\rm LDA}$ を使うこの手順は、
実際には最終的な数値結果にごく弱い依存しか与えない。
ただし、形式的には結果は LDA/GGA の交換相関汎関数に依存しない。

## ${\bf q}$ と ${\bf G}$ ベクトルの生成。qg4gw

理論は [kotani_quasiparticle_2014] の 3.2 節を見よ。
${\bf k}$ 点は BZ での積分に使う正則メッシュの点である。正則メッシュの点は次で与える。
$$\begin{aligned}
{\bf k}_{i_1,i2,i3}=(i1/N1)*{\bf P}_1+(i2/N2)*{\bf P}_2+(i3/N3)*{\bf P}_3
\end{aligned}
\tag{I21}
$$
ここで ${\bf P}_i$ は基本ベクトル $\frac{2\pi}{\rm alat} \times$ `qlat(:,i)`、すなわち逆格子ベクトル ${\bf Q}_i$ である（`[gw] chi_RegQbz = true` が既定で、${\bf k}=0$ を含む。`false` ではメッシュを半ステップずらす）。
（BZ で周期的な）量 $f({\bf k})$ は、正則メッシュの点での和として積分できる。
ただし $f({\bf k}=0)$ が発散するときは、注意深い扱いをする。

`mkqg.f90` の `mkQG2` が **qg4gw**（`main_qg4gw.f90`）の主要部である。ジョブ `--job=N`（`iq0pin`）は、1（offset-Γ 点 Q0P を生成する。`gwsc`）、2（Q0P を読む。`job_eps` は `--job=2 --n1n2n3eps` を使う）、3（epsPP）、4（χ+-。ngp = ngc = 0）、201（Q0P に加えて x、y、z 方向の Q0P）、10002（Γ セルを飛ばす）、20002（Γ セルだけ）、40001（マグノン用の `QforEPSL`）である。`--n1n2n3eps` を付けると、`[gw] n1n2n3eps` があればそのメッシュを使う。
目的は、正則メッシュの点と offset-Γ 点で必要な q と G ベクトルを生成することである。
BZ のデータ `__BZDATA` と、各 q について $|{\bf q}+{\bf G}|$ が `[gw] QpGcut_psi` より小さい G ベクトル（固有関数用、`__QGpsi`）と `QpGcut_cou` より小さい G ベクトル（MPB の IPW 用、`__QGcou`）を書く。

### $G$ ベクトルを作る: getgv2

${\bf G}$ ベクトルを得るには `getgv2.f90` の `getgv2` を使う。その先頭は次のとおりである。

```fortran
subroutine getgv2(alat,plat,qlat,q, QpGcut,job, ng,ngvec) ! Set up a list of recip vectors within cutoff |q+G| < QpGcut a.u.
  ! job==1 -> return ng (number of G ) only
  ! job==2 -> return ng and ngvec
  ! True G is given by  G(1:3,1:ng) = 2*pi/alat * matmul(qlat * ngvec(1:3,1:ng))
```

これで、与えた $q$ に対する $\{ {\bf G}\}$ が得られる。
`getgv2` のアルゴリズムは少し込み入っている。
まず下限と上限 `n1min` $\le n_1 \le$ `n1max` を（1 の余裕を付けて）与える。
ここで ${\bf G}= n_1 Q_1 + n_2 Q_2 + n_3 Q_3$ であり、$n_2$ と $n_3$ も同様である。

${\bf x}= {\bf q}+{\bf G}$ の三次元空間を考える。
与えた ${\bf q}$ に対し、許される ${\bf G}$ は正則メッシュの点の集合 $\{ {\bf q}+{\bf G}\}$ をなす。
`getgv2` の目的は、これらのメッシュ点のうち
$|{\bf q}+{\bf G}|\lt $`QpGcut` を満たすものだけを拾うことである。

まず、$|{\bf q}+{\bf G}|$ の最大値（=`QpGcut`）に対して $n_1$ の許される範囲を計算できる。
$|{\bf x}|=|{\bf q}+{\bf G}|$=`QpGcut`
は球を与えることに注意する。この球の中のメッシュ点を拾えばよい。
$n_1$ を指定すると平面が一つ決まる
（$n_2,n_3$ は任意の値をとってよい）。
範囲は、球
$|{\bf q}+{\bf G}|$=`QpGcut` が $n_1$ で指定した平面と交わるという条件で決まる
（正確には、そのような $n_1$ は実数である）。
平面の法線ベクトルは外積 $Q_2 \times Q_3$ である。

$n_1$ の範囲（$n_2,n_3$ も同様）が得られたら、
$(n_1,n_2,n_3)$ に対する ${\bf q}+{\bf G}$ が許されるかどうかを単純に調べる。

### $G$ ベクトルを作る: shortn3

$\{Q_i\}$ を法として最短のベクトルを見つける。すなわち ${\bf q}$ を第 1 BZ に引き戻す（`m_shortn3.f90` の `shortn3`）。
注意: ${\bf q}$ が BZ の境界にあるときは一意でないことがある。
そのときは、すべての ${\bf q}$ と縮重度を知る必要がある。

## 混合積基底 {#mixed-product-basis}

混合積基底は二種類の基底、すなわち
積基底と IPW からなる。
$\{M^{\bf k}_I({\bf r}) \}\equiv 
\{ P^{\bf k}_{\bf G}({\bf r}), B_{ {\bf R}\mu}^{\bf k}({\bf r})\}$
であり、添字 $I\equiv \{ {\bf G},{\bf R}\mu\}$
が基底の要素を区別する。
PB $B^{\bf k}_{ {\bf R}\mu}({\bf r})$ は次で定義する。
$$\begin{aligned}
  B^{\bf k}_{ {\bf R}\mu}({\bf r}) & =  
          \sum_{\bf T} B_{ {\bf R}\mu}({\bf r-R}-{\bf T}) e^{i {\bf k\cdot T} },
\end{aligned}
\tag{I22}
$$
ここで $B_{ {\bf R}\mu}({\bf r})$ は動径関数の積から作る。
$B_{ {\bf R}\mu}({\bf r})$ は実数で、
MT の外 $|{\bf r}| \gt  R$ ではゼロである（Sec.([§](#product-basis-hbasfp0)) を見よ）。
$\{B_{ {\bf R}\mu}({\bf r})\}$ は正規直交になるように作る。
$$\begin{aligned}
  \int_{|{\bf r}| < R} B_{ {\bf R}\mu}({\bf r}) B_{ {\bf R}\mu'}({\bf r}) d^3r =
  \delta_{\mu \mu'}.
  \end{aligned}
\tag{I23}
$$
加えて、$\{B_{ {\bf R}\mu}({\bf r})\}$ と
$\{ P^{\bf k}_{\bf G}({\bf r})\}$ が直交することは自明である。
したがって重なり行列
$$\begin{aligned}
 O^{\bf k}_{IJ}= 
  \int_\Omega \{M_I^{\bf k}({\bf r})\}^* {M}_J^{\bf k}({\bf r}) d^3r
\end{aligned}
\tag{I24}
$$
は PB のブロックでは単位行列であり、IPW のブロックでだけ対角でない。IPW $\{ P^{\bf k}_{\bf G}({\bf r})\}$ が直交しないからである。
そこで、式 (I13) と同じやり方で ${M}^{\bf k}_{I}({\bf r})$ の双対
$\tilde{M}^{\bf k}_{I}({\bf r})$ を定義すると便利である。

固有関数の積から作った関数は、
${M}^{\bf k}_{I}({\bf r})$ の基底で事実上完全に
次のように展開できる。
$$\begin{aligned}
\begin{cases}
  \displaystyle 
  F^{\bf k}({\bf r}) = \sum_I M_I^{\bf k} ({\bf r}) F_I({\bf k}) & \cr
  \displaystyle
  F_I({\bf k}) = 
  \int_\Omega \{\tilde{M}_I^{\bf k}({\bf r})\}^* 
  F^{\bf k}({\bf r}) d^3r.
\end{cases}
  \end{aligned}
\tag{I25}
$$

GW の部分は $M^{\bf k}_I$ と、式 (I14) の Coulomb 正規直交化した $E^{\bf k}_\mu$ の両方を使う。どの段でどちらを使うかは式 (I14a) を見よ。

### 積基底（hbasfp0） {#product-basis-hbasfp0}

原子 $a$ の動径関数を次のように書く。
$$\begin{aligned}
   u_{apl\sigma}(r) = r \phi_{apl\sigma}(r),
\end{aligned}
\tag{I26}
$$
ここで添字 $p$ は、$\phi$ なら 1、$\dot{\phi}$ なら 2、
局所軌道なら 3 などをとる。動径関数には $m$ 依存性
（$m$ は $L=(l,m)$ の $m$）を許さない。

加えて、$p$ はコア関数の添字もとる。コアと価電子の動径関数は、`main_hbasfp0.f90` で `__PHIVC`
（`lmf --jobgw=1`、`sugw.f90` が書く）からまとめて読む。
各原子と各スピンについて、まずコア（`ncore`。その $l$ と $n$ は `lcindx`、`ncindx`）、次に価電子の動径関数（`nrad`、`lindx_r`、`nindx_r`）が並ぶ。
したがって価電子の $n$ は、その $l$ のコアの数（`nc_max`）だけずれる。
正規直交化した動径関数 $u_{apl\sigma}(r)$ は `phitoto` に、生のものは `phitotr` に入る（コアも同じ）。

*本当の* 動径関数は
$\phi_{apl\sigma}(r) = u_{apl\sigma}(r)/r$ であることに注意する。
規格化は
$1=\int_0^{R_a} \{ u_{apl\sigma}(r) \}^2 dr = 
   \int_0^{R_a} \{ \phi_{apl\sigma}(r) \}^2 r^2dr$ である。
関数 $u_{apl\sigma}(r)$ は `phitoto` に入る。

積関数を作るときは、
次で与えるスピン平均した関数 `phiav` を使う。
$$\begin{aligned}
   u_{apl}(r) = \frac{1}{N_{\rm spin} }\sum_{\sigma}
   u_{apl\sigma}(r).
\end{aligned}
\tag{I27}
$$
（`basnfp.f90` の `subroutine basnfp_v2` を見よ。生の `phitotr` を、下向きスピンの関数の符号を上向きスピンのものに揃えて平均する。）
これらから積関数 `rprod` を作る。
$$\begin{aligned}
  \tilde{b}_{al\nu}(r) = \frac{1}{r}{u_{a p l }(r) u_{a p' l'}(r)}
  = r \phi_{a p l}(r) \phi_{a p' l'}(r),
\end{aligned}
\tag{I28}
$$
ここで添字 $l$ は $|l-l'| \le l \le |l+l'|$ を動き、
$\nu$ は組み合わせ $(p,p')$ の添字である。
*本当の* 積関数は次で与えられることに注意する。
$$\begin{aligned}
  \tilde{B}_{al\nu}(r) = \frac{1}{r}\tilde{b}_{al\nu}(r).
\end{aligned}
\tag{I29}
$$
この関係は $\phi_{apl}(r)=u_{apl}(r)/r$ と同じである。

次に重なり行列 `ovmt` を計算する。
$$\begin{aligned}
  O_{\nu_1\nu_2} &= \int_0^{R_a} 
  \tilde{B}_{al\nu_1}(r) \tilde{B}_{al\nu_2}(r) r^2 dr \\
  &= \int_0^{R_a}  \phi_{a p_1 l_1}(r) \phi_{a p'_1 l'_1}(r)
                  \phi_{a p_2 l_2}(r) \phi_{a p'_2 l'_2}(r) r^2 dr
\end{aligned}
\tag{I30}
$$
そして重なり行列の固有値問題
$Oz_{\nu}=\epsilon_{\nu} z_{\nu}$ を `call rss(..)` で解く。
（`basnfp_v2` を見よ。）

固有値が $\epsilon_{\nu} \lt  {\rm tolerance} \sim 10^{-4}$ の固有ベクトル
$z_{\nu}$ を捨てると（tolerance は `[product_basis] pb_tolerance`（[GW input](./gwinput#product-basis)）で $l$ ごとに与える。それより大きい $l$ には最後の値を使う）、
最終的に、積関数の一次結合として
最適な積関数が次のように得られる。
$$\begin{aligned}
   b_{al\nu}(r) = \frac{1}{\sqrt{\epsilon_{\nu} }}
   \sum_{\nu'}\tilde{b}_{al\nu'}(r)z_{\nu'\nu},
\end{aligned}
\tag{I31}
$$
これは `rprodx` に入り、`__BASFP<ic>` に書かれて、
続く Coulomb 行列のルーチン `hvccfp0` で使われる。
もちろん、*本当の* 積関数は
$B_{al\nu}(r) = b_{al\nu}(r)/r$ である。

最適な積関数の規格化は
標準出力で確かめる
（`gwsc` が書く `lbasC` と `lbas` を見よ）。

    ...
    Use rs diagonalization for real symmetric 
    Diag ibx ovv=  1 0.9999999999999930D+00 eb=   0.2716113799D-01 nod=   2
    Diag ibx ovv=  2 0.9999999999999980D+00 eb=   0.4993303381D-01 nod=   3
    Diag ibx ovv=  3 0.1000000000000001D+01 eb=   0.1467546915D+00 nod=   3
    Diag ibx ovv=  4 0.9999999999999996D+00 eb=   0.4415639258D+01 nod=   0

    ...

`basnfp` では、必要な動径積分
$\lt \phi \phi B\gt $`=ppbrd` をすべて計算する。
$$\begin{aligned}
   \langle\phi \phi B\rangle &= \int_0^{R_a}
   \phi_{ap_1 l_1}(r) \phi_{ap_2 l_2}(r) B_{al\nu}(r) r^2 dr \\
   &= \int_0^{R_a} \frac{1}{r}
     u_{ap_1l_1}(r) u_{ap_2l_2}(r) b_{al\nu}(r) dr,
\end{aligned}
\tag{I32}
$$
これは `__PPBRD_V2_<ic>` に入る。
`rdpp.f90` の `call rdpp(...`（module `m_rdpp`。`m_zmel` が使う）で、
`ppbrd` を確保して読む。

加えて、\"回転した Clebsch-Gordan 係数\"
$C(L,L_1,L_2,g)$ `cgr(lm,lm1,lm2,ng)` を読む。
ここで $g$ は空間群の
係数（点群の対称操作で回転したもの）の添字である。

## Coulomb 行列（hvccfp0）

式 (I13) のすぐ後に現れる
$v_{IJ}^{\bf k}= \langle M^{\bf k}_{I} |v|  M^{\bf k}_J \rangle$ を計算しなければならない。
`hvccfp0`（`main_hvccfp0.f90`）は
$\frac{\exp(- |\kappa| |{\bf r}_1-{\bf r}_2|)}{|{\bf r}_1-{\bf r}_2|}$ の場合を扱える。
エネルギー変数 $E$ は $E=\kappa^2$ で与える。
負の $E$ では $\kappa$ の虚部を正にとる。
この $E=$`eee=screenfac()` は `main_hvccfp0.f90` にある。`switch.f90` の `screenfac()` は $-$`TFscreen`$^2$ を返す。
既定は `[gw] TFscreen` $=\sqrt{10^{-5}}$、すなわち $E=-10^{-5}$ で、ほぼ裸の Coulomb 相互作用である（小さすぎる値は $q=0$ で Coulomb 行列に負の固有値を生む）。
例えば必要なら、`ctrlg.<sname>.toml` の `[gw]` に `TFscreen = 0.1` と書けば $|\kappa|=0.1$（a.u.）にできる。

MPB は IPWcou（$|{\bf q}+{\bf G}|$`<QpGcut_cou`）と PB からなる。
すなわち $\{M^{\bf k}_J\}=\{P_{\bf G}^{\bf k}({\bf r}), B^{\bf k}_{ {\bf R}{liL}\mu}({\bf r})\}$ である。
与えた ${\bf k}$（簡単のため明示しない場合がある）に対し、
IPWcou は ${\bf G}$ ベクトルだけで指定される。
PB は各 $l$ の動径関数だけで指定できる
（$m$ 依存性もスピン依存性もないからである）。
PB は Sec.([§](#mixed-product-basis)) のやり方で作る。
$GW$ 計算には Coulomb 行列要素
$\langle B|v|B \rangle,\langle B|v|P_{\bf G}\rangle,$ と $\langle P|v|P \rangle$ が要る。
これらの計算は、`main_hvccfp0.f90` から呼ぶ `vcoulq_4`（`mkjp.f90`）が受け持つ。

### 球 Bessel 関数と関連する関数

$X_L({\bf r})=X_l(r) Y_L(\hat{r})$ という記法を使う。
動径部分は $l$ だけに依存することに注意する。

与えたエネルギー $E$ に対し、Bessel 関数
$J_L(E,{\bf r}),H_L(E,{\bf r})$ は次のとおりである
（ここで $\sqrt{E}=\kappa$ を導入する。
$E\lt 0$ では $\sqrt{E}=\kappa=i|\kappa|$ である）。
$$\begin{aligned}
 &{J}_L(E,{\bf r})=
j_l(\sqrt{E} r)Y_L(\hat{r}) =
\begin{cases}
j_l(i|\kappa| r)Y_L(\hat{r}) \ {\rm for} \ E<0, \\
j_l(\kappa r)Y_L(\hat{r})    \ {\rm for \ E>0},
\end{cases}
 \\
 &{H}_L(E,{\bf r})=h_l(i|\kappa| r)Y_L(\hat{r})
\ {\rm for \ E<0}, 
\end{aligned}
\tag{I33}
$$
ここで $j_l(z)$ と $h_l(z)$ は通常の球 Bessel 関数と球 Hankel 関数で、
$z \sim 0$ で
$$\begin{aligned}
 &j_l(z) \sim \frac{z^l}{(2l+1)!!}  \\
 &h_l(z) \sim \frac{-i(2l-1)!!}{z^{l+1} }
\end{aligned}
\tag{I34}
$$
のようにふるまう。
便利のため、
**Methfessel の Bessel 関数（の規約）**
$\bar{J}_L=\bar{J}_lY_L$ と $\bar{H}_L=\bar{H}_lY_L$ を定義する。ここで
$$\begin{aligned}
 &\bar{J}_l(E,{\bf r})=j_l(i|\kappa| r) / (i |\kappa|)^l  \\
 &\bar{H}_l(E,{\bf r})=h_l(i|\kappa| r) i (i |\kappa|)^{l+1}. 
%&\bar{J}_l(E,\bfr)=(2l+1)!! j_l(i|\kappa| r) / (i |\kappa|)^l  \\
%&\bar{H}_l(E,\bfr)=\frac{1}{(2l-1)!!} h_l(i|\kappa| r) i (i |\kappa|)^{l+1}. 
\end{aligned}
\tag{I35}
$$
（メモ: 例えば $7!!=7\cdot5\cdot3\cdot1$ `;function fac2m(i)`）。
これらは $E\le 0$ で実関数である。
`mkjp.f90` の `genjh` の注を見よ。
$E=0$ では次に帰着する。
$$\begin{aligned}
 &\bar{J}_l(E=0,{\bf r})=r^l/(2l+1)!!,  \\
 &\bar{H}_l(E=0,{\bf r})=r^{-l-1} (2l-1)!!.
\end{aligned}
\tag{I36}
$$

（Coulomb 行列を作る `hvccfp0` では、数値的な問題を避けるため、
既定で $E=0$ の代わりにごく小さな負の $E=-10^{-5}$ を使う
（`switch.f90` の既定の `screenfac` を見よ）。）

bessl を定義するソースコードは、いくつかの規約が混ざっているため
少し紛らわしい\...（`besslr.f90` の冒頭の注）。

主なものは `mkjp.f90`（`genjh`）の `call bessl(ex2,lx,phi(0:lx),psi(0:lx))` である。
これは `besslr.f90` の `bessl(ex2,lmax,phi,psi)` で定義される（module `m_bessl` にも同じ `bessl` があり、`vcoulq_4` で使う）。
`lx` は $l$ の上限である。
与えた `ex2`$=E\times x^2$ に対し、
これは球 Bessel 関数 $j_l(\kappa x)$ を
次のように返す。


$$\begin{aligned}
&E<0,\ \kappa=i|\kappa|,\ |\kappa|=\sqrt{|E|}: \\
&\qquad \bar{J}_l(E,{\bf r})=j_l(i|\kappa|r)/(i|\kappa|)^l={\tt phi(l)}\,r^l, \\
&\qquad \bar{H}_l(E,{\bf r})=h_l(i|\kappa|r)\,i (i|\kappa|)^{l+1}={\tt psi(l)}/r^{l+1}, \\
&E>0,\ \kappa=\sqrt{E}: \\
&\qquad \bar{J}_l(E,{\bf r})=j_l(\kappa r)/\kappa^l = {\tt phi(l)}\,r^l, \\
&\qquad \bar{N}_l(E,{\bf r}) = {\tt psi(l)}/r^{l+1}\quad (\text{Neumann function})
\end{aligned}
\tag{I37}
$$

ここで $0\le l \le{\tt lmax}$ である。この定義から、$E \to 0$ で `phi`$=1/(2l+1)!!$、
`psi`$=(2l-1)!!$ となる。
ここで $n_l(r)$ は球 Neumann 関数である。
すなわち `psi` は、負の $E$ では Hankel 関数、
正の $E$ では Neumann 関数に対するものである。
例えば `m_hansr.f90` の `hansmr` を見よ。`psi` はあまり使わない。

便利のため、次をときどき使う。
$$\begin{aligned}
 &R^{+l}=\bar{J}_l(E,r)(2l+1)!!={\tt phi(l)}r^l(2l+1)!!,   \\
 &R^{-l}=\bar{H}_l(E,r)/(2l-1)!!={\tt psi(l)}r^{-l-1}(2l-1)!!.
\end{aligned}
\tag{I38}
$$

ここで $E\to 0$ で $R^{+k} \rightarrow r^k$、
$E\to 0$ で $R^{-k} \rightarrow r^{-k-1}$ となることに注意する。
これらの $rR^{+k},rR^{-k}$ は、サブルーチン `genjh` で作る `rkpr,rkmr` である。

### Green 関数

（この小節は読み飛ばしてよい。）
ここでは自由空間の Green 関数
$G({\bf r}-{\bf r}',E)$ を説明する。
次を満たす $G({\bf r}-{\bf r}',E)$ から始める。
$$\begin{aligned}
(E+i\delta+\nabla^2)G({\bf r}-{\bf r}',E)=\delta({\bf r}-{\bf r}').
\end{aligned}
\tag{I41}
$$
（$+i\delta$ は時間軸に沿った境界条件を指定するためのものである。これにより遅延 Green 関数が選ばれる。）
大まかに言えば、これは $(\omega -H) G=1$ である。
その Fourier 変換は簡単に
$G({\bf k},E)=1/(E-|{\bf k}|^2+ i \delta)$ と書ける。
これを Fourier 逆変換すると $G({\bf r}-{\bf r}',E)$ が得られる。
それは
$$\begin{aligned}
G({\bf r}-{\bf r}',E)=-\frac{1}{4 \pi} 
\frac{e^{i \kappa|{\bf r}-{\bf r}'|} }{|{\bf r}-{\bf r}'|}, \end{aligned}
\tag{I42}
$$
である。ここで $\kappa \equiv \sqrt{E}$ であり、$E\lt 0$ では $\kappa$ の虚部は正である。
すなわち負の $E$ では $\kappa = i |\kappa|= i \sqrt{-E}$ である。
式 (I42) は Helmholtz
微分方程式の解にほかならない。$E=0$ では通常の Poisson 方程式に帰着する。
$E\lt 0$ では Thomas-Fermi 型の関数になる。
式 (I42) の分子は $\exp(i \kappa |{\bf r}-{\bf r}'|)=\exp(- |\kappa| |{\bf r}-{\bf r}'|)$ である。

### 使う公式 {#used-formulas}

$$\begin{aligned}
\frac{1}{|{\bf r}-{\bf r}'|}
=4 \pi \sum_K \frac{r^k_<}{r^{k+1}_>}\frac{1}{2k+1}
Y^*_K(\hat{ {\bf r} })Y_K(\hat{ {\bf r} }') \end{aligned}
\tag{I43}
$$
文献 [svane_evaluation_1986] の Appendix A を見よ。
これは次のように一般化される。
$$\begin{aligned}
\frac{e^{-|\kappa| |{\bf r}-{\bf r}'|} }{|{\bf r}-{\bf r}'|}
&=4 \pi \sum_K R^{+k}(r_<) R^{-k}(r_>)\frac{1}{2k+1}
Y^*_K(\hat{ {\bf r} })Y_K(\hat{ {\bf r} }') \\
&=4\pi \sum_K 
\bar{J}_k(E,r_<) \bar{H}_k(E,r_>)\, Y^*_K(\hat{ {\bf r} })Y_K(\hat{ {\bf r} }'),
\end{aligned}
\tag{I44}
$$
$\bar{J},R^{+l}$ の定義は式 (I35) と式 (I38) にある。
ここで $E\lt 0$、$\kappa= i|\kappa|=i\sqrt{-E}$ である。

$$\begin{aligned}
 &\exp( i {\bf k} {\bf r}) = 4 \pi \sum_L i^l 
j_l(|{\bf k}|r)  Y^*_L( {\bf \widehat{ k} }) Y_L(\widehat{ {\bf r} }) 
\\
 &\frac{2l+1}{4 \pi} P_l(\cos\Theta) = \sum_m Y^*_L( {\bf \widehat{r} }_1) 
Y_L({\bf \widehat{r} }_2) \ \ \ \ \ \ 
[\cos \Theta = {\bf \widehat{r} }_1 \cdot {\bf \widehat{r} }_2].
\end{aligned}
\tag{I45}
$$

IPW の重なり（`pwmat.f90` の `mkppovl2` と `matgg2`）は ${\bf G}'-{\bf G}$ だけに依存し、${\bf k}$ には依存しない。
$$\begin{aligned}
\langle P^{\bf k}_{\bf G}|P^{\bf k}_{ {\bf G}'} \rangle
&= \Omega \delta_{ {\bf G},{\bf G}'} \\ &\quad - \sum_{a} e^{ i ({\bf G}'-{\bf G})\cdot {\bf R}_a}\,\frac{4\pi}{g^3}\big(\sin(gR_a)-gR_a\cos(gR_a)\big),
\end{aligned}
\tag{I46}
$$
ここで $g=|{\bf G}'-{\bf G}|$ であり、$g=0$ では $a$ の項は $4\pi R_a^3/3$ である。

### Hankel 関数と構造定数

$v({\bf r},{\bf r}')=\frac{e^{-|\kappa||{\bf r}-{\bf r}'|} }{|{\bf r}-{\bf r}'|}$
は一中心展開の式 (I44) で展開できる。

式 (I44) の Hankel 関数には、
Hankel 関数の中心外展開定理の式 (I47) を使う。すなわち、
中心が ${\bf X}\equiv{\bf R}+{\bf T}$ にある Hankel 関数 $H_L({\bf r}-{\bf X})$ は、
中心が ${\bf X}'$ にある Bessel 関数で展開できる。
$$\begin{aligned}
\bar{H}_L(E,{\bf r}-{\bf X})=\frac{1}{4\pi} \sum_{L'}\bar{J}_{L'}(E,{\bf r}-{\bf X}') S_{ {\bf X}'{L'},{\bf X}{L} }\ ,
\end{aligned}
\tag{I47}
$$
ここで $E\lt 0$ である（負のエネルギーの Hankel 関数）。
$S_{ {\bf X}'{L},{\bf X}{L'} }$ の $E$ 依存性は明示しない。
$\bar{J}_L$ と $J_L$ の違いに注意する（$\bar{H}_L$ も同様）。

したがって $({\bf R}',{\bf T}')\neq({\bf R},{\bf T})$ では、
二中心展開が得られる。
$$
\begin{aligned}
&\frac{e^{-|\kappa||{\bf r}+{\bf R}+{\bf T}-({\bf r}'+{\bf R}'+{\bf T}')|} }
{|{\bf r}+{\bf R}+{\bf T}-({\bf r}'+{\bf R}'+{\bf T}')|} \\
&\qquad = \sum_L \sum_{L'} \bar{J}_L(E,{\bf r})\,
S_{ {\bf R}+{\bf T}{L},{\bf R}'+{\bf T}'{L'} }\, \bar{J}_{L'}(E,{\bf r}')
\end{aligned}
\tag{I48}
$$
これは $({\bf R}',{\bf T}')\neq({\bf R},{\bf T})$ で成り立つ。

$S_{ {L},{\bf X}{L'} }$ の Bloch 和は、
${\bf k}$ の構造定数を次のように与える。
$$\begin{aligned}
S^{\bf k}_{ {\bf R}{L},{\bf R}'{L'} }=
\sum_{ {\bf T}'} S_{ {\bf R}{L},{\bf R}'+{\bf T}'{L'} }\exp(i {\bf k}{\bf T}').
\end{aligned}
\tag{I49}
$$
既定では $E=-10^{-5}$（`screenfac`）を使う。ほぼ裸の Coulomb 相互作用である。

`m_hvccfp0_util.f90` で定義する `strxq` の先頭を見よ（`main_hvccfp0.f90` から呼ぶ。$E\lt 0$ では `strxq_all` ですべての組を一度に計算する）。
このルーチンは、通常の Bloch 和をとった Hankel 関数の一中心展開のためのものである（smooth Hankel 用ではない）。
この結果は最終的に、
式 (I50) で使う構造定数の Bloch 和
$S^{ {\bf k} }_{ {\bf R}{L},{\bf R}'{L'} }$ に変換される。

`strxq` は因子 $4\pi$ を除いた構造定数を返す（その展開定理は $H_{RL}(E,r)=\sum J_{R'L'}(E,r)\,S_{R'L',RL}$ であり、その記法では $S_{R'L',RL}=4\pi\sum_{L''}C_{LL'L''}(-1)^l(-E)^{(l+l'-l'')/2}H_{L''}(E,R-R')$ である。ここで $R$ は $T$ を含む）。
`main_hvccfp0.f90` はこれに $4\pi$ を掛け（`strx = fpi*s`）、`strx` が式 (I48)–(I50) の $S^{\bf k}_{{\bf R}L,{\bf R}'L'}$ になる。

### $\langle B|v|B \rangle$ の部分 {#b-v-b-part}

$\langle B|v|B \rangle$ の部分から始める。
この計算には、構造定数といくつかの種類の
動径積分が要る。
Bloch 和をとった構造定数 `strx`
$=S^{\bf k}_{RL,R'L'}$ を使うと、
$$
\begin{aligned}
v^{\bf k}({\bf r},{\bf r}')
&=\sum_{ {\bf T}'} v({\bf r},{\bf r}'+{\bf T}') e^{i {\bf k}{\bf T}'} \\
&=e^2 \sum_{L,L'} \bar{J}_L(E,{\bf r}-{\bf R})\,
S^{ {\bf k} }_{ {\bf R}{L},{\bf R}'{L'} }\, \bar{J}_{L'}(E,{\bf r}'-{\bf R}'),
\end{aligned}
\tag{I50}
$$
となる。これは ${\bf R}\ne {\bf R}'$ のときである。
`hvccfp0` では $e^2=1$（a.u.）を使うことに注意する。
（`mkjp.f90` の冒頭の注を見よ。また `hvccfp0` の最後の規格化の確認にも注意する。
$\langle \exp(i {\bf q}{\bf r}) |v |\exp(i {\bf q}{\bf r})\rangle=4\pi \Omega/|{\bf q}|^2$ であり、$\Omega$ はセルの体積である。）

この `strx` は `call vcoulq_4` で使う。
`nlx1` は $R$
`(=ibas1)` の `(l+1)**2` を意味する。`nlx2` も同様である。

$({\bf R},{\bf T})=({\bf R}',{\bf T}')$ の寄与を除けば、
$\langle B^{\bf k}_{ {\bf R}{L}\mu}({\bf r})|v^{\bf k}({\bf r},{\bf r}')|B^{\bf k}_{ {\bf R}'{L'}\mu'}({\bf r}') \rangle$ は\
`rojb` 積分 $\rho^l(B_{ {\bf R}{l}\mu})$ から次のように評価できる。
$$
\begin{aligned}
\rho^l(B_{ {\bf R}{l}\mu})
&=\int_0^{R} r\bar{J}_l(E,r) \ r B_{ {\bf R}{l}\mu}(r)  dr \\
&= \frac{1}{(2l+1)!!}\int_0^{R} {\tt rkpr(r) \ rprodx(r)} \ dr
\end{aligned}
\tag{I51}
$$
ここで $B_{ {\bf R}{l}\mu}(r)$ は $B_{ {\bf R}{L}\mu}({\bf r})$ の動径部分である。
この `rojb` 積分は、
`mkjp.f90` のサブルーチン `mkjb_4` で計算する。
動径関数として
`rprodx` $=r B_{ {\bf R}{l}\mu}(r)$、
`rkpr`$=r\bar{J}_l(E,r)(2l+1)!!$、
`rkmr`$=r\bar{H}_l(E,r)/(2l-1)!!$ を使う。
（ここで `rkpr` と `rkmr` は、$E =0$ では
$r^l$ と $r^{-l-1}$ に比例する。）

$({\bf R},{\bf T})=({\bf R}',{\bf T}')$ からの寄与を
加えなければならない。これは ${\bf k}$ に依存せず、\
`sgbb` 積分で与えられる。これも `mkjb_4` で計算する。
$$
\begin{aligned}
\sigma^l(B_{ {\bf R}{l}\mu},B_{ {\bf R}{l}\nu})
&={4\pi}\int_0^{R}\!\!\int_0^{R} r_<\bar{J}_l(E,r_<)\, r_>\bar{H}_l(E,r_>) \\
&\qquad\times r B_{ {\bf R}{l}\mu}(r)\, r' B_{ {\bf R}{l}\nu}(r')\, dr dr'  \\
&=\frac{4\pi}{2l+1}\int_0^{R}\!\!\int_0^{R} {\tt rkpr(r_<) \ rkmr(r_>)} \\
&\qquad\times {\tt rprodx(n1,r)}\ {\tt rprodx(n2, r')}\, dr dr'.
\end{aligned}
\tag{I52}
$$

積分 `rojb` と `sgbb` を使って、`vcoulq_4` で
$\langle B|v|B \rangle$ を次のように計算する
（`nbloch` は PB の総数を意味する。`vcoulq_4` で PB の組についてループする）。

```fortran
vcoul(ibl1,ibl2) = rojb(n1, l1, ibas1) *strx(lm1,ibas1,lm2,ibas2) *rojb(n2, l2, ibas2)    ! offsite Coulomb
if(ibas1==ibas2.AND.lm1==lm2) vcoul(ibl1,ibl2) = vcoul(ibl1,ibl2) + sgbb(n1,n2,l1, ibas1) ! sigma-type onsite parts
```

### $|P^{\bf k}_{ {\bf G}'} \rangle$ の $RL$ 展開

$\langle P^{\bf k}_{\bf G}|v| P^{\bf k}_{ {\bf G}'} \rangle$ を評価するには、
十分大きな $l_{\rm Pmax}$ を使う限り、
${P}^{\bf k}_{ {\bf G} }$ の代わりに
$$\begin{aligned}
\bar{P}^{\bf k}_{ {\bf G} } \equiv 
\left( 1- \sum_{ {\bf R}{L} }^{l\le l_{\rm Pmax} } {P}_{ {\bf R}{L} } \right) 
e^{i ({\bf k+ G}){\bf r} }, \end{aligned}
\tag{I53}
$$
を使える。
ここで ${P}_{ {\bf R}{L} }$ は、
${ {\bf R}{L} }$ の寄与の成分を取り出す射影演算子である。
実際、十分大きな $l_{\rm Pmax}$ を使う。
$l_{\rm Pmax} = 2 \times l_{\rm max}$`=2*LMXA`（`basnfp.f90` と `hvccfp0` の `lxx = 2*lmxax`）であり、
ここで $l_{\rm max}$ は MT の中での固有関数の展開の
最大角運動量（式 (I1) の $\alpha^{ {{\bf k} }n}_{ {\bf R}u}$ の
$l$ の最大の打ち切り）である。
既定の設定では $l_{\rm max}=4$ を使うので、
$l_{\rm Pmax}=8$ を使う。

行列要素 $\langle P_1 |v| P_2 \rangle$ は
次のように計算できる。
$$
\begin{aligned}
\langle P_1(\phi\phi) |v| P_2(\phi\phi) \rangle
&=\sum_{\bf G_{1'} G_{1''} G_{2'} G_{2''} } 
\langle P_1(\phi\phi) | {P_{1'} } \rangle  
\langle {P}_{1'}| {P}_{1''} \rangle^{-1} \\
&\qquad\times\langle {P}_{1''} |v| {P_{2''} } \rangle  
\langle {P}_{2''}| {P}_{2'} \rangle^{-1} 
\langle {P_{2'} } | P_2(\phi\phi) \rangle,
\end{aligned}
\tag{I54}
$$
ここで $1\equiv ({\bf k},{\bf G_1})$ などである。
ここで $P_1(\phi\phi)$ は、
IPW の積から作った IPW を
表す。行列要素
$\langle {P}_{2''}| {P}_{2'} \rangle^{-1}$
（`__PPOvlp` に入る。逆行列は使う所でとる）
と $\langle {P_{2'} } | P_2(\phi\phi) \rangle$
（`__PPOvlpGG` に入る）
は `rdata4gw`（`sugw.f90` の中のブロック、すなわち `lmf --jobgw=1`）で与える。

### IPW の重なり行列（`__PPOvlp*` ファイル）

IPW の重なり行列は `lmf --jobgw=1`（`sugw.f90` のブロック `rdata4gw`）が書き、module `m_read_ppovl`（`rppovl.f90`）の `getppx2` が読む。`getppx2` は、行列要素 `<IPW psi|psi>` の IPW の部分のために `build_zmel`（`m_zmel.f90`）から呼ばれる。
二種類ある。

- MPB の IPW（`QpGcut_cou`）の $\langle P^{\bf k}_{\bf G}|P^{\bf k}_{{\bf G}'}\rangle$。IBZ の各 $q$ と offset-Γ 点（`qibze`）について、
  `__PPOvlp`（直接アクセス、`ppovl(ngc,ngc)`）、`__PPOvlpG`（G ベクトル `ngvecc`）、`__PPOvlp.info`（大きさ）。
  行列は ${\bf G}'-{\bf G}$ だけに依存する（式 (I46)）が、G の集合（`ngc`）が $q$ に依存するので、$q$ ごとに格納する。
- `__PPOvlpGG`: ${\bf G}'-{\bf G}$ の関数としての $\langle P_{\bf G}|P_{{\bf G}'}\rangle$（`ggg`、集合 `nvggg` の上）。$\langle G({\rm eigenfunction})\,G({\rm eigenfunction})|G({\rm Coulomb})\rangle$ に使う
  （二つの固有関数の IPW の積を MPB の IPW で展開する）。$k=0$ だけである。範囲は
  `QpGcutggg = (2+0.01)*QpGcut_psi + QpGcut_cou + 2π/alat*dQpG` と `QpGcutgcgp = (1+0.01)*QpGcut_psi + QpGcut_cou + 2*2π/alat*dQQ` で決める
  （三角不等式。`dQpG`、`dQQ` は offset-Γ 点のために範囲を広げる）。

### $\langle {P}^{\bf k}_{\bf G}|v| {P}^{\bf k}_{ {\bf G}'} \rangle$ {#p-v-p-part}

式 (I54) を評価するには、
その主要部 $\langle {P}^{\bf k}_{\bf G}|v| {P}^{\bf k}_{ {\bf G}'} \rangle$ が要る。それは次のように書ける。
$$
\begin{aligned}
\langle {P}^{\bf k}_{\bf G}|v| {P}^{\bf k}_{ {\bf G}'} \rangle
&\approx \langle \bar{P}^{\bf k}_{\bf G}|v| \bar{P}^{\bf k}_{ {\bf G}'} \rangle \\
&= \langle e^{i({\bf k}+{\bf G}){\bf r}} |v| e^{i({\bf k}+{\bf G}'){\bf r}} \rangle
-\sum_{ {\bf R}{L} } \langle P^{ {\bf k}+{\bf G} }_{ {\bf R}{L} }|v| e^{i({\bf k}+{\bf G}'){\bf r}} \rangle  \\
&\quad - \sum_{ {\bf R}'{L'} } \langle e^{i({\bf k}+{\bf G}){\bf r}} |v| P^{ {\bf k}+{\bf G}'}_{ {\bf R}'{L'} }\rangle \\
&\quad +\sum_{ {\bf R}{L} } \sum_{ {\bf R}'{L'} } 
\langle P^{ {\bf k}+{\bf G} }_{ {\bf R}{L} }|v|P^{ {\bf k}+{\bf G}'}_{ {\bf R}'{L'} } \rangle,
\end{aligned}
\tag{I55}
$$
ここで $P^{ {\bf k}+{\bf G} }_{ {\bf R}{L} }$ は
PW の ${ {\bf R}{L} }$ への射影を表す。
すなわち
$P^{ {\bf k}+{\bf G} }_{ {\bf R}{L} } \equiv 
{P}_{ {\bf R}{L} } e^{i ({\bf k+ G}){\bf r} }$ である。

##### 第 1 項

式 (I55) の右辺の第 1 項は
$$\begin{aligned}
\langle \exp(i({\bf k}+{\bf G}){\bf r}) |v| \exp(i({\bf k}+{\bf G}'){\bf r}) \rangle
=\frac{4 \pi \Omega}{|{\bf k}+{\bf G}|^2+|E|} \delta_{ {\bf G}{\bf G}'},
\end{aligned}
\tag{I56}
$$
である。IPW には単に $\exp(i({\bf k}+{\bf G}){\bf r})$
（規格化の前因子なし）を使うからである。
ここで $E$`=eee` は負（または（ほぼ）ゼロ）である。
これは `vcoulq_4`（`mkjp.f90`）の次の一行で書かれている。

    vcoul(nbloch + ig1, nbloch + ig1) = vcoul(nbloch + ig1, nbloch + ig1) + fpivol/(absqg2(ig1) - eee)

##### 第 2 項と第 3 項

第 2 項では、
$v$ は $|\exp(i({\bf k}+{\bf G}'){\bf r}) \rangle$ に対して対角なので、$v$ を $\frac{4 \pi \Omega}{|{\bf k}+{\bf G}'|^2+|E|}$ で置き換えられる。第 3 項も同様である。
$v$ を除くと、
$$
\begin{aligned}
\langle e^{i({\bf k}+{\bf G}){\bf r}}| P^{ {\bf k}+{\bf G}'}_{ {\bf R}'{L'} }\rangle
&=\sum_{ {\bf R}{L} } \langle P^{ {\bf k}+{\bf G} }_{ {\bf R}{L} }| P^{ {\bf k}+{\bf G}'}_{ {\bf R}'{L'} }\rangle  \\
&=\sum_{ {\bf R}{L} } \big({\tt pjyl\_}({\bf k}+{\bf G},L)\, e^{i({\bf k}+{\bf G}){\bf R}}\big)^* \\
&\qquad\times R^{JJ}(|{\bf k}+{\bf G}|,|{\bf k}+{\bf G}'|,l) \\
&\qquad\times {\tt pjyl\_}({\bf k}+{\bf G}',L)\, e^{i({\bf k}+{\bf G}'){\bf R}},
\end{aligned}
\tag{I57}
$$
となる。ここで次を使った。
$$\begin{aligned}
 &P^{ {\bf k}+{\bf G} }_{ {\bf R}{L} }({\bf r})=4 \pi i^l j_l(|{ {\bf k}+{\bf G} }|r) 
Y_L(\widehat{ {\bf k}+{\bf G} }) Y_L(\hat{ {\bf r} })\exp(i({\bf k}+{\bf G}){\bf R}) \\
 &= {\tt pjyl\_}({\bf k}+{\bf G},L) \bar{J}_l(|{ {\bf k}+{\bf G} }|r) Y_L(\hat{ {\bf r} })
\exp(i({\bf k}+{\bf G}){\bf R}), 
\end{aligned}
\tag{I58}
$$
ここで ${\bf r}$ は中心 ${\bf R}$ から測る。
ここで `pjyl_` は次で定義する。
$$\begin{aligned}
{\tt pjyl\_}({\bf k}+{\bf G},L) 
=4 \pi i^l |{\bf k}+{\bf G}|^l Y_L(\widehat{ {\bf k}+{\bf G} }) 
\end{aligned}
\tag{I59}
$$
（$\bar{J}_l$ の定義を思い出すこと。
コードでは `cy(lm)*yl(lm)`
$=Y_L(\widehat{ {\bf k}+{\bf G} })$ である。）
`mkjp.f90` で `pjyl_` を検索せよ。
ここで Bessel 関数は PW の展開に現れる。式 (I45) を見よ。
$R^{JJ}(|{\bf k}+{\bf G}|,|{\bf k}+{\bf G}'|,l)$
は次で与えられる。
$$\begin{aligned}
R^{JJ}(|{\bf k}+{\bf G}|,|{\bf k}+{\bf G}'|,l)
= \int_0^R  r^2 \bar{J}_l(|{ {\bf k}+{\bf G} }|r)
\bar{J}_l(|{ {\bf k}+{\bf G}'}|r) dr ,
\end{aligned}
\tag{I60}
$$
これはロンスキアン（`wronkj`）で、次の式により計算できる。
$$
\begin{aligned}
R^{JJ}(\kappa_{\rm A},\kappa_{\rm B},l)
&=\int_0^R \!\!r^2 \bar{J}_l(\kappa_{\rm A} r) \bar{J}_l(\kappa_{\rm B} r) dr \\
&=\frac{R^2}{\kappa_{\rm A}^2-\kappa_{\rm B}^2}\Big(
\bar{J}_l(\kappa_{\rm A} r) \frac{d \bar{J}_l(\kappa_{\rm B} r)}{dr} \\
&\qquad -\frac{d \bar{J}_l(\kappa_{\rm A} r)}{dr}\bar{J}_l(\kappa_{\rm B} r)
\Big)\Big|_{r=R} = {\tt -fjj}
\end{aligned}
\tag{I61}
$$

（`mkjp.f90` では `-fjj` を使う。）
`vcoulq_4` では、式 (I55) の第 2 項、第 3 項、第 4 項のオンサイト部分を、同じ動径メッシュをもつ原子の組ごとにまとめて足す（`eee = 0` の場合。`eee < 0` では `vcoul_termb` から来る）。

```fortran
vcoul(nbloch+ig1,nbloch+ig2) = vcoul(nbloch+ig1,nbloch+ig2) &
  + cPhi * sum( dconjg(pjyl_(1:lm2x,ig1)) * pjyl_(1:lm2x,ig2) &
    * ((fpi/(absqg2(ig1)-eee)+fpi/(absqg2(ig2)-eee))*fjj(llx(1:lm2x)) + radsig(llx(1:lm2x))) )
```

ここで `cPhi` は位相因子 $e^{-i({\bf k}+{\bf G}_1){\bf R}}e^{i({\bf k}+{\bf G}_2){\bf R}}$ の組の原子についての和、
`fpi` は $4\pi$、`absqg2(ig)` は $|{\bf k}+{\bf G}|^2$、`-eee` は $|E|$、`radsig` は式 (I62) の $\sigma$ 型の積分（`sigintpp`）である。

##### 第 4 項

式 (I55) の最後の項は、本質的に
$\langle B| v| B\rangle$ と同じやり方で計算できる。ただし
Sec.([§](#b-v-b-part)) に現れた $B_{ {\bf R}{l}\mu}(r)$ の代わりに Bessel 関数を使う。
そこで積分
`rojp` と `sgpp` を次のように定義する（`mkjp.f90` の中。`sgpp` はいまは配列として格納しない。上の `radsig` を見よ）。
$$
\begin{aligned}
{\tt rojp}({\bf R}L)&= {\tt pjyl\_}\, \exp(i ({\bf q}+{\bf G}) {\bf R})\,  \rho^l(\bar{J}_l), \\
{\tt sgpp}({\bf R}L,{\bf G},{\bf G}')&= {\tt pjyl\_^*(ig1)} \exp(-i ({\bf q}+{\bf G}) {\bf R}) \\
&\qquad\times {\tt pjyl\_(ig2)} \exp(i ({\bf q}+{\bf G}') {\bf R})\, {\tt radsig}, \\
{\tt radsig}&=\sigma^l(\bar{J}_l(|{\bf q}+{\bf G}|r),\bar{J}_l(|{\bf q}+{\bf G}'|r)).
\end{aligned}
\tag{I62}
$$
このように、`rojp` と `sgpp` は展開の係数と
動径積分からなる。
第 4 項のオフサイト部分は $\sum$ `rojp`$^*$ `strx` `rojp`（`vcoulq_4` の `zmm`）である。

### $\langle P|B \rangle$ の部分

$$
\begin{aligned}
\langle {P}^{\bf k}_{ {\bf G} } |v| B^{\bf k}_{ {\bf R}{L}\mu} \rangle
&=\langle  e^{i({\bf k}+{\bf G}){\bf r}} |v|B^{\bf k}_{ {\bf R}{L}\mu} \rangle \\
&\quad -\sum_{ {\bf R}'{L'} } \langle P^{ {\bf k}+{\bf G} }_{ {\bf R}'{L'} }|v| B^{\bf k}_{ {\bf R}{L}\mu} \rangle  \\
&={\tt fouvb}({\bf G},{\bf R}L \mu) - {\tt sgpb}({\bf G},{\bf R}L \mu) \\
&\quad -\sum_{ {\bf R}'L'} {\tt rojp}^*\, {\tt strx}\, {\tt rojb}
\end{aligned}
\tag{I63}
$$
第 1 項は `fouvb(ngc,nxx,nlxx,nbas)` に入る。これは `mkjp_4`（`mkjp.f90`）が `rojp`、`sgpb` と一緒に作り、`vcoulq_4` で使う。
$v$ は PW に対して対角なので、
これは式 (I57) と同様のやり方で評価できる。
第 2 項は、前節と同じやり方で
`rojb` と `sgpb` から評価できる。オンサイト部分 ${\bf R}'={\bf R}$ が `sgpb` であり、オフサイト部分は $\sum$ `rojp`$^*$ `strx` `rojb`（`vcoulq_4` のブロック `PvB2` の `zmm`）である。

## Offset-$\Gamma$ 法：$\Gamma$ セルで平均した W(k=0) {#offset-gamma}

offset-$\Gamma$ 法は、もともと Ref.[kotani_all-electron_2002] のために Kotani が考案したもので
（Ref.[kotani_quasiparticle_2007] に記述がある）、
我々の論文で精度のよい $GW$
計算を行う鍵であった。
これは Eqs.(I15) と (I20) の ${\bf k}$ 積分のためのものである。
これらの被積分関数は ${\bf k}\to 0$ で発散する。
もとの offset $\Gamma$ 法は、対称性の高い系ではよく働く。
しかし、対称性の低い系に適用すると問題となりうる。
被積分関数の異方的な発散を
正確に扱えない場合があるからである。

ここでは、$W({\bf k},\omega)$ の異方性を正確に扱う、改良した offset-$\Gamma$ 法を示す。
以下では、ここで問題にするのは ${\bf k}$ 積分なので、
簡単のため $W_{\mu \nu}({\bf k},\omega)$ の代わりに
$W({\bf k})$ という表記を用いる（添字と $\omega$ を省く）。

${\bf k}$ メッシュ上の離散和を用いて
$\int_{\rm BZ} f({\bf k}) d^3k$ を計算する式を与える。
ここで $f({\bf k})=G({\bf q}-{\bf k}) \times W({\bf k})$ である。
${\bf k}$ メッシュには次を用いる。
$$
{\bf k}(i_1,i_2,i_3) = 2 \pi \Big(\frac{i_1}{N_1} {\bf b}_1 
+ \frac{i_2}{N_2} {\bf b}_2 + \frac{i_3}{N_3} {\bf b}_3\Big),
\tag{I64}
$$
ここで ${\bf b}_1, {\bf b}_2$, ${\bf b}_3$ は基本逆格子
ベクトルである（Ref.[kotani_quasiparticle_2007] の Eq.(47) と同じ）。
第 1 BZ は $N=N_1 \times N_2 \times N_3$
個の微小セルに分割される（$i_1=0,1,... N_1-1$。$i_2$ と $i_3$ も同様）。
$\Gamma$ 点を含む微小セルを
$\Gamma$ セルと呼ぶ [freysoldt_dielectric_2007]。
主な問題は、$\Gamma$ セルの寄与をどう評価するかである。
$f({\bf k})$ の発散部分は $\approx$ （${\bf k}$ の解析関数）
$/({\bf k}^{\rm T}{\bf L}{\bf k})$ のようにふるまう。
ここで ${\bf k}^{\rm T}$ は ${\bf k}$ の転置を表し、${\bf L}$
は $3\times3$ の Hermite 行列である [friedrich_efficient_2010]。
上の（${\bf k}$ の解析関数）のうち ${\bf k}$ の奇の部分は無視する。
${\bf k}=0$ のまわりの積分に寄与しないからである。
したがって、発散部分が ${\bf k}\to 0$ で
$f({\bf k})=\sum_L \frac{f_L Y_L(\widehat{ {\bf k} })}{|{\bf k}|^2}$ のようにふるまう $f({\bf k})$ の積分を考えれば十分である。
ここで $L\equiv(l,m)$ の $l$ は偶数に限られる。
積分を次の式で評価する。
$$\begin{aligned}
\int_{\rm BZ} f({\bf k}) d^3k \approx \frac{1}{N}\sum^{ {\bf k}\ne 0} f({\bf k})
 + \sum_L f_L w_L + \frac{1}{N} \tilde{f},
\end{aligned}
\tag{I65}
$$
この式は Ref.[freysoldt_dielectric_2007] で導入された。
ここで重み $w_L$ は、$\Gamma$ セルの中での ${\bf k}\to 0$ における $f({\bf k})$ の
発散部分の寄与を取り入れるように、
以下のやり方で決める。
$\tilde{f}$ は ${\bf k}\to 0$ での $f({\bf k})$ の定数部分である。

$w_L$ を決めるには、Ref.[freysoldt_dielectric_2007] で与えられた手順の代わりに
次の手順を用いることができる。
まず補助関数を導入する。
$$\begin{aligned}
F_L({\bf k}) =\sum_{ {\bf G} } \frac{\exp(- \alpha
 |{\bf k}-{\bf G}|^2)Y_L(\widehat{ {\bf k}-{\bf G} })}{|{\bf k}-{\bf G}|^2}.
\end{aligned}
\tag{I66}
$$
これは offset-$\Gamma$ 法で用いた補助関数の
一般化である（当時は $F_{00}$ だけを用いた [kotani_quasiparticle_2007]）。
通常は $\alpha \to 0$ の極限をとるか、代わりに十分小さな $\alpha$ を用いる（`[gw] alpha_offG`、`mkqg.f90` で読む。`m_q0p.f90`）。
Eq. (I65) を $F_L({\bf k})$ に適用する。
すると Eq. (I65) の左辺は
厳密に評価できる（厳密な値は $L=(0,0)$ 以外ではゼロである）。
一方、Eq. (I65) の右辺の第 1 項と
第 3 項は数値的に評価できる。
さらに、$F_L({\bf k})$ に対する $f_{L'}$ は $L'=L$ のとき 1、
それ以外ではゼロであることがわかっている。
したがって、任意の $L$ について $F_L({\bf k})$ に対して Eq. (I65) が厳密に成り立つように、
Eq. (I65) の $w_L$ を決めることができる。

Eq. (I65) を $f({\bf k})=G({\bf q}-{\bf k}) \times W({\bf k})$ に適用する。
そして、$W({\bf k})$ の解析的な部分に加えて、最も強く発散する項だけをとる
近似を行う。すなわち、${\bf k}\to 0$ で
$$\begin{aligned}
W_{\mu \nu}({\bf k}) \sim \widetilde{W}_{\mu \nu}({\bf 0}) +
 \frac{4 \pi}{ {\bf k}^{\rm T} {\bf L} {\bf k} }
 \delta_{1\mu}\delta_{1\nu}
\end{aligned}
\tag{I67}
$$
を用いる。
$\mu=1$ または $\nu=1$ のとき $\widetilde{W}_{\mu \nu}({\bf 0})=0$ である。
Eq. (I67) の近似で何が無視されているかは、
Ref.[friedrich_efficient_2010] の Eq.(36) を参照。

こうして最終的に次を得る。
$$\begin{aligned}
\int_{\rm BZ} d^3k G({\bf q}-{\bf k}) W({\bf k}) 
\approx \overline{\sum G({\bf q}-{\bf k}) W({\bf k})},\end{aligned}
\tag{I68}
$$
ここで右辺は次のように定義する。
$$\begin{aligned}
 &\overline{\sum G({\bf q}-{\bf k}) W({\bf k})}  \\
 &\equiv 
\frac{1}{N} \sum_{ {\bf k}\ne 0} G({\bf q}-{\bf k}) W({\bf k}) + \frac{1}{N} G({\bf q})\overline{W}({\bf 0}),\\
 &\overline{W}({\bf 0})\equiv N\sum w_L W_L + \widetilde{W}({\bf 0}).
\end{aligned}
\tag{I69}
$$
ここで $\overline{W}({\bf 0})$ は
$\Gamma$ セルでの $W$ の平均である。
この $\overline{W}({\bf 0})$ を用いれば、
積分を離散的な ${\bf k}$ メッシュ上の和と同じように評価できる。
行列 ${\bf L}$ が与えられると
（${\bf L}$ の計算法は次の段落で述べる）、
非解析的な（ただし発散しない）関数
${ {\bf k}^{\rm T}{\bf L}{\bf k} }/|{\bf k}^2|$
を球面調和関数で展開する。そして、
与えられた ${\bf L}$ に対して Ref.[friedrich_efficient_2010] に示されたやり方で $W_L$ を計算する。
Eq. (I67) の近似と組み合わせた離散的な ${\bf k}$ メッシュによる積分の精度は、
${\bf k}$ メッシュの大きさを変えて計算することで
評価できる。\

コードでは $L=(0,0)$ の $w_L$ だけを用いる（$L\ne(0,0)$ では $w_L=0$。`m_q0p.f90` の `lxklm = 0`。`m_q0p.f90` は重み `wklm` を書き出す）。
すべての $w_L$ を用いると、$W(\omega)$ の実部が単調でなくなった（因果律に反する）。これは $L\ne(0,0)$ の $w_L$ から来ていた。
金属の $\omega=0$ では ${\bf L}$ は非常に大きく（事実上無限大）、$1/|{\bf k}|^2$ のふるまいはない。

残る問題は、Eq. (I67) の行列 ${\bf L}$ を
どう計算するかである。二つの方法がありうる。
一つは Ref.[friedrich_efficient_2010] で用いられた ${\bf k}\cdot {\bf p}$ 法（摂動）である。
もう一つは、${\bf k}=0$ の近くのいくつかの
${\bf k}$ 点で ${\bf L}$ を数値的に計算する方法である。ここでは後者を用いる。
系の点群対称性により、${\bf L}$ は単位胞の対称性に対する
不変テンソル $\mu_{ij}^g$ の線形結合で表せる。
$$\begin{aligned}
L_{ij}(\omega)=  \sum_{g=1}^{N_g} a_g(\omega) \mu_{ij}^g,
\end{aligned}
\tag{I70}
$$
ここで $g$ は不変テンソルの添字である。
$g$ の数 $N_g$ は、
1（立方対称）から 6（対称性なし）までとりうる。
係数 $a_g(\omega)$ は、
${\bf k}=0$ のまわりの $\{ {\bf k}_{0i}\}$ 点で計算した誘電関数 ${\hat{ {\bf k} }_{0i} }^{\rm T}{\bf L}\hat{ {\bf k} }_{0i}$
から決めることができる。ここで $\{ {\bf k}_{0i}; i=1,N_g\}$ は offset-$\Gamma$ 点の組である。
offset-$\Gamma$ 点は、$\hat{ {\bf k} }_{0i}^{\rm T}{\bf L}(\omega)\hat{ {\bf k} }_{0i}$ から $a_g(\omega)$ への
変換行列が数値的に縮退しないように選ぶ。長さ $|{\bf k}_{0i}|$（`[gw] Q0Pchoice`、または直接与える `deltaq_scale`）
は十分小さく選べるが、$\Gamma$ セルでの $W({\bf k})$ の平均として
数値誤差が出ないようにする。
ここで示した改良 offset-$\Gamma$ 法は、
$\hat{ {\bf k} }_{0i}^{\rm T}{\bf L}(\omega)\hat{ {\bf k} }_{0i}$ がバンド内遷移の寄与を
含む限り、金属の場合にも適用できる。

## ecalj における tetrahedron method のデータ構造 {#tetrahedron-data-structure}

分極関数の tetrahedron の重み、Eqs. (I18a) は、module `m_tetwt`（`m_tetwt.f90`）の `gettetwt` が作る。

```fortran
call gettetwt(q, iq, isp_k, isp_kq, ekxx1, ekxx2, nband=nband, ikbz_in=k_lo, fkbz_in=k_hi)   ! in x0kf_v4h.f90
```

- `q`, `iq`: ベクトル ${\bf q}$ とその添字。
- `isp_k`, `isp_kq` ($\alpha$, $\beta$): ${\bf k}$ と ${\bf k}+{\bf q}$ のスピン（横スピン感受率では異なる）。
- `ekxx1(nband,nqbz)`, `ekxx2(nband,nqbz)`: BZ のすべての ${\bf k}$ についての、${\bf k}$ と ${\bf k}+{\bf q}$ での固有値（`readeval`）。
- `ikbz_in`, `fkbz_in`: 扱う ${\bf k}$ の範囲（${\bf k}$ のループはメモリと並列化のために分割する）。
- そのほかの入力は module から来る。`m_freq` のヒストグラムのビン `frhis(1:nwhis+1)` と `npm`（`npm = 2` では負の振動数の部分も足し込む。時間反転対称性がないときに要る）、
  `m_read_bzdata` の四面体（`ntetf`, `idtetf`, `ib1bz`）、`m_core_state` のコア準位、`m_ReadEfermi` の $E_{\rm F}$、`m_readgwinput` の `ebmx`, `nbmx`, `mtet` である。
  `[gw] t_tetrakbt` では、有限温度の占有の Fermi 準位を用いる（[kBT](./kBT)）。

出力は module 変数に入る。

- `nbnb(k,jpm)`: ${\bf k}$ で重みがゼロでない組 $(n,n')$ の数。`n1b(ibib,k,jpm)`（占有、コアは価電子の後）と `n2b(ibib,k,jpm)`（非占有）がそのバンドの添字である。
- `ihw(ibib,k,jpm)`, `nhw(ibib,k,jpm)`: その組について、重みがゼロでない最初のヒストグラムのビンとビンの数。
  `jhw(ibib,k,jpm)` は重み `whw(jhw : jhw+nhw-1)` を指す。

これらは `x0kf_v4h.f90` で次のように使われる（簡略化したもの）。

```fortran
do jpm = 1, npm
  do ibib = 1, nbnb(k,jpm)
    it  = n1b(ibib,k,jpm)                      ! n  at k
    itp = n2b(ibib,k,jpm)                      ! n' at k+q
    ! M(I,J) = conjg(zmel(I,it,itp)) * zmel(J,it,itp)   (zmel from build_zmel)
    do iw = ihw(ibib,k,jpm), ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1
      imagweight = whw(jhw(ibib,k,jpm)+iw-ihw(ibib,k,jpm))
      ! Im Pi(:,:,iw,jpm) += M(:,:)*imagweight
    enddo
  enddo
enddo
```

これが Eq. (I18a) の $n,n'$ と ${\bf k}$ についての和である。

## 非局所的な量の Fourier 変換

MTO 基底での静的な自己エネルギーは `m_rdsigm2`（`rdsigm2.f90`）で Fourier 変換する。`rdsigm2` が `sigm.<sname>` を全 BZ に展開し、`fftz3` が実空間の行列 `hrr` を与え、`bloch2` がそれを任意の ${\bf k}$ に戻す。
これは 3 次元の FFT のためのものである。
通常の FT は次で与えられる。
$$\begin{aligned}
f({\bf T})=\sum_{\bf k}f({\bf k}) \exp(i {\bf k}{\bf T}),
\end{aligned}
\tag{I71}
$$
ここで $\{ {\bf k}\}$ は規則的なメッシュ点上にある。
その要素の総数は $N_1 \times N_2 \times N_3$ である
（この数は $\{ {\bf T}\}$ の数と同じ）。
${\bf k}$ 点にも ${\bf T}$ 点にも周期性があることに
注意する。周期性のため、$\{ {\bf k}\}$ の範囲は
一意でなく、$\{ {\bf T}\}$ の範囲も同様である。

${\bf T}-{\bf T}'$ に依存する非局所的な量を考える。すると
$$\begin{aligned}
f({\bf R}{\bf T},{\bf R}'{\bf T}')=\sum_{\bf k}
f_{ {\bf R}{\bf R}'}({\bf k}) \exp(i {\bf k}({\bf T}-{\bf T}')),
\end{aligned}
\tag{I72}
$$
となる。
実際の計算（`bloch2` で扱う自己エネルギーの静的な版）では、
まず規則的なメッシュ点の ${\bf k}$ 上で $f_{ {\bf R}{\bf R}'}({\bf k})$ を計算する。次に
その実空間表示
$f({\bf R}{\bf T},{\bf R}'{\bf T}')$ を得る必要がある。周期性のため、
ありうる $|{\bf T}-{\bf T}'|$ の選び方に
任意性がある。$\bar{\bf T}={\bf T}-{\bf T}'$ を導入すると、
並進対称性により Eq. (I72) は
$f({\bf R}\bar{ {\bf T} },{\bf R}'0)=\sum_{\bf k}
f_{ {\bf R}{\bf R}'}({\bf k}) \exp(i {\bf k}\bar{ {\bf T} })$ と書ける。

妥当な選び方は、
$|{\bf R}-{\bf R}'+\bar{ {\bf T} }|\le\eta_{\rm FTmax}$ を満たす $\bar{ {\bf T} }$ を許すことである。
ここで $\eta_{\rm FTmax}$ は、
許される $\{\bar{\bf T}\}$ の数が $N_1 \times N_2 \times N_3$ になるように選ぶべきである。
しかし、縮退のためにそれができない場合がある。
許される $\{\bar{\bf T}\}$ の中で $|{\bf R}-{\bf R}'+\bar{ {\bf T} }|$ が最大の値をとるものについて、
いくつかの $\bar{ {\bf T} }$ がありうる（縮退しているという）。
その場合は、そのような $\bar{ {\bf T} }$ に分数の重みを与える必要がある。

コードでは $\bar{\bf T}$ のリストは自動的に作られる。サイトの各組 $({\bf R},{\bf R}')$ と $N_1\times N_2\times N_3$ セルの各 ${\bf T}$ について、
`m_gennlat_sig`（`m_gennlat.f90`）が最も短い等価なベクトル $|{\bf R}-{\bf R}'+\bar{\bf T}|$ を保持し（`nlat`、範囲は `nlatS:nlatE`）、
`bloch2` が縮退したものについて位相因子を平均する（等しい分数の重み）。

## Brillouin ゾーンでの自己エネルギーの内挿 {#sigma-interpolation}

ここでは、規則的なメッシュ点 ${\bf k}(i_1,i_2,i_3)$ でだけ計算した $V^{\rm xc}_{ {\bf k} }$ から、
任意の ${\bf k}$ での $V^{\rm xc}_{ {\bf k} }$ を与える
内挿の手順を示す。
この内挿は、$\{ {\bf k}_{0i}\}$ での $W(\omega)$ を必要とする
offset-$\Gamma$ 法で使われる。
この $W(\omega)$ を計算するには、
規則的なメッシュ点 ${\bf k}(i_1,i_2,i_3)$ だけでなく
${\bf k}(i_1,i_2,i_3)+{\bf k}_{0i}$ でも
固有関数と固有値が要る。
この内挿はエネルギーバンドを描くのにも
役に立つ。
内挿の要点は、次のように $V^{\rm xc}$
を局在性の高い MTO で実空間展開することである。

QSGW の反復の終わり（`hqpe_sc`）で、規則的な ${\bf k}$ メッシュ点上の行列要素
$\langle {\Psi_{ {\bf k}n} }| {\it \Delta} V^{\rm xc}_{\bf k}| {\Psi_{ {\bf k}m} }\rangle$ を得る。ここで
${\it \Delta} V^{\rm xc}_{\bf k}=V^{\rm xc}_{\bf k}-V^{\rm xc,LDA}_{\bf k}$ である。
次に、これを APW と MTO 基底での表示に変換する。
$$\begin{aligned}
\langle \chi^{ {\bf k} }_a| {\it \Delta} V^{\rm xc}_{\bf k}| \chi^{ {\bf k} }_b \rangle
= \sum_{n,m} \left(z^{-1}\right)^*_{an}\langle {\Psi_{ {\bf k}n} }|{\it \Delta} V^{\rm xc}_{\bf k}|
{\Psi_{ {\bf k}m} }\rangle z_{bm}^{-1},  \\
\end{aligned}
\tag{I73}
$$
ここで基底を指定する簡略化した基底の添字 $a$ を用いる
（MTO なら ${\bf R}L{j}$、APW なら ${ {\bf G} }$）。
したがって $\chi^{ {\bf k} }_a$ は Eq. (I1) の APW または MTO を表す。
$z_{na}$（簡単のため ${\bf k}$ を省く）
は ${\bf k}$ での固有関数の係数、すなわち
Eq. (I1) の $z^{ {{\bf k} }n}_{ {\bf R}L{j} }$ と $z^{ {\bf k}n}_{ {\bf G} }$ をまとめたものを表す。
この $z_{an}$ は、固有関数（バンドの添字 $n$）と APW・MTO 基底（基底の添字 $a$）を
結ぶ変換行列とみなせる。

${\it \Delta} V^{\rm xc}$ の実空間表示を得るには、
BZ の ${\bf k}$ について周期的な、Bloch 和をとった局在軌道からなる基底で
展開した表示が要る。しかし、Eq. (I73) の APW はそうなっていない。
この問題を克服するため、
MTO に関係する行列要素、すなわち
$a$ と $b$ が MTO を指定する要素
$\langle \chi^{ {\bf k} }_a| {\it \Delta} V^{\rm xc}_{\bf k}| \chi^{ {\bf k} }_b \rangle$
だけをとる近似を用いる。
これは、${\it \Delta} V^{\rm xc}$ のうち APW に関係する部分を MTO の基底に射影することを意味する。
精度を確かめるには数値的な試験が要るが、
${\it \Delta} V^{\rm xc}$ の主要部分が MTO でよく展開できる限り、
この近似は妥当でありうる。
こうして、$\langle \chi^{ {\bf k} }_a| {\it \Delta} V^{\rm xc}_{\bf k}| \chi^{ {\bf k} }_b \rangle$ の MTO 部分から
Fourier 変換によって、MTO で展開した ${\it \Delta} V^{\rm xc}$ の
実空間表示を得る。
そして、逆 Fourier 変換によって任意の ${\bf k}$ で内挿した ${\it \Delta} V^{\rm xc}$ が得られる。
局在性の高い MTO を用いるので、
この内挿は FP-LMTO-QSGW での以前のもの [kotani_quasiparticle_2007]
より安定である。
Ref.[kotani_quasiparticle_2007] の Sec.II-G で与えた
複雑な内挿手順は
もはや必要ない。

計算時間を減らすため、
行列要素 $\langle {\Psi_{ {\bf k}n} }|{\it \Delta} V^{\rm xc}_{\bf k}| {\Psi_{ {\bf k}m} }\rangle$ は、
固有値が $E_{\rm MAX}^\Sigma$（`[gw] emax_sigm`）より小さい状態までしか計算しない。
行列要素の高エネルギー部分は
対角であると仮定し、その値は計算した対角要素の
エネルギー重み付き平均（`hqpe_sc` の `ESEAVR`）で与える。

## GW 計算のための結晶の対称性

**空間群の操作。** `lmf` は空間群の操作を求め（既定では spglib を用いる。[lmf: spglib による対称性](./lmf#spglib-symmetry)）、
`--jobgw` のとき、点群の部分 `symops(3,3,ig), ig=1,ngrp` と並進 `ag` を `__HAMindex` に書き出す（`m_hamindex.f90` の `m_hamindexW_init`、`main_lmf.f90` で呼ぶ）。
たとえば GaAs では操作は 24 個である。

GW のプログラムは module `m_hamindex`（`symops`, `ngrp`, `invgx`）を通してこれらを読む。
`m_zmel.f90` の `mptauof_zmel` は `mptauof`（`m_symderive.f90`）を呼ぶ。`mptauof` は、各操作について、各原子が写される原子（`miat`）と格子並進（`tiat`）を与える。
`mptauof` の冒頭のコメントを参照。

操作は事前に `lmchk` で確かめられる。${\bf k}$ 点については、時間反転対称性のために反転が加えられる
（実の局所ポテンシャルでは $|\psi_{-{\bf k}}({\bf r})|^2=|\psi_{\bf k}({\bf r})|^2$）。軌道モーメントがあるとき（LDA+U、SOC）は注意が要る。

対称性は `ctrlg.<sname>.toml` の最上位の `symgrp` で制御する。`"find"`（既定）は結晶構造の対称性を用いる。
生成元を与えて、意図的に低い対称性を用いることもできる（たとえば NiO は実際には (111) 方向に歪んでいるが、歪みなしで計算してもよい）。
**反強磁性の対称性**: `[[site]]` のサイトの組に `af = ±k` を書き、`symgrpaf = "find"` とする。上向きスピンのサイトを下向きスピンのサイトに写す操作
（スピン反転を伴う）が加えられる。`pwmode = 11` と一緒に使う。

## ecalj における球面調和関数と実調和関数

ecalj では、fortran のコードは主に、通常の球面（複素）調和関数 $Y_{lm}(\hat{\bf r})$ ではなく、
実調和関数 $y_{lm}(\hat{\bf r})$ に基づく。
固有関数の係数などは、たとえば $l=2$ では
$(m=-2, m=-1, m=0, m=1,m=2)$ の順に並ぶ。
たとえば `LMXA=4` では `(4+1)**2=25` 個の調和関数があり、
次の順に並ぶ。
$y_{00},y_{-11},y_{01},y_{11},y_{-22},y_{-12},... ,y_{22},y_{-33},...y_{33},
y_{-44},...y_{44}$.
実調和関数（$r^l$ を掛けたもの）は `sylm`（`sylm.f90`）が `sylmnc` の規格化 `cy` で作る。$l=1$ では順序は $y, z, x$、$l=2$ では $xy, yz, 3z^2-r^2, xz, x^2-y^2$ である（規格化を除いて）。$m\gt 0$ は $\cos m\phi$ 型、$m\lt 0$ は $\sin |m|\phi$ 型である。

$y_{lm}(\hat{\bf r})$ は $Y_{lm}(\hat{\bf r})$ から定義する。
ここで $\hat{\bf r}=(\theta, \phi)$ である。

$$\begin{aligned}
 y_{l0}(\hat{\bf r}) 
  & \equiv  Y_{l0}(\hat{\bf r}). \\
 y_{lm}(\hat{\bf r}) 
  & \equiv  \frac{1}{\sqrt{2} }
           [ (-1)^m Y_{lm}(\hat{\bf r}) + Y_{l-m}(\hat{\bf r}) ]. \\
 y_{l-m}(\hat{\bf r})
  & \equiv  \frac{1}{\sqrt{2}i}
           [ (-1)^m Y_{lm}(\hat{\bf r}) - Y_{l-m}(\hat{\bf r}) ].
\end{aligned}
\tag{I87}
$$
、ここで $m\gt 0$ である。同値な形では、
$$\begin{aligned}
 Y_{l0}(\hat{\bf r}) 
  & \equiv  y_{l0}(\hat{\bf r}). \\
 Y_{lm}(\hat{\bf r}) 
  & \equiv  \frac{(-1)^m}{\sqrt{2} }
           [ y_{lm}(\hat{\bf r}) + iy_{l-m}(\hat{\bf r}) ]. \\
 Y_{l-m}(\hat{\bf r})
  & \equiv  \frac{1}{\sqrt{2} }
           [ y_{lm}(\hat{\bf r}) - iy_{l-m}(\hat{\bf r}) ].
\end{aligned}
\tag{I88}
$$

$Y_{lm}(\hat{\bf r})$ の定義は次のとおりである。
$$\begin{aligned}
 &Y_{lm}(\theta, \phi)
=(-1)^m \left[ \frac{(2l+1)(l-m)!}{4 \pi (l+m)!} \right]^{\frac{1}{2} } P^m_l(\cos(\theta)) e^{i m \phi}, \\
 &P^m_l(x) = \frac{(1-x^2)^{m/2} }{2^l l!}\frac{d^{l+m} \ \ }{dx^{l+m} } (x^2-1)^l
\end{aligned}
\tag{I89}
$$

これらの定義は
(1)A.R.Edmonds, Angular Momentum in quantum Mechanics,
Princeton University Press, 1960,
(2)M.E.Rose, Elementary Theory of angular Momentum,
John Wiley & Sons, INC. 1957 からとった。球面調和関数の定義はこれらの本で同じである。

## 空間群の操作による固有関数と MPB の回転 {#spacegroup-rotation}

### PMT 固有関数の空間群による回転

固有関数を空間群の操作で写す必要がときどきある。

空間群の操作は $(g,{\it \Delta})$ で指定する。これは
$3 \times 3$ の空間回転行列 $g$ と並進ベクトル
${\it \Delta}$ からなる。共型（symmorphic）の場合は ${\it \Delta}= 0$ であることを思い出そう。
この $(g,{\it \Delta})$ は回転 ${\bf r}\to {\bf r}'$ を
$$\begin{aligned}
 &{\bf r}'=g({\bf r})+ {\it \Delta},
\end{aligned}
\tag{I96}
$$
のように行う。これは次と同値である。
$$\begin{aligned}
 &{\bf r}=g^{-1}({\bf r}')+ {\it \Delta}^{-1},\\
 &{\it \Delta}^{-1}=-g^{-1}({\it \Delta})
\end{aligned}
\tag{I97}
$$
さらに、与えられた $(g,{\it \Delta})$ について、どの原子サイトがどの
原子サイトに写されるかを計算する。それは次で与えられる。
$$\begin{aligned}
 &g({\bf R})+{\it \Delta}={\bf R}'+{\it \Delta}{\bf T}_{\bf R}.
\end{aligned}
\tag{I98}
$$
`ecalj` では、配列 `miat(ipos=1,npos)` が基本セル内の ${\bf R}'$ を
与える。
`tiat(3,ipos)` は基本セルベクトルの整数倍の並進
${\it \Delta}{\bf T}_{\bf R}$ を与える。ここで ${\bf R}$ は単位胞の中の原子である。

$(g,{\it \Delta})$ による関数の回転 $F \to F'$ は、
$F'({\bf r}')=F({\bf r})$ となるように定義する。ここで $F'({\bf r}')$ を $g[F]({\bf r}')$ と書く。
したがって
$$\begin{aligned}
g[F]({\bf r})=F(g^{-1}({\bf r})+{\it \Delta}^{-1})
\end{aligned}
\tag{I99}
$$
となる。

これを Eq. (I3) の Bloch 和の関数
$A^{\bf k}_{ {\bf R}u L}({\bf r}) \equiv \sum_{\bf T} 
A_{ {\bf R}u L}({\bf r-R-T}) \exp(i {\bf kT})$ に適用する。
（ここで角運動量の添字 $L$ をとる。$u$ は動径部分の添字である。）
すると
$$
\begin{aligned}
g[A^{\bf k}_{ {\bf R}u L}]({\bf r})&=A^{\bf k}_{ {\bf R}u L}(g^{-1}({\bf r})+{\it \Delta}^{-1}) \\
&=\sum_{\bf T}A_{ {\bf R}u L}(g^{-1}({\bf r})+{\it \Delta}^{-1} -{\bf R}-{\bf T}) \exp(i {\bf kT})  \\
&=\sum_{\bf T}A_{ {\bf R}u L}\big(g^{-1}({\bf r}+g({\it \Delta}^{-1}) -g({\bf R})-g({\bf T}))\big) \\
&\qquad\times \exp(i {\bf kT})  \\
&=\sum_{\bf T}\sum_{L'}  A_{ {\bf R}u L'} ({\bf r}+g({\it \Delta}^{-1}) -g({\bf R})-g({\bf T})) \\
&\qquad\times D_{L'L}(g) \exp(i {\bf kT})  \\
&=\sum_{\bf T}\sum_{L'}  A_{ {\bf R}' u L'}({\bf r}-{\bf R}'-{\it \Delta}{\bf T}_{\bf R}-g({\bf T})) \\
&\qquad\times D_{L'L}(g)\exp(i {\bf kT})  \\
&=\sum_{ {\bf T}''} \sum_{L'}  A_{ {\bf R}' u L'}({\bf r}-{\bf R}'-{\bf T}'')\, D_{L'L}(g) \\
&\qquad\times \exp(i g({\bf k}){ {\bf T} }'') \\
&\qquad\times \exp(-i g({\bf k}){\it \Delta}{\bf T}_{\bf R})  \\
&= \sum_{L'}  A^{g({\bf k})}_{ {\bf R}' u L'}({\bf r})\, D_{L'L}(g) \exp(-i g({\bf k}){\it \Delta}{\bf T}_{\bf R}),
\end{aligned}
\tag{I100}
$$
となる。ここで $R$ と $R'$ のサイトは等価なサイトのはずなので、$A_{ {\bf R}' u L}({\bf r})=A_{ {\bf R}u L}({\bf r})$ を用いた。
また打ち消し ${\it \Delta}=-g({\it \Delta}^{-1})$ を用いた。
$D_{L'L}(g)$ は $g$ に対する実球面調和関数の回転行列である。

### コード: `rotevec`

ルーチンは `rotwave.f90` にある。`rotevec` は PMT 基底（MTO と APW の部分）での固有関数の係数を回転し、`rotmatMTO` は MTO 部分の回転行列を与える。

この Eq. (I100) により、固有関数の MTO 部分の係数が
空間群の操作 $g$ で回転される（写される）。
これは `subroutine rotevec`（固有関数を回転するルーチン）で次のようにコードされている。

    ...
     phase = [(exp(-img2pi*sum(qtarget*tiat(:,ibas,igg))),ibas=1,nbas)]
     do iorb=1,norbmto !orbital-blocks are specified by ibas, l, and k.
       ibas = ibastab(iorb)
       l   = ltab(iorb)
       k   = ktab(iorb)
       init1 = offl(iorb)+1
       iend1 = offl(iorb)+2*l+1
       init2 = offlrev(miat(ibas,igg),l,k)+1
       iend2 = offlrev(miat(ibas,igg),l,k)+2*l+1
       evecout(init2:iend2,:)= matmul(dlmm(-l:l,-l:l,l,igg),evec(init1:iend1,:))*phase(ibas)
     enddo
    ...

ここで

> `ndimh = nlmto + napw`: PMT ハミルトニアンのサイズ\
> `evec(1:nlmto,1:nband)` : MTO 上の固有関数の係数\
> `igg`: 空間群の操作 $g$ を指定する添字\
> `qtarget`: qtarget=$g({\bf k})$\
> `ibas`: 基本セル内のサイトの添字。\
> `iorb` は MTO 部分のブロックの添字（$m=-l,-l+1,... ,l$ の $2l+1$ 個の要素）。\
> `l` は角運動量の添字。L の l。\
> `k` は動径の添字。\
> `miat`: サイトの添字 ibas は `miat(ibas,igg)` に写される\
> `tiat`: ${\it \Delta}{\bf T}_{\bf R}$\
> `offl`: `iorb` ブロックのオフセット\
> `offlrev`: (ibas,l,k) ブロックのオフセット\
> ブロック `init1:iend1` はブロック `init2:iend2` に写される。

`rotevec` を呼ぶには、あらかじめ `readhamindex()` を呼んで、`m_hamindex` module の変数
`ltab,ktab,offl,offlrev,...` を設定しておく必要がある。

軌道の添字の表（`ibas_tbl`, `l_tbl`, `k_tbl`, ...）は、ファイル `@MNLA_CPHI`
（`sugw.f90` が書く）から `readmnla_cphi`（`readeigen.f90`）が読む。
`@MNLA_CPHI` の最後の列は軌道ブロックの添字である。
軌道ブロックは $l$ と $k$ で指定され、そのサイズは $2l+1$ である。
$i$ の数はちょうど `MNLA_CPHI` の行数である。\

`subroutine rotevec`（固有関数を回転するルーチン）には APW の部分もある。
式は次のとおりである。
$$
\begin{aligned}
g[\exp(i ({\bf k}+{\bf G}){\bf r})]  
&=\exp\big(i ({\bf k}+{\bf G}) (g^{-1}({\bf r})+{\it \Delta}^{-1})\big) \\
&=\exp\big(i (g({\bf k})+g({\bf G})) ({\bf r}+g({\it \Delta}^{-1}))\big)  \\
&=\exp\big(i g({\bf k}+{\bf G}) {\bf r}\big) \exp\big(i ({\bf k}+{\bf G}) {\it \Delta}^{-1}\big)  \\
&= \exp\Big({i \big(\bar{g}({\bf k}) + g({\bf G})+ (g({\bf k})-\bar{g}({\bf k}))\big) {\bf r} }\Big) \\
&\qquad\times\exp\big(-i\, g({\bf k}+ {\bf G}){\it \Delta} \big).
\end{aligned}
\tag{I101}
$$

（我々の PW (IPW) の定義には「規格化因子がない」。）
ここで $g({\bf k})$ は第 1 BZ の $\bar{g}({\bf k})$ に引き戻される。
その差は $g({\bf k})-\bar{g}({\bf k})$ である。
`barG` を $\bar{G}=g({\bf G})+ (g({\bf k})-\bar{g}({\bf k}))$ と定義する。
したがって次のようにコードされている\
（ここでは ${\bf k}$ の代わりに ${\bf q}$ を用いている。紛らわしくて申し訳ない）。

    ...
      igg: index for space group operation.
    nlmto: the size of MTO part.

      ikt  = getikt(q)       !index for q 
      ikt2 = getikt(qtarget) !index for bargq = \bar{g}(q) 
      do ig = 1,napw_in !number of APW 
        qpg  = q + matmul( qlat(:,:),igv2(:,ig,ikt))  ! q+G
        qpgr = matmul(symops(:,:,igg),qpg)            !g(q+G)
        nnn= nint(matmul(platt,qpgr-qtarget)) !integer sets for barG= g(G)+ g(q)-bargq
        ig2 = ngvecprev(nnn(1),nnn(2),nnn(3))   ! index for barG
        evecout(nlmto+ig2,:)= evec(nlmto+ig,:) * exp( -img2pi*sum(qpgr*shtvg(:,igg)) )
      enddo
    ...

### 固有関数の MTIPW 展開の空間群による回転

ecalj の GW 計算では、固有関数を MTIPW の形、すなわち
Eq. (I7) のように MT 部分 + IPW 部分で展開し直す。
MTIPW の形で、固有関数を空間群の操作で写す必要がある。
`rotevec` と同様に、`rotmto` と `rotipw`（`rotwave.f90`）がある。
これらは MT 部分と IPW 部分の回転のためのものである。`readeigen.f90` で使われ、保存されていない ${\bf k}$ 点での `__CPHI` と `__GEIG` を、保存されている等価な点から与える。

### Mixed Product Basis の空間群による回転

与えられた ${\bf q},{\bf k},n,n'$ に対して、すべての行列要素 $\langle \Psi_{ {\bf q+k}n'}| \Psi_{ {\bf k}n} E^{\bf q}_\nu \rangle$ が要る。
ここで $E^{\bf q}_\nu$ は Coulomb 相互作用を対角化した mixed product basis である
[friedrich_efficient_2010]。

これは EIBZ の方式についての覚え書きで、EIBZ は今は使っていない（[IBZ と EIBZ](#ibz-and-eibz-scheme)）。`m_rotMPB2` と `x0kf_v4hz_symmetrize` は今のコードにはなく、`rotmto2` と `rotipw2` は `rotwave.f90` に残っている。
EIBZ の手順（Eq. (I95) のあたりを参照）では、MPB を回転する必要がある。
module `module m_rotMPB2 (in m_rotMPB.F)` は回転行列 `zrotm` を作るためのものである。
`zrotm` を用いて `call X0kf_v4hz_symmetrize` で MPB を回転する。そこでは
MPB から $\{E^{\bf q}_\nu\}$ への変換行列 `zcousq` を次のように掛ける\
`  rcxq(:,:,iw,jpm) = matmul(zcousqc,matmul(rcxq_core,zcousq))`（`x0kf_v4hz_symmetrize` のコードを参照）。

`zrotm` は `rotmto2,rotipw2` を呼んで作る。
これらは本質的に `rotmto,rotipw` と似ているが、MT 領域の基底と
格子間領域の基底が固有関数の場合とは異なる。
そこで `ibas_tbl` などを次のように与える。

          use m_pbindex,only: norbt, ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,
         &     max_ibas_tbl,max_l_tbl,max_k_tbl,max_offset_tbl

（`module m_rotMPB2` の中。）

## MLO のための SOC 行列

Wannier 関数の経路（`hsocmat`）は削除した。SOC 行列は今は MLO のために使う（`job_mlo_soc`。[MLO §4](./mlo#mlo-soc)）。

1. `lmf --writeham --socmatrix`（スカラー相対論、`[ham] so = 0`）は、PMT 基底でのハミルトニアン（`__HamiltonianPMT`）と、それとは別に
   SOC 行列 $\langle F_{{\bf k}i}|H_{\rm SOC}|F_{{\bf k}j}\rangle$（`aughsoc` が作る `hammhso(ndimh,ndimh,3)`。三つは $\uparrow\uparrow$, $\downarrow\downarrow$, $\uparrow\downarrow$ のブロックで、スピノルで倍にせず軌道ごと）を `__HamiltonianPMTsoc` に書き出す（`m_bandcal.f90`）。
   ここで $F_{{\bf k}i}$ は PMT 法の基底関数で、$\Psi_{{\bf k}n\sigma}=\sum_i F_{{\bf k}i}\alpha^{\bf k}(i,n,\sigma)$ である。
2. `mlo --socmatrix` はそれをハミルトニアンと一緒に MLO 基底に変換し（`HamRsMLO` の `hammhsor`、`m_mlo_ham.f90` で読む）、
   MLO ハミルトニアンに摂動として加える（スピノルの次元 `nspinor = 2`）。
3. Fermi エネルギーは full-SOC の `lmf` の計算（`[ham] so = 1`、`efermi_soc`）からとる。

MLO バンドの SOC 分裂は、SOC を入れた完全な計算のバンドと比べて確かめるべきである（`job_mlo_soc` は両方を描く）。

## Harris-Foulkner エネルギーと Kohn-Sham エネルギー

LDA/GGA では、自己無撞着に至る途中では、入力密度と出力
密度は同じでない（自己無撞着でない）。
そこで、与えられた入力密度 $n^{\rm in}$ に対して二つの全エネルギー、
Harris-Foulkner エネルギー $E_{\rm Harris}$ と Hohenberg-Kohn エネルギー $E_{\rm HoKohn}$ を定義する。
$$
\begin{aligned}
E_{\rm Harris}& =  E_{\rm k}^{\rm core} + E_{\rm B} \\
 &\quad - V_{\rm es}[n^{\rm Zc}+n^{\rm in},{\bf R}_a] \cdot n^{\rm in}
- V_{\rm xc}[n^{\rm c}+n^{\rm in}] \cdot n^{\rm in} \\
 &\quad + E_{\rm es}[n^{\rm Zc}+n^{\rm in},{\bf R}_a] + E_{\rm xc}[n^{\rm c}+n^{\rm in}],\\
E_{\rm HoKohn}& =  E_{\rm k}^{\rm core} + E_{\rm B} \\
 &\quad - V_{\rm es}[n^{\rm Zc}+n^{\rm in},{\bf R}_a] \cdot n^{\rm out}
- V_{\rm xc}[n^{\rm c}+n^{\rm in}] \cdot n^{\rm out} \\
 &\quad + E_{\rm es}[n^{\rm Zc}+n^{\rm out},{\bf R}_a] + E_{\rm xc}[n^{\rm c}+n^{\rm out}], \\
E_{\rm B} & =  \sum_p^{\rm occupied} \alpha_{p}^{i*} \langle F_i|H^{\rm in}|F_j \rangle \alpha_p^j,
\end{aligned}
\tag{I80}
$$
Ref.[kotani_formulation_2015] の Eqs.(B.1),(B.2) を参照。
$E_{\rm xc}$ についての Eq.(31) が誤っているので、これらは少し誤っている。$E_{\rm xc}$ は
全電荷（原子核を含む）$n^{\rm Zc}+n^{\rm in}$ の汎関数ではなく、
全電子密度 $n^{\rm c}+n^{\rm in}$ の汎関数である。
Ref.[kotani_formulation_2015] の 3 成分形式に基づくと、
$E_{\rm es}$ は原子位置 ${\bf R}_a$ に陽に依存するが、$E_{\rm xc}$ は依存しない。
その依存性は Eq.(14) の多重極変換を通じたものである。
`m_bndfp.f90` で `ham_ehf,ham_ehk` を探すこと。
$E_{\rm Harris}$ は `call m_mkehkf_etot1(sev, eharris)` で与えられ、
$E_{\rm HoKohn}$ は
`call m_mkehkf_etot2(ekinval, eksham)`（`mkehkf.f90`）で与えられる。
$E_{\rm Harris}$ と $E_{\rm HoKohn}$ は `save.*` ファイルに示される。

バンドエネルギー $E_{\rm B}$ は `m_bndfp.f90` の module 変数 `sev` に格納される。
これは `save.*` に `sev` として示される。$E_{\rm HoKohn}$ に含まれる
運動エネルギー
$E_{\rm k}^{\rm core}+E_{\rm B} 
- V_{\rm es}[n^{\rm Zc}+n^{\rm in},{\bf R}_a] \cdot n^{\rm out}
- V_{\rm xc}[n^{\rm c}+n^{\rm in}] \cdot n^{\rm out}$ は、
`mkekin`（`m_bndfp_util.f90`）が `ekinval` として計算する。

`lmf` による LDA/GGA の計算では、`save.*` ファイルは反復ごとに 1 行を含む（`m_lmfp.f90` で書く。エネルギーは、原子についての種の参照エネルギー `eref` の和に対する相対値である。`m_lmfinit.f90` の `ham_seref`）。

    c ehf(eV)=-15730.0982239 ehk(eV)=-15730.0982222 sev(eV)=-15.7737393

これは $E_{\rm Harris}=$-15730.0982239 eV と $E_{\rm HoKohn}=$-15730.0982222 eV、および価電子の
バンドエネルギー $=-15.7737393$ eV を示す。
原理的には、収束したとき $E_{\rm Harris}$ と $E_{\rm HoKohn}$ は厳密に一致するはずである。
その差は数値誤差である。行頭の `c` は「収束した」を意味する。
${\tt h}$ は `atm.*` ファイル（原子密度の重ね合わせ）からの最初の反復を意味する。これらの文字は `m_lmfp.f90` の `svflg='cxhiC'` である（`x`: 反復の最大回数に達した、`i`: 途中の反復、`C`: 緩和が収束した）。

## IBZ と EIBZ の方式 {#ibz-and-eibz-scheme}

以下の IBZ の部分は、コードが使っているものである。自己エネルギーは IBZ の ${\bf k}$ と回転について和をとる（`m_sxcf_sc.f90` の `kx` と `irot`、表は `irkip`）。EIBZ の部分は覚え書きで、今は使っていない（`eibzgen.f90` はあるが呼ばれない。`[gw] EIBZmode` は読むが使わない。`m_eibz`、`x0kf_v4hz_symmetrize`、`zsecsym`、`rotwvigg` は今のコードにない）。再び使うことがあるかもしれない。

Ref.[friedrich_efficient_2010] の III.G. で与えられた EIBZ の手順を、分極関数と自己エネルギーについて述べる。
[IBZ]{.underline}\
まず IBZ を思い出そう。空間群の操作の集合
$S_A=\{A_i |i=1,... N_A\}$ がある（逆格子ベクトルによる並進を除く）。
$A_i$ は空間回転と並進で表される。
この $S_A$ は時間反転操作を含んでよい。${\bf q}'=A_i({\bf q})$ は BZ に属する。
（${\bf q}$ に回転行列を単に掛けたものは、（第 1）BZ への引き戻しが要る場合と要らない場合がある。）
「BZ の中の $\{ {\bf q}\}$」というメッシュ点の集合に対して、$A_i$ はそれらの間の一対一の写像を与える。
「IBZ の中の $\{ {\bf q}\}$」は第 1 BZ の部分集合である。

${\bf q}$ に依存する関数の集合 $\{g_{ {\bf q} }({\bf r}, n)\}$ を考える。
$F[g_{\bf q}]$ の和として定義される $X$ を次のように計算できる。
$$\begin{aligned}
X=\sum_{ {\bf q}\in {\rm BZ} } F[g_{\bf q}]
= \sum_{A_i \in S_A } 
\sum_{ {\bf q}\in {\rm IBZ} } \frac{N_{\bf q} }{N_A} F[g_{A_i({\bf q})}] \end{aligned}
\tag{I91}
$$
ここで IBZ は、BZ のすべてのメッシュ点を生成できる ${\bf q}$ の集合である。
$N_{\bf q}$ は $S_A$ によって ${\bf q}$ から生成されるメッシュ点の数（集合 ${\bf q}^*$ の数）である。
${\bf q}=0$ では $N_{\bf q}=1$ である。一般の ${\bf q}$（最も対称性の低い点）では $N_{\bf q}=N_A$ である。

$A_i[g_{\bf q}]= g_{A_i({\bf q})}({\bf r},n)$ が成り立てば、
$$
\begin{aligned}
X=\sum_{ {\bf q}\in {\rm BZ} } F[g_{\bf q}]
&= \sum_{A_i \in S_A } \sum_{ {\bf q}\in {\rm IBZ} } \frac{N_{\bf q} }{N_A} F[g_{A_i({\bf q})}] \\
&= \sum_{A_i \in S_A } \sum_{ {\bf q}\in {\rm IBZ} } \frac{N_{\bf q} }{N_A} F[A_i[g_{\bf q}]] .
\end{aligned}
\tag{I92}
$$
となる。
これは、$\{g_{ {\bf q} }({\bf r}, n)|{\bf q}\in {\rm IBZ}\}$ だけから $X$ を評価する方法を示している。

（注: 対角化で計算した固有関数は、位相因子を除いて
$A_i[g_{\bf q}]= g_{A_i({\bf q})}({\bf r},n)$ を満たす。位相因子が関係しない場合を考える。）

 \
[EIBZ]{.underline}\
与えられた ${\bf k}$ を保つ $S_A$ の部分群を考える。
この部分群を $S^{\bf k}_A=\{A^{\bf k}_i |i=1,2,...,N_A^{\bf k}\}$ と書く。
この部分群に IBZ と同じ論理を適用すると、次を得る。
$$\begin{aligned}
X=\sum_{ {\bf q}\in {\rm BZ} } F[g_{\bf q}] = 
\sum_{A^{\bf k}_i \in S_A^{\bf k} } \sum_{ {\bf q}\in {\rm EIBZ({\bf k})} }
\frac{N_{\bf q}^{\bf k} }{N_A^{\bf k} } F[A^{\bf k}_i[g_{\bf q}]] \end{aligned}
\tag{I93}
$$
ここで EIBZ$({\bf k})$ は、BZ のすべてのメッシュ点を生成できる ${\bf q}$ の集合である。
$N^{\bf k}_{\bf q}$ は、$S^{\bf k}_A$ によって $N_{\bf q}$ から生成される
メッシュ点の数である。Eq. (I93) では ${\bf k}$ は固定されていることに注意。
$S_A^{\bf k}\subset S_A$ なので、IBZ $\subset$ EIBZ(${\bf k}$) である。すなわち、IBZ で等価な $q$ 点が
EIBZ(${\bf k}$) では区別されることがある。

Ref.[friedrich_efficient_2010] の Eq.(50) を参照。
`call Seteibz -->call eibzgen` は $S_A^{\bf k}$ を返し、module `m_eibz` に格納する。

規格化は次のとおりであることに注意。
$$\begin{aligned}
N=\sum_{\bf q}1 = 
\frac{1}{N_A^{\bf k} } \sum_{A^{\bf k}_i \in S_A^{\bf k} } 
\sum_{ {\bf q}\in {\rm EIBZ}({\bf k})}N_{\bf q}^{\bf k}=
\sum_{ {\bf q}\in {\rm EIBZ}({\bf k})}N_{\bf q}^{\bf k}\end{aligned}
\tag{I94}
$$

[分極関数について]{.underline}\
Eq. (I93) を $F^{\bf k}[g_{\bf q}]$ に適用すると、次を得る。
$$\begin{aligned}
X({\bf k})
=\frac{1}{N_A^{\bf k} } \sum_{A^{\bf k}_i \in S_A^{\bf k} } \left( 
\sum_{ {\bf q}\in {\rm EIBZ}({\bf k})}N_{\bf q}^{\bf k}F^{\bf k}[A^{\bf k}_i[g_{\bf q}]] \right).
\end{aligned}
\tag{I95}
$$
これは、$g_{\bf q}\in {\rm EIBZ}({\bf k})$ についてだけ重み付き平均を計算し、
最後に対称化すればよいことを意味する（対称化は
`hx0fp0.sc.m.F,hx0fp0.m.F` または新たに開発中の `hrcxq` の中の
`call x0kf_v4hz_symmetrize` にある）。

Eq. (I95) を Ref.[friedrich_efficient_2010] の Eq.(22) に適用すると、
Ref.[friedrich_efficient_2010] の Eq.(51) を得る。
そこで、与えられた空間群の回転 $A^{\bf k}_i$ に対して MPB を回転する必要がある。
$A^{\bf k}_i$ は MPB の ${\bf k}$ を保つことに注意。これは次の節で述べる。

`call eibzgen` は $S_A^{\bf q}$（`eibzsym`）だけでなく、
重み $N_{\bf q}^{\bf k}$（`nwgt`）も module `m_eibz` に格納することに注意。
（`m_eibz` では ${\bf q}$ と ${\bf k}$ の役割が入れ替わっているかもしれない。）

 \
[自己エネルギーについて]{.underline}\
`hsfp0.sc.m.F`（`hsfp0.m.F` ではない）では、
Ref.[friedrich_efficient_2010] の Eq.(52) のような対称化を行う。
ただし用いる式は少し異なる。

Eq. (I94) をそのまま用いる。
分極関数と同じように、最初に `eibzgen` を
呼ぶ。${\bf q}$ での自己エネルギーの重み付き和を得るため、
${\bf k}\in {\rm EIBZ({\bf q})}$ について重み付き和を計算する。
（和は `nrkip` で制御されることに注意。`nrkip` は `nwgt` から写した重みである
（MPI 並列化を用いる）。`nrkip=0` のときは、
対応する do ループを飛ばす。）

対称化を行うのは `call zsecsym` で、`hsfp0.sc.m.F` の主プログラムから呼ばれる。
その中心部分は次のとおりである。

           zsect(ii1:ie1,ii2:ie2)= zsect(ii1:ie1,ii2:ie2) 
         &   + matmul( dconjg(transpose(rmatjj(1:ne1,1:ne1,iblk1))), 
         &     matmul(zsec(ii1:ie1,ii2:ie2,iqxx),
         &     rmatjj(1:ne2,1:ne2,iblk2)) )

。これは回転行列 `rmatjj` で回転した自己エネルギーである。
`rmatjj` はその中の `call rotwvigg` で作られる。
（注: 対称性をよく保つため、縮退したバンドを注意深く扱う。
これにより自己エネルギー行列が（縮退による）ブロックに分解される。
ブロックは `zsect(ii1:ie1,ii2:ie2)` で、バンドの添字
`ii1:ie1` は縮退したバンドである。）

## MPB での電子格子結合（中断） {#electron-phonon}

この作業は中断した。フォノンのための静的な $\epsilon^{-1}$ は、`w4pmode` がオンのとき、すなわち $x,y,z$ 方向の offset-Γ 点を `qg4gw --job=201` で加えたとき（`m_q0p.f90` の `ixyz`。`m_llw.f90`）に、`hgw` が `__W4PHONON.<iq>` に書き出す。

電子格子（EP）結合のための微分を計算するには、次が要る。
$$
\begin{aligned}
&\frac{\partial}{\partial {\bf r}_i} \left(\bar{J}_l(E,{\bf r}) Y_L(\hat{r})\right)\Big|_{ {\bf r}=0}=
\frac{1}{3} \sqrt{\frac{3}{4\pi} } \\
&\qquad (\text{only for } l=1;\ \text{zero otherwise}).
\end{aligned}
\tag{I39}
$$
実球面調和関数では、場合は次のとおりである。
$i=y$ と $L=(1,-1)$、$i=z$ と $L=(1,0)$、$i=x$ と $L=(1,1)$。
MPB を用いるので、原子核と
$\langle B\rangle$ の間の Coulomb 相互作用の裸の行列要素を次のように評価する必要がある。
$$
\langle B|\frac{\partial v({\bf r}-{\bf R})}{\partial {\bf R}_i}\rangle
\tag{I40}
$$
