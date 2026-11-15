# `lmf,lmfa,lmchk`

入力ファイルは `ctrlg.<sname>.toml`（古い形式では `ctrl.<sname>`）です。
`<sname>` はプログラムに渡す**位置引数**で（例えば `lmf si` は `ctrlg.si.toml` を読みます）、このページでは、
`ctrlgenToml.py` や `Legacy2toml.py` で決めた物質の名前を `<sname>` と書きます。

## lmfa: 球対称な原子の計算
例:
```
lmfa si
```
数秒で終わります。
```
lmfa si |grep conf
```
で電子配置が分かります。`lmfa` は、単位胞に含まれる原子について、球対称な原子の電子密度を計算します。
その重ね合わせから、`lmf` の電子密度の初期値を作ります。
また `lmfa` は、MT 半径での動径関数の対数微分を `atmpnu.*.<sname>` に書きます（[`[ham].readp`](#readp-option)）。

## lmchk: 結晶の対称性を確かめる
`lmchk` は、`lmf` が認識した空間群の対称性を表示します。MT 半径も決めます。空間群の操作の書き方は [SYMGRP](#symgrp) にあります。

### lmchk で決める MT 半径
理想的な球の半径は、MT の中では球対称、外では平らなポテンシャルをいちばんよく近似するものです。
`lmchk` には、MT 半径の妥当な初期値を決める方法が一つ入っています。
ある原子から比較的近くの原子へのすべての方向について、自由原子の密度を重ねて作った（静電）ポテンシャルを計算し、
どの方向でも最初に現れるポテンシャルの極大までを MT 半径とします。
中のポテンシャルがほぼ球対称になるという意味で、かなり妥当な見積もりです。
完全に対称な結合では、ポテンシャルの極大はちょうど結合の中点に来るので、二つの球はちょうど接し、ポテンシャルも等しくなります。
この半径を求めるには `lmchk --getwsr` とします。`ctrlgenToml.py` はこれを自動で行います。
* `ctrlgenToml.py` は、一部の原子の半径に上限も設けます（`r_upper_limit`）。上限の値とその理由（偽の準位）は [MT 半径](./mtradius) にあります。

* 分子（特に二量体）を扱うときは、とても小さな MT 半径を使う必要があります。そのような系でも PMT がうまく働くことは確かめていますが、もう少し調べる必要があります。

## lmf: 一体問題を解く
LDA/GGA と LDA+U のバンド計算をします。QSGW のポテンシャルの項も `sigm.<sname>` で加えられます。

例:
```
mpirun -np 4 lmf si [options] > llmf
```
これは `ctrlg.si.toml` のある場合です。

### オプション
* 途中で止めるオプション: `--quit=show`、`--quit=ham`、`--quit=mkpot`、`--quit=dmat`、`--quit=band`（[コマンドラインのオプション](./cmdopts) を参照）。
* `--tdos`: 全状態密度の計算。ふつうは `job_tdos` を使います。
* `--ctrlg:<節.キー>=<値>`: 一回の実行だけキーを変えます。[コマンドラインのオプション](./cmdopts#run-time-ctrlg-overrides) を参照。

電子密度の図や BoltzTraP などのためのオプションもあります（まだ説明していません）。

### 入力ファイル

- `syml.<sname>`: バンドプロットの対称ライン（ふつうは `job_band` が中で `lmf` を呼びます）
- `sigm.<sname>`: QSGW の計算の一部
- `atmpnu.{1,2,3}.<sname>`: `lmfa` が作ります（`[ham] readp = true` のとき）

#### `[ham].readp` と `[ham].pnufix` {#readp-option}
`readp = true` のとき、`lmf` は `lmfa` が `atmpnu.*` に書いた、MT 半径での対数微分（価電子の $P$ とローカル軌道の $P_z$）を読みます。
`pnufix = true` のとき、それを反復の間固定します。どちらも `ctrlgenToml.py` が書く既定です。そのため反復の間、基底は変わりません。
詳しくは [GW で用いる波動関数の表現](../theory/gw_basis) の「動径方程式について」。

### 出力ファイル
* `rst.<sname>`
  再開のためのファイルで、電子密度とそれに関わる情報が入っています。
  `ctrlg.<sname>.toml`、`atmpnu*`、`sigm*`、`rst*` があれば、`lmf` の計算を再開できます
  （LDA+U では `dmats*` も使います）。

* `save.<sname>`
  `lmf`、`lmfa`、`lmchk` の最小限の履歴で、`lmf` の自己無撞着の反復ごとの全エネルギーが入ります（[下](#save-sname)）。

* `__mixm.<sname>`
  電子密度の混合のファイルです。入力と出力の密度の組を前の反復の分まで持ち、自己無撞着への収束を速める Anderson や Broyden の混合に使います。
  `lmf` は、前の実行のこのファイルを始めに捨てます（オプション `--keepmixm` で残せます）。格子定数を変えたときのような新しい計算に使わないためです。
  磁性の無い解で終わった実行の履歴が残っていると、新しいスピン分極の計算が、最初の Broyden の一歩で磁性の無い解に引き込まれます（Fe: 2.13 → 0.02 μB）。
  引き継ぐのは `rst.<sname>` だけです。すでにある `rst` の上で `nspin` を変えるときは、先に `rst.<sname>` を消してください。

**主なソース**: [`SRC/subroutines/main_lmf.f90`](https://github.com/tkotani/ecalj/blob/main/SRC/subroutines/main_lmf.f90)

## ctrlg.&lt;sname&gt;.toml — いまの入力の形式（2026-05〜）

`lmf`、`lmfa`、`lmchk` は、2026-05 から **TOML だけ**を読みます。
一つのファイル `ctrlg.<sname>.toml` に、以前 `ctrl.<sname>` にあったものすべてと、以前 `GWinput` にあった GW の設定と積基底の節が入っています。
ファイルの最後には、原子ごとの積基底の表があります。古い作業ディレクトリでは `Legacy2toml.py <sname>` で作ります。

ファイルの最後の `[product_basis]` の表は GW のためのもので、`gwinit` が書き、ふつうは手で直しません。[GW の設定: 積基底](./gwinput#product-basis) を参照。

#### 画面の出力 {#reading-the-output}
画面の出力は、いまは主にデバッグのためのものです。
それでも、いくつかの値は画面の出力から読む必要があります（直していくつもりです）。
収束の様子、バンドのエネルギー、フェルミエネルギー、`sigm` や `rst` を正しく読んだかどうかを、画面の出力で確かめられます。

#### `save.<sname>` と `log.<sname>` {#save-sname}
* `save.<sname>` には、`lmf`、`lmchk`、`lmfa` を起動した履歴が残ります。さらに `lmf` の反復ごとに 1 行、全エネルギーを書きます。
  各行の最初の文字は状態を表します。c: 収束、x: 収束しないまま反復の上限に達した、i: 反復の途中、
  h: 最初の反復（自由原子の密度の重ね合わせからの Harris のエネルギー）、C: 構造の緩和の収束（`m_lmfp.f90` の `svflg='cxhiC'`）。
* 使った `--ctrlg:` の上書きも記録します。

* Kohn-Sham と Harris-Foulkes の二つの全エネルギーが出ます。[二つはほぼ同じになるはず](#q-should-the-harris-foulkes-and-hohenberg-kohn-sham-functionals-agree-at-self-consistency)ですが、
  大きな系ではいくらか違うことがあります。どちらか一つを使ってください。
* QSGW のモード（`sigm` のファイルを使う `lmf`）では、全エネルギーには意味がありません。`gwsc` の中の `lmf` の反復でも全エネルギーは出ますが、収束の目安としてだけ見てください。

* ログのファイル
  `lmf`、`lmchk`、`lmfa` が作る `log.<sname>` は、いまはデバッグのためのものです。

### ファイルの構成（節） {#file-structure-sections}

```toml
# Top-level scalars
symgrp  = "find"        # space-group spec: 'find' = auto-detect
verbose = 35            # console verbosity
time    = [0, 0]        # CPU timing log: [depth, on-the-fly]

[struc]   alat / plat                              (no nbas / nspec: the counts are those of the tables below)
[[site]]  atom / pos / xpos / af / relax           (one table per atom)
[[spec]]  atom / z / r / rsmh / eh / lmx / ...     (one table per species)
[bz]      nkabc / metal / tetra / npts / ...
[iter]    nit / mix / b / conv / convc / umix / tolu
[ham]     nspin / rel / so / phispinsym / xcfun / gmax / pwmode / pwemax / oveps / ...
[esm]     boundary / origin / shiftmode / zb / potential / field   (slabs only)
[gw]      n1n2n3 / QpGcut_psi / HistBin_dw / iSigMode / niw / t_sigmaw / t_tetrakbt / ...
          QforEPS / QforGW (multi-line q lists)
[mlo]     mlo_method / mlo_delta / mlo_w / mlo_nkabc (required for job_mlo) / mlo_lm (formerly Worb)
[blocks]  QPNT, QforEPSL, hrotr (raw multi-line blocks with no better home)
[product_basis]   pb_tolerance / pb_lcutmx + nlx / valence / core   (always the last section)
```

`[struc]` に `nbas` / `nspec` を書かないでください。[入力ファイル](./toml_migration#an-older-ctrlg-sname-toml-ctrlg-update-py) を参照。

### 例: fcc の Cu（1 原子、磁性なし）

`Samples/EPS/EPS_Cu/` の `ctrlg.cu.toml` の全体を下に示します。行末の `#` のコメントが、それぞれのキーの役割です。
fcc の単体の金属の出発点として写し、`[spec]`、`[bz]`、`[ham]` をそれぞれの場合に合わせてください。
このファイルはサンプルの古い入力から変換したものなので、いまの `ctrlgenToml.py` が書く値や既定値と違うものがあります
（例えば、ここでは `pwmode = 1` だが既定は 11。`[iter] mix = "A2"`、`nit = 30`。`[gw] n1n2n3`、`QpGcut_*`、`t_sigmaw = 262`、`t_tetrakbt = 0`）。
いまの雛形の値は、[GW の設定](./gwinput) と、新しく作った `ctrlg.<sname>.toml` のコメントを見てください。

```toml
# Top-level scalars
symgrp  = "find"  # 'find' = auto-detect space group from lattice + species
verbose = 35      # 31 default, 35 verbose, 41+ debug
time    = [0, 0]  # CPU timing log: [depth, on-the-fly]

[struc]
alat  = 6.798     # lattice constant (a.u.)
plat  = [[0.0, 0.5, 0.5],   # primitive vectors in units of alat
         [0.5, 0.0, 0.5],
         [0.5, 0.5, 0.0]]

[[site]]
atom = "Cu"
pos  = [0.0, 0.0, 0.0]

[[spec]]
atom   = "Cu"
z      = 29
r      = 2.33    # MT radius (a.u.)
rsmh   = [1.17, 1.17, 1.17, 1.17]   # smoothing radii per l
eh     = [-1.0, -1.0, -1.0, -1.0]   # tail energies
rsmh2  = [1.17, 1.17, 1.17]         # second basis group
eh2    = [-2.0, -2.0, -2.0]
lmx    = 3       # basis l-cutoff
lmxa   = 4       # augmentation l-cutoff
nmcore = 1       # 1: spin-averaged core; 0: spin-polarized core
kmxa   = 5       # radial expansion order at tail sites

[bz]
nkabc  = [12, 12, 12]   # k-mesh divisions; 1 entry => isotropic
metal  = 3              # 0 insulator, 3 tetra+broadening (safe metal)
savdos = true           # write dos.tot.* file
npts   = 2001
dosmax = 1.5

[iter]
nit   = 30        # max self-consistency iterations
mix   = "A2"      # B3 = Broyden hist=3 (default); A2 = Anderson hist=2
conv  = 1e-05     # max dE between iterations (Ry)
convc = 1e-05     # max d(rho_out - rho_in) (Ry)
b     = 0.2       # mixing ratio

[ham]
nspin       = 1   # 1 nonmag, 2 spin-polarized (then spec.mmom seeds)
rel         = true
so          = 0   # 0 none, 1 L.S, 2 Lz.Sz
gmax        = 12.0   # real-space mesh G-cutoff (Ry); alt: ftmesh
frzwf       = false
xcfun       = 1   # 1=VWN, 2=Barth-Hedin, 103=PBE-GGA
forces      = 0
scaledsigma = 1.0   # QSGW mixing: 1.0 full, 0.8 = QSGW80
oveps       = 0.0
pwmode      = 1   # 0=MTO 1=PMT 2=APW 11=PMT(|q+G|, default) 12=LAPW
pwemax      = 3.0  # APW cutoff (Ry)

# === GW driver (used by lmf --jobgw=N and gwsc) ===
[gw]
n1n2n3        = [10, 10, 10]   # GW BZ mesh (typically smaller than [bz].nkabc)
QpGcut_psi    = 3.0     # |q+G| cutoff for eigenfunctions
QpGcut_cou    = 2.0     # |q+G| cutoff for Coulomb / W
unit_2pioa    = false   # false: a.u.;  true: 2*pi/alat
HistBin_dw    = 1e-05   # bin width on real-axis at omega=0
HistBin_ratio = 1.03
iSigMode      = 3       # self-energy mode (3 = standard QSGW)
niw           = 10      # # of imag-axis frequencies
deltaw        = 0.02
t_sigmaw      = 262     # (K) Fermi-Dirac width of the Sigma levels
t_tetrakbt    = 0       # (K) required. 0: no smearing of chi0; T > 0: finite-T tetrahedron; -T: Gaussian on Im chi0
QforEPS = """           # q points for eps (a.u. if QforEPSau = true)
 0 0 0.00050
 0 0 0.00100
 0 0 0.00200
"""
QforGW = """            # q points for the one-shot GW driver (gw_lmfh)
 0.0 0.0 0.0
 0.1 0.0 0.0
"""

# === MLO (muffin-tin based localized orbitals; see manual/mlo) ===
[mlo]
mlo_method = 4          # theta = sigma((eps - ecut_j)/mlo_w), ecut_j = max(CBM + mlo_delta, eps^MTO_j)
mlo_delta  = 2.0        # (eV) how far above the band edge the model must be accurate
mlo_w      = 2.0        # (eV) width of the fall-off; the knob if the residual is too large
mlo_lm     = """        # which lm channels of which atom make the model (1=s, 2-4=p, 5-9=d, 10-16=f)
  1 Cu   5 6 7 8 9
"""

[product_basis]
pb_tolerance = [0.001]  # drop near-linear-dep products (default 1e-3)
pb_lcutmx    = [4]      # max l-cutoff per atom
# ----------------------------------------------------------------
# nlx [iatom, l, nnvv, nnc] ... (per-atom tables written by gwinit; see manual/gwinput)
nlx = [
  [1, 0, 2, 3],
  [1, 1, 2, 2],
  [1, 2, 2, 0],
  [1, 3, 2, 0],
  [1, 4, 2, 0],
]
# ----------------------------------------------------------------
valence = [
  [1, 0, 1, 1, 1],   # 4s_phi
  [1, 0, 2, 0, 0],   # 4s_phidot
  # ...
]
# ----------------------------------------------------------------
core = [
  [1, 0, 1, 0, 0, 0, 0],   # 1S
  # ...
]
```

（q 点の一覧 `QforEPS` / `QforGW` は、複数行の文字列として `[gw]` の最後に置きます。上の例にあります。）入力はこれで全部です。

生成されるファイルの体裁（空行、`# ----` の罫線、節の順）は [入力ファイル](./toml_migration#file-layout) にあります。

#### 物理ごとの差分（Cu を出発点として）

- **スピン分極**: `[ham].nspin = 2` とし、各 `[[spec]]` に `mmom` を書きます。[その他の説明](./UsageDetailed#spin-polarized-case-without-soc) を参照。きれいなディレクトリから始めてください（上の `__mixm.<sname>`）。
- **SOC**: `[ham].nspin = 2` で、`[ham].so = 1`（L·S）または `2`（Lz·Sz だけ）。[その他の説明](./UsageDetailed#spin-orbit-coupling) と `Samples/MLOsamples/FeSoc/` を参照。
- **絶縁体・半導体**（Si、GaAs）: `[bz].metal = 0`、`[bz].tetra = true`。`[bz].nkabc` は小さめでよく、`[ham].pwemax` は 3 くらいで足ります。
- **LDA+U**: 下の [LDA+U](#lda-u-toml) を参照。
- **原子位置の緩和**: [力と構造の緩和](./relax) を参照。
- **GW の BZ のメッシュ**: `[gw].n1n2n3` ≤ `[bz].nkabc`。[GW の設定](./gwinput) を参照。

### 実行時の `--ctrlg:` による上書き

一回の実行だけキーを変えるには `--ctrlg:<節.キー>=<値>` を付けます（例 `lmf si --ctrlg:bz.nkabc=[8,8,8]`）。ファイルは書き換えません。書き方の細部、ファイルに無いキーの扱い、やめた古い書き方（`-v<NAME>=<VAL>` など。使うと止まります）は [コマンドラインのオプション](./cmdopts#run-time-ctrlg-overrides) にあります。

古い `ctrl.<sname>` のキーと TOML のパスの対応は [入力ファイル](./toml_migration#legacy-toml-map) にあります。
ESM（真空のあるスラブ）の `[esm]` の節は [ESM](./esm) のページにあります。

### lmf の収束

`lmf` 自身が収束しないときは、`[iter] b` を 0.015、0.01、0.005 のように小さくします。`gwsc` での自己エネルギーの混合（`[gw].mixbeta`、`mixpriorit`）は [gwsc](./gwsc#mixbeta-and-mixpriorit) にあります。

### 使っている MTO を確かめる
画面の出力の始めの方に、使っている MTO が表示されます（GaAs の場合）。
 ```
 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-06
 spec      l    rsm    eh     gmax    last term    cutoff
  Ga       0*   1.13  -1.00   6.579    1.19E-06    1459
  Ga       1*   1.13  -1.00   7.028    1.26E-06    1807
  Ga       2*   1.13  -1.00   7.475    1.09E-06    2109
  Ga       3    1.13  -1.00   7.920    1.06E-06    2637
  Ga       0*   1.13  -2.00   6.579    1.19E-06    1459
  Ga       1*   1.13  -2.00   7.028    1.26E-06    1807
  Ga       2    1.13  -2.00   7.475    1.09E-06    2109
  As       0*   1.18  -1.00   6.300    2.13E-06    1243
  As       1*   1.18  -1.00   6.720    1.26E-06    1471
  As       2*   1.18  -1.00   7.140    1.37E-06    1837
  As       3    1.18  -1.00   7.558    1.05E-06    2229
  As       0*   1.18  -2.00   6.300    2.13E-06    1243
  As       1*   1.18  -2.00   6.720    1.26E-06    1471
  As       2    1.18  -2.00   7.140    1.37E-06    1837
 ```

#### ローカル軌道（`pz`）と `ecore>evalence` のメッセージ
コアの電子のはみ出し（spillout）が小さくない（0.05 くらいより大きい）ときは、`pz`（ローカル軌道）を使い、そのコアを価電子として扱う方がよいです。
Bi の 4d がその例で、`pz = [0.0, 0.0, 4.9]` などとします。

メッセージ `ecore>evalence` は、コアの準位が価電子の準位より高いことを意味します。安全のため、これは許していません。
ECORE のファイルと、ログや画面の出力にある価電子の準位を比べてください。

### 分子
`ctrlgenToml.py <sname> --systype=molecule` で、箱の中の分子の雛形を作ります（Γ 点だけ、`[bz] nkabc = [1, 1, 1]`）。
分子では `[bz]` に次を書きます。
```toml
tetra = false
n = -1        # Negative is the Fermi distribution function. w gives temperature.
w = 0.001     # w=0.001 corresponds to T=157K as shown in console
```
磁気モーメントがあるときは、`fsmom`（n_up − n_down）も要ります（`fsmommethod = 1`）。
サンプル: `Samples/AtomDimer/N2`（箱の中の N₂ と N 原子。結合長と結合エネルギー）。
`ctrlgenToml.py --systype=molecule` が書く雛形は `tetra = true`、`w = 0.01` で、上の値（サンプル `Samples/AtomDimer/N2` のもの）と違います。
どちらがよいかはまだ確かめていません（2026-11-14）。
H₂ から Kr₂ までの分子に PMT 法を使った例: T. Kotani and H. Kino, J. Phys. Soc. Jpn. **82**, 124714 (2013),
[doi:10.7566/JPSJ.82.124714](https://doi.org/10.7566/JPSJ.82.124714)
（[pdf](../presentations/KotaniKino2013PMTMolecule.pdf)）。

---

## LDA+U {#lda-u-toml}

LDA+U は、種ごとに、各 `[[spec]]` の表に `idu`、`uh`、`jh` のベクトルを書いて選びます。添字は `[s, p, d, f]` です。

```toml
[[spec]]
atom = "Er"
z    = 68
r    = 3.00
# ... rsmh, eh, ... ...
idu  = [0, 0, 2, 2]          # 0 = none, 1 = AMF, 2 = FLL.  +10 disables when sigm.* exists
uh   = [0.0, 0.0, 0.10, 0.632]    # Hubbard U (Ry)
jh   = [0.0, 0.0, 0.00, 0.055]    # exchange J (Ry)
```

`idu[l]` の値:
- `0` — この l には U を入れない
- `1` — **AMF**（Around-Mean-Field）、Petukhov, PRB 67, 153106 (2003)
- `2` — **FLL**（Fully Localized Limit）、Liechtenstein, PRB 52, R5467 (1995)
- `12` / `21` / ... — 同じ値に 10 を足すと、「`sigm.<sname>` があるときは LDA+U を使わない」の意味になります
  （`sigm` のある QSGW の計算は、ただの QSGW になります）。

`dmats.<sname>`（密度行列）は自動で読み書きします。LDA+U での `rst.<sname>` に当たるものです。U なしの最初の LDA の計算で作られます。
`occnum.<sname>`（スピンごとに `2l+1` 個の数を 1 行）を手で書くと、初期の占有を指定できます。ErAs の例とスピン配置の落とし穴は
[古い LDA+U の覚え書き](./lmf_legacy_ctrl#lda-u) にあります。

密度行列の混合は、電荷の混合とは**別に**、線形の混合で行います。

```toml
[iter]
umix = 0.2          # dmat mixing beta (0..1)
tolu = 0.0          # rms tol where dmat mixing stops (0 = always mix)
```

`lmf <sname> --quit=dmat` で、密度行列を一度確かめてください。初期の設定の後で止まり、再開のファイルを書かずに行列を表示します。

## コアホール {#core-hole-toml}

コアの占有を部分的に変えることで、Slater の遷移状態の考え方を実現します（XAS や、4f を外すときに使えます）。種ごとに設定します。

```toml
[[spec]]
atom    = "Nd"
z       = 60
r       = 3.00
p       = [0.0, 0.0, 0.0, 5.2]    # treat 4f as core, 5f as valence
c-hole  = "4f"                     # core channel to perturb
c-hq    = [-11.0]                  # excess electron charge in that channel
# c-hq    = [-7.0, 7.0]            # optional 2nd entry = spin moment
```

電荷の勘定:
```
Q(spin1) = full_core_count + c-hq[1]/2 + c-hq[2]/2
Q(spin2) = full_core_count + c-hq[1]/2 - c-hq[2]/2
```
`c-hq[1]` が負なら電子を取り除きます（穴）。2 番目の `c-hq[2]` はコアを分極させます
（例えば Gd の 4f の技法: `c-hq = [-7, 7]` で、4f のコアを上向き 7、下向き 0 に固定し、5f は価電子のままにします）。

中性になっていることを `lmfa <sname> | grep -i 'add core hole'` で確かめてください。Gd、Nd、CrN の例は
[古いコアホールの説明](./lmf_legacy_ctrl#core-hole-options) にあります。

## ローカル軌道 {#local-orbitals-toml}

ローカル軌道 `pz` は、連続な主量子数 `p` に対応するものです。どちらも `[[spec]]` の中の、種ごと・`l` ごとのベクトルです。

```toml
[[spec]]
atom = "Ga"
# p[0..3] = s,p,d,f valence principal quantum numbers (continuous);
#          fractional part 0.5 = mid-band, ->0 = bottom, ->1 = top.
p   = [0.0, 0.0, 4.0, 0.0]       # 4d "valence-like"
# pz[0..3] = semicore (lower n) or high-lying (higher n) local-orbital prinq;
#           fractional part ~0.9 = semicore, ~0.3 = high-lying.
pz  = [0.0, 0.0, 3.9, 0.0]       # 3d treated as semicore local orbital
# add 10 to pz[l] to attach a smooth-Hankel tail (extended local orbital).
```

注:
- `pz` と `p` の整数部分が、動径関数の主量子数を決めます（semicore なら `int(pz) = int(p) − 1`、高い準位なら `int(pz) = int(p) + 1`）。
- 小数部分は、semicore で `0.9`、高い準位で `0.3` がふつうです。小数部分が小さいほど、線形化のエネルギーは低くなります。
- `pz[l] + 10`（つまり `pz = [0, 0, 13.9, 0]` と書く）で、スムーズ Hankel の尾を付けます（拡張ローカル軌道）。深い状態が格子間領域にはみ出すときに使います。
- `lmfa <sname> | grep -i conf` で配置を確かめてください。どの状態をコア、価電子、ローカル軌道として扱っているかが分かります。

よくある落とし穴:
- ふつうの（拡張でない）ローカル軌道を semicore に使うと、MT 球の外の表現は価電子の包絡関数に頼ることになります。その l の `rsmh` を少なくとも一つ小さく（1 くらい）、`eh` をかなり負に（−1 くらい）しておいてください。
- `pz` の小数部分が 1.0 に近すぎると、拡張ローカル軌道の接続が失敗することがあります（"`mtchre: failed to match phi'/phi`"）。小さくしてください（例えば `0.98 → 0.85`）。
- 主量子数の大きい方の `p` は自動で固定され、いちばん深い状態だけが動けます（`idmod` で制御します）。

スムーズ Hankel の `mtchre` の昔の議論は、[古いローカル軌道の覚え書き](./lmf_legacy_ctrl#old-note-for-local-orbitals) にあります。

## SYMGRP / symgrp — 空間群の指定 {#symgrp}

TOML では、空間群の指定は最上位の文字列 `symgrp` です。

```toml
symgrp     = "find"           # auto-detect from lattice + species (default)
symgrpaf   = "find"           # AF symmetry (optional): spglib finds it from the af labels; generators only for the old search
```

`symgrp` には、文字列 `"find"`（自動で見つける）か、空白で区切った生成元の並びを書きます。生成元は、古い `SYMGRP` と同じです。

- `R<n>X`、`R<n>Y`、`R<n>Z`、`R<n>D` — X、Y、Z、(1,1,1) のまわりの `n` 回回転
- `M<axis>` — 鏡映
- `I` — 反転
- `<axis> = X | Y | Z | D | (nx,ny,nz)` — 直交座標の軸か、方向を直接指定
- `*` で積をとる。例 `I*R4X`
- 並進: `r6z:(0,0,0.4778973)`（直交座標）、`r6z::(0,0,1/3)`（分率座標、`plat` の倍数）

例（立方晶）:
```toml
symgrp = "R4X MX R3D"        # 4-fold around X, mirror in X, 3-fold around (1,1,1) ⇒ 48 symops
```

`find` と一緒に生成元を書く（古い形で、いまも受け付けます）と、spglib が群全体を見つけ、生成元はその操作でなければいけません（違えば `lmf` は止まります）。生成元は群を変えません。
```toml
symgrp = "R4Z*I MX*I R3D find"   # zinc blende: the same 24 operations as symgrp = "find"
```

対称性をすべて切るには:
```toml
symgrp = "e"
```

`lmchk` は、認識した対称操作を画面の出力の始めに表示します。対称性が意図どおりに読まれたか、これで確かめてください。
ポテンシャルが局所的なとき、DFT の計算では反転は暗に含まれます。画面の出力の冒頭の、反転についての注意を見てください。

### spglib による対称性: `symmetry.<sname>.json`（2026-10） {#spglib-symmetry}

`symgrp = "find"`（既定）のとき、`lmf`、`lmchk`、MLO のプログラムは spglib で空間群を見つけます。spglib は ecalj に組み込んであります（`SRC/external/spglib`、版 2.6.0。Python は要りません）。
操作 $x' = Rx + t$（`plat` の分率座標で、$R$ は整数）を、空間群の番号と記号、それが属する構造と一緒に、作業ディレクトリの `symmetry.<sname>.json` に書きます。
次の実行はこのファイルを読みます。構造が違えば（alat、plat、種、`af` のラベル、位置。`--ctrlg:` の上書きを含む）、操作を見つけ直してファイルを書き直します。
どちらだったかは画面の出力で分かります: `mksym: N operations, read from ...` または `found by spglib, <symbol> (<number>), written to ...`。

* 超格子（基本単位胞より大きい単位胞。例えば Si の 8 原子の立方晶の単位胞）の純粋な並進も操作に含まれます。Si8 では 192 個です。全エネルギー、MLO のバンド、QSGW の準粒子エネルギーは、基本単位胞の操作のときと一致します。
* 反強磁性: 反強磁性の対のサイトに `af = ±k` を書くと、spglib の磁気対称性が操作を与えます。時間反転を伴う操作は上と下のスピンを入れ替えます。`symgrpaf = "find"` で、これらの操作で反強磁性のモードに入ります（2026-10-02）。
  `symgrpaf` の生成元（古い形。例えば NiO の `"i:( 1 1 1 )"`）も受け付け、それはこの反強磁性の操作の中になければいけません。`pwmode = 11` と一緒に使ってください。
  NiO などでの手順は [その他の説明](./UsageDetailed#antiferro-symmetry-without-soc) にあります。
* 生成元を書いた `symgrp`（わざと対称性を下げるとき。例えば軌道秩序のある 4f の状態での `"r4z"`）では、その生成元が生成する群を前と同じように使い、ファイルは書きません。このときの `symgrpaf` には生成元が要ります（`"find"` は使えません）。
* ecalj の古い探し方（格子の操作と、結晶を不変に保つ操作を探す。`ECALJ_SYMFIND=ecalj`）は 2026-10-02 に外しました。探すのは spglib だけです。上の入力の形はどれもいまも受け付けます。
* `symfind.py <sname>`（Python の spglib）で、`lmf` の外でも同じファイルを作れます。

## Q: 自己無撞着のとき、Harris-Foulkes と Hohenberg-Kohn-Sham の汎関数は一致するべきか？ {#q-should-the-harris-foulkes-and-hohenberg-kohn-sham-functionals-agree-at-self-consistency}
（Mark van Schilfgaarde による）
A: はい。一致しないなら、自己無撞着に達していないか、行列要素の評価の許容誤差が緩すぎます。よくある原因は、格子間領域の積分のメッシュが細かくないことです。
次は、格子間の密度に 14×14×14 の分割のメッシュを使ったときの出力の一部です。

```
 GVLIST: gmax = 4.702 a.u. created 973 vectors of 2744 (35%)
         mesh has 14 x 14 x 14 divisions; length 0.588, 0.588, 0.668
 SGVSYM: 170 symmetry stars found for 973 reciprocal lattice vectors
...

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol= 1.00E-05
 spec      l    rsm    eh     gmax    last term    cutoff
  Te       0*   1.30  -0.10   5.220    9.30E-05     973*
  Te       1*   1.30  -0.10   5.597    4.36E-04     973*
  Te       2    1.30  -0.10   5.976    2.04E-03     973*
  Te       0*   1.40  -1.00   4.847    2.11E-05     973*
  Te       1    1.40  -1.00   5.182    9.90E-05     973*
  E        0*   1.50  -0.30   4.524    1.06E-05     847
  E        1    1.50  -0.30   4.823    2.01E-05     973*
```

この表は、格子間のポテンシャルの行列要素を作るのに使う平面波展開で、各基底の軌道がどれだけ収束しているかを示します。
指定した許容誤差（$10^{-5}$）より下まで収束しているのは E の s 軌道だけです。Te の d 軌道の収束は特に悪いですが、Te の d は結合にほとんど関わらないので、エネルギーにはあまり影響しません。
この場合（Te 3 原子）、HF と HK の全エネルギーの差は 10⁻⁴ Ry くらいです。分割の数に対する HF のエネルギーの収束と、HK と HF のエネルギーの差を次の表に示します。
```
  n       ehf         ehf-ehk
  14 -14561.140624    0.000089
  15 -14561.133966    0.000017
  16 -14561.132309    0.000005
  20 -14561.131170    0.000000
```

## 古い文書

昔の `Category_Token` の説明の大部分と、物理の古い覚え書き（APW のテスト、スムーズ Hankel の `mtchre` の議論、古い SPEC の覚え書き、古い四面体法のメモ）は、別のページに移しました。

- [古い `ctrl.<sname>` の説明（Category_Token）](./lmf_legacy_ctrl)

古いキーの*意味*はそのページで引き、[古いキーと TOML のパスの対応](./toml_migration#legacy-toml-map) で TOML に直してください。
