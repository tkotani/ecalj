# gwsc: a script to run QSGW calculation

> ⚠️ `gwsc` reads `ctrlg.<sname>.toml` only, and stops when the file is missing. A directory with the legacy `ctrl.<sname>` + `GWinput` is converted by `Legacy2toml.py <sname>` (run once per directory). See [TOML migration](./toml_migration) and the [Samples/TestInstall/*_gwsc/](https://github.com/tkotani/ecalj/tree/main/Samples/TestInstall) directories as templates.

gwscがQSGW計算実行スクリプトである。
QSGW計算は，複数のfortran実行ファイルを呼び出して実行される．
入力は `ctrlg.<sname>.toml` から読み込まれる(GW の設定は `[gw]`。キーの説明は [gwinput](./gwinput.md))。

## Usage
**Usage**: ` gwsc -np NP [-np2 NP2] [--gpu] [--prec=tf32|fp32|fp64] [--mp] [--fp32] nloop extension [Options]`

### `-np NP`
MPI並列数を指定する。

### `-np2 NP2`
GPU版を使用する場合のみ指定する。GPUで実行される計算のMPI並列数を指定する。
通常は使用できるGPU数を指定する。

### `--gpu`
GPU版を使用する場合のみ指定する。

### `--prec=tf32|fp32|fp64` (2026-09)
GW 部分の精度。`fp64` は倍精度版、`fp32` は混合精度版で単精度の積を FP32 で（`--mp --fp32` と同じ）、`tf32` は `fp32` のうち
Σc の最後の積だけ入力の仮数を 10 ビット（TF32 か FP16、GPU ごとの表が選ぶ）にする（LiTi₂O₄ 6³ で E_F ±1 eV の Re Σc は倍精度と 1.3 meV 以内、RTX 5090 の `hgw` は fp32 の約半分の時間: 173 対 367 秒、2026-09-27 夜の版）。旧来の `--mp` だけの指定は `tf32` と同じ
（2026-09-27 まではすべての積が TF32 だった。誤差が 6 倍で、時間は 1 割しか縮まない）。
GPU でどの方法で積や逆行列を計算するかは、インストール時に測った表（`<bindir>/ecalj_linalg_policy.toml`）で自動的に決まる。
詳細は [GPU 版マニュアル § 精度の選び方と行列演算の自動選択](./ecaljgpu#精度の選び方と行列演算の自動選択-2026-09)。
`--mp`・`--fp32` と一緒に書くときは矛盾しないこと（矛盾すると止まる）。

### `--prec-final=<prec>[:N]` (2026-09)
その回の最後の N 反復（省略時 1）だけ別の精度にする。例: `gwsc 10 --prec=tf32 --prec-final=fp32:2 ...` は 8 反復を tf32、最後の 2 反復を fp32 で回す。

### `--no-tetwt-helper` (2026-09)
`--gpu` のとき、四面体の重みを CPU で先に計算する補助の実行（`hgw --tetwt_write`）を使わない。
既定では使う（[GPU 版マニュアル § 四面体の重み](./ecaljgpu#四面体の重みを-cpu-で先に計算する)）。

### `--mp`
混合精度版を使う。2026-09-27 から `--prec=tf32` と同じ意味（Σc の最後の積だけ入力の仮数 10 ビット、ほかの単精度の積は FP32）。
それまでは GPU のすべての単精度の積が TF32 で、悪条件な誘電行列（重元素 + 分子アニオン NO3 / N3 / ClO など）では TF32 の誤差が
$W$ / $\Sigma^c$ で増幅された。精度の選び方は上の `--prec` と [GPU 版マニュアル § 精度の選び方](./ecaljgpu#精度の選び方と行列演算の自動選択-2026-09)。

### `--fp32` (2026-06)
`--mp` と併用し、`--prec=fp32` と同じ意味（Σc の積も FP32）。ecalj_auto の GW1500 バッチ（`gwscconv --gpu --mp --fp32 --conv-tol 0.1`）は
これを既定で使う。TF32 で壊れた事例の再現は [`Samples/mptf32problem/`](https://github.com/tkotani/ecalj/tree/main/Samples/mptf32problem)。


### `nloop`
QSGWのイテレーション数を指定する。

### `extension`
`ctrlg.<sname>.toml` の `<sname>` を指定する。

### `Options`
追加のオプションを指定する。
追加オプションは, 全ての実行ファイルの実行時引数となる。以下追加のオプションのリスト。
またlmfへのTOML override (`--ctrlg:ham.so=1` など、`--ctrlg:<toml-path>=<value>` 形式) もここに書く。

------
#### メモリの使い方と MPI の分け方
遮蔽クーロン相互作用 $W$ は `hgw`（`--gpu` では `hgw_gpu`、`hgw_mp_gpu`）の中でメモリに持ち、ファイルには書かない
（解析のために `__WVR.<iq>`・`__WVI.<iq>` に書き出すのは `--dumpW`）。MPI ランクの分け方（$q$ 点の組と、その中の $\omega$・$k$ の並列）は
使えるメモリの量から自動で決まり、`lgw` の `MPI layout:` に出る。調整はコマンドラインのオプションではなく、`ctrlg.<sname>.toml` の `[gw]` のキーで行う。

- `KeepWV`: 自己エネルギー(相関部分)を計算する際に, $W$ を GPU のメモリに載せておく。GPU 版の既定は `true`、CPU 版は `false`。`false` では $W$ を周波数ごとに主記憶から送る（GPU のメモリは減るが遅くなる）
- `zmel_batch_gb`: 行列要素を一度に作るバッチの大きさ（[gwinput](./gwinput.md)）
- `mpi_worker_exch`、`mpi_worker_corr`: 交換・相関の計算で、1 つの $q$ 点の組に割り当てるランクの数。0（既定）は自動。全ランク数の約数を書く

`--keepwv`、`--nwpara=X` というオプションは無い（書くとプログラムが止まる）。`--nb=X` は受け付けるが、`gwsc` の計算では使われない。

#### `--tetwtk`
指定すると, 分極関数を計算する際に, 結合状態間の四面体重みをメモリ上に保持しない。$k$点が多い計算でメモリ不足になる場合に使用する。

#### `--skipGS`
`lmf --jobgw=1` で使用される。
GW計算ではDFT計算(`lmf`)で得られた波動関数を`lmf`とは異なる基底関数で展開しなおす。再展開後の波動関数についてGram-Schmidt正規直行化をしている。その規格直交化をスキップする場合に指定する。

通常は指定する必要はないが、`lmf --jobgw=1`計算が遅い場合には指定することによって計算の高速化が期待される。
* ecaljでは有限のqで誘電関数が計算できるーこのとき分母分子のキャンセレーションが起こるため、波動関数の直行性が正確である必要があり、そのときには--skipGSを使うべきでない。

#### `--normcheck`
`lmf --jobgw=1` で使用される。 GW計算で使用される波動関数の規格直交性を確認したいときに使用する。
normchk.fobar は
```
> head -20 normchk.si
#       IPW         IPW(diag)   Onsite(tot)   Onsite(phi)      Total
      0.436015      0.805123      0.563972      0.562573      0.999988
      0.339134      0.620353      0.660515      0.656881      0.999649
      0.339133      0.620353      0.660516      0.656882      0.999649
      0.339133      0.620353      0.660516      0.656882      0.999649
      0.507738      0.648515      0.492040      0.487673      0.999778
...
```
などとなる。右端の値が、1になっているべきであるが、展開し直しているため従来では高いエネルギー(下の方。ここでは見えてない）でかなり小さくなっていた(0.8などのおおきさ)。最近デフォルトでは、Gram-Schmidt正規直行化をかけているので、正規直行化は8桁程度以上は保たれている。

The first line (corresponding to 1st band of 1st q point) means that total normalization almost unity = 0.999988 = 0.436015 + 0.563972. 

### `--ntqxx`
This fix the number of bands to calculate self energy at the first iteration for each $\bf q$ point in the IBZ.
In principle, the number is determined by

## Cautions
* QPU.[number]runをチェックして、number回のQSGWイテレーションが終了している、と認識する。
(初期状態から実行したいときはすべての`*run*`ディレクトリ、ファイルを消すこと）。

* QSGW.[number]runディレクトリには、QSGWのnumber回目の結果 `rst`, `sigm` (加えて `atmpnu`, `ctrlg.<sname>.toml`, `QPU`, `efermi.lmf`, その反復のログ) が格納されており、これを用いてバンドプロットなどができる。

* 途中経過とメモリの使用量は、各段のログ（`lgw`、`lsxC`、`lvcc` など）で見る。GW のプログラムはランク 0 の画面出力をそこに書き、
  ほかのランクの出力は捨てる（2026-09-27 から）。ランクごとの出力が要るときは `--fullstdo`（`STDOUT/stdout.<rank>.<prog>` になる）。
  どのランクのエラーも、ランク番号付きで `gwsc` の標準エラーに出る。



---

# Other scripts 

`cleargw`: clean up temporary files

`gwscconv`: `gwsc` を自動収束まで回すラッパー。`--conv-tol <eV>` (例 `0.1`) で
QSGW 反復の停止条件を指定する。**収束判定 (2026-06 改定)**: 直近 3 イテレーションの
ギャップ振れが連続して 2 回 `conv-tol` 以下になったら停止 (一発の偶然収束で止まる
のを防ぐ)。**金属 (2026-09-30)**: `lmf` がギャップを出さない反復は、固有値の変化で判定する。
`QPU.<n>run` と `QPU.<n-1>run` の `eQP` (E_F 基準) のうち |e| < 5 eV の状態について、変化の最大が
`--conv-qp <eV>` (既定 0.03) 未満の反復が 2 回続いたら停止。以前のようにギャップが無いところで止める (exit 3) には
`--no-metal`。`--gpu --mp --fp32` を渡せばそのまま `gwsc` に流される。
ecalj_auto の GW1500 バッチもこの script を呼んでいる (`MD/auto.md`)。

`qsgw_status` (2026-06): 走行中の `gwsc` (または gwscconv) 反復の進行
ダッシュボード。実行ディレクトリで以下のように使う:
```bash
qsgw_status               # 一度だけ表示
qsgw_status --last 21     # 最終 iter 番号を渡すと完了 ETA も計算
watch -n 30 qsgw_status   # 常駐モニタ
```
`QPU.*run` / `sigm.iter*` の mtime から反復ごとの経過時間・ペース(中央値)・
ETA、`gwsc` ログ末尾から現在の step (`hgw` 等) と step 内経過時間、
`save.<target>` の sev 差分と `lqpe` の `rmsdel` (σ 収束指標) を表示。
対象名・ログ・反復アンカーは自動検出。run の妨げにならないよう成果物の
読取のみで重い処理はしない。

`gw_lmfh`: The one-shot \GW calculation. Lifetime(impact ionization rate) of QPs.

`job_eps`: dielectric function. No local field corrections ([optical](./optical))
 
(`job_eps --lcf` : Dielectric function with local-field corrections. computationally expensive.)


<!-- \item
{\bf epsPP\_lmfh\_chipm} : non-interacting spin susceptibility. 
One-degree of freedom like Rigid moment approx.
After it ends, you need to do \verb#calj_nlfc_metal# and/or \verb#calj_summary_mat#
to get the full spin susceptibility. -->

`job_mloW` : the MLO model and the matrix elements of v, W and the cRPA W (`--crpa`) in it ([MLO](./mlo) section 6). (`genMLWFx`, the Wannier functions, was removed 2026-10-02.)


# Files used in gwsc
Temporary files are with `__*`. Thus we can delete __* (or use `cleargw`) after you finish `gwsc/job_eps` and so on.

To repeat a small test for gwsc:

```bash
cd ~/ecalj/Samples/TestInstall
testecalj si_gwsc -np 8
```

This is one of the install tests. `testecalj` creates
`si_gwsc_work/` next to `si_gwsc/`, copies inputs in, and runs there.
After that, `ls -rlt si_gwsc_work/` roughly shows which step generates
which files.

The console output of `gwsc 1 si -np 2` for the same input is
```text
===== Ititial band structure ====== 
--> No sigm. LDA caculation for eigenfunctions 
00:00:00.007   mpirun -np 1 /home/takao/bin/lmfa si stdout='llmfa'  Elap. 1.5s
00:00:01.481   mpirun -np 2 /home/takao/bin/lmf si --ctrlg:iter.b=0.5 stdout='llmf_lda'  Elap. 6.1s
lmf successful with b=0.5
Removing __mixm.si
00:00:07.563   mpirun -np 2 /home/takao/bin/lmf si --jobgw=0 stdout='llmfgw00'  Elap. 1.8s
00:00:09.335   mpirun -np 1 /home/takao/bin/qg4gw si --job=1 stdout='lqg4gw'  Elap. 1.0s
===== QSGW iteration start iter 1 === (GW precision fp64)
00:00:10.381   mpirun -np 2 /home/takao/bin/lmf si --jobgw=1 stdout='llmfgw01'  Elap. 3.7s
00:00:14.042   mpirun -np 1 /home/takao/bin/heftet si --job=1 stdout='leftet'  Elap. 1.1s
00:00:15.193   mpirun -np 1 /home/takao/bin/hbasfp0 si --job=3 stdout='lbasC'  Elap. 1.4s
00:00:16.570   mpirun -np 2 /home/takao/bin/hvccfp0 si --job=3 stdout='lvccC'  Elap. 1.4s
00:00:17.966   mpirun -np 2 /home/takao/bin/hsfp0_sc si --job=3 stdout='lsxC'  Elap. 1.9s
00:00:19.862   mpirun -np 1 /home/takao/bin/hbasfp0 si --job=0 stdout='lbas'  Elap. 2.1s
00:00:21.973   mpirun -np 2 /home/takao/bin/hvccfp0 si --job=0 stdout='lvcc'  Elap. 2.5s
00:00:24.472   mpirun -np 2 /home/takao/bin/hgw si --jobgw=1 stdout='lgw'  Elap. 4.8s
00:00:29.276   mpirun -np 1 /home/takao/bin/hqpe_sc si stdout='lqpe'  Elap. 1.5s
Using cached b-value 0.5 for /home/takao/bin/lmf
00:00:30.814   mpirun -np 2 /home/takao/bin/lmf si --ctrlg:iter.b=0.5 stdout='llmf'  Elap. 8.0s
lmf successful with b=0.5
===== QSGW iteration end   iter 1 ===
```
`lmf --jobgw=0` and `qg4gw` are run once, before the iterations; the other steps are repeated in every iteration.

## `lmfa`
### atmpnu*
Logalithmic derivative at MT boundaries generated by `lmfa`. This is
used by `lmf` when `[ham].readp = true` in `ctrlg.<sname>.toml`.

### `__atm.<sname>`
This contains electron densities of spherical atoms.

## `lmf`
### `rst.<sname>`
This contains self-consistent electron density revised by each iteration of lmf

### QPLIST*.chk
QPLIST.jobgw1.chk (jobgw1 means --jobgw=1 for lmf) is for human, containing q for irr=1.
QPLIST.lmf.chk  (no jobgw option) is for human. Irreducible q points for lmf self-consistent calculations.

### __HAMindex
q points table and so on for generating Hamiltonian

### `__mixm.<sname>`
mixing file for lda

## `lmf --jobgw=0`
###  NLAindx.chk
This is for human. It shows index to expand eigenfunctions in MTs.

### __HAMindex0
Generated by `m_hamindex0_init` called in main_lmf.f90.
Index of MTOs, space-group symmetries and so on.

### QBZ.chk
q points mesh. Just for human reading.

## `qg4gw`

### __QGpsi, __QGcou
q+G of the interstitial plane wave (IPW). Type lqg4gw, which shows
```
--- Max number of G for psi, G for Cou= 116 36
 iq=       1 q=  0.000000  0.000000  0.000000 ngp ngc=    111    29 irr.= 1 <--R
 iq=       2 q= -0.250000 -0.250000  0.750000 ngp ngc=    106    34 irr.= 1 <--R
 iq=       3 q= -0.250000  0.750000 -0.250000 ngp ngc=    106    34 irr.= 0 <--R
 iq=       4 q= -0.500000  0.500000  0.500000 ngp ngc=    116    36 irr.= 1 <--R
 iq=       5 q=  0.750000 -0.250000 -0.250000 ngp ngc=    106    34 irr.= 0 <--R
 iq=       6 q=  0.500000 -0.500000  0.500000 ngp ngc=    116    36 irr.= 0 <--R
 iq=       7 q=  0.500000  0.500000 -0.500000 ngp ngc=    116    36 irr.= 0 <--R
 iq=       8 q=  0.250000  0.250000  0.250000 ngp ngc=    116    36 irr.= 1 <--R
 iq=       9 q= -0.012500 -0.012500  0.037500 ngp ngc=    111    29 irr.= 1 <--Q0P
 iq=      10 q= -0.262500 -0.262500  0.787500 ngp ngc=    106    34 irr.= 1 <--Q0P+R
 iq=      11 q= -0.262500  0.737500 -0.212500 ngp ngc=    106    34 irr.= 1 <--Q0P+R
 iq=      12 q= -0.512500  0.487500  0.537500 ngp ngc=    116    36 irr.= 1 <--Q0P+R
 iq=      13 q=  0.737500 -0.262500 -0.212500 ngp ngc=    106    34 irr.= 0 <--Q0P+R
 iq=      14 q=  0.487500 -0.512500  0.537500 ngp ngc=    116    36 irr.= 0 <--Q0P+R
 iq=      15 q=  0.487500  0.487500 -0.462500 ngp ngc=    116    36 irr.= 1 <--Q0P+R
 iq=      16 q=  0.237500  0.237500  0.287500 ngp ngc=    116    36 irr.= 1 <--Q0P+R
 iq=      17 q= -0.012500  0.012500  0.012500 ngp ngc=    111    29 irr.= 1 <--Q0P
 iq=      18 q= -0.262500 -0.237500  0.762500 ngp ngc=    106    34 irr.= 1 <--Q0P+R
 iq=      19 q= -0.262500  0.762500 -0.237500 ngp ngc=    106    34 irr.= 0 <--Q0P+R
 iq=      20 q= -0.512500  0.512500  0.512500 ngp ngc=    116    36 irr.= 1 <--Q0P+R
 iq=      21 q=  0.737500 -0.237500 -0.237500 ngp ngc=    106    34 irr.= 1 <--Q0P+R
 iq=      22 q=  0.487500 -0.487500  0.512500 ngp ngc=    116    36 irr.= 1 <--Q0P+R
 iq=      23 q=  0.487500  0.512500 -0.487500 ngp ngc=    116    36 irr.= 0 <--Q0P+R
 iq=      24 q=  0.237500  0.262500  0.262500 ngp ngc=    116    36 irr.= 1 <--Q0P+R
  OK! End of qg4gw 
```
Here, we have regular mesh points specified by `<--R`. Q0P is the offset Gamma points shown in Q0P. 
irr=1 shows the irreducible q points at which we have to generate eigenfunctions.
ngp is the number of IPW for the expansion of eigenfunctions (controlled by QpGcut_phi). 
ngc is for IPW for the MPB (controlled by QpGcut_cou).

### __BZDATA 
__BZDATA contains info on regular mesh points, and offset Gamma (Q0P). Data for tetrahedron division. 

## `lmf --jobgw=1`
Generate all the following data to perform GW.
See at subroutines/sugw.f90.

### @MNLA_core.chk 
human readable: core index

### @MNLA_CPHI 
human readable (but program use this):
Eigenfunctions expanded within MT

### hbe.d.chk 
human readable: size file for check

###  __PHIVC
radial functions.

### __MTOindex
MTO index

### __VxcEvec*
 Coefficients of eigenfunctions and  eigenvalues for $\langle F_i|H^0|\psi F_j$ in the basis of PMT$\{F_i\}$.

### __GEIG,__CPHI,__EValue
__GEIG: Coefficients of eigenfunctions. IPW part
__CPHI: Coefficients of eigenfunctions. MT  part
__EValue: eigenvalue

### __PPOvlp, __PPOvlpG, __PPOvlpGG
overlap matrix of IPW.

### __VXCFP
XC term in LDA (only diagonal part). A part of vxcevec
Used at hsfp0 but not essential (only for convenience of presentantion).

## `heftet --job=1`

### EFERMI
The Fermi energy in the tetrahedron method

### EFERMI_kbt
Written when `[gw] t_tetrakbt` > 0: the Fermi energy at that temperature, used for $\chi_0$ and $\Sigma$ ([kBT](./kBT) §2.3).


## hbasfp0 (we have --job=3 for core and --job=0 for valence)
### __BASFP* 
 Product basis functions
### __PPBRD*
 Radial integrals on each MT, symbolically written as  $\int \phi(r) \phi(r) B(r) dr$

## `hvccfp0` 
We call two times one for core, and the other for valence.

### __Vcoud*
the Coulomb matrix (eigenvalues and eigenfunctions of the Coulomb matrix in the expansion of MPB)

## `hgw`
`hgw --jobgw=1` calculates the valence exchange, the screened Coulomb interaction $W$ and the valence correlation in one run.
$W-v$ in the expansion of mixed product basis (along the real axis and along the imag axis) is kept in memory and handed to the
calculation of the correlation $q$ by $q$; it is not written to files.
Note that we use a technique to define $W-v$ at ${\bf q}=0$ as an average of the Gamma cell.

### __WV.d
Size of the dielectric function

### __WVR.*, __WVI.*
Only with `gwsc --dumpW` (for analysis): $W-v$ along the real axis and along the imag axis, a file for each $q$.

### freq_r
human readable: energy mesh to accumrate imaginary parts of W-v.


## `hsfp0_sc --job=3` (core exchange) and `hgw` (valence exchange and valence correlation)
These are moved to SEBK at the end of gwsc iteration cycle.

### SEXcoreU,SEXcore2U :  
core exchange by `hsfp0_sc --job=3` . SEXcoreU is diagonal part for human check but not used so much. 
We have *D for down spin (isp=2)as well.
### SEXU,SEX2U 
valence exchange by `hgw` .       SEXU is diagonal part for human check but not used so much. 
### SECU,SEC2U:  
valence correlation by `hgw` .     SECU is diagonal part for human check but not used so much. 
### XCU
LDA exchange correlation
## `lqpe`

### QPU, QPD
human readable format. 
decomposion of self-energy for diagonal elements.QPD is for isp=2.

An example of one-shot GW by `gw_lmfh si` (small size calculation in TestInstall) is:
```
...
           q               state  SEx   SExcore SEc    vxc    dSE  dSEnoZ  eLDA    eQP  eQPnoZ   eHF  Z    FWHM=2Z*Simg  ReS(elda)
  0.00000  0.00000  0.00000  1  -16.91  -1.85   6.62 -12.47   0.22   0.33 -12.24 -12.03 -11.92 -18.54 0.66   1.25708    -12.14031
  0.00000  0.00000  0.00000  2  -13.87  -1.96   2.81 -13.61   0.47   0.59  -0.31   0.16   0.28  -2.53 0.80   0.00000    -13.02308
  0.00000  0.00000  0.00000  3  -13.87  -1.96   2.81 -13.61   0.47   0.59  -0.31   0.16   0.28  -2.53 0.80   0.00000    -13.02308
  0.00000  0.00000  0.00000  4  -13.87  -1.96   2.81 -13.61   0.47   0.59  -0.31   0.16   0.28  -2.53 0.80   0.00000    -13.02308
  0.00000  0.00000  0.00000  5   -4.60  -1.41  -4.27 -11.81   1.19   1.52   2.23   3.42   3.75   8.03 0.78  -0.02515    -10.28546
  0.00000  0.00000  0.00000  6   -4.60  -1.41  -4.27 -11.81   1.19   1.52   2.23   3.42   3.75   8.03 0.78  -0.02515    -10.28546
  0.00000  0.00000  0.00000  7   -4.60  -1.41  -4.27 -11.81   1.19   1.52   2.23   3.42   3.75   8.03 0.78  -0.02515    -10.28546
  0.00000  0.00000  0.00000  8   -5.11  -3.79  -5.14 -15.23   0.91   1.20   2.95   3.86   4.14   9.28 0.76  -0.07397    -14.03341

  0.50000  0.00000  0.00000  1  -16.80  -1.91   5.97 -12.64  -0.06  -0.10 -11.18 -11.24 -11.27 -17.25 0.67   0.84187    -12.73584
  0.50000  0.00000  0.00000  2  -13.77  -2.37   2.91 -13.84   0.47   0.60  -3.84  -3.37  -3.24  -6.14 0.78   0.04107    -13.24015
  0.50000  0.00000  0.00000  3  -13.57  -1.74   3.01 -12.76   0.37   0.47  -2.20  -1.84  -1.74  -4.75 0.78   0.00000    -12.29249
  0.50000  0.00000  0.00000  4  -13.57  -1.74   3.01 -12.76   0.37   0.47  -2.20  -1.84  -1.74  -4.75 0.78   0.00000    -12.29249
  0.50000  0.00000  0.00000  5   -4.27  -1.19  -4.08 -10.97   1.17   1.43   0.76   1.92   2.19   6.27 0.81  -0.00000     -9.53378
  0.50000  0.00000  0.00000  6   -3.83  -1.10  -4.38 -10.98   1.33   1.68   2.98   4.31   4.65   9.03 0.79  -0.01203     -9.30307
  0.50000  0.00000  0.00000  7   -4.73  -1.66  -4.56 -12.65   1.34   1.70   5.45   6.79   7.14  11.70 0.79  -0.08743    -10.95431
  0.50000  0.00000  0.00000  8   -4.73  -1.66  -4.56 -12.65   1.34   1.70   5.45   6.79   7.14  11.70 0.79  -0.08743    -10.95431

  1.00000  0.00000  0.00000  1  -15.60  -2.13   4.43 -13.20  -0.08  -0.11  -8.13  -8.21  -8.24 -12.67 0.78   0.61936    -13.30329
  1.00000  0.00000  0.00000  2  -15.60  -2.13   4.43 -13.20  -0.08  -0.11  -8.13  -8.21  -8.24 -12.67 0.78   0.61936    -13.30329
  1.00000  0.00000  0.00000  3  -13.66  -1.70   3.17 -12.58   0.30   0.39  -3.16  -2.86  -2.77  -5.94 0.77   0.07961    -12.19216
  1.00000  0.00000  0.00000  4  -13.66  -1.70   3.17 -12.58   0.30   0.39  -3.16  -2.86  -2.77  -5.94 0.77   0.07961    -12.19216
  1.00000  0.00000  0.00000  5   -3.97  -0.91  -3.95 -10.33   1.22   1.50   0.31   1.53   1.81   5.76 0.81  -0.00000     -8.82976
  1.00000  0.00000  0.00000  6   -3.97  -0.91  -3.95 -10.33   1.22   1.50   0.31   1.53   1.81   5.76 0.81  -0.00000     -8.82976
  1.00000  0.00000  0.00000  7   -3.59  -2.37  -5.91 -13.53   1.20   1.66   9.81  11.01  11.47  17.37 0.72  -0.40499    -11.86796
  1.00000  0.00000  0.00000  8   -3.59  -2.37  -5.91 -13.53   1.20   1.66   9.81  11.01  11.47  17.37 0.72  -0.40499    -11.86796
```
All of the unit of energy is in eV. Detailed value of  eLDA} is in {TOTE.UP}.
For insulators, the Fermi energy is at the middle of band. For metals, one shot GW can be problematic
if we consider the self-consistency of the Fermi energy.

q  : q vector

state:  Band index n for valence.

SEx: $= \langle\Psi_{{\bf k}n}|\Sigma_{\rm x}^{\rm valence}({\bf r},{\bf r}^{\prime})|\Psi_{{\bf k}n}\rangle$

SExcore: $= \langle\Psi_{{\bf k}n}|\Sigma_{\rm x}^{\rm core}({\bf r},{\bf r}^{\prime})|\Psi_{{\bf k}n}\rangle$

SEc: $= \langle\Psi_{{\bf k}n}|\Sigma_{\rm c}^{\rm valence}({\bf r},{\bf r}^{\prime}, \epsilon_n({\bf k})|\Psi_{{\bf k}n}\rangle$

vxc: LDA exchange correlation energy.$\langle \Psi_{{\bf k}n}|V_{\rm xc}^{\rm LDA}([n_{\rm total}],{\bf r})|\Psi_{{\bf k}n}\rangle$

dSE: $=Z_{n{\bf k}}\times$ dSEnoZ

dSEnoZ: $=\langle\Psi_{{\bf k}n}|\Sigma_{\rm x}^{\rm core}({\bf r},{\bf r}^{\prime})+\Sigma_{\rm xc}^{\rm valence}({\bf r},{\bf r}^{\prime},\epsilon_n({\bf k}))|\Psi_{{\bf k}n}\rangle - \langle\Psi_{{\bf k}n}|V_{\rm xc}^{\rm LDA}([n_{\rm total}],{\bf r})|\Psi_{{\bf k}n}\rangle$= SEx + SExcore + SEc - vxc

eLDA: LDA eigenvalues. $\epsilon_n({\bf k})$

eQP: QP energy.  $\epsilon_n({\bf k})$+dSE

eQPnoZ: QP energy without $Z$. $\epsilon_n({\bf k})$+dSEnoZ

eHF: HF energy of 1st iteration. $\epsilon_n({\bf k})$+SEx + SExcore -vxc
Z: Z factor. $Z_{n{\bf k}}$

FWHM=2Z*S_{img}: Quasi-particle life time. $2 Z_{n{\bf k}} \times {\rm Im}\langle\Psi_{{\bf k}n}|\Sigma_{\rm c}^{\rm valence}({\bf r},{\bf r}^{\prime},\epsilon_n({\bf k}))|\Psi_{{\bf k}n}\rangle$  

ReS(elda):${\rm Re} \langle\Psi_{{\bf k}n}|\Sigma_{\rm x}^{\rm core}({\bf r},{\bf r}^{\prime})+\Sigma_{\rm xc}^{\rm valence}({\bf r},{\bf r}^{\prime},\epsilon_n({\bf k}))|\Psi_{{\bf k}n}\rangle$ 

* NOTE: QPU for `gwsc` is a little different. No Z and no life time shown.
  Shown eQP is just the eigenvalues of starting point of lmf. 

### TOTE.UP
numerical detailed values of QPU. TOTE.DN for QPD

In one-shot GW `gw_lmfh`, TOTE.UP contains LDA and QP energies. 
It contains two kind of QP energies {\tt QP QPnoZ}.

### __mixsig
mixing file for `sigm.<sname>`

### `sigm.<sname>`
self-energy file in the expansion of eigenfunctions of $H^0$.


# Product basis

The product basis (the union of MT-local products and interstitial
plane waves) is the heart of how ecalj represents the screened Coulomb
interaction `W`. It is specified in the **`[product_basis]` section of
`ctrlg.<sname>.toml`**, always the last section of the file:

| key | what |
|---|---|
| `pb_tolerance`, `pb_lcutmx` | global cut-offs — the part you may tune |
| `nlx`, `valence`, `core` | per-atom radial-function tables, written by `gwinit` (via `ctrlgenToml.py` / `Legacy2toml.py`) from `[[spec]]`; almost never hand-edited |

(Until 2026-09 the three tables were a separate `PB.<sname>.toml`. The
binaries do not read such a file: they abort and point at
`ctrlg_absorb.py <sname>`, which appends its tables to `[product_basis]`
and keeps the original as `.bk`.)

> Product basis is originally given by F.Aryasetiawan and O.Gunnarsson
> ([PRB 49, 16214](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.49.16214)).
> The mixed product basis
> ([Solid State Commun. 2002](https://doi.org/10.1016/S0038-1098(02)00028-5))
> is the successor of original.

## `[product_basis]` cut-offs (in `ctrlg.<sname>.toml`)

```toml
[product_basis]
pb_tolerance = [0.001]   # cut off of linear dependency (default 1e-3)
                         #   larger ⇒ fewer product basis. 1d-2 is OK
                         #   when you want a cheaper calculation.
pb_lcutmx    = [4, 4]    # per-site l-cutoff for the product basis.
                         #   4 needed when valence 3d exists (Ni, Ga, ...);
                         #   2 OK for oxygen; 6 needed for 4f / 5f atoms.
                         #   site order = [[site]] order; see lmchk SiteInfo.
```

After editing, re-run `lmf --jobgw=1` (or the next `gwsc` step) — the
`hbasfp0` driver consumes these values and reports the resulting basis
sizes in `lbas` / `lbasC`.

## The per-atom tables (`nlx`, `valence`, `core`)

They follow the cut-offs in the same `[product_basis]` table; every row
starts with the site index `iatom` (in `[[site]]` order). `gwinit` writes
them from `rsmh`/`eh`/`pz` of `[[spec]]` and adds a label per row; the
explanation below is written into the file itself as comments, and the
three tables are divided by `# ----` rules.

```toml
[product_basis]
pb_tolerance = [1e-3]
pb_lcutmx    = [4]
# ----------------------------------------------------------------
# nlx [iatom, l, nnvv, nnc]: one row per (iatom, l) ...
nlx = [
  [1, 0, 2, 3],
  [1, 1, 2, 2],
  [1, 2, 2, 0],
  [1, 3, 2, 0],
  [1, 4, 2, 0],
]

# ----------------------------------------------------------------
# valence [iatom, l, n, occ, unocc]: one row per valence radial function ...
valence = [
  [1, 0, 1, 1, 1],   # 4s_phi
  [1, 0, 2, 0, 0],   # 4s_phidot
  [1, 1, 1, 1, 1],   # 4p_phi
  [1, 1, 2, 0, 0],   # 4p_phidot
  [1, 2, 1, 1, 1],   # 3d_phi
  [1, 2, 2, 0, 0],   # 3d_phidot
  [1, 3, 1, 0, 1],   # 4f_phi
  [1, 3, 2, 0, 0],   # 4f_phidot
  [1, 4, 1, 0, 0],   # 5g_phi
  [1, 4, 2, 0, 0],   # 5g_phidot
]

# ----------------------------------------------------------------
# core [iatom, l, n, occ, unocc, forX0, forSxc]: one row per core level ...
core = [
  [1, 0, 1, 0, 0, 0, 0],   # 1S
  [1, 0, 2, 0, 0, 0, 0],   # 2S
  [1, 0, 3, 0, 0, 0, 0],   # 3S
  [1, 1, 1, 0, 0, 0, 0],   # 2P
  [1, 1, 2, 0, 0, 0, 0],   # 3P
]
```

Semantics (same as the legacy `<PRODUCT_BASIS>` block — only the
container changed):

- `nlx` row `[iatom, l, nnvv, nnc]` — how many radial functions ecalj has
  on each `l` channel of that site: `nnvv = 2` is `{phi, phidot}`,
  `nnvv = 3` adds the local orbital `phiz`; `nnc` is the number of core
  levels at that `l`. Derived, do not edit.
- `valence` row `[iatom, l, n, occ, unocc]` — one per radial function,
  `n = 1` phi, `2` phidot, `3` phiz. `occ` / `unocc` are 0/1 switches that
  put the function into the "occupied" / "unoccupied" group; the product
  basis is formed from products (occupied group) × (unoccupied group) and
  pruned by `pb_tolerance`. `phidot` is usually left out (`0 0`) for speed;
  including it gives a more accurate basis.
- `core` row `[iatom, l, n, occ, unocc, forX0, forSxc]` — one per core
  level (`n = 1` deepest). `forX0` / `forSxc` came from the CORE1/CORE2
  distinction of PRB 76, 165106 (2007). Historical: modern runs leave
  every core row all-zero.

The site index `iatom` follows the order of `[[site]]` in `ctrlg.<sname>.toml`;
`lmchk` prints the order at the top of its console output.

> **Legacy input.** `Legacy2toml.py` translates the pre-TOML
> `<PRODUCT_BASIS>` block of `GWinput` into the `[product_basis]` section
> above.


# MEMO
* We need to explain how to set Gamma-cell averaged $\tilde{W}({\bf q}=0,\omega)$.
* 旧い開発者向けの理論メモ ecaljdetails（2015〜2022、offset-Γ 法・クーロン行列・EIBZ など）は 2026-10-02 に外した。今も通じる要点は ecalj の `MD/past_log.md` §14、原文は ecaljdoc のコミット `3284d0b` の `ecaljdetails/ecaljdetails.tex`。
