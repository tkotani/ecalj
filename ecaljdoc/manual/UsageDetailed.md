# その他の説明 {#usagedetailed}

下に名前の出るサンプルは [サンプル](./samples) にあります。

## SOC なしのスピン分極 {#spin-polarized-case-without-soc}
磁性のある系を扱うには、`ctrlg.<sname>.toml` で `[ham].nspin = 2` とし、さらに各 `[[spec]]` に `mmom`（スピン磁気モーメントの初期値）を書きます。
スピン磁気モーメントは、スピン `isp=1` と `isp=2` の電子の数の差で指定します。

例えば `ecalj/Samples/TestInstall/nio_gwsc/` の `ctrlg.nio.toml` を新しいディレクトリに写し、
```
lmfa nio
mpirun -np 8 lmf nio
```
を走らせます。`cat save.nio` とすると、`lmf` の反復の履歴が見えます。`save.nio` の最後の行は
`c mmom= 0.0000 ehf(eV)=-86704.349038 ehk(eV)=-86704.348942 sev(eV)=-142.498970` のようになります（エネルギーは設定によります）。`c ` は収束を意味します。
ehf と ehk は [二つの式による全エネルギー](./lmf#q-should-the-harris-foulkes-and-hohenberg-kohn-sham-functionals-agree-at-self-consistency) です。
QSGW では全エネルギーに意味がないことに注意してください（[`save.<sname>`](./lmf#save-sname) を参照）。

`ctrlg.nio.toml` を見ると、次のようになっています。
```toml
[[site]]
atom = "Niup"                 # ← Niup site
pos  = [0.0, 0.0, 0.0]
[[site]]
atom = "Nidn"                 # ← Nidn site
pos  = [1.0, 1.0, 1.0]
[[site]]
atom = "O"
pos  = [0.5, 0.5, 0.5]
[[site]]
atom = "O"
pos  = [1.5, 1.5, 1.5]

[[spec]]
atom   = "Niup"
z      = 28
r      = 2.12
mmom   = [0.0, 0.0, 1.2, 0.0]      # ← initial spin
# rsmh, eh, rsmh2, eh2, lmx, lmxa, nmcore, kmxa follow
[[spec]]
atom   = "Nidn"
z      = 28
r      = 2.12
mmom   = [0.0, 0.0, -1.2, 0.0]     # ← initial spin
```
Niup と Nidn という名前の二つのサイトがあります。`mmom = [0.0, 0.0, 1.2, 0.0]` は MT の中のスピンモーメントの初期値で、`mmom = [s, p, d, f]` の各値が、原子のサイトでの $l$ ごとのスピンモーメントです。
さらに `[ham]` で `nspin = 2` とします。

計算した MT の中のスピンモーメントは、画面の出力の 'true mm' の列と、`mmom.nio.chk` に
```
# Qtrue    MagMom(up-dn) Rmt   MT
1 8.527587  1.200085  2.120000 Niup    
2 8.527604 -1.200081  2.120000 Nidn    
3 5.380352 -0.000002  1.700000 O       
4 5.380352 -0.000002  1.700000 O    
```
のように出ます。反復ごとに上書きされます。原子のサイトの番号は `SiteInfo.lmchk` にあります。全体の磁気モーメントは
```
 Magnetic moment=     2.241805 !this is a case of bulk Fe
```
のように出ます。

GW の設定（`ctrlg.<sname>.toml` の `[gw]` / `[mlo]` / `[blocks]` / `[product_basis]` の節。[GW の設定](./gwinput)）があれば、
```
gwsc -np 32 5 nio 
```
のように `gwsc` を回せます。32 MPI で QSGW の 1 反復に 275 秒かかりました。次の qsub のスクリプトを使いました。
```
#!/bin/sh
#$ -N NioTest
#$ -pe smp 32
#$ -cwd
#$ -M <your e-mail>
#$ -m be
#$ -V
#$ -S /usr/bin/bash
export OMP_NUM_THREADS=1
EXEPATH=/home/takao/bin/
$EXEPATH/gwsc -np 32 5 nio >& lgwsc
```
計算の後は、`rm __*` か `cleargw .` で一時的なファイルを消します。
`job_band` なども、nspin=1 の DFT と同じように走らせられます。

* `gwsc` で `ctrlg.<sname>.toml` の `[bz].fsmommethod` を 0 以外にしていても、`job_band` を走らせるときは `fsmommethod = 0` にする（またはキーを消す）必要があります。
  `fsmommethod` を 0 以外のまま `job_band` を走らせると、ずれ（磁場のバイアスを足したもの）が入ります。
* [反強磁性の対称性](#antiferro-symmetry-without-soc) を使って計算を効率よくできます（so=0、2 のときだけ）。

## 軌道磁気モーメント
so=1 のとき、MT の中の軌道磁気モーメントが、画面の出力（`IORBTM:` の行から）とファイル `orbitalmom.chk` に次のように出ます。
```
IORBTM:  orbital moments :
 ibas  Spec        spin   Moment decomposed by l ...
    1    Pr           1    0.000000   -0.011483   -0.010535   -4.881144    0.000234
    1    Pr           2    0.000000    0.012332    0.003604    0.000373   -0.000090
 total orbital moment   1:   -4.886708
    2    N            1    0.000000   -0.004174    0.001584    0.000291    0.000129
    2    N            2    0.000000    0.004622    0.000236    0.000045    0.000006
 total orbital moment   2:    0.002738
```

## SOC なしの反強磁性の対称性 {#antiferro-symmetry-without-soc}
反強磁性の対称性を使って、LDA や QSGW の計算ができます（so=1 は除く）。
上向きのスピンだけを計算し、電荷密度と固有値を反強磁性の対称性で対称化します。

反強磁性の対称性のサンプル（NiO と NiSe）を見てください。キー `symgrpaf` と `af` の意味は [lmf](./lmf#symgrp) にあります。手順は次のとおりです。

1. `ctrlg.<sname>.toml` の最上位に `symgrpaf = "find"` と書きます（spglib が操作を見つける。2026-10-02）。
   生成元を書いた `symgrp`（対称性を下げるとき）では、`"find"` は使えないので、反強磁性の生成元を書きます。
   NiO の例: `symgrpaf = "i:( 1 1 1 )"`（`g = i:(1,1,1) * $i \sigma_y$` が磁気空間群を生成する。`i:(1,1,1)` は並進を伴う反転）。
2. 反強磁性の対を `af` で指定します。
   ```toml
   [[site]]
   atom = "Niup"
   pos  = [0.0, 0.0, 0.0]
   af   = 1                  # Antiferro pair
   [[site]]
   atom = "Nidn"
   pos  = [1.0, 1.0, 1.0]
   af   = -1
   ```
3. `pwmode = 11`（平面波の |q+G|）を使います。スピン 2 の関数はスピン 1 の関数を回して作るので、各 q の G の組が要ります。
   `pwmode = 1` では、いくつかの k のメッシュ（NiO の 4x4x4。3x3x3 はたまたま通る）で `rotwave: q+G rotation error` で止まります。

NiSe の例もあります（`symgrpaf = "find"`。古い探し方では `"r2z:( 0 0 0.5 )"`）。(0, 0, 0) と (0, 0, 0.5) の二つの Ni のサイトに `af = 1` と `af = -1` を付けます。

### 反強磁性のモーメントを決まった値に保つ（`mmtarget.aftest`）
反強磁性の対の二つのサイトに逆向きの場を掛け、反復ごとに調整して、モーメント $(m_1-m_2)/2$ を目標の値に保ちます。
QSGW でモーメントが不安定な反強磁性の金属や、全エネルギーをモーメントの関数 $E(m)$ として見るときに使えます。
例: ecalj の `Samples/AFfixMMOM`（NiO、目標 1.6 $\mu_B$。`symgrpaf` あり・なしの両方で、$10^{-6}$ eV 以内で一致）。

* 目標のモーメント（$\mu_B$）の数一つを書いたファイル `mmtarget.aftest` を作業ディレクトリに置きます。このファイルがあるとモードに入ります。
* 場は LDA+U のポテンシャル（`vorb`）に入れるので、対の二つのサイトに LDA+U のブロックが要ります。`idu = 1` で `uh = jh = 0` とすれば、LDA の結果を変えずにブロックを作れます。`sigm`（QSGW）があってもブロックは残るので、QSGW でも使えます。
* 対はサイト 1 と 2（LDA+U のブロック 1 と 2）です。ブロックが 6 つ（NiSe）なら、ブロック 1〜3 と 4〜6 です。
* `lmf` は `mmagfield.aftest`（場 `uhx` とモーメント）と `mixmag.aftest`（場の混合の履歴）を書きます。やり直すときは両方を消します。
  `lmf` を始めるたびに（`job_band` でも）場が一歩進み、`mmagfield.aftest` を書き直します。
* 出力の `mmaftest:` の行に反復、`uhx`（Ry）、$m_1$、$m_2$、... が、`mmhist:` の行に履歴が出ます。
* `ehk` は拘束のもとでの全エネルギーで、$dE/dm_d = 2\,u_{hx}$ です（$m_d$: `dmats` の d のモーメント）。
* モーメントを下げる（NiO で 1.28 から 1.0）と、収束が遅くなります（70 反復くらい。場の更新の利得は固定）。
* `symgrpaf` と使うときは `pwmode = 11` にします（上を参照）。

## スピン軌道相互作用 {#spin-orbit-coupling}
`ctrlg.<sname>.toml` の `[ham]` に `so` のスイッチがあります。

* LDA/GGA では、`nspin = 2` と `so = 1` にすると SOC を入れた計算ができます。so=1 は SOC を入れる、so=2 は LxSx+LySy を無視する LzSz のモードです。
* GaAs のような半導体では、価電子帯の上端のバンド構造を見るのに so=1 が要ります。
* いまは so=1 で QSGW はできません。まず so=0 か 2 で `gwsc` を回し、`sigm` のファイルができてから、摂動として `--ctrlg:ham.so=1` で `lmf` を走らせます（nit=1 でもよい）。
  80% QSGW と組み合わせる手順は [gwsc](./gwsc) の「80% QSGW + 20% LDA と SO の設定」にあります。

* 扱えるのは共線的なスピンだけです。スピンの軸は既定で (0,0,1) で、`[ham] socaxis` でほかの方向も選べます。SOC のサンプル（FePt と MnGa の磁気異方性）を見てください。完全には確かめていませんが、うまく働くようです。

### スピン軌道相互作用を入れたバンドプロット
* 方法 1: バンドプロットにだけ SOC を入れる
  ```bash
  job_band mp-2534 -np 8 --ctrlg:ham.so=1 --ctrlg:ham.nspin=2
  ```
  注意: nspin=2 にすると、rst の大きさは 2 倍になります。nspin=1 の rst には戻せないので、rst を残しておく必要があるかもしれません。

* 方法 2: `[ham].so = 1` で 1 回だけ反復する
  ```bash
  mpirun -np 8 lmf mp-2534 --ctrlg:ham.so=1 --ctrlg:ham.nspin=2 --ctrlg:iter.nit=1
  ```
  `rst.<sname>` が更新されます。そのあと `job_band mp-2534 -np 8 --ctrlg:ham.so=1 --ctrlg:ham.nspin=2` を走らせます。

* 方法 3: `[ham].so = 1` で収束まで反復する
  ```bash
  mpirun -np 8 lmf mp-2534 --ctrlg:ham.so=1 --ctrlg:ham.nspin=2
  ```
  そのあと `job_band mp-2534 -np 8 --ctrlg:ham.so=1 --ctrlg:ham.nspin=2` を走らせます。

有効質量: サンプル `Samples/EffectiveMass/GaAs`（SOC を入れた QSGW のバンド、`syml.<sname>` の質量のモード）。

## フェルミ面
フェルミ面のサンプル（Cu）を見てください。
`lmf cu` が収束した後で `job_fermisurface cu -np 4 -vnk1=10 -vnk2=10 -vnk3=10` を走らせます（`-vnk1` などはスクリプト `job_fermisurface` のオプション）。
BZ の 10x10x10 のすべてのバンドのエネルギー（Ef のまわり）を書き出し、`xcrysden --bxsf fermiup.bxsf` のように xcrysden で見られます。
`job_fermisurface` に `--allband` を付けると、すべてのバンドが出るので、どのエネルギーの等エネルギー面も見られます。

## PROCAR のモード {#procar-mode}
サンプルは ecalj/Samples/PROCAR/MgO_PROCAR にあります。
固有関数の 2 乗の大きさを入れた PROCAR のファイルを作れます。
サンプル（`testecalj MgO_PROCAR`。手順はその `test.py` にある）は、O 2p の成分のファットバンドを示す pdf を作ります。
手順: `job_band mgo` でバンドを描いた後、`lmf --mkprocar --band mgo` で PROCAR（vasp の形式）をランクごとに書き（`PROCAR.UP.<rank>`）、それを k の順につないで `PROCAR.UP` にし、
サンプルのディレクトリにあるスクリプト `BandWeight.py` で重みを取り出して（`bw.dat`）、`gnuplot bnds.gnu.mgoW` で描きます。

## LDA+U
LDA+U のサンプル `Samples/LDAU/ReN`（GdN と PrN。4f に U、`so = 2`、初期の占有 `occnum.<sname>`）と、試験のターゲット `TestInstall/gdn` があります。
`[[spec]]` のキー `idu`、`uh`、`jh` の意味は [lmf の LDA+U](./lmf#lda-u-toml) にあります。例（PrN の Pr）: `idu = [0, 0, 0, 2]`（FLL）、`uh = [0.0, 0.0, 0.0, 0.535]`、`jh = [0.0, 0.0, 0.0, 0.0692]`（Ry）。

## スピンゆらぎ
最近、スピンゆらぎの計算を進展させています。問い合わせてください。
`~/ecalj/Samples/Magnon`（試験のターゲット `Fe_mlo_magnon`。そこの `README.md` に `job_mlo_magnon` の手順がある）。
モデルは [MLO](./mlo) のモデル（窓も反復もなし）です。サンプルは、それを 2026-10-02 より前の版（git のタグ `last-wannier`）の最局在 Wannier 関数と比べます。
次は、[Okumura2021](../presentations/okumura2021.pdf) の Fe のスピンゆらぎの図です（LDA の上で）。
![Fe のスピンゆらぎ](image-6.png)

## 4f と 5f の原子
4f や 5f が部分的に詰まった原子の扱いには、特別な注意が要ります。
4f の原子の atomlist の既定値は 2025-08-21 に更新しました（5f はまだ）。`ctrlgenToml.py` も同じ既定値を使います（表 `SRC/exec/pylib/atomlist.py`）。
試験はまだ限られているので、自分で確かめてください。

ErN のような岩塩構造の希土類の窒化物 ReN の設定は次のとおりです（以前の計算のもの。いま試験しているのは、同じ考え方の `Samples/LDAU/ReN` の GdN と PrN）。
1. `[ham].nspin = 2`、`[ham].so = 2` にします（いまの QSGW は `so=1` を扱えないため）。必要なら、`gwsc` が終わってから `so=1` にします。
   `ctrlgenToml.py gd2pdo4 --so=2 --nspin=2` を使うか、作った後で `ctrlg.<sname>.toml` を直します。
2. `symgrp = "r4z"` にします。格子より低い対称性になるので必要です。スピンと軌道モーメントの向きを固定します。
3. `ctrlg.ern.toml` の設定は次のとおりです。
    ```toml
    [[spec]]
    atom   = "Er"
    z      = 68
    r      = 2.58
    pz     = [0.0, 5.0, 0.0, 5.0]
    idu    = [0, 0, 0, 12]
    uh     = [0.0, 0.0, 0.0, 0.773]
    jh     = [0.0, 0.0, 0.0, 0.0945]
    mmom   = [0.0, 0.0, 0.0, 3.0]
    rsmh   = [1.29, 1.29, 1.29, 1.29]
    eh     = [-1.0, -1.0, -1.0, -1.0]
    rsmh2  = [1.29, 1.29, 1.29, 1.29]
    eh2    = [-2.0, -2.0, -2.0, -2.0]
    lmx    = 3
    lmxa   = 6
    nmcore = 1
    kmxa   = 5
    # ([ham] pwemax = 2, [bz] nkabc = [8, 8, 8] was used)

    [gw]
    n1n2n3 = [6, 6, 6]

    [product_basis]
    pb_tolerance = [1e-2]
    pb_lcutmx    = [6, 2]
    ```
    `pz` はローカル軌道の主量子数で、ここでは 5p と 5f をローカル軌道にしています。原子の配置は `lmfa ern|grep conf` で確かめます。
    * `idu`、`uh`、`jh` は LDA+U のパラメータです。`idu` が 10+2 のときは、`sigm` があれば U の効果を使いません。
    * `mmom` は磁気モーメントの初期値です。
    * augmented な波の展開に `lmxa = 6` が要ります。
    * `[ham] readp = true` であることを確かめてください。

**注意**: `idu` > 10 のとき、LDA+U は QSGW の初期条件にだけ使います。そのため QSGW の結果は `jh` と `uh` によりません。

この設定で、4f の ReN の一連の計算ができます。入力のファイルは [ここ](../data/RENctrlGWinput.tar.gz) にあります。使いたいときはお問い合わせください。
