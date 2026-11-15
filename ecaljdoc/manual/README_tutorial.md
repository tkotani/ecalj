# ecalj メインドキュメント {#ecalj-maindocument}
**ecaljdoc のメインドキュメントです。ecalj の特徴、QSGW の概観、最初の計算（GetStarted）を書いています。ほかのページへは左のメニューからたどってください。**

* ここでは [GetStarted](#getstarted) を、インストールと QSGW の概観とあわせて説明します。
* [UsageDetails](./UsageDetailed.md) は別のファイルにあります。
* [Qiita の日本語の記事](https://qiita.com/takaokotani/items/9bdf5f1551000771dc48) も参考になるかもしれませんが、ほとんどはここに書いてあります。
* リンクが切れていたら、https://github.com/tkotani/ecalj からダウンロードしてください。
* インストールと案内は AI にさせてください（[GetStarted with AI](./getstartedAI) を AI に読ませる）。インストールの手順の文書は [install](../install/install.md)、[ISSP でのインストール](../install/installISSP.md)。
* ライセンス: AGPLv3（AI にも例外なく適用）。[ライセンス](./license) のページを見てください。

## ecalj パッケージの特徴 {#features-of-ecalj-package}

1. **全電子フルポテンシャルの PMT 法** 
   PMT 法とは、2 種類の補強された波（augmented waves）、つまり APW+MTO の混合基底の方法です。
   言い換えると、PMT 法 = 線形化した（APW+MTO）法です。ecalj と同じ起源をもつ [Questaal](https://www.questaal.org/) を除けば、ほかにない方法です。MTO と APW は、波動関数の表現において相補的に機能します。それぞれ固有関数の局在した性質と広がった性質に対応します。

   図 1 を見てください。<br>
   MTO と APW では、MT の外の包絡関数（図の破線）の曲率が逆になっています。曲率の符号は、その場所の $V-E$ の符号で決まります（$\nabla^2\psi/\psi = (2m/\hbar^2)(V-E)$）。$E>V$ では振動し（APW）、$E<V$ では減衰します（MTO）。3d のような局在した状態の固有値は、MT の外の有効ポテンシャル（遠心力の項 $l(l+1)/r^2$ を含む）の下にきます。一方、Na の 3s のような広がった状態の固有値は、その上にきます。要するに、固体中の波動関数の曲率は、MT の外の領域でも、平面波的（振動する）にも局在的（減衰する）にもなるということです。これらに対応した基底関数が APW と MTO です。
   局在した MTO（減衰因子 $\exp(-\kappa r)$、$\kappa \sim 1\,{\rm bohr}^{-1}$。最近接の原子までしか届かないということです）と APW（カットオフは $\approx 3$ Ry）を合わせると、十分な収束が得られます。APW を含むので、散乱状態もとてもよく記述できます。これは Gaussian などでは困難な点です。一方で、極端に少ない APW でも、3d などのコアの形に近い原子の波動関数を記述できるのは、MTO のおかげです。GGA/LDA のレベルで原子の位置の緩和もできます。 
   ![MTO と APW](image-1.png)

   *図 1*. MTO と APW（M. Methfessel と M. van Schilfgaarde による [nfp のマニュアル](../presentations/nfpmanual.pdf) から）。

   いまの PMT の定式化は次にあります。

   [1][KotaniKinoAkai2015, PMT formalism](../presentations/KotaniKinoAkai2015FormulationPMT.pdf)
   
   [2][KotaniKino2013, PMT applied to diatomic molecules](../presentations/KotaniKino2013PMTMolecule.pdf).

   ただ、セッティングや細かい改良はいろいろと加えられています。とくに、MTO の設定が自動化されている点が重要です。通常の LMTO 法では、MTO のパラメータの設定が難しく、かなり使いづらいものでした。
   基底関数のパラメータは自動で設定するので、パラメータの設定に悩む必要はありません。MTO に APW を加えたことで、MTO のパラメータへの依存がかなり減ったのです。従来の MTO だけの方法では、非常に長距離の（1/r 型に近い）基底を入れる必要もありました。計算に要るのは結晶構造（POSCAR）だけです。 
   <!-- In principle, it is possible to perform reasonable calculations just from crystal structures and very minimum setting.  -->
    * PMT のはじまりは https://journals.aps.org/prb/abstract/10.1103/PhysRevB.81.125117 です。
    * PMT は [E. Bott、M. Methfessel、W. Krabs、P. C. Schmidt による smooth Hankel 関数の論文](../presentations/Bott1988.pdf) にある smooth Hankel 関数を使います。これは [B][M. Methfessel、M. van Schilfgaarde、R.A. Casali による nfp の論文](../presentations/nfp.pdf) で使われたもので、私たちの PMT はその上に作られています。
   *  PMT の基底に加えて、局所軌道（local orbital）も合わせて使います。
   * カットオフ $\approx 3$ Ry は驚くほど小さく聞こえるかもしれません。しかし、平面波の方法で高いエネルギーの部分は、動径関数の曲率が負になる、固有関数の内殻のような形を表すためだけに使われていることに注意してください。さらに、平面波には周期ゲージ（periodic gauge）を満たさないという難しさがあります（[arXiv:2010.03959](https://arxiv.org/pdf/2010.03959)）。ecalj では MTO のおかげで、最大局在ワニエ法に対応する方法（[MLO モデル化](./mlo)）をシンプルな方法で実行できます。一方、LMTO やガウス関数のような局在基底の方法は、真空の領域を扱うのが苦手です。表面の扱いなどでは問題になりえます。

2. **PMT-QSGW 法** 
   
   QSGW 法とは、準粒子バンドを的確に与える有効一体ハミルトニアンを計算する方法です。
   PMT-QSGW とは、`the Quasiparticle self-consistent GW method (QSGW) based on the PMT method`、つまり PMT 法にもとづく準粒子セルフコンシステント GW 法（QSGW）のことです。PMT の上で QSGW を行うための鍵は [3] の開発です。磁性金属も扱えます。ecalj は GPU に実装したので、4 枚の GPU で約 40 原子を扱えます。
   [3][Kotani2014, Formulation of PMT-QSGW method](../presentations/Kotani2014QSGWinPMT.pdf)
   [4][D.Deguchi PMT-QSGW applied to a variety of insulators/semiconductors](../presentations/deguchi2016.pdf)
   [5][M.Obata GPU implementation](https://arxiv.org/abs/2506.03477)。ここでは Type II の GaSb/InAs（40 原子）を 4 枚の GPU で扱っています。
   Wannier などの内挿をしなくても、バンドの図を簡単に描けます。自己エネルギーの内挿の仕組みが内部に組み込まれているからです。
   これは、QSGW の基底状態の上で固有関数から物理量を計算するのにとても大事です。
   （MLO をもとにした QSGW 法（[MLO-gwsc](./mlo_gwsc)）も動きますが、まだ試験できていません。）

   すべての電子を扱うので、擬ポテンシャルのあいまいさに悩まされません。

3. **誘電関数と磁化率**
    誘電関数、グリーン関数のスペクトル関数、磁気ゆらぎなど、GW に関係する量を計算できます。QSGW は任意の k 点で固有関数と固有値を作れるからです。

4. **MLO によるモデルハミルトニアン** 
   MLO（muffin-tin based localized orbitals、マフィンティンにもとづく局在軌道）で有効模型を作ります。模型のハミルトニアンと、MLO の間の有効相互作用 v、W、cRPA の W です
   （`job_mlo`、`job_mloW --crpa`。[MLO](./mlo)）。
   * 2026-10-02 まで、ecalj には cRPA つきの最局在 Wannier 関数（Juelich のグループの定式化）がありました。最局在 Wannier 関数の部分は、
     Miyake、Sakuma、Kino らによるコードから出発したものです（これを使った論文: [Phys. Rev. B 108, 035141](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.108.035141)）。
     MLO がこれに置き換わりました。もとのものは git のタグ `last-wannier` にあります。 

5. **Fortran の singleton パターン**
   コードの主な部分は Fortran です。コードは singleton のデザインパターンにもとづくオブジェクト指向の構造をもっています。この構造は、AI にとっても理解しやすく、間違いにくいという利点があります。データは一つの module にだけ置かれ、どこで何を使うかは `use ..., only:` で明示されるので、データの流れを追いやすく、写しの食い違いも起きません。

## QSGW の概観 {#overview-of-qsgw}

* バンド計算（LDA のレベル）はプログラム `lmf` で行います。最初の設定ファイルは `ctrlg.<sname>.toml` です（`<sname>` は利用者が選ぶ物質の拡張子。古い形は `ctrl.<sname>`）。`lmf` の前に `lmfa` を走らせる必要があります。これは球対称の原子の計算で、電子密度の初期条件を決めます（`lmfa` は数秒で終わります）。
* ファイル `sigm.<sname>` が QSGW 計算の鍵です。`sigm.<sname>` には非局所ポテンシャル $\Delta V_{\rm xc}=V_{\rm xc}^{\rm QSGW}-V_{\rm xc}^{\rm LDA}$ が入っています。`lmf` が行う通常の LDA 計算にこのポテンシャルの項を加えることで、QSGW 計算ができます。下の図を見てください。
* 問題は $V_{\rm xc}^{\rm QSGW}({\bf r},{\bf r}')$ をどう作るかです。これは GW 近似で計算した自己エネルギー $\Sigma({\bf r},{\bf r}',\omega)$ から計算します。大まかに言えば、$\Sigma({\bf r},{\bf r}',\omega)$ の omega 依存性を取り除いて $V_{\rm xc}^{\rm QSGW}({\bf r},{\bf r}')$ を得ます。
* したがって、$V_{\rm xc}^{\rm QSGW}$ の計算が QSGW のサイクルの主な部分で、二重のループで計算します。つまり、`lmf` の内側のループと、`lmf` が与える固有関数を使って $V_{\rm xc}^{\rm QSGW}$ を計算する外側のループがあります。外側のループは gwsc という python のスクリプト（Fortran のプログラムを走らせます）で実行できます。QSGW の計算時間は LDA の計算よりずっと長くなります。 
目安として、QSGW の反復 1 回に、20 原子でおよそ 10 時間かかります 
（電子の数によります。セルあたり 10 原子くらいより大きい系では GPU を使うことを勧めます）。ですが、GPU の混合精度の計算なら 30 分程度以下だと思います。
* 次は https://arxiv.org/abs/2506.03477 の Figure 1 に示した QSGW のサイクルです。MPB は、固有関数の積を展開するための混合積基底（mixed product basis）のことです。 
![alt text](image-4.png) 
* [QSGW の GPU による高速化](https://arxiv.org/abs/2506.03477) があります。この論文は QSGW の基本も説明しています。QSGW のアルゴリズムは GPU の計算にとてもよく合います。4 枚の GPU で、セルあたり 40 原子の系を計算できます（lmf の部分では、まだ GPU を効率よく使えていません）。PMT ではスラブ模型の大きな真空の領域も扱えます。
* 手でパラメータを設定することなく、ほぼ自動で QSGW を行えます。そのため、ecalj は GW/QSGW を行うのに最も簡単なものの一つだと考えています。QSGW のバンドのデータベースは [DOSnpSupplement](https://github.com/tkotani/DOSnpSupplement) を見てください。

![alt text](image-2.png) この図は [4][D.Deguchi](../presentations/deguchi2016.pdf) から取りました。

LDA と比べると、QSGW では次の違いが見られます。
* バンドギャップ。QSGW は実験よりやや大きくなる傾向があります。上の図のように系統的に見えます。
* バンド幅。ふつう sp バンドは広がります（Na のようなとても密度の低い場合を除く）。
    これは一様電子ガスの場合と整合しています。
    3d 電子のような局在したバンドは狭くなります。
* バンドの相対的な位置。例えば O(2p) と Ni(3d)。
より局在したバンドほど深くなる傾向があります。
上向きと下向きの間の交換分裂は大きくなります（LDA+U のように）。
NiO のような場合、磁気モーメントは大きくなり、実験値に近づきます。
* 3d バンドとほかのバンドの混成。QSGW は固有関数を局在させる傾向があります。

ただし実際は複雑で、場合によってはそう単純ではありません。

<!-- 
### 基礎知識check 
* LDA計算の流れ
* MT division of space
* 基底関数 MTOとAPW　（AugmentationとMT内の原子的な波動関数). 3つの成分をもつ。envelope + true - counter
* バンド計算法PMT=LMTO+LAPW法. ３つの成分で電荷も表現する。MTOの張る空間の不満足な部分を3Ry以下程度のAPWでサポートする。
* 独立粒子近似とは
* Hartree-Fock(HF)近似で水素Hがきちんと解ける。一様ガスではフェルミ面で状態密度がゼロになる。
* LDAでは一様ガスがきちんと解けるがHのバンドギャップはかなり小さい。
* LDAとHFは両極にある。HFではSiのバンドギャップは10eV以上（なぜか？）　ハイブリッド法HSE。自己相互作用。
* 交換項の重要性。H2のbonging-antibondingを区別するには、占有状態へのプロジェクタの役割をするFock項が重要
* GW近似. Fock項にくわえて相関項。分極媒質中を走る荷電粒子の感じる時間依存ポテンシャル。スクリーンHF+クーロンホール
　バンドギャップ、dバンドの位置、Uの効果
* QSGW法。GW近似における自己エネルギーからいくらか強引に時間依存性をとりのぞいて特殊な交換相関項をつくる。
* （できたら線形応答理論の基礎、ステップ関数のフーリエ変換、古典的な減衰振動系でのグリーン関数（線形応答関数）。）

### ecalj 何ができるか？

* LDA計算 
VWN,PBE-GGA,LDA+Uがえらべる。
構造緩和（格子変化は手動）、Colinear,  SOC(軸を選べる）, AFの対称性を入れることができる。
ESM法でのスラブ計算. 最適化しきれてないので現状では構造緩和などでは優位性がない。

* QSGW計算：
self-consistent GW法。特殊な交換相関項を作る計算であるといえる。
自己エネルギースペクトルプロット。
インパクトイオン化率（オージェによる寿命）。
GPU化、自動化セッテイング。
GW法でのバンドプロットが直接に可能。
QSGWでは現状全エネルギー計算ができない。バンド構造（固有値、波動関数）のみ。
* 線形応答などの計算。
RPAでの誘電率計算、スピンゆらぎ計算 （改良の余地。金属でもできる。ドルーデウエイト(q→0)）
MLO による模型化（模型のハミルトニアン、v・W・cRPA。2026-10-02 に最局在 Wannier から置き換えた）。 自動化がすこしできてないところがある。
その他の物理量についても応用できるはず。
* かなりの部分で自動計算が可能。QSGW法ではMaterial Projectから1500個程度の構造ファイルを持ってきて自動化でQSGW計算しているが
ほぼ問題なく可能。個別にセッティングを手動でいじらなくても良い（4f,5fについては自動化がまだ設定できてないが基本的に可能）
バンドプロットも対称ラインも含め自動化してある。データはgnuplotなどでプロットするので読みやすい。結晶構造についてはPOSCARとの相互コンバータあり。
複数のPOSCARを一括計算するecalj_autoも梱包してある（整備中）。 -->



## 始めよう（GetStarted） {#getstarted}
ここでは ecalj による DFT/QSGW の計算を説明します。そのあとでバンドの図の描き方を説明します。簡単のため、常磁性（nsp=1）で、4f なし、SOC なしの場合を扱います。
一歩ずつ説明します。

> **実例** — 以下のすべての Step で **GaAs**（物質の名前 sname は `gaas`）を例に使います。
> 最小のシード（`ctrls.gaas` と、生成した `ctrlg.gaas.toml`）は
> [`Samples/GetStarted/GaAs/`](https://github.com/tkotani/ecalj/tree/main/Samples/GetStarted/GaAs) にあります。
> Step 0 のとおり作業の場所にコピーして、Step 1-6 に従ってください。一度に再現するコマンドの並びは、その
> [README](https://github.com/tkotani/ecalj/blob/main/Samples/GetStarted/GaAs/README.md)
> にあります。

詳しくは [UsageDetailed](./UsageDetailed.md) で説明します。

### Step 0. POSCAR を用意する {#step-0-get-poscar}
GaAs の実例では、サンプルのディレクトリを ecalj のクローンの**外**の作業の場所にコピーして、そこで計算します（計算で作られるファイルがクローンの中に散らばらないようにするためです）。
```bash
mkdir -p ~/work
cp -r ecalj/Samples/GetStarted/GaAs ~/work/gaas    # ecalj はクローンの場所
cd ~/work/gaas
```
ここには `ctrls.gaas`（Step 1 で作る構造のファイル）と `ctrlg.gaas.toml`（Step 2 で作る入力ファイル）がもう入っています。

自分の物質では、POSCAR（VASP の形式の結晶構造）から始めます。
Materials Project（MP）の POSCAR の例が ecalj/ecalj_auto/INPUT/testSGA/POSCARALL にあります（物質の一覧は同じ場所の `joblist`）。
```bash
cat ecalj/ecalj_auto/INPUT/testSGA/joblist
mkdir -p ~/work/mp-2534
cp ecalj/ecalj_auto/INPUT/testSGA/POSCARALL/POSCAR.mp-2534 ~/work/mp-2534/
```

例えば、mp-2534 の GaAs の POSCAR は次のとおりです。
```
Ga1 As1
1.0
   3.5212530000000002    0.0000000000000000    2.0329969999999999
   1.1737510000000000    3.3198690000000002    2.0329969999999999
   0.0000000000000000    0.0000000000000000    4.0659929999999997
Ga As
1 1
direct
   0.0000000000000000    0.0000000000000000    0.0000000000000000 Ga
   0.2500000000000000    0.2500000000000000    0.2500000000000000 As

```
この格子定数は 5.75 Å で、実験値の 5.65 Å より大きくなっています。MP は PBE で構造を緩和するためです（GGA は格子定数を大きめに出します）。
実例の `ctrls.gaas` は実験値 a=5.65325 Å を使います。

cif を `POSCAR` に変換したいときは、次のようにします。
`cif2cell foobar.cif -p vasp --vasp-cartesian --vasp-format=5`


### Step 1. POSCAR を ctrls に変換する {#step-1-convert-poscar-to-ctrls}
次に vasp2ctrl で POSCAR を ctrls に変換します。
ctrls は ecalj で使う構造のファイルです。GaAs の実例では `ctrls.gaas` がすでにあるので、この Step は要りません。

`POSCAR.foo` から始めるときは次のようにします。`foo` が物質の名前（sname）になり、この後のコマンドで使います。
```
vasp2ctrl POSCAR.foo
mv ctrls.POSCAR.foo.vasp2ctrl ctrls.foo
```

例えば mp-2534 では `vasp2ctrl POSCAR.mp-2534` のあと `mv ctrls.POSCAR.mp-2534.vasp2ctrl ctrls.mp-2534` とします。
ctrls.mp-2534 には POSCAR と同じ結晶構造が入っています。
```
STRUC
     ALAT=1.8897268777743552
     PLAT=       3.52125300000       0.00000000000       2.03299700000  
                 1.17375100000       3.31986900000       2.03299700000 
                 0.00000000000       0.00000000000       4.06599300000 
SITE
     ATOM=Ga POS=     0.00000000000       0.00000000000       0.00000000000 
     ATOM=As POS=     1.17375100000       0.82996725000       2.03299675000 
```
STRUC と SITE のような区分のブロックを示すために、字下げが必要です。
実例の `ctrls.gaas`（[`Samples/GetStarted/GaAs/ctrls.gaas`](https://github.com/tkotani/ecalj/blob/main/Samples/GetStarted/GaAs/ctrls.gaas)）は、手で書いた ctrls の例です。

```
%const bohr=0.529177 a=5.65325/bohr
STRUC
     ALAT={a}
     PLAT=0 0.5 0.5 0.5 0 0.5 0.5 0.5 0
SITE
    ATOM=Ga POS=0.0 0.0 0.0
    ATOM=As POS=0.25 0.25 0.25
```
ctrls.srtio3:
```
%const au=0.529177
%const a0=2*1.95/au v=a0^3 a1=v^(1/3)
STRUC ALAT={a1}
   PLAT=1 0 0 0 1 0 0 0 1
SITE
   ATOM=Sr POS=1/2 1/2 1/2
   ATOM=Ti POS= 0 0 0
   ATOM=O POS=1/2 0 0
   ATOM=O POS= 0 1/2 0
   ATOM=O POS= 0 0 1/2
```

`%const` の行は、`ctrlgenToml.py` が展開する変数を定義します。1 行に複数の定義を書いて
かまいません。演算子（`+ - * / ** ^`）と `sin / cos / log / sqrt / pi`
は評価されます。先に定義したものは後の定義から見えます（そのため、
`au` が前の `%const` の行で定義されているので
`a0=2*1.95/au` が使えます）。`1/3` は浮動小数点の割り算（Python 3）なので、
`v^(1/3) ≡ v**(1/3)` は立方根です。

> **落とし穴** — `%const x= y=...` のような、値の空のプレースホルダは
> **使わないでください**。右辺のない項目は次の語とくっついてしまい、
> 行の残りが黙って壊れます（そのため `y`、`z`、... は定義されません）。
> 空の項目を消すか、実際の値を入れてください
> （`%const x=0 y=...`）。

原子の種の既定の名前は `ctrlgenToml.py --showatomlist` で表示されます。
既定の記号の代わりに、自分の記号を次のように使えます。
```
SITE
   ATOM=M1 POS=1/2 1/2 1/2
   ATOM=M2 POS= 0 0 0
   ATOM=O POS=1/2 0 0
   ATOM=O POS= 0 1/2 0
   ATOM=O POS= 0 0 1/2
SPEC
   ATOM=M1 Z=38
   ATOM=M2 Z=22
   ATOM=O Z=8
```
この場合は SPEC が必要です。これを使うと、反強磁性の NiO は次のようになります。
```
%const bohr=0.529177 a=7.88
STRUC 
   ALAT={a} PLAT= 0.5 0.5 1.0 0.5 1.0 0.5 1.0 0.5 0.5
SITE 
   ATOM=Niup POS= .0 .0 .0
   ATOM=Nidn POS= 1.0 1.0 1.0
   ATOM=O POS= .5 .5 .5
   ATOM=O POS= 1.5 1.5 1.5
SPEC
   ATOM=Niup Z=28 MMOM=0 0 1.2 0
   ATOM=Nidn Z=28 MMOM=0 0 -1.2 0
   ATOM=O Z=8 MMOM=0 0 0 0
```
反強磁性の対称性を使いたいとき（上向きのバンドだけを計算し、下向きのバンドを作る）は、
[AFsymmetry](./UsageDetailed.md#antiferro-symmetry-without-soc) を見てください。MMOM は各サイトの磁気モーメントの初期条件です。
数値（ここでは 1.2）はあまり正確でなくてかまいません（整数くらいで十分です）。
例は [Samples/MATERIALS](#samples-materials-lda-and-mlo-samples-of-65-materials) にあります。まずは簡単な物質から始めてください。ただ、
自分の POSCAR から VASP や MP の結果を再現するのは、それほど難しくありません。


- メモ: 
    - `ctrl2vasp ctrls.gaas` で VASP のファイルに戻せます。VESTA で確かめてください。viewvesta（変換して VESTA を起動します）も使えます。
    - 使わないファイルがたくさん作られます（気にしなくてかまいません）。
    - サイトには [Niup などの](./UsageDetailed.md#antiferro-symmetry-without-soc) 好きな名前を使えます。その場合は SPEC の区分も設定する必要があります。[Z は整数でなくてもかまいません。](./Memo_bgcharge.md#background-charge-and-fractional-z) これは後で学んでください。
    - 古い ctrls には NL、NBAS、NSPEC があるかもしれませんが、いまは要りません。

### Step 2. ctrls から `ctrlg.<sname>.toml` を作る {#step-2-get-ctrlg-sname-toml-from-ctrls}

`ctrlg.<sname>.toml` は、2026-05 から Fortran のプログラム（`lmf`、`lmfa`、`lmchk`、`gwsc`、...）
が読むただ一つの入力ファイルです。Step 1 の ctrls から **`ctrlgenToml.py`** で作ります。

```bash
ctrlgenToml.py gaas
```

このコマンドは次のことをします。

1. 最上位（`symgrp` / `verbose` / `time`）と `[struc] / [[site]] / [[spec]]` を、周期表の
   既定値（`atomlist` の表、`SRC/exec/pylib/atomlist.py`）で埋めます。
2. 内部で `lmfa → lmf --jobgw=0 → gwinit` を走らせて、
   `[gw] / [mlo] / [blocks] / [product_basis]`（最後のものは
   原子ごとの表 `nlx` / `valence` / `core` を含みます）を埋めます。
3. `ctrlg.<sname>.toml` を書いて終わります。

うまく終わると、最後は次のようになります。

```text
ctrlgenToml: wrote ctrlg.gaas.toml (2 spec, 2 sites)
ctrlgenToml: running lmfa -> lmf --jobgw=0 -> gwinit  to fill GW sections
ctrlgenToml: done. ctrlg.gaas.toml has top-level/[struc]/[[site]]/[[spec]]/...
             plus [gw]/[mlo]/[blocks]/[product_basis] (incl. nlx/valence/core).
```

実例のディレクトリにある
[`ctrlg.gaas.toml`](https://github.com/tkotani/ecalj/blob/main/Samples/GetStarted/GaAs/ctrlg.gaas.toml)
は、こうして作ったものです。走らせ直すとファイル全体を ctrls から書き直します（前のファイルは
`ctrlg.<sname>.toml.bakup` に移しますが、控えは一世代だけです）。

**勧める手順**: `ctrlgenToml.py <sname>` を**ほかのフラグなしで**走らせて
既定値を得て、**そのあとで `ctrlg.<sname>.toml` を編集します**。
ファイルには説明のコメントと TOML の型が付いているので、エディタでその場で直すのが自然なやり方です。
`ctrlgenToml.py --help` に出るフラグは、こうした編集をあらかじめ焼き込むためのもので（多くの物質をまとめて作るときに便利です）、必須ではありません。
`[product_basis]` の末尾の `nlx` / `valence` / `core` の表は GW の計算だけで使います。ふつうは手で編集しません。

`ctrlg.<sname>.toml` でよく編集するものです。

1. `[bz].nkabc` — k 点のメッシュ（LDA）。
2. `[ham].nspin` — 磁性なら 2、非磁性なら 1。
3. `[ham].so` — 0 / 1（完全な L·S）/ 2（Lz·Sz だけ）。`nspin = 2` が必要です。
   `so=1` はまだ QSGW に対応していません。QSGW + SOC では、`so=0`
   か `so=2` で走らせて `sigm.<sname>` を得てから、`so=1` にします。
4. `[ham].xcfun` — 1 = VWN（既定）、2 = Barth-Hedin、103 = PBE-GGA。
5. `[ham].scaledsigma` — QSGW の混合。`1.0` = 完全な QSGW、`0.8` = QSGW80。
   `sigm.<sname>` のファイルがあれば、`(scaledsigma) × (V_xc^QSGW − V_xc^LDA)` が lmf の
   ハミルトニアンに加わります。
6. `[gw].n1n2n3` — GW の k 点のメッシュ。`[bz].nkabc` より小さくします（1/2 か 2/3）。
7. LDA+U: `[[spec]]` の中に配列 `idu / uh / jh` を加えます。
   [`lmf.md` の LDA+U の節](./lmf) を見てください。

[すべてのキーの詳しい説明](./lmf) があります。
GW を使わないときの `--skipgw`、手で直したところを残したままあとから GW の区分を足す `--addgw`、
一回の実行だけキーを変える `--ctrlg:<section.key>=<value>`（例 `lmf gaas --ctrlg:bz.nkabc=[8,8,8]`。古い `-v` の形ではプログラムが止まります）は、
[入力ファイル](./toml_migration) を見てください。

* `lmchk gaas --ctrlg:verbose=60` で、認識された対称性を確かめられます（古い `lmchk --pr60 gaas` は廃止しました）。

ここで、次のものを目で確かめられます。

* `SiteInfo.lmchk` — MT 半径、原子の位置
* `PlatQlat.chk` — 基本格子ベクトル（plat）と逆格子（qlat）

#### Step 2-Migration. すでにある `ctrl.<sname>`（と、あれば `GWinput`） {#step-2-migration-existing-ctrl-sname-and-optional-gwinput}

手で編集した `ctrl.<sname>`（と、あれば `GWinput`）のある 2026-05 より前のディレクトリは、ctrls から作り直さず、`Legacy2toml.py <sname>` でその場で `ctrlg.<sname>.toml` に変換します。
診断の出力とキーの対応は [入力ファイル](./toml_migration#migration-of-an-existing-legacy-directory) を見てください。

#### VESTA をインストールする {#install-vesta}
結晶構造を見るには VESTA が便利です。
ecalj/StructureTool/ に 'viewvesta' のコマンドがあります。
次を試してください。
```
viewvesta ctrlg.gaas.toml
viewvesta ctrls.gaas
```
これで構造を VESTA で見られます。変換のプログラム `vasp2ctrl` / `ctrl2vasp`（これも `/StructureTool` にあります）も、同じく TOML / ctrls / ctrl の形の入力を受け付けます。
<small>（~/ecalj/GetSyml/README.md もありますが、利用者は読まなくてかまいません。）</small>



### Step 3. LDA の計算 {#step-3-lda-calculation}
この Step に載せた出力（`conf` の表、反復のエネルギー、`bndfp` の行）は、mp-2534（格子定数 5.75 Å）で走らせたときのものです。gaas では数値が少し違います。

1. まず lmfa を走らせます。結晶に含まれる、球対称の原子の電子密度を求めるものです。lmfa は数秒で終わります。
   ```bash
   lmfa gaas
   ```
   で、初期化のための球対称の原子の計算をします。`lmfa` は
   `ctrlg.gaas.toml` を読み（cwd に `ctrlg.*.toml` がちょうど一つだけあれば、
   名前を省いても自動で見つけます）、下の `lmf` に要る
   ファイルを作ります。
   コンソールの出力の `conf ` の部分を、次のようにして確かめてください。
   ```bash
   lmfa gaas |grep conf
   ```
   これで原子の配置が表示されます（`lmfa` を何度走らせても副作用はありません）。
   ```
   conf:------------------------------------------------------
   conf:SPEC_ATOM= Ga : --- Table for atomic configuration --- (n and nz are Pr. Quantum numbers for valence) 
   conf:  isp  l    n nz    Qval     Qcore   CoreConf 
   conf:    1  0    4  0    2.000    6.000 => 1,2,3,
   conf:    1  1    4  0    1.000   12.000 => 2,3,
   conf:    1  2    4  3   10.000    0.000 => 
   conf:    1  3    4  0    0.000    0.000 => 
   conf:    1  4    5  0    0.000    0.000 => 
   conf: Species  Ga       Z=  31.00 Qc=  18.000 R=  2.230000 Q=  0.000000 nsp= 1 mom=  0.000000
   conf: rmt rmax a=  2.230000  49.011231  0.015000 nrmt nr= 661 867
   conf: Core rhoc(rmt)= 0.001940 spillout= 0.001839
   conf:------------------------------------------------------
   conf:SPEC_ATOM= As : --- Table for atomic configuration --- (n and nz are Pr. Quantum numbers for valence) 
   conf:  isp  l    n nz    Qval     Qcore   CoreConf 
   conf:    1  0    4  0    2.000    6.000 => 1,2,3,
   conf:    1  1    4  0    3.000   12.000 => 2,3,
   conf:    1  2    4  0    0.000   10.000 => 3,
   conf:    1  3    4  0    0.000    0.000 => 
   conf:    1  4    5  0    0.000    0.000 => 
   conf: Species  As       Z=  33.00 Qc=  28.000 R=  2.330000 Q=  0.000000 nsp= 1 mom=  0.000000
   conf: rmt rmax a=  2.330000  49.695135  0.015000 nrmt nr= 675 879
   conf: Core rhoc(rmt)= 0.011286 spillout= 0.017467
   ```
   この表から、Ga の初期条件の価電子の数は 4s、4p、4d+3d（局所軌道を使う）でそれぞれ 2、1、10 だとわかります。 
   内殻の配置と電荷も表示されています。Q=0 は電荷の中性を満たしていることを意味します。 
   * lmfa は球対称の原子を（実質的に $r \to \infty$ まで）計算し、MT の境界での対数微分を計算して `atmpnu.*` に保存します。 
   * lmf の初期の電子密度は、lmfa が与える球対称の原子の密度の重ね合わせです。 
   * [`[ham].readp = true`](./lmf.md#readp-option)（`ctrlgenToml.py` が既定で設定します。古い `READP=T` も同じでした）のときは、MT の境界での動径波動関数の対数微分を固定したままにします。

2. `lmfa` のあと、次のようにして LDA の計算を走らせます（`tee` でコンソールの出力を `llmf` にも書きます）。

   ```bash
   mpirun -np 8 lmf gaas | tee llmf
   ```
   数秒で終わります。すると次のように表示されます。   
   ```bash
   it 10  of 80    ehf=   -8398.499465   ehk=   -8398.499466
   From last iter    ehf=   -8398.499464 ehk=   -8398.499465
   diffe(q)= -0.000001 (0.000004)    tol= 0.000010 (0.000010)   more=F
   c ehf(eV)=-114268.304015 ehk(eV)=-114268.304037 sev(eV)
   ```
   ここで `c ` は、前の反復からのエネルギーと電荷の変化を示す `diffe(q)` が基準より小さい、つまり収束したという印です。収束までに 10 回かかったことがわかります。[二つのエネルギー](./lmf.md#q-should-the-harris-foulkes-and-hohenberg-kohn-sham-functionals-agree-at-self-consistency) を表示していますが、これらは同じになるはずです（ただし少し違います）。
   * `llmf` には反復の情報が入っています。固有値、フェルミエネルギー、バンドギャップを確かめてください。
   * `rst.gaas` が作られます。セルフコンシステントな電荷が入っています。
   * `save.gaas` には lmfa と lmf の履歴が入ります。プログラムを始めるたびにコマンドライン（オプションを含む）を `Start` の行に書き、lmf は反復ごとに 1 行を書きます。
     反復の行の先頭の文字は状態を示します。`c` は収束、`x` は収束しないまま反復の上限に達した、`i` は反復の途中、`h` は 1 回目の反復（自由原子の密度の重ね合わせからの Harris のエネルギー）、`C` は構造の緩和の収束です。
   * lmf をもう一度走らせると、rst のファイルを自動で読みます。そのため、もう収束しているはずです。それを確かめるために、2 回の反復で止まります。
   * 格子定数は `ctrlg.gaas.toml` の `[struc].alat`（a.u.）を編集して
     変えられます。TOML の値の中の演算子は
     評価されません。計算した値を書く（例えば `alat = 10.6818`）か、
     実行時に `--ctrlg:struc.alat=10.6818` を使ってください。ctrls や古い `ctrl.<sname>` の
     `%const` の計算（`1.8897268777743552*5.65/5.75` など）は、
     `ctrlgenToml.py` や `Legacy2toml.py` が値にして書き込みます。
   * 古い版にあった中間のファイル `ctrlp.<sname>`（`ctrl2ctrlp.py` が作っていたもの）は、いまはありません。Fortran が `ctrlg.<sname>.toml` を直接読みます。
   * SiteInfo.lmchk : サイトの情報
   * PlatQlat.chk : 格子の情報
   * efermi.lmf: フェルミエネルギー 
   * estaticpot.dat : なめらかな部分の静電ポテンシャル。
   * __mixm* は反復の履歴を入れた混合のファイルです。

一体ポテンシャルのゼロの基準の定義に注意してください。MT の境界での静電ポテンシャルの平均をゼロにしています。 
`__*` は一時的なファイルで、消してかまいません。 
コンソールの出力には、次のようなものが見られます。
```
...
 bndfp: kpt     1 of    29 k=  0.0000  0.0000  0.0000 ndimh = nmto+napw =    82   55   27
 bndfp: kpt     2 of    29 k=  0.0355 -0.0126  0.0000 ndimh = nmto+napw =    79   55   24
 bndfp: kpt     3 of    29 k=  0.0710 -0.0251  0.0000 ndimh = nmto+napw =    82   55   27
 bndfp: kpt     4 of    29 k=  0.1065 -0.0377  0.0000 ndimh = nmto+napw =    83   55   28
 bndfp: kpt     5 of    29 k=  0.1420 -0.0502  0.0000 ndimh = nmto+napw =    89   55   34
 bndfp: kpt     6 of    29 k=  0.0355  0.0251  0.0000 ndimh = nmto+napw =    78   55   23
 bndfp: kpt     7 of    29 k=  0.0710  0.0126  0.0000 ndimh = nmto+napw =    80   55   25
...
```
ここでハミルトニアンの次元がわかります。`[ham].pwmode = 11`（既定）のときは、APW のエネルギーのカットオフが `|q+G|^2 < pwemax`（`[ham].pwemax`、Ry）で、q に依存します（コードでは Rydberg 単位 $2m=\hbar=e^2/2=1$ を使います）。 

* 注意: ctrlgenToml.py の既定値では、すべての計算が次の設定です。  
   >空の球は置きません。 
   >EH=-1,EH=-2、MT 半径は接しないように -3% です。
   >RSMH=RSMH2=R/2

   この設定が最良だと主張するわけではありません。しかし、非線形の最適化はとても厄介です。一つずつ最適化することはしたくありません。
   この設定を変えるには、大がかりな系統的な研究が必要かもしれません。

#### job_tdos, job_pdos {#job-tdos-job-pdos}
全 DOS は次のようにして計算できます。
```
job_tdos gaas -np 8
gnuplot -persist tdos.gaas.glt
```
PDOS は次のようにします。
```
job_pdos gaas -np 8
```
`job_pdos` は対称性を使わないので、少し重い計算です。終わったら、`gnuplot -p pdos.site001.gaas.glt` で、実の球面調和関数の $lm$ のチャネルに分けた PDOS が表示されます
（チャネルは `job_pdos` の最後に表示されます）。
次の図は mp-2534 の As のサイトでの PDOS です。
![PDOS at the As site of mp-2534](image-7.png)
k 点の数に注意してください。正確な DOS / PDOS には、メッシュを
大きくする必要があるかもしれません。コマンドラインで、TOML のパスの書き方で
渡してください。
`job_tdos gaas -np 8 --ctrlg:bz.nkabc=[10,10,10]`
`grep '^bz.nkabc\|^nkabc' llmf_tdos` で確かめ、`save.gaas` も見てください。

### Step 4. バンドの図のための k の経路とブリルアンゾーンを作る {#step-4-create-k-path-and-brillouin-zone-for-band-plot}
`lmf` が収束したら、Step 5 の `job_band` でバンドの図を描きます。バンドの計算が正常かどうかは、バンドの図で確かめられます（磁性の系では、全体の磁気モーメントとサイトごとの磁気モーメントも確かめてください）。

`job_band` の前に、対称線のファイル [`syml.<sname>`](./syml) を `getsyml` で作ります（`gaas` は `ctrlg.<sname>.toml` の `<sname>` です）。
```
getsyml gaas
```
足りないパッケージは pip で入れてください。Togo さんの spglib と seekpath を使っています。終わったら BZ.html を見てください。BZ の中の k の経路が表示されます。
plotly で書いた対話的な図なので、座標の値を読み取れます。
getsyml による BZ.html の例は https://ecalj.sakura.ne.jp/BZgetsyml/ で見られます（ba2pdo2cl2 の図は[下の事例](#a-case-study-of-ba2pdo2cl2)にあります）。


### Step 5. バンドの図 {#step-5-band-plot}
```
job_band gaas -np 8
```

gnuplot のスクリプト `bandplot.isp1.glt`（スピンが二つなら `bandplot.isp2.glt` も）が作られます。必要なら編集して `gnuplot -persist bandplot.isp1.glt` で描き直してください。`job_band` の前に `syml.gaas` を編集すれば、対称線とメッシュの大きさを調整できます。
対称点の名前は BZ.html と `syml.<sname>` で確かめられます。0 eV がフェルミエネルギーです。

* `job_band` の最後には lmf のための上書きを付けられます。例えば `--ctrlg:ham.so=1 --ctrlg:ham.nspin=2` で、SOC を摂動として入れたバンドになります。

* ecalj では、`lmf` の k 点のメッシュ（`ctrlg.<sname>.toml` の `[bz].nkabc`）と GW の k 点のメッシュ（`[gw].n1n2n3`）は別にできます。前者は LDA の実時間にあまり効きませんが、後者は大きく効きます（そのため `[gw].n1n2n3` は小さくしたい）。

* ecalj のバンドの図のモードでは、BZ の端で対称性のために理論上縮退するバンドが縮退しません。APW の基底関数の数が限られているためです。そこで、pwemax=4 など（`job_band gaas -np 8 --ctrlg:ham.pwemax=4`）でバンドの図を描いてください（一時的な解決策です。自動化したいと考えています）。


#### job_tdos, job_fermisurface, job_pdos {#job-tdos-job-fermisurface-job-pdos}
`job_pdos` は PDOS を、`job_tdos` は全 DOS を計算し、`job_fermisurface` は Xcrysden でフェルミ面を描きます。
`job_fermisurface` は、Si の CBM の底の形を楕円体として描くのにも使えます。


### Step 6. QSGW の計算 {#step-6-qsgw-calculation}

ここで QSGW の計算を走らせます。QSGW は計算がとても重いので、
まず小さな系から走らせることを勧めます。

GW の設定は、同じ `ctrlg.<sname>.toml` の `[gw]`、`[mlo]`、`[blocks]`、`[product_basis]` の区分にあり、`ctrlgenToml.py`（Step 2）が既定値を書きます。
古い入力からの変換や、GW の区分をあとから足す方法は [入力ファイル](./toml_migration) を見てください。

`[gw]` には、ケルビン単位の二つの温度があります。`t_sigmaw`（自己エネルギーの中の準位の
フェルミ・ディラックの幅）と `t_tetrakbt`（$\chi_0$ をどう
ぼかすか。**必須**で、ないとプログラムは止まります）です。生成されたファイルでは
どちらも 300 K で、絶縁体・半導体ではそのままでかまいません（[kBT](./kBT) §7）。金属は [kBT](./kBT) を見てください。
実例の `ctrlg.gaas.toml` は古い入力から変換したもので、`t_tetrakbt = 0`（T=0 の四面体法）、`t_sigmaw = 262` です。

QSGW を始める前に最もよく調整するキーは、GW の
k 点のメッシュ `[gw].n1n2n3` です。実時間を大きく左右するので、
`[bz].nkabc` より小さく（ふつうはその 1/2 か 2/3）しておきます。

```toml
[bz]
nkabc = [12, 12, 12]      # LDA k-mesh

[gw]
n1n2n3 = [6, 6, 6]        # GW k-mesh (smaller -> faster)
```

実例の `ctrlg.gaas.toml` は `nkabc = [8, 8, 8]`、`n1n2n3 = [4, 4, 4]` です。
Si では `n1n2n3 = [6, 6, 6]` で十分です。非磁性の半導体では、`[gw]` のほかのキーを
触る必要はほとんどありません。`[gw]` のキーは [GW の設定](./gwinput) で説明しています。
次は Ref.[3] から取った、GaAs のバンドギャップの収束の様子です。
![convGaAs](./gaasnk.png) 
この図から、~0.1eV の精度でよければ 4 4 4 でも悪くないと言えそうです。6 6 6 なら十分です。8 8 8 より大きいのは数値の確認のためです（実用上は得るものがありません）。


#### スクリプト gwsc による QSGW の計算の流れ {#flow-of-qsgw-calculation-with-the-script-gwsc}
QSGW の計算は gwsc で走らせます。半導体では、QSGW の反復を数回行えば十分で、最終の結果に十分近くなります。
QSGW で得るのはバンド構造（一体のハミルトニアン）で、全エネルギーはまだ得られません。

`QPU` のファイルには GW の計算の対角成分が入っています。
混合積基底（`Mixed Product basis`）は GW の計算の鍵となる技術です。
```
gwsc -np NP [-np2 NP2] [--gpu] [--prec=tf32|fp32|fp64] nloop extension
```
（`--gpu` は GPU のプログラムを使い、`--prec` はその精度を選びます。[gwsc](./gwsc) を見てください。
磁性体で上向きと下向きに同じ基底を使いたいときは、
`ctrlg.<sname>.toml` で `[ham] phispinsym = true` にするか、実行時に
`--ctrlg:ham.phispinsym=true` を渡してください。）

GaAs の実例では次のようにします（1 は QSGW の反復の回数です）。
```
gwsc -np 8 1 gaas
```

`gwsc` のコンソールの出力は次のようになります（Si の例です）。
```text
--- Start gwsc ---
option= []
gwsc: using ctrlg.si.toml (TOML mode)
### START gwsc: ITERADD= 1, MPI size=  2, 2 TARGET= si
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
OK! ==== All calclation finished for  gwsc ====
```
... 

コンソールの出力はログのファイル `l*` に書かれます。`lsxC` は内殻による交換の自己エネルギーです。`lgw` はプログラム `hgw` のもので、`hgw` は価電子の交換の自己エネルギー、遮蔽されたクーロン相互作用 $W$、相関の自己エネルギーを一度に計算します（$W$ はメモリに保ちます）。`lvcc` はクーロン行列のものです。
この計算では `gwsc -np 2 1 si` を走らせています。1 は QSGW の反復の回数です。
gwsc を繰り返すと、前の計算の上に QSGW の反復が追加されます。

##### La2CuO4 の場合 {#a-case-of-la2cuo4}
La2CuO4 では次のようになりました（2025-06 のログをそのまま載せます。下に記録された古い `-vssig=0.8`
の形では、今はプログラムが止まります。同じ上書きは
`--ctrlg:ham.scaledsigma=0.8` です。また、三つの段 `hsfp0_sc --job=1`、`hrcxq`、
`hsfp0_sc --job=2` は、今は一つの段 `hgw --jobgw=1` です）。
```
2025-06-27 19:09:01.465241   mpirun -np 1 echo --- Start gwsc ---
--- Start gwsc ---
option=   -vssig=0.8 
### START gwsc: ITERADD= 5, MPI size=  32, 32 TARGET= lcuo
===== Ititial band structure ====== 
---> We use existing sigm file
0:00:00.041902   mpirun -np 32 /home/takao/bin/lmf lcuo   -vssig=0.8  >llmf_start
We found QPU.5 -->start to generate QPU.6...
===== QSGW iteration start iter 6 ===
0:00:18.111042   mpirun -np 1 /home/takao/bin/lmf lcuo   -vssig=0.8  --jobgw=0 >llmfgw00
0:00:18.233440   mpirun -np 1 /home/takao/bin/qg4gw   -vssig=0.8 --job=1 > lqg4gw
0:00:18.488197   mpirun -np 32 /home/takao/bin/lmf lcuo   -vssig=0.8  --jobgw=1 >llmfgw01
0:00:40.910973   mpirun -np 1 /home/takao/bin/heftet --job=1   -vssig=0.8 > leftet
0:00:41.052760   mpirun -np 1 /home/takao/bin/hbasfp0 --job=3   -vssig=0.8 >lbasC
0:00:43.290214   mpirun -np 32 /home/takao/bin/hvccfp0 --job=3   -vssig=0.8 > lvccC
0:01:00.327823   mpirun -np 32 /home/takao/bin/hsfp0_sc --job=3   -vssig=0.8 >lsxC
0:01:45.547149   mpirun -np 1 /home/takao/bin/hbasfp0 --job=0   -vssig=0.8 > lbas
0:01:46.806858   mpirun -np 32 /home/takao/bin/hvccfp0 --job=0   -vssig=0.8 > lvcc
0:01:59.034735   mpirun -np 32 /home/takao/bin/hsfp0_sc --job=1   -vssig=0.8 >lsx
0:02:45.526614   mpirun -np 32 /home/takao/bin/hrcxq  -vssig=0.8 > lrcxq
0:06:51.996781   mpirun -np 32 /home/takao/bin/hsfp0_sc --job=2   -vssig=0.8 > lsc
0:23:31.022636   mpirun -np 1 /home/takao/bin/hqpe_sc   -vssig=0.8 > lqpe
0:23:33.062303   mpirun -np 32 /home/takao/bin/lmf lcuo   -vssig=0.8  >llmf
===== QSGW iteration end   iter 6 ======== QSGW iteration start iter 7 ===
0:24:12.969096   mpirun -np 1 /home/takao/bin/lmf lcuo   -vssig=0.8  --jobgw=0 >llmfgw00
...
```
これは GPU なしです。QSGW の 1 回の反復に 24 分かかっています（各行の先頭に開始の時刻が表示されています）。
5 回目の QSGW の反復がすでに終わっていた（QPU.5run があることで確かめます）ので、6 回目の反復から始まっています。


##### ba2pdo2cl2 の事例 {#a-case-study-of-ba2pdo2cl2}

ba2pdo2cl2 の計算の記録です（下の時間、図、数値は記録したときのもので、今のコードで走らせ直したものではありません）。
POSCAR は次のとおりです。これを `POSCAR.ba2pdo2cl2` に書きます。
```
POSCAR_ba2pdo2cl2
1.0
-2.06443 2.06443 8.40383 
2.06443 -2.06443 8.40383 
2.06443 2.06443 -8.40383 
Ba Pd O Cl 
2 1 2 2 
Cartesian
0.0 0.0 6.5153213224
0.0 0.0 10.2923386776
0.0 0.0 0.0
0.0 2.06443 0.0
2.06443 0.0 0.0
0.0 0.0 3.1625293056
0.0 0.0 13.6451306944
```

Step 1-5 と同じ流れで LDA のバンドまで進みます。
```bash
vasp2ctrl POSCAR.ba2pdo2cl2
mv ctrls.POSCAR.ba2pdo2cl2.vasp2ctrl ctrls.ba2pdo2cl2
ctrlgenToml.py ba2pdo2cl2
lmfa ba2pdo2cl2
mpirun -np 8 lmf ba2pdo2cl2 | tee llmf
getsyml ba2pdo2cl2
job_band ba2pdo2cl2 -np 8
```

* getsyml が作る BZ.html は次のとおりです。
![image.png](https://qiita-image-store.s3.ap-northeast-1.amazonaws.com/0/482779/6af0e96e-c573-ede7-bee2-821134910c27.png)

* 次の図は、既定の計算での LDA のバンドです（対称点の名前は BZ.html で確かめられます。`syml.ba2pdo2cl2` も見てください）。0 eV がフェルミエネルギーです。金属なので、バンドギャップはありません。
![image.png](https://qiita-image-store.s3.ap-northeast-1.amazonaws.com/0/482779/8169cff0-78e6-ea50-0084-66a4cefc4cda.png)

QSGW を走らせる前にふつう見る必要があるキーは、
`t_tetrakbt`（上を見てください）と GW の k 点のメッシュです。

```toml
[gw]
n1n2n3 = [4, 4, 4]   # smaller than [bz].nkabc; 1/2 or 2/3 of it works
```

非磁性の半導体では、`[gw]` のほかのキーに
注意を払う必要はほとんどありません。

次のようにして QSGW の計算を走らせます。 
```
gwsc -np 32 1 ba2pdo2cl2
```
ここで 1 は QSGW の反復の回数です。QSGW の反復にはかなり時間がかかります。`gwsc` が表示するヘルプは最小限です（オプションは [gwsc](./gwsc) を見てください）。

反復の結果は `rst.<sname>`（電子密度）と `sigm.*`（vxcqsgw）に保たれます。
（最初からやり直したいときは、これらのファイルに加えて *run のファイルやディレクトリも消してください。）

* QSGW の 1 回の反復に 53 分かかります。

* job_band ba2pdo2cl2 -np 32 で次の図が得られます。QSGW の 1 回目で、Ef の近くのバンド構造が LDA のものから変わります。ただしまだ金属で、バンドギャップはありません。
![image.png](https://qiita-image-store.s3.ap-northeast-1.amazonaws.com/0/482779/0ddfb2d0-e63c-1cdb-ea16-d73a972565b1.png)

* QSGW の反復を続けるには、次を走らせます。
```
gwsc -np 32 nx ba2pdo2cl2
```
すでに 1 回行っているので、1+nx 回の QSGW の反復の結果が得られます。 

* ba2pdo2cl2 で 8 回反復すると、バンドギャップは 2.1 eV になりました。4 回目の反復のあとでバンドギャップが現れました。
![image.png](https://qiita-image-store.s3.ap-northeast-1.amazonaws.com/0/482779/4aebc98b-a050-e198-1769-d2d4236f0bf8.png)
* 実験と比べるには、QSGW80 を使うことを勧めます
（`ctrlg.<sname>.toml` で `[ham].scaledsigma = 0.8` と設定します）。

* スピン軌道相互作用。`sigm.<sname>` を得たあと、`[ham].so = 1`（L·S の方式。`[ham].nspin = 2` が要ります）にして lmf を走らせます。これで SOC の効果を入れられます。`so = 1`
は gwsc のサイクル全体には実装していないので、最後の段でだけ SOC を入れる必要があります（VxcQSGW を固定したあとで SOC を入れます）。

* 次を走らせると 
```
gwsc -np 32 5 ba2pdo2cl2 --ctrlg:ham.scaledsigma=0.8
```
lmf の計算で `ctrlg.ba2pdo2cl2.toml` の `[ham] scaledsigma` が上書きされます（save.ba2pdo2cl2 と `llmf` の先頭で確かめてください）。

##### KTaO3 の場合 {#a-case-of-ktao3}
KTaO3（ペロブスカイト、mp-3614）の QSGW の例
![image.png](https://qiita-image-store.s3.ap-northeast-1.amazonaws.com/0/482779/2e785d2a-9418-acc3-93f9-788daa2bd19b.png)


### lmchk {#lmchk}
  ```
  lmchk gaas
  ```
  で結晶の対称性を確かめます。さらに MT 半径を決め、MT の重なりを確かめます。既定の設定は -3% の重なり（重なりなし）です。

  - 対称性
  - MT の重なり
磁性の系で、格子の対称性より低い対称性をもつときは、
結晶の対称性を手で設定する必要があります。

これは `ctrlg.<sname>.toml` の `symgrp` に（`find` の代わりに）空間群の対称操作の生成元を書くことでできます。
SOC の場合はこの点に注意が必要です。

#### 空間群の操作の書き方 {#how-to-write-space-group-operation}
例えば r3x は x 軸のまわりの 3 回軸を意味します。ecalj での空間群の操作の書き方は [lmf.md の SYMGRP の節](lmf.md#symgrp) で説明しています。


### 計算をやり直す方法 {#how-to-start-over-calcualtions}
rst* を消します（前の実行の混合のファイル `__mixm.<sname>` は、`lmf` が始めるときに捨てます）。
MT 半径を変えたときは、lmfa からやり直します（atm* のファイルを消します）。

- 収束していれば問題はありません。 
- Ef に 3d のスパゲッティのように絡み合ったバンドがあるときは、注意が必要です。


## Samples/MATERIALS: 65 物質の LDA と MLO の例 {#samples-materials-lda-and-mlo-samples-of-65-materials}
ecalj の `Samples/MATERIALS/` には物質ごとに一つのディレクトリがあります。62 の簡単な物質（Si、GaAs、ZnO、MgO、SrTiO3、NiO、EuO、Bi2Te3、...）
と、La2CuO4、16 原子の InAs/GaSb 超格子、BaTiO3 です。それぞれに `ctrls.<sname>` と `ctrlg.<sname>.toml` があります。
これらは LDA の計算と MLO の模型の例です。65 物質すべての LDA のバンドと MLO の模型を
2026-10-01 に計算して比べました（MLO の模型の選び方と結果: [mlo](./mlo) §9。物質の一覧: ecalj の
[Samples/MATERIALS/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/MATERIALS/README.md)）。試験ではありません。

例えば Si を走らせます（まずディレクトリをコピーします）。
  ```
  cp -r ~/ecalj/Samples/MATERIALS/Si ~/work/ && cd ~/work/Si
  mpirun -np 1 lmfa si
  mpirun -np 8 lmf si
  ```
（2026-10-01 までは、これらの入力は最上位の `MATERIALS/` にある `Materials.ctrls.database` から `jobmaterials.py` がその場で作っていました。
このディレクトリは今はありません。中身は ecalj の `MD/past_log.md` §9 に書いてあります。）

* 主な入力ファイルは
```ctrls.si,ctrlg.si.toml```
です。下の節を見てください。```rst.si``` にはセルフコンシステントな電子密度が入っています。反復は出力ファイル `save.si` で確かめてください。lmf のコンソールの出力は llmf にあります。コンソールの出力をすべて知る必要はありません。 

* QSGW の前に、LDA のレベルの計算に問題がないことを確かめておくのがよいです。確かめるにはバンドの図が便利です。
バンドの図には `syml.si` のような対称線が必要で、これは
次のようにして作れます。
    ```
    getsyml si
    ```
  次に、
  ```
  job_band si -np 8
  ```
  を走らせると、gnuplot でバンドの図が表示されます。


-----

GetStarted はここまでです。[UsageDetailed](./UsageDetailed.md) に進んでください。
