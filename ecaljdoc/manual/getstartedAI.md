# GetStarted with AI

このページを AI に読ませれば、ecalj をガイドしてくれます。人間がこのページを読まなくていいです。

::: tip このページの使い方
**1. AI のチャットで（インストール無し）。** Web ページを読める AI アシスタント（例えば Claude）を開き、次のように頼みます。

`https://ecalj.github.io/ecaljdoc/manual/getstartedAI を読んで、ecalj を日本語で案内してください。`

**2. 手元の計算機の Claude Code で（一緒に計算を回す。おすすめ）。**

```bash
git clone https://github.com/tkotani/ecalj
cd ecalj
claude        # then say: "ecaljdoc/manual/getstartedAI.md を読んで、ecalj を案内して"
```

AI がやりたいことを尋ね、ecalj を短く説明し、インストール、最初の計算（GaAs）、試験まで、あなたの言葉で案内します。
クラスタ（`qsub` / `sbatch`）では、Claude Code がジョブを投入し、走っているジョブを見張ります（ログを見る、エラーに対応する、失敗したら再開する）。
ecalj の文書と道具（`ctrlgenToml.py`、`Legacy2toml.py`、`testecalj`、`ecalj_auto` のスケジューラなど）も Claude Code と一緒に作っています。
:::

**AI へ**: このページの残りはあなたに向けて書いています。何をどの順に読み、人をどう案内するかを書きます。
このサイトとリポジトリには日本語だけ、英語だけのところがあります。説明するときは、その人の言葉に訳してください（このサイトには各ページの英語版もあります）。
（2026-10-03 に書いた。このページとコードが食い違うときは、コードとリポジトリの `Changes.md` が正しい。）

## ライセンス

ライセンス（AI に関する注意を含む）は [ライセンス](./license) のページにあります。最初に読み、従ってください。

## ecalj を 5 行で

1. 第一原理の電子状態計算のパッケージ。PMT 基底（補強平面波 + マフィンティン軌道）による DFT（LDA/GGA、LDA+U）と、
   **QSGW**（準粒子セルフコンシステント GW）。QSGW は、調整のパラメータ無しに、実験に近いバンド構造・バンドギャップ・有効質量を与える。
2. QSGW の上に: 誘電関数、スピンゆらぎとマグノン、自己エネルギーのスペクトル、衝突電離、MLO（マフィンティン軌道に基づく局在軌道。
   U、J、cRPA も）による模型ハミルトニアン。
3. 入力: 物質ごとに TOML のファイル一つ、`ctrlg.<sname>.toml`。POSCAR（VASP の形式）から `vasp2ctrl` と `ctrlgenToml.py` で作る。
   設定のほとんどは自動。
4. CPU（MPI）と NVIDIA の GPU（`--gpu`。GW の行列積の精度は `--prec=tf32|fp32|fp64` で選ぶ）で走る。
5. オープンソース（AGPLv3）。小谷岳生と共同研究者が開発し、その多くを Claude Code と一緒に進めている。

## 何をどの順に読むか

その人の質問に要る深さまでだけ読んでください。このサイトのページは名前で書きます。このページの隣の `./<名前>.md` として開きます。

**段階 1: このサイト（マニュアル）**

| 順 | ページ | なぜ |
| --- | --- | --- |
| 1 | [What's new](./whatsnew.md) | 最近の変更と大きな話題。何が今のものか |
| 2 | [メインドキュメント](./README_tutorial.md) | ecalj で何ができるか、QSGW の概観、そして **GetStarted**（GaAs で一歩ずつ: POSCAR → ctrlg → LDA → バンドの図 → QSGW） |
| 3 | [Install](../install/install.md) | `InstallAll.py` とインストールの試験 |
| 4 | [Samples](./samples.md) | どのサンプルが何を示すか（サンプルのディレクトリにはそれぞれ README がある） |
| 5 | [入力ファイル](./toml_migration.md) | 入力ファイル、上書き `--ctrlg:<path>=<value>`、古い `ctrl` + `GWinput` からの変更 |
| 6 | 話題ごと | [lmf](./lmf.md)（DFT、対称性）、[gwsc](./gwsc.md)（QSGW）、[GPU 版](./ecaljgpu.md)、[kBT](./kBT.md)（温度と均し）、[MLO](./mlo.md)（模型、U/J/cRPA、マグノン）、[MLO-gwsc](./mlo_gwsc.md)、[誘電関数](./optical.md)、[スペクトル](./spectrum.md)、[UsageDetailed](./UsageDetailed.md)（反強磁性、SOC、LDA+U）、[コマンドラインのオプション](./cmdopts.md) |

**段階 2: リポジトリ**（その人がクローンしたとき、または GitHub で）

- `README.md`（リンク 3 つ）、次に `CLAUDE.md`: コードを扱う AI の入口。規則と文書の地図を読み込む。
  `Changes.md`: すべての変更、新しいものが上（正本。[What's new](./whatsnew.md) はその要約）。
- `Samples/<dir>/README.md`: 各サンプルが何か、回し方、試験が何を比べるか。
- 開発の記録はリポジトリの `MD/` にある（このサイトには無い）。`MD/README.md` から始める: 読む順、話題ごとのジャンプ
  （有限温度、MLO、マグノン、GPU、試験）。その人が ecalj を開発するとき、なぜそうなっているかを尋ねたときだけ読む。
- コードの構造（主プログラム → module）: `MD/module_map.md`。

## 人の案内の仕方

その人の言葉で説明する。一度に一歩ずつ。コマンドを示し、その人に走らせて結果を知らせてもらう。

1. **まず尋ねる**（短い質問を 3 つ）:
   - 何をしたいか（ecalj を試す / ある物質のバンドギャップやバンド構造 / 誘電関数 / マグノン / 模型ハミルトニアン /
     有限温度 / 大規模な自動の計算 / コードの開発）
   - 背景: DFT は? GW は?（理論をどこまで説明するかを決める）
   - 計算機: OS、コアの数、NVIDIA の GPU の有無、コンパイラ（gfortran、ifx、nvfortran）。ecalj はもうインストールしてあるか。
2. **概観を数行で**（§2）、答えに合わせて。初めての人には、QSGW が DFT より何を得るか（調整のパラメータ無しのギャップ）と、
   何がかかるか（数秒ではなく数分〜数時間）を言う。
3. **目的で道を選び**（表 1）、一緒に進める。
4. **結果を確かめる**（その人と一緒に）: lmf が出すギャップ（`llmf` の `gap = ... eV`）、バンドの図のファイル、QSGW の反復の履歴。
   サンプルの README や既知の値と比べる。
5. **失敗したとき**: ログの最後（`llmf`、`lgwsc`、失敗した段の `l*` のファイル）を読み、いちばん近いサンプルの README と比べ、
   エラーをリポジトリの `Changes.md` で探す。

*表 1*. 目的 → どこから始めるか

| 目的 | 始め | 次に |
| --- | --- | --- |
| ecalj を試す | [Install](../install/install.md)、次に GaAs で [GetStarted](./README_tutorial.md) | GaAs の QSGW を `gwsc` で（[gwsc](./gwsc.md)） |
| ある物質のギャップ / バンド | その人の POSCAR で GetStarted | QSGW は `gwsc`。QSGW80 は `ctrlgenToml.py` に `--ssig=0.8` |
| GPU | [GPU 版](./ecaljgpu.md) | ギャップには `--prec=tf32` で速く十分に正確 |
| 誘電関数 | [誘電関数](./optical.md)、サンプル `Samples/EPS` | |
| マグノン、U/J、cRPA、模型 | [MLO](./mlo.md) §6、サンプル `Samples/Magnon/Fe_mlo_magnon` | |
| 有限温度 | [kBT](./kBT.md): キーは `t_tetrakbt` と `t_sigmaw` の 2 つだけ | サンプル `Samples/kBT/scanT` |
| 反強磁性、SOC、LDA+U | [UsageDetailed](./UsageDetailed.md)、[lmf](./lmf.md) | サンプル `Samples/AFsymmetry`、`Samples/SOC`、`Samples/LDAU` |
| 多くの物質を自動で | リポジトリの `ecalj_auto/`（そこの `README.md`） | |
| コードの開発 | リポジトリの `CLAUDE.md` と `MD/README.md` | |

## 試験

`InstallAll.py` のあと（最小の試験はこれが回す）、全部の試験は `testecalj` で、組ごとにはリポジトリの `TOOLS/samples_tests.sh` で回す。

```bash
cd ecalj/Samples/TestInstall && testecalj -np 8 --all          # the install tests (CPU)
ecalj/TOOLS/samples_tests.sh -np 8                               # every group of the samples (CPU)
ecalj/TOOLS/samples_tests.sh --gpu -np 8 install mlo heavy       # some groups, GPU
```

組ごとに一行 `<group> PASSED <n> checks passed, 0 failed <seconds> s` を出す。まとめはログのディレクトリの `summary.txt`。
件数は最後のまとめだけで数える（ログには前の行がくり返し出る）。

*表 2*. 試験の組ごとの時間（開発者の計算機、`-np 8`、通った run のうち最も速いもの。t14 と mic は CPU、kt1 と kr7 は `--gpu`。
t14 は 2026-10-03、ほかは 2026-09-30 の記録で、組の中身が今と違うものがある（例えば magnon はサンプルが 1 つでなく 4 つ）。✗: 通った run が無い。bench の参照が古い）

| 試験の組 | t14 gfortran（CPU、16 コア） | mic ifx（CPU、32 コア） | kt1（RTX 5090 ×2） | kr7（RTX 5090） |
| --- | --- | --- | --- | --- |
| inputs | 2 min | 56 s | 58 s | 56 s |
| install (TestInstall, all) | 6 min | 4 min | 4 min | 4 min |
| eps | 2 min | 15 s | 53 s | 12 s |
| procar | 47 s | 19 s | 22 s | 20 s |
| mlo | 11 min | 5 min | 4 min | 4 min |
| mloqsgw | 2 min | 54 s | 54 s | 47 s |
| afsym | 55 s | 12 s | 9 s | 7 s |
| affix | 2 min | — | — | — |
| FermiSurface | 9 s | 6 s | 4 s | 5 s |
| Doping | 22 s | 14 s | 9 s | 10 s |
| HomoGas | 5 s | 3 s | 3 s | 3 s |
| BoltzTraP | 6 s | 4 s | 3 s | 3 s |
| SLAB | 32 s | 18 s | 11 s | 12 s |
| SOC | 4 min | 2 min | 1 min | 1 min |
| LDAU | 2 min | 1 min | 40 s | 42 s |
| EffectiveMass | 2 min | 42 s | 33 s | 26 s |
| Relax | 4 min | 2 min | 2 min | 2 min |
| DOS | 57 s | 32 s | 22 s | 21 s |
| IIR | 45 s | 35 s | 21 s | 19 s |
| AtomDimer | 27 min | 11 min | 7 min | 7 min |
| kBT_scanT | 28 s | 15 s | 27 s | 16 s |
| magnon | 2 min | 11 min | 8 min | 8 min |
| heavy（大きな GW の試験 4 つ） | — | 25 min | 7 min | 14 min |
| bench（BenchmarkTest） | — | 25.9 h ✗ | — | 2.8 h ✗ |
| **bench と heavy を除いた合計** | **71 min** | **41 min** | **32 min** | **31 min** |

つまり、bench と heavy 以外は 30〜70 分。bench は GPU で 3 時間ほど足される。ほかのジョブで混んだ計算機では試験はずっと遅くなる。

## 説明するときに間違えないこと

- QSGW が与えるのはバンド構造とギャップ（固有値と波動関数）で、**全エネルギーではない**。
- 入力は `ctrlg.<sname>.toml`（2026-05 から）。古い `ctrl.<sname>` + `GWinput` は読まない。`Legacy2toml.py` で変換する。
- 温度と均しは **2 つのキーだけ**で決める。`t_tetrakbt`（χ0）と `t_sigmaw`（Σ）（[kBT](./kBT.md)）。
- `ctrlgenToml.py` は既定で `scaledsigma = 1.0`（完全な QSGW）を書く。QSGW80（`--ssig=0.8`）は QSGW の自己エネルギーを 80 % 混ぜ、
  QSGW のギャップの平均的な過大評価を補正する。
- lmf が出すギャップは k のメッシュの上で取ったもの。バンドの端がメッシュから外れるとき（層状や六方晶の物質、K 点）は、
  バンドの経路に沿ったギャップのほうがよい見積もり。
- 対称性は spglib で見つける（`symgrp = "find"`、既定）。`symgrp` に生成元を書くと対称性を下げられる。
- MLO は直交化（Löwdin）した局在軌道。模型、U/J、cRPA、マグノンはこれを使う（[MLO](./mlo.md) §6）。
- GPU の精度: `--prec=tf32` は fp32 に比べてギャップを 1 meV ほど変える。古いオプション `--mp`（2026-09 より前、すべての積を TF32）は
  ギャップを 0.1〜0.3 eV ずらすことがあった。

## クローンせずにリポジトリを見る

ecalj のファイルはすべて GitHub で読める。Web ページを取ってくる AI には、生のテキストが次にある。

```text
https://raw.githubusercontent.com/tkotani/ecalj/main/<path>        e.g. .../main/Changes.md, .../main/Samples/README.md
https://github.com/tkotani/ecalj/tree/main/<directory>             the listing of a directory
```

*表 3*. リポジトリの最上位

| パス | 何があるか |
| --- | --- |
| `README.md`、`CLAUDE.md`、`Changes.md` | 入口（リンク 3 つと公開の手順）、コードを扱う AI の入口、すべての変更（新しいものが上） |
| `InstallAll.py` | ビルドとインストール（`--fc` に `gfortran`、`ifx`、`nvfortran`。`--gpu`。`--bindir`） |
| `SRC/subroutines/` | Fortran のライブラリ `libecaljF.so`（module `m_*.f90`。プログラムは薄い入口） |
| `SRC/main/` | Fortran のプログラムの入口（表 4） |
| `SRC/exec/` | Python とシェルの駆動スクリプト（`gwsc`、`job_*`、`ctrlgenToml.py`、`testecalj` など）と `pylib/` |
| `Samples/` | 計算例と試験。話題ごとに一つのディレクトリ、それぞれ README 付き（`Samples/README.md` が一覧） |
| `ecaljdoc/` | このサイトのソース |
| `ecalj_auto/` | 多くの物質の自動の計算（GW1500）、データベースを作る `gw1500db_build.py` |
| `TOOLS/` | 保守のスクリプト: `samples_tests.sh`（組ごとの試験）、`publish_*.sh`、`doclinks.py`、`module_map.py` |
| `MD/` | 開発の記録（このサイトには無い。まず `MD/README.md`） |
| `GetSyml/`、`StructureTool/` | バンドの図の k の経路、構造の変換 |

## プログラムと駆動スクリプト

*表 4*. 何が何を走らせるか（`InstallAll.py` が bin のディレクトリに入れる）

| コマンド | すること |
| --- | --- |
| `vasp2ctrl POSCAR` | POSCAR → `ctrls.POSCAR.vasp2ctrl`（ecalj の形の構造。`ctrls.<sname>` に写す） |
| `ctrlgenToml.py <sname>` | `ctrls.<sname>` → 入力 `ctrlg.<sname>.toml`（基底、k のメッシュ、GW の設定。QSGW80 は `--ssig=0.8`） |
| `lmfa <sname>` | 自由原子（出発の密度と原子のファイル） |
| `lmf <sname>` | DFT のセルフコンシステント計算（`sigm` があれば QSGW の一体の部分も）。`lmchk <sname>` は構造を確かめる |
| `getsyml <sname>` | バンドの図の k の経路 |
| `job_band <sname> -np N` | バンドの図（ファイル `bnd*.spin*`、`bandplot.isp1.glt`）。`job_tdos`、`job_pdos`: 状態密度 |
| `gwsc <n> -np N <sname>` | QSGW を `n` 回反復: lmf → GW の入力 → W と Σ（`hgw`）→ `hqpe_sc`（QP のエネルギー、`sigm`）→ lmf |
| `gwscconv -np N <sname> --conv-tol 0.1` | ギャップ（金属では E_F 近くの準位）が変わらなくなるまで QSGW をくり返す |
| `job_eps <sname>` | 誘電関数（RPA） |
| `job_mlo <sname>`、`job_mlo_soc`、`job_mloW <sname> [--crpa]`、`job_mlo_magnon` | MLO の模型、SOC 付き、U/J/cRPA、マグノン |
| `testecalj [--all] [targets]` | 試験（`Samples/TestInstall` で、またはサンプルのディレクトリで） |
| `ctrlg_update.py`、`Legacy2toml.py` | 古い `ctrlg` を今の形に、古い `ctrl` + `GWinput` を変換 |

どのプログラムも `--help` でオプションを出す。[コマンドラインのオプション](./cmdopts.md) が一覧。

## 入力ファイルと出力ファイル

`ctrlg.<sname>.toml` の節は `[struc]`（格子）、`[[site]]`（原子）、`[[spec]]`（種: 基底、半径、LDA+U）、`[bz]`（k のメッシュ、金属かどうか）、
`[iter]`（混合）、`[ham]`（汎関数、スピン、SOC、`scaledsigma`、平面波）、`[options]`、`[gw]`（GW のメッシュとカットオフ、`t_tetrakbt`、`t_sigmaw`）、
`[mlo]`、最後に `[product_basis]`（自動で作る）。どのキーもコマンドラインで上書きできる: `lmf <sname> --ctrlg:bz.nkabc=[8,8,8]`
（[入力ファイル](./toml_migration.md)）。

*表 5*. 見る出力ファイル

| ファイル | 何か |
| --- | --- |
| `llmf`（gwsc の中の lmf の標準出力）、`log.<sname>` | 計算の経過。ギャップの行 `VBmax= ... CBmin= ... gap = ... eV`。`efermi.lmf` はフェルミ準位 |
| `rst.<sname>` | 密度（再開用） |
| `sigm`、`sigm.<sname>` | lmf が使う QSGW の自己エネルギー（静的、エルミート） |
| `QPU`、`QPD` | 最後の反復の QP のエネルギー（上向き、下向きスピン） |
| `bnd*.spin1` | `job_band` のバンドの線 |
| `dos.tot.<sname>*` | 全状態密度 |
| `QSGW.<n>run/` | 各反復のファイル（`gwsc`、`gwscconv`） |

## 用語

PMT: 基底。補強平面波とマフィンティン軌道。QSGW: 準粒子セルフコンシステント GW。Σ の静的でエルミートな部分が交換相関ポテンシャルに代わり、
セルフコンシステントになるまで反復する。QSGW80: QSGW の Σ を 80 %、LDA を 20 %（`scaledsigma = 0.8`）。MLO: マフィンティン軌道から作り、
直交化（Löwdin）した局在軌道。模型に使う。cRPA: 遮蔽された U のための制限 RPA。ESM: スラブのための有効遮蔽媒質。`<sname>`: 物質の名前で、
ファイルの拡張子（`ctrlg.gaas.toml` → `gaas`）。

## データと、うまくいかないときに見るところ

- QSGW のバンドのデータベース: [DOSnpSupplement](https://github.com/tkotani/DOSnpSupplement)。データベースは逐次更新されている（過去の版もそこに残る）。
  版、各値の条件、利用条件はそこの README を見る。**新しい計算と比べる前に、その値をどの版・どの条件の ecalj が作ったかを README で確かめる。**
  版と条件が違えば値も違いうる。
- 古い `ctrl.<sname>` や `GWinput`: もう読まない。`Legacy2toml.py` で変換する。
- 古い `ctrlg` で GW のプログラムが止まる（「`[gw] … is required`」、`SmearX0`、`nbas`）: `ctrlg_update.py ctrlg.<sname>.toml` を走らせる。
  `[gw]` の必須のキー（`t_tetrakbt`、`t_sigmaw`、`HistBin_dw`、`HistBin_ratio`、`QforEPSau`、`EMAXforGW`）が無ければ、前のコードの既定値で足すので、結果は変わらない。
  古い `esmr`・`t_sigmakbt` はもう読まない（`ctrlg_update.py` が `t_sigmaw` に換算する）。詳しくは [GW の設定](./gwinput)。
- README と数値が違う: まず k のメッシュ、`scaledsigma`、反復の回数を確かめる。次に、README を書いたあとの既定値の変更を `Changes.md` で探す。
- 質問と報告: GitHub のリポジトリ https://github.com/tkotani/ecalj
