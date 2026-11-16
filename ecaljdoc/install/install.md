# ecalj のインストール {#ecalj-install}

> ⚠️ **TOML への移行（2026-05）** — Fortran のバイナリは `ctrlg.<sname>.toml` だけを読みます。古い `ctrl.<sname>` / `GWinput` は `Legacy2toml.py <sname>` で変換してください。手引きの全体と移行済みのサンプルの一覧は [入力ファイル](../manual/toml_migration) を見てください。

ecalj をインストールするには、ecalj が使うソフトウェアと Python のモジュールを入れる必要があります。
さらに、便利のためにいくつかのソフトリンクを作ります。ほとんどの手順はほぼ自動です。以下の手順に従ってください。



## 必要なソフトウェアとネイティブのツール {#_1-required-software-and-native-tools}

ecalj の中心部はモダンな Fortran（Fortran 90/2003 の書き方）で書かれています。Fortran のソースファイルは次の場所にあります。

```text
ecalj/SRC/main/*.f90
ecalj/SRC/subroutines/*.f90
```

GPU 対応（OpenACC/CuBLAS）の cpp スイッチは、下で述べるコンパイルのときに自動で入/切できます。

推奨または必要なネイティブのツールとライブラリです（Ubuntu/Debian の例を後に示します）。

* `git`（ecalj のリポジトリの clone と更新に使います）。履歴を目で見て調べるには `gitk` が便利です。
* Fortran コンパイラ: `gfortran`（GNU）、`ifort`（Intel）、`nvfortran`（NVIDIA）のうち、お使いの計算機で使えるものを選びます。`nvfortran` は GPU 用です。
* BLAS/LAPACK/FFT のライブラリ。使えるなら Intel MKL を推奨します。OpenBLAS も使えます。
* MPI の実装: Open MPI または MPICH。（Ubuntu 24.04 では Open MPI 4.x で動きました。）
* ビルドのツール: `cmake`、`make`、`bash`。
* 描画とユーティリティ: `gnuplot`（描画のスクリプトと試験に使います。必須ではありません）。

Debian/Ubuntu では、次のコマンドでこれらのツールを入れます。

```bash
sudo apt update
sudo apt install -y git build-essential cmake make gfortran libopenmpi-dev openmpi-bin \
  libblas-dev liblapack-dev libfftw3-dev libopenblas-dev gnuplot
```


注意とヒント:

* Intel MKL は Intel のコンパイラでも gfortran でもうまく動きます。
<!-- * For GPU builds use `nvfortran` and pass `--gpu` to the installer; CUDA toolkit and drivers must be present. (See 3.Install below) -->
* 最近は、こうしたツールをホームディレクトリの下に入れる必要があるかもしれません。


例: 私の場合（ThinkPad T14、Ubuntu 24.04）は次を使っています。

* Open MPI 4.1.6
* CMake 3.28.x
* make 4.3
* gfortran 13.x
* Intel MKL 2020.x（使える場合）


`git` は更新や手元の変更の追跡（例えば `ecalj/` で `git diff`）に役立ちます。



## Python と Python のモジュール（pyenv を使う） {#_2-python-and-python-modules-using-pyenv}

> [!TIP]
> ecalj で使う Python の版は `pyenv` で管理することを勧めます。システムの Python に手を付けずに済み、環境を再現しやすくなります。複数の環境を管理したいのでなければ `pyenv-virtualenv` は要りません。ecalj の環境が一つだけなら `pyenv` だけで足ります。
>
> **注意:** 最近の Ubuntu などのシステムでは、プロジェクトの Python の環境をシステムの Python と分けておくのが最もうまくいきます。

Python 3.9 以上が必要です。以下に、`pyenv` を入れ、シェルを設定し、Python の版を入れ、ecalj に必要な Python のパッケージを入れる推奨の手順を示します。


1. 前提となるものを入れます（Debian/Ubuntu の例）。

```bash
sudo apt update
sudo apt install -y build-essential curl libssl-dev zlib1g-dev libbz2-dev \
  libreadline-dev libsqlite3-dev wget llvm libncurses5-dev libncursesw5-dev \
  xz-utils tk-dev libffi-dev liblzma-dev python-openssl git
```

1. `pyenv` を入れます（推奨のインストール方法）。

```bash
curl https://pyenv.run | bash
```

1. シェルの起動ファイル（`~/.bashrc` または `~/.profile`）に `pyenv` の初期化を加えます。

```bash
# pyenv
export PATH="$HOME/.pyenv/bin:$PATH"
eval "$(pyenv init --path)"
eval "$(pyenv init -)"
```

そしてシェルを読み込み直します。

```bash
source ~/.bashrc
```

1. ecalj で使う Python の版を入れます。例: 最新の 3.11 を入れて local の版に設定します。

```bash
pyenv install 3.11.6    # choose a concrete, available patch version
pyenv local 3.11.6      # or pyenv global 3.11.6 for system-wide in that account
```

1. pip を更新し、必要なパッケージを入れます。

```bash
pip install --upgrade pip
pip install numpy pandas seekpath spglib pymatgen mp-api scipy plotly gnuplot cif2cell
```


注意:
* conda/anaconda を好むなら、そちらの手順に置き換えてもかまいません。ただし `pyenv` はシステムの Python に手を付けず、複数のプロジェクトの環境ともうまく共存します。
* システムによっては、Python のビルドに要るものをパッケージマネージャで入れるほうがよいかもしれません（上の例は Debian/Ubuntu 用です）。


必要に応じて、上のパッケージを並べた `requirements.txt` をプロジェクトに作っておくと、同じインストールを再現できます。

```text
numpy
pandas
seekpath
spglib
pymatgen
mp-api
scipy
plotly
gnuplot
cif2cell
```


そして次のように入れます。

```bash
pip install -r requirements.txt
```



## インストールとインストールの試験 {#_3-install-and-installtest}

**ISSP の ohtaka と kugui（GPU 機）では、この節を飛ばして [installISSP.md](./installISSP.md) を見てください。**


### ecalj のパッケージを取得する {#get-ecalj-package}

ecalj のパッケージは [https://github.com/tkotani/ecalj/](https://github.com/tkotani/ecalj/) にあります。
master ブランチをチェックアウトしてください。

ISSP の計算機については [https://github.com/msobt/ecaljdoc](https://github.com/msobt/ecaljdoc) を見てください。

ecalj を clone するには次のようにします。

```bash
git clone https://github.com/tkotani/ecalj.git
```


clone すると `ecalj/` というディレクトリができます。
ecalj の履歴は次のように見られます。

```bash
cd ecalj
gitk --all
```


### InstallAll.py {#installall-py}

インストールには `InstallAll.py [Options]` を使います。ヘルプは次で出ます。

```bash
./InstallAll.py --help
```


ヘルプの出力の例:

```text
usage: InstallAll [-h] [-np NP] [--clean] [--gpu] [--bindir BINDIR] --fc FC [--notest] [--verbose]

Install ecalj and run tests.

options:
  -h, --help       show this help message and exit
  -np NP           number of mpi cores for install test
                   default: 8
                   specify the number of MPI parallelization in test calculation
  --clean          Clean CMakeCache CMakeFiles before make
                   default: none
                   delete the cache files before compiling
  --gpu            nvfortran for GPU
                   default: none
                   compile the GPU and GPU-MP version
  --bindir BINDIR  ecalj binaries and scripts
  --fc FC          fortran compiler gfortran/ifort/ifx/nvfortran
  --notest         no test. only compile
  --verbose        verbose on for debug
```


`InstallAll.py` は ecalj のバイナリとスクリプトをすべて、`--bindir <path>/` で指定したディレクトリ（既定は `$HOME/bin/`）にコピーします。`<path>/` は好きに決めてください。`--fc` は必須です。
GPU では `--fc nvfortran` と `--gpu` を一緒に使います。



試験は `ecalj/Samples/TestInstall/` にあります。試験を走らせるのは
`ecalj/SRC/exec/` の下の `testecalj` スクリプトです（バイナリと一緒に
`BINDIR` に入ります）。


## ecalj をインストールする {#install-ecalj}

`InstallAll.py` は ecalj をコンパイルしてインストールし、そのあと
`ecalj/Samples/TestInstall/` のインストールの試験を走らせます（走らせるのは `testecalj`、比較の補助は
`SRC/exec/comp.py`）。
成功すると最後に 'OK! ALL PASSED!' と出ます。かかる時間は 3〜10 分です。

ecalj をインストールするには、まずヘルプを見てください。

```bash
./InstallAll.py --help
```

典型的なインストールのコマンド:

```bash
./InstallAll.py --fc gfortran
```

GPU 用のビルドでは次のようにします。

```bash
./InstallAll.py --fc nvfortran --gpu
```

一から始めるには `--clean` を加えます。

これで Fortran のソースをコンパイル、リンクし、すべてのプログラムとスクリプトを bin にコピーし、`ecalj/Samples/TestInstall/` の最小限のインストールの試験を走らせます（3〜10 分）。

*ヒント:* バイナリをどこに入れるかの詳細は `./InstallAll.py` を見てください（既定は `BINDIR=~/bin`）。`getsyml`、`viewvesta`、`vasp2ctrl` などのツールはソフトリンクになります。

すべての試験が通ると、次のような出力になります。

```text
PASSED! copt TEST 1 out.lmf.copt
PASSED! te TEST 1 out.lmf.te
...
PASSED! nio_gwsc /home/<you>/ecalj/Samples/TestInstall/nio_gwsc_work/QPU
PASSED! nio_gwsc /home/<you>/ecalj/Samples/TestInstall/nio_gwsc_work/log.nio
PASSED! fe_gwsc  /home/<you>/ecalj/Samples/TestInstall/fe_gwsc_work/QPU
PASSED! fe_gwsc  /home/<you>/ecalj/Samples/TestInstall/fe_gwsc_work/log.fe
PASSED! ni_crpa  /home/<you>/ecalj/Samples/TestInstall/ni_crpa_work/Screening_W-v.h
PASSED! srvo3_crpa /home/<you>/ecalj/Samples/TestInstall/srvo3_crpa_work/Screening_W-v.h
OK! ALL PASSED ===

real    4m49.712s
user    19m38.952s
sys     3m38.583s
```

`ecalj/Samples/TestInstall/` で、個々の試験を次のように走らせることもできます。

```bash
cd ecalj/Samples/TestInstall
testecalj si_gwsc -np 8
```

* 試験ごとにディレクトリがあります（例: `si_gwsc`）。
* `testecalj` は隣に `si_gwsc_work/` を作ってそこで計算を走らせ、
  結果を `ecalj/Samples/TestInstall/si_gwsc/` の参照データと比べます。
* 試験の手順は各ディレクトリの `test.py`（小さなファイル）にあります。
  比較の補助は `ecalj/SRC/exec/comp.py` にあります。自分の試験を
  加えるのは簡単です。

インストールに失敗したら、きれいな状態からやり直す必要があるかもしれません。`ecalj/SRC/exec/` の `CMakeFiles` と `CMakeCache.txt` を消してください。`SRC/` の下に残った古い `*.mod` ファイルが問題を起こすこともあります。たいていは `./InstallAll.py --clean` を走らせれば直ります。

並列の make（`make -j`）が失敗するときは、make を手で走らせてください。

```bash
cd ecalj/SRC/exec/
rm -rf CMakeFiles CMakeCache.txt
FC=nvfortran cmake .
make
```


インストールのあと、コマンドのパス（BINDIR）を設定します。例えば ecalj のバイナリをすべて `~/bin` に移したなら、`.bashrc` に次を加えます。

```bash
PATH="~/bin/:$PATH"
```




## アンインストール {#_4-uninstall}

[アンインストール](./uninstall) を見てください。コマンドは三つです: `python3 uninstall.py`、
新しいシェルを開く、そして `rm -rf <ecalj-checkout>`。



## VESTA と getsyml のインストール {#install-vesta-and-getsyml}

構造を見るには VESTA が便利です。
（例: 私は Ubuntu 24 に `VESTA-gtk3.tar.bz2`（ver. 3.5.8、2022 年 8 月 11 日のビルド、23.8MB）を入れました。）

`ecalj/StructureTool/` に `viewvesta` コマンドがあります。次を試してください。

```bash
viewvesta ctrlg.si.toml      # canonical (post-2026-05)
viewvesta ctrls.si           # legacy structure-only seed file
```

これでビューアで構造を確かめられます。`/StructureTool` には変換の道具
`vasp2ctrl` と `ctrl2vasp` があります。これらで POSCAR との間で構造を変換できます。
`ctrl2vasp` は `viewvesta` と同じ三つの形の入力を受け付けます。

バンドの図の対称線を作る `getsyml` は、`InstallAll.py` が bin に入れます（ソフトリンク）。
作られた `syml.*` ファイルは ecalj のバンドの図に使います。
`spglib` と `seekpath` があれば、ふつうは追加の手順は要りません。
詳しくは上流のメモ [ecalj/GetSyml/README.md](https://github.com/tkotani/ecalj/blob/main/GetSyml/README.md) を見てください。



## 補足のメモ {#additional-memo}

* `InstallAll.py` が終わると、すべてのバイナリとシェルスクリプトが `--bindir <path>`（既定は `~/bin/`）にコピーされています。

* 後片付け:

  何かがうまくいかないときは、`--clean` を付けて `InstallAll.py` をもう一度走らせてください。

  次のことが必要になるかもしれません。

  ```bash
  rm -f CMakefiles CMakeCache.txt
  ```

  これを `ecalj/SRC/exec/` で行います（さらに `SRC/` の下の `*.mod` ファイルをすべて消します）。ふつうは必要ありません。

* Fortran だけをコンパイルする:

  Fortran のソースだけをコンパイルするには、`ecalj/SRC/exec/` に移って次を走らせます。

  ```bash
  FC=gfortran cmake . -D CMAKE_BUILD_TYPE=Debug  # or FC=ifort / FC=nvfortran
  make
  ```

  `CMakeLists.txt` を見てみてもよいでしょう。きれいにビルドしたいときは、先に `CMakeCache.txt` と `CMakeFiles/` を消してください。

* コンパイラのバグ: コンパイラが `-O2` の最適化でうまくいかないことがあります。`CMakeLists.txt` に書いてあるように、ソースファイルに応じて低い最適化のオプション（`-O1` や `-O0`）にすると、たいていそうしたバグを避けられます。ふつう `-O3` は使いません。こうした条件付きのコンパイルの設定をしてもかまいません。`CMakeLists.txt` を見てください。

* ソースコードと make の仕組みは `SRC/` の下にあり、試験の例は
  `Samples/` に移しました。

  ```text
  ecalj/
  ├── SRC/
  │   ├── exec        : CMakeLists.txt, build scripts, testecalj driver, comp.py
  │   ├── main        : All main *.f90
  │   ├── subroutines : All subroutines *.f90
  │   └── external    : bundled third-party (toml-f, ...)
  └── Samples/
      ├── TestInstall : Root of install tests (testecalj --all)
      └── EPS / PROCAR / MLOsamples / ... : topical samples (the old Legacy/ was rebuilt and removed on 2026-09-30)
  ```

  Fortran のコードはすべて `SRC/main/` と `SRC/subroutines/` にあります。
  `SRC/exec/CMakeLists.txt` が cmake で Makefile を作ります。
  （昔の名残: `SRC/TestInstall_is_moved_to_under_ecaljSamples`
  が TestInstall の古い場所を示しています。）

* 試験の仕組み: インストールの試験は `ecalj/Samples/TestInstall/` にあります。

  各サブディレクトリに自分の `test.py`（その試験の手順）があります。
  共通の駆動は `ecalj/SRC/exec/testecalj`（`BINDIR` に
  `testecalj` として入ります）、共通の比較の補助は
  `ecalj/SRC/exec/comp.py` にあります。自分の試験を加えるのは難しくありません。

  * まず何かを計算します。
  * 入力と最小限の結果を `Samples/TestInstall/` の下（または `Samples/` の下のどこか）の
    新しいディレクトリに置きます。
  * 比べ方をそのディレクトリの `test.py` に書きます。

  インストールの試験を走らせるには次のようにします。

  ```bash
  cd ecalj/Samples/TestInstall
  testecalj --all -np 8                # all install tests
  testecalj --gwall -np 8               # GW-only subset
  testecalj si_gwsc nio_gwsc -np 8      # specific tests
  ```

* 注意: Ubuntu 22 では、gfortran + openmpi がときどき失敗するのを見ました。問題が出たら mpich を使ってください。（今の状況は不明です。）
