
# ISSPのkugui（GPU）とohtaka（CPU）でのインストール

> 📝 **注意** — 東大物性研（ISSP）の計算機に ecalj をインストールするための 2025 年時点のメモです。
> `kugui` / `ohtaka` のモジュールの名前、ソフトウェアの版、ファイルシステムの構成は変わっていきます。
> うまくいかない段があれば、ISSP の最新の説明（[SCC のページ](https://www.issp.u-tokyo.ac.jp/supercom/)）を見て、
> モジュールの版やパスを合わせてください。下の設定は著者のふだんの使い方で動くことを確かめたものですが、
> ISSP の環境が変わっても自動では更新されません。

## kuguiでのインストール

### ログイン

### [pyenvなどを用いてpythonとtoolをインストール](./install.md)

### `~/.bashrc`の書き換え

`~/.bashrc`の最後に以下を挿入します。

```bash
ulimit -s unlimited
export PATH=$HOME/bin:$PATH
module purge
module load nvhpc-nompi/24.7 openmpi_nvhpc compiler-rt tbb mkl
if which nvidia-cuda-mps-control > /dev/null 2>&1 ; then
  export CUDA_MPS_PIPE_DIRECTORY=$(pwd)/nvidia-mps-$(hostname)
  export CUDA_MPS_LOG_DIRECTORY=$(pwd)/nvidia-log-$(hostname)
  echo "start nvidia-cuda-mps-control at" $(hostname)
  nvidia-cuda-mps-control -d
fi
```

> [!TIP]
> 必要な`module`のロードとMPS起動のスクリプト例です。
> OneAPI_MPIはやや速いですが安定性はやや劣ります。
> Intel MKLはnvfortran/gfortran両方で利用可能です。

ログインや mpirun のときに次のメッセージが出ることがありますが, 害はありません（§1.4）。

>--------------------------------------------------------------------------
>The library attempted to open the following supporting CUDA libraries,
>but each of them failed.  CUDA-aware support is disabled.
>libcuda.so.1: cannot open shared object file: No such file or directory
>libcuda.dylib: cannot open shared object file: No such file or directory
>/usr/lib64/libcuda.so.1: cannot open shared object file: No such file or directory
>/usr/lib64/libcuda.dylib: cannot open shared object file: No such file or directory
>If you are not interested in CUDA-aware support, then run with
>--mca opal_warn_on_missing_libcuda 0 to suppress this message.  If you are interested
>in CUDA-aware support, then try setting LD_LIBRARY_PATH to the location
>of libcuda.so.1 to get passed this issue.
>

### ecaljのインストールとテスト

```bash
git clone https://github.com/tkotani/ecalj.git
cd ecalj
```

そのあとインストール＋テストのコマンドを実行します。

`jobinstall_kugui.sh` はリポジトリから外した（2026-10-01）ので, 次の内容のジョブスクリプト（例えば `jobinstall_kugui.sh`）を `ecalj` に作って投入します。フロントエンドでは MPI の試験を走らせられないためです。

```bash
#!/bin/sh
#PBS -q i2cpu
#PBS -l select=1:ncpus=16:mpiprocs=8:ompthreads=1
#PBS -N install
ulimit -s unlimited
module purge
module load nvhpc-nompi/24.7 openmpi_nvhpc compiler-rt tbb mkl
export PATH=$HOME/bin:$PATH
./InstallAll.py --fc nvfortran --gpu
```

```bash
qsub jobinstall_kugui.sh
```

これを実行すると、コンパイルのあとテスト計算が進みます。`jobinstall_kugui.sh` ではCPU版とGPU版の両方がインストールされます。その中で呼ばれる
InstallAll.pyは`$HOME/bin`ディレクトリを作り、そこへecaljのバイナリやスクリプトをコピーします。
もし困る場合は、現在のコピーのバックアップをとったうえでecaljに上書きさせる、あるいはInstallAll.py に `--bindir <ディレクトリ>` を渡してください（既定は `$HOME/bin`）。そのディレクトリにパスを通しておく必要があります。

テストはOKを出しながら進みます。以下のようなエラーっぽいメッセージが出ても無視して大丈夫です（正常動作）。

```text
[cpu121:54969] 7 more processes have sent help message help-mpi-common-cuda.txt / dlopen failed
[cpu121:54969] Set MCA parameter "orte_base_help_aggregate" to 0 to see all help / error messages
--------------------------------------------------------------------------
The library attempted to open the following supporting CUDA libraries,
but each of them failed.  CUDA-aware support is disabled.
libcuda.so.1: cannot open shared object file: No such file or directory
libcuda.dylib: cannot open shared object file: No such file or directory
/usr/lib64/libcuda.so.1: cannot open shared object file: No such file or directory
/usr/lib64/libcuda.dylib: cannot open shared object file: No such file or directory
If you are not interested in CUDA-aware support, then run with
--mca opal_warn_on_missing_libcuda 0 to suppress this message.  If you are interested
in CUDA-aware support, then try setting LD_LIBRARY_PATH to the location
of libcuda.so.1 to get passed this issue.
```

テストの最後は

```text
PASSED! ni_crpa/Screening_W-v_crpa.h
PASSED! srvo3_crpa/Screening_W-v.h
PASSED! srvo3_crpa/Screening_W-v_crpa.h
OK! ALL PASSED ===
  See work/summary.txt
Elapsed time for make        : 51 seconds
Elapsed time for testecalj.py: 409 seconds
```

という感じで終了します。これでインストールとテスト完了です（合計10分程度）。

テスト中は

ジョブの標準出力（PBS では `install.o<ジョブ番号>`）で進行状況を確認できます。各ターゲットの結果は `ecalj/Samples/TestInstall/<ターゲット>_work/summary.txt` に書かれます（上の出力例の `See work/summary.txt` は古い版の表示）。最後に`PASSED! srvo3_crpa/Screening_W-v_crpa.h`と表示され、FAILやエラーがなければOKです。
nvfortranでもテストはCPU実行です。

### GPU計算テスト

`ecalj/Samples/BenchmarkTest/inas2gasb2` を作業用のディレクトリに写し（作業ファイルが大きい）、そこに
`ecalj/Samples/BenchmarkTest/ISSP/job_kugui.sh` を置いて

```bash
qsub job_kugui.sh
```

を実行。`cat lgwsc` してみて

```
===== QSGW iteration end   iter 2 ===
OK! ==== All calclation finished for  gwsc ====
```

で終了していればOKです。サンプルに前の結果 `QPU.1run` が入っているので、`gwsc 1` の 1 回の反復は 2 番と数えられ、結果は `QPU.2run` に書かれます。これをサンプルの `QPU.1run` と比べてください。スクリプトは `-c kugui_gpu`（`SRC/exec/clusters.toml` の GPU ノードの設定）を渡します。

### 本番の計算

本番の計算のキューは `F1accs` などです（`i1accs` は 30 分までの試験用）。例:

```bash job.sh
#!/bin/sh
#PBS -q F1accs
#PBS -l select=1:ncpus=64:mpiprocs=64:ompthreads=1
#
ulimit -s unlimited
id=inas6gasb6
gwsc 5 -np 64 -np2 4 --gpu -c kugui_gpu $id > lgwsc
```

- `-np` は CPU で走る部分の MPI の並列数、`-np2` は GPU で走る部分の並列数（ふつうは GPU の数）。最初の数（5）は反復の回数です。詳しくは [GPU 版](../manual/ecaljgpu)。
- GPU 版は小さい系ではあまり効きませんが、4 枚の GPU で 40 原子くらいまで扱えます（計算の条件によります。計算量は k 点の数の 2 乗で増えます）。

## ohtakaでのインストール

ほぼ同様です。`~/.bashrc`の最後に以下を追加します。

```bash
ulimit -s unlimited
export PATH=$HOME/bin:$PATH
module purge
module load openmpi/4.1.5-oneapi-2023.0.0-classic
```

`jobinstall_ohtaka.sh` もリポジトリから外したので, 次の内容のスクリプトを `ecalj` に作って投入します。

```bash
#!/bin/sh
#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -c 1
#SBATCH --job-name=install
#SBATCH --ntasks-per-node=8
./InstallAll.py --fc ifort
```

```bash
sbatch jobinstall_ohtaka.sh
```

です。この段階で最後に`OK! ALL PASSED!`が見れない場合はインストールできていません。

それを確認後、計算テストは `ecalj/Samples/BenchmarkTest/inas2gasb2` を作業用のディレクトリに写し、そこに
`ecalj/Samples/BenchmarkTest/ISSP/job_ohtaka.sh` を置いて

```bash
sbatch job_ohtaka.sh
```

を行ってください。ここでも結果は `QPU.2run` に書かれるので、サンプルの `QPU.1run` と比べます。スクリプトは `-c ohtaka`（srun で起動）を渡します。16 原子の `inas4gasb4` も同じ手順です（スクリプトの `id` を書き換える）。
スクリプトは 2026-09-30 に今の入力（`ctrlg.<sname>.toml`）とオプションに合わせたもので、それ以後 ISSP の機械では実行していません。

- `InstallAll.py` は中で `mpirun` を使いますが、1 ノードなので動きます。
- 大きな計算の並列の指定の例（16 ノード、512 プロセス）。メモリが足りなければ、ノードあたりのコアを減らして、コアあたりのメモリを増やします。

```bash
#SBATCH -p F16cpu
#SBATCH -N 16
#SBATCH -n 512
#SBATCH -c 4
#SBATCH --ntasks-per-node=32
#SBATCH -t 24:00:00
```

- 実行中に MKL のエラーが出たら `-c 1` にします。

## ISSPシステムでのカスタマイズについて

ecaljソースコードの `SRC/exec/clusters.toml`（`gwsc` などの `-c <名前>` で選ぶ）と `SRC/exec/pylib/MachineDependence.py` にマシン依存性を記述しています。
ohtakaではmpirunのかわりにsrunが使われ, kuguiのGPUノードでは `mpirun --bind-to none` を使うなど、gwsc等を走らせるときにこれが必要です。
