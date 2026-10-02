
# ISSPのkuguiとohtaka(GPU)でのインストール

> 📝 **Note** — This is a memo to install ecalj to the machines in
> ISSP@u-tokyo as of 2025. Module names, software versions, and
> filesystem layouts on `kugui` / `ohtaka` change over time; if any
> step fails, check the latest ISSP documentation
> ([SCC のページ](https://www.issp.u-tokyo.ac.jp/supercom/)) and adjust
> module versions / paths accordingly. The configuration below has
> been verified to work for the author's typical workflow but is not
> automatically updated when ISSP's environment changes.

## kuguiでのインストール

### 1. ログイン

### 2. [pyenvなどを用いてpythonとtoolをインストール](./install.md)

### 3. `~/.bashrc`の書き換え

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

### 4. ecaljのインストールとテスト

```bash
git clone https://github.com/tkotani/ecalj.git
cd ecalj
```

そのあとインストール＋テストのコマンドを実行します。

```bash
qsub jobinstall_kugui.sh
```

これを実行すると、コンパイルのあとテスト計算が進みます。`jobinstall_kugui.sh` ではCPU版とGPU版の両方がインストールされます。その中で呼ばれる
InstallAll.pyは`$HOME/bin`ディレクトリを作り、そこへecaljのバイナリやスクリプトをコピーします。
もし困る場合は、現在のコピーのバックアップをとったうえでecaljに上書きさせる、あるいはInstallAll.pyの中で--bidir binを検索し書き直してください。binにパスを通しておく必要があります。

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

```bash
tail -f ecalj/Samples/TestInstall/summary.txt
```

で進行状況を確認できます。最後に`PASSED! srvo3_crpa/Screening_W-v_crpa.h`と表示され、FAILやエラーがなければOKです。
nvfortranでもテストはCPU実行です。


### 5. GPU計算テスト

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

で終了していればOKです。


### 6. 使い方メモ

[使い方メモ](../manual/UsageISSP.md)


## ohtakaでのインストール

ほぼ同様です。`~/.bashrc`の最後に以下を追加します。

```bash
ulimit -s unlimited
export PATH=$HOME/bin:$PATH
module purge
module load openmpi/4.1.5-oneapi-2023.0.0-classic
```

インストールは

```bash
sbatch jobinstall_ohtaka.sh
```

です。この段階で最後に`OK! ALL PASSED!`が見れない場合はインストールできていません。

それを確認後、計算テストは `ecalj/Samples/BenchmarkTest/inas2gasb2` を作業用のディレクトリに写し、そこに
`ecalj/Samples/BenchmarkTest/ISSP/job_ohtaka.sh` を置いて

```bash
sbatch job_ohtaka.sh
```

を行ってください。16 原子の `inas4gasb4` も同じ手順です（スクリプトの `id` を書き換える）。
スクリプトは 2026-09-30 に今の入力（`ctrlg.<sname>.toml`）とオプションに合わせたもので、それ以後 ISSP の機械では実行していません。


kuguiではmpirunのかわりにsrunがつかわれるなどのためgwscなどを走らせるときにこれが必要です。

## ISSPシステムでのカスタマイズについて

ecaljソースコードの `SRC/exec/MachineDependence.py` にマシン依存性を記述しています。
kuguiではmpirunのかわりにsrunが使われるなど、gwsc等を走らせるときにこれが必要です。
