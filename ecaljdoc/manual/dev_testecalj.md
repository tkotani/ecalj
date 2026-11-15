# 試験の仕組み（testecalj） {#dev-testecalj}

ecalj の試験は `testecalj` で回します。`InstallAll.py` がインストールのときに最小の試験を回すのも、これです。

## 置き場 {#where}

- [`SRC/exec/testecalj`](https://github.com/tkotani/ecalj/blob/main/SRC/exec/testecalj)（本体）と
  [`SRC/exec/comp.py`](https://github.com/tkotani/ecalj/blob/main/SRC/exec/comp.py)（ファイルを比べる関数）、`pylib/diffnum0.py`。
  インストールすると、バイナリと同じディレクトリ（BINDIR）に入ります。
- `testecalj` が使うプログラム（`lmfa`、`lmf`、`gwsc` など）は、`testecalj` と同じディレクトリにあるものです。
- インストールの試験は [`Samples/TestInstall`](https://github.com/tkotani/ecalj/tree/main/Samples/TestInstall) にあります。
  ほかのサンプル（MLOsamples、EPS、Magnon など）も、`test.py` のあるディレクトリならどこでも試験にできます。

## 使い方 {#usage}

```bash
cd Samples/TestInstall
testecalj -np 8 --all                                              # all targets on CPU (26 targets, 66 checks since 2026-10-02)
testecalj -np 8 --gwall                                            # the GW targets only
testecalj -np 8 copt si_gwsc                                       # some targets
testecalj -np 8 -np2 2 --gwall --gpu [--mp] [--run-args=--prec=fp32]   # GW on the GPU (fp64 / tf32 / fp32)
cd Samples/MLOsamples; testecalj -np 8 Si666gwsc Cu NiO           # MLO samples (--all works only in TestInstall)
```

- ターゲットも `--all` も `--gwall` も無ければ、何も走りません。`testecalj --help` で選択肢が出ます。
- ターゲット `foobar` の試験は、写し `foobar_work/` の中で走ります。`testecalj` は走らせる前に `foobar_work/` を消して作り直します
  （`--checkonly` のときは作り直さず、そこにあるファイルを比べるだけです）。
- `testecalj` は、ターゲットを一つ終えるたびに、それまでの全部の要約を出します。件数は**最後の要約だけ**で数えます（ログの PASSED の行は延べ数です）。

## Samples の組ごとの試験 {#sample-groups}

[`TOOLS/samples_tests.sh`](https://github.com/tkotani/ecalj/blob/main/TOOLS/samples_tests.sh) が、Samples の組（inputs、TestInstall、EPS、PROCAR、MLOsamples、
MLOQSGW、AFsymmetry、BenchmarkTest、重い GW、Magnon など）を順に回し、一つの要約を書きます。

```bash
TOOLS/samples_tests.sh -np 8                      # every group (CPU)
TOOLS/samples_tests.sh --gpu -np 8 install mlo    # some groups, GPU
```

- 使うバイナリのディレクトリは環境変数 `BIN`（既定は `~/bin`）、ログの置き場は `LOG` で変えられます。インストールしたものに触らずに、別に作ったバイナリで試験できます。
- 途中で止まったターゲットは STOPPED と数えます。
- 組ごとの時間の目安は [GetStarted with AI](./getstartedAI) の §5 の表 2 にあります。
- 試験に要る Python のモジュール: numpy、pandas、matplotlib、pymatgen、seekpath、spglib、scipy、plotly。gnuplot は無くてもかまいません（図は描かずに `.glt` だけを書きます）。

## 試験を比べるときの注意 {#cautions}

- GPU で `--mp` のときは、比べるファイルが少なくなります（QPU の厳密な比べ方を省きます）。CPU のほうが比べ方は厳しいです。
- GPU では、同じ入力でも GPU 0 と 1 で下の桁が違うことがあります。比べるときは `CUDA_VISIBLE_DEVICES` を固定します。
- **2 原子の試験では、GPU の非同期の誤り（競合）は出ません。** GPU での実行の順番を変えたときは、実際の物質（LiTi₂O₄ の 6³ など）で Σ を前の版と桁まで比べます。
- QPU を比べるときは dSEnoZ の列を見ます（出発の固有値の列ではありません）。

## 仕組みと試験の足し方 {#how-it-works}

各ターゲットのディレクトリには、入力（`ctrlg.<sname>.toml`）と、比べるための計算済みのファイル（`out.foobar`、`EPS*`、`QPU`、`QPD`、`log.*` など）と、
`test.py` を置きます。`test.py` は、計算の手順と、ファイルの比べ方を書くものです（`test.py` という名前は `testecalj` が決めています）。例えば `si_gwsc/test.py` は次のとおりです。

```python
from comp import runprogs,diffnum,dqpu
def test(args,bindir,testdir,workdir):
        gwsc0= bindir + f'/gwsc 0 -np {args.np} '
        tall=''
        runprogs([
                 "rm log.si QPU",
                 gwsc0+ " si"
        ])
        dfile="QPU"
        for outfile in dfile.split():  #this loop is for dfile="QPU QPD"
                tall+=dqpu(testdir+'/'+outfile, workdir+'/'+outfile)
        outfile='log.si'
        tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=3e-3,comparekeys=['fp evl'])
        return tall
```

1. 走らせるコマンドを決めます。
2. 計算の手順を `runprogs` に並べます。
3. ファイルを比べます（`diffnum` は、指定した文字列を含む行の数値を比べます。`dqpu` は QPU を比べます）。

新しい試験を足すときは、計算を一度走らせて比べるファイルを残し、ほかの `Samples/TestInstall/*/test.py` にならって `test.py` を書きます。

数値だけを比べるときは `comp.py` の `test2_check` が使えます。

```python
from comp import test2_check
testdir = "<path>/Samples/Magnon/Fe_mlo_magnon/"
workdir = "<path>/Samples/Magnon/Fe_mlo_magnon_work/"
dat = 'MagSuscep.syml001'
tall = test2_check(testdir+'/'+dat, workdir+'/'+dat, abs_tol=0.01)
```
