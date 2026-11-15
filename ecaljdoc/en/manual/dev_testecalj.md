# Tests (testecalj) {#dev-testecalj}

> English version of the page [manual/dev_testecalj](../../manual/dev_testecalj.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/dev_testecalj.md of commit c1c4537df (source blob f43b00823148).

The tests of ecalj are run with `testecalj`. The minimum tests that `InstallAll.py` runs at installation also use it.

## Where it is {#where}

- [`SRC/exec/testecalj`](https://github.com/tkotani/ecalj/blob/main/SRC/exec/testecalj) (the program),
  [`SRC/exec/comp.py`](https://github.com/tkotani/ecalj/blob/main/SRC/exec/comp.py) (functions that compare files) and `pylib/diffnum0.py`.
  Installation puts them in the directory of the binaries (BINDIR).
- The programs that `testecalj` runs (`lmfa`, `lmf`, `gwsc`, ...) are those in the same directory as `testecalj`.
- The installation tests are in [`Samples/TestInstall`](https://github.com/tkotani/ecalj/tree/main/Samples/TestInstall).
  Any directory with a `test.py` (MLOsamples, EPS, Magnon, ...) can be a test.

## Usage {#usage}

```bash
cd Samples/TestInstall
testecalj -np 8 --all                                              # all targets on CPU (26 targets, 66 checks since 2026-10-02)
testecalj -np 8 --gwall                                            # the GW targets only
testecalj -np 8 copt si_gwsc                                       # some targets
testecalj -np 8 -np2 2 --gwall --gpu [--mp] [--run-args=--prec=fp32]   # GW on the GPU (fp64 / tf32 / fp32)
cd Samples/MLOsamples; testecalj -np 8 Si666gwsc Cu NiO           # MLO samples (--all works only in TestInstall)
```

- Without targets, `--all` or `--gwall`, nothing runs. `testecalj --help` shows the options.
- The test of a target `foobar` runs in its copy `foobar_work/`, which `testecalj` removes and makes again before the run
  (not with `--checkonly`, which only compares the files that are there).
- After each target, `testecalj` prints the summary of all targets so far. Count the checks in **the last summary only** (the PASSED lines of the log are cumulative).

## Test groups of Samples {#sample-groups}

[`TOOLS/samples_tests.sh`](https://github.com/tkotani/ecalj/blob/main/TOOLS/samples_tests.sh) runs the groups of Samples in turn (inputs, TestInstall, EPS, PROCAR, MLOsamples,
MLOQSGW, AFsymmetry, BenchmarkTest, the heavy GW tests, Magnon, ...) and writes one summary.

```bash
TOOLS/samples_tests.sh -np 8                      # every group (CPU)
TOOLS/samples_tests.sh --gpu -np 8 install mlo    # some groups, GPU
```

- The directory of the binaries is set by the environment variable `BIN` (default `~/bin`), the place of the logs by `LOG`. So binaries built elsewhere can be tested without touching the installed ones.
- A target that stops on the way is counted as STOPPED.
- The time of each group is in table 2 of section 5 of [GetStarted with AI](./getstartedAI).
- Python modules needed by the tests: numpy, pandas, matplotlib, pymatgen, seekpath, spglib, scipy, plotly. gnuplot is not needed (without it, only `.glt` files are written).

## Cautions when comparing {#cautions}

- On the GPU with `--mp`, fewer files are compared (the strict comparison of QPU is skipped). The CPU comparison is stricter.
- On the GPU, the same input can give different last digits on GPU 0 and 1. Fix `CUDA_VISIBLE_DEVICES` when comparing.
- **The two-atom tests do not reveal asynchronous errors (races) on the GPU.** After changing the order of execution on the GPU, compare Σ with the previous version to the last digit on a real material (LiTi₂O₄ 6³ etc.).
- When comparing QPU, look at the dSEnoZ column (not the column of the starting eigenvalues).

## How it works and how to add a test {#how-it-works}

Each target directory holds the input (`ctrlg.<sname>.toml`), finished files to compare with (`out.foobar`, `EPS*`, `QPU`, `QPD`, `log.*`, ...),
and `test.py`, which gives the steps of the calculation and how the files are compared (the name `test.py` is fixed by `testecalj`). For example, `si_gwsc/test.py`:

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

1. Decide the commands to run.
2. List the steps of the calculation in `runprogs`.
3. Compare the files (`diffnum` compares the numbers in the lines that contain the given words; `dqpu` compares QPU).

To add a test, run the calculation once, keep the files to compare, and write a `test.py` like the other `Samples/TestInstall/*/test.py`.

To compare only numbers, `test2_check` of `comp.py` can be used:

```python
from comp import test2_check
testdir = "<path>/Samples/Magnon/Fe_mlo_magnon/"
workdir = "<path>/Samples/Magnon/Fe_mlo_magnon_work/"
dat = 'MagSuscep.syml001'
tall = test2_check(testdir+'/'+dat, workdir+'/'+dat, abs_tol=0.01)
```
