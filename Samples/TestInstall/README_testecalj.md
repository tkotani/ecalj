# 2025-10-6 test system in python

The commands for CPU and GPU (`--gpu`, `--mp`, `--run-args=--prec=fp32`), and how the checks are counted, are in
the developer's guide: ecaljdoc [manual/ForDevelopers](https://ecalj.github.io/ecaljdoc/manual/ForDevelopers) (sections 4, 5 and 12).
The groups of sample tests (TestInstall, EPS, PROCAR, MLOsamples, ...) are run by `ecalj/TOOLS/samples_tests.sh`.

`testecalj` is installed in your ecalj binary directory BINDIR.

## Usage
`testecalj` uses `comp.py` and `pylib/diffnum0.py` internally; they are in your bin together with `testecalj`.
(for developers: we can use `ecalj/SRC/exec/testecalj`. Then we use binaries at  `ecalj/SRC/exec`.)

>testecalj [-np mpi_size] [list of tests]

Run `testecalj --help`


To run only copt and si_gwsc test with mpi_size=8, run
>testecalj -np 8 copt si_gwsc
at ecalj/Samples/TestInstall

To run all tests (26 targets, 66 checks since 2026-10-02: the cRPA tests compare 3 files each), or the GW tests only,
>testecalj -np 8 --all
>testecalj -np 8 --gwall

Without a list of tests, `--all` or `--gwall`, nothing runs.

The test runs in `si_gwsc_work`, a copy of `si_gwsc`; testecalj removes and recreates it before each run
(not with `--checkonly`, which only compares the files that are there).
testecalj prints the summary of all tests so far after each test: count the checks in the last summary only.

-----------
## How the testecalj work?
For each test, we need finished calculations. Keep files for comparison 
{such as `out.foobar, EPS*,QPU,QPD,log.*`} in the test directory such as si_gwsc/.
In the manner described in si_gwsc/test.py, we reproduce such calculations, and compare these filed. The name 'test.py' is special name used in ecalj/SRC/exec/testecalj.

## How to add your test
If you like to add new test, you do run calculations
and keep some files for comparison. Then you add your own test as in si_gwsc/test.py.


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
This is called from `testecalj`. As this shows, we need to 
1. Speficy commnands for test.
2. Write steps of computation in runprogs.
3. Comparison (diffnum is for numerical comparison for lines including 'fp' and 'eval') in this case.

See other Samples/TestInstall/*/test.py as examples.