# For developer

> ⚠️ Test directories ship `ctrlg.<sname>.toml` as the input (see [Samples/EPS/](https://github.com/tkotani/ecalj/tree/main/Samples/EPS), [Samples/PROCAR/](https://github.com/tkotani/ecalj/tree/main/Samples/PROCAR), [Samples/MLOsamples/](https://github.com/tkotani/ecalj/tree/main/Samples/MLOsamples), [Samples/TestInstall/](https://github.com/tkotani/ecalj/tree/main/Samples/TestInstall)). (The older examples of `Samples/Legacy/` were rebuilt as samples with `ctrlg.<sname>.toml` or removed on 2026-09-30.) See [TOML migration](./toml_migration) and [Samples](./samples). The present procedure of the tests (`testecalj --all` with 26 targets and 66 checks, `TOOLS/samples_tests.sh` for the groups of samples) is in [ForDevelopers](./ForDevelopers) §5.

## Test system `testecalj` (2025-10-8).

Our new test system is made from two files 
```
SRC/exec
├── comp.py    utitities (functions to show difference of files) called from testecalj
└── testecalj   main program for test
```

* Install test described in InstallAll.py is performed at `ecalj/Samples/TestInstall`.
* We can run testecalj as `>testecalj foobar`, where `foobar/` is the name of a test directory.
`testecalj` creates `foobar_work/` directory and do test in it. 
* `foobar/` should contain initial settings (`ctrlg.<sname>.toml`) and files to be compared. In addition, we have to write `test.py` which describe schedules to run programs and to compare files. `ecalj/Samples/TestInstall/foobar` contains samples of test directories.
* For any test directory `foobar`, we can perform the test by `testecalj foobar`. For example, `Fe_mlo_magnon` lives at `ecalj/Samples/Magnon/Fe_mlo_magnon` with its `test.py`.
* We can write your own `test.py` easily. 

* `testecalj` use binaries such as `lmfa,lmf,qg4gw...` in the directory containing the `testecalj`. 

* We can do only numerical comparison of files by
```
from comp import test2_check
testdir= "/home/takao/ecalj/Samples/Magnon/Fe_mlo_magnon/"
workdir= "/home/takao/ecalj/Samples/Magnon/Fe_mlo_magnon_work/"
dat='MagSuscep.syml001'
tall=test2_check(testdir+'/'+dat, workdir+'/'+dat,abs_tol=0.01)
```
where, we need to import `ecalj/SRC/exec/comp.py`.

---

* Samples/TestInstallにテストを移した。InstallAll.pyを実行するとテストではそこへ飛びます。個別実行もSamples/TestInstallで行います。どんなところにあるfoobar/に対してもその横にfoobar_work/をつくってそこでテストするようにしました。foobar/にはtest.pyをいれること
* testecaljというコマンドを作った。短いプログラム。これはバイナリのあるディレクトリBINDIRに入ります。使うバイナリはその同じディレクトリにあるものです。
test.py を置いたディレクトリならどこでもテストできる。
 
