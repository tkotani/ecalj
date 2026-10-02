# Test system `testecalj`
 (2025-10-8).

Our new test system is made from two files 
```
SRC/exec
├── comp.py    utitities (functions to show difference of files) called from testecalj
└── testecalj   main program for test
```

* Install test described in InstallAll.py is performed at `ecalj/Samples/TestInstall`.
* We can run testecalj as `>testecalj foobar`, where `foobar/` is the name of a test directory.
`testecalj` creates `foobar_work/` directory and do test in it. 
* `foobar/` should contain initial settings such as `ctrl,GWinput` and files to be compared. In addition, we have to write `test.py` which describe schedules to run programs and to compare files. `ecalj/Samples/TestInstall/foobar` contains samples of test directories.
* For any test directories `foobar`, we can perform test by `testecalj foobar`. In fact, we can run test of `Fe_magnon` at `ecalj/Samples/Magnon` since we added `ecalj/Samples/Magnon/test.py`. 
* We can write your own `test.py` easily. 

* `testecalj` use binaries such as `lmfa,lmf,qg4gw...` in the directory containing the `testecalj`. 

* We can do only numerical comparison of files by
```
from comp import test2_check
testdir= "/home/takao/ecalj/Samples/Magnon/Ni_magnon/"
workdir= "/home/takao/ecalj/Samples/Magnon/Ni_magnon_work/"
dat='wan_ChiPMz.mat.syml1'
tall=test2_check(testdir+'/'+dat, workdir+'/'+dat,tol=0.01)
```
where, we need to import `ecalj/SRC/exec/comp.py`.

---

* Samples/TestInstallにテストを移した。InstallAll.pyを実行するとテストではそこへ飛びます。個別実行もSamples/TestInstallで行います。どんなところにあるfoobar/に対してもその横にfoobar_work/をつくってそこでテストするようにしました。foobar/にはtest.pyをいれること
* testecaljというコマンドを作った。51行の短いプログラム。これはバイナリのあるディレクトリBINDIRに入ります。使うバイナリはその同じディレクトリにあるものです。
たとえばSamples/AHC/xxxにでもtest.pyをおけばテストできるはず。
 
