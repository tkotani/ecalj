# ecalj 過去ログ — 片付けたものに含まれていたノウハウと、役に立つかもしれない情報

リポジトリから外したもの（各計算機の `ecalj/trash/` に移した）の中にあったノウハウや経緯を、話題ごとにまとめる。
入口は [CLAUDE.md](CLAUDE.md)（中身は [ecaljclaude.md](ecaljclaude.md)）と [README.md](README.md)。
利用者向けの説明は ecaljdoc が軸なので、ここには開発の側の細かいことだけを書く。

外したファイルの中身は git の履歴から取り出せる（表 1 の「外す前のコミット」を使う）:

```bash
git show <コミット>:<パス>                 # 1 ファイルを見る
git checkout <コミット> -- <パス>          # 作業ツリーに戻す（戻したら追跡し直すかを決める）
```

`trash/` は `.gitignore` に入っていて、空にする（本当に消す）のはメンテナが決める。

---

## 1. ビルドとインストール

### 1.1 `InstallAll.py` の打ち込み（最上位の `g`・`gg`・`i`・`ii`・`n`、2025-10〜2026-02）

1 行だけのシェルの控えだった。いまの使い方は README.md と ecaljdoc の ForDevelopers §4。

| ファイル | 中身 |
| --- | --- |
| `g` | `./InstallAll.py --fc gfortran` |
| `gg` | `./InstallAll.py --fc gfortran --notest` |
| `i` | `./InstallAll.py --fc ifx` |
| `ii` | `./InstallAll.py --fc ifort` |
| `n` | `./InstallAll.py --fc nvfortran --gpu` |

### 1.2 nvfortran が `gaugm.f90` で内部エラーになった回避（`build_nvfortran.sh`、2026-03-25）

- 当時の nvfortran は `gaugm.f90` のコンパイルで内部エラー（ICE）を出した。回避は、4 つの版（無印・`-D__MP`・`-D__GPU`・両方）ごとに
  `mpifort -O1 -fpic -Mbackslash -cpp`（`-acc`・`-Mvect` なし）で `gaugm.f90` だけを手でコンパイルし、`gaugm.f90` のタイムスタンプを
  古くしてから `gmake` を続け、最後にタイムスタンプを戻す、というもの
- 2026-09 の nvfortran 26.1 では、この回避なしで `InstallAll.py` が通る（kt1・kr7 の試験）。`InstallAll.py` にも CMake にも回避は入っていない
- nvfortran の別の症状（ecaljclaude.md「ビルド環境」）: 毎回違うファイルで signal 11 になることがある。同じコマンドをやり直せば通る

### 1.3 ISSP のスーパーコンピュータでのインストール（`jobinstall_kugui.sh`・`jobinstall_ohtaka.sh`）

インストール（コンパイルと試験）をジョブとして投げるもの。フロントエンドでは MPI を走らせられないため。詳しい手順は ecaljdoc の `UsageISSP.md`。
ベンチマークの投入スクリプトは `Samples/BenchmarkTest/ISSP/`。

- kugui（システム C）: `#PBS -q i2cpu`、`select=1:ncpus=16:mpiprocs=8:ompthreads=1`。
  `ulimit -s unlimited; module purge; module load nvhpc-nompi/24.7 openmpi_nvhpc compiler-rt tbb mkl`、`./InstallAll.py --fc nvfortran --gpu`
  （GPU ノードなら `-q i1accs`、`ncpus=64:mpiprocs=64`）
- ohtaka（システム B）: `#SBATCH -p i8cpu -N 1 -n 8 -c 1 --ntasks-per-node=8`、`./InstallAll.py --fc ifort`

### 1.4 CMake 以前のビルドと古いソースの控え（`SRC/exec/BK`、〜2026-09）

- `Make.inc.{gfortran,gfortran_mpik,ifort,ifort_mpik}`、`Makefile`・`Makefile.issp`・`Makefile.mac`: Makefile でビルドしていた頃の設定
- `*.F`（`asars.F`、`hambls.F`、`m_rdctrl.F`、`potpus.F` など 26 本）: 固定形式の頃のソースの控え。今のソースは `SRC/subroutines/*.f90`
- `gwsc.py`・`gwsc_org`・`gwsc_sym`・`gwscempty`、`run_arg2`・`run_argloc`: `gwsc` と `run_arg` の古い版

### 1.5 SRC の中のビルドの生成物（2026-10-01 に外した）

`SRC/exec/build`・`SRC/exec/build_gf14`（CMake の作業ツリー）、`SRC/exec_gfortran`・`SRC/exec_gfortran-14`・`SRC/exec_ifx`（`libecaljF.so` の写し）の
1853 ファイル・78 MB が、2026-09-24 のコミット `a0c7a7300` で誤って追跡されていた（公開のリポジトリへの push の前に見つけた）。
`.gitignore` に `/SRC/exec/build*/`・`/SRC/exec_*/` を足した。ビルドの置き場は `SRC/build_<コンパイラ>`（既に無視されている）。

## 2. GPU と計算機

### 2.1 GPU のロック（`gpu_wait_run.sh`、2026-04-29）

- `/tmp/gpu.lock` を `flock` で取ってから引数のコマンドを実行する、1 本の GPU を順番に使うための包み
- いまは `SRC/exec/pylib/gpu_lock.py`（GPU ごとに `/tmp/ecalj_res/gpu<N>.lock`）を `pylib/run_cmd.py` が取る。`CUDA_VISIBLE_DEVICES` の GPU だけ、
  `ECALJ_GPU_LOCK=0` で切れる

## 3. リファクタの記録

### 3.1 `m_sxcf_sc` の状態の整理の事前調査（`.refactor_notes/`、2026-04-30）

`m_itq` のリファクタ（`5cebb798`）の後、streaming（Phase 1-C）をやり直すために、`m_sxcf_sc.f90` の module 変数を
どこまで安全に動かせるかを調べたもの。

- OpenACC の directive がかかる変数（`zsecall`、`wvi_upper`・`wvr_upper`、`idx_i`・`idx_j`）は動かすと危ない。
  directive の無い、相関の側だけで使う変数（`omega`、`keepwv`、`nt0p`・`nt0m`、`ifrcw`・`ifrcwi`、ストップウォッチ）は安全。
  交換と相関の両方で使う変数（`ekc`、`eq`、`ntqxx`、`wkkr`）は両方を同時に書き換える必要がある
- ループの添字（`kx`、`irot`、`ip`、`isp`）は、BLOCK の host scope を確かめたうえで、各サブルーチンのローカルにして安全
- 回帰の道具 `snapshot_baseline.sh`: 各 `*_gwsc_work` の `SECU`・`SEC2U`・`SEXU`・`SEX2U` の md5 と、`lcombined` の診断行
  （`ntq`、`nbandmx`、`ncount` など）を段階ごとに記録して比べる（`snapshots/pre_step2`〜`stepWB1`）
- 調査では derived type にまとめる段階（Step 2.1〜2.5）を提案したが、最終的には type を使わず module の singleton にした
  （WB.3、ecaljclaude.md「実例: WB.3 リファクタ」）

## 4. サンプルと試験

### 4.1 MLOsamples の古い試行（2026-10-01 に外した）

`Al2O3_Cr/test1`〜`test7`、`RuO2/temp`、`Si666gwsc/temp`、`SrTiO3/temp`、`Al2O3_Cr/ctrl.tmp`、`MLOsamples/GWinput.tmp`、各サンプルの
`bandplot.isp*.glt.bk`。旧形式の入力（`GWinput`）で MLO のパラメータを試したときの跡で、試験（`test.py`）からも README からも参照されていなかった。

### 4.2 `Samples/TestInstall/TESTunused`（cdte、er、gaslc、srtio3、tio2、2025-10）

`testecalj` の対象から外した試験の入力。どの試験の道具からも参照されていなかった。

---

## 表 1. 片付けたもの（trash に移したもの）

| 日 | もの | 元の場所 | 外す前のコミット | 過去ログの節 |
| --- | --- | --- | --- | --- |
| 2026-10-01 | Materials Project の API キーを含む設定の写し | `ecalj_auto/OUTPUT/*/config.ini` の `apikey` 行 | `290397b34` の前 | ecalj_auto/README.md |
| 2026-10-01 | 古い写し | `SRC/BK`、`SRC/execgfortran`、`SRC/execAHC` | `116254e1d` の前（`4f9332d98`） | — |
| 2026-10-01 | GW1500 の古いスクリプト・表 | `ecalj_auto/run_gw1500_addrun*.sh`、`jobgw1500.sh`、`gw1500_recheck_*.sh`、`qpu_change.py`、`gw1500_rerun_20260930.tsv` | `116254e1d` の前 | ecalj_auto/GW1500_status.md |
| 2026-10-01 | ビルドの生成物 | `SRC/exec/build`・`build_gf14`、`SRC/exec_gfortran`・`exec_gfortran-14`・`exec_ifx` | `3f0771f2f` の前 | §1.5 |
| 2026-10-01 | 打ち込み用・古いスクリプト | `g`・`gg`・`i`・`ii`・`n`、`build_nvfortran.sh`、`gpu_wait_run.sh`、`jobinstall_kugui.sh`・`jobinstall_ohtaka.sh` | `3f0771f2f` | §1.1〜1.3、§2.1 |
| 2026-10-01 | 古いビルドとソースの控え | `SRC/exec/BK` | `3f0771f2f` | §1.4 |
| 2026-10-01 | リファクタの調査と回帰の控え | `.refactor_notes/` | `3f0771f2f` | §3.1 |
| 2026-10-01 | サンプルの古い試行と控え | 表の下の注 | `3f0771f2f` | §4.1、§4.2 |

注: サンプルの古い試行と控えは、MLOsamples の `test*`・`temp`・`*.bk`・`*.tmp`（§4.1）、`Samples/TestInstall/TESTunused`、
`Samples/TestInstall/eras/occnum.eras.bk`、`TOOLS/FparserTools/f_calltree.py.bk*`、`TOOLS/SrcFragments/f_calltree.py.bk*`・`ANALYZEnotusednow/analyze_temp~`、
`ecalj_auto/INPUT/testSGA/joblist.bk`。
