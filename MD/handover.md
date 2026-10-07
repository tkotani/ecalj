# handover.md — 引き継ぎ（Claude の個人メモリから、今も有効なものを写したもの）

2026-10-01 に、t14 の Claude の個人メモリ（`~/.claude/projects/-home-takao-ecalj/memory/`、パッケージに入らない）から、別の機械や
記憶の無いセッションに引き継ぐ価値があり、今も正しいものだけを写した（user「メモリの内容はパッケージに入らないので MD/ に。ただし混乱を招くものは良くない」）。
その時点の状況（push の状況、夜間作業の途中経過）、古くなったもの、研究ログや TODO と重なるものは写していない。各項目の終わりの（ ）は確かめた日と元のメモの名前。

読む順: [../CLAUDE.md](../CLAUDE.md) → [ecaljclaude.md](ecaljclaude.md)（方針、記録の方針） → このファイル → [TODOandQuestion.md](TODOandQuestion.md)（未決と実行中） →
[research_kotani_log.md](research_kotani_log.md) の先頭（最新の経過）。片付けたもののノウハウは [past_log.md](past_log.md)。人向けの手引きは ecaljdoc（`manual/ForDevelopers.md` など）。

---

## 1. user との取り決め

ecaljclaude.md「記録の方針」「Claude の作業の進め方」にあるもの（コメント・コミット・研究ログ・文書・図の書き方、着手前に計画を言う、エージェントの指摘は確かめる）は
そちらが正本。ここにはそれ以外を書く。

- **push**: commit はその場でよい。push は dev（`git@github.com:tkotani/ecaljdeveloper.git`、`main`、private）も含めて user の指示を待つ。
  rel（`git@github.com:tkotani/ecalj.git`、**`main`**、公開）は、dev で寝かせてから、3 環境（t14 gfortran、kt1 nvfortran GPU、ifx の mic か ucgw）の試験を同じコミットで通してから。
  早送りだけ（`git merge-base --is-ancestor rel/main main`）、force-push しない。壊れた rel を直す緊急の修正だけは例外（2026-10-01 確認: rel の既定ブランチは `main`、`master` は無い。
  feedback_push_timing・feedback_rel_push_policy）
- **待つときに `sleep` を書かない**: 完了のファイルやログの行が出たら終わる背景タスクか Monitor で受ける。状態は待たずに 1 回だけ問い合わせる（2026-09-25、feedback_no_sleep）
- **壊れた run は報告に持ち出さない**: 報告・図・表は有効な run だけで組み、無効の経緯は研究ログに短く残す（2026-09-29、feedback_skip_broken_runs）
- **片付け**: 消さずに各計算機の [`ecalj/trash/`](../trash) へ（リポジトリのものは追跡も外す）。ノウハウは past_log.md、直すべき点は直さずに TODOandQuestion.md。
  [`MD/`](.) は Claude が読むもの、人向けは ecaljdoc（2026-10-01、feedback_trash・feedback_doc_axis）
- **バンドの比較の図**（MLO-QSGW）: MLO バンドだけで比べる（sigm の内挿のバンドは出さない）。複数の run は重ねず横に並べる（[`Samples/kBT/LiTi2O4/plot_bands_side.py`](../Samples/kBT/LiTi2O4/plot_bands_side.py)）、
  反復の推移は縦に並べる（`plot_runs_iter_grid.py` など）。重ねるときも線は細く（0.6〜0.7 pt）。枠の題に計算の日時（`llmf.<N>run` の時刻）。
  図の数値を npz で隣に残す（2026-09-28、feedback_band_plots）
- **API キーなどの秘密を画面やファイルに出さない**: Materials Project のキーは `<ecalj>/MaterialProject.key`（`.gitignore`）か `MP_API_KEY`。確かめるときは長さや件数だけ（2026-10-01、reference_mp_api_key）
- 外付けドライブ `/media/takao/TAKAOMINI` の `ecalj_20261001`・`ecaljdoc_20261001` はバックアップ専用。そこで作業しない（2026-10-01、reference_backup_takaomini）

## 2. 計算機と環境

計算機の表の正本はここ（2026-10-02 20:51、user「handover に取りまとめて」。ecaljdoc の ForDevelopers §7 にあった表をまとめた）。

原本は t14 の `~/ecalj`（文書の ecaljdoc も 2026-10-02 から `~/ecalj/ecaljdoc/`。`~/ecaljdoc` は古いクローンで使わない。[ecaljdoc_publish.md](ecaljdoc_publish.md)）。ほかの機械へは `TOOLS/sync_ecalj_src.sh <host> [<dir>]` で送る（git archive + `SRC/.ecalj_rev` の刻印、送り先に git は無い。
`--samples` で Samples ごと、`--check-all` で各機の版。SRC に未コミットの変更があると止まる → `ALLOW_DIRTY=1`）。
数ファイルだけの差分なら `scp` で置いて `SRC/.ecalj_rev` を打ち直すと増分ビルドで済む（全部送ると mtime が変わって全部ビルドし直す）（reference_build）。

| 機械 | 中身 | 使い方の要点 |
| --- | --- | --- |
| t14（手元） | 16 コア、メモリ 30 GB（普段 20 GB 以上使用中）、gfortran-14（`mpif90` は Intel MPI のラッパー）、Python 3.13（mise） | [`SRC/build_gfortran`](../SRC/build_gfortran) → `~/bin`（シンボリックリンク）。確かめのビルドは別の木（`~/work/<x>tree`）で（試験中は `~/bin` の先を変えない）。`~/bin` には [`SRC/exec`](../SRC/exec) から退かせたファイルを指す行き先の無いリンクが残ることがある（`InstallAll.py` の次の実行で消える） |
| kt1 | RTX 5090 ×2（32 GB）、コア 0〜63（SMT 無し）、メモリ 251 GB、nvfortran 26.1（NVIDIA HPC SDK 2026）・CUDA 13、HPC-X の OpenMPI。gfortran 13/14 はあるが gfortran 用の MPI は無い | 本番 `~/ecalj` → `~/bin` は触らない。開発 `~/ecalj_dev` → `~/bin_dev`。何時間もかかる計算は実体をコピーした `~/bin_frozen_<rev>`。作業と試験用のツリーは `/mnt/data1`（`/mnt/data1/ecalj_test*`）。環境 `source /mnt/data1/LiTi2O4_kbt_runs/liti_src_full9/env.sh`。GPU は 2 枚とも使ってよい（2026-09-30 から。それまでは GPU 1 を user のジョブ用に空けていた）。GW のジョブと並べる試験は `-np2 1`。HPC-X の `mpirun` 無しで `lmfa` は動かない（`mpirun -np 1 lmfa`）。`lmf --listcmdopt` は正常でも exit 1 |
| kr7 | Ubuntu 24.04、Zen 5 の 16 スレッド、RTX 5090 32 GB、メモリ 30 GB | 試験用（2026-09-28 から）。sudo は user が打つ。kr5 と同じ構成: `source ~/nvenv.sh`（HPC SDK 26.1 を `~/opt/nvhpc`、MKL を `~/opt/intel`、uv の venv `~/venv`）。`python3 InstallAll.py --fc nvfortran --gpu --bindir ~/bin --notest --no-bashrc`、`sync_ecalj_src.sh --samples kr7` で Samples ごと送る、`-np 8`。試験に gnuplot が要る |
| kr5 | Ubuntu 24.04、Ryzen 7 9700X（**8 物理コア**／16 スレッド）、メモリ 30 GB、RTX 5090 32 GB、sudo 不可。kr7 と同じ構成（kr7 の MKL は kr5 から写した） | user が GPU を使う機械。2026-09-25 に LiTi₂O₄ を回せるようにして停止中（`~/liti/run_kr5.sh`、`-np 8 -np2 1`。`nproc` の 16 で `-np 14` にすると OpenMPI がスロット不足で落ちる）。ソースは `sync_ecalj_src.sh` で `SRC` と `InstallAll.py` だけ。最後に送った版は 2026-09-25（`sync_ecalj_src.sh --check-all`） |
| mic | RHEL 8.7、ifx・ifort（oneAPI 2023。ifx 2026 は `~/.local/bin/mpiifx` の包み）、ホーム `/home/cp/tkotani` | sudo なし、Python は uv の 3.12（`~/.local/bin` を先に）。試験用のツリーは `~/ecalj_test0928`、環境は `~/ecalj_test0928_env.sh`（ifx 2026 と venv `~/venv_ecalj`） |
| ucgw | SGE のクラスタ（ログインは 4 コア）、ifx 2024.2 | `module load intel/2024.2 intelMKL/2024.2 intelMPI/2021.13`。ジョブは `-pe x32 64`（x56・x64 は使えない）、`#$ -l mem_free=150G`（93 GB の機械で swap で死ぬ）、`I_MPI_HYDRA_BOOTSTRAP=sge`、ジョブの中で `module load`。`qsub` には `PATH` に `/usr/sge/bin/linux-x64`。`.git` が大きく `git push` が固まる（下の「肥大した `.git`」）。試験は `sync_ecalj_src.sh --samples ucgw ~/ecalj_test<日付>` で送った新しい木で、ジョブの中でビルドして `TOOLS/samples_tests.sh`（1 ノード `-pe smp 32`、優先度を下げるなら `#$ -p -100`。例 `~/ecalj_test1007/job_test1007d.sh`）。Python はジョブの中で venv `~/venv_ecalj`（pyenv 3.12.13 から。anaconda3 は bin の OpenMPI の `mpirun` が Intel MPI の代わりに走るので使わない。2026-10-08） |
| ISSP（kugui・ohtaka） | | past_log.md §1.3、ecaljdoc `UsageISSP.md`、[`Samples/BenchmarkTest/ISSP/`](../Samples/BenchmarkTest/ISSP/README.md) |

**kt1 で複数のジョブを並べるとき**（2026-09-30、reference_kt1_taskset）: ジョブごとに `taskset -c 16-27 ...` で別のコアの組を与え、**`export OMPI_MCA_hwloc_base_binding_policy=none` を必ず付ける**。
HPC-X の OpenMPI は起動した側の `taskset` を見ずにランク i を機械のコア i に固定するので、付けないと全部がコア 0〜11 に重なる（一度、コア 0 で NVIDIA ドライバのロック待ちが回り続けて
kt1 全体が 17 分止まった）。GPU のプログラムは [`SRC/exec/pylib/gpu_lock.py`](../SRC/exec/pylib/gpu_lock.py) の GPU ごとのロックで順番に GPU を使う。GW のジョブと GPU を分け合う試験は `-np2 1`（`-np2 2` は 2 枚同時に空くのを待って進まない）。
コアは 0〜63 だけ（`taskset -c 64-` は失敗する）。何かを足す前と後に、GPU のジョブが進んでいるかを確かめる。

**リモートの操作**（2026-09-30、feedback_remote_scripts）: ssh の 1 行に `pkill -f <pattern>` を書くと、その `bash -c` の行自身に当たって接続ごと落ちる。
kill や起動はスクリプトのファイルにして `scp` し、`ssh host 'bash /path/script.sh'` で実行する。パターンは `name.s[h]` のように最後の文字を括る。
手元でも同じ（2026-10-01: Claude のシェルで `pkill -f "lmf eute"` がそのコマンド行自身に当たり、後続が走らなかった）。止めるときは `pgrep` で PID を見てから `kill <PID>`。
待つときは親の PID で `while kill -0 <PID>; do sleep 15; done` のようにし、`pgrep -f` のパターンで待たない（待つ側の行にも一致する）。
背景ジョブは `nohup setsid ... < /dev/null &`。走っているスクリプトのファイルを書き換えない（bash は少しずつ読む）。

**肥大した `.git` を持つリモートへ main を入れる**（2026-06-05、project_ecalj_install_transfer）: `git push` が相手側の `index-pack --fix-thin` で固まる。
t14 で `printf '<main>\n^<相手のHEAD>\n' | git pack-objects --revs --stdout > main.pack`（`--thin` なし）→ `scp` → 相手で `git unpack-objects < main.pack`、
`git update-ref refs/heads/main <main>`。今は sync_ecalj_src.sh（git を使わない）を使うのが普通。

## 3. ビルドと試験

- **ビルドと試験の手順**は [ForDevelopers.md](ForDevelopers.md) §4・§5 が正本（2026-10-02 21:06 に重複を外した）。ここは落とし穴と、手順書に無い事実だけ
- **投入の前の 4 点**: (1) `sync_ecalj_src.sh --check-all` で版が HEAD と同じか、(2) `.so` に新しい印があるか、(3) 使う bindir の `gwsc` などがリンクか実体のコピーか（実体だと古いまま）、
  (4) 走り出したらログに期待する印が出ているか。2026-09-25 に三つとも見た目は正常なまま無効な計算を回した（feedback_verify_before_run）
- **凍結した bindir の作り方**（2026-10-01）: `cp -rL <bindir> ~/bin_frozen_<rev>`。開発ツリーに壊れたリンク（エディタの `.#name`）があると `cp -rL` が途中で失敗するので、
  失敗を見てから足りないものを入れる。Python の道具（`vasp2ctrl`・`getsyml`）は写すと隣のモジュール（`convctrl` など）を読めないので、元のツリーの bin から使う。
  `gwscconv` は自分と同じディレクトリの `gwsc` を呼ぶ
- **HPC-X の OpenMPI では、プログラムを `mpirun -np 1` なしで起動すると `MPI_Init` で止まる**。スクリプトから 1 プロセスで呼ぶときも `mpirun -np 1` を付ける（2026-09-30、`run_arg`・`hx0ahc.py`）

- **試験の最中に手元で `libecaljF.so` を作り直さない**（2026-10-02）: 走っている試験の次のプログラムが作り直し中のライブラリを読み、`heftet` が空の出力で止まった（`fe_kbt`）。試験が終わるのを待つか、別のビルド場所で
- **試験の最中に手元で別のものをビルドする**（2026-10-02 06:27）: `rsync -a --exclude 'build_*' ~/ecalj/SRC ~/work/<名前>tree/` で写し、その中で
  `mkdir build_gfortran && cd build_gfortran && FC=gfortran cmake .. && make -j4`。実行ファイルは同じ場所の `libecaljF.so` を読む（`$ORIGIN`）ので試験と混ざらない。
  cmake を設定し直すときは `FC` が要る（`FC=gfortran cmake .`。無いと「Fortran compiler must be set via FC」で止まり、古い .so が残る）
- **`sync_ecalj_src.sh`** は送り先の [`SRC/subroutines`](../SRC/subroutines)・`main`・`exec` にある HEAD に無いファイルを `trash/` へ移す（2026-10-02 から）。古い `.f90` が残ると CMake の GLOB が拾う。送り先の `bin` に残る古いプログラム（外したもの）は自動では消えない

- **t14 の trash は外付けの TAKAOMINI にある**（2026-10-02、ディスクが満杯になったため）: `/media/takao/TAKAOMINI/trash/`（`ecalj_trash/`、`work_20261001/`、`home_trash/`）。t14 のディスクは 468 GB で、試験の作業ディレクトリ（`Samples/*/*_work`、合わせて約 8 GB）や `~/work` の GW の途中のファイルでいっぱいになる。満杯になると Claude Code の Bash の出力も受け取れなくなる（`/tmp` も同じディスク）

- **対称操作は spglib から**（2026-10-02 08:34、[`MD/symmetry_spglib.md`](symmetry_spglib.md) §4.7g）: `symgrp = "find"` のとき lmf・lmchk などは同梱の spglib（[`SRC/external/spglib`](../SRC/external/spglib)、C、
  静的ライブラリ `symspg`）で操作を求め、作業ディレクトリに `symmetry.<sname>.json` を書く。次からは読み、構造が違えば作り直す。`symgrp` に生成元だけを
  書いた入力（対称性を下げる）は、その生成元から群を閉じる（`gensym`）。古い探し方（`symlat`・`symcry`、`ECALJ_SYMFIND=ecalj`）は 2026-10-02 21:33 に外した。ビルドには C コンパイラが要る（CMake の `project(... Fortran C)`）

- **MLO は Löwdin で直交化した関数**（2026-10-02、user の判断。`d10a63716`）: 実空間で規格化した MLO を各 k で O(k)^(−1/2) で直交化したもの（射影 Wannier）が
  ecalj の MLO。模型は H̃(R) だけ（O(R) = δ）、__cmlo（U・J・cRPA・マグノン）、MLO-QSGW の Σ、SOC もこの基底。対称性と軌道の名前を保ち、窓によらない。
  メッシュの外のバンドは生の (H, O) の模型より E_F 近くで良い（ecaljdoc mlo §6 表 M8）。`job_mlo --mlo_raw` は比べるためだけの以前の模型。
  HamRsMLO に直交化の印があり、それが無い古い HamRsMLO は読まずに止まる（job_mlo を回し直す）。
  最大局在化（MV、Sakuma の拘束つき）は試作の道具（`TOOLS/gadget/mlo_maxloc.py --sym`、2026-10-02 に SRC/exec から移した。bindir には入らない）として残すがメインにしない。拘束なしの MV は結合の方向にずれた混成軌道になる

## 4. コンパイラと実行時の落とし穴

- **S6 より前（2026-10-02 08:30 以前）に作った GW・MLO の中間ファイル**（`__HAMindex`、`__cmlo.*`、`__WV*` など）を今の版のプログラムで読むと、対称操作の並びが違って
  `rotwvigg: qtarget is not a star of q` や `rppovl: qi is not found` で止まる（2026-10-02 10:15、Fe のマグノンの作業場所で）。作り直すか、`ECALJ_SYMFIND=ecalj` で回す（それでも合わない所がある）
- **Fortran の `.and.` は短絡しない**: `if (present(x) .and. x)` は x が無いときも x を読んで落ちる。入れ子にする（past_log.md §3.2）
- **確保しただけの配列の中身**: gfortran は 0 になっていることが多く、ifort・ifx はごみが残る。未初期化に頼る誤りは Intel で初めて出る
  （2026-05-09 `iors.f90` の readrho: 基底を広げたとき新しい部分を 0 にしていなかった、`-1.32e124` の負の密度から `RSEQ : nit gt 80 and bad nodes`）（project_iors_nlm0_bug）
- **nvfortran**: 組み込みの `findloc` が一部の型で未実装（`FINDLOC: unimplemented for data type`）→ `m_nvfortran` の findloc を使う。同じファイルを 2 回 open すると
  `FIO-F-207` で止まる（gfortran・ifx は通る）。`newunit` の番号は負なので「`ifile < 0` なら未 open」という判定は使えない（2026-09-30）
- **gfortran 13・14 の誤コンパイル**（`m_HamPMT.f90`）: 旧形式の読み込みの到達しない `else` の本体を消すと、TOML の側の `lmindex` の配列が誤ってコンパイルされ、
  MLO で `zhev_tk2: nev /=nevx` になる。`aaa = trim(aaa) // ' '` の 1 行（「compiler bait」）を残してある。消さない（2026-05-07 発見、2026-10-01 確認、project_gfortran14_bait）
- **Intel の名前の衝突**: ローカル変数 `mpi_info` が Intel MPI の `MPI_INFO` と衝突（#6401）。module の `public` に type を二重に並べると ifort が #6251（project_mic_ifort_install）
- **OpenMPI の MPI-IO が固まる**: SIGKILL されたプロセスが `/dev/shm/sem.OMPIO_<file>` を残すと、次の同じ名前の `mpi_file_open` が永久に待つ。
  固まったら `/proc/<pid>/wchan` が `futex_do_wait`、`ls /dev/shm/sem.OMPIO_*` を消す（2026-05-08、project_mpiio_hang）
- **`nvfortran` と `mpi_bcast`**: `m_blas` の 2 関数だけ `include 'mpif.h'` のまま（generic の `mpi_bcast` が complex(4)+integer(8) を解決できない）（project_main_validation_20260603）
- **64 プロセスで hgw（W を作る部分、以前の hrcxq）がヒープ破壊**: `m_llw` の `WVIllwI` の `iw > niw` のガードが無かった（MPI の集団通信の後に置く。修正済み、project_hrcxq_fix）
- **部分 DOS（`job_pdos`）がプロセス数で変わった**（`co` の TEST 2 が `-np 6` で 0.14 ずれ）: `bandcal` は (k, スピン) の組をランクに配るので、1 つの k の 2 スピンが別ランクに
  分かれる。重みのファイル `__DWGT` に k ごとに両スピンを書き戻して相手のスピンを古い値で上書きしていた。(k, スピン) ごとのレコードに自分のスピンだけ書くよう直した
  （`3fa93489d`）。**(k, スピン) を配る処理で、k ごとのファイルに全スピンを書くと同じことが起きる**
- **GPU の Σc が NaN**: 非同期のカーネルと cuBLAS が別のストリーム（ecaljclaude.md「OpenACC 一般注意」）

## 5. 計算の中身の落とし穴（物理と数値）

- **Σ の準位の均しは、実軸の極の項と虚軸の積分で同じ核にする**: 虚軸積分の Gaussian の正則化（`sig = esmr/2`）は Gaussian の準位の均しそのもの。極の項だけ Fermi–Dirac にすると
  段差が打ち消さず、NiO 2³ の O 2s の対で 0.6 eV ずれた（2026-09-20、project_tsigmaw。ecaljdoc kBT.md）
- **金属の QSGW が反復で荒れる正体**: 小さい q の W_c(ω) の鋭いプラズモン極を Σ_c の実軸の極の項が 1 点で拾う（LiTi₂O₄ で 1.8 eV、幅 0.1 eV）。`wcsmear = true`（既定）
  （2026-09-18〜19、project_kbt_plasmon_pole。ecaljdoc kBT.md §3.5）
- **GPU の精度**: TF32（旧 `--mp` だけ）は条件の悪い誘電行列を壊す。既定は fp32（`--prec=fp32`、`--prec=tf32` は Σ_c の最後の積だけ TF32）。GEMMul8 は小さい積（辺 < 64 か mnk < 1e8）に使わない。
  OpenACC の非同期キュー 1 は既定ストリームの cuBLAS を待たないので、非同期にするなら m_blas の積も同じストリームに（project_hgw_speedup_20260927、ecaljclaude.md「OpenACC 一般注意」）
- **APW の G の選び方**（`m_igv2x.f90`）: `pwmode` の 10 の位が 0（`pwmode = 1`）なら |G| < √pwemax で q によらない基底、1（`pwmode = 11`）なら |q+G|。
  QSGW（Σ の実空間の内挿）は |q+G| で、バンドの図は |G| でないと基底の入れ替わりが跳びとして見える。両立しない（2026-09-24、project_pwmode_qgcut。[`MD/mlo_notes/sigma_mlo_design.md`](mlo_notes/sigma_mlo_design.md) §3.4）
- **MLO-QSGW のバンドの描き方**: [`Samples/kBT/LiTi2O4/draw_mloband.sh`](../Samples/kBT/LiTi2O4/draw_mloband.sh)（SCF の H を保存した χ̃ で書き出し、`SigRsMLO` を退けて凍結なしの `mlo --mlo`）。
  `job_mlo --mlofreeze`（凍結した `HamRsMLO` と今の `SigRsMLO` で基底が混ざり 0.3〜0.5 eV ずれる）と `lmf --band --mlo`（Γ 以外無効）は使わない（project_mlo_qsgw_sigma）
- **混合の履歴 `__mixm` の継承**: 前の run の履歴を継ぐと、Broyden の 1 歩目で非磁性の解へ落ちることがある（Fe 2.13 → 0.02 μB）。lmf は起動時に捨てる（`--keepmixm` で従来どおり）。
  混合のパラメタを疑う前に `__mixm`・`rst` の残骸を疑う（2026-09-18、project_mlo_night_20260918）
- **スラブの E_F が数 eV ずれたら ESM を疑う**: ESM は `ctrlg.<sname>.toml` の `[esm]`。無いと `esmsmves: ESM is off` と出るだけで止まらず、E_F が 4.4 eV ずれた例がある（FeMgO、2026-09-15、project_femgo_esm）
- **LDA+U**: `idu = 10 + mode` は `sigm` が無ければ mode の LDA+U、あれば U を切る（2026-09-30 に直した）。LDA+U では ehf が U の寄与を含まないので ehf と ehk は合わない（`m_lmfp`、仕様）
- **`lmf --jobgw=1` のメモリ**: IPW の重なり行列が 1 ランク 32·ngp² バイト（ngp ≈ V·Q³/6π²）。疎な構造（1 原子 500 Å³）で 1 ランク 35 GB。並列数に比例（2026-10-01、TODOandQuestion.md）

- **MLO の模型（基準 1・2・3、2026-10-01、user と決めた）**: 説明の正本は ecaljdoc mlo §9（§1 に入力の書き方）。
  - 基準 1（既定）: `mlo_lm` の EH ＋ 半内殻の局所軌道。局所軌道は**入力にキーが無く自動**（`m_HamPMT` の ShallowLO: 同じ種類の原子の LO をまとめた
    部分空間への重みが 1/2 を超える占有状態の上端で決める。E_F − 8 eV より上は EH と**入れ替え**（窓の中、Ni 3d・Zn 3d）、−17〜−8 eV は EH に**加える**
    （Ga 3d・Eu 5p・La 5p）、それより下は外す）。判定は `lmlo` の `local orbital atom ...` の行で見る。入れ替えにしないと、帯 1 本にシード 2 本が付き、
    部分バンドの模型（[`Samples/MLOsamples/NiO666lda`](../Samples/MLOsamples/NiO666lda)、Ni d ＋ O p、O 2p と Ni 3d の帯の試験）が `nskip` の検査で止まる（2026-10-01 19:2x）。
    自動にした理由は、判定に SCF のバンドの位置が要り、gwinit（SCF の前）では書けないため（user 19:1x に同意）。閾値によらず入れるときは `mlo_lm3`（その原子は書いたとおり）
  - 基準 2: gwinit が `mlo_lm2` に陽イオン（N O F P S Cl As Se Br Sb Te I 以外、遷移金属 Sc–Cu・Y–Ag・La–Au と Ac 以降を除く、Zn・Cd・Hg は含む）の s,p を
    **`!` 付き**で書く。**何もしなければ基準 1、`!` を外せば基準 2**（`job_mlo` だけ回し直す）。遷移金属・4f に EH2 を入れると、生の MLO の模型では壊れた（Cu・Ni は特定の k で崩れた）。
    Löwdin の模型（2026-10-02 から）では Cu・Ni も通る（重なりを内挿しないため。ecaljdoc mlo §9）。EuO は `Hreduction` の規格化の確かめで止まる（同じ原子の EH と EH2 の一次従属）
  - 基準 3: 空隙に空格子球（SiO₂ が要る。`[[site]]`・`[[spec]]`・`mlo_lm`、`lmfa` から）。自動の置き方は未完（[`SRC/exec/ctrlg_addes.py`](../SRC/exec/ctrlg_addes.py) は空隙を探す道具。ctrlg の生成への組み込みは TODO）
  - 選び方: 基準 1 で不満足なら基準 2、空隙があれば基準 3（2 と 3 を一緒も）。`mlo_delta`・`mlo_w` の細かい調整はしない
  - 評価は `mlo_bandcheck.py`: 窓 [VBM − 8, CBM + mlo_delta] の中の固有値のずれ（rms を 2 つの向き）と |Δgap| の最大。2026-10-01 17:54 より前の tsv・表は別の窓（CBM + 3 / + 1）
  - **スピン軌道の MLO（`job_mlo_soc`）はスピン軌道ありの DFT と比べる**。`job_mlo_soc` の段 1b が同じ条件で `lmf --band`（so=1）を回して `bnd*.spin1` を書く
    （2026-10-01）。SOC なしの DFT と比べると Δ_SO/3 ずれて見える（GaAsSoc で −0.1 eV に見えた）。`mlo_bandcheck.py` が WARNING を出す。
    MLO の窓の基準（E_F と CBM）は `--efermi=` のファイル（`efermi.lmf`、SOC は `efermi_soc`）から取る（2026-10-01 の夜まで E_F だけ `qplist.dat` の 1 行目だった）。
    `qplist.dat` の値は図の 0 点だけ
  - **模型の検査**: `mlo_bandcheck.py`（`job_mlo`・`job_mlo_soc` が最後に回す）が CHECK PASS/FAIL を出す。FAIL は 1 点のずれ > 0.1 eV、帯ごとのずれの隣の k との跳び > 0.1 eV、
    模型のバンドの帯のとげ（2 階差分、窓 [VBM − 3, CBM + 2] で 2 eV 超。線形独立性の崩れ。2026-10-02 から。それまでは重なりの最小固有値で見ていたが、
    模型が直交化した MLO で O = 1 なので替えた。ecaljdoc mlo §9 式 (12)）

- **MLO は実空間で規格化する**（2026-10-02、ecaljdoc mlo 式 (7a)）: 各 MLO を実空間の 2 乗積分 O_ii(R=0) の平方根で割る（k によらない定数）。`mlo` が `HamRsMLO` を作るときに決めて末尾に書き、以後の `Hreduction`（`__cmlo`、sugw、`m_sigmlo`）は同じ値で割る。k ごとに割る（旧 `--mlo_diagnorm`）と内挿バンドが動くので使わない。規格化の後に Löwdin で直交化するのが標準（上の「MLO は Löwdin で直交化した関数」）。古い `HamRsMLO` は読めないので `job_mlo` を回し直す
- **Wannier 関数・AHC・lmfham2 は外した**（2026-10-02）: git のタグ `last-wannier`（`dbcd6e51d`）に残る。cRPA は `job_mloW --crpa`、広がりは `mlo_spread.py`、マグノンは `job_mlo_magnon`。MLO と Wannier の違いは [`MD/wannier_vs_mlo.md`](wannier_vs_mlo.md)

## 6. 道具の癖

- **VSCode のチャットのリンク**: md は相対でも絶対でも開くが、画像（png）は開かない。図を見せるときは、図を埋め込んだ md へリンクし、パスは素のテキストで添える。
  md では画像をリンクで包む `[![名前](パス.png)](パス.png)`。md の表の中の数式で `|` を生で書かない（`\vert` か `\mid`）（reference_vscode_links）
- **ecaljdoc の図が ecalj のワークスペースのプレビューに出なかった**: Markdown Preview Enhanced がワークスペースを文字列の前方一致で判定したため（`/home/takao/ecalj` と
  `/home/takao/ecaljdoc`）。2026-10-02 に ecaljdoc を [`ecalj/ecaljdoc/`](../ecaljdoc/README.md) に入れたので、この問題は無くなったはず（未確認）
- **Materials Project の API**: 新しい ID（`mp-aaacpdie` の形、古い番号の 26 進）を返し、手元の `mp_api` の検証が止まる。読むだけなら requests で REST を直接
  （`x-api-key` と User-Agent。urllib の既定は 403）（2026-10-01、reference_mp_api_key、TODOandQuestion.md）
- **Python から libecaljF を呼ぶ**（構想、2026-03）: ctypes でシンボルを直接呼べる（nvfortran の名前は `<module>_<routine>_`、例 `m_bndfp_bndfp_`）。
  `bind(C)` は要らない（project_python_integration）

- **履歴の書き換え**（2026-10-02 12:59、user の指示）: 未公開の範囲の 657 コミットを書き換え、旧 `a0c7a7300`（新 `6e2903731`）に誤って入っていたビルドの生成物 3197 本を
  履歴から除いた。main のツリーは変わらない。dev・rel にあるコミットは変わらない。ecalj・ecaljdoc の文書とコミットメッセージの中のハッシュは新しいものに直した。
  古いハッシュ（他の機械の `SRC/.ecalj_rev`、作業場所のメモ、古いアーティファクト）は [`MD/commit_map_20261002.txt`](commit_map_20261002.txt) で引く。書き換え前の丸ごとの控えは
  `/media/takao/TAKAOMINI/ecalj_mirror_before_filter_20261002.git`
