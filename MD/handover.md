# handover.md — 引き継ぎ（Claude の個人メモリから、今も有効なものを写したもの）

2026-10-01 に、t14 の Claude の個人メモリ（`~/.claude/projects/-home-takao-ecalj/memory/`、パッケージに入らない）から、別の機械や
記憶の無いセッションに引き継ぐ価値があり、今も正しいものだけを写した（user「メモリの内容はパッケージに入らないので MD/ に。ただし混乱を招くものは良くない」）。
その時点の状況（push の状況、夜間作業の途中経過）、古くなったもの、研究ログや TODO と重なるものは写していない。各項目の終わりの（ ）は確かめた日と元のメモの名前。

読む順: [../CLAUDE.md](../CLAUDE.md) → [ecaljclaude.md](ecaljclaude.md)（方針、記録の方針） → このファイル → [TODOandQuestion.md](TODOandQuestion.md)（未決と実行中） →
[research_log.md](research_log.md) の先頭（最新の経過）。片付けたもののノウハウは [past_log.md](past_log.md)。人向けの手引きは ecaljdoc（`manual/ForDevelopers.md` など）。

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
- **片付け**: 消さずに各計算機の `ecalj/trash/` へ（リポジトリのものは追跡も外す）。ノウハウは past_log.md、直すべき点は直さずに TODOandQuestion.md。
  `MD/` は Claude が読むもの、人向けは ecaljdoc（2026-10-01、feedback_trash・feedback_doc_axis）
- **バンドの比較の図**（MLO-QSGW）: MLO バンドだけで比べる（sigm の内挿のバンドは出さない）。複数の run は重ねず横に並べる（`Samples/kBT/LiTi2O4/plot_bands_side.py`）、
  反復の推移は縦に並べる（`plot_runs_iter_grid.py` など）。重ねるときも線は細く（0.6〜0.7 pt）。枠の題に計算の日時（`llmf.<N>run` の時刻）。
  図の数値を npz で隣に残す（2026-09-28、feedback_band_plots）
- **API キーなどの秘密を画面やファイルに出さない**: Materials Project のキーは `<ecalj>/MaterialProject.key`（`.gitignore`）か `MP_API_KEY`。確かめるときは長さや件数だけ（2026-10-01、reference_mp_api_key）
- 外付けドライブ `/media/takao/TAKAOMINI` の `ecalj_20261001`・`ecaljdoc_20261001` はバックアップ専用。そこで作業しない（2026-10-01、reference_backup_takaomini）

## 2. 計算機と環境

原本は t14 の `~/ecalj`・`~/ecaljdoc`。ほかの機械へは `TOOLS/sync_ecalj_src.sh <host> [<dir>]` で送る（git archive + `SRC/.ecalj_rev` の刻印、送り先に git は無い。
`--samples` で Samples ごと、`--check-all` で各機の版。SRC に未コミットの変更があると止まる → `ALLOW_DIRTY=1`）。
数ファイルだけの差分なら `scp` で置いて `SRC/.ecalj_rev` を打ち直すと増分ビルドで済む（全部送ると mtime が変わって全部ビルドし直す）（reference_build）。

| 機械 | 中身 | 使い方の要点 |
| --- | --- | --- |
| t14（手元） | gfortran-14、Python 3.13（mise） | `SRC/build_gfortran` → `~/bin`（シンボリックリンク）。`~/bin` には `SRC/exec` から退かせたファイルを指す行き先の無いリンクが 58 本ある（TODO） |
| kt1 | RTX 5090 ×2、コア 0〜63（SMT 無し）、メモリ 251 GB、nvfortran 26.1・CUDA 13、HPC-X の OpenMPI | 本番 `~/ecalj` → `~/bin` は触らない。開発 `~/ecalj_dev` → `~/bin_dev`。何時間もかかる計算は実体をコピーした `~/bin_frozen_<rev>`。試験用のツリーは `/mnt/data1/ecalj_test*`。環境 `source /mnt/data1/LiTi2O4_kbt_runs/liti_src_full9/env.sh`。GPU は 2 枚とも使ってよい（2026-09-30 から） |
| kr7 | RTX 5090 ×1、16 スレッド、メモリ 30 GB | sudo なし。`source ~/nvenv.sh`（HPC SDK 26.1 を `~/opt/nvhpc`、MKL、uv の venv）。`python3 InstallAll.py --fc nvfortran --gpu --bindir ~/bin --notest --no-bashrc`、`-np 8`。試験に gnuplot が要る |
| kr5 | kr7 と同じ sudo なしの構成（kr7 の MKL は kr5 から写した） | 最後に送った版は 2026-09-25（`sync_ecalj_src.sh --check-all`） |
| mic | RHEL 8.7、ifx・ifort（oneAPI）、ホーム `/home/cp/tkotani` | sudo なし、Python は uv（3.12）。ifx の MPI は `~/.local/bin/mpiifx` の包み。試験のときの環境は `~/ecalj_test0928_env.sh` |
| ucgw | SGE のクラスタ（ログインは 4 コア） | ifx: `module load intel/2024.2 intelMKL/2024.2 intelMPI/2021.13`。ジョブは `-pe x32 64`（x56・x64 は使えない）、`#$ -l mem_free=150G`（93 GB の機械で swap で死ぬ）、`I_MPI_HYDRA_BOOTSTRAP=sge`、ジョブの中で `module load` |
| ISSP（kugui・ohtaka） | | past_log.md §1.3、ecaljdoc `UsageISSP.md`、`Samples/BenchmarkTest/ISSP/` |

**kt1 で複数のジョブを並べるとき**（2026-09-30、reference_kt1_taskset）: ジョブごとに `taskset -c 16-27 ...` で別のコアの組を与え、**`export OMPI_MCA_hwloc_base_binding_policy=none` を必ず付ける**。
HPC-X の OpenMPI は起動した側の `taskset` を見ずにランク i を機械のコア i に固定するので、付けないと全部がコア 0〜11 に重なる（一度、コア 0 で NVIDIA ドライバのロック待ちが回り続けて
kt1 全体が 17 分止まった）。GPU のプログラムは `SRC/exec/pylib/gpu_lock.py` の GPU ごとのロックで順番に GPU を使う。GW のジョブと GPU を分け合う試験は `-np2 1`（`-np2 2` は 2 枚同時に空くのを待って進まない）。
コアは 0〜63 だけ（`taskset -c 64-` は失敗する）。何かを足す前と後に、GPU のジョブが進んでいるかを確かめる。

**リモートの操作**（2026-09-30、feedback_remote_scripts）: ssh の 1 行に `pkill -f <pattern>` を書くと、その `bash -c` の行自身に当たって接続ごと落ちる。
kill や起動はスクリプトのファイルにして `scp` し、`ssh host 'bash /path/script.sh'` で実行する。パターンは `name.s[h]` のように最後の文字を括る。
背景ジョブは `nohup setsid ... < /dev/null &`。走っているスクリプトのファイルを書き換えない（bash は少しずつ読む）。

**肥大した `.git` を持つリモートへ main を入れる**（2026-06-05、project_ecalj_install_transfer）: `git push` が相手側の `index-pack --fix-thin` で固まる。
t14 で `printf '<main>\n^<相手のHEAD>\n' | git pack-objects --revs --stdout > main.pack`（`--thin` なし）→ `scp` → 相手で `git unpack-objects < main.pack`、
`git update-ref refs/heads/main <main>`。今は sync_ecalj_src.sh（git を使わない）を使うのが普通。

## 3. ビルドと試験

- **ビルド**: `python3 InstallAll.py --fc <gfortran|ifx|ifort|nvfortran> [--gpu] --bindir <dir> [--notest] [--no-bashrc]`。`--fc` は裸のコンパイラ名
  （`FC=<path>/mpifort` は CMake の `MATCHES "ifort"` に当たって Intel のフラグになる）。ビルドの置き場は `SRC/build_<fc>`、並列度は 8 まで（nvfortran の ICE を避ける）。
  Python 3.11 以上が要る（`tomllib`・`contextlib.chdir`）。新しい `.f90` を足したら cmake の configure からやり直す（`file(GLOB)`）（reference_build、project_master_validation_20260509）
- **nvfortran の signal 11**: 毎回違うファイルで落ちる。同じコマンドをやり直せば通る。ループで回し、`libecaljF.so`・`_mp`・`_gpu`・`_mp_gpu` の 4 本すべてに新しい印が
  入ったかを `strings ... | grep -c <文字列>` で確かめる（実体は `.so`、実行ファイルは薄い包み。GEMMul8 の包みは `libgemmul8wrap.so`）（reference_build）
- **投入の前の 4 点**: (1) `sync_ecalj_src.sh --check-all` で版が HEAD と同じか、(2) `.so` に新しい印があるか、(3) 使う bindir の `gwsc` などがリンクか実体のコピーか（実体だと古いまま）、
  (4) 走り出したらログに期待する印が出ているか。2026-09-25 に三つとも見た目は正常なまま無効な計算を回した（feedback_verify_before_run）
- **凍結した bindir の作り方**（2026-10-01）: `cp -rL <bindir> ~/bin_frozen_<rev>`。開発ツリーに壊れたリンク（エディタの `.#name`）があると `cp -rL` が途中で失敗するので、
  失敗を見てから足りないものを入れる。Python の道具（`vasp2ctrl`・`getsyml`）は写すと隣のモジュール（`convctrl` など）を読めないので、元のツリーの bin から使う。
  `gwscconv` は自分と同じディレクトリの `gwsc` を呼ぶ
- **試験**: `testecalj` の前に `*_work` を消す（testecalj も各ターゲットの `_work` を作り直すが念のため）。件数は最後の要約で数える（ログの PASSED 行は延べ数）。
  組ごとの試験は `TOOLS/samples_tests.sh [--gpu] [-np N] [-np2 M] <組>`（組: inputs install gwall eps procar mlo mloqsgw afsym samples bench heavy magnon）（feedback_test_cleanup、project_samples_sweep）
- **HPC-X の OpenMPI では、プログラムを `mpirun -np 1` なしで起動すると `MPI_Init` で止まる**。スクリプトから 1 プロセスで呼ぶときも `mpirun -np 1` を付ける（2026-09-30、`run_arg`・`hx0ahc.py`）

## 4. コンパイラと実行時の落とし穴

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

## 5. 計算の中身の落とし穴（物理と数値）

- **Σ の準位の均しは、実軸の極の項と虚軸の積分で同じ核にする**: 虚軸積分の Gaussian の正則化（`sig = esmr/2`）は Gaussian の準位の均しそのもの。極の項だけ Fermi–Dirac にすると
  段差が打ち消さず、NiO 2³ の O 2s の対で 0.6 eV ずれた（2026-09-20、project_tsigmaw。ecaljdoc kBT.md）
- **金属の QSGW が反復で荒れる正体**: 小さい q の W_c(ω) の鋭いプラズモン極を Σ_c の実軸の極の項が 1 点で拾う（LiTi₂O₄ で 1.8 eV、幅 0.1 eV）。`wcsmear = true`（既定）
  （2026-09-18〜19、project_kbt_plasmon_pole。ecaljdoc kBT.md §3.5）
- **GPU の精度**: TF32（旧 `--mp` だけ）は条件の悪い誘電行列を壊す。既定は fp32（`--prec=fp32`、`--prec=tf32` は Σ_c の最後の積だけ TF32）。GEMMul8 は小さい積（辺 < 64 か mnk < 1e8）に使わない。
  OpenACC の非同期キュー 1 は既定ストリームの cuBLAS を待たないので、非同期にするなら m_blas の積も同じストリームに（project_hgw_speedup_20260927、ecaljclaude.md「OpenACC 一般注意」）
- **APW の G の選び方**（`m_igv2x.f90`）: `pwmode` の 10 の位が 0（`pwmode = 1`）なら |G| < √pwemax で q によらない基底、1（`pwmode = 11`）なら |q+G|。
  QSGW（Σ の実空間の内挿）は |q+G| で、バンドの図は |G| でないと基底の入れ替わりが跳びとして見える。両立しない（2026-09-24、project_pwmode_qgcut。`Samples/kBT/sigma_mlo_design.md` §3.4）
- **MLO-QSGW のバンドの描き方**: `Samples/kBT/LiTi2O4/draw_mloband.sh`（SCF の H を保存した χ̃ で書き出し、`SigRsMLO` を退けて凍結なしの `mlo --mlo`）。
  `job_mlo --mlofreeze`（凍結した `HamRsMLO` と今の `SigRsMLO` で基底が混ざり 0.3〜0.5 eV ずれる）と `lmf --band --mlo`（Γ 以外無効）は使わない（project_mlo_qsgw_sigma）
- **混合の履歴 `__mixm` の継承**: 前の run の履歴を継ぐと、Broyden の 1 歩目で非磁性の解へ落ちることがある（Fe 2.13 → 0.02 μB）。lmf は起動時に捨てる（`--keepmixm` で従来どおり）。
  混合のパラメタを疑う前に `__mixm`・`rst` の残骸を疑う（2026-09-18、project_mlo_night_20260918）
- **スラブの E_F が数 eV ずれたら ESM を疑う**: ESM は `ctrlg.<sname>.toml` の `[esm]`。無いと `esmsmves: ESM is off` と出るだけで止まらず、E_F が 4.4 eV ずれた例がある（FeMgO、2026-09-15、project_femgo_esm）
- **LDA+U**: `idu = 10 + mode` は `sigm` が無ければ mode の LDA+U、あれば U を切る（2026-09-30 に直した）。LDA+U では ehf が U の寄与を含まないので ehf と ehk は合わない（`m_lmfp`、仕様）
- **`lmf --jobgw=1` のメモリ**: IPW の重なり行列が 1 ランク 32·ngp² バイト（ngp ≈ V·Q³/6π²）。疎な構造（1 原子 500 Å³）で 1 ランク 35 GB。並列数に比例（2026-10-01、TODOandQuestion.md）

## 6. 道具の癖

- **VSCode のチャットのリンク**: md は相対でも絶対でも開くが、画像（png）は開かない。図を見せるときは、図を埋め込んだ md へリンクし、パスは素のテキストで添える。
  md では画像をリンクで包む `[![名前](パス.png)](パス.png)`。md の表の中の数式で `|` を生で書かない（`\vert` か `\mid`）（reference_vscode_links）
- **ecaljdoc の図が ecalj のワークスペースのプレビューに出ない**: Markdown Preview Enhanced がワークスペースを文字列の前方一致で判定するため（`/home/takao/ecalj` と `/home/takao/ecaljdoc`）。
  パスは正しいので書き換えない。`code -n /home/takao/ecaljdoc <md>` で別ウィンドウに開いてもらう（reference_vscode_links）
- **Materials Project の API**: 新しい ID（`mp-aaacpdie` の形、古い番号の 26 進）を返し、手元の `mp_api` の検証が止まる。読むだけなら requests で REST を直接
  （`x-api-key` と User-Agent。urllib の既定は 403）（2026-10-01、reference_mp_api_key、TODOandQuestion.md）
- **Python から libecaljF を呼ぶ**（構想、2026-03）: ctypes でシンボルを直接呼べる（nvfortran の名前は `<module>_<routine>_`、例 `m_bndfp_bndfp_`）。
  `bind(C)` は要らない（project_python_integration）
