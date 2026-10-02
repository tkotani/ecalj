# ForDevelopers — ecalj 開発の引き継ぎメモ

> 読者: ecalj の開発を引き継ぐ人と、次に作業する Claude のセッション（前の会話の記憶が無くても、ecalj と ecaljdoc だけで研究とジョブの投入ができることを目指す）。
> 約束ごと・ビルドとテストの手順・計算機ごとの注意・ジョブの投入・落とし穴・研究の現在地と経過を 1 枚にまとめ、詳しくはリンク先に任せる。
> 2026-09-28 時点。ここに書いたことが古くなったら、このページを直す。

**初めて読むときの順**: ecalj の `CLAUDE.md`（`ecaljdoc/MD/ecaljclaude.md` を読み込む。「記録の方針」が書き方の正本）→ このページの §1（push の約束）・§2 →
§9（研究の現在地）と §13（研究ログの要約）→ 作業に応じて §4・§5（ビルドとテスト）、§7・§12（計算機とジョブの投入）→ 各テーマの ecaljdoc のページ。
研究の細部は ecalj の研究ログ（§13 の索引から日付と時刻で引く）。

## 0. どこに何があるか

| 文書 | 中身 |
| --- | --- |
| ecalj `ecaljdoc/MD/README.md` | 2026-05〜09 の日付ごとの新機能（ログ）と TOML の流れのクイックスタート（英語）。2026-10-01 に最上位から移した（ecalj は Claude が読み、人間は Claude を通して情報を取る構造） |
| ecalj `Changes.txt` | 変更の記録（日本語）。上が新しい |
| ecalj `ecaljdoc/MD/handover.md` | 引き継ぎ: user との取り決め、計算機の癖、ビルド・試験・コンパイラ・数値の落とし穴（Claude の個人メモリから今も正しいものを写したもの、2026-10-01） |
| ecalj `ecaljdoc/MD/TODOandQuestion.md` | 直すべき点、メンテナへの質問、実行中のこと、やったこと（2026-10-01〜） |
| ecalj `ecaljdoc/MD/past_log.md` | リポジトリから外したもの（`trash/`）にあったノウハウと経緯。2026-06〜09 の要点（§6）、一発 GW と任意の k 線の落とし穴（§5）、hgw の統合の設計判断（§3.2） |
| ecalj `ecaljdoc/MD/ecaljclaude.md` | 設計方針（singleton module）、GPU 開発の教訓、hgw の構成（節の名前は hgw_combined）、GW1500 量産 |
| ecalj `ecalj_auto/README_slot_scheduler.md` | GW1500 量産（kt1、スロットスケジューラ、NaN 監視） |
| ecalj `Samples/TestInstall/README_testecalj.md`、[developer](./developer) | テストの仕組み（`testecalj`、`test.py`） |
| ecalj `ecaljdoc/MD/research_log.md` | 研究ログ（有限温度、MLO-QSGW、GPU 高速化、GW1500、片付け。2026-10-01 に `Samples/kBT/kBT_research.md` から移した）。上が新しく、時刻付き |
| ecaljdoc [ForDevelopers_research](./ForDevelopers_research) | 研究ログのテーマ別の要約（結論、訂正、未解決）、日付の索引、計算の置き場所（§13） |
| ecalj `Samples/kBT/gpu_fp32_report.md` | GW の GPU 高速化の報告（2026-09-27、§11 に要約） |
| ecalj `Samples/kBT/sigma_mlo_design.md` | MLO-QSGW の設計 |
| ecaljdoc [ecaljgpu](./ecaljgpu)、[gwsc](./gwsc)、[cmdopts](./cmdopts) | GPU 版の使い方、gwsc のオプション、コマンドラインの一覧 |
| ecaljdoc [mlo](./mlo)、[mlo_gwsc](./mlo_gwsc)、[kBT](./kBT)、[toml_migration](./toml_migration) | MLO、MLO-QSGW、有限温度、TOML 入力への移行 |
| ecalj `Samples/kBT/LiTi2O4/`（`input/qmlo`、`run_gwsc10.sh`、`cmp_gwsc10.py`、描画のスクリプト） | LiTi₂O₄ の MLO-QSGW の入力と、投入・比較・描画の道具（§12.2） |
| Claude のメモリ（t14 の `~/.claude/projects/-home-takao-ecalj/memory/`、索引は `MEMORY.md`） | 補助。約束・計算機・落とし穴で要るものはこのページと `ecaljdoc/MD/ecaljclaude.md` に移してある（§10） |

## 1. リポジトリと push

| リモート | 場所 | 使い方 |
| --- | --- | --- |
| `dev` | github.com/tkotani/ecaljdeveloper | 開発。まずここへ |
| `rel` | github.com/tkotani/ecalj（`main`） | 公開。dev で寝かせてから |
| ecaljdoc `origin` | github.com/ecalj/ecaljdoc | GitHub Pages で公開される文書 |
| ecaljdoc `dev` | github.com/tkotani/ecaljdoc | 文書の開発用 |

- **commit はその場で、push は指示があってから**（dev も）。まとめて `git push dev main` で出す
- **rel は慎重に**: dev で寝かせ、同じコミットで 3 環境（手元 gfortran、kt1 nvfortran GPU、ucgw か mic の ifx）の `testecalj --all` が通ってから
  `git push rel main`（fast-forward のみ、force は使わない）。破壊的変更は dev で警告だけの期間を置いてから止めるようにする。
  rel が壊れているときの修正だけは例外
- commit メッセージは**英語**（ecalj・ecaljdoc とも）。`Changes.txt`、ecaljdoc、研究ログの中身は日本語
- コミットには `Changes.txt`（日本語）と、使い方が変わるなら `README.md` の追記を含める
- `Samples/` の下の計算の出力（`*_work` など）はコミットしない
- local の main が dev よりどれだけ先かは `git rev-list --count dev/main..main`

## 2. 書き方の約束

コードのコメント、コミット、研究ログ、文書の方針は ecalj の `ecaljdoc/MD/ecaljclaude.md` の「記録の方針」が正本（ecalj の `CLAUDE.md` から
Claude Code に毎回読み込まれる。2026-09-27）。要点: コメントの注記には日付時刻、バグ修正は書く、将来混乱を招く数値は書かない、
コミットメッセージは英語、研究ログは最新が上で `### HH:MM`、数式に番号、利用者向けノートに苦労話を書かない。ここには ecaljdoc の書式だけ:

- Markdown の表のセルの中で `|` を生で書かない（`\vert`、`\mid` を使う）
- 図は md に埋め込んでその md を示す（VSCode のチャットでは png へのリンクは開かない）

## 3. コードの方針（詳しくは `ecaljdoc/MD/ecaljclaude.md`）

- **Fortran の module を singleton として使う**: 1 module = 1 責務 = 1 状態。状態は `protected, public` の module 変数、import は `use m_foo, only: ...`、
  subroutine の引数は指示的な値（`ef`, `qp`, `isp` など）だけで、配列データは `use` で取る。`type`（derived type）で状態を持たない
  （例外は `m_struc_def` の ragged array 用コンテナ）。module 間の `use` は DAG
- OpenACC: structured `!$acc data` の中で `BLOCK` を使えない（→ `enter/exit data`）。`goto` の経路でも `!$acc wait` が要る
- 入力は `ctrlg.<sname>.toml` だけ（`PB.<sname>.toml`、`esm_input.dat`、`GWinput`、`GWinput.toml` は読まない。`PB.<sname>.toml` と `esm_input.dat` は残っていれば止まるので `ctrlg_absorb.py <sname>`）。
  上書きは `--ctrlg:<path>=<value>`

## 4. ビルド

```bash
python3 InstallAll.py --fc gfortran --bindir ~/bin              # CPU。Python 3.11 以上
python3 InstallAll.py --fc nvfortran --gpu --bindir ~/bin       # GPU（kt1）。最後に linalgtune_gpu が表を書く
```

- `--fc` は**コンパイラの裸の名前**（`mpifort` のパスを渡すと `CMakeLists.txt` の判定を誤る）。ビルドは `SRC/build_$FC`、並列は 8 まで。
  `--notest`、`--clean`、`--notune`（計測を省く）
- 実体は `libecaljF.so` / `_mp` / `_gpu` / `_mp_gpu` の **4 本**で、実行ファイルは薄い入口。新しい機能が入ったかは
  `strings <build>/libecaljF*.so | grep -c <新しい文字列>` で 4 本とも確かめる
- `<bindir>` の実行ファイルとスクリプト（`gwsc` など）は build と `SRC/exec` への symlink。実体のコピーだと古いまま走る
- 新しい `.f90` を足したら cmake を configure し直す（`file(GLOB)`）
- nvfortran の `signal 11`（コンパイラの間欠的な落ち）は同じコマンドの再試行で通る
- 数ファイルだけ変えたときは、そのファイルを送って `make -j8`（全部送ると mtime が変わってフルビルド約 15 分）。
  リモートのソースの版は `SRC/.ecalj_rev`、`TOOLS/sync_ecalj_src.sh --check-all`
- GPU の行列演算の表（`<bindir>/ecalj_linalg_policy.toml`）は GPU・ドライバ・CUDA を替えたら `linalgtune_gpu` で作り直す（[ecaljgpu](./ecaljgpu)）
- **`m_HamPMT.f90` の `aaa = trim(aaa) // ' '` の行を消さない**（gfortran 13/14 の誤コンパイルを避ける「えさ」。消すと MLO が `zhev_tk2 nev` で落ちる）
- **確認用のビルドは別の worktree と別の bindir で**（`<bindir>` は build への symlink なので、作業ツリーで作り直すと使っている `~/bin` が変わる）:
  ```bash
  git worktree add --detach temp/ecalj_check <commit>      # 既にあれば: git -C temp/ecalj_check checkout --detach <commit>
                                                           # コミット前の変更を試す: git diff | git -C temp/ecalj_check apply（後で checkout -- SRC で戻す）
  cd temp/ecalj_check && python3 InstallAll.py --fc gfortran --bindir ~/ecalj/temp/bin_check --no-bashrc --notest
  cd Samples/TestInstall && ~/ecalj/temp/bin_check/testecalj -np 6 --all
  ```
  メモリが少ない機械では、worktree の中の `InstallAll.py` の並列数（`jobs = min(os.cpu_count(), 8)`）を下げてから。
  2026-09-27 に t14 でこの手順（並列 3）で確認: ビルド約 7 分。`5c35edb55` では `co` の部分 DOS だけ `-np 6` で失敗（§8）、
  修正（`3fa93489d`）を入れて `testecalj -np 6 --all` は全部合格

## 5. テスト

```bash
cd Samples/TestInstall; for d in *_work; do rm -rf "$d"; done     # 毎回、前の *_work を消す
testecalj -np 8 --all                                              # CPU 全部（26 ターゲット。有限温度は fe_kbt）
testecalj -np 8 -np2 2 --gwall --gpu [--mp] [--run-args=--prec=fp32]   # GW を GPU で（fp64 / tf32 / fp32）
cd Samples/MLOsamples; testecalj -np 8 Si666gwsc Cu NiO ...        # MLO（--all は TestInstall の中だけ。対象は README.md の表）
cd Samples/MLOQSGW;    testecalj -np 8 -np2 2 --gpu --mp GaAs NiO  # MLO-QSGW
```

- `testecalj` は走らせる各ターゲットの `<target>_work` を消してから作り直す（2025-10-09 から）。以前、前の失敗で残った `rst` による偽の失敗があったので、
  手で `*_work` を消す手順も残している（害は無い）。`-np` は手元（16 コア・30 GB）で 6〜12、kt1 で 8
- GPU＋MP のテストは照合するファイルが少ない（`--mp` では QPU の厳密照合を省く）。CPU のほうが照合は厳しい
- GPU では同じ入力でも GPU 0 と 1 で下の桁が違うことがある。比べるときは `CUDA_VISIBLE_DEVICES` を固定
- **2 原子のテストでは GPU の非同期のバグ（競合）は出ない**。GPU の実行順を触ったら、実物質（LiTi2O4 6³ など）で Σ を前の版と桁まで比べる
- QPU を比べるときは dSEnoZ の列（出発の固有値の列ではない）
- Samples の組ごとの試験は ecalj の `TOOLS/samples_tests.sh [--gpu] -np 8 [組 ...]`（inputs、TestInstall、EPS、PROCAR、MLOsamples、MLOQSGW、
  AFsymmetry、BenchmarkTest、重い GW、Magnon）。件数は testecalj の最後の要約で数え、途中で止まったターゲットは STOPPED にする
  （ログの PASSED 行は延べ数。2026-09-28）
- 試験に要る Python: numpy、pandas、matplotlib、pymatgen、seekpath（getsyml）、spglib、scipy、plotly。gnuplot は無くてもよい
  （図は飛ばして .glt だけ書く、2026-09-28 から）

## 6. GPU 開発の注意

- **行列積は `m_blas` の `gemm` を通す**。固定の行列には `key=` を付ける（変換した形を使い回す）。どの方法で計算するかは `m_linalg_policy` の表
  （精度 × 演算 × 大きさ → cublas / realsgemm / realhgemm / gemmul8 / 逆行列 lu64・mixed1・mixed2）。上書きは `--linalg=`
- **ストリームと非同期、`m_stopwatch` の同期、`acc routine` の中の自動配列**は ecalj の `ecaljdoc/MD/ecaljclaude.md`「OpenACC 一般注意」
  （非同期の区間の作り方は `m_sxcf_sc` の `sigma_stream_begin/end`）
- `!$acc update host` を Σc のループの中で使うと GPU 版の Σc が壊れた。hgw を同じ GPU で 2 本同時に走らせると W-build で segfault（未追跡）
- GEMMul8（Ozaki、INT8）は小さな積（辺 < 64 か m·n·k < 1e8）には使わない。計測は 64 の倍数を避けた寸法で（1024 などは GEMMul8 に有利に出る）
- 精度の性質: TF32・FP16 の誤差の大半はテンソルコアの**和の途中の切り捨て**による一様な縮み（入力の丸めは最近接）。
  W を作る側（χ0、誘電行列の逆）の誤差は $(1-v\chi_0)^{-1}$ で拡大されるので、精度を下げるのは Σc の最後の積だけにしてある
- nvfortran の癖: `findloc` は長さの違う文字列を空白で埋めずに比べる。`use mpi` の総称 `mpi_bcast` は complex(4)＋integer(8) を解決できない
  （`m_blas` の 2 関数は `include "mpif.h"`）。derived type の allocatable 成分の `present` は OpenACC で黙って失敗する
- プロファイル: rank 0 だけ `nsys profile -t cuda,openacc,osrt --duration=150` で包み（`OMPI_COMM_WORLD_RANK` で分岐する小さなラッパー）、
  `nsys stats --report cuda_api_sum`。同期（`cudaDeviceSynchronize`、`cuStreamSynchronize`）の回数と時間を最初に見る

## 7. 計算機

| 機械 | 中身 | 注意 |
| --- | --- | --- |
| t14（手元） | 16 コア、メモリ 30 GB（普段 20 GB 以上使用中）、gfortran（`mpif90` は Intel MPI のラッパー） | `~/bin` は作業ツリーの `SRC/build_gfortran` への symlink。確認用のビルドは §4 のとおり別の worktree で |
| kt1 | RTX 5090 ×2（32 GB）、64 コア、メモリ 251 GB、nvfortran 26.1（NVIDIA HPC SDK 2026）、CUDA 13。gfortran 13/14 はあるが gfortran 用の MPI は無い（HPCX は nvfortran 用） | GPU は 2 枚とも使ってよい（2026-09-30 から。それまでは GPU 1 を user のジョブ用に空けていた）。GW のジョブと並べる試験は `-np2 1`。本番 `~/ecalj`→`~/bin` は触らず、開発は `~/ecalj_dev`→`~/bin_dev`。作業は `/mnt/data1`。HPCX の `mpirun` 無しで `lmfa` は動かない（`mpirun -np 1 lmfa`）。`lmf --listcmdopt` は正常でも exit 1 |
| kr5 | Ubuntu 24.04、Ryzen 7 9700X（**8 物理コア**／16 スレッド）、メモリ 30 GB、RTX 5090 32 GB、`sudo` 不可。NVIDIA HPC SDK 26.1 を `~/opt/nvhpc` に、環境は `~/nvenv.sh`、MKL は kt1 から `~/opt/intel` に（研究ログ 2026-09-25 14:45） | user が GPU を使う機械。2026-09-25 に LiTi₂O₄ を回せるようにして停止中（`~/liti/run_kr5.sh`、`-np 8 -np2 1`。`nproc` の 16 で `-np 14` にすると OpenMPI がスロット不足で落ちる）。ソースは `sync_ecalj_src.sh` で `SRC` と `InstallAll.py` だけ |
| kr7 | Ubuntu 24.04、Zen 5 の 16 スレッド、RTX 5090 32 GB、`sudo` は user が打つ。kr5 と同じ構成（`~/opt/nvhpc`、`~/opt/intel`、`~/nvenv.sh`）、Python は uv の `~/venv`（研究ログ 2026-09-28 22:58） | 試験用（2026-09-28 から）。`sync_ecalj_src.sh --samples kr7` で Samples ごと送り、`-np 8` |
| mic | RHEL 8.7、ifort/ifx 2023（ifx 2026 は `~/.local/bin/mpiifx`） | Python は uv の 3.12（`~/.local/bin` を先に）。試験用のツリーは `~/ecalj_test0928`、環境は `~/ecalj_test0928_env.sh`（ifx 2026 と venv `~/venv_ecalj`） |
| ucgw | SGE クラスタ、ifx 2024.2（`module load intel/2024.2 intelMKL/2024.2 intelMPI/2021.13`） | 64 コアは `-pe x32 64`、大きな GW は `#$ -l mem_free=150G`、ジョブで `module load` と `I_MPI_HYDRA_BOOTSTRAP=sge`。`qsub` には `PATH` に `/usr/sge/bin/linux-x64`。`.git` が大きく `git push` が固まる → `git pack-objects --revs --stdout`（thin にしない）を送って `git unpack-objects` |

- `pkill -f <パターン>` は ssh 越しだと自分のコマンド行にも当たる。PID で止めるか `pgrep -f "pat[t]ern"` の形で
- 走っているシェル脚本を書き換えない（bash は脚本を逐次読む）。変えるなら止めて起動し直す
- 長いビルドは ssh を張ったまま走らせる（`nohup setsid` でも切断で死ぬことがある）

## 8. 既知の落とし穴

| 症状 | 原因・対処 |
| --- | --- |
| スラブで E_F が数 eV ずれる（止まらない） | ESM が効いていない。`ctrlg` の `[esm]`（`llmf` で `effective screening medium` の行を確認） |
| 磁性が消える・非磁性に落ちる | 前の run の混合履歴 `__mixm` を継いでいた。lmf は起動時に捨てる（続けるときだけ `--keepmixm`） |
| `mpi_file_open` で固まる | `/dev/shm/sem.OMPIO_*` の残骸（kill したジョブの）。消す |
| ifort/ifx だけ巨大な負の密度で落ちる | 未初期化の allocatable（gfortran は 0 で隠す）。例: `iors.f90` の基底拡張時の `nlm0`（修正済み） |
| `testecalj` が偽の失敗 | 前の `*_work` の `rst` |
| MLO-QSGW のバンドがおかしい | `job_band --mlo`、`job_mlo --mlofreeze` は MLO-QSGW のバンドにならない。[mlo_gwsc](./mlo_gwsc) §1.7 の方法で |
| QSGW とバンド描画で平面波の打ち切りが合わない | QSGW は `pwmode = 11`（$\vert q+G\vert$）、バンドは $\vert G\vert$。両立しない |
| 64 プロセスで hgw（W を作る部分、以前の hrcxq）がヒープ破壊 | `m_llw` の `WVIllwI` の `iw > niw` ガード（MPI の集団通信の後に置く。修正済み） |
| GPU の Σc が NaN | 非同期のカーネルと cuBLAS が別ストリーム（§6） |
| 部分 DOS（`job_pdos`）がプロセス数で変わる（`co` の TEST 2 が `-np 6` で 0.14 ずれ、`-np 8` は合格） | `bandcal` は (k, スピン) の組をランクに配るので、1 つの k の 2 スピンが別ランクに分かれる。重みのファイル `__DWGT` に k ごとに両スピンを書き戻して相手のスピンを古い値で上書きしていた。(k, スピン) ごとのレコードに自分のスピンだけ書くよう直した（`3fa93489d`）。**(k, スピン) を配る処理で、k ごとのファイルに全スピンを書くと同じことが起きる** |

## 9. 研究の現在地（2026-09-28）

テーマごとの結論・経過・未解決と研究ログの索引は §13。ここは一覧だけ。

- **MLO**（局在軌道の模型）: `mlo_method = 4` が既定、パラメタは `mlo_delta`・`mlo_w`。`Samples/MLOsamples` の系を `testecalj` で回帰。[mlo](./mlo)
- **MLO-QSGW**（Σ を MLO で持って k 空間で内挿）: LiTi₂O₄ 6³・9³ で 10 反復で収束し、メッシュ点の間のこぶが消えた。`gwsc N --mlo`
  （Σ^MLO は既定で `[gw] mixbeta` で混合する。`ECALJ_MLO_MIX=0` で切る）。6³ は tf32 と fp32 で LDA から 10 反復したバンドが rms 0.3 meV で一致（2026-09-27）。9³ の tf32 も LDA から 10 反復（3.5 時間）し、
  09-26 の fp32 の鎖と MLO バンドが rms 5.7 meV（窓の基準の違いを含む）、なめらかさは同じ（2026-09-28）。
  反強磁性は未確認。[mlo_gwsc](./mlo_gwsc)、ecalj の `Samples/kBT/sigma_mlo_design.md`
- **有限温度**: χ0 側（`t_tetrakbt`）と Σ 側の準位の幅（`t_sigmaw`、既定 1000 K）。金属の QSGW が反復で荒れる正体は、第一殻の q（offset-Γ ではない）の W のプラズモン極を
  Σc の実軸極項が踏むことで、`wcsmear`（既定 true）で均す。2026-09-27 に CoreEx の 2 つのバグと熱の核の積分（4 区間 × GL5）を直し、
  回帰テスト `fe_kbt` を足した。残る課題は [kBT](./kBT) §9
- **GPU の高速化**（§11）: 精度は `gwsc --prec=tf32|fp32|fp64` の 1 つで選び、方法は表で自動。LiTi₂O₄ 6³ の QSGW 1 反復は 342 → 212 秒（tf32）、
  `hgw` は 833 → 173 秒。hgw の間は GPU 2 枚とも 90〜95%・電力の上限（500 W）で回っている
- push していない仕事が多い（2026-09-28 に ecalj は dev より約 480、ecaljdoc も約 100 先）。push の前に §1 のゲートを通す

## 10. Claude と作業するとき

- 会話は日本語。**着手前に何をするかを一〜三行で言う**（長い作業ほど。途中で方針を変えられるように）
- **待つときに `sleep` を書かない**。完了の印（`steps.log` の `done`、ログの行、ファイル）が出たら終わるループを背景で走らせ、終わりの知らせで受ける（§12.2 に例）
- 投入の前の確認は §12.1、結果の残し方は §12.5。自律で作業したら、研究ログ（時刻付き）と報告書に残し、`Changes.txt` を更新してコミット（push は待つ）
- user の取り決め: kt1 の GPU 1 は user のジョブ用、kt1 の本番 `~/ecalj`・`~/bin` は触らない（開発は `~/ecalj_dev`・`~/bin_dev`）、push は指示があってから
- Claude のメモリ（t14 の `~/.claude/projects/-home-takao-ecalj/memory/`、索引 `MEMORY.md`）は補助。そこにしか無い約束や手順を見つけたら、
  `ecaljdoc/MD/ecaljclaude.md` かこのページに移す（記憶の無いセッションや別の機械のセッションが困らないように）

## 11. GW の GPU 高速化と QSGW 1 反復の短縮（2026-09-27）

報告は ecalj の `Samples/kBT/gpu_fp32_report.md`（§1〜10、FP16 経路の式は §2.1）、経過は `ecaljdoc/MD/research_log.md`、
変更の一覧は ecalj の `Changes.txt`（2026-09-27 (1)〜(3)）、利用者向けは [ecaljgpu](./ecaljgpu)。数値はすべて kt1（RTX 5090 ×2、電力上限 500 W）。

### 11.1 結果

| LiTi₂O₄ MLO-QSGW の `hgw` 1 回（`-np2 2`） | 09-26 までのコード | fp32 | tf32 |
| --- | --- | --- | --- |
| 6³ | 833 秒 | **367 秒**（2.3 倍） | **173 秒**（4.8 倍） |
| 9³ | 5004〜5128 秒 | 2990 秒（1.7 倍） | 1975 秒（2.6 倍） |

- 9³ の fp32・tf32 の値は 09-27 午前の版（午後の非同期化・FP16・k まとめなどの前）。**夜の版の tf32 では 9³ の `hgw` は 1164 秒（4.3 倍）**
  （2026-09-28、`gwsc 10` の反復 3〜9 の平均。中身は Σc 76%〔うち実軸の極項 6 割〕、W の構築 16%、Σx 6%。研究ログ 2026-09-28 03:05・03:11）

| LiTi₂O₄ 6³ の QSGW 1 反復（`gwsc 1`、tf32、MLO、GPU 2 枚・CPU 60 本） | 秒 |
| --- | --- |
| 15:50 のコード | 342 |
| hvccfp0 の GPU 化、hgw の午後の改良（18:34、`a74c9d118`） | 237 |
| hgw 以外の段（lmf --jobgw=1、mlo、hqpe_sc、hsfp0_sc、起動）（21:34、`e2cb402ef`） | 212 |

いまの内訳（212.8 秒、22:43、`b4bcb7adf`）: hgw 171（80%。うち Σc 115、W の構築 40）、lmf の SCF 13、lmf --jobgw=1 10、hvccfp0 ×2 8、
hsfp0_sc 4、mlo 3.5、hqpe_sc 1.5、heftet・hbasfp0 1。60 ランクのプログラムには 1 本 1.5〜2 秒の起動が含まれる。

### 11.2 精度

| `--prec` | 何を下げるか | 倍精度との差（E_F ±1 eV の Re Σc、最大） |
| --- | --- | --- |
| `fp64` | — | 基準 |
| `fp32` | 単精度の積を FP32 で | 0.02 meV |
| `tf32`（`--mp` だけも同じ） | そのうえ Σc の最後の積だけ入力の仮数を 10 ビットに（RTX 5090 の表では FP16） | 1.3 meV（TF32 の積なら 0.8） |
| （参考）すべての積を TF32 | 2026-09-27 までの `--mp` | 5 meV、NiO のギャップが 0.1 eV ずれる。tf32 より 1 割速いだけなのでやめた |

- tf32 の誤差の大半は Σc の一様な縮み（Re Σc の −5×10⁻⁵〜−1×10⁻⁴ 倍）。テンソルコアの**和の途中の切り捨て**が原因で、入力の丸めは最近接。
  打ち消し合う和の途中の値が大きいので拡大されて見える。安くは消せない
- W を作る側（χ0、誘電行列の逆）の誤差は $(1-v\chi_0)^{-1}$ で拡大されるので、精度を下げるのは Σc の最後の積だけにしてある
- 最後の反復だけ fp32: `gwsc 10 --prec=tf32 --prec-final=fp32:2`。tf32（FP16）は途中の反復に使い、最後を fp32 で締めるのが前提
  （FP16 は TF32 より精度が 1.6 倍悪い代わりに hgw が約 12% 速い。誤差は乱雑ではなく一様な縮みなので最後の fp32 の反復で消えるはず。報告書 §2.1）。
  **2026-09-28 の時点で `--prec-final` を使った計算はまだ無い**（6³ は tf32 だけの 10 反復でもバンドが fp32 と rms 0.3 meV で一致）

### 11.3 何が効いたか

`hgw` 1 回（6³、tf32。括弧は fp32）:

| 手 | 秒 |
| --- | --- |
| Hilbert 変換（dpsion）を GPU へ | 833 → 722 |
| 方法の振り分け表: 複素の単精度の積を実数 SGEMM 1 回に（realsgemm）、誘電行列の逆を FP32 の LU ＋ Newton（mixed1）、倍精度の積を GEMMul8 14 分解に | → 494 |
| 四面体の重みを空いている CPU コアで並走し、ファイルで渡す（`hgw --tetwt_write`、ビット一致） | → 438 |
| readeigen の MT 係数の小さな回転（1 回の hgw で 115 万回の倍精度 GEMM）を 1 カーネルに | → 424 |
| Σc のバッチの非同期化（積とカーネルを 1 本のストリームに、区間の中のタイマーの同期を外す、極の重みを GPU の裏で CPU） | → 412（fp32） |
| tf32: Σc の最後の積を 10 ビットの入力で（TF32 なら 274 秒） | → 242（FP16） |
| q の割り振り（LPT）の重みに W の構築分（2 ランクの終わりの差 14 → 5 秒） | → 237 |
| Σc の重み付けカーネルを平たい並列ループに（`ed946023e`） | → 233 |
| `build_zmel` の平面波の積を全状態で実数 SGEMM 1 回に（`d49bc8d78`） | → 218（388） |
| FP16 経路の B を重み付けのカーネルが直接 FP16 で書く（B の走査と変換をやめた、`cc9113900`） | → 197 |
| readeigen: 回転した固有関数を k 点ごとにデバイスに取り置く（`9bfea3432`、`ECALJ_WF_CACHE_GB`） | → 187（381） |
| χ0 を最大 8 k 点まとめてビンごとに積 1 回（`33c7b01b7`） | → 173（367） |

hgw 以外の段（QSGW 1 反復の中、秒）:

| 段 | 前 → 後 | 何をしたか |
| --- | --- | --- |
| hvccfp0 ×2（GPU） | 40.7 → 8.0 | Bessel 表・Wronskian・全原子対の Ewald 和（`strxq_all`）・CG 和を GPU で、原子群ごとの表は持ち回さない |
| lmf --jobgw=1（CPU） | 20.8 → 9.8 | H と xc 抜きの H を 1 回で（`hambl2`）、hsibl は全サイトの PW 係数を 1 回、pwmat は MTO の列だけ（FFT の相関）、Cholesky QR、mkpot 1 回 |
| mlo --mlofreeze | 10.1 → 3.5 | Σ^MLO だけを実空間へ、配列は主ランクへの reduce |
| hqpe_sc | 4.1 → 1.5 | Σ^MLO(R) の Bloch 和の位相を原子対ごとに（lmf の getsenex も同じ関数） |
| hsfp0_sc --job=3（GPU） | 6.6 → 4.4 | コア状態の M→E 変換を原子ブロックだけで、交換の積を split-K に（CPU 60 本では 153 → 57） |
| heftet・hbasfp0 ×2 | 1.8 → 1.1 | 1 ランクの CPU プログラムは GPU を隠して起動 |

- fp32 の誤差の主因は `build_zmel` の平面波の積（opA=N）の cuBLAS cgemm だった。realsgemm にして倍精度との差が 1/5 になった
- GEMMul8（Ozaki、INT8）が RTX 5090 で速いのは倍精度の積だけ（ネイティブ FP64 の 7〜13 倍）。単精度では同じ精度の FP32/TF32/FP16 に勝たない
  （INT8 と FP32 のピーク比が 8 倍しかない。比の大きい GPU では違いうる。報告書 §6.1）

### 11.4 仕組み

- **行列積は `m_blas` の `gemm` を通す**。固定の行列には `key=` を付ける（W(iω)・W(ω)、M2E の基底）。backend が変換した形（realsgemm の A'、GEMMul8 の分解）を
  キーごとに保持して使い回す（`ECALJ_LA_CACHE_GB`、既定 4 GB。q ごとに捨てる）。GEMMul8 の取り置きは m・k・分解数まで一致したときだけ使う
- **`m_linalg_policy` の表**: 精度（tf32/fp32/fp64）× 演算（cgemm/zgemm/dgemm/epsinv）× 大きさ（small/large）→ 方法
  （cublas / realsgemm / realhgemm / gemmul8[:分解数] / lu64・mixed1・mixed2）。優先順は 既定（全部 cuBLAS）< `<bindir>/ecalj_linalg_policy.toml`
  （GPU の名前が合うときだけ）< `--use_gemmul8` < `--linalg=`。使った表は標準出力に 1 回出る。`--sigma_tf32` では Σc の積（`policy=BACKEND_SIGMA`）だけが tf32 の行を使う
- **realsgemm**: 複素の C = A B を、B を実数の (2k×n) として読み、A だけを (2k×2m) の実数行列 A'（opA=N のときは A''）に組み替えて SGEMM 1 回。
  cuBLAS の cgemm（約 31 TFLOPS）より速く（50〜54）、和の取り方のおかげで精度もよい
- **realhgemm（tf32 の FP16 経路）**: realsgemm の FP16 版（入力 FP16、和と出力 FP32、`cublasGemmEx` の `CUBLAS_COMPUTE_32F`）。A' と B を 2 のべきで
  最大要素 [2^13, 2^14) にスケール（FP16 は 65504 で溢れる）。W 側の A' は key で取り置き、スケールは作るときに 1 回。Σc の積では重み付けのカーネルが
  FP16 の B を直接書き、スケールは上限 max\|w\|·max\|zmel\| から決める（一般の呼び出しは `cublasIsamax` で B を走査）。係数 1/(f_A f_B) は alpha に入れ、
  ホストは待たない。仮数は TF32 と同じ 10 ビットで、RTX 5090 では TF32 の 1.75〜1.9 倍速い。BF16 は同じ速さで誤差 8 倍。式と数値は報告書 §2.1
- **split-K**（`cmm_d`/`zmm_d` の `splitk=`）: 出力が小さく k が長い積（コアとの交換 158×158、k = 46 状態 × 788）を、k を等分した strided-batched の積と
  足し合わせのカーネルにする（1 回の積では数 SM しか使わない）。表・`policy`・`key` は見ない
- **GEMMul8** は小さな積（辺 < 64 か m·n·k < 1e8）には表に関係なく使わない（分解の固定費が回収できない）
- **計測ツール `linalgtune_gpu`**（`InstallAll.py --gpu` の最後、GPU が空いていれば。`--notune` で省く）: hgw の形 6 つ（寸法は 64 の倍数を避ける: 1037、8191 など）と
  逆行列 2 つで各方法の誤差と時間を測り、誤差の約束（単精度の複素積は同じ形での cuBLAS の誤差の 4 倍以内かつ tf32 5e-3・fp32 1e-5、倍精度 1e-12、
  逆行列 1e-6 / 1e-13）を満たし、hgw での時間の割合で重みを付けた幾何平均で cuBLAS より 1 割以上速く、重み 0.1 以上の形で 1 割以上遅くならないものを選ぶ。
  表は形ごとではなく組（small/large）ごと。tf32 の行は Σc の形（虚軸・実軸・zsec）だけで重み付け。表はファイルに一時名で書いてから rename
- **Σc の非同期化**（`m_sxcf_sc`）: バッチごとに `sigma_stream_begin`（`cudaDeviceSynchronize`、m_blas の積をキュー 1 のストリームへ）→ カーネルは全部 `async(1)`
  → 最後に `!$acc wait(1)` と `sigma_stream_end`。バッファの確保・解放はバッチの外。虚軸の重みは GPU のカーネル、実軸の極の重みは虚軸の積を投げた後に CPU で
- **χ0 の k まとめ**（`x0kf_v4h` の `accumulate_chi0_group`）: zmel を 1 回で作れる k 点を最大 8 個まとめ、ビンごとに gather と重みのカーネル 1 回、積 1 回。
  正の周波数のビン（npm=1）だけ扱う。npm=2（時間反転が破れた場合、今は `x0kf_zxq` が止める）は k ごとの `accumulate_chi0` が jpm を扱う
- **readeigen の取り置き**: `build_zmel` が同じ k 点を何度も求めるので、回転した固有関数を (k, スピン) ごとにデバイスに置く（`ECALJ_WF_CACHE_GB`、既定 2 GB、
  mpi_mode では使わない）。固有関数が変わらない間だけ有効（リセットしない）
- **コアとの交換**（`hsfp0_sc --job=3`）: コア状態の zmel はその原子の積基底ブロックの行にしか値が無いので、E 基底への変換を原子ごとにそのブロックだけで行う
  （k = ngb → nblocha）。交換の積は split-K
- **hvccfp0 の GPU 化**: `mkjp_4` が原子群（同じ動径メッシュと lx）ごとに Bessel 表 ajr と核を掛けた表 a1r を作り、`vcoul_termb` がその群のオンサイト項を足す
  （群の表だけを持つ）。`m_bessl` の作業配列は固定長（`acc routine` の中の自動配列はデバイスのヒープ確保になる。radkj・wronkj は lmax ≤ 21）。
  `strxq_all`（e<0）は全原子対の Ewald 和（q 空間は ZGEMM、実空間は (T, 対) ごと）と CG 和を GPU で作り、strx はデバイスにだけ置く
- **lmf --jobgw=1（ホスト）**: `hambl2` が H と xc 抜きの H を 1 回で作る（augmbl・smhsbl・hsibl がもう一方のポテンシャルの分も同じ所で足す）。
  xc 抜きのポテンシャル（spotx、oppix）は全体と同じ mkpot で（smves の後の smpot を写し、oppix は v1es・v2es で potpus〜gaugm を 1 回足す）。
  hsibl は全サイトの基底関数の PW 係数を 1 回作り（256 MB まで）、同じ打ち切りのサイトの並びごとに積 1 回。pwmat は APW の列を表引き、MTO の列は FFT の相関。
  Gram-Schmidt は Cholesky QR。GPU 版は hambl 2 回と従来の積のまま
- **mlo --mlofreeze**: HamRsMLO を書き直さない回なので、Σ^MLO(q) → QMLO_SigRs だけを行う（MLO の添字は HamRsMLO の末尾から読み、`__HamiltonianPMT` は読まない）。
  実空間の配列は主ランクへの MPI_Reduce。`sigmlo_sigq`（hqpe_sc と lmf の getsenex）は位相を原子対ごとに作る
- **四面体の重みの補助**: `gwsc --gpu` は `hvccfp0`・`hgw` の間、同じ `hgw` を GPU なしで `-np` の CPU コアで `--tetwt_write` として走らせ、
  `__TETWT.<iq>.<isp>` を書く。`hgw` は合うファイルがあれば読む（無い・合わないときは自分で計算）。MPMD にしなかったのは `MPI_COMM_WORLD` を直接使う箇所が多いため
- **GPU ロック**（`pylib/gpu_lock.py`）: `gwsc` は GPU を使う実行ファイルの前に `/tmp/ecalj_res/gpu<N>.lock` を必要な枚数だけ取り（全部取れるまで待つ）、
  取れた GPU を `CUDA_VISIBLE_DEVICES` で渡す。プロセスが終われば外れる。1 ランクの CPU プログラム（heftet、hbasfp0、hqpe_sc）は GPU を隠して起動する
  （`run_cmd`、MPI_Init が GPU を調べる分。60 ランクでは効かなかった）
- **q の割り振り**: `mpi_assign_qtask_lpt` の LPT。重み = 星の大きさ（Σc の量）＋ W の構築分（平均の Σc の 1/3）
- **画面出力とエラー**（2026-09-27 22:40）: GW のプログラムはランク 0 がそのままログ（`lgw` など、途中経過と `Memused`）に書き、ほかのランクは /dev/null。
  `--fullstdo` でランクごとの `stdout.<rank>.<prog>`。エラー終了（`rx` 系）は `rx_stop` にまとめ、stdout と stderr（ランク番号付き）に出して MPI_Abort。
  正常終了の `rx0s` も MPI_Finalize してから終わる

### 11.5 使い方

```bash
python3 InstallAll.py --fc nvfortran --gpu --gemmul8 --bindir ~/bin      # 最後に linalgtune_gpu が表を書く
gwsc 10 -np 60 -np2 2 --gpu --prec=tf32 <sname>                          # 最速
gwsc 10 -np 60 -np2 2 --gpu --prec=tf32 --prec-final=fp32:2 <sname>      # 最後の 2 反復だけ fp32
gwsc 10 -np 60 -np2 2 --gpu --prec=fp32 <sname>                          # 精度重視（= --mp --fp32）
```

| 指定 | 意味 |
| --- | --- |
| `--linalg=fp32.cgemm.large=realsgemm,...` | 表の行を上書き（`gwsc` に書けば全部に渡る） |
| `ECALJ_LINALG_POLICY=<file>` | 表のファイルを替える（試すとき） |
| `ECALJ_LA_CACHE_GB=<GB>` | キーで保持する行列の上限（既定 4、FP16 の分はその半分） |
| `ECALJ_WF_CACHE_GB=<GB>` | readeigen がデバイスに取り置く固有関数の上限（既定 2） |
| `ECALJ_GPU_WAIT=<秒>`、`ECALJ_GPU_LOCK=0` | ロックを待つのをやめる秒数、ロックしない |
| `--no-tetwt-helper` | 四面体の重みの補助を使わない |
| `--fullstdo` | ランクごとの画面出力を `stdout.<rank>.<prog>` に（調べるとき） |

GPU・ドライバ・CUDA を替えたら `linalgtune_gpu` を実行し直す。表が無いと全部 cuBLAS（以前と同じ動き）。

### 11.6 検証

- 2026-09-27 22:20（`d409f1958`、点検の反映の後）: TestInstall の `--gwall` を GPU の fp64・tf32・fp32 と CPU で回して全部合格。`Samples/MLOQSGW`（GaAs、NiO）も
  tf32・fp32 で合格（NiO の `log.nio` の差は tf32 0.0154、fp32 0.0008、許容 0.02）。手元の gfortran の `--all` は全部合格。
  LiTi₂O₄ 6³ の 1 反復は点検の前後と、画面出力の変更の前後で QPU がバイト単位で一致
- 2026-09-28 03:29（`24e3b6068`、有限温度の 2 つの修正・GL の区間分け・`fe_kbt` の入った版）: 手元の `--all` と、kt1 の TestInstall `--gwall`
  （`fe_kbt` を含む）を GPU の fp64・tf32・fp32 と CPU で、`Samples/MLOQSGW` を tf32・fp32 で、どれも合格
- 実行順を変える改良（非同期化など）は、LiTi₂O₄ 6³ の `hgw` で Σ が前の版と表示の桁（0.001 meV）まで一致し、2 回回しても一致することを確かめた
- 2026-09-28 22:38（`b4b64be66`、`t_tetrakbt` の必須化）: 手元の gfortran の `--all` の照合 64 件がすべて合格。件数は testecalj の最後の要約で数える（testecalj は
  ターゲットごとにそれまでの要約を出し直すので、ログ全体の PASSED 行は延べ数になる）。Samples の組ごとの試験は
  ecalj の `TOOLS/samples_tests.sh`（最後の要約だけを数える）
- 収束テスト（2026-09-27 23:35）: LiTi₂O₄ 6³ を LDA から tf32 で `gwsc 10`（`qmlo_k6_tf32n`、39 分）。fp32 の 10 反復（`qmlo_k6_gwsc10`）との差は
  10 反復目の後のバンドで rms 0.3 meV・最大 0.8 meV（MLO バンド、±3 eV）、QP シフトは毎反復 5 meV 以内で膨らまない。ehf は −22〜−91 meV
  （Σc の一様な縮みが占有状態の和に効く分）。表と図は ecalj の `ecaljdoc/MD/research_log.md` の 23:35

### 11.7 測り方（kt1）

- 1 反復: `/mnt/data1/LiTi2O4_kbt_runs/qmlo_k6_tf32h`（5 反復済み）をコピーして `run1.sh`（`gwsc 1 -np 60 -np2 2 --gpu --prec=tf32 --ntqxx --mlo`。Σ^MLO の混合は既定）
- `hgw` だけ: `/mnt/data1/LiTi2O4_kbt_runs/bench_hgw666`（6³）、`bench_hgw999`（9³）。10 反復後の状態で `hgw` を 1 回だけ回す（道具の写しと使い方は
  ecalj の `Samples/kBT/bench_hgw/`。`bench3.sh` の既定は GPU 0 の 1 枚、2026-09-27 の計測は `GPUS=0,1`）。
  `BENCH_BIN=~/bin_dev BENCH_PREC="--use_fp32 --sigma_tf32" ./bench3.sh <tag> 2`、その前に `hgw --tetwt_write` で重みを書いておく。
  倍精度の基準は `res_fp64_oz`、比較は `cmpef.py`（E_F からの窓ごとの最大と rms）、`cmpse.py`
- 行列積だけ: `/mnt/data1/tf32test/hgemmbench.cu`（FP32/TF32/FP16/BF16）、`TOOLS/ozbench`
- 各段の中身: GW のプログラムはランク 0 のタイマー（`Time:`）がログにそのまま出る。nsys は §6

### 11.8 残り

| 候補 | 見積もり | 注意 |
| --- | --- | --- |
| q 点ごと・原子ごとのファイル（`__Vcoud.<iq>`、`__PPBRD_V2_<ic>`、`__BASFP<ic>`）を 1 本に | 速さはほぼ同じ、ファイル数が減る | 2026-09-28 に調べて見送り（報告書 §10.1 の箇条: 書き手と読み手の場所、得が小さい理由）。ランクごとの stdout と gwsc の `PROCAR.UP.<rank>` は止めた。`__TETWT.<iq>.<isp>` は並走の仕組みなのでこのまま |
| hsfp0_sc（コア）の zmel のカーネル起動と同期 | GPU 1 枚で約 2 秒 | 原子ごとの小さな OpenACC 領域（1 回に約 130 起動、100 同期） |
| hqpe_sc の q ごとの matmul と混合履歴の読み書き | 約 0.7 秒 | 履歴ファイルの形を変えると走行中のチェーンの履歴が読めない |
| χ0 の積算（x0_gemm）を FP16 に | 6³ で数 % | W で拡大されるので Σ で確かめてから |
| zmel の構築（Σc 側）の非同期化 | 数 % | 小さな同期カーネルが多い |
| 表の行を用途（Σc、χ0、zmel）ごとに、単精度の候補に gemmul8:8 | RTX 5090 では 1〜2% | 別の GPU で効きうる。5090 だけでは確かめられない |
| 別の GPU（A100、H100 など）での表 | — | FP64 が速い GPU では倍精度の積は cuBLAS に、FP16 と TF32 の比、INT8 と FP32 の比も GPU で違う |

push はしていない（ecalj・ecaljdoc とも。push の前に §1 のゲート）。

## 12. ジョブの投入（2026-09-28）

記憶のないセッションでも同じ手順で投入できるように、実際に使っている形を書く。計算機の性格は §7、ビルドとテストは §4・§5。

### 12.1 投入の前に

1. **版**: リモートの `SRC/.ecalj_rev` が手元の HEAD と一致するか（`TOOLS/sync_ecalj_src.sh --check-all`。kt1 は開発の `~/ecalj_dev` を見る。
   送るのも `sync_ecalj_src.sh kt1` で `~/ecalj_dev` へ。2026-09-28 までは既定が本番の `~/ecalj` だった）。送るのは SRC と InstallAll.py だけで、
   `Samples`（新しい試験など）は `git archive HEAD Samples/TestInstall | ssh kt1 'tar xf - -C ~/ecalj_dev'` のように別に送る。
   変えたファイルだけを scp したときは `SRC/.ecalj_rev` を手で書き直し、全部が同じかはチェックサムで確かめる:
   `git ls-files SRC InstallAll.py | xargs md5sum | sort -k2` を手元と kt1 の `~/ecalj_dev` で取って diff
2. **実体**: `strings <build>/libecaljF*.so | grep -c <今回の変更で入った文字列>` が 4 本とも 1 以上（実行ファイルは薄い入口で、中身は `.so`）
3. **スクリプト**: `ls -la <bindir>/gwsc` が symlink か（実体のコピーは古いまま走る）
4. **走り始めの印**: ログに期待する行が出たか（MLO-QSGW なら `llmf` の `MLO Sigma interpolation ON`、`lmlo_sigr` の `promoted`）。
   出なければ、その計算は測りたいものを測っていない
5. **GPU**: kt1 は GPU 0 だけ（`CUDA_VISIBLE_DEVICES=0`、`-np2 1`）。2 枚は user が良いと言ったときだけ。
   先客は `nvidia-smi --query-compute-apps=pid,process_name --format=csv` で見る（瞬間の使用率ではどれが user のものか分からない）。
   このページや他のページの例の `-np2 2` は、2 枚を使ってよいときの形。`InstallAll.py --gpu` の最後の計測（`linalgtune_gpu`）も
   空いている GPU を使うので、kt1 では `CUDA_VISIBLE_DEVICES=0` を付けて走らせる。`/tmp/slot_scheduler.sock` があると `run_cmd` は
   量産用のスロットの GPU を使う（GW1500 の仕組み、ecaljdoc/MD/ecaljclaude.md。普段は無い）
6. **同じ機械の 2 本目**: GPU のロックは段ごとにしか取らないので、2 本の計算は GPU を交互に取り合い、CPU の段（`lmf -np 60` など）は
   ぶつかる。`hgw` は 1 ノードに 1 本（2 本だと /dev/shm の MPI の共有窓が足りなくなる、研究ログ 2026-09-19 21:15）。先客があれば終わるのを待つか、順番を user に聞く
7. **何時間もかかる計算**は、ビルドしても変わらないバイナリで（`<bindir>` はビルドディレクトリへの symlink なので、途中でビルドすると後の段が新しい版で走る）:
   ```bash
   cp -rL ~/bin_dev ~/bin_frozen_<rev>          # 実体のコピー。実行ファイルは RUNPATH の $ORIGIN で同じ所の .so を読む
   echo "<rev>: copy of ~/bin_dev, $(date)" > ~/bin_frozen_<rev>/FROZEN_REV
   ```

### 12.2 kt1: 長い GPU の計算（例: LiTi₂O₄ の MLO-QSGW を LDA から 10 反復）

入力は ecalj の `Samples/kBT/LiTi2O4/input/qmlo`（6³。9³ はメッシュの 3 行を変える。README に env.sh の作り方）、スクリプトは
`Samples/kBT/LiTi2O4/run_gwsc10.sh`。kt1 では作業場所 `/mnt/data1/LiTi2O4_kbt_runs`（`$S`）にスクリプトのコピーと入力
（`liti_src_full9`、`liti_src_full9_k9`）がある。**走っているスクリプトは書き換えない**（bash は逐次読む）。

```bash
ssh kt1
S=/mnt/data1/LiTi2O4_kbt_runs; cd $S
nohup env PREC=tf32 GPUS=0 bash $S/run_gwsc10.sh <tag> $S/liti_src_full9 ~/bin_frozen_<rev> > $S/<tag>.nohup.log 2>&1 < /dev/null &
                                      # GPUS=0,1 で 2 枚（許可があるとき）。kt1 の $S のコピーは 2026-09-28 03:20 にリポジトリの版にした
cat $S/<tag>/steps.log                # 版、.so の印、設定、精度、先客の GPU の数が最初に出る
```

- 進み: `$S/<tag>/gwsc10.log`（段ごとの経過時間と `QSGW iteration end iter N`）と `llmf.<N>run` の時刻。`steps.log` の反復ごとの行（`mloON=` と ehf）は
  `gwsc 10` が終わってから書かれる。反復ごとの中身は、`gwsc` が各反復の終わりに
  `QSGW.<N>run/` へ写す `lgw`・`lqpe`・`llmfgw01`・`lmlo_sigr`（2026-09-28 から）と `llmf.<N>run` を、[mlo_gwsc](./mlo_gwsc) §2.3 の表のとおりに見る
  （6³・9³ の正常値は `input/qmlo/README.md`）。`check_iter.sh` は `run_snap.sh` の鎖（`gwsc 1` を繰り返す形）用
- 終わり: `steps.log` の最後が `done`（失敗なら `ABORT`）。待つ間は自分のセッションを `sleep` で止めず、条件が成り立ったら終わるループ
  （中の `sleep` は見に行く間隔）を背景のタスクで走らせ、終わりの知らせで受ける:
  ```bash
  ssh kt1 'L=/mnt/data1/LiTi2O4_kbt_runs/<tag>/steps.log; while ! grep -q -E " done$|ABORT" $L; do
           pgrep -f "run_gwsc10.s[h] <tag>" > /dev/null || { echo "runner gone"; break; }; sleep 120; done; cat $L'
  ```
  （`[h]` は pgrep がこのコマンド自身の行に当たらないようにするため。`run_gwsc10.sh <tag>` のままだと ssh の先の bash 自身に当たって "runner gone" が出ない）
- 反復ごとのバンド（2026-09-28 から）: `run_gwsc10.sh` は `gwsc 10` が走っている間に、反復が終わるたびに MLO バンドを描く（`bnd_mlo_iter<N>.dat`、1 反復 約 1 分、
  `steps.log` に `bands of iter N` の行）。gwsc が `QSGW.<N>run/` に efermi.lmf・QMLO_SigRs・QMLO_z も残すので（`e1cf8e50f`、9³ で 1 反復 600 MB）、
  `draw_iter_bands.sh <tag> <from> <to> <bindir>` で後からも描ける。従来の sigm バンドは `SIGM_BAND=1` のときだけ
- 続き: `cont_gwsc.sh <tag> <from> <to> <bindir>`（`gwsc 1` を 1 反復ずつ、反復ごとに `snap/iter<N>` と MLO バンド）。LDA からも回せる（`$S/<tag>/` に ctrlg と env.sh を置いて `<from>` = 1）。
  **途中で止める**: `$S/<tag>/stop_after` に N を書くと、反復 N の後で止まる（走っている鎖の脚本を書き換えずに短くできる）
- 比べる: `python3 cmp_gwsc10.py <new> <old>`（反復ごとの ehf と QP、最後の MLO バンドと sigm バンド）。**バンドの図は MLO バンドで、複数の計算は重ねずに横に並べる**
  （user の指定、2026-09-28）: `mlo_rows.py`（[mlo_gwsc](./mlo_gwsc) の図と同じ描き方。`<root>/<列>/bnd_iter<N>.dat`、行 = 反復、`HILITE=1 MESHCOLS=… ROWS=… ROWH=…`）、
  `plot_bands_side.py`（占有の 2 本の拡大つき）。2 本だけ重ねるなら `plot_band_pair.py`（細い線）
- かかる時間（2026-09-27 夜のコード、GPU 2 枚）: 6³ tf32 39 分、9³ tf32 3.5 時間（1 反復 21 分、うち `hgw` 19.4 分）。fp32 の `hgw` は tf32 の約 2 倍（6³ で 367 対 173 秒、
  9³ で 43.5 分、1 反復 47 分。2026-09-28）。
  GPU 1 枚は今のコードでは測っていない（2026-09 の旧コードの 9³ fp32 で 1 反復 1.8 → 4.3 時間、2.4 倍。9³ を 1 ランクで回したときのメモリも確かめていない）
- 2 つの計算の精度だけを比べるときは、同じバイナリで回す（`qmlo_k9_tf32n` は `~/bin_frozen_9e881`）。9e881 の後の SRC の変更（有限温度の修正、GL の区間分け、
  gwsc の反復ごとのログ、PROCAR）は `t_tetrakbt = 0` の LiTi₂O₄ の結果を変えないので、新しい凍結コピーでもよい。反復ごとのログが `QSGW.<N>run/` に残るのは
  `daf41c19d` 以降の gwsc
- 1 反復だけ、`hgw` だけのベンチは §11.7

### 12.3 ucgw（SGE）

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N <name>
#$ -pe x32 64                 # 32 スロット × 2 ノード（x56・x64 は載るノードが無い）
#$ -l mem_free=150G           # 大きな GW。93 GB のノード（ucs11-17）ではスワップで止まる
#$ -j y
#$ -o job.log
source /usr/share/Modules/init/bash
module load intel/2024.2 intelMKL/2024.2 intelMPI/2021.13     # ジョブの中で読む（-V に頼らない）
export I_MPI_HYDRA_BOOTSTRAP=sge
export PATH=$HOME/bin:$PATH
gwsc 10 -np 64 <sname> > gwsc.log 2>&1
```

- 投入と確認: ssh 越しでは `export PATH=$PATH:/usr/sge/bin/linux-x64` の後に `qsub job.sh`、`qstat`、`qdel <id>`
- ノード: ucs11-17 は 32 スロット・62〜93 GB、ucs18-28 は 32 スロット・187 GB、ucs29-40 は 40 スロット・376 GB
- ビルド: 上の `module load` の後に `python3 InstallAll.py --clean --fc ifx --notest`
- ecalj の更新: ucgw の `.git` は大きく（約 10 GB、NFS）、`git push` は向こうの `index-pack --fix-thin` で固まる。thin にしない pack を送る（2026-06 に kt1・ucgw で使った形）:
  ```bash
  # 手元: ucgw の HEAD（<old>）から手元の main（<new>）までを 1 本の pack に
  printf '<new>\n^<old>\n' | git pack-objects --revs --stdout > /tmp/main.pack
  scp /tmp/main.pack ucgw:/tmp/
  # ucgw の ecalj で（固まった gc や receive-pack があれば先に止め、*.lock を消し、git config gc.auto 0）
  git unpack-objects < /tmp/main.pack && git update-ref refs/heads/main <new> && git checkout main
  ```
  チェックアウトを邪魔する未追跡のファイルがあれば退避してから。版は `SRC/.ecalj_rev` ではなく `git log -1`（ucgw は `sync_ecalj_src.sh` の対象に入っていない）
- ジョブの中の Python は 3.11 以上が要る（gwsc が tomllib を使う）。ucgw の計算ノードでどれを使うかはここには記録が無いので、最初に
  `qsub` の小さなジョブで `python3 --version` を確かめる
- `ecalj_auto/jobtemplate.ucgw`（GW1500 の量産に使った雛形）は `-pe smp`（1 ノード）と `-V`、`-q` でノードを並べる別の形。上の形は多ノードの 1 本用
- 本番の NiO の入力（メッシュ、時間とメモリの目安）はまだ無い。試験の `Samples/TestInstall/nio_gwsc444`（AF II、4³）が出発点になる

### 12.4 mic と手元

- mic: ifort/ifx 2023（ifx 2026 は `~/.local/bin/mpiifx`）。Python はシステムが 3.6 なので uv の 3.12 を `~/.local/bin` に置き、`PATH` で先にする。
  2026-09-17 に ifx 2026 で `testecalj --all` が通っている
- 手元（t14）: 確認用のビルドとテストは §4・§5。小さな有限温度の試験は `Samples/kBT/Fe`（bcc Fe 3000 K、GW 5³、`gwsc 1 -np 8 fe` が 40 秒）

### 12.5 結果を残す

- 研究ログ（ecalj の `ecaljdoc/MD/research_log.md`）に `### HH:MM` で、最新を上に。表と図には番号（*表 23:35-1* など）、図は png を
  `Samples/kBT/<系>/` に置いて md に埋め込み、描いたスクリプトも同じ所に
- 結論が固まったら ecaljdoc の該当ページ（[kBT](./kBT)、[mlo_gwsc](./mlo_gwsc)、[ecaljgpu](./ecaljgpu) など）と、このページの §9・§13 に移す
- コードを変えたら `Changes.txt`。commit はその場で、push は指示を待つ（§1）

## 13. 研究ログの要約と索引（2026-09-28）

研究ログ（ecalj の `ecaljdoc/MD/research_log.md`、2026-09-16〜、最新が上）のテーマ別の要約・日付の索引・計算の置き場所は
付録 [ForDevelopers_research](./ForDevelopers_research) にまとめた。ここはその見出しだけ。

| テーマ | いま成り立っていること（要点） | 残っていること |
| --- | --- | --- |
| 有限温度と金属 QSGW の荒れ | 荒れの正体は第一殻 q の $W_c$ のプラズモン極を Σc の実軸極項が踏むこと（09-19）。`wcsmear`（既定）と `t_sigmaw` で均す。2 準位の間に極が挟まる非対角の針は静的 QSGW の限界（09-20）。LiTi₂O₄ の実用設定 P（χ0 T=0 で Im χ0 を Gaussian で均す `t_tetrakbt = -992.4` + `t_sigmaw = 1000`）で 6³ が収束（09-22） | 300 K で Σ 側が 0.1 eV 動く件の再評価、絶縁体の `EFERMI_kbt`、Bose 項（[kBT](./kBT) §9） |
| 従来 QSGW の Σ(k) 内挿 | E_F 直上の凸凹はメッシュ点の間の内挿のはみ出し。Σ(R) がセル端まで減らないのが本体（09-23、09-25） | Σ(R) の窓掛けは未実施（MLO-QSGW へ移った） |
| MLO 模型 | `mlo_method = 4`、`mlo_nkabc` は Σ のメッシュに合わせる、(サイト, lm) あたり EH 1 枚、LiTi₂O₄ は全原子 s+p+d の 126 軌道（09-23〜25） | 76 と 126 軌道の比較、損失関数（MLO_v6_report §9） |
| MLO-QSGW | Σ を MLO の実空間 `QMLO_SigRs` で持つ。Σ^MLO も混合が要る（β=0.5）。6³・9³ で 10 反復、MLO バンドでメッシュ点の間のこぶが消える（09-26）。tf32 と fp32 は rms 0.3 meV（09-27） | 反強磁性、pwmode の違いの切り分け、`sigma_mlo_design.md` §13 |
| GPU の高速化 | `--prec` 1 つと方法の表。tf32 の誤差は Σc の一様な縮み（テンソルコアの和）。1 反復 342 → 212 秒、hgw の間は GPU 2 枚とも上限の 500 W（09-27〜28） | FP16 の縮みが TF32 の 2 倍な理由（報告 §2.1 の (b)(c)）、§11.8 の候補 |
| 運用 | 投入前の 4 点、長い計算は凍結したバイナリ、比べる前に作業ディレクトリを掃除、`/mnt/data1` に大物 | — |
