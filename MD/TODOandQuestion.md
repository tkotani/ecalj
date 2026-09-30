# TODO と質問、やったこと

直すべき点は見つけてもその場では直さず、ここに書く（user 2026-10-01）。メンテナに決めてほしいことも、やったことも、ここに書く。
入口は [CLAUDE.md](../CLAUDE.md)（中身は [ecaljclaude.md](ecaljclaude.md)）。片付けたものの中のノウハウは [past_log.md](past_log.md)。
細かい経緯と数値は研究ログ [research_log.md](research_log.md)。

書き方: 各項目に見つけた日を付ける。済んだら「やったこと」へ移し、済んだ日とコミットを書く。

---

## 1. TODO（直すべき点）

### コード

- **`sugw`（`lmf --jobgw=1`）のメモリ**（2026-10-01）: `GEIGpart` が IPW の重なり行列 `ppovl(ngp,ngp)` と LU 用の写し `ppovlLU` を
  各ランクで持つ。1 ランクあたり 32·ngp² バイト、全体は並列数に比例。ngp ≈ V·Q³/6π²（Q = `QpGcut_psi`）なので胞の体積の 2 乗で増える
  （Rb8 の 4410 Å³ で 1 ランク 35 GB）。案: (a) 投入前に ngp から並列数を決める、(b) O は正定値なので Cholesky の因子だけ持つ（半分）、
  (c) O·x を FFT で作り反復法で解く（ngp に比例）
- **`auto_mpquery.py` が今の Materials Project で動かない**（2026-10-01）: MP の API が新しい ID（`mp-aaacpdie` の形、古い番号の 26 進）を返し、
  手元の `mp_api` の検証（pydantic）が止まる。`mp_api` を上げるか、REST を直接読む（requests、`x-api-key`、User-Agent が要る。urllib の既定は 403）
- **GW1500 の選定に構造の確かめを入れる**（2026-10-01）: 別の化合物から副格子を抜き出した MP の項目（ICSD の備考が「〜 part」、
  1 原子の体積が同じ組成の 2 倍以上など）が PBE のギャップ > 0 で選ばれていた（`ecalj_auto/GW1500_status.md` 表 6 の `INVALID_STRUCTURE`）。
  `auto_mpquery.py` で弾くか、印を付ける
- **`m_bndfp` が `m_clsmode_finalize` に渡す `ndimh`**（2026-09-30）: module の状態のまま。`--cls` を `pwmode = 11` で使うときに確かめる
- **GPU の build で module が循環している**（2026-10-01、`MD/module_map.md` を作って見つけた）: `m_mpi`（`MPI__AutoSetup` の中で `use m_gpu`）→
  `m_gpu`（`__GPU` のとき `use m_blas`）→ `m_blas`（`use m_gemmul8`）→ `m_gemmul8`（`gemmul8_init` の中で `use m_mpi, only: ipr`）→ `m_mpi`。
  まっさらな状態からは順序が決まらない。いまは `SRC/CMakeLists.txt` が GPU の変種を基本の変種の後に作り、`-I` で基本の `mod/` を見せる回避策で通っている
  （そのコメントは「CMake の走査が m_gemmul8 → m_mpi を見落とす」と書くが、見落としではなく循環）。案: `m_gemmul8` の `ipr` を下の層から取る、
  または `MPI__AutoSetup` の GPU メモリの問い合わせを引数で受ける。直したら GPU（kt1・kr7）で clean build と試験
- **AFTEST（`mmtarget.aftest`）の残り**（2026-10-01、研究ログ 2026-10-01 朝 06:46 と表 06:46-1）: 場をスピン 1 にしか入れていなかった誤りはブランチ `aftest-fix`
  で直した（下の質問）。残り: (a) モーメントを下げる向きで更新が行き過ぎる（利得 2 固定と 0.1(m − m_t)² の項。NiO 1.28 → 1.0 は 70 反復で、途中で符号が反転）。
  応答を割線で見積もる更新にするか。(b) 対がサイト 1・2、ブロック 1・2 の決め打ち（NiSe だけ 6 ブロック）。`AF=` の印や種の名前から対を決める。
  (c) `m_ldau_init` が lmf の起動のたびに場を一歩更新して `mmagfield.aftest` を書き直す（`job_band` でも）。(d) ecaljdoc の `UsageDetailed.md` の「直す必要がある」を、
  使い方（LDA+U のブロックが要る、U = 0 でよい、`SYMGRPAF`、目標のモーメントの定義）に書き直す
- **`Samples/MATERIALS/Database` の 62 物質を LDA で一度回して確かめる**（2026-10-01、急がない）: 入力は `lmchk` で読めることだけ確かめた

- **`hgw` の残り**（2026-04 の統合の残タスク、2026-10-01 に past_log.md §3.2 から）: ノード内の W の共有を `MPI_Win_allocate_shared` で（メモリの重複を減らす）。
  `hsfp0_sc` の Sx（`--job=1`）・core の交換（`--job=3`）も `hgw` に入れる（時間は小さいので優先度は低い）

- **`m_tetrakbt.f90` の使われていないルーチン**（2026-10-01）: 中点分解の `tetrakbt`・`eaf_triangle`・`eafww`・`factri`・`factri0`・`integral_1t`・
  `integral_0t`・`funcgx`・`numintall`・`funcgx2`・`ndiv_tt`・`int_simpson`（約 300 行）はどこからも呼ばれない（`tetwt5.f90` は `use` の only に
  `tetrakbt` を挙げるだけ）。今も使うのは `tetrakbt_init`・`kbt`・`integtetn`。2026-06 に「参照のため残す」としたもの。消すなら kBT の試験で確かめる（past_log.md §13）

- **MLO の自動の模型の既定**（2026-10-01、`Samples/MATERIALS/Database` の試験。研究ログ 2026-10-01）: 既定（`mlo_lm` が Ne まで s,p・Na から s,p,d、Δ = w = 2 eV）で
  多くの物質はバンドの rms が 0.01 eV 前後だが、次の 2 つの型で外れる。どちらも追加の動径関数で直る（*表 1*、`~/work/mlocheck_variants`）。
  (a) 陽イオンの半内殻 d と陰イオンの 2p の混成（GaN、InN、両相）: `m_HamPMT` の「浅い局所軌道」の判定（局所軌道が主の状態の上端が E_F − 10 eV より上）で
  Ga 3d（−11.8 eV）・In 4d（−12.5 eV）が「深い」とされ、模型から落ちて VBM が 0.33〜0.6 eV 下がる。`mlo_lm3` に d を足すと rms 0.002〜0.008 eV。
  合っている GaAs・InP・InAs は −13.9〜−14.8 eV。閾値を −13 eV にするか、陰イオンの p との近さで決めるか。
  (b) 空隙の大きい構造（MgS・MgSe・MgTe・CdTe・ZnTe・AlN・SiO₂ クリストバライト）: 伝導帯の底が高すぎる（+0.04〜+0.29 eV、SiO₂ は +4.1 eV で伝導帯が丸ごと無い）。
  `mlo_lm2` に全原子の s,p（EH2）を足すと 0.00〜0.05 eV。gwinit が書くコメント「同じ (atom, lm) に 2 本目（mlo_lm2）を足さない。重なりが特異に近くなる」は、
  これらの系では当たらなかった（`m_hreduction` の規格化の確かめも警告なし）。ただし常に足すのは不可: Cu は s,p の EH2 を足すと**止まらずに**模型が壊れる
  （rms 0.012 → 0.66 eV、最大 8 eV）。EuO は `Hreduction: PMT completeness loss too large` で止まる。SrTiO₃ は Sr 4p の半内殻を足すと 0.045 → 0.003 eV、
  La₂CuO₄ は La 5p で 0.094 → 0.001 eV（(a) の型。La₂CuO₄ では判定の表示が `********` で、局所軌道が主の状態が見つかっていない）。
  足す条件（陰イオンの 2p があるときの半内殻、空隙の大きい構造の EH2）を決めるか、模型を作った後にバンドを確かめて足す仕組み（`mlo_bandcheck.py`）にするか。
  模型が黙って壊れる場合（Cu）の検知も要る
- **`mlo` が so = 1 の入力で NaN**（2026-10-01、GaAs_so）: `job_mlo` を so = 1 の ctrlg で回すと `zhev_tk4: the Hamiltonian contains NaN`。
  SOC は `job_mlo_soc`（摂動）で扱う設計なので、`mlo`（か `job_mlo`）が so = 1 を見て止まり、`job_mlo_soc` を案内するのがよい
- **`job_mlo_soc` が空の `band_MLO_spin2.dat` を書く**（2026-10-01）: 2N のスピノルの帯は `band_MLO_spin1.dat` に入る。空の spin2 があると
  `mlo_bandplot.py` が空の枠を描く（`bandplot_MLO.isp2.glt` も残る）
- **`Samples/AFsymmetry/NiO` は `pwmode = 1` で `symgrpaf`**（2026-10-01）: k 点を 4³ にすると `rotwave: q+G rotation error (We have to set PWmode=11 for symgrpAF)`
  で止まる。試験の 3³ では通っているだけ。`pwmode = 11` にして参照を作り直すか

- **`InstallAll.py` が bindir の行き先の無いリンクを消さない**（2026-10-01）: t14 の `~/bin` に 58 本（`uutest`、`genMLWFmod`、`FLEX_interaction.py` など、
  `SRC/exec` から退かせたファイルと、エディタの一時ファイル `job_mlo~`・`#ctrlgenM1.py#`、`TAGS`）。インストールのときに、自分が作ったリンクのうち
  行き先の無いものを消すか（install manifest と照らす）

### 試験と入力

- `MLOsamples/RuO2` の保存してある `rst`・`dmats` は、`pwmode = 11` の LDA+U の誤り（2026-03-30〜09-30、`2498e5283` で修正）の時期に作ったもの。
  作り直すか（GdCo5・SmP は 2026-09-30 に作り直した）
- ecaljdoc に `QforGW`（メッシュの外の q の一発 GW）の落とし穴が無い: `EMAXforGW` が必須（1d20 の詰め物）、窓を変えたら交換からやり直す、
  `epsWVR` の行で小さい q の誘電関数が見られる（past_log.md §5）。spectrum.md か gwsc.md に書くか
- ecaljdoc に、TOML の流れの最短の手順（新規: `ctrls` → `ctrlgenToml.py` → `lmfa`・`lmf`・`gwsc`、旧い作業ディレクトリの移行、`--ctrlg:` の上書き）を
  ある程度まとめる（2026-10-01、user「ecaljdoc にある程度は書く」）。元の英語のクイックスタートは `MD/README.md`
- 試験の入力の温度（`t_tetrakbt = 262`、`t_sigmaw = 0`）をテンプレート（300/300）に揃えるか。揃えると gas_gwsc・fe_gwsc などの参照が動く

## 2. 質問（メンテナに決めてほしいこと）

- **ecaljdoc の古い文書**（`BackUp/`、`ecaljdetails/` の LaTeX・PS、2019 年以前）は、trash に移した `TOOLS/checkmodule`・`TOOLS/ModuleCodingSample` などを参照している。ecaljdoc の側も同じ決まり（trash へ、要点は過去ログへ）で片付けるか（2026-10-01）
- **GW1500 の `INVALID_STRUCTURE`（12）・`SUSPECT_STRUCTURE`（4）を集合から外すか**。いまは注記だけ
- **ビルドの生成物が入ったコミット `a0c7a7300`（78 MB）を、push の前に履歴から消すか**。消すと以後 674 コミットのハッシュが変わり、
  研究ログなどに書いたハッシュが合わなくなる。消さなければ、公開のリポジトリの履歴に 78 MB が残る
- **Materials Project の API キーの作り直し**（user の作業）: 2025-04〜2026-09 の公開リポジトリの履歴に入っていた
- **trash を空にするか**: t14 の `~/ecalj/trash`（旧ツリー 58 GB ほか）、kt1 の `~/ecalj/trash`（255 GB）。kt1 は空にするまで `/` の空きが 49 GB のまま
- 試験に使った古いツリー（mic の `~/ecalj_test0928`、kt1 の `/mnt/data1/ecalj_test0930b`・`0930c`）も trash に入れるか
- ブランチ `fix-idu10`（main にマージ済み）を消すか
- **ブランチ `aftest-fix`（`6d6f3451b`）をマージするか**: AFTEST の場をスピン 2 にも逆符号で入れる。afsym では SCF の中身は変わらず ehk だけが制約の下の全エネルギーになり、afsym なしでも正しくなる（NiO で afsym の結果と 10⁻⁶ eV で一致）。AFTEST を使う試験は今は無い
- push: dev・rel とも、t14 の main より 674 コミット遅れ（2026-10-01）

## 3. 実行中

- kt1 `/mnt/data1/gw1500_rerun/run3`: GW1500 の NOTCONV の残り（fp32、5 本）。2026-10-01 04:40 に約 17 時間の見積もり

---

## 4. やったこと（新しい順）

### 2026-10-01

- AFTEST を調べた（研究ログ 2026-10-01 朝 06:46）: afsym ではモーメントを目標に保てるが、表示の ehk に −uhx·m_d、ehf に −2·uhx·m_d が残る。afsym なしでは誤り（サイトの電荷が分かれる）。修正はブランチ `aftest-fix`
- GW1500: kr7 で fp32 と TF32 を同じバイナリ・同じ入力で比べた（8 物質、02:03〜06:24）。最終のギャップの差は 1 meV 未満。5 月との 0.29〜1.48 eV の差は精度ではなく、5 月の振動と設定の違い（`GW1500_status.md` §5.1 の表 7）。GOOD 1120 を精度の理由で見直す必要は無い
- `SRC/subroutines/m_tetrakbt_BUGREPORT.md`（2026-06、直し済みの不具合の報告）の要点を past_log.md §13 に移して trash へ。空の道標 `SRC/TestInstall_is_moved_to_under_ecaljSamples` も trash へ
- `.gitignore` に `/build/`（VSCode の CMake 拡張が最上位に作る）
- `MD/module_map.md`（生成物）と `TOOLS/module_map.py`: 主プログラム → 入口の module、module の階層（226 module、最大 29 段、`m_lmf` が頂点）、
  依存と被依存の数、module の冒頭のコメント。`ecaljclaude.md` から参照。GPU の build の module の循環が 1 つ見つかった（上の TODO）
- `GetSyml/README.md`・`StructureTool/README.md`・`Samples/EPS/EPS_GaAs/README_eps.md` を今の形に書き直した（入力は `ctrlg.<sname>.toml`、
  `getsyml --nobzview`、StructureTool の各スクリプトの向きと出力のファイル名、EPS は `[gw] QforEPS`・`QforEPSau`・`n1n2n3`・`[product_basis] pb_lcutmx`、
  EPS の出力の列）。`getsyml` の使い方の表示も（`-nobzview` → `--nobzview`、`ctrl.nio` → `ctrlg.nio.toml`）
- `ctrlgenToml.py`: nspin=1 のとき原子表の IDU/UH/JH をコメントにして書く（LDA+U は nspin=2 が要る）。`Database/Ce`（nspin=1、idu=12）が
  09-30 の idu の修正（`78475ea3d`）から lmf で止まっていた（`LDA+U must be spin-polarized!`）ので同じく直した（`2cbeb93d5`）
- `Samples/MATERIALS` を仕分けた: 構造のデータベース（62 物質）を `Database/` に展開（`ctrls` と、GW・MLO の節つきの `ctrlg`、`lmchk` で全部読める）。ほかの旧形式の 30 項目は trash（624 ファイル）。残したのは `La2CuO4`・`InAsGaSb`・`BaTiO3`。拾ったノウハウは past_log.md §9
- `README.txt` と `.org` を Markdown に（user「README.txt とあるのは md 形式に。org もそう」）: `StructureTool/README.txt`、`Samples/EPS/EPS_GaAs/README_eps.org`、`Samples/MATERIALS/{LaGaO3_relax,MLOsamples,NiSe_aftest}/README.org`、`Samples/MATERIALS/Si_doping_sample/Memo_bgcharge.org`（`git mv`）。pandoc は `_` を下付きに、`--` をダッシュに変えるので使わず、原文を保つ変換（見出し・`#+TITLE`・`#+begin_src`・コマンド行と設定の断片をコードブロックに）
- `GetSyml/`・`StructureTool/` の古い例（旧形式の `ctrl.*` 86 本、`syml.*`、鉱物の POSCAR 135 本）と使わないスクリプトを trash へ（257 ファイル）。本体（`getsyml`、`vasp2ctrl`・`ctrl2vasp`、`viewvesta`、`refineposcar.py`、`superlattice/`）は残し、動くことを確かめた。past_log.md §12
- `Doxygen/` を trash へ（user「そうしよう」）。コメントを Doxygen 形式に揃えることはしない。作り直し方は past_log.md §11
- `TOOLS/` の古い道具（約 60 項目、632 ファイル）を trash へ。残したのは `samples_tests.sh`・`sync_ecalj_src.sh`・`ozbench/`。中身は past_log.md §10。`diffnum` は試験で今も使うが、使うのは `SRC/exec/pylib/diffnum0.py`（TOOLS の版は古い）
- Claude の個人メモリから、引き継ぐ価値があり今も正しいものを [handover.md](handover.md) に写した（user「メモリの内容はパッケージに入らないので MD/ に。混乱を招くものは良くない」）。
  写す前に 5 点をコードと照らした（rel の既定ブランチは `main`、`master` は無い、など）
- `README.md` を `MD/README.md` へ（user「Claude に読ませて、人間は Claude から情報を取る構造にする。README を人間に読ませるのは好ましくない」）。
  最上位の `README.md` は数行の案内だけ（GitHub の表紙が空にならないように）。`ecaljclaude.md` の方針の文言を直した
- 開発の文書を `MD/` に（`ecaljclaude.md`、`TODOandQuestion.md`、`past_log.md`、研究ログ `research_log.md`、`bf8985687`）
- ecalj の片付けの 3 回目: `PHASE1B_REFACTOR.md`、`HIGHLIGHTS_2026-06_09.md`、`FiniteT_and_QPE_HOWTO.md`、`ecaljdoc_drafts/`、`jobauto/`、
  `SRC/exec_legacy/` を trash へ。中身は past_log.md §3.2・§4.3・§5〜§8 に整理し、README と ecaljdoc（ForDevelopers・kBT・README_tutorial）の参照を直した。
  最上位の `MATERIALS/` は trash ではなく `Samples/MATERIALS/` の下へ移した（user「いったん Samples の下へ」、624 ファイル、`git mv`）
- ecalj の片付け（user「いらないものは trash へ。ノウハウは過去ログへ。重複は整理してから trash へ」）:
  ビルドの生成物 1853 ファイル（`3f0771f2f`）、最上位の打ち込み用と古いスクリプト、`.refactor_notes/`、`SRC/exec/BK`、
  MLOsamples の古い試行、`TestInstall/TESTunused`、`*.bk` などを trash へ。ノウハウは [past_log.md](past_log.md)
- 最上位の `SRC/BK`・`SRC/execgfortran`・`SRC/execAHC`（181 MB）と、GW1500 の古いスクリプトを追跡から外して trash へ（`116254e1d`）
- kt1 の GW1500 の作業ファイルと退避（255 GB）を `~/ecalj/trash` へ（`4f9332d98`）
- ecalj・ecaljdoc の控えのクローンを外付けドライブ `/media/takao/TAKAOMINI`（`ecalj_20261001`・`ecaljdoc_20261001`）に
- Materials Project の API キーを追跡から外した。`<ecalj>/MaterialProject.key`（無視される）と見本 `.example`、`pylib/mpkey.py`（`290397b34`）
- GW1500: 物質ごとの注記の表 `ecalj_auto/gw1500_notes_20261001.tsv`、分類と条件の違い（`GW1500_status.md` §5、`40ffd89e8`）。
  5 月に収束した GOOD はやり直さない（run3 のキューを NOTCONV だけに）。Rb8（mp-1179832）は不正な構造と分かり止めた（`e616f55fb`）
- `gwscconv`: ギャップの無い反復（金属）は固有値の変化で収束を判定（`ecfac2c6b`）
- kt1・kr7（nvfortran GPU）・mic（ifx）で試験の組 inputs・mlo・afsym・install・samples がすべて PASS（HEAD `74ba72dad`）

### 2026-09-30

- `fix-idu10` をマージ（`78475ea3d`）: `idu = 10 + mode` が `sigm` の無いとき二重計数なしの LDA+U になっていた（2023-09 から）。
  MLOsamples の SmP・GdCo5 の LDA+U を作り直した（`e547b81e0`、`a7b19be75`）
- GW1500: FAILED の 147 物質を fp32 で回し直し、144 が収束（`d00b9778e`）
- Samples/Legacy を片付けて組み直し（AtomDimer/N2 など）、ecalj と ecaljdoc を clone し直した
