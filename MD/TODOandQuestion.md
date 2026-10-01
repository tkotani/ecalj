# TODO と質問、やったこと

直すべき点は見つけてもその場では直さず、ここに書く（user 2026-10-01）。メンテナに決めてほしいことも、やったことも、ここに書く。
入口は [CLAUDE.md](../CLAUDE.md)（中身は [ecaljclaude.md](ecaljclaude.md)）。片付けたものの中のノウハウは [past_log.md](past_log.md)。
細かい経緯と数値は研究ログ [research_log.md](research_log.md)。

書き方: 各項目に見つけた日を付ける。済んだら「やったこと」へ移し、済んだ日とコミットを書く。

---

## 1. TODO（直すべき点）

### コード

- **MLO: EH2 s,p を単体の金属（Cu・Ni）に足すと、Γ–X の 2〜3 点の k だけで模型が崩れる**（2026-10-01、`~/work/mlocheck_eh2cat/{Cu,Ni}`、ページ https://claude.ai/artifact/VCWpe5W8Fei6umatXGeZGH の図 4）: その k で MLO の帯が E_F + 0.3 eV に集まり、DFT の帯が抜ける（Cu 0.012 → 0.802 eV）。同じ原子の EH と EH2 がほぼ一次従属になって `Hreduction` の重なりが特異に近い、と疑っている（未確認: その k の重なり行列の固有値を見る）。gwinit の `!` 付きの `mlo_lm2` の行は Cu・Ni には書かれない（既定は基準 1 で EH2 なし）ので、既定には影響しない。EuO（Eu に EH2 s,p）は模型を作る所で止まる: `m_hreduction` の NormalizationCheck でシードのノルムの減りが 1 % を超えた（band 26、−1.3 %。EH と EH2 の両方をシードにすると `zhev_tk4` が一次従属に近い向きを落とすため）。Cu・Ni も同じ原因で、減りが 1 % 未満なので規格化し直して進み、特定の k で崩れているのかもしれない（未確認）
- **`sugw`（`lmf --jobgw=1`）のメモリ**（2026-10-01）: `GEIGpart` が IPW の重なり行列 `ppovl(ngp,ngp)` と LU 用の写し `ppovlLU` を
  各ランクで持つ。1 ランクあたり 32·ngp² バイト、全体は並列数に比例。ngp ≈ V·Q³/6π²（Q = `QpGcut_psi`）なので胞の体積の 2 乗で増える
  （Rb8 の 4410 Å³ で 1 ランク 35 GB）。案: (a) 投入前に ngp から並列数を決める、(b) O は正定値なので Cholesky の因子だけ持つ（半分）、
  (c) O·x を FFT で作り反復法で解く（ngp に比例）
- **`auto_mpquery.py` が今の Materials Project で動かない**（2026-10-01）: MP の API が新しい ID（`mp-aaacpdie` の形、古い番号の 26 進）を返し、
  手元の `mp_api` の検証（pydantic）が止まる。`mp_api` を上げるか、REST を直接読む（requests、`x-api-key`、User-Agent が要る。urllib の既定は 403）
- **GW1500 の選定に構造の確かめを入れる**（2026-10-01）: 別の化合物から副格子を抜き出した MP の項目（ICSD の備考が「〜 part」、
  1 原子の体積が同じ組成の 2 倍以上など）が PBE のギャップ > 0 で選ばれていた（`ecalj_auto/GW1500_status.md` 表 6 の `INVALID_STRUCTURE`）。
  `auto_mpquery.py` で弾くか、印を付ける
- **GPU の build で module が循環している**（2026-10-01、`MD/module_map.md` を作って見つけた）: `m_mpi`（`MPI__AutoSetup` の中で `use m_gpu`）→
  `m_gpu`（`__GPU` のとき `use m_blas`）→ `m_blas`（`use m_gemmul8`）→ `m_gemmul8`（`gemmul8_init` の中で `use m_mpi, only: ipr`）→ `m_mpi`。
  まっさらな状態からは順序が決まらない。いまは `SRC/CMakeLists.txt` が GPU の変種を基本の変種の後に作り、`-I` で基本の `mod/` を見せる回避策で通っている
  （そのコメントは「CMake の走査が m_gemmul8 → m_mpi を見落とす」と書くが、見落としではなく循環）。案: `m_gemmul8` の `ipr` は
  `gemmul8_init` の中で既定の moduli を一度表示する判定にしか使っていない（2026-10-01 に確かめた）ので、`use mpi` で rank 0 を見れば切れる。
  または `MPI__AutoSetup` の GPU メモリの問い合わせを引数で受ける。直したら GPU（kt1・kr7）で clean build と試験。CMake の回避策はその後で外せる
- **AFTEST（`mmtarget.aftest`、例は `Samples/AFfixMMOM`）の残り**（2026-10-01、研究ログ 2026-10-01 朝 06:46 と表 06:46-1）: 場をスピン 1 にしか入れていなかった誤りは直した
  （`6eff2df53`）。残り: (a) モーメントを下げる向きで更新が行き過ぎる（利得 2 固定と 0.1(m − m_t)² の項。NiO 1.28 → 1.0 は 70 反復で、途中で符号が反転）。
  応答を割線で見積もる更新にするか。(b) 対がサイト 1・2、ブロック 1・2 の決め打ち（NiSe だけ 6 ブロック）。`AF=` の印や種の名前から対を決める。
  (c) `m_ldau_init` が lmf の起動のたびに場を一歩更新して `mmagfield.aftest` を書き直す（`job_band` でも）。(d) ecaljdoc の `UsageDetailed.md` の「直す必要がある」を、
  使い方（LDA+U のブロックが要る、U = 0 でよい、`SYMGRPAF`、目標のモーメントの定義）に書き直す

- **`hgw` の残り**（2026-04 の統合の残タスク、2026-10-01 に past_log.md §3.2 から）: ノード内の W の共有を `MPI_Win_allocate_shared` で（メモリの重複を減らす）。
  `hsfp0_sc` の Sx（`--job=1`）・core の交換（`--job=3`）も `hgw` に入れる（時間は小さいので優先度は低い）

- **`m_tetrakbt.f90` の使われていないルーチン**（2026-10-01）: 中点分解の `tetrakbt`・`eaf_triangle`・`eafww`・`factri`・`factri0`・`integral_1t`・
  `integral_0t`・`funcgx`・`numintall`・`funcgx2`・`ndiv_tt`・`int_simpson`（約 300 行）はどこからも呼ばれない（`tetwt5.f90` は `use` の only に
  `tetrakbt` を挙げるだけ）。今も使うのは `tetrakbt_init`・`kbt`・`integtetn`。2026-06 に「参照のため残す」としたもの。消すなら kBT の試験で確かめる（past_log.md §13）

- **MLO の模型の残り**（2026-10-01。既定は基準 1・2・3 の形に決まった: ecaljdoc mlo §1・§9、`MD/handover.md` §5）:
  (a) 基準 3（空格子球）の自動化。案は空隙の半径（最も近い MT 球の表面までの距離）が 3 a.u. を超える所に置く `SRC/exec/ctrlg_addes.py`（未コミット・未検証）。
  Bi₂Te₃ の空隙の半径は、手早い見積もり 3.56 a.u. と道具の 2.61 a.u. が食い違ったまま（菱面体の長い c で、像の数 ±1 と格子の点の数が足りないのを疑う）。
  SiO₂ は手で置いて 0.001 eV（研究ログ 2026-10-01）。
  (e) 検査（`mlo_bandcheck.py` の CHECK）で、基準 1 の Sn・AlSb・Bi₂Te₃・InSb・Cu が 1 点だけ 0.1 eV を少し超えて FAIL（rms は 0.01〜0.03）。
  線形独立性の崩れではない（重なりの最小固有値は 3e-3 程度）。帯の交差の付近の形が違う。直すか、目安を見直すか
- **`job_mlo_soc` が空の `band_MLO_spin2.dat` を書く**（2026-10-01）: 2N のスピノルの帯は `band_MLO_spin1.dat` に入る。空の spin2 があると
  `mlo_bandplot.py` が空の枠を描いていた（描く側は 2026-10-01 に空のファイルを飛ばすようにした）。`bandplot_MLO.isp2.glt` も残る。
  `job_mlo_soc` は始めに古い spin2 を消すが、`mlo` がまた空のものを書く（`mlo_bandcheck.py` は SOC のとき spin1 だけを見る）
- **`Samples/AFsymmetry/NiO` は `pwmode = 1` で `symgrpaf`**（2026-10-01）: k 点を 4³ にすると `rotwave: q+G rotation error (We have to set PWmode=11 for symgrpAF)`
  で止まる。試験の 3³ では通っているだけ。`pwmode = 11` にして参照を作り直すか

### 試験と入力

- `MLOsamples/RuO2` の保存してある `rst`・`dmats` は、`pwmode = 11` の LDA+U の誤り（2026-03-30〜09-30、`2498e5283` で修正）の時期に作ったもの。
  作り直すか（GdCo5・SmP は 2026-09-30 に作り直した）
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
- push: dev・rel とも、t14 の main より 674 コミット遅れ（2026-10-01）

## 3. 実行中

- kt1 `/mnt/data1/gw1500_rerun/run3`: GW1500 の NOTCONV の残り（fp32、5 本）。2026-10-01 04:40 に約 17 時間の見積もり

---

## 4. やったこと（新しい順）

### 2026-10-01

- MLO の模型の検査（`mlo_bandcheck.py` の CHECK、`job_mlo`・`job_mlo_soc` が最後に回す、`mlo` が `MLO_ovlpmin.dat` を書く）、
  窓の基準の E_F を `--efermi=` のファイルから取る直し、`Samples/MATERIALS` の 63 の ctrlg の `[mlo]` を今の gwinit で書き直し（研究ログ 2026-10-01 21 時）。
  試験: mlo 45、MLO-QSGW 5、install 64、inputs 176 が t14 で PASSED
- MLO の模型を基準 1・2・3 に整理した（user と決めた。ecaljdoc mlo §1・§4・§9、`MD/handover.md` §5、研究ログ 2026-10-01 16〜20 時）:
  半内殻の局所軌道は帯の上端で自動（E_F − 8 eV より上は EH と入れ替え、−17〜−8 eV は加える、`d759110e4`）、gwinit は f と `!` 付きの `mlo_lm2`（基準 2）を書く、
  MLO-QSGW の凍結した模型の本数の誤り（`2053bb5ea`）、`mlo_bandcheck.py` の窓を [VBM − 8, CBM + mlo_delta] に（`bd1b893ec`）、
  `job_mlo_soc` がスピン軌道ありの DFT のバンドも描く（GaAsSoc の −0.1 eV は比べ方の誤り、`dc66439c5`）。Samples/MATERIALS の結果のページ
  https://claude.ai/artifact/VCWpe5W8Fei6umatXGeZGH（`Samples/MATERIALS/mlocheck/gallery.py`）。試験: inputs 176、mlo 25 試料、MLO-QSGW、Fe_mlo_magnon が t14 で PASSED
- AFTEST の修正（`aftest-fix`）を main にマージ（`6eff2df53`、user「マージできるよね」→「良い方を壊さないように」）。AF の対称性ありの計算は ehf・sev・モーメント・場がマージ前と同じで、ehk だけが制約の下の全エネルギーになった。afsym 4・install 64 が PASS。例を `Samples/AFfixMMOM`（user の命名、`NiO_afsym`・`NiO_noafsym`、試験の組 `affix` 12 件 PASS）に置いた（`ca1e4c7fa`）
- ecaljdoc の TOML の流れの最短の手順は、`manual/README_tutorial.md` の GetStarted（Step 0〜6: POSCAR → ctrls → `ctrlgenToml.py` → `lmfa`・`lmf` → バンド → `gwsc`、Step 2-Migration、`--ctrlg:` の上書き）に既にあった。TODO から外した。ecaljdoc `manual/mlo.md` に既定の模型が外れる 2 つの型と表 M1 を書いた（`422dc92`、未 push）
- `Samples/MATERIALS` の 65 物質で LDA と MLO の自動の模型を回した（InAs/GaSb n10 の 40 原子は user の判断で回さない）。結果の表は `Samples/MATERIALS/MLOcheck_20261001*.tsv`、一覧のページ https://claude.ai/artifact/VCWpe5W8Fei6umatXGeZGH、研究ログ 2026-10-01 朝 07:16
- `job_mlo`: ctrlg の `[ham] so = 1` なら `job_mlo_soc` を案内して止まる（`mlo` がハミルトニアンの NaN で止まっていた。GaAs_so）。`--ctrlg:ham.so=` の上書きは尊重。`mlo_bandplot.py` は空の spin2 を描かない
- `--cls`: `m_clsmode_finalize` に `ndimh`（`m_igv2x` の最後の k 点の値）でなく `nbandmx` を渡す（`vcdmel` の重みの並びと `dostet` の読み方を揃える）。CrN（`Samples/TestInstall/crn`）で APW なし・pwemax 2 と 5・4×4×4（ndimh が 182〜188 と変わる）のどれも `dos-vcdmel.crn` が一致（ずれるのは DOS の窓より上の帯だけだった）。ブランチ `cls-nbandmx` をマージ（`7ca6f1caa`）。試験は `~/work/clstest_20261001`
- `InstallAll.py`: build の後に、bindir の中で「この ecalj の木を指していて行き先の無いリンク」だけを消す（`remove_dangling_links`）。`SRC/exec` の行き先の無いリンク（消した `SRC/exec/build/` を指す 12 本）は張らず、trash に移した。これが `~/bin/mlo` などを一時的に死んだリンクで上書きしていた（CMake の deliver が build の後に戻していた）。SRC/exec のエディタの一時ファイル（`~` で終わる、`#`・`.#` で始まる）はリンクしない。t14 の `~/bin` の 58 本は次のインストールで消える（一時の bindir で試験）
- ecaljdoc `manual/spectrum.md` に、メッシュの外の q（`QforGW`）の注意 3 点（`EMAXforGW` が必須、窓を変えたら交換から、`epsWVR` の行）を書いた（past_log.md §5 から、ecaljdoc の未 push のコミット）
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
