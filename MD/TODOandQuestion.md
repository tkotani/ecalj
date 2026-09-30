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
- **module の依存と主プログラムの流れの一覧を `MD/` に置く**（2026-10-01、Doxygen をやめた代わり）: `use m_foo, only:` を拾って、module の DAG と、主プログラム（`SRC/main/*.f90`）から各 module への流れを機械的に書き出す小さなスクリプト。人が大局をつかむ入口として、Claude が説明に使う
- **VSCode の CMake 拡張が最上位に `build/` を作る**（2026-10-01）: `.vscode/settings.json` の設定か、`.gitignore` に `/build/` を入れる

- **`hgw` の残り**（2026-04 の統合の残タスク、2026-10-01 に past_log.md §3.2 から）: ノード内の W の共有を `MPI_Win_allocate_shared` で（メモリの重複を減らす）。
  `hsfp0_sc` の Sx（`--job=1`）・core の交換（`--job=3`）も `hgw` に入れる（時間は小さいので優先度は低い）

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
- **`Samples/MATERIALS` の旧形式の古いサンプル**（2026-10-01 に最上位から移した 30 項目）を、今の形にするか、trash に移すか。
  中身は past_log.md §9。ecaljdoc の `README_tutorial.md` の `jobmaterials.py` の節もこれを前提にした古い説明
- **GW1500 の `INVALID_STRUCTURE`（12）・`SUSPECT_STRUCTURE`（4）を集合から外すか**。いまは注記だけ
- **5 月（TF32）と fp32 の差**: kr7 の比較（同じバイナリ・同じ入力で fp32 と TF32）の結果しだいで、GOOD の 1120 を見直すか
- **ビルドの生成物が入ったコミット `a0c7a7300`（78 MB）を、push の前に履歴から消すか**。消すと以後 674 コミットのハッシュが変わり、
  研究ログなどに書いたハッシュが合わなくなる。消さなければ、公開のリポジトリの履歴に 78 MB が残る
- **Materials Project の API キーの作り直し**（user の作業）: 2025-04〜2026-09 の公開リポジトリの履歴に入っていた
- **trash を空にするか**: t14 の `~/ecalj/trash`（旧ツリー 58 GB ほか）、kt1 の `~/ecalj/trash`（255 GB）。kt1 は空にするまで `/` の空きが 49 GB のまま
- 試験に使った古いツリー（mic の `~/ecalj_test0928`、kt1 の `/mnt/data1/ecalj_test0930b`・`0930c`）も trash に入れるか
- ブランチ `fix-idu10`（main にマージ済み）を消すか
- push: dev・rel とも、t14 の main より 674 コミット遅れ（2026-10-01）

## 3. 実行中

- kt1 `/mnt/data1/gw1500_rerun/run3`: GW1500 の NOTCONV の残り（fp32、5 本）。2026-10-01 04:40 に約 17 時間の見積もり
- kr7 `~/gw1500ab`: fp32 と TF32 の比較（8 物質）

---

## 4. やったこと（新しい順）

### 2026-10-01

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
