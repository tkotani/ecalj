# ecalj Claude Development Guide

ecalj を Claude (LLM) で開発する際のガイドライン。
コーディング規約、ビルドの罠、GPU 開発の教訓、アーキテクチャ方針をまとめる。

## 文書の地図（まずここから。2026-10-02）

[`CLAUDE.md`](../CLAUDE.md)（この文書を読み込む）から次の順にたどる。人が MD/ を開いたときの入口は [README.md](README.md)（読む順だけ）。パスはリポジトリの最上位から。リンクの切れは `python3 TOOLS/doclinks.py` で確かめる。

| 何を知りたいか | 文書 |
| --- | --- |
| いまの状況、記憶の無いセッションへの引き継ぎ（計算機、取り決め、落とし穴） | [handover.md](handover.md) |
| まだやっていないこと（やりかけ・未着手・判断待ち）と、やったこと | [TODOandQuestion.md](TODOandQuestion.md) |
| 経緯、計測、判断の時刻（最新が上） | [research_log.md](research_log.md) |
| 片付けたものに入っていたノウハウ、trash に移したものの表 | [past_log.md](past_log.md) |
| 新機能と変更（利用者向け、正本） | `Changes.txt`（最上位） |
| 主プログラム → module の階層（生成物） | [module_map.md](module_map.md) |
| 個別の設計・手順: 対称性 [symmetry_spglib.md](symmetry_spglib.md)、Wannier と MLO の比較 [wannier_vs_mlo.md](wannier_vs_mlo.md)、ecaljdoc の管理と公開 [ecaljdoc_publish.md](ecaljdoc_publish.md)、履歴の書き換えのハッシュの対応 `commit_map_20261002.txt` | [`MD/`](.) |
| 使い方と理論の本体（公開サイト https://ecalj.github.io/ecaljdoc/ の元） | [`ecaljdoc/manual/`](../ecaljdoc/manual)（入口 [../index.md](../ecaljdoc/index.md)、サンプルの一覧 [../manual/samples.md](../ecaljdoc/manual/samples.md)）、[`ecaljdoc/theory/`](../ecaljdoc/theory) |
| 変更の記録（利用者向け） | `Changes.txt` |
| サンプルと試験 | `Samples/<dir>/README.md`（各 README の冒頭に本体の節へのリンク）、試験の組は [`TOOLS/samples_tests.sh`](../TOOLS/samples_tests.sh) |
| 標準の流れに入れない試作の道具 | [`TOOLS/gadget/README.md`](../TOOLS/gadget/README.md) |

**開発者向けの文書を [`MD/`](.) に集めた**（2026-10-02 20:59、user「開発者用で ecaljdoc からたどらなくてもいいもの（余計に混乱を招くもの）は ecaljdoc/MD へ」。
ファイルはまとめずにそのまま移した。中身の重複と整合の整理はこの後）。元の場所は次のとおり（`git log --follow` で前の履歴もたどれる）:

| いまの場所（[`MD/`](.) の下） | 元の場所 | 中身 |
| --- | --- | --- |
| [ForDevelopers.md](ForDevelopers.md) | [`ecaljdoc/manual/`](../ecaljdoc/manual)（サイトにあった） | 開発の引き継ぎの手引き: リポジトリと push、ビルド、試験、GPU、ジョブの投入、研究の現在地 |
| [ForDevelopers_research.md](ForDevelopers_research.md) | [`ecaljdoc/manual/`](../ecaljdoc/manual) | 研究ログの要約と日付の索引 |
| [developer.md](developer.md) | [`ecaljdoc/manual/`](../ecaljdoc/manual) | 試験の仕組み `testecalj`（2025-10） |
| [auto.md](auto.md) | [`ecaljdoc/manual/`](../ecaljdoc/manual) | `ecalj_auto`（GW1500 の自動の流れ） |
| [kBT_history.md](kBT_history.md) | [`ecaljdoc/manual/`](../ecaljdoc/manual) | 有限温度の 2026-06〜09 の経緯と誤り（今の仕様は [`ecaljdoc/manual/kBT.md`](../ecaljdoc/manual/kBT.md)） |
| [mlo_backup.md](mlo_backup.md) | [`ecaljdoc/manual/`](../ecaljdoc/manual) | MLO の経緯と作業記録（今の仕様は [`ecaljdoc/manual/mlo.md`](../ecaljdoc/manual/mlo.md)） |
| [`implementation/`](implementation/README.md)（hsfp0、hx0fp、issue） | [`ecaljdoc/implementation/`](../ecaljdoc/implementation) | 自己エネルギーと W の実装の覚え書き |
| [MemoCode.md](MemoCode.md) | [`ecaljdoc/theory/`](../ecaljdoc/theory) | m_zmel などのコードの覚え書き |
| [`kBT/`](kBT/README.md)（sigma_mlo_design、gpu_fp32_plan、gpu_fp32_report、finiteT_202606、LiTi2O4_finiteT_202606。案内は [kBT/README.md](kBT/README.md)） | [`Samples/kBT/`](../Samples/kBT/README.md)、[`Samples/kBT/LiTi2O4/`](../Samples/kBT/LiTi2O4/README.md)（finiteT は `README_202606_finiteT.md` だった。six_patterns は 2026-10-02 22:03 にサンプルのそば [`Samples/kBT/LiTi2O4/six_patterns.md`](../Samples/kBT/LiTi2O4/six_patterns.md) に戻した。図と表がそこにあり、README から上のディレクトリへのリンクが VSCode で開けなかったため） | MLO-QSGW の設計書、GW の GPU 高速化の計画と報告、有限温度と LiTi₂O₄ の計算の記録 |
| `mlo_notes/`（[README.md](mlo_notes/README.md) ほか、図も） | [`Samples/MLOsamples/BackUp_notes/`](../Samples/MLOsamples/BackUp_notes)、`Samples/MLOsamples/README_SOC.md`（→ `mlo_soc_memo.md`） | MLO の理論の整理、自動窓の探索の報告、最適化の覚え書き、nskip の問題、MP の比較、SOC の覚え書き |
| [testecalj_2025.md](testecalj_2025.md) | `Samples/TestInstall/README_testecalj.md` | 試験の仕組み（2025-10） |
| [GW1500_failures.md](GW1500_failures.md) | [`Samples/mptf32problem/`](../Samples/mptf32problem/README.md) | GW1500 の失敗と fp32 での回し直しの記録 |

## 記録の方針: コードのコメント、コミット、研究ログ、文書（2026-09-27）

Claude（や人）が後でコードや記録を読むときの手がかりになるように、何をどこに書くかを決めておく。
**この節が正本。** 開発者の Claude の個人メモリ（`~/.claude/projects/<project>/memory/feedback_*.md`）に
同じ方針を写してある場合はそちらを要約とし、方針を変えるときはまずこの節を直してからメモリを合わせる。

### コードのコメント

- 書くこと: 何を計算するか（式・添字・配列の並び）、なぜその形か（正しさの条件・不変量）、壊さないための注意
- 変更・判断・計測・バグ修正の注記には**日付時刻**を付ける（例: `(2026-09-27 21:26)`）。
  日付があれば、古くなった注記も「その時点の記録」として読める
- 経緯や計測値は、設計の理由になるものだけ短く（例: `(2026-09-27, RTX 5090: 2.1 ms as one product, 0.46 ms split)`）。
  すぐ古くなって**将来混乱を招くもの**（性能値の羅列、他のルーチンの内部動作の説明）は書かない
- **バグ修正は書いてよい**: 何が起きていたか、いつ直したか
  （例: `Bug fixed 2026-09-27 16:07: the partial DOS depended on the number of ranks`）。後戻りを防ぐ
- 式の定義のような時間で変わらない説明には日付は要らない
- コードの中のコメントは英語。詳しい経緯は研究ログへ

### 整合性を保つためのメモ

- コードを直したら、近くのコメント（特に構造や流れを説明するもの）を読み直す。誤りになったものは直すか、日付付きで追記する。
  古い日付の注記を最新に合わせ続ける必要はない（日付で古さが分かる）。**誤りだけは残さない**
- ルーチン・配列・ファイルの名前を消したり変えたりしたら、`grep` でコメントや文書の中の言及も探して直す
- コメントに行番号（`mkjp.f90:296` など）を書かない。すぐずれる。ルーチン名やブロック名で指す
- 「使われていない」と判断したコードは、呼び出し元を `grep` で確かめてから消す（別の流れ、例えば job_mloW・job_mlo_magnon だけが使う段がある）
- AI エージェントにレビューさせた指摘は、必ずコードで確かめてから反映する。誤りが混じる
  （2026-09-27: 「npm=2 で χ0 が壊れる」と報告されたが、実際は x0kf_zxq の冒頭で止まっていた）

### コミットと文書の言語

- git のコミットメッセージは英語（ecalj、ecaljdoc とも）。`Changes.txt`、`MD/*.md`、`manual/*.md` の中身は日本語でよい
- コミットはこまめに、論理的に分ける。使い方が変わるときは `Changes.txt` と ecaljdoc も直す

### 研究ログ（[research_log.md](research_log.md)）

- 最新を上に。日付ごとに `## YYYY-MM-DD — 要約`、その中は `### HH:MM 見出し` を新しい順に積む
- 投入・完了・停止・判断（誰の判断か）の時刻を書く。古い仮説は消さず、後の判定を添える
- 経緯・失敗談・細かい計測はここに書く（コードの注記や利用者向け文書には要点と日付だけ）
- 時刻は推さない。`date` の出力を変数で差し込んで書く（`T=$(date +%H:%M) python3 - <<'EOF'` など）。「05:3x」のようにぼかした時刻も書かない
  （2026-10-02: 一晩で十数回、見ずに書いた時刻を直した）

### TODO・質問・やったこと、過去ログ、片付け（2026-10-01、user の指示）

- 直すべき点は見つけてもその場では直さず、[TODOandQuestion.md](TODOandQuestion.md) に書く。メンテナに決めてほしいこと、実行中のこと、
  やったこと（日付とコミット）も同じファイルに書く
- いらないものは消さずに各計算機の [`ecalj/trash/`](../trash) に移す（リポジトリのものは `git rm --cached` で追跡も外す。`/trash/` は `.gitignore`）。
  trash は適宜減らしていく（user 2026-10-02。作り直せるもの、要点を past_log に写したものから消す。消したものは past_log の表 1 に書く）
- 移すものの中にあるノウハウや役に立ちそうな情報は、移す前に読んで [past_log.md](past_log.md) に話題ごとに整理して書く。
  情報が重複しているものは、まとめてから trash へ。past_log.md の表 1 に、元の場所と外す前のコミットを書く（`git show <コミット>:<パス>` で取り出せる）

### 図と数値（2026-09-28）

- 図には、描いた数値をファイルで添える。図を作るスクリプトが npz や dat を図の隣に書き、元のファイルの場所・日時も入れる
  （例: [`Samples/kBT/LiTi2O4/mlo_rows.py`](../Samples/kBT/LiTi2O4/mlo_rows.py) の `DATAOUT`）。手元の作業場所（/tmp など）は消える
- 反復ごとの結果（バンドなど）は、その反復のうちに描いてファイルに残す。次の反復で上書きされる作業ファイル（QMLO_SigRs など）に頼らない
  （2026-09-28: `gwsc 10` で回した tf32 の反復 1〜9 の MLO バンドは、あとから作れなかった）

### 利用者向け文書（ecaljdoc、各 README）

- 必要な理由と注意だけ。どこで迷ったかの経緯（苦労話）は書かない（研究ログへ）
- 数式には必ず式番号（`$$ ... \tag{1}$$`）を付け、本文から番号で参照する。図表にも番号を付ける
- 説明の軸は ecaljdoc。ecalj の Changes.txt・各サンプルの README には要点と ecaljdoc への参照だけを書き、同じ説明を二か所に書かない。
  ecaljdoc の正本は 2026-10-02 から [`ecalj/ecaljdoc/`](../ecaljdoc/README.md)（公開は GitHub の `ecalj/ecaljdoc` へ `TOOLS/publish_ecaljdoc.sh` で写して出す、手順は [ecaljdoc_publish.md](ecaljdoc_publish.md)）。
  古いクローン `~/ecaljdoc` では作業しない
- 文書の置き場（2026-10-02、user）:
  - `ecaljdoc/`（公開サイト）: 使い方と理論の本体。サイトは `ecaljdoc/` だけから作られる
  - [`MD/`](.)（最上位）: 開発の記録（Claude が読む）。公開しない（2026-10-02、user「公開リポジトリに MD は送らない」。一度 `ecaljdoc/MD/` に置いたが同じ日に戻した）
  - `Samples/<dir>/README.md`: そのサンプルのそばに要るもの（何の例か、回し方、試験が比べるもの、そのサンプルの数値と図）。冒頭に本体の節へのリンクを置く
  - リンクの向き: ecaljdoc の本体（manual・theory）から外（Samples、ソース）へは GitHub の URL（`https://github.com/tkotani/ecalj/tree/main/...`。サイトで
    `../` の相対リンクは切れる）。外から ecaljdoc へはリポジトリの中の相対パス（[`../../ecaljdoc/manual/mlo.md`](../ecaljdoc/manual/mlo.md)）
  - ecaljdoc の本文で `<...>` や `{{` をコードの印（バッククォート）の外に書かない。VitePress のビルドが止まる（2026-10-02 に 3 か所あった）
  - `MD/` の文書の中のファイルのパスは、開けるように相対リンクで書く（`[`Samples/.../README.md`](../../Samples/.../README.md)`）。
    コードの印で書いたパスは `python3 TOOLS/mdlinkify.py <md>` がリンクに直す（実在するものだけ。生成物の `module_map.md` には使わない。2026-10-02、user）
  - 文書を直したら `python3 TOOLS/doclinks.py` で確かめる（切れた相対リンク、サイトから `../` で出るリンク、本体へのリンクの無い Samples の README、
    入口から名指しされない MD のファイル。2026-10-02）
- ecalj は「Claude が読み込んで人間に説明できる」パッケージにする。人間は Claude を通して情報を取る（user 2026-10-01）。最上位は [`CLAUDE.md`](../CLAUDE.md)（入口）と数行の [`README.md`](../README.md)、`Changes.txt`。
  開発の文書は [`MD/`](.)（新機能の正本は `Changes.txt`。2026-05〜09 の英語のログとクイックスタートだった `MD/README.md` は 2026-10-02 に trash へ）。
  コードの大局（主プログラム → 入口の module、module の DAG と階層、各 module の依存）は [module_map.md](module_map.md)。
  `python3 TOOLS/module_map.py` で作り直す生成物で、手では直さない（2026-10-01、Doxygen の代わり）
- 引き継ぎ: Claude の個人メモリ（`~/.claude/.../memory/`）はパッケージに入らないので、別の機械や記憶の無いセッションに要ること
  （user との取り決め、計算機の癖、ビルドと実行の落とし穴）は [handover.md](handover.md) にも書く。その時点の状況や古くなったことは写さない（混乱の元）
- [`MD/`](.) は基本的に Claude が読むもの（user 2026-10-01）。人に読ませる体裁より、Claude があとで正確に引けることを優先する:
  日付時刻、コミット、ファイルのパス、数値、式番号・表番号を省かない。人への説明は Claude がそこから組み立てる
  （メンテしにくくなる。2026-09-28）。ecaljdoc の中でも、同じ式や表を二度書かずに番号で参照する
- 理論は、式とアルゴリズムが追えるように正確に書く: 何を計算するか、離散化と規格化、その結果をどこで使うか。
  キーの名前や定義を変えたら、変遷の表（期間ごとのキーと意味）を残す。古い入力や記録を読むときに混乱しないように（2026-09-28）

### 計算を回す前の確認

- リモート機のソースが手元と同じか（`SRC/.ecalj_rev`、ファイルのチェックサム）
- 新しい機能が `.so` に入っているか（`nm -D` や `strings libecaljF*.so`。実行ファイルは薄い wrapper）。GEMMul8 の包み（`gemmul8_wrapper.cu`）は
  別の `libgemmul8wrap.so` に入る（2026-09-29: libecaljF だけを見て「入っていない」と取り違えた）
- インストール先の `gwsc` などが symlink か実体コピーか（実体だと古いスクリプトのまま走る）
- 走り始めたらログに期待する印が出ているか確かめる
- `testecalj` の前に `*_work` を消す（前回の失敗の rst が残ると偽の失敗になった。いまは testecalj も走らせる各ターゲットの `_work` を
  作り直す〔2025-10-09 から〕ので念のため）
- 試験の最中は、`~/bin` が symlink で指すもの（[`SRC/exec`](../SRC/exec) のスクリプト、`libecaljF.so`）を変えない。別の木でビルドする
  （2026-10-02: 試験の途中で `job_mloW` を書き換え、後半の組が別の版で走って偽の FAIL になった）
- コミットはファイルを名指しで（`git commit -a` を使わない。2026-10-02: 確かめる前の参照のファイルを別のコミットに混ぜた）
- `testecalj` の件数は最後の要約だけで数える。ターゲットごとにそれまでの要約を出し直すので、ログ全体の PASSED 行は延べ数
  （2026-09-28: `--all` の 64 件を 832 件と数えていた）。Samples の組ごとの試験は [`TOOLS/samples_tests.sh`](../TOOLS/samples_tests.sh)

### Claude の作業の進め方

- 着手前に、何をどうするかを一〜三行で述べる（途中で方針を変えられるように）
- commit はその場でしてよい。公開リポジトリへの push はメンテナの指示があってから（ecaljdoc の公開は ecalj の push と同じときに
  `TOOLS/publish_ecaljdoc.sh --push`。[ecaljdoc_publish.md](ecaljdoc_publish.md) §3）

## コーディング規約: Fortran Singleton Pattern (kotani 設計)

ecalj は **Fortran module を singleton class として使う** 設計パターンを採用している。
長年の Fortran 開発で type (derived type) ベースの設計が引き起こす問題
(データ多重化、サイズ不整合、OpenACC 非互換、可読性劣化) と格闘した経験から、
module の言語特性を活かした引き算の設計に到達したもの。

Fortran 2003/2008 の機能は選択的に取り込む:
- **使う**: `allocatable`, `module`, `protected`, `use only`, `associate`, `block`
- **使わない**: `class`, `type-bound procedure`, `select type`, `inheritance`

判断基準は「Fortran の配列計算言語としての強みを活かすか、OOP の後付け機能で
複雑さを持ち込むか」。前者を取り、後者を捨てる。

### なぜ singleton か

Fortran の module は言語仕様上、プログラム内で唯一のインスタンスを持つ。
これは OOP の singleton と同じ性質であり、科学計算の多くのサブシステム
(ハミルトニアン、波動関数ストレージ、自己エネルギー計算器など) は
複数インスタンスが不要なため、module = singleton が自然に適合する。

C++ や Python では singleton pattern を「わざわざ」実装する必要があるが、
Fortran では module 自体がそれを提供する。ecalj はこの言語特性を
意識的・体系的に活用している。

### 設計原則
- **1 module = 1 責務 = 1 状態セット**: 複数インスタンス不要なら `type` を被せない
- **状態は module 変数で公開**: `protected, public` で読み取り公開、書き換えは module 内 subroutine 経由のみ
- **import は明示**: 必ず `use m_foo, only: bar, baz` で必要な symbol だけ取る
- **subroutine 規約**: 入力=引数、出力=module 変数で返す
- **`type` (derived type) は積極的に避ける**: Fortran の type には以下の罠がある
  - `intent(out)` で allocatable components が暗黙に deallocate される
  - `present(state%alloc_member)` が nvfortran + OpenACC で silent fail
  - 合成 type の copy semantics が不明瞭
  - `use only` でデータフローが見えない (type の中身は呼出側から不透明)

### 階層化設計: module が module を use する DAG

singleton module 間のデータ受け渡しは `use` の階層で表現する。
各 module は自分より下位の module だけを `use` し、循環依存を禁止する。
これにより依存関係が DAG (有向非巡回グラフ) になり、データの流れが一方向に追える。

```
m_struct_from_lmf (構造: plat, alat, pos, z)
    ↑ use
m_gw_product_basis (積基底: ndima, nlnx)
    ↑ use
m_sxcf_sc (自己エネルギー計算)
    ↑ use
main_hrcxq / hgw_combined (メインプログラム)
```

**データフローの規約:**
- **大量データの入力**: `use m_foo, only: bar` で下位 module の `protected` 変数を読む。subroutine 引数では渡さない
- **subroutine 引数**: フェルミエネルギー `ef`、q 点 `qp`、スピン `isp`、反復回数 `niter` など、**指示的・局所的な値に限る**。配列データを引数で渡し回さない
- **出力**: subroutine が計算した結果は module 変数に格納し `protected` で公開。caller は `use only` で読む

この規約により:
- `type(state)` を引数で渡し回す必要がない (type を避けられる)
- どの module がどのデータを提供し、誰が消費するかが `use only` で完全に見える
- `grep 'use m_foo' SRC/subroutines/*.f90` で依存関係を即座にトレースできる
- module の差し替え (例: FILE → MEMORY backend) が caller に影響しない

### なぜ Fortran OOP (type-bound procedure, class) を使わないか

Fortran 2003 で `class`, `type-bound procedure`, `select type`, inheritance が導入されたが、
これらは C++/Java のために設計された OOP 概念を配列計算特化言語に後付けしたもので、
実用上の問題が多い:

- `select type` による polymorphism は冗長で可読性が低い
- `type-bound procedure` は module procedure への syntax sugar に過ぎない
- inheritance は Fortran のメモリモデル (allocatable component) と相性が悪い
- コンパイラのバグが多い (特に nvfortran + OpenACC + derived type の組み合わせ)

ecalj の singleton + 階層 use は、Fortran の元々の強み (module, allocatable, 配列操作) だけで
設計を完結させる **言語に逆らわない** アプローチ。

### singleton の利点

1. **OpenACC/GPU との相性**: device 変数の管理が単純。module 変数は全 subroutine から直接アクセスでき、`!$acc` ディレクティブで素直に扱える。type のメンバを device に置くと nvfortran の制約に次々ぶつかる
2. **データフローが grep で追える**: `use m_foo, only: bar` を grep すれば誰が何を使っているか一目瞭然。type だと `state%bar` を追う必要があり、複数の type が絡むとトレースが困難
3. **re-entrant 化が楽**: module 変数を dealloc-realloc すればいい。type のネストした copy semantics を考えなくていい
4. **backend 切替が caller 透過**: module 内部で FILE / MEMORY_3D 等の実装を切り替えても、caller は同じ module 変数を `use` するだけ
5. **LLM にも優しい**: module のヘッダ (変数宣言部) を見れば状態が全部わかる。type の定義を辿って中身を理解する必要がない

### 核心的利点: 変数宣言の唯一性

従来の subroutine 引数渡し方式では、同じ変数を caller と callee の両方で宣言する
必要がある (二重記述):

```fortran
! 呼ぶ側
call calc_sigma(ndim, nqibz, nspinmx, qibz, wk, ekc, ...)

! 呼ばれる側 — 同じ変数を再宣言、サイズも再指定
subroutine calc_sigma(ndim, nqibz, nspinmx, qibz, wk, ekc, ...)
  integer, intent(in) :: ndim, nqibz, nspinmx
  real(8), intent(in) :: qibz(3,nqibz), wk(nqibz), ekc(ndim)
```

引数が 10〜20 個になると caller と callee の対応がずれてバグになる。
配列サイズの整合性も人間が保証しなければならない。二重記述は百害あって一利なし。

singleton + `use only` なら変数宣言は module 内の一箇所だけ:

```fortran
! 呼ばれる側 — use で取るだけ、再宣言不要
subroutine calc_sigma(ef, qp, isp)  ! 指示値のみ引数
  use m_read_bzdata, only: qibz, wk, nqibz
  use m_struct_from_lmf, only: nspinmx
```

- 変数の生成元が唯一 — 二重記述ゼロ
- サイズ不整合が原理的に起きない
- 引数リストが短い (指示値のみ) ので caller-callee の対応ミスが起きない

**深い呼び出し階層での追加データ問題**: 5段の call 階層の末端で `natom` が1つ必要になった場合:
- **引数渡し**: 途中の全 subroutine (5箇所) の引数リストに `natom` を追加。大量の変更、リグレッションリスク
- **type**: `state%natom` として渡すと、データが本来の出自 (生成元 module) から切り離されコピーが増える。実際に起きるのは、たまたま 5 段目まで通っている type に `nbas` を便乗で放り込んでしまうこと。やがて `crystal_state%nbas`、`gw_state%nbas`、`basis_state%nbas` と複数の type に同じ値のコピーが散在し、どれが正でいつ同期されるか誰にもわからなくなる — **データの唯一性が簡単に消失する**
- **singleton use**: `use m_struct_from_lmf, only: natom` を末端の1箇所に追加するだけ。途中の subroutine は変更不要、データの出自は `m_struct_from_lmf` と明確

### 弱点と対策

- **subroutine の出力がシグネチャに現れない**: 出力は module 変数に書かれるため、関数定義だけ見ても「何を返すか」がわからない。ただし caller 側の `use m_foo, only: bar` を見れば出力は特定でき、LLM なら module 宣言部と caller の use only を同時に見て追跡できるため、実質的な問題は小さい
- **module 状態の一括クリーンが面倒**: type なら `deallocate(state)` で一発だが、module 変数は個別に deallocate する subroutine が必要。ただし LLM は module ヘッダの変数宣言を見て漏れなく dealloc subroutine を書けるため、LLM 時代にはこの弱点の実害が縮小している

### 実例: WB.3 リファクタ

`type(wv_storage)` / `type(sxcf_state)` を捨て、module-level singleton 化した結果:
- **行数減**: type 定義、constructor、accessor が不要に
- **OpenACC 制約回避**: structured data region 内の BLOCK 禁止等を自然に回避
- **FILE ↔ MEMORY_3D backend 切替が caller 透過に**: m_wv_storage 内部のフラグ切替だけで、m_sxcf_sc や main_hrcxq は変更不要
- **hgw_combined 統合が容易に**: hrcxq + hsfp0_sc の状態が module 変数として共有できるため、2つのプログラムを1プロセスに統合する際の状態受け渡しが自然
- `comm` のような caller ごとに異なる値は subroutine 引数の optional として直交化

### 例外: m_struc_def のポインタ配列型

Fortran には「allocatable 配列の配列」(ragged array) を直接書く構文がない。
`real(8), allocatable :: a(:)` の配列 `a(n)` は作れない。

`m_struc_def.f90` はこの言語制約を回避するための最小限の wrapper type を定義している:
```fortran
type s_rv1
   real(8),allocatable:: v(:)
end type s_rv1
type s_rv2
   real(8),allocatable:: v(:,:)
end type s_rv2
! ... s_rv4, s_rv5, s_cv3, s_cv4, s_cv5, s_nv2, s_sblock
```

使用例: `type(s_rv2) :: basis(natom)` で、`basis(1)%v` は (10,5)、`basis(2)%v` は (8,3) のように各原子で異なるサイズの配列を持てる (C の `double **a` に相当)。

これは singleton 規約の例外ではなく、**状態を持たない純粋なコンテナ型**。
振る舞い (subroutine) を持たず、module 間のデータフローを隠蔽しない。
singleton が禁じている「状態管理 type」とは明確に区別される。

### 適用ガイドライン
新規モジュールを設計する時:
- まず module 変数 (state) と公開する subroutine を決める
- `type` を導入したくなったら一度立ち止まる — 本当に複数インスタンス必要か?
- caller に `type(...) :: x` を持たせるくらいなら module-level singleton にする
- subroutine 引数は scalar/array of intrinsic types に限定 (state は use で渡す)
- comm のような「caller ごとに異なる値」は subroutine 引数の optional として直交化
- ポインタ配列 (ragged array) が必要な場合のみ m_struc_def のコンテナ型を使用
- **k 点ごとのデータは引数で明示的に渡す** (2026-03-30, addrbl): 「今の k 点」を module 状態に置く形 (`m_Igv2x_setiq`) は、
  GPU で複数の k 点を同時に処理すると使えない。`napw, ndimh, igapw` のような k ごとのデータは `m_Igv2x_getiq` で取り、
  引数で渡す。module 状態に置くのは全 k 共通のものだけ

### Python/PyTorch への移植性

singleton module 設計は Python module と 1:1 対応する:

| Fortran | Python |
|---------|--------|
| `module m_foo` | `m_foo.py` |
| `use m_foo, only: bar` | `from m_foo import bar` |
| `protected :: bar` | module-level 変数 (慣習的に `_` prefix や getter) |
| `subroutine init(ef)` | `def init(ef):` + `global bar` |

type ベースの設計を Python に移植すると class 設計をゼロから考え直す必要があるが、
singleton なら Fortran module → Python module にほぼ機械的に変換できる。
GPU 計算も `torch.tensor` + `.cuda()` で module 変数を device に置けるため、
OpenACC の `!$acc` ディレクティブと構造が対応する。

将来的に考えている方向:
- **Python ラッパー**: 各 Fortran singleton module を Python module でラップし、
  外部 (ユーザー、ワークフロー管理、機械学習パイプライン) には Python インターフェースを見せる。
  Fortran は module 内部の高速計算カーネルを担う
- **段階的移行**: module 単位で Fortran → Python/PyTorch に置換可能。
  singleton 間の依存が `use only` / `import` で明示されてい���ため、
  1 module ずつ差し替えても全体が壊れない
- **既にライブラリ化済み**: ecalj は `libecaljF.so` (共有ライブラリ) + 薄い main
  プログラムの構成。各 main (lmf, hsfp0_sc, hgw_combined 等) は数十行で、
  singleton module の subroutine を呼ぶだけ。Python から `ctypes` / `f2py` で
  `libecaljF.so` を直接呼べるため、main を Python で書き直すのは容易
- **gwsc のライブラリ直接呼び出し化**: 現在の gwsc は subprocess で Fortran バイナリ
  を起動しているが、将来は `libecaljF.so` の subroutine を Python から直接呼ぶ構成に
  移行可能。プロセス起動オーバーヘッド (MPI_Init, 入力ファイル再読み込み、共有ライブラリ
  再ロード) が消え、singleton module 変数がプロセス内で共有されるため、ステップ間の
  状態受け渡しがゼロコストになる。hgw_combined が hrcxq + hsfp0_sc を 1 プロセスに
  統合したのと同じ原理を、gwsc 全体に拡張するもの

## ビルド環境

手順は [ForDevelopers.md](ForDevelopers.md) §4、計算機ごとの中身と落とし穴は [handover.md](handover.md) §2〜§4（2026-10-02 21:06 に重複を外した）。

## GPU 開発の教訓

### WB.4 async 事件 (2026-05-01 → 05-04)

`m_sxcf_sc.f90` に `!$acc async(1)` を追加して GPU kernel を非同期実行する最適化を入れたところ、hsfp0_sc_mp_gpu が生成する sigm が完全に壊れ、全物質の lmf が `zhev_tk4: nev/=nevout` で即死した。

**教訓:**
1. **testecalj PASS は GPU 信頼性の証明にならない** — 2 原子系ではタイミング依存バグが顕在化しない。実物質 (6-8 原子系) でのテスト必須
2. **GPU async は同期ポイントの設計が極めて難しい** — `!$acc wait` が 1 箇所でも不足すると sigm が完全破壊 (中間的劣化ではなく NaN/overflow)
3. **ベンチマーク未実施で production に入ってはいけない** — 性能ゲイン未計測のままリスクを取った
4. **async 再導入時の必須テスト**: GW1500 実物質で sigm の bit-wise 一致テスト、`!$acc wait` の網羅性検証 (全 data region exit 前、全 host_data 前、全 CPU-GPU 境界)

### OpenACC 一般注意
- nvfortran は `!$acc data` (structured) スコープ内の Fortran `BLOCK` を許可しない → `!$acc enter/exit data` (unstructured) で回避
- `attributes(device)` 変数は Fortran スコープで自動解放される → スコープを跨ぐ場合は外部スコープに移動
- `!$acc wait(1)` は goto パス (エラーハンドリング) でも必要
- `!$acc routine seq` の中の自動配列は、呼ぶたびにデバイスのヒープ確保になる (2026-09-27: m_bessl の Bessel 表が 1 回 101 ms → 固定長にして 15 ms)。
  作業配列は固定長 (`lmaxb` のような上限) にし、上限を超えない条件をコメントに書く
- OpenACC の非同期キュー 1 は non-blocking のストリームで、既定ストリームの cuBLAS の積を待たない (2026-09-27, m_sxcf_sc で NaN)。
  非同期にするなら m_blas の積も `cublas_set_stream` で同じストリームに揃え、ストリームを替える所で `cudaDeviceSynchronize`
- `m_stopwatch` の start/pause はデバイスを同期させる。非同期の区間の中のタイマーは `stopwatch_init(..., hostonly=.true.)` にする (2026-09-27)
- デバイス配列の確保・解放もデバイスを待たせる。非同期の区間の外で行う

## 入力ファイル: TOML 方式 (2026-05 〜)

### 新方式
Fortran バイナリは以下のみを読む:
- `ctrlg.<sname>.toml` — ctrl + GW driver sections ([gw] [mlo] [blocks]) + 末尾の
  [product_basis] (cut-offs と per-atom tables nlx / valence / core) の統合 TOML
- 旧 `PB.<sname>.toml` / `esm_input.dat` / `GWinput.toml` は読まない。前二者が残って
  いれば abort → `ctrlg_absorb.py <sname>` で ctrlg に取り込む (2026-09-17)

### 生成
```bash
# POSCAR から新規生成
ctrlgenToml.py <sname>        # writes ctrlg.<sname>.toml

# 旧形式 (ctrl + GWinput) からの移行
Legacy2toml.py <sname>        # writes ctrlg.<sname>.toml
```

### 旧 -v オーバーライドの新形式
```
OLD:  lmf si -vnk=8 -vmetal=3
NEW:  lmf si --ctrlg:bz.nkabc=[8,8,8] --ctrlg:bz.metal=3
```

### ファイル命名規約
- `ctrlg.foobar.toml` (foobar = 物質ID, mp-1234 等)
- ディレクトリで識別可能でも**ファイル名に ID を残す** — HTP バッチで取り違え防止、LLM のクロスリファレンス誤認防止

## hgw_combined (in-memory W)

（今のプログラム名は `hgw`（`hgw_mp_gpu` など、エントリは main_hgw.f90）。hgw_combined・hrcxq はこの統合をした 2026-05 の時期の名前）

hrcxq (screened interaction W) + hsfp0_sc (exchange Sx + correlation Sc) を 1 MPI プロセスに統合。
`__WVR/__WVI` ファイル中継 (30-90 GB/物質) を排除し、module-level buffer で in-memory 受け渡し。

gwsc での呼び出し:
```python
# 旧 (3ステップ、ファイル中継)
hsfp0_sc --job=1   # Sx
hrcxq --jobgw=1    # W  →  __WVR/__WVI 書き出し
hsfp0_sc --job=2   # Sc ←  __WVR/__WVI 読み込み

# 新 (1ステップ、in-memory)
hgw_combined --jobgw=1   # Sx + W + Sc (メモリ内)
```

### 開発の経緯 (WB.1 → WB.3e)

段階的に singleton リファクタ + ストレージ抽象化 + streaming 統合を進めた:

| Step | 内容 | 主要ファイル |
|------|------|-------------|
| WB.1 | m_wv_storage: FILE/MEMORY_3D 両 backend のストレージ抽象化 | m_wv_storage.f90 |
| WB.2a | m_sxcf_sc: sxcf_scz_correlation を init/step_kx/finalize に分割 | m_sxcf_sc.f90 |
| WB.2b | Phase 1-B 4D in-memory モード削除 (WB.1 の MEMORY_3D に統一) | m_llw, m_w0w0i, main_hrcxq |
| WB.3a | m_wv_storage: type(wv_storage) → module-level singleton 化 | m_wv_storage.f90 |
| WB.3b | m_sxcf_sc: type(sxcf_state) → module-level singleton 化 | m_sxcf_sc.f90 |
| WB.3c | m_hsfp0_sc: hsfp0_sc を setup/consume/writeout 三相に分割 | main_hsfp0.sc.f90 |
| WB.3d | m_hgw_iq_loop: iq production loop を main_hrcxq から抽出 | m_hgw_iq_loop.f90 |
| WB.3e | streaming 統合: per-iq W 生成と per-kx sxcf 消費をインターリーブ | main_hrcxq.f90 |

### アーキテクチャ

```
main_hrcxq.f90 (= hgw_combined のエントリポイント)
  │
  ├── m_hgw_iq_loop: iq ループ (W の生成、各 iq で)
  │     ├── m_llw: W(iq) を計算 → m_wv_storage に書く
  │     └── m_w0w0i: W(iq=1) の修正 (effective W(0))
  │
  ├── m_sxcf_sc: kx ループ (self-energy の消費、各 kx で)
  │     ├── sxcf_correlation_init()
  │     ├── sxcf_correlation_step_kx(kx)  ← m_wv_storage から W を読む
  │     └── sxcf_correlation_finalize()
  │
  └── m_hsfp0_sc: writeout (SEC/SEX ファイル書き出し)

m_wv_storage: W データの橋渡し
  ├── FILE backend: __WVR/__WVI ファイル (旧方式、fallback)
  └── MEMORY_3D backend: module-level 3D buffer (hgw_combined 用)
```

### singleton が可能にしたこと

この統合は singleton 設計なしには困難だった:
- **m_wv_storage が singleton** だから、hrcxq (producer) と sxcf (consumer) が同じ module 変数を共有できる。type で渡すなら巨大な W バッファの所有権管理が必要
- **m_sxcf_sc が singleton** だから、init/step_kx/finalize の3相に分割しても状態が module 変数に保持される。type なら caller が state を管理しなければならない
- **m_hsfp0_sc が singleton** だから、hrcxq の後に skip_init で呼べる。独立バイナリだった hsfp0_sc を「関数呼び出し」に変換できた
- **backend 切替が caller 透過**: m_wv_storage 内部で FILE → MEMORY_3D に切り替えても、m_sxcf_sc と main_hrcxq のコードは変更不要

### WB.4 (async) の失敗と教訓

WB.3e の上に GPU async overlap (stream 1 で imagaxis を非同期実行) を追加した WB.4 は、
testecalj では PASS したが実物質で sigm を破壊した。詳細は「GPU 開発の教訓」セクション参照。
現在は WB.3e ベース (同期実行) に戻して production 稼働中。
（2026-09-28 の注: 2026-09-27 に Σc のバッチの非同期化を入れ直した。原因だったのは、OpenACC のキュー 1 が non-blocking のストリームで
既定ストリームの cuBLAS の積を待たないこと。m_blas の積を同じストリームに揃え、LiTi₂O₄ 6³ の Σ が前の版と表示の桁まで一致し、
2 回回しても一致することを確かめた。上の「OpenACC 一般注意」と m_sxcf_sc の `sigma_stream_begin/end`）

## GW1500 量産インフラ

量産の仕組み（スロットスケジューラ、NaN の監視）は [`ecalj_auto/README_slot_scheduler.md`](../ecalj_auto/README_slot_scheduler.md)、結果と失敗の分類は [`ecalj_auto/GW1500_status.md`](../ecalj_auto/GW1500_status.md)
（2026-10-02 21:06 にこの節の写しを外した。5 月の NaN は「物質に固有」と書いていたが、旧 `--mp` の TF32 の精度が原因だった。GW1500_status §3・§5.2）。
