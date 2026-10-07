# ecaljdoc を ecalj の中で書き、公開リポジトリへ写して出す手順（2026-10-02）

文書（ecaljdoc）の正本は ecalj のディレクトリ [`ecaljdoc/`](../ecaljdoc/README.md)。公開サイト（GitHub Pages、https://ecalj.github.io/ecaljdoc/）は
GitHub の **`ecalj/ecaljdoc`**（organization の ecalj）の `main` から出る（2026-10-02 21:05 に `gh api repos/ecalj/ecaljdoc/pages` で確かめた。
`tkotani/ecaljdoc` は以前の開発用の写しで、Pages は無い）。公開するときは、ecalj の `ecaljdoc/` を公開リポジトリの clone へ**写して**コミットし、push する
（[`TOOLS/publish_ecaljdoc.sh`](../TOOLS/publish_ecaljdoc.sh)）。URL と文書の中のリンクは変わらない。開発の記録 [`MD/`](ecaljclaude.md) は ecaljdoc の外にあり、公開しない。

## 1. 経緯（2026-10-02 の一日のうちに）

1. ecaljdoc を `git subtree add --squash --prefix=ecaljdoc`（`~/ecaljdoc` の `e6bb431` の時点）で ecalj に取り込んだ（ecalj `b02f59e99`、user「一体管理したほうが便利」）。
   ecaljdoc の古い履歴（353 コミット）は ecalj に入れていない。GitHub の `ecalj/ecaljdoc` と、手元の古いクローン `~/ecaljdoc`（使わない。目印 `MOVED_TO_ecalj_ecaljdoc.txt`）にある
2. 開発の記録 `MD/` を一度 `ecaljdoc/MD/` に移し、VitePress の `srcExclude` でサイトから外した（user「MD ごと移して vitepress が無視すればいい」）
3. 公開を `git subtree push` にする手順を書いたが、user「subtree にしなくても、写して出せば公開のときに困らない」「公開リポジトリに MD は送らない」
   →「それなら ecaljdoc の下にする必要はない」で、`MD/` を最上位に戻し、公開は写して出す形にした（この文書。もとは `ecaljdoc_publish.md`）
4. 途中で、公開のリポジトリを `tkotani/ecaljdoc` と書いていたのを `ecalj/ecaljdoc` に訂正した（手元の古いクローンの remote がそれだけだったための思い込み）

`ecaljdoc/` は ecalj の普通のディレクトリ。subtree の印は 1 の取り込みのコミットのメッセージ（`git-subtree-dir`、`git-subtree-split`）に残るだけで、使わない。

## 2. ふだんの作業

- 文書は [`ecaljdoc/`](../ecaljdoc/README.md) の下を直す。コードと文書を同じコミットで直してよい。コミットメッセージは英語（ecalj の決まり）
- **`~/ecaljdoc`（古いクローン）では作業しない**。公開リポジトリ（`ecalj/ecaljdoc`）を GitHub の画面などで直接直すこともしない（次の公開で上書きされる）
- 手元で見る（Node が要る）:

```bash
cd ~/ecalj/ecaljdoc
npm install            # 初回だけ（node_modules は ecaljdoc/.gitignore で無視される）
npm run docs:dev       # http://localhost:5173/ecaljdoc/
npx vitepress build    # ビルドが通るか（本文の <...> や {{ をコードの印の外に書くと止まる）
```

- `ecaljdoc/.github/workflows/deploy.yml` は ecalj の中では働かない（GitHub は最上位の `.github` だけを見る）。写した先の `ecalj/ecaljdoc` で働く
- サイトのページから ecaljdoc の外（Samples、ソース、MD）へは GitHub の URL で書く（`../` の相対リンクはサイトで切れる）。`MD/` へはリンクを張らない。
  確かめは `python3 TOOLS/doclinks.py`

## 3. 公開（push はメンテナの指示があってから）

```bash
cd ~/ecalj
TOOLS/publish_ecaljdoc.sh          # HEAD の ecaljdoc/ を ~/work/ecaljdoc_publish（ecalj/ecaljdoc の clone）へ写してコミット。差分を見せる
TOOLS/publish_ecaljdoc.sh --push   # そのうえで push（GitHub Actions がサイトを作り直す）
```

- 写すのは HEAD のコミット済みの `ecaljdoc/`（`git archive`）。`ecaljdoc/` に未コミットの変更があると止まる
- 公開側の clone は写す前に `origin/main` にそろえる（そこでの直接の編集は上書きされる）
- 公開側の履歴は「`Sync from ecalj <ハッシュ>`」のまとまったコミットになる。文書の細かい履歴は ecalj の側にある
- 公開のタイミングは ecalj の push（dev・rel）と独立に決めてよい（user 2026-10-02）。その時点の HEAD の文書が全部出る
  （まだ push していないコードの機能を説明した文書も出る）。途中の版を出すときは、そのコミットを checkout した ecalj で回す
- 最初の公開（2026-10-02 以後の最初）では、サイトから外したページ（ForDevelopers など、今は `MD/`）と trash に移したもの（`BackUp/`、`ecaljdetails/` など）が
  公開側から消える

### 3.1 三つを順に出す（2026-10-03、user）

公開するリポジトリは三つで、順番がある（サイトとデータベースは GitHub の ecalj のファイルを指す）: ① ecalj の push（dev、試験のあと rel）、
② ecaljdoc の公開（この文書）、③ GW1500 のデータベース（`ecalj_auto/gw1500db_snapshot.py` で日付つきのスナップショット `QSGW80_<yyyymmdd>` を作り、`TOOLS/publish_gw1500db.sh` で
github.com/tkotani/DOSnpSupplement へ。木には最新の版だけを置き、前の版はタグ `QSGW80_<yyyymmdd>` と差分の過去ログ `changes_from_<yyyymmdd>.md` で残す（user 2026-10-07）。
2025 年の補遺は `DOSnp2025.md` と元のファイル）。一括は [`TOOLS/publish_all.sh`](../TOOLS/publish_all.sh)（`--push` で push、`--rel` で rel も）。説明は最上位の [`README.md`](../README.md)

## 4. ほかの機械

- [`TOOLS/sync_ecalj_src.sh`](../TOOLS/sync_ecalj_src.sh) `--samples <host>` は追跡しているツリーを全部送るので、`ecaljdoc/` と `MD/` も送られる。`--samples` なしは SRC と InstallAll.py だけ
- ほかの機械では文書を直さない（直すのは t14 の ecalj だけ）

## 5. 文書の中の書き方

- ecalj の文書（`MD/`、`Changes.md`、各 README）から ecaljdoc を指すときは、`ecaljdoc/manual/mlo.md` のようにリポジトリの中の相対パスで書く。
  公開サイトの URL（`https://ecalj.github.io/ecaljdoc/manual/mlo`）も使える
- 2026-10-02 より前の記録にある「ecaljdoc のコミット `xxxxxxx`」は、`ecalj/ecaljdoc`・`tkotani/ecaljdoc`（と `~/ecaljdoc`）の履歴のハッシュ。ecalj の中には無い
  （ecalj の中の ecaljdoc の最初は、`e6bb431` をまとめたコミット）
