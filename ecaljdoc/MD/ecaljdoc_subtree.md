# ecaljdoc を ecalj の中で管理する手順（git subtree、2026-10-02）

2026-10-02 に、文書のリポジトリ ecaljdoc を ecalj の中のディレクトリ [`ecalj/ecaljdoc/`](../README.md) に取り込んだ（user「一体管理したほうが便利」「いまの公開サイトは維持」）。
文書の正本は [`ecalj/ecaljdoc/`](../README.md)。公開サイト（GitHub Pages、https://ecalj.github.io/ecaljdoc/）は今までどおり GitHub の **[`ecalj/ecaljdoc`](../README.md)**（organization の ecalj）の
`main` から出る（2026-10-02 21:05 に `gh api repos/ecalj/ecaljdoc/pages` で確かめた。`tkotani/ecaljdoc` は開発用の写しで Pages は無い）。そこへは ecalj から
`git subtree push` で文書の部分だけを送る。URL と文書の中のリンクは変わらない。
（2026-10-02 21:05 訂正: 最初の版はこの節と §3・§4 で公開を `tkotani/ecaljdoc` としていた。手元の古いクローンの remote がそれだけだったための思い込み）

## 1. 何をしたか

```bash
cd ~/ecalj
git subtree add --prefix=ecaljdoc /home/takao/ecaljdoc main --squash
```

- ecalj に二つのコミットが入った: ecaljdoc の中身（ecaljdoc のコミット `e6bb431` の時点）を一つにまとめたコミットと、それを `ecaljdoc/` に取り込むマージのコミット
  （ecalj の `b02f59e99`）。ecaljdoc の古い履歴（353 コミット）は ecalj に入れていない。GitHub の `tkotani/ecaljdoc` と、手元の `~/ecaljdoc`（使わない控え）にある
- 取り込んだ大きさは約 40 MB（発表の PDF [`ecaljdoc/presentations/`](../presentations) 26 MB、図 [`ecaljdoc/public/`](../public) 9 MB）
- 確かめ: `git subtree split --prefix=ecaljdoc` で切り出したツリーが ecaljdoc の `e6bb431` と同じ（d0f9cdb46b3a）。その後 ecalj の側で直した文書を
  切り出すと、`e6bb431` を親に持つコミットになる。[`ecalj/ecaljdoc`](../README.md) の main（`2012b30`）も `tkotani/ecaljdoc` の main（`2482c4b`）も `e6bb431` の祖先なので（2026-10-02 21:05 に確かめた）、どちらへも早送りで送れる

## 1a. 開発の記録 `MD/` も ecaljdoc の中に（2026-10-02 の後半）

`ecalj/MD/` を [`ecalj/ecaljdoc/MD/`](.) に移した（user「MD ごと ecaljdoc に、vitepress は無視してくれたらいい」）。[`ecaljdoc/.vitepress/config.mts`](../.vitepress/config.mts) の
`srcExclude: ['MD/**', 'trash/**']` でサイトには出さない（MD の本文は `<sname>` などを Vue がタグと読んでビルドが止まる。2026-10-02 に確かめた）。
`git subtree push` で [`ecalj/ecaljdoc`](../README.md) にも送られる（ecalj の公開リポジトリにも入っているので同じ）。

## 2. ふだんの作業

- 文書は [`ecalj/ecaljdoc/`](../README.md) の下を直す。コードと文書を同じコミットで直してよい。コミットメッセージは英語（ecalj の決まり）
- **`~/ecaljdoc`（古いクローン）では作業しない**。直接 [`ecalj/ecaljdoc`](../README.md)・`tkotani/ecaljdoc` を直すこともしない（した場合は §4）
- 手元で見る（Node が要る）:

```bash
cd ~/ecalj/ecaljdoc
npm install            # 初回だけ（node_modules は ecaljdoc/.gitignore で無視される）
npm run docs:dev       # http://localhost:5173/ecaljdoc/
```

- [`ecaljdoc/.github/workflows/deploy.yml`](../.github/workflows/deploy.yml) は ecalj の中では働かない（GitHub は最上位の `.github` だけを見る）。送った先の [`ecalj/ecaljdoc`](../README.md) で働く

## 3. 公開（push。メンテナの指示があってから）

文書の部分を [`ecalj/ecaljdoc`](../README.md)（公開）へ送る。開発用の写し `tkotani/ecaljdoc` へも同じに送ってよい。

```bash
cd ~/ecalj
git remote add ecaljdoc git@github.com:ecalj/ecaljdoc.git        # 初回だけ（公開）
git remote add ecaljdoc-dev git@github.com:tkotani/ecaljdoc.git  # 初回だけ（開発用の写し）
git subtree push --prefix=ecaljdoc ecaljdoc main          # 公開
git subtree push --prefix=ecaljdoc ecaljdoc-dev main      # 開発用の写し（任意）
```

- `subtree push` は ecalj の履歴から `ecaljdoc/` に触れたコミットだけを切り出し（`git subtree split`）、送り先の main に送る。
  送った先では、ファイルは今までどおり最上位にある。push で GitHub Actions がサイトを作って公開する
- コードと文書を両方直したコミットは、文書の部分だけが送られる（メッセージはそのまま）
- 送る前に確かめるとき:

```bash
S=$(git subtree split --prefix=ecaljdoc)       # 送られるコミット
git fetch ecaljdoc main
git merge-base --is-ancestor ecaljdoc/main $S && echo fast-forward   # 早送りで送れる
git log --oneline ecaljdoc/main..$S                                   # 送られるコミットの一覧
```

- 公開のタイミングは ecalj の push（dev・rel）と独立に決めてよい（user 2026-10-02）。`subtree push` はその時点の ecalj の HEAD までの文書の変更を
  全部送る（まだ push していないコードの機能を説明した文書も出る）。途中の時点までを出すときは、そのコミットで切り出して送る:
  `git push ecaljdoc $(git subtree split --prefix=ecaljdoc <ecalj のコミット>):main`
- バックアップは ecalj だけでよい（文書と開発の記録も入っている）
- `subtree split` は ecalj の全履歴をたどるので、履歴が長くなると遅くなる。遅いときは `git subtree split --prefix=ecaljdoc --rejoin` で
  切り出した点を ecalj に記録しておく（次からそこより後だけをたどる）

## 4. [`ecalj/ecaljdoc`](../README.md)（や `tkotani/ecaljdoc`）の側で直接直された場合

GitHub の画面などで送り先が直接直されたら、ecalj に取り込んでから作業を続ける（取り込まずに `subtree push` すると早送りにならず拒まれる）。

```bash
cd ~/ecalj
git subtree pull --prefix=ecaljdoc ecaljdoc main --squash
```

## 5. ほかの機械

- `TOOLS/sync_ecalj_src.sh --samples <host>` は追跡しているツリーを全部送るので、`ecaljdoc/` も送られる（約 40 MB 増える）。`--samples` なしは SRC と InstallAll.py だけ
- ほかの機械では文書を直さない（直すのは t14 の ecalj だけ）

## 6. 文書の中の書き方

- ecalj の文書（[`ecaljdoc/MD/`](.)、`Changes.txt`、各 README）から ecaljdoc を指すときは、[`ecaljdoc/manual/mlo.md`](../manual/mlo.md) のようにリポジトリの中のパスで書けばよい。
  公開サイトの URL（`https://ecalj.github.io/ecaljdoc/manual/mlo`）も今までどおり使える
- 2026-10-02 より前の記録にある「ecaljdoc のコミット `xxxxxxx`」は、[`ecalj/ecaljdoc`](../README.md)・`tkotani/ecaljdoc`（と `~/ecaljdoc`）の履歴のハッシュ。ecalj の中には無い
  （ecalj の中の ecaljdoc の最初は、`e6bb431` をまとめたコミット）
