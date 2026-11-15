# ecalj_knowledge.md — コードとサンプルで確かめた事実、文書を直す手順（AI 向け）

> **人間向けの文書ではない。** AI（Claude）が作業の前に読み、事実を引くためのメモ。読みやすさより、引ける正確さ（ファイル、ルーチン、日付、コミット）を優先している。
> 人が読む説明は ecaljdoc（公開サイト）にある。

2026-11-15 00:04 に書いた。2026-11-14 のセッション（ecaljdoc の整理、日本語化、二重化の片付け、朱の確定）で、
コード・サンプル・git で確かめたことをまとめる。Claude が後で引くためのもの。人への説明はここから組み立てる。
入口は [`CLAUDE.md`](../CLAUDE.md)（Claude 以外の AI は [`AGENTS.md`](../AGENTS.md)）。取り決めと計算機は [handover.md](handover.md)、まだのことは [TODOandQuestion.md](TODOandQuestion.md) の「ecaljdoc（文書）」、
ecaljdoc の書き方の決まり（公開）は [`ecaljdoc/manual/writing_docs.md`](../ecaljdoc/manual/writing_docs.md)。ここではそれらを繰り返さない。

**使い方**: 書いてあることも、使う前にコードで確かめる（日付の時点の事実）。違っていたら直し、日付を付ける。

---

## 1. ecaljdoc のどこに何があるか（正本の割り当て、2026-11-14）

同じ説明を二か所に書かない。ほかのページはリンクだけにする。

| 話題 | 正本のページ |
|---|---|
| 入力ファイル、`Legacy2toml.py`、古いキーと TOML のパスの対応（表 1 `#legacy-toml-map`）、生成されるファイルの体裁（`#file-layout`）、`ctrlg_update.py` | `manual/toml_migration.md`（メニュー「入力ファイル」） |
| `ctrlg.<sname>.toml` の構成、lmf のキー、LDA+U（`#lda-u-toml`）、芯の穴、ローカル軌道、`symgrp`/`symgrpaf`（`#symgrp`）、`save.<sname>`（`#save-sname`）、出力ファイル、`readp`/`pnufix`（`#readp-option`） | `manual/lmf.md` |
| コマンドラインのオプションの一覧、`--ctrlg:` の仕様（`#run-time-ctrlg-overrides`）、`--fullstdo` | `manual/cmdopts.md` |
| `[gw]` のキー、積基底 `[product_basis]`（`#product-basis`、2026-11-14 に gwsc から移した） | `manual/gwinput.md` |
| `gwsc` のオプション、反復の混合 `mixbeta`（`#mixbeta-and-mixpriorit`、lmf から移した）、gwsc の使うファイル、QPU の列 | `manual/gwsc.md` |
| GPU の精度（`--prec`、`--mp`、`--linalg`、GEMMul8）、GPU のロック | `manual/ecaljgpu.md` |
| ESM | `manual/esm.md` |
| `[dyn]`（緩和） | `manual/relax.md` |
| サンプルの一覧 | `manual/samples.md`（ecalj の `Samples/README.md` の表 1 が元） |
| スピン分極・反強磁性・SOC のバンド・PROCAR・4f の手順 | `manual/UsageDetailed.md`（メニュー「その他の説明」） |
| 誘電関数 | `manual/optical.md`（2026-11-14 に `Samples/EPS/EPS_GaAs/README_eps.md` を土台に書き直した） |
| 試験の仕組み | `manual/dev_testecalj.md` |
| 有限温度 | `manual/kBT.md` |

メインドキュメント（`manual/README_tutorial.md`）は通しのチュートリアルなので、画面の出力の例などを持ってよい（gwsc はそこへリンクする）。

## 2. 文書を直す手順（Claude が実際に使ったもの）

- **英語版の印**: `ecaljdoc/en/X.md` の 3 行目が「…made from X.md of commit <c> (source blob <b>)」。日本語の側をコミットした後、
  `C=$(git rev-parse --short=9 HEAD); B=$(git rev-parse HEAD:ecaljdoc/X.md | cut -c1-12)` を 3 行目に書く（sed で
  `of the tree of …` と `of commit …` の両方の形を置き換える）。`python3 TOOLS/ecaljdoc_en.py status` が空ならそろっている
- **メニュー**: `.vitepress/config.mts` の `sidebarRoot`。日本語の見出しには `enText` に英語を足す。config を変えると開発サーバが再起動するので、確かめる前に `curl --retry` で待つ
- **見出しのアンカー**: 日本語の見出しにすると自動のアンカーが変わる。ほかのページからリンクされている見出しには `{#…}` を明示する
  （訳す前に `grep -rno "page#[A-Za-z0-9_-]*"` でリンクされているアンカーを集める）
- **数式のはみ出し**: 番号付きの式は CSS（`custom.css`）で `overflow-x: visible` なので、スクロールバーは出ずに右の欄へはみ出す。
  検査は playwright で、各 `mjx-container[display="true"]` の `svg` の右端と `.vp-doc` の右端を比べる（1300px と 960px）。
  `scrollWidth > clientWidth` だけでは見落とす。直し方は `aligned` で折り返す。`\braket{…}` は行をまたげないので、`\langle … \rangle` に分ける
- **doclinks の誤検出**: 数式の中の `[n](r)` のような並びを Markdown のリンクと見る（2026-11-14 に `V_\mathrm{eff}[n](r)` で起きた）。書き方を変えて避ける
- **朱の確定**: `del.clfix` を消し、`span.clfix` と `div.clfix` の中身を本文にし、`div.clfix-old` を消す。確定のとき「…から移した」「もとの式は…」のような経緯だけの注も消す。
  確かめられないものは残す（user 2026-11-14「わからないものは残す」。残したものは TODO にページごとに書いた）。誤記は印を付けずに直してよい（同）
- **開発サーバ**: `NODE_OPTIONS=--max-old-space-size=16384 npm run docs:dev -- --host 127.0.0.1 --port 5173`。編集を重ねるとメモリが増えて落ちる（終了 134）ので、立て直す。
  ビルドの確かめは `npx vitepress build .`（約 90 秒）
- **サブエージェントの報告**: そのまま採らず、両方のページとコードで確かめてから直す（2026-11-14 の重複の報告は 3 件の食い違いが本当、正本の見立ては 2 件が要検討だった）

## 3. コードで確かめた事実（2026-11-14）

### 入力と既定値
- `[gw]` の必須のキー: `t_sigmaw`、`HistBin_dw`、`HistBin_ratio`、`QforEPSau`、`EMAXforGW`（2026-11-14 から。`m_GWinput.f90` の `need_gw_key`）、`t_tetrakbt`（2026-09-28 から）。
  無い入力は止まる。`ctrlg_update.py`（`pylib/ctrlg_rules.py` の規則 3）が以前の既定値を書き足す
- gwinit の雛形（`main_gwinit.f90`）: `t_tetrakbt = 300`（2026-06-13 から書き、09-17 にコメントにし、09-28 の b4b64be66 でふたたび書く）、
  `HistBin_dw = 1e-4`、`HistBin_ratio = 1.05`、`QforEPSau = true`
- 古い `esmr` は読まない。`ctrlg_update.py` が `t_sigmaw` に換算して書く。換算は π/√3 を 1.81 に丸める（`ctrlg_rules.py`。0.003 Ry → 262 K）
- 実軸のメッシュ: $\omega_i = b(e^{a(i-1)}-1)$、$a$ = `HistBin_ratio` − 1、$b$ = `HistBin_dw`/$a$（`m_freq.f90` の `Getfreq`）
- `[dyn]` の既定（`m_lmfinit.f90`）: `mode` 0（4 CG、5 Fletcher-Powell、6 Broyden）、`nit` 1、`hess` true、`xtol` 1e-3、`gtol` 0、`step` 0.015、`nkill` 0。
  緩和のモードでは `AtomPos.<sname>` があると位置を上書きする（`main_lmf.f90` の `ReadAtomPos`）
- ドーピングのキー: `[[spec]].q`、`[bz].zbak`（`ctrl_schema.py`）
- `[[spec]]` の `idu` は 0/1(AMF)/2(FLL)、+10 で `sigm` があれば U を使わない

### 動径関数と基底
- `[ham].readp = true`、`pnufix = true` が `ctrlgenToml.py` の既定。`lmfa` が `atmpnu.*` に書いた MT 半径での対数微分（価電子の $P$ とローカル軌道の $P_z$）を
  `lmf` が読み（`m_lmfinit`）、反復の間固定する（`m_bndfp` で `pnunew` を呼ばない）。動径関数は LDA のポテンシャルの球対称部分で解く
- 占有状態の重心（部分状態密度の中心）で解くのは `pnufix = false` のときの `pnunew` だけ（ローカル軌道の無い $l$ の価電子）。以前の方法で、不安定だった
- GW の波動関数の表現: $\phi_{Ru}$ は MT の中だけの関数（MTO そのものではない）、IPW は MT の中で 0。そのため内積の交差項は消える。IPW の重なりは $\Omega\delta_{GG'}$ − MT の積分
- 「2 種の MPB」は M（もとの混合積基底 $M_I$、直交しない）と E（Coulomb を対角にする直交基底 $E_\mu=\sum_I z_{\mu I}M_I$、係数は `zcousq`）のこと（user 2026-11-15）。
  `hgw` は χ0 の足し込みと Σc を M で（`build_zmel(..., is_m_basis=.true.)`）、ε の逆行列と W^c を E で（`m_llw` で ngb×ngb を変換）、Σx を E で行う。W^c は M でしまう。
  `hx0fp0` は E だけ。論文 Obata 2025 は E だけで書いている。詳しくは [`ecaljdoc/manual/dev_impl.md`](../ecaljdoc/manual/dev_impl.md) の (I14a)
- コアの交換項 Σ^Cx 用の MPB_C（M_C と E_C）は別にあり、使っている: `gwsc` で `hbasfp0 --job=3`（M_C、コア×価電子の動径関数、`incwfin=-2`）→ `hvccfp0 --job=3`（Coulomb 行列と E_C、`__Vcoud.*` はあとの `--job=0` で上書き）→ `hsfp0_sc --job=3`（E_C で Σ^Cx）。dev_impl の (I15) の下

### 出力とファイル
- `save.<sname>` の行の先頭の文字は `svflg='cxhiC'`（`m_lmfp.f90`）: c 収束、x 上限に達した、h 最初（Harris）、i 途中、C 緩和の収束
- 軌道磁気モーメントは画面（`IORBTM:`）とファイル `orbitalmom.chk`（`m_bndfp_util.f90`）、スピンは `mmom.<sname>.chk`（`mkrout.f90`）
- 反強磁性のモーメントを保つモード（`m_ldau.f90`）: 入力 `mmtarget.aftest`、出力 `mmagfield.aftest`、`mixmag.aftest`、画面の `mmaftest:`・`mmhist:`
- `pwmode = 1` で `symgrpaf` を使うと `rotwave: q+G rotation error`（`rotwave.f90`）。`pwmode = 11` が要る
- ESM: `readesm`（`m_lmfinit.f90`）は `[esm]` を読み、`esm_input.dat` が残っていれば止まる

### job_eps（`SRC/exec/job_eps`）
- 段: `lmf --jobgw=0 --novxc`、`qg4gw --job=2 --n1n2n3eps`、`lmf --jobgw=1`、`heftet`、`hbasfp0`（`--lcf` で `--job=0`、なしで `--job=4`）、`hvccfp0 --job=202`、`hx0fp0 --job=202|203`
- k 点は `[gw] n1n2n3eps` があればそれ、無ければ `n1n2n3`（`mkqg.f90`）
- 出力 `EPS000n.nlfc.dat`（`--lcf` で `EPS000n.dat` も、`--decompose` で `.interbandonly`/`.intrabandonly`）。列: q（2π/alat）、ω（Ry）、ε、1/ε（`main_hx0fp0.f90`）
- gnuplot: `eps_without_lcf_<sname>.glt`、`eps_with_lcf_<sname>.glt`、`eps_{interbandonly,intrabandonly,total}_<sname>.glt`。参照のファイルは指数が `D` のことがある（`set datafile fortran`）
- `[blocks]` の `QforEPS` は `[gw]` に無いときだけ読む

### 廃止したオプション（与えると止まる。`m_cmdopt_registry.f90` の `check_retired`）
- `--nb=`、`--nww=`、`--cutuu=`、`--Wtype=`、`--ahc`、`--mloahc`、`--AHCMAT`、`--UUMAT`、`--cmlo`、`--q2q1test`、`--mlo_diagnorm` など（2026-10-02）
- `--zmel0` は登録だけ残っていて何もしない

### サンプル（ecalj の `Samples/`、2026-11-14）
- `AFfixMMOM/`（`NiO_afsym`、`NiO_noafsym`、組 `affix`）がある
- `TestInstall/yh3fcc_gwsc666` は走らなかった（Y の拡張局所軌道で動径方程式が失敗）ので 2026-11-15 に外した（MD/past_log.md 表 1）
- LDA+U の試験済みの例は `LDAU/ReN`（GdN、PrN。`symgrp = "r4z"`、`so = 2`、`occnum`）。UsageDetailed の ErN の設定は以前の計算
- PROCAR: `BandWeight.py` はサンプルのディレクトリにある（`SRC/exec` には無い）

## 4. 理論のページの式番号の接頭辞（2026-11-14）

衝突を避けるため、ページごとに接頭辞を決めた。新しい式は同じ接頭辞で続ける。

| ページ | 接頭辞 |
|---|---|
| theory/basic.md | GF（Green 関数とは）、Gd・Ge・Gh（定義）、G（Gell-Mann–Low）、N（相互作用のない系）、P（分極関数） |
| theory/basic2.md | K |
| theory/gw.md | W、H（Hedin） |
| theory/qsgw.md | Q |
| theory/gw_basis.md | B、Bq（動径方程式の囲み） |
| theory/gw_product_basis.md | PB |
| theory/gw_zmel.md | Z |
| theory/optical_properties.md | O |
| theory/appendix.md | A（ブロック逆行列、ダウンフォールディング、因果律。2026-11-15 に実装メモから移した） |
| manual/dev_impl.md | I（実装メモ。番号は写したときのまま、抜けがある） |
| manual/gwsc.md | M（mixbeta） |

## 5. ドキュメントとコードの書き方（要点。正本はリンク先）

決まりの本文はリンク先にある。ここは作業の前に見落としを防ぐための要点と、2026-11-14 に user から受けた書き方の指示（正本にまだ無いもの）。
決まりを変えるときは正本を直し、ここは合わせる。

### 5.1 コード（正本: [ecaljclaude.md](ecaljclaude.md) の「コーディング規約」「記録の方針」）
- Fortran の module を singleton として使う。状態は module 変数（`protected, public`）、`use m_foo, only: …` で明示して取る。
  subroutine の引数は指示的な値（`ef`、`qp`、`isp` など）と、k 点ごとのデータに限る。`type`、`class`、type-bound procedure は使わない（例外は `m_struc_def` の入れ物の型）
- module の依存は DAG（下位だけを `use`）。大局は [module_map.md](module_map.md)（`python3 TOOLS/module_map.py` が作る。手で直さない）
- コメントは英語。何を計算するか、なぜその形か、壊さないための注意を書く。変更・判断・計測・バグ修正の注記には日付時刻（`date` の出力）を付ける。
  すぐ古くなる性能値や行番号は書かない。名前を変えたら `grep` でコメントと文書の言及も直す
- 「使われていない」と判断する前に、呼び出し元を `grep` で確かめる（2026-11-14: Numerical Recipes の `bessjy` は `smves → ropbes → sphbes → bessjy` で使われていた）
- 消したコードは消さずに `trash/` へ移し、[past_log.md](past_log.md) の表 1 と ecaljdoc の「外したコード」（`manual/past_code.md`）に書く
- コミットは英語で、ファイルを名指しする（`git commit -a` は使わない）。末尾にセッションの帰属の 2 行。push は user の指示があってから。`ecaljdoc/package-lock.json` はコミットしない
- 試験の間はインストール先の bin が指すものを変えない。別の木でビルドし、試験用の bin を別に作る。bin を丸ごと写すときは、ecalj のファイルだけを写す（bin にはほかの個人のファイルがあることがある）
- ライセンス（AGPLv3）の決まりは [`ecaljdoc/manual/license.md`](../ecaljdoc/manual/license.md)。ecalj から派生したコード（移植や、構造や手順をなぞった書き直しを含む）は AGPLv3。
  ecalj の引用は必須、spglib と（バンドの図の）seekpath の引用はお願い。ecalj のアルゴリズムの出どころを残すのはお願い

### 5.2 研究ログ・TODO・MD の文書（正本: [ecaljclaude.md](ecaljclaude.md) の「記録の方針」）
- 研究ログは最新を上に、`## 日付` の下に `### HH:MM`。時刻は推さずに `date` で入れる
- 直すべき点は見つけても、その場では直さずに TODO に書く（文書の誤りで、確かめられるものはその場で直してよい。2026-11-14 の作業はそうした）
- MD/ の文書の中のパスは相対リンクにする（`python3 TOOLS/mdlinkify.py <md>`）。確かめは `python3 TOOLS/doclinks.py`

### 5.3 ecaljdoc（正本: [`ecaljdoc/manual/writing_docs.md`](../ecaljdoc/manual/writing_docs.md)）
- 日本語のページが正本、英語版は `en/`。メニュー（`config.mts` の `sidebarRoot`）がただ一つの目次で、すべてのページを載せる。一つの話題は一つのページに書き、ほかはリンク
- 必要な理由と注意だけを書く。迷った経緯（苦労話）は書かない。数式には `\tag`、図と表にも番号
- サイトのページから `MD/` にリンクしない（言及するときは「開発記録（リポジトリの MD/…、サイトには無い）」）。ecaljdoc の外へは GitHub の URL
- 直したところを確かめてもらうときは朱（`clfix`）。確定は user が確かめてから、または確かめられたものだけ

### 5.4 2026-11-14 に user から受けた書き方の指示（ecaljdoc）
- 英語の混ざった日本語でよい。用語は無理に訳さない（例: 「準粒子自己無撞着 GW」でなく QSGW、quasiparticle self-consistent $GW$）
- 数式は横のスクロールバーを出さない。折り返してよい（全ページ同じ）
- メニューの見出し、ページの題も日本語にする（英語版のメニューは `enText`）
- Claude の助けは README とトップに一行だけ（「2026-02 から Claude Code と開発」）。誰が書いたか・誰の指示かは書き分けない
- 人が読む文書でも半分は AI が読む。釣り合いのとれた短さにし、AI のための長い断り書きは入れない
- 理論のページ（Basic I）は Fetter-Walecka に書いてあることを、見方を変えて書く: バネのインパルス応答から Green 関数へ、電子と正孔、
  T 積は摂動論で自然に出る、Gell-Mann–Low と分離したダイヤグラム、Feynman の見方（引用の囲み）、多体摂動論も通常の摂動論の書き直し
- Samples の README とサンプル（`test.py`）の方が文書より新しく正しいことが多い（user 2026-11-14）。文書を直すときは、まずそちらとコードに合わせる
- 実装メモ（`ecaljdoc/manual/dev_impl.md`）は、いまのコードの姿だけを書く。経緯（いつ何を変えたか、日付の注）はほとんど要らない（user 2026-11-15）。
  コードとの対応（どのファイル・ルーチン・キーか）は書く。古い名前を見つけたら、注を足すのでなく本文をいまの名前で書き直す
