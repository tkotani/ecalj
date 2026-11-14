# 入力ファイル {#toml-migration-2026-05}

ecalj は、一体計算の部分（DFT パート）と GW パートからなります。

今の ecalj は、TOML 形式のファイル `ctrlg.<sname>.toml`（以下 ctrlg）だけを基本の入力ファイルとして読みます。すなわち、あらゆる計算はこの ctrlg のみから出発します。
古い形式（2026-05 以前）では `ctrl.<sname>` と `GWinput` が必要でした。後方互換性のため、これらを ctrlg に変換するコンバータ
`Legacy2toml.py` は残してあります。過去の計算を再現するには、あらかじめ変換して ctrlg を作ってください（[下の節](#migration-of-an-existing-legacy-directory)）。

## 基本となる入力ファイルは TOML ファイル一つ {#the-one-toml-file}

| ファイル | 役割 |
|---|---|
| `ctrlg.<sname>.toml` | ecalj の基本入力ファイル。すべての計算はこのファイルから出発する。`[struc]`、`[[site]]`、`[[spec]]`、`[bz]`、`[iter]`、`[ham]`、`[gw]`、`[mlo]`、`[blocks]`、`[product_basis]` などの節からなる（ほかに `[options]`、`[str]`、`[ewald]`、`[esm]`、`[dyn]` など） |

基本となる入力ファイルは ctrlg 一つです。ctrlg は結晶構造から自動生成し（`ctrlgenToml.py`）、それを編集して使うのが基本です。
もしくは、生成するときにオプションを付けて生成します。
自動生成は、POSCAR（VASP の形式）を `vasp2ctrl` で構造だけのファイル `ctrls.<sname>` に変換し、それから `ctrlgenToml.py <sname>` で ctrlg を作ります。
AI に任せるなら、POSCAR が入力だと思ってかまいません。
キーの詳しい説明は、ctrl の側（`[struc]`、`[[site]]`、`[[spec]]`、`[bz]`、`[ham]` など）が [lmf](./lmf)、GW の側（`[gw]` など）が [GW の設定](./gwinput) にあります。

**最後にある原子ごとの表**

- GW 計算で使います。ふつうは手で直しません。
- `gwsc`、eps の道具、`mlo` の W の側が使う混合積基底を作るプログラム（`hbasfp0` / `hvccfp0` など）に渡すものです。
- `gwinit` が書きます（`ctrlgenToml.py` か `Legacy2toml.py` を通して）。作られたままにしておいてください。
- 2026-09 までは別のファイル `PB.<sname>.toml` でした。このファイルはもう読みません。
  残っているとバイナリは止まり、`ctrlg_absorb.py <sname>` を使うように示します。このスクリプトが表を `ctrlg` に移します。

**GW の節が要らないとき（`--skipgw`）**

- DFT だけの使い方（`lmf` / `lmfa` / バンドの図だけ）なら、`ctrlgenToml.py` に `--skipgw` を渡してください。
- `lmfa → lmf --jobgw=0 → gwinit` の段を飛ばし、出力に `[gw]` / `[mlo]` / `[blocks]` / `[product_basis]` を書きません。

**あとから GW の節を足すとき（`--addgw`、手で直したところを残す）**

- `ctrlgenToml.py <sname> --addgw` を使います。
- `[gw]` / `[mlo]` / `[blocks]` / `[product_basis]` を、ctrl の側のキーを作り直さずに書き足します。
  `[bz]`、`[ham]`、`[[spec]]` などに入れた変更はそのまま残ります。
- `[gw]` がすでにあるときは走りません（知らないうちに重複するのを防ぐため）。
- この目的で `ctrlgenToml.py <sname>` をただ走らせ直しては**いけません**。ふつうの流れでは `ctrls.<sname>` から
  `ctrlg.<sname>.toml` 全体を書き直します（前のファイルは `ctrlg.<sname>.toml.bakup` に移しますが、控えは一世代だけです）。

`<sname>` は利用者が決める物質の拡張子です（`si`、`cu`、`gaas`、`fe` など）。一つの作業ディレクトリには、ふつう `ctrlg.<sname>.toml` が
ちょうど一つあります。そのディレクトリの中で走らせれば、`lmf` は引数無しでもそのファイルを見つけます。

## 古い形式のディレクトリの移行 {#migration-of-an-existing-legacy-directory}

```bash
# inside an old directory containing ctrl.<sname> [+ GWinput]
Legacy2toml.py <sname>
# produces ctrlg.<sname>.toml
# resume normal workflow (lmf, gwsc, job_eps, ...)
```

`Legacy2toml.py` は `ecalj/SRC/exec/Legacy2toml.py` にあります（`InstallAll.py` がバイナリと一緒に `~/bin/` に入れます）。
何度走らせても同じ結果になります。変換で生き残らないコマンドラインの `-v` の上書きがあれば、`[INFO] / [WARN] / [ERROR]` の診断を出します。

### `esm_input.dat`（スラブ） {#esm-input-dat-slabs}

別ファイルだった `esm_input.dat` は 2026-09-16 にやめ、`ctrlg.<sname>.toml` の `[esm]` 節にしました。Fortran はこのファイルを読みません。
残っていると `lmf` は止まり、変換のスクリプト `ctrlg_absorb.py <sname>` の名前を示します。このスクリプトは、ctrlg に `[esm]` がまだ無ければ足し、
元のファイルを `esm_input.dat.bk` に移して、設定をどこへ移したかを頭に書きます。`Legacy2toml.py` も古い形式の変換のときに同じことをします。
同じスクリプトは、残っている `PB.<sname>.toml`（原子ごとの積基底の表。2026-09 までは別ファイル）も `[product_basis]` に取り込みます。
GW のバイナリはこのファイルがあっても同じように止まります。[ESM](./esm) を見てください。

### 古い `ctrlg.<sname>.toml`: `ctrlg_update.py` {#an-older-ctrlg-sname-toml-ctrlg-update-py}

前の版の変換スクリプトが書いた `ctrlg.<sname>.toml` には、プログラムが止まる原因になるものが入っていることがあります。

- `[struc]` に `nbas` / `nspec` があってはいけない（サイトと種の数は `[[site]]` / `[[spec]]` の表の数で決まる）。
- `[gw]` には `t_tetrakbt`（K）が要り、`SmearX0` / `SmearX0q0` があってはいけない。
- `[gw]` には `t_sigmaw`、`HistBin_dw`、`HistBin_ratio`、`QforEPSau`、`EMAXforGW` が要る（2026-11-14 から）。

```bash
ctrlg_update.py ctrlg.<sname>.toml [...]    # in place; the original is kept as <file>.bak_update
```

これは `nbas` / `nspec` を消し（`nbas` が `[[site]]` の表の数より少なかったときは、それより後の表をコメントにする）、
`SmearX0 = s` を同じ幅の `t_tetrakbt = -T` に書き換え、キーが無いところには `t_tetrakbt = 0` を足し、ほかの必須のキーもそれまでのコードの既定値で足して、計算が変わらないようにします（[GW の設定](./gwinput)）。
`--no-backup` を付けると元のファイルの控えを残しません。規則は [lmf.md](./lmf#file-structure-sections) と [kBT](./kBT) §2 に、
`[gw]` の古いキーは [GW の設定](./gwinput) の表 1 にあります。

## 実行時の `--ctrlg:` による上書き {#run-time-ctrlg-overrides}

古い `ctrl.<sname>` に埋め込まれた `%const` の定数は、変換のときに `ctrlg.<sname>.toml` に書き込まれます。実行時の上書きは
`--ctrlg:<ドットで区切ったパス>=<値>` で行い、`m_toml_override.f90` がメモリの中で処理します（ファイルは書き換えません）。

```bash
# OLD (the programs stop with a message that names the new form)
lmf si -vnk=8 -vmetal=3
lmf si -v[bz.nkabc]=[8,8,8] -v[bz.metal]=3
lmf si --toml.bz.nkabc=[8,8,8] --toml.bz.metal=3

# NEW
lmf si --ctrlg:bz.nkabc=[8,8,8] --ctrlg:bz.metal=3
```

`--ctrlg:` の後ろは、ドットで区切った TOML のパス（`節.キー`）です。値は TOML として読みます（整数のベクトルは `[8,8,8]`、
真偽値は小文字の `true`/`false`、文字列は引用符で囲む、など）。上書きの層は、裸の識別子の値を自動で引用符で囲むので（`find` → `"find"`）、
シェルが引用符を取ってしまっても困りません。

> ⚠️ **やめた上書きの書き方（どれもヒントを出して止まる）。**
> 移行の途中の版では、いくつもの書き方を受け付けていました（`-v<NAME>=<VAL>`、`-v[<path>]=<VAL>`、`--[<path>]=<VAL>`、
> `--<a.b.c>=<VAL>`、`--toml.<path>=<VAL>`、`--pr=N`、`--time=N,M`、`--phispinsym`）。これらはすべて無くなりました。
> コマンドラインにどれかがあると、正しい書き方 `--ctrlg:<path>=<value>` を示す一行を出して止まります。
>
> - `--ctrlg:bz.nkabc=[8,8,8]` — この実行の間だけ、`ctrlg.<sname>.toml` の TOML のパス `[bz].nkabc` の値をメモリの中で直接置き換えます。
>   `%const` を介さず、ファイルも書き換えません。パスはスキーマにある本当のキーでなければいけません。
>   使った上書きはランク 0 がログに書くので、`llmf` などに出ます。
>
> 実際の使い方で変わること:
>
> - `-vmetal=3`、`-vnspin=2 -vso=0` のような古い癖は、移行のヒントを出して止まります。
>   `--ctrlg:bz.metal=3`、`--ctrlg:ham.nspin=2 --ctrlg:ham.so=0` を使ってください。
> - `Legacy2toml.py` は、変換のときに見つけた `%const` の上書きについて `[INFO] / [WARN] / [ERROR]` の診断を出し、
>   知らないうちに壊れるのを見つけやすくします。

## 移行済みのサンプル（手本に使う） {#migrated-samples-use-these-as-templates}

`ecalj/Samples/` の下の次のディレクトリは `ctrlg.<sname>.toml` を持ち、`testecalj` を通ります（全体の一覧は [Samples](./samples)）。

| ディレクトリ | 役割 |
|---|---|
| [Samples/MLOsamples/](https://github.com/tkotani/ecalj/tree/main/Samples/MLOsamples) | MuffinTin Localized Orbitals（Wannier の代わり）。試験のターゲットは 25 — [MLOsamples/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/MLOsamples/README.md) を参照 |
| [Samples/TestInstall/](https://github.com/tkotani/ecalj/tree/main/Samples/TestInstall) | インストールの検証の組。26 ターゲット（66 件の確認）を `testecalj --all` で回す |
| [Samples/EPS/](https://github.com/tkotani/ecalj/tree/main/Samples/EPS) | 誘電関数 ε(q,ω) — `EPS_Cu`、`EPS_GaAs`、`EPS_Ag`（`job_eps`） |
| [Samples/PROCAR/](https://github.com/tkotani/ecalj/tree/main/Samples/PROCAR) | ファットバンドの重み — `MgO_PROCAR`、`Ni2MnGa_L21_PROCAR` |

サンプルの木全体の索引は [Samples/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/README.md) です。

## 古い形式のサンプル {#legacy-samples}

`Samples/Legacy/` にあった古い例は、`Samples/` の下にそれぞれのディレクトリとして作り直し（反強磁性の対称性、マグノンのスペクトル、
フェルミ面、ドーピング、磁気異方性、LDA+U、有効質量、緩和、DOS など）、元のものは 2026-09-30 に消しました。一覧は ecalj の
`Samples/README.md` の表 1 にあります。古い形式の自分の入力は、上に書いたとおりに変換してください。

## 結果のまとめ（計算済みの参照値） {#result-summaries-worked-out-reference-values}

- [MLOsamples/README.md § Reference value: bcc Fe Fe-3d on-site W](https://github.com/tkotani/ecalj/blob/main/Samples/MLOsamples/README.md#reference-value-bcc-fe-fe-3d-on-site-w-with---mlo_feb4) — bcc Fe の Fe-3d のオンサイトの遮蔽された W はおよそ 1.5 eV（job_mloW での再現、BZ は 2×2×2 と 4×4×4）
- [MLOsamples/README.md § Regression test (testecalj Fe)](https://github.com/tkotani/ecalj/blob/main/Samples/MLOsamples/README.md#regression-test-testecalj-fe) — 対角の V/W-V を、中に書いた参照値と自動で比べる試験
- MLOsamples/README_SOC.md（`MD/mlo_notes/mlo_soc_memo.md`）— SOC を摂動として扱う版（`job_mlo_soc`）

## 自動でジョブを回す道具（`ecalj_auto`） {#auto-job-runner-ecalj-auto}

多くの物質で QSGW / GW のジョブをまとめて回すときは、
[ecalj_auto/README.md](https://github.com/tkotani/ecalj/blob/main/ecalj_auto/README.md) と
[ecalj_auto/README_slot_scheduler.md](https://github.com/tkotani/ecalj/blob/main/ecalj_auto/README_slot_scheduler.md) を見てください。
どの駆動スクリプトが今のプログラムで動くかは `MD/auto.md` に書いてあります。`jobtemplate.{kugui,ohtaka,ucgw,...}` は SLURM/PBS のクラスタ用の雛形です。

## よくある落とし穴 {#common-gotchas}

- **やめた上書きの書き方** — `lmf` は `-vnspin=2 -vso=0`、`-v[ham.so]=1`、`--toml.ham.so=1`、`--pr=50`、`--time=5,5`、
  `--phispinsym` で止まります。止まるときには、それぞれ正しい書き方 `--ctrlg:<path>=<value>` を示します。
- **古い `Worb` ブロック** — `gwinput2toml.py`（`Legacy2toml.py` の GWinput を扱う部分）は、`<Worb>` ブロックが重複していると
  以前は最後のものを採っていました。今は最初のものを採ります。
- **`mlo_emax`** — 古い `mlo_method = 0` と `1` だけが読みます。既定の `mlo_method = 4` は使いません。
- **`t_tetrakbt` などの必須のキーが無い、`SmearX0`、`nbas` / `nspec`** — プログラムが止まります。`ctrlg_update.py`（上）でファイルを今の規則に合わせます。
- **gfortran 13.3 / 14.2 のコード生成** — 知られているコード生成のバグを `m_HamPMT.f90` で回避しています（一行のおとり
  `aaa = trim(aaa) // ' '`）。master のコミット `35ee9f225` と `f693848ce` を見てください。

nvfortran 26.1 に限った話として、論理型の配列に対する `findloc` が実行時に落ちる件は、影響を受ける 7 つのソースファイルに
`use m_nvfortran, only: findloc` を足して直しました（コミット `f693848ce`）。
