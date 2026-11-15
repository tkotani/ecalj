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

### 生成されるファイルの体裁 {#file-layout}

`ctrlgenToml.py`、`gwinit`、`Legacy2toml.py` と整形の道具 `pylib/toml_tidy.py` は、どれも同じ体裁のファイルを作ります。
空行で区切るのは大きな見出し（`# === X ===` のコメントの並びと、一つだけの `[section]` の見出し）だけです。
節の中（`[[site]]` や `[[spec]]` の表どうし、複数行のブロックとその前後、積基底の三つの表）は `# ----` の罫線で区切ります。
`[blocks]` の節は、ほかに置き場のない古い `<...>` のタグ（`QPNT`、`QforEPSL`、`hrotr`）があるときだけ現れ、複数行の文字列として持ち、GW のプログラムがそのまま読みます。
節の順は `[gw]` `[mlo]` `[blocks]` `[product_basis]` で、`[product_basis]` がいつも最後です（`toml_tidy.py` がこの順にそろえます）。

## 古い形式のディレクトリの移行 {#migration-of-an-existing-legacy-directory}

```bash
# inside an old directory containing ctrl.<sname> [+ GWinput]
Legacy2toml.py <sname>
# produces ctrlg.<sname>.toml
# resume normal workflow (lmf, gwsc, job_eps, ...)
```

`Legacy2toml.py` は `ecalj/SRC/exec/Legacy2toml.py` にあります（`InstallAll.py` がバイナリと一緒に `~/bin/` に入れます）。
何度走らせても同じ結果になります。変換で生き残らないコマンドラインの `-v` の上書きがあれば、`[INFO] / [WARN] / [ERROR]` の診断を出します。

### 古いキーと TOML のパスの対応 {#legacy-toml-map}

古い `Category_Token` のキーの多くは、小文字にすると `[節.キー]` にそのまま対応します。名前や形が変わったものは次のとおりです（表 1）。
キーの意味は [古い ctrl の文書](./lmf_legacy_ctrl) で引けます。

表 1. 古いキーと TOML のパス

| 古い `Category_Token` | TOML の場所 | 注 |
|---|---|---|
| `STRUC_ALAT` / `STRUC_PLAT` | `[struc].alat` / `.plat` | |
| `SITE_ATOM=`, `POS=`, `XPOS=` | `[[site]] atom = ...; pos = [...]`（または `xpos`） | 原子ごとに `[[site]]` の表を一つ |
| `SPEC_ATOM=`, `Z=`, `R=`, `RSMH=`, `EH=`, ... | `[[spec]] atom = ...; z = ...; r = ...; rsmh = [...]` | 種ごとに `[[spec]]` の表を一つ。`R/W=` と `R/A=` は引用符付きのキー `"r/a"` |
| `BZ_NKABC` / `BZ_METAL` / `BZ_TETRA` | `[bz].nkabc` / `.metal` / `.tetra` | |
| `HAM_NSPIN` / `HAM_SO` / `HAM_XCFUN` / `HAM_GMAX` / `HAM_PWMODE` / `HAM_PWEMAX` | `[ham].nspin` など | |
| `HAM_FORCES` / `HAM_OVEPS` / `HAM_FRZWF` / `HAM_REL` | `[ham].forces` など | |
| `HAM_ScaledSigma` | `[ham].scaledsigma` | 小文字 |
| `HAM_READP` / `HAM_PNUFIX` | `[ham].readp` / `.pnufix` | `HAM_V0FIX` はコマンドラインのオプション `--v0fix` |
| `esm_input.dat`（別のファイル） | **`[esm]`** の節 | 2026-09-16 にやめた。読まないので `ctrlg_absorb.py <sname>` で変換する（[ESM](./esm)） |
| `ITER_NIT` / `ITER_MIX` / `ITER_CONV` / `ITER_CONVC` / `ITER_UMIX` / `ITER_TOLU` | `[iter].nit` / `.mix` / ... | |
| `SYMGRP` / `SYMGRPAF` | 最上位の `symgrp = "..."`（と `symgrpaf`） | 節の中には入れない |
| `EWALD_TOL` | `[ewald].tol` | ほとんど触らない |
| `DYN_MODE` / `DYN_NIT` / `DYN_HESS` / ... | `[dyn].mode` など | |
| `IO_VERBOS` / `IO_TIM` | 最上位の `verbose` / `time` | `[io]` から出した（節はなくした） |
| `<Worb>...</Worb>` のブロック | `[mlo].mlo_lm = """ ... """` | 複数行の文字列 |
| `<QforEPS>...</QforEPS>` のブロック | `[gw].QforEPS = """ ... """`（`QforGW` も同じ） | |
| `n1n2n3`（GWinput） | `[gw].n1n2n3 = [k1,k2,k3]` | 整数のベクトル |
| `HistBin_dw` / `HistBin_ratio` / `niw`（GWinput） | `[gw].HistBin_dw` など | 名前は GWinput と同じ（大文字・小文字も） |
| `esmr`, `SmearX0` / `GaussianFilterX0`, `GaussSmear`, `delta`, `dw`, `omg_c`, `WgtQ0P`（GWinput） | `[gw].t_sigmaw` と `[gw].t_tetrakbt`（K。`t_tetrakbt` は必須） | 古いキーとその扱い: [GW の設定](./gwinput) の表 1、[kBT](./kBT) §3.5 の表 2 |
| 積基底の上限（`<PRODUCT_BASIS>` の `tolerance`、`lcutmx`） | `[product_basis].pb_tolerance` / `.pb_lcutmx` | |
| 積基底の原子ごとの表（`nlx`、`valence`、`core`） | `[product_basis].nlx` / `.valence` / `.core`（ファイルの最後。2026-09 より前は別の `PB.<sname>.toml`） | |

全体のスキーマは [`SRC/exec/pylib/ctrl_schema.py`](https://github.com/tkotani/ecalj/blob/main/SRC/exec/pylib/ctrl_schema.py) にあります（変換の正本）。
注釈と行末のコメントは [`SRC/exec/pylib/toml_comments.py`](https://github.com/tkotani/ecalj/blob/main/SRC/exec/pylib/toml_comments.py) が付けます。

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
`--no-backup` を付けると元のファイルの控えを残しません。規則は [lmf.md](./lmf#file-structure-sections) と [kBT](./kBT) §3 に、
`[gw]` の古いキーは [GW の設定](./gwinput) の表 1 にあります。

## 実行時の `--ctrlg:` による上書き {#run-time-ctrlg-overrides}

一回の実行だけキーを変えるには `--ctrlg:<節.キー>=<値>` を付けます（例 `lmf si --ctrlg:bz.nkabc=[8,8,8]`）。ファイルは書き換えません。書き方の細部、ファイルに無いキーの扱い、やめた古い書き方（`-v<NAME>=<VAL>` など。使うと止まります）は [コマンドラインのオプション](./cmdopts#run-time-ctrlg-overrides) にあります。

## サンプル {#migrated-samples-use-these-as-templates}

手本になる `ctrlg.<sname>.toml` のサンプルの一覧は [Samples](./samples) にあります。古い `Samples/Legacy/` の例も、そこに書いたディレクトリとして作り直してあります。古い形式の自分の入力は、上に書いたとおりに変換してください。

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

- **古い `Worb` ブロック** — `gwinput2toml.py`（`Legacy2toml.py` の GWinput を扱う部分）は、`<Worb>` ブロックが重複していると
  以前は最後のものを採っていました。今は最初のものを採ります。
- **`mlo_emax`** — 古い `mlo_method = 0` と `1` だけが読みます。既定の `mlo_method = 4` は使いません。
- **`t_tetrakbt` などの必須のキーが無い、`SmearX0`、`nbas` / `nspec`** — プログラムが止まります。`ctrlg_update.py`（上）でファイルを今の規則に合わせます。
- **gfortran 13.3 / 14.2 のコード生成** — 知られているコード生成のバグを `m_HamPMT.f90` で回避しています（一行のおとり
  `aaa = trim(aaa) // ' '`）。master のコミット `35ee9f225` と `f693848ce` を見てください。

nvfortran 26.1 に限った話として、論理型の配列に対する `findloc` が実行時に落ちる件は、影響を受ける 7 つのソースファイルに
`use m_nvfortran, only: findloc` を足して直しました（コミット `f693848ce`）。
