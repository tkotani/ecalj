# 結晶の対称性と spglib（2026-10-02 の検討）

user の依頼（2026-10-01〜02）: ecalj の対称性の作り方を、構造（与えた胞と原子の並び）を保ったまま spglib に置き換えたい。k 点の生成にも使われ、
反強磁性の対称性（SYMGRPAF）もある。与えた胞より小さい並進がある（超格子）とこける。この文書は、今の仕組み、超格子の問題と 2026-10-02 の直し、
spglib を入れる案と順番をまとめる。コードの場所はルーチン名で書く（行番号は書かない）。

## 1. 今の仕組み

**表 1**. 対称操作を作る流れ（`SRC/subroutines/m_mksym.f90`、`m_mksym_util.f90`。2026-10-02 の S1 で `m_symfind.f90`・`m_symderive.f90`・`m_symop_util.f90` に分けた。§4.7b）

| 段 | ルーチン | 何をするか |
| --- | --- | --- |
| 入力 | `m_lmfinit`（`SYMGRP`、`SYMGRPAF`）、`--nosym`、`--noinv`、`--afsym` | `SYMGRP` は既定 `find`（自動）か生成元の文字列。`SYMGRPAF` は反強磁性の追加の生成元 |
| 1 | `m_mksym_init` → `mksym` | 普通の対称性。`faithful=.true.`（2026-10-02） |
| 2 | `gensym` → `symlat` | 格子の点群（立方なら 48） |
| 3 | `symcry` | 各回転に、結晶を不変にする並進を 1 つ付ける（48 個の (R, t)） |
| 4 | `groupgblock`（`gensym` の中） | 群を一番大きくする生成元を貪欲に選ぶ。`sgroup` が生成元から群を作る（上限 `ngmx = 48`） |
| 5 | `symtbl`、`splcls` | 原子の置換表 `istab`、原子をクラスに分ける |
| 6 | `m_mksym`（AF） | `SYMGRPAF` があると、上下のスピンの原子を同じ種類と見なして `mksym` をもう一度呼び、「格子 + AF」の群を作る。普通の群に無い操作が AF の操作（`ngrpAF` 個） |
| 7 | `mptauof`、`rotdlmm` | 原子の写し先と格子の並進（`miat`、`tiat`）、D^l 行列 |

**使う側**（`use m_mksym` は 21 ファイル、`symops` は 40 ファイル）: k 点（`mkqp`、`bzmesh`。回転だけを使い、時間反転のために反転を足す）、
密度の対称化（`symrho` の `symprj`。**群として閉じていることを前提にする**）、力（`totfrc`）、LDA+U（`m_ldau_util`）、
GW（`m_hamindex` が `symops`・`ag` をファイルに書き、`rotwave` の `rotmatMTO` などが固有関数を回す）、MLO（`m_HamPMT`、SOC は `so3_to_su2`）。

## 2. 超格子でこける問題と、2026-10-02 の直し

- 再現: Si を立方晶の慣用胞（8 原子、fcc の中心の並進が 3 つ余分にある）で与えると、`symcry` は 48 個の (R, t) を見つけるが、`sgroup` が生成元の積で
  「同じ回転で、並進が中心の分だけ違う」要素を作り、格子の並進だけで同一視するので 48 を超え、`SGROUP: ng=6 > 48 probably bad translation` で止まる
  （ctrlgenToml の `lmchk --getwsr` の段階で）。2×1×1 の超格子は、余分の並進が回転と組まず、止まらない
- 最初の直し（同じ回転なら同じ要素とする）は**誤り**だった: lmchk は通るが、`symprj`（原子ごとの密度の対称化）が群の作用（各原子を最初の原子へ移す
  操作の数がそろう）を前提にするので、8 原子の電荷がばらばらになり NaN（user「対称化操作として足りない」）
- **採った直し**（`m_mksym_util` の `lfaithful`〔S1 から `m_symfind` の `gensym(...,faithful)`〕、`puretrans`、`distinctrot`、`sgroup` の `overflow`）: 結晶の純粋な並進 τ を探し、あるときだけ、
  生成元の並進を t + τ から選び直して、**回転がすべて異なる閉じた群**のうち大きいものを選ぶ（最初の 2 つの生成元は総当たり、最大の組から最大 20 組を
  貪欲に伸ばす）。純粋な並進そのものは対称操作として使わない（対称性を低めに使うので正しい）。AF の 2 回目の呼び出しには使わない
  （そこでは上下を入れ替える純粋な並進が AF の操作そのもの）
- 結果: Si の慣用胞で 24 個の閉じた群（48 個の Pn-3m までは届かない。2.5 秒）。lmf は対称性なし（`--nosym`）と全エネルギーが 3 μeV/胞で一致。
  純粋な並進の無い入力は処理が働かず、14 入力で lmchk の出力（操作の一覧を含む）が前後で同じ
- 残り: 純粋な並進も使う（user の決定、4.6）。超格子の GW は 2×1×1 の Si で動く（kt1 で 1 反復 3 分半）

## 3. まず Fortran の中で対称性の部分を分ける（user 2026-10-02 未明「fortran code で対称部分をしっかり分離して抜き出すのがいる」）

spglib に替える前に、**見つける部分**と、それから**導く部分**、**使う部分**を分け、受け渡す中身（契約）を一つにする。見つける部分だけを差し替えられるように。

**表 2**. 三つの層（2026-10-02 の調べ。`use m_mksym, only:` で多いもの: ngrp 15、symops 14、miat・tiat・shtvg・dlmm・ag 各 5、ngrpAF 4）

| 層 | 今のルーチン | 入力 → 出力 | 備考 |
| --- | --- | --- | --- |
| 見つける | `gensym`（`symlat`、`symcry`、`groupgblock`、`sgroup`、2026-10-02 の `puretrans`・`distinctrot`）、`psymop`（生成元の文字列）、`m_mksym` の AF の 2 回目の `mksym` | 格子、原子の位置と種類、`SYMGRP`・`SYMGRPAF` → 操作 (g, ag)、AF の操作の数 `ngrpAF` | spglib に差し替える所。純粋な並進の扱い（閉じた部分群を選ぶ）もここ |
| 導く | `symtbl`（`istab`）、`splcls`（クラス `iclasst`、`oics`）、`mptauof`（`miat`、`tiat`、`invgx`、`shtvg`）、`rotdlmm`（`dlmm`）、`grpgen`（k 点用に反転を足す → `npgrp`） | 操作 + 構造 → 置換表、クラス、写し先、D 行列 | 純粋な計算。`mptauof` は GW 側（`m_hamindex0`、`m_zmel`、`main_hwmatK`）でも別に呼ばれ、`rotdlmm` は `m_procar`、`rotcg` でも呼ばれる → 一か所で作り、`__HAMindex` で渡すのに揃える |
| 使う | lmf（`symrho`、`totfrc`、`mkqp`、`m_ldau_util`、`m_bandcal`、`smves`、`supot` など 21 ファイル）、GW（`m_hamindex` が `__HAMindex` から `symops`・`ag` を読む）、MLO（`m_HamPMT`、`m_mlo_ovlppair`） | `m_mksym` の `protected` 変数 | 変えない |

**手順**（どれも等価変換で、Samples の lmchk の出力（操作の一覧）と試験の組が前後で同じことで確かめる）:

1. 見つける部分を新しい module `m_symfind`（singleton、`protected` で g、ag、ng、ngAF を出す）に移す。`m_mksym_util` には導く部分だけを残す
2. 導く部分を `m_mksym` の初期化の中で一度だけ計算し、GW 側は `__HAMindex` から読む（`mptauof` の重複を外す）
3. `m_symfind` に「外から操作を受け取る」口を足す（ctrlg の `[symmetry]` か、Python のラッパーから配列で）。spglib の操作を Python で
   作り（純粋な並進があれば閉じた部分群を選ぶ）、その口から入れる。ecalj の自前の方法と結果を比べられるように、両方を残す

## 4. 設計（2026-10-02、user「対称化部分をクリーンにモジュール化して、中身は spglib に置き換える（インターフェースが要る）。
対称性の記述法は ecalj 独自のものにこだわらず、spglib の標準記法に。でたらめにならないように作戦・設計をしっかりやってから。commit とテストを重ねて安全に」）

### 4.1 目的と、しないこと

- 目的: (1) 対称性を「見つける」部分を一つの module に閉じ込め、spglib に替えられるようにする。(2) 入力と出力の対称性の書き方を spglib の標準
  （整数の回転行列と分数座標の並進の組、空間群の番号・国際記号・Hall 番号、`symprec`）にする。(3) 超格子（純粋な並進がある胞）と反強磁性も、同じ仕組みで扱う
- しないこと: 対称操作を「使う」側（21 ファイル）の計算の中身は変えない。受け取る配列の形と意味も、移行の間は今と同じにする（3 節の表 2 の「使う」層）

### 4.2 module の形（singleton、`MD/ecaljclaude.md` のコーディング規約に従う。4.5 も見る）

| module | 役目 | 公開するもの（`protected`） | 入力 |
| --- | --- | --- | --- |
| `m_symfind`（新） | 見つける。backend を一つ選ぶ: `spglib`（外で作った操作を受け取る）、`ecalj`（今の `gensym` を移したもの）、`input`（`[symmetry] ops` に書いた操作） | `nsym`、`rotf(3,3,nsym)`（整数、格子の基底での回転）、`transf(3,nsym)`（分数座標）、AF の操作の数と印、空間群の番号・記号、`symprec`、どの backend か | 格子、原子の位置と種類、（AF なら）スピンの向き |
| `m_symderive`（`m_mksym_util` から分ける） | 導く。純粋な計算 | なし（ルーチンだけ: `symtbl`、`splcls`、`mptauof`、`rotdlmm`、反転の追加） | 操作と構造 |
| `m_mksym`（今のまま） | 使う側への窓口。`m_symfind` の結果をデカルト座標の `symops`・`ag` に直し、`m_symderive` で導いた `istab`・`miat`・`tiat`・`dlmm`・クラスを持つ | 今の公開変数（`ngrp`、`symops`、`ag`、`ngrpAF` ほか）はそのまま | — |

GW 側は `__HAMindex`（`m_hamindex0` が書く）から読む今の形のまま。`mptauof` を GW 側で別に呼んでいる 3 か所は、`__HAMindex` に入れて読む形に揃える（3 節の手順 2）。

### 4.3 spglib との口: 別ファイル `symmetry.json`（user 2026-10-02 の決定）

- Python の道具 `symfind.py <sname>`（新。`ctrlgenToml.py` と、lmf を回すスクリプトの前段から呼ぶ）が spglib（`get_symmetry_dataset`、AF は
  `get_magnetic_symmetry_dataset`）を回し、**`symmetry.json`** に書く。ctrlg には書かない（user「TOML は人間用、計算機には JSON」）
- `symmetry.json` には、**ctrlg の構造に関わる部分（`[struc]`、`[[site]]`、種類とスピンの向き、`symprec`）のハッシュ**を入れる。Fortran は読むときに
  今の ctrlg から同じハッシュを作って照らし、違えば止める（「`symfind.py <sname>` をもう一度」と言う）。このやり方は lmfa などの初期条件にも使う
  （メモリ `feedback_sidecar_json_hash`）
- ハッシュは Python と Fortran で同じ値になるように、正規化した文字列（例えば分数座標を `%.10f` で並べ、種類の名前、`symprec`）の SHA-256 にする。
  Fortran 側の SHA-256 は小さなルーチンを書くか、Python が正規化した文字列そのものも JSON に入れて Fortran は文字列を比べる（簡単。こちらを先に）
- ビルドには何も足さない。spglib の C を結ぶ案（旧案 II）は、要るとわかるまでしない

### 4.4 `symmetry.json` の形（単位はキーで明記。alat は使わない）

```json
{
  "format": "ecalj-symmetry-1",
  "made_by": "symfind.py (spglib 2.6.0)", "date": "2026-10-02T05:10",
  "ctrlg": "ctrlg.si8.toml",
  "structure_key": "<正規化した構造の文字列>", "structure_sha256": "...",
  "symprec": 1e-5, "symprec_unit": "angstrom",
  "spacegroup": {"number": 227, "international": "Fd-3m", "hall_number": 525},
  "lattice_unit": "fractional (columns of plat)",
  "n_operations": 192, "n_pure_translations": 4,
  "operations": [
     {"rotation": [[1,0,0],[0,1,0],[0,0,1]], "translation": [0.0,0.0,0.0], "time_reversal": false},
     ...
  ],
  "rotation_unit": "integer matrix in the basis of plat (r_frac' = R r_frac + t)",
  "translation_unit": "fractional coordinates of plat",
  "equivalent_atoms": [0,0,0,0,0,0,0,0]
}
```

- 操作の並びは spglib の順（恒等が先頭）。Fortran はこの並びのまま使う
- AF（磁性）: 各原子のスピンを spglib に渡し、`time_reversal: true` の操作を AF の操作として使う（今の `SYMGRPAF` を置き換える）
- 旧 `SYMGRP`（`find` や生成元の文字列）・`SYMGRPAF` は移行の後に廃止し、旧書式から `symmetry.json` を作る変換を付ける

### 4.5 Fortran 側（すり替えの口だけを持つ。高い機能は Python 側）

- `m_symfind` は backend を一つ選ぶだけ: `json`（`symmetry.json` を読み、ハッシュを照らす）、`ecalj`（今の `gensym` を移したもの。比べ役）、`none`（恒等だけ）。
  空間群の判定、超格子の扱い、AF の磁気対称性は持たない（Python 側）
- 公開するもの: `nsym`、`rotf(3,3,nsym)`（整数）、`transf(3,nsym)`（分数座標）、`trev(nsym)`（AF の印）、`spgnum`。`m_mksym` がデカルト座標の
  `symops`・`ag` と、導く量を作る（使う側は変えない）

### 4.6 純粋な並進も対称操作として使う（user 2026-10-02 の決定）

超格子では spglib は |点群| × |純粋な並進| 個の操作を返し、同じ回転の操作が複数ある。これをそのまま使うために、Fortran の使う側で次を確かめて直す。

**表 4**. 確かめる所（2026-10-02 の調べ。上限 48 を固定で持つのは `m_mksym`・`m_mksym_util`（S1 から `m_symfind`）の `ngmx` だけで、使う側の配列は `ngrp` で取る）

| 所 | 何を確かめるか |
| --- | --- |
| `m_mksym`、`m_symfind`（旧 `m_mksym_util`） | `ngmx = 48` を外す（操作の数で配列を取る）。印字の `48i3` などの書式 |
| k 点（`mkqp`、`bzmesh`） | 同じ回転が重なっても、既約な k と重みが正しいか（回転だけを使うので重なりを除いてよいはず） |
| 密度の対称化（`symrho` の `symprj`、`totfrc`） | 群なので正しいはず（2 節で分かった前提を満たす）。純粋な並進で原子の密度が平均される |
| LDA+U（`m_ldau_util`）、SOC（`so3_to_su2`） | 同じ回転の操作が複数あるときの D 行列とスピンの回転 |
| GW（`m_hamindex0` の `__HAMindex`、`rotwave`、`eibzgen`、`m_zmel`） | 既約な q から全体への写し、EIBZ の小群に純粋な並進が入ったときの位相 |
| MLO（`m_HamPMT` の q を保つ操作での平均、`m_sigmlo`） | 純粋な並進での H(k) の平均（並進の位相の付き方） |

**S4 の要点**（2026-10-02 05:40 の調べ）: `mptauof`（`m_symderive`）は渡された並進 `ag` を使わず、回転 `am` ごとに原子の対 (ib1, ib2) から
並進 `tran = bas(ib2) − am·bas(ib1)` を試し、全原子が同じクラスの原子に（格子ベクトル ±3 の範囲で）移る最初のものを `delta`（= `shtvg`）にする。
逆操作 `invgx` は回転だけ（転置）で探す。GW（`__HAMindex0` の `shtvg`、`m_zmel` の `miat`・`tiat`）と MLO はこの `delta` を使う。
- 基本格子（純粋な並進が無い）では、回転ごとに正しい並進は格子を法に一つなので、`delta` と `ag` は格子ベクトルだけ違い、問題ない
- 今の超格子（暫定の直しで回転が重ならない 24 操作）では、`delta` は `ag` と純粋な並進の分だけ違いうる。どれも結晶の対称操作なので、一つずつの
  波動関数の回転としては正しい（2×1×1 の Si の GW が通ったのと合う）。ただし `{(am, delta)}` が群として閉じている保証は無い
- 同じ回転の操作が複数あると、`mptauof` は同じ `delta`・`miat` を返し、`invgx` も最初の一つを指す。S4 では `mptauof` に `ag` を渡してそれを
  使い（探すのは格子ベクトルのずれ `tiat` だけ）、`invgx` は回転と並進の両方で探すように直す。基本格子でも `delta` が格子ベクトルだけ変わりうるので、
  GW の出力は位相の約束がずれる（物理は同じ）。試験の組で確かめる

確かめ方: Si の慣用胞（8 原子、192 操作）と 2×1×1 の超格子で、対称性なしの計算（`--nosym`）と全エネルギー・バンド・QSGW の準粒子のエネルギーが
丸めの範囲で一致し、基本格子の Si の値（1 原子あたり）とも合うこと。

### 4.7 手順（一段ごとにコミットし、前の段と比べて同じであることを確かめる）

| 段 | 中身 | 確かめ（通らなければ次に進まない） |
| --- | --- | --- |
| S0 | `symfind.py`（書くだけ）と `symcheck.py`: Samples の全部の ctrlg（172）で spglib を回し、lmchk の操作の一覧と集合として比べる。Fortran は変えない | 違う入力を一つずつ調べる（許容、超格子、AF）。ここで `symprec` を決める |
| S1 | Fortran の分離（等価変換）: `m_symfind`（今の `gensym` 一式を `ecalj` の backend として）と `m_symderive` | 172 入力で lmchk の出力（操作の一覧）がビット単位で同じ。試験の組（inputs、install、afsym、affix、mlo、mloqsgw） |
| S2 | `mptauof` などの重複を外し、GW 側は `__HAMindex` から読む | 同上と gwall がビット単位で同じ |
| S3 | `m_symfind` に `json` の backend（読む、ハッシュを照らす）を足す。`symmetry.json` が無ければ今までどおり | S0 で同じ操作になる入力で、`json` と `ecalj` の全部の試験の組が丸めの範囲で同じ |
| S4 | 純粋な並進を使う（表 4）。`ngmx` を外す | Si の慣用胞・2×1×1 で表 4 の確かめ。基本格子の入力は変わらない |
| S5 | AF を spglib の磁気対称性で（`time_reversal`） | afsym、affix、`Samples/AFsymmetry`。NiO・Fe2O3 で今の `SYMGRPAF` と比べる |
| S6 | 既定を `json`（`symfind.py` を lmf の前段で自動）に。旧 `SYMGRP`・`SYMGRPAF` を廃止（変換の道具）。ecaljdoc。3 台（t14、kt1、kr7）で試験 | 3 台で全部の試験の組 |

### 4.7a S0 の結果（2026-10-02 05:15 終了、`TOOLS/symcheck_samples.sh ~/work/symcheck_s0`、symprec 1e-5 Å）

Samples の 172 入力: **166 で ecalj と spglib の操作が集合として同じ、6 で ecalj が spglib の部分群、食い違い 0**。lmchk が止まった入力は無い。
部分群の 6 つ（`LDAU/ReN` の cgdn・cprn、`MLOsamples/GdION`、`MLOsamples/SmP`、`TestInstall/eras`、`TestInstall/felz`）は、どれも入力の `SYMGRP` で
対称性をわざと下げている（SOC の磁化の向き、4f の占有、LDA+U の秩序）。spglib は構造だけから全部を出すので、**入力で対称性を下げる口**が要る:
磁化の向き（非共線のスピンを spglib の磁気対称性に渡す）か、部分群の指定（生成元、または操作の番号の一覧）。S5 で AF と一緒に設計する。

### 4.7b S1 の分け方（2026-10-02 05:22、等価変換）

`m_mksym_util.f90` を三つに分け、元のファイルは trash へ。DAG は `m_symop_util` ← `m_symderive` ← `m_symfind` ← `m_mksym`。

**表 5**. S1 の module

| module | ルーチン | 役目 |
| --- | --- | --- |
| `m_symop_util` | `asymop`（`parsvc2`）、`spgcop`、`spgprd`、`spgeql`、`grpeql`、`latvec` | 一つの操作の道具（積、等しいか、格子ベクトルか、記号） |
| `m_symderive` | `symtbl`、`splcls`、`grpgen`、`mptauof`、`rotdlmm`、`iclbsjx` | 操作から導く（状態を持たない）。GW 側（`m_zmel`、`m_hamindex0`、`main_hwmatK`）と `rotcg`・`m_procar` もここを `use` する |
| `m_symfind` | `gensym`、`sgroup`、`psymop`、`parsop`、`parsvc`、`skipbl`、`symlat`、`csymop`、`symcry`、`puretrans`、`distinctrot` | 見つける（ecalj の backend）。`faithful` は module 変数をやめて `gensym` の引数に |
| `m_mksym` | `mksym`（private に移した）、`m_mksym_init` | 窓口 |

### 4.7c S2 の判断（2026-10-02 05:31）

`mptauof` を呼ぶ所は 4 つ: `m_mksym_init`（群全体と AF の操作）、`m_hamindex0_init`、GW 側の `m_zmel` の `mptauof_zmel`、`main_hwmatK`。
`m_hamindex0_init` は `m_mksym_init` と同じ引数（操作 1:ngrp、同じ位置とクラス）で作り直していたので、`m_mksym` の `miat`・`tiat`・`invgx`・`shtvg` の写しに替えた。
GW 側の二つは `__HAMindex0` から読んだ操作のうち、呼ぶ側が決めた一部（`hgw` の W の段は単位元だけ、`ng=1`）で呼ぶので、そのまま残す。
どれも同じ `mptauof`（`m_symderive`）を通るので、S4 で同じ回転の操作が重なっても、直すのは `mptauof` の一か所で済む。

### 4.7d S3 の形（2026-10-02 05:39。§4.3・§4.4 の案から変えた所も）

- ファイル名は **`symmetry.<sname>.json`**（§4.3 の `symmetry.json` から変えた）。一つのディレクトリに ctrlg が二つあることがある（`Samples/LDAU/ReN`）。
  ファイル名に ID を残す決まり（`MD/ecaljclaude.md`「ファイル命名規約」）にも合う
- 構造の照らし方: JSON に構造を数値で入れ（`structure`: `alat_bohr`、`plat_alat`、`species`、`frac`）、Fortran（`m_symfind` の `symfind_json`）は
  今の ctrlg と数値で比べる（alat・plat は相対 1e-6、分数座標は 1 を法に 1e-6、種の名前は一致）。§4.3 の「正規化した文字列を比べる」は、Fortran で
  Python と同じ書式・丸め（`%.10f` と `f0.10` の先頭の 0 など）を作るのが脆いのでやめた。文字列とその SHA-256 は人と Python のために残す
- JSON の読み手は toml-f の試験用の `json_lexer.f90`・`json_parser.f90` を `SRC/subroutines/tjson_*.f90` に写したもの（`json_load` で toml-f の表になり、
  ctrlg と同じ `get_value` で読める。`SRC/external/toml-f/VENDORED_FROM.txt`）。lmfa の初期条件などの JSON にも使える
- 使う条件（`m_mksym` の `mksym`）: ファイルがある、`symgrp` が既定の `find`、結晶の群（AF の二回目の呼び出しではない）、`SYMGRPAF` が無い。
  使えないとき（純粋な並進がある → S4、時間反転の操作 → S5）は理由を印字して gensym。環境変数 `ECALJ_SYMFIND=ecalj` で gensym に固定（比べるため）
- `SYMGRPAF` の入力で使わない理由: gensym は見つけた生成元を `ssymgr` に書いて返し、`m_mksym_init` の AF の呼び出しはそれと AF の操作から群を作る。
  json のときは `find` が残り、AF の群を別の道（`find` と AF の操作）で探すことになって並びが変わる（NiO で確かめた）。S5 で AF ごと spglib にする
- 操作の並びは spglib の順（恒等を先頭に置き直す）。gensym とは並びが違うだけで集合は同じ（Fe・NiO の lmchk で確かめた）。並びに依る出力（既約な k の代表など）は
  変わりうるので、S3 の確かめは試験の組（数値の比べ）で行う
- 古い JSON（ctrlg の位置を 0.001 動かした）は「does not match the ctrlg (a site position). Run symfind.py <sname> again」で止まる。
  超格子の Si8（192 操作、純粋な並進 4）は「not used: the cell has pure translations」で gensym に戻る

### 4.7e S3・S4 の確かめで分かったこと（2026-10-02 06:06）

**表 6**. 試験の組を json の口で（kt1、`91f9bdcaf`、CPU、lmf などの前に `symfind.py` を走らせる包み、05:40〜05:59）

| 組 | 結果 | 備考 |
| --- | --- | --- |
| inputs 172、afsym 4、affix 12、install 66、mlo 45、mloqsgw 5、procar 5、FermiSurface〜kBT_scanT | PASS | json が使われたことはログの `space group from symmetry` で確かめた |
| BoltzTraP（Si） | TEST 2 FAIL | `si.struct.boltztrap` は操作の一覧。整数の行列として集合は同じで、並びだけが違う。S6 で json を既定にするとき参照を作り直す |
| AtomDimer（N₂） | STOPPED | 試験が `--ctrlg:site.n.pos=` で位置を上書きしていた。照らし合わせが止めたのは正しい動き。`symfind.py` が同じ上書きを受け取るようにした（`3572a0859`）。S6 で lmf の前段から呼ぶときも引数をそのまま渡す |

- Si8（Si の 8 原子の立方胞、[bz] 8×8×8、gfortran、t14）: json の口で 192 操作（純粋な並進 4）と、gensym の 24 操作（`ECALJ_SYMFIND=ecalj`）で、
  収束した ehf・ehk・sev が表示の桁まで一致（ehf −62922.978555 eV、sev −64.768559 eV）。既約な k は両方 35。対称性なし（`--nosym`、512 k 点、
  06:15 終了）は ehf −62922.978553 eV（差 2 μeV）、sev −64.768595 eV
  （以前に見た −62923.0390 eV は収束していない反復の値だった）
- S4a で LiTi₂O₄（`Samples/kBT/LiTi2O4`）の lmchk が止まった: gensym の `ag` は格子ベクトル何個分も長いことがあり、`mptauof` の格子のずれの探索
  （各方向 −3〜3）が届かなかった。分数座標を丸めて直接求めるように直した（一致すれば同じ格子ベクトル）

### 4.7f S5 の形（2026-10-02 06:16）

- `symfind.py`: `[[site]] af` が 0 でない原子があれば、AF の対（+k と −k）を同じ種類にまとめ、磁気モーメント ±1（ほかは 0）で
  `spglib.get_magnetic_symmetry_dataset` を回す。時間反転なしの操作（スピンを保つ）が結晶の群で、`get_symmetry_dataset`（上下を別の種類として）の群と
  同じであることを確かめてから書く。時間反転つきの操作（上下を入れ替える）が AF の操作。並びは恒等、結晶の群の残り、AF の操作。`structure.af` も書く
- Fortran（`m_mksym` の `mksym`）: 一回目（結晶の群）は時間反転なしの操作だけを読む。それを使ったら（`jsonused`）、二回目（AF の対をまとめた種類の群）は
  全部の操作を読む。AF の操作を選ぶのは今までどおり `m_mksym_init`（二回目の群のうち一回目に無いもの）。構造の照らし合わせに `af` の印も加えた
- `SYMGRPAF` は AF の型を入れる印として残る（何を書いても json の操作が使われる）。S6 で `af` の印から決めるようにして廃止
- NiO（`Samples/AFsymmetry/NiO`）: 結晶の群 12、AF をまとめた群 24 が、gensym の群と集合として一致（lmchk、格子を法として）。
  spglib の磁気空間群は UNI 1332、型 4（反ユニタリの操作が並進を伴う黒白群）

### 4.8 決めたこと（user 2026-10-02）

1. 口は別ファイル `symmetry.json`（ctrlg に書かない）。ctrlg のハッシュを入れて整合性を保つ
2. 単位はキーで明記。alat は使わない
3. 純粋な並進も対称操作として使う
4. Fortran は backend をすり替えられるようにする。高い機能は Python 側に置き、Fortran には持たせない

## 5. 未確認・決めていないこと

- （済み 2026-10-02 05:11）2×1×1 の Si の超格子の GW（`hgw`）が手元（t14、4 コア）で 30 分で終わらなかったのは遅いだけ: kt1 の 32 コアでは 1 反復が 3 分半で終わった（`/mnt/data1/symtest_si2`、今のコード）
- spglib の許容（`symprec`）と ecalj の `toll = 1e-4` の合わせ方
- AF で、spglib の磁気空間群の操作と `SYMGRPAF` の操作が一致するか（NiO、Fe2O3 などで比べる）
