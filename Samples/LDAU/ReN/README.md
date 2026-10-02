# ReN: 希土類窒化物 GdN, PrN の LDA+U とスピン軌道相互作用

> 説明の本体: ecaljdoc の [manual/UsageDetailed.md](../../../ecaljdoc/manual/UsageDetailed.md)「LDA+U」（サイト https://ecalj.github.io/ecaljdoc/manual/UsageDetailed）。サンプルの一覧は [manual/samples.md](../../../ecaljdoc/manual/samples.md)。

岩塩構造の希土類窒化物 GdN（`cgdn`）と PrN（`cprn`）で、4f 殻に LDA+U をかけ、スピン軌道相互作用を
$L_z S_z$ の形（`so = 2`）で入れて自己無撞着に解く。LDA+U は 4f の占有の初期値によって違う解に落ちるので、
初期値を `occnum.<sname>` でフント則に合わせて与える。`Samples/TestInstall/gdn` が LDA+U だけを試すのに対し、
ここでは LDA+U・スピン軌道相互作用・`occnum` を合わせて使い、スピンと軌道のモーメントを見る。

## 入力の要点

`ctrlg.cgdn.toml`, `ctrlg.cprn.toml` は `ctrlgenToml.py <sname> --skipgw` の出力で次を変えたもの。

- `[[spec]] idu = [0, 0, 0, 2]`: 4f に LDA+U（2 = FLL、fully localized limit）。`uh`, `jh` は Ry 単位（表 1）。
- `[ham] nspin = 2`, `so = 2`。
- `symgrp = "r4z"`: z 軸まわりの 4 回回転だけを使う。`occnum` の占有（$L_z$ の固有状態）は立方対称ではないため。
- `[ham] pwmode = 1`, `[bz] nkabc = [6, 6, 6]`（旧サンプルは 8×8×8）。

`occnum.<sname>` は 4f の占有の初期値で、`#` の行は注釈、続く 2 行が spin 1 と spin 2、
各行は $m=-3,\dots,3$（球面調和関数）の占有である。`dmats.<sname>` があるとそちらが優先して読まれる
（計算をやり直すときは `dmats.<sname>` と `rst.<sname>` を消す）。

表 1: U, J と占有の初期値。

| 物質 | 4f | uh (Ry) | jh (Ry) | spin 1 の占有 ($m=-3..3$) | 初期値の $L_z$ |
|---|---|---|---|---|---|
| GdN | f⁷ | 0.676 | 0.0882 | 1 1 1 1 1 1 1 | 0 |
| PrN | f² | 0.535 | 0.0692 | 1 1 0 0 0 0 0 | −5（スピンと逆向き） |

## 実行

```bash
lmfa cgdn > llmfa.cgdn
mpirun -np 4 lmf cgdn > llmf.cgdn
python3 ldau_result.py cgdn      # result.cgdn.txt

lmfa cprn > llmfa.cprn
mpirun -np 4 lmf cprn > llmf.cprn
python3 ldau_result.py cprn      # result.cprn.txt
```

試験は `Samples/LDAU` で `testecalj ReN -np 4`。手で実行するときは、このディレクトリを別の場所に写してその中で行う
（`testecalj` はこのディレクトリをまるごと写して使うので、`rst.*` や `dmats.*` をここに残さない）。

## 結果の見方

`ldau_result.py` は `save.<sname>`、`mmom.<sname>.chk`、`orbitalmom.chk`、`dmats.<sname>` から
`result.<sname>.txt` を作る。`orbitalmom.chk` は名前に物質名が付かず次の計算で上書きされるので、`ldau_result.py` は
最初に呼ばれたとき `orbitalmom.<sname>.chk` に写す（lmf の直後に実行する）。

表 2: 得られた値（6×6×6、gfortran、MPI 4 並列、2026-09-30）。モーメントは μB、MT は MT 球の中の値。

| 量 | GdN | PrN |
|---|---|---|
| 全スピンモーメント mmom | 6.9999 | 2.0122 |
| 希土類の MT 内スピンモーメント | 7.0564 | 2.0273 |
| 希土類の軌道モーメント | −0.0063 | −4.8911 |
| N の MT 内スピンモーメント | −0.1294 | −0.0726 |
| 4f の占有（spin 1 / spin 2 の和） | 6.9930 / 0.0246 | 2.0422 / 0.0448 |
| 全エネルギー ehk (eV) | −307996.977114 | −252676.330746 |
| 反復の数 | 12 | 27 |

- GdN は 4f が半分詰まり軌道モーメントはほぼ 0、PrN は $m=-3,-2$ の二つが詰まったままで軌道モーメントが −4.89 μB
  （スピン 2.03 μB と逆向き）になる。初期値の占有がそのまま残っていることを表 2 の 4f の占有で確かめられる。
- PrN は密度行列の収束が遅く、終わりのほうでも ehk が 1 反復あたり 1e-4 eV ほど動く。
- 旧サンプル（8×8×8）の PrN は mmom = 2.0016、Pr の軌道モーメント −4.892 μB。

## 注意

- `idu` に 10 を足した値（旧サンプルの `idu = 12`。`sigm` があるとき U を切る指定）は、2026-09-30 の版では `sigm` が無いときに
  二重計数の補正が入らない。PrN では mmom = 1.43、軌道モーメント −2.78 μB の別の解になった。`idu = 2` と書く。
- `pwmode = 1` にしてあるのは、2026-03-30 から 2026-09-30 までの版が `pwmode = 11`（基底の数が k ごとに違う）で LDA+U の密度行列と
  軌道モーメントを誤って計算するため（全エネルギーが桁違いになり、並列数で値が変わる）。2026-09-30 に直した版は `pwmode = 11` でも計算できる。

## 計算時間

MPI 4 並列（16 コアの計算機）で GdN が約 45 秒、PrN が約 1.8 分、合わせて 2.6〜3.0 分。

## 試験の内容（`test.py`）

`result.cgdn.txt` と `result.cprn.txt`（表 2 の量と 4f の軌道ごとの占有）を参照と比べる（許容 1e-3）。

## ほかの希土類窒化物

`INIT/` に、15 種の希土類窒化物（CeN〜LuN と LaN）の構造 `ctrls.<sname>` と、フント則に合わせた 4f の占有の初期値 `occnum.<sname>`、
バンドを描く経路 `syml.ren` がある。ほかの窒化物は、GdN・PrN と同じ手順で計算する: `ctrls.<sname>` から `ctrlgenToml.py <sname> --skipgw` で
入力を作り、「入力の要点」の項目（`idu`・`uh`・`jh`、`nspin = 2`、`so = 2`、`symgrp = "r4z"`、`pwmode = 1`）を書き、`occnum.<sname>` を同じ場所に置く。

## 由来

`Samples/Legacy/ReNcub`（15 種の希土類窒化物を `job` で順に計算するサンプル）から GdN と PrN を取り出し、
今の入力形式と命令に直したもの。
