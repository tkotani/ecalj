# Si のドーピング（分数の核電荷、一様な背景電荷、固定スピンモーメント）

単位胞の電子数を整数からずらす方法の例です。方法は二つあり、組み合わせることもできます。
一つは原子の核電荷 `[[spec]] z` を分数にする方法（仮想結晶）、もう一つは一様な背景電荷 `[bz] zbak` を置く方法です。
全体はいつも電気的に中性で、価電子の数は式 (1) で決まります。

$$ N_{\rm val} = \sum_{\rm atoms} Z - N_{\rm core} - z_{\rm bak} \tag{1}$$

Si（単位胞に 2 原子、$N_{\rm core} = 20$）で、表 1 の四つの場合を計算します。
場合 D は、場合 C の正孔を固定スピンモーメント `[bz] fsmom` でスピン分極させたものです。

## 実行

```bash
# 場合 A: Z = 14.2（ctrlg.si.toml のまま）
lmfa si > llmfa
mpirun -np 4 lmf si > llmf

# 場合 B: Z = 14.2 と背景電荷 0.3
rm -f rst.si __mixm.si
lmfa si --ctrlg:bz.zbak=0.3 > llmfa
mpirun -np 4 lmf si --ctrlg:bz.zbak=0.3 > llmf

# 場合 C: Z = 14 と背景電荷 0.2
rm -f rst.si __mixm.si
OPT="--ctrlg:spec.1.z=14 --ctrlg:spec.1.q=[2.0,2.0] --ctrlg:bz.zbak=0.2"
lmfa si $OPT > llmfa
mpirun -np 4 lmf si $OPT > llmf

# 場合 D: 場合 C をスピン分極させ、全スピンモーメントを 0.2 に固定
rm -f rst.si __mixm.si
lmfa si $OPT --ctrlg:ham.nspin=2 --ctrlg:bz.fsmom=0.2 > llmfa
mpirun -np 4 lmf si $OPT --ctrlg:ham.nspin=2 --ctrlg:bz.fsmom=0.2 > llmf
```

## 見るところ

電荷は次で確かめます。分数の Z が使われていれば `Z= 14.20`、`z= 14.200`、`nucleii -28.40000` と出ます。

```bash
grep 'conf: Species' llmfa        # 原子の計算: Z, 内殻 Qc, 原子の電荷 Q（0 なら中性）
grep ' z=' llmf | head -2         # lmf が使う核電荷
grep -A1 'Charges:' llmf | tail -2
```

場合 A の出力:

```
conf: Species  Si       Z=  14.20 Qc=  10.000 R=  2.160000 Q=  0.000000 nsp= 1 mom=  0.000000
 site  1  z= 14.200  rmt= 2.16000  nr=555   a=0.015  nlml=25  rg=0.540  Vfloat=T
   Charges:  valence     8.40000   cores    20.00000   nucleii   -28.40000
   hom background     0.00000   deviation from neutrality:      0.00000
```

表 1: 四つの場合（2026-09-30、gfortran、4 並列、k メッシュ 6×6×6）。$E_{\rm F} - E_{\Gamma}$ は、Γ 点の三重に縮退した価電子帯の頂上から測った
フェルミエネルギー（`log.si` の最後の `fp evl` と `bzmet` の行から）。

| 場合 | Z | zbak | 価電子数 | 全エネルギー ehf (eV) | $E_{\rm F}$ (Ry) | $E_{\rm F} - E_{\Gamma}$ (eV) |
| --- | --- | --- | --- | --- | --- | --- |
| A | 14.2 | 0 | 8.4 | −16274.921209 | 0.319350 | +1.91 |
| B | 14.2 | 0.3 | 8.1 | −16277.374785 | 0.293646 | +1.49 |
| C | 14 | 0.2 | 7.8 | −15732.101789 | 0.163924 | −0.88 |
| D | 14 | 0.2 | 7.8 | −15732.054415 | 0.189198 | （スピン分極） |

A と B は電子を加えた場合で、フェルミエネルギーは伝導帯の中にあります。C と D は電子を減らした場合で、価電子帯の中にあります。
場合 D では、スピンモーメントが 0.2001 $\mu_{\rm B}$ になるように、スピンごとに逆向きのポテンシャルのずれが加わります。
その大きさ $-V_\uparrow + V_\downarrow$ = 0.082841 Ry は `efermi.lmf` の一行目の二つ目の数です。

## 注意

- 分数の Z を使うときは、原子が中性になるように価電子の配置 `[[spec]] q`（s, p, d, f の順）も合わせます。
  Z = 14.2 なら `q = [2.0, 2.2]`。内殻は Z の整数部分の原子のものです。
- `lmfa` にも `lmf` と同じ `--ctrlg:` を付けます。場合を変えるときは `rst.si` と `__mixm.si` を消します。
- `--ctrlg:<節.キー>=値` が書き換えるのは、`ctrlg.si.toml` に書いてあるキーだけです。無いキーは読み飛ばされ、
  出力の先頭に `(note) --ctrlg: key "..." not in target section; override skipped` と出ます。
  そのため `zbak` と `fsmom` を使わないときの値で（`zbak = 0.0`、`fsmom = -99999.0`）、`mmom` を場合 D の初期値で書いてあります。
- 古い入力を変換した `ctrlg.*.toml` では、Z が整数に丸められていることがあります。分数の Z は手で書きます。

## 時間

`testecalj Si -np 4`（四つの場合を続けて計算）は 24〜28 秒（16 コアの計算機、4 並列）。

## テスト

`testecalj Si -np 4`（一つ上のディレクトリで）。場合ごとに、`save.si.<場合>` の全エネルギーとスピンモーメント、
`log.si.<場合>` の Z・価電子数・背景電荷・フェルミエネルギーを比べます。

もとは `Samples/Legacy/Si_doping_sample` です。説明は ecaljdoc の `manual/Memo_bgcharge.md` にもあります。
