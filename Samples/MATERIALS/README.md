# Samples/MATERIALS — LDA の計算と MLO の計算のサンプル（65 物質）と、大きめの系の QSGW の入力と結果

**LDA の計算と MLO の計算のサンプル**。62 物質（1 節）と La₂CuO₄・InAs/GaSb 超格子・BaTiO₃（3 節）の入力があり、その 65 物質で
LDA（Eu の化合物は LDA+U）と MLO の模型を回して、DFT のバンドと比べた（2 節）。MLO の模型の選び方（基準 1・2・3）、誤差の測り方、
結果は ecaljdoc の [mlo §9](https://ecalj.github.io/ecaljdoc/manual/mlo)。

試験（`testecalj`）にはしない（`TOOLS/samples_tests.sh` の `inputs` の組が、ctrlg が読めることだけ確かめる）。
どれも `ctrlg.<sname>.toml` が入力で、旧形式（`ctrl.<sname>`、`GWinput`）は読まれない。
2026-10-01 まで 62 物質は `Database/` の下にあったが、1 段にした（user「2 段にする必要はない」）。

- 1 節: 62 物質の入力（物質ごとのディレクトリ、表 1）
- 2 節: LDA と MLO の計算（2026-10-01、1 物質の回し方、`mlocheck/`）
- 3 節: 大きめの系の QSGW の入力と結果（`La2CuO4/`、`InAsGaSb/`、`BaTiO3/`、表 2）

## 1. 62 物質の入力（`Materials.ctrls.database` を展開したもの、2026-10-01）

以前は `job_materials.py` が `Materials.ctrls.database`（構造の雛形 17 種と物質ごとの格子定数・指定）から、その場でディレクトリと入力を作って
LDA/GGA を回していた。これを物質ごとのディレクトリに展開した（user「展開してしまうのが便利かも」）。元の `job_materials.py` とデータベースは
`trash/`（中身は ecalj の `MD/past_log.md` §9）。

各ディレクトリ `<物質>/` には 2 つだけ置く:

- `ctrls.<sname>`: 構造だけ（`job_materials.py --noexec` がデータベースの雛形から作ったもの。`%const` の式と `{a}` の置き換えを含む）
- `ctrlg.<sname>.toml`: `ctrlgenToml.py <sname>` で作った入力。データベースの指定を次のように写した: `--nk1..3` → LDA の k 点（`[bz] nkabc`）、
  `--nspin=2`・`lmf-vnspin=2` → `--nspin=2`、`lmf-vso=1` → `--so=1`、`mkGW-a,b,c` → `[gw] n1n2n3`、`ctrls` の `MMOM` → `[[spec]] mmom`
  （`ctrlgenToml.py` は `ctrls` の `MMOM` を写さないので後から書き足した。原子の表に既定の `MMOM` がある Eu では、その既定より `ctrls` の値を優先する。
  2026-10-01 に直すまで EuS・EuSe・EuTe の `Nidn` が +6 のままだった）。GW の節（`[gw]` `[mlo]` `[blocks]` `[product_basis]`）も入れてあり、
  QSGW と MLO の自動の模型（`mlo_method = 4`、`mlo_lm`、`mlo_nkabc`）がそのまま回せる

**LDA と MLO の模型は 2026-10-01 に全物質で回した**（下の「MLO の試験」）。GW は回していない。回すときは、物質のディレクトリを写してその中で:

```bash
lmfa <sname>; mpirun -np 4 lmf <sname> > llmf         # LDA/GGA
getsyml <sname>; job_band <sname> -np 4               # バンド（syml.<sname> を作ってから）
gwsc 5 <sname> -np 8                                  # QSGW（ecaljdoc manual/gwsc）
job_mlo <sname> -np 4                                 # MLO の自動の模型（lmf と job_band の後。ecaljdoc manual/mlo）
```

注意:

- Eu の化合物は、`ctrlgenToml.py` の原子の表の既定で 4f に `idu = 12`・`uh`・`jh` が入る（`sigm` が無いときは FLL の LDA+U、QSGW では U を切る）。
  Ce は nspin = 1 なので、その LDA+U はコメントにしてある（LDA+U は nspin = 2 が要る。2026-10-01）
- AF II の雛形（NiO、MnO、EuS、EuSe、EuTe）は、2 つの磁性の位置を種 `Niup`・`Nidn` として持つ（元素は `z` で決まる。名前は NiO の雛形のまま）
- `job_materials.py --all` は LaGaO3・4hSiC・Bi2Te3 を重いので外していた。Bi2Te3 はデータベースの注に「混合をゆっくり（beta 0.3）」とある
- 格子定数は、`Å` と書いたものはデータベースで Å から換算したもの、`bohr` はデータベースが bohr で与えたもの

**表 1**. 62 物質（ディレクトリ名、`sname`、構造、格子定数、LDA の k 点、GW の `n1n2n3`、スピン・SO・初期の磁気モーメント、データベースの注）

| 物質 | sname | 構造 | 格子定数 | k 点 | n1n2n3 | スピンなど | 注 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2hSiC | `2hsic` | wurtzite | a=3.076 Å c=5.048 Å | 6×6×4 | 4×4×2 |  | 2H-SiC COD 9008875.cif |
| 3cSiC | `3csic` | zinc blende | a=4.348 Å | 8×8×8 | 6×6×6 |  | ZincBlend COD 9008856.cif 3C-SiC |
| 4hSiC | `4hsic` | 4H-SiC |  | 6×6×2 | 4×4×2 |  |  |
| AlAs | `alas` | zinc blende | a=5.6611 Å | 8×8×8 | 6×6×6 |  |  |
| AlN | `aln` | wurtzite | a=3.112 Å c=4.982 Å | 6×6×4 | 4×4×2 |  |  |
| AlNzb | `alnzb` | zinc blende | a=4.38 Å | 8×8×8 | 6×6×6 |  |  |
| AlP | `alp` | zinc blende | a=5.4672 Å | 8×8×8 | 6×6×6 |  |  |
| AlSb | `alsb` | zinc blende | a=6.1355 Å | 8×8×8 | 6×6×6 |  |  |
| Bi2Te3 | `bi2te3` | Bi2Te3 (rhombohedral) | a=4.7825489 coa=4.0154392 bohr | 5×5×5 | 3×3×3 | nspin 2, so 1 | slower mixing needed (the database note: beta 0.3) |
| C | `c` | diamond | a= | 8×8×8 | 6×6×6 |  |  |
| CdO | `cdo` | rock salt | a=8.92 bohr | 8×8×8 | 6×6×6 |  |  |
| CdS | `cds` | zinc blende | a=11.01 bohr | 8×8×8 | 6×6×6 |  |  |
| CdSe | `cdse` | zinc blende | a=11.44 bohr | 8×8×8 | 6×6×6 |  |  |
| CdTe | `cdte` | zinc blende | a=12.25 bohr | 8×8×8 | 6×6×6 |  |  |
| Ce | `ce` | fcc | a=9.052 bohr | 12×12×12 | 12×12×12 | mmom (1 spec) |  |
| Cu | `cu` | fcc | a=6.798 bohr | 12×12×12 | 12×12×12 |  |  |
| EuO | `euo` | rock salt | a=5.14 Å | 8×8×8 | 6×6×6 | nspin 2, mmom (1 spec) |  |
| EuS | `eus` | rock salt, AF II | a=5.968 Å | 6×6×6 | 4×4×4 | nspin 2, mmom (2 spec) | species `Niup`/`Nidn` = the two magnetic sites (template of NiO) |
| EuSe | `euse` | rock salt, AF II | a=6.195 Å | 6×6×6 | 4×4×4 | nspin 2, mmom (2 spec) | species `Niup`/`Nidn` = the two magnetic sites (template of NiO) |
| EuTe | `eute` | rock salt, AF II | a=6.598 Å | 6×6×6 | 4×4×4 | nspin 2, mmom (2 spec) | species `Niup`/`Nidn` = the two magnetic sites (template of NiO) |
| Fe | `fe` | bcc | a=2.8665 Å | 12×12×12 | 12×12×12 | nspin 2, mmom (1 spec) |  |
| GaAs | `gaas` | zinc blende | a=5.65325 Å | 8×8×8 | 6×6×6 |  |  |
| GaAs_so | `gaas_so` | zinc blende | a=5.65325 Å | 8×8×8 | 6×6×6 | nspin 2, so 1 |  |
| GaN | `gan` | wurtzite | a=3.189 Å c=5.189 Å | 6×6×4 | 4×4×2 |  |  |
| GaNzb | `ganzb` | zinc blende | a=4.50 Å | 8×8×8 | 6×6×6 |  |  |
| GaP | `gap` | zinc blende | a=5.4505 Å | 8×8×8 | 6×6×6 |  |  |
| GaSb | `gasb` | zinc blende | a=6.0959 Å | 8×8×8 | 6×6×6 |  |  |
| Ge | `ge` | diamond | a= | 8×8×8 | 6×6×6 |  |  |
| HfO2 | `hfo2` | HfO2 | a=6.700 bohr | 6×6×4 | 4×4×2 |  |  |
| HgO | `hgo` | HgO (COD 9012530) |  | 3×3×6 | 2×2×4 |  |  |
| HgS | `hgs` | zinc blende | a=5.8517 Å | 8×8×8 | 6×6×6 |  | beta-HgS |
| HgSe | `hgse` | zinc blende | a=11.5 bohr | 8×8×8 | 6×6×6 |  |  |
| HgTe | `hgte` | zinc blende | a=12.21 bohr | 8×8×8 | 6×6×6 |  |  |
| InAs | `inas` | zinc blende | a=6.0583 Å | 8×8×8 | 6×6×6 |  |  |
| InN | `inn` | wurtzite | a=3.545 Å c=5.703 Å | 6×6×4 | 4×4×2 |  |  |
| InNzb | `innzb` | zinc blende | a=4.98 Å | 8×8×8 | 6×6×6 |  |  |
| InP | `inp` | zinc blende | a=5.8697 Å | 8×8×8 | 6×6×6 |  |  |
| InSb | `insb` | zinc blende | a=6.4794 Å | 8×8×8 | 6×6×6 |  |  |
| LaGaO3 | `lagao3` | LAGAO3 | a=5.485 Å c=5.523 Å | 3×2×3 | 2×2×2 |  |  |
| Li | `li` | bcc | a=6.35 bohr | 12×12×12 | 12×12×12 |  |  |
| MgO | `mgo` | rock salt | a=7.96 bohr | 8×8×8 | 6×6×6 |  |  |
| MgS | `mgs` | zinc blende | a=10.6198035 bohr | 8×8×8 | 6×6×6 |  |  |
| MgSe | `mgse` | zinc blende | a=11.1678005 bohr | 8×8×8 | 6×6×6 |  |  |
| MgTe | `mgte` | zinc blende | a=12.1315193 bohr | 8×8×8 | 6×6×6 |  |  |
| MnO | `mno` | rock salt, AF II | a=8.40 bohr | 6×6×6 | 4×4×4 | nspin 2, mmom (2 spec) | species `Niup`/`Nidn` = the two magnetic sites (template of NiO) |
| Ni | `ni` | fcc | a=3.524 Å | 12×12×12 | 12×12×12 | nspin 2, mmom (1 spec) |  |
| NiO | `nio` | rock salt, AF II | a=7.88 bohr | 6×6×6 | 4×4×4 | nspin 2, mmom (2 spec) | species `Niup`/`Nidn` = the two magnetic sites (template of NiO) |
| PbS | `pbs` | rock salt | a=5.9362 Å | 8×8×8 | 6×6×6 |  |  |
| PbTe | `pbte` | rock salt | a=6.4620 Å | 8×8×8 | 6×6×6 |  |  |
| Si | `si` | diamond | a= | 8×8×8 | 6×6×6 |  |  |
| SiO2c | `sio2c` | ideal beta cristobalite | a=13.54 bohr | 6×6×6 | 4×4×4 |  | sc-cubic-cristobalite 13.45485515 |
| Sn | `sn` | diamond | a= | 8×8×8 | 6×6×6 |  | gray tin alpha-tin |
| SrTiO3 | `srtio3` | perovskite | a=7.37 bohr | 6×6×6 | 4×4×4 |  |  |
| SrVO3 | `srvo3` | perovskite | a=7.2644 bohr | 8×8×8 | 6×6×6 |  | metal |
| YMn2 | `ymn2` | C15 | a=14.51 bohr | 7×7×7 | 5×5×5 |  |  |
| ZnO | `zno` | wurtzite | a=6.14933837 bohr | 6×6×4 | 4×4×2 |  |  |
| ZnS | `zns` | zinc blende | a=10.23 bohr | 8×8×8 | 6×6×6 |  |  |
| ZnSe | `znse` | zinc blende | a=10.71 bohr | 8×8×8 | 6×6×6 |  |  |
| ZnTe | `znte` | zinc blende | a=11.53 bohr | 8×8×8 | 6×6×6 |  |  |
| ZrO2 | `zro2` | ZrO2 | a=6.726 bohr | 6×6×4 | 4×4×2 |  |  |
| wCdS | `wcds` | wurtzite | a=4.160 Å c=6.756 Å | 6×6×4 | 4×4×2 |  |  |
| wZnS | `wzns` | wurtzite | a=3.82 Å c=6.26 Å | 6×6×4 | 4×4×2 |  |  |

## 2. LDA と MLO の計算（2026-10-01）

全物質（InAs/GaSb の n10 を除く 65）で、LDA のあと MLO の模型を作り、対称線の上で DFT のバンドと比べた。
模型は 3 通り: 基準 1（半内殻の局所軌道を自動で入れる既定）、基準 2（基準 1 + 陽イオンに EH2(s,p) シード）、
基準 3（基準 1 + 空格子球、SiO₂ だけ）。65 物質のうち基準 1 で 53 が good（誤差の最大値 0.02 eV 以下）、物質ごとに一番良い模型で 63 が good、
残る 2 つ（C、Bi₂Te₃）も fair（0.05 eV 以下）。表と図と選び方は ecaljdoc の mlo §9。

1 物質の回し方（GaAs の例。`-np` は手元の数に）:

```bash
cp -r Samples/MATERIALS/GaAs ~/work/ && cd ~/work/GaAs
mpirun -np 1 lmfa gaas
mpirun -np 8 lmf gaas                 # LDA
getsyml gaas --nobzview               # syml.gaas（対称線）
job_band gaas -np 8                   # DFT のバンド
job_mlo gaas -np 8                    # MLO の模型とそのバンド（so = 1 の物質は job_mlo_soc）
mlo_bandcheck.py .                    # 誤差（ecaljdoc mlo の式 (9)〜(11)）
mlo_bandplot.py . -o mlo_GaAs.png     # DFT（灰）と MLO（赤）の図
```

注意:
- ここの ctrlg の `[mlo]` は 2026-10-01 19:3x より前の `gwinit` が書いたもの（H〜Ne は s,p、Na 以降は s,p,d、`mlo_lm2` の行なし、コメントも古い）。
  4f が価電子の原子（Ce、Eu の化合物、HfO₂ の Hf）は f（10〜16）を手で足すか、今の `gwinit` で `[mlo]` を書き直す。
  2026-10-01 の計算では `mlocheck/prep.py` が全物質の `[mlo]` を f 入りの形に書き直し、`rdsig = 0`（DFT のみ）にした
- 基準 2 は `mlo_lm2` に陽イオン（遷移金属・4f・5f を除く。Zn・Cd・Hg は含む）の s,p の行（`<番号> <原子>   1 2 3 4`）を書く。
  今の `gwinit` はこの行を `!` 付きで書くので `!` を外すだけ。ここの古い ctrlg には行が無いので手で書く
- スピン軌道の物質（GaAs_so、Bi₂Te₃）は `job_mlo` でなく `job_mlo_soc`。比べる相手のスピン軌道ありの DFT のバンドも `job_mlo_soc` が描く

ファイル:
- `mlocheck/`: 回したスクリプト（`README.md`）と、結果のページを作る `gallery.py`（図に描いた数値は `page_data/*.npz`）
- `MLOcheck_20261001.tsv`: 最初の版（旧既定）の物質ごとの結果。`MLOcheck_20261001_variants.tsv`: 物質ごとに動径関数を足した変種。
  どちらも 2026-10-01 17:54 までの窓で評価したもの（MLO → DFT が [VBM − 8, CBM + 3]、DFT → MLO が [VBM − 8, CBM + 1] eV）。
  今の `mlo_bandcheck.py` は両方とも [VBM − 8, CBM + mlo_delta]（`mlo_delta` は模型を合わせる範囲の上端、既定 2 eV）
- 経緯は ecalj の `MD/research_log.md`（2026-10-01）

## 3. 大きめの系の QSGW の入力と結果

試験（`testecalj`）にはしない、時間のかかる計算の入力と、以前の計算の結果。`Samples/Legacy` から 2026-09-30 に移した。
どれも `ctrlg.<sname>.toml` が入力で、旧形式（`ctrl.<sname>`、`GWinput`）は読まれない。

**表 2**. 大きめの系のディレクトリ

| ディレクトリ | 系 | 入力 | 残してある結果 | 時間の目安 |
| --- | --- | --- | --- | --- |
| `La2CuO4/` | La₂CuO₄（7 原子、`nkabc` 8³、`n1n2n3` 4³）。QSGW 20 反復（`sigm` なし、`-vssig=1.0` = `scaledsigma` 1.0） | `ctrlg.la2cuo4swj.toml`、`ctrls.la2cuo4swj` | `QPU`、バンド `bnd00*.spin1`、`bandplot.isp1.glt`、ログ `out`（2025-08） | 1 反復 35 分（CPU 32 コア） |
| `InAsGaSb/n10/` | (InAs)₁₀/(GaSb)₁₀ 超格子（40 原子） | `ctrlg.inas10gasb10.toml` | なし | 未計算 |
| `InAsGaSb/n4/` | (InAs)₄/(GaSb)₄（16 原子）。`BenchmarkTest/inas4gasb4` と同じ系で `nkabc` だけ違う（8×8×2） | `ctrlg.inas4gasb4.toml` | なし | `BenchmarkTest` の表 1 |
| `BaTiO3/` | BaTiO₃（5 原子、`nkabc` 4³、`n1n2n3` 3³）。`LDA/`、`QSGW0run/`（1 反復）、`QSGW5run/`（5 反復） | `QSGW5run/ctrlg.batio3.toml`（2026-09-30 に `Legacy2toml.py` で作った）。`LDA/` と `QSGW0run/` は旧形式のまま | 各段のバンド、`rst`、`sigm` | 5 反復で数時間（CPU、2019 の記録） |

回し方（La₂CuO₄）:

```bash
cd La2CuO4
lmfa la2cuo4swj > llmfa
gwsc 10 -np 32 la2cuo4swj --ctrlg:ham.scaledsigma=1.0 > lgwsc      # 旧 jobgw.sh の -vssig=1.0
job_band la2cuo4swj -np 32
```

- `jobgw.sh`（La₂CuO₄）と `job`（BaTiO₃）は当時の投入スクリプトで、旧形式のオプション（`-vssig=1.0`）を含む。上の書き方に読み替える
- 保存してある `rst.*` は当時の版のもので、今の `lmf` では読めないことがある（`EffectiveMass/GaAs` の例）。`sigm` は読める

最上位にあった旧形式の古いサンプル（2026-10-01 にここへ移した 30 項目）は、同じ日に仕分けて trash に移した。中身と拾ったノウハウは ecalj の
[MD/past_log.md](../../MD/past_log.md) §9。
