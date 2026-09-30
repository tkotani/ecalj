# Samples/MATERIALS/Database — 62 物質の入力（`Materials.ctrls.database` を展開したもの、2026-10-01）

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

**LDA と MLO の自動の模型は 2026-10-01 に全物質で回した**（結果は研究ログ `MD/research_log.md` の 2026-10-01）。GW は回していない。回すときは、このディレクトリを写してその中で:

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

**表 1**. 物質（ディレクトリ名、`sname`、構造、格子定数、LDA の k 点、GW の `n1n2n3`、スピン・SO・初期の磁気モーメント、データベースの注）

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
