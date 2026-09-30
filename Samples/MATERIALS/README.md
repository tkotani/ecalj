# Samples/MATERIALS — 大きめの系の QSGW の入力と結果

試験（`testecalj`）にはしない、時間のかかる計算の入力と、以前の計算の結果。`Samples/Legacy` から 2026-09-30 に移した。
どれも `ctrlg.<sname>.toml` が入力で、旧形式（`ctrl.<sname>`、`GWinput`）は読まれない。

**表 1**. ディレクトリ

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

## 旧形式の古いサンプル（2026-10-01 に最上位の `MATERIALS/` から移した）

最上位にあった `MATERIALS/` を、いったんここへ移した（user 2026-10-01）。旧形式（`ctrl.<sname>`、`GWinput`）の入力と結果で、今のプログラムはそのままでは読まない。
今の入力にするには `Legacy2toml.py <sname>`（`ctrl.<sname>` と `GWinput` から `ctrlg.<sname>.toml`）。各ディレクトリの中身と、近い内容の今のサンプルは
ecalj の [past_log.md](../../past_log.md) §9。整理（今の形にするか、trash に移すか）は ecalj の [TODOandQuestion.md](../../TODOandQuestion.md)。

`CTRLsample`、`EPS_Ag`、`FeMTOHAM`、`Fe_5`、`Fe_HamMTO`、`GdNldau`、`H_atom`、`LaGaO3_relax`、`LiDOS_Discrete`、`Li_atom`、`MLOsamples`、`NiMnSb_magnon`、`NiO`、
`NiSe_aftest`、`SYMLsamples`、`SiSigma`、`SiSigmaAny`、`Si_HamMTO`、`Si_doping_sample`、`TESTsamples`、`TETRAHEDRON_HomoGas`、`TestHomoDimerAtom`、
`cugase2_gwsc222`、`erasldau`、`mass_fit_test`、`pdo_gwsc443`、`yh3fcc_gwsc666`、`Materials.ctrls.database`・`job_materials.py`・`job_check`・`job_vbm`・`sortvbm.py`
