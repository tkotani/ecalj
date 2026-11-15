# サンプル — `ecalj/Samples/` のどこに何があるか

このページは、[`ecalj/Samples/`](https://github.com/tkotani/ecalj/tree/main/Samples) の木の簡単な索引です。
各行は、そのディレクトリの説明の正本である README を指します。役割と「何から始めるか」の列は、すべての枝を開かなくても入口を選べるようにするためのものです。

## 最上位の構成（2026-05 から）

```
ecalj/Samples/
├── README.md          (https://github.com/tkotani/ecalj/blob/main/Samples/README.md)
├── GetStarted/        minimal seeds for the README_tutorial walk-through
├── TestInstall/       testecalj-runnable (install validation)
├── EPS/               testecalj-runnable
├── PROCAR/            testecalj-runnable
├── MLOsamples/        testecalj-runnable
├── MLOQSGW/           testecalj-runnable (QSGW with gwsc --mlo)
├── AFsymmetry/        testecalj-runnable (antiferromagnetic symmetry, LDA and QSGW)
├── AFfixMMOM/         testecalj-runnable (AF moment held at a given value, mmtarget.aftest)
├── Magnon/            testecalj-runnable (magnon spectra, job_mlo_magnon)
├── kBT/               electron temperature in GW: scanT/ (testecalj-runnable), LiTi2O4, research log
├── FermiSurface/ Doping/ HomoGas/ BoltzTraP/ SLAB/ SOC/ LDAU/ EffectiveMass/ Relax/ DOS/ IIR/
│                      short examples rebuilt from Legacy/, each testecalj-runnable
├── MATERIALS/         LDA and MLO samples of 65 materials; inputs and results of QSGW of larger systems: not tests
├── BenchmarkTest/     testecalj-runnable (QSGW of InAs/GaSb superlattices, for GPU)
├── mptf32problem/     GW1500 failure corpus: TF32 vs FP32 on AgNO3 etc.
└── AtomDimer/         a molecule and an atom in a box (N2): bond length and binding energy
```

どのディレクトリの入力も `ctrlg.<sname>.toml` で、`lmf`、`gwsc`、`job_eps` などがそのまま読みます。
試験のターゲットの組は、ecalj の `TOOLS/samples_tests.sh` が順に回します（[試験の仕組み](./dev_testecalj)）。
各ディレクトリのターゲットと、それを回す `samples_tests.sh` の組は、ecalj の
[`Samples/README.md`](https://github.com/tkotani/ecalj/blob/main/Samples/README.md) の表 1 にあります。
古い `ctrl.<sname>` + `GWinput` の自分のディレクトリは、[入力ファイル](./toml_migration) の手順で変換してください。

## 試験になるサンプル（testecalj で回せる）

| ディレクトリ | 物理・役割 | ディレクトリの README | 説明のページ |
|---|---|---|---|
| **GetStarted/GaAs** | [チュートリアル](./README_tutorial#getstarted) の最小の例。閃亜鉛鉱の GaAs、2 原子。`ctrls.gaas` と `ctrlg.gaas.toml` を含む | [GetStarted/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/GetStarted/README.md) / [GaAs README](https://github.com/tkotani/ecalj/blob/main/Samples/GetStarted/GaAs/README.md) | [メインドキュメント](./README_tutorial) |
| **EPS/EPS_Cu** | 誘電関数 ε(q,ω)、`job_eps`、fcc の金属 | [EPS_Cu/](https://github.com/tkotani/ecalj/tree/main/Samples/EPS/EPS_Cu)（README なし。test.py を参照） | [誘電関数](./optical) |
| **EPS/EPS_GaAs** | ε(q,ω)、閃亜鉛鉱の半導体 GaAs | [EPS_GaAs/](https://github.com/tkotani/ecalj/tree/main/Samples/EPS/EPS_GaAs)（README なし。test.py を参照） | [誘電関数](./optical) |
| **EPS/EPS_Ag** | ε(q,ω)、Ag（fcc、局所場補正なし） | [EPS_Ag/](https://github.com/tkotani/ecalj/tree/main/Samples/EPS/EPS_Ag)（README なし。test.py を参照） | [誘電関数](./optical) |
| **MLOsamples/**（試験のターゲット 25） | MTO の局在軌道（Wannier の代わり）。SOC の例、`job_mloW` によるオンサイトの W | [MLOsamples/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/MLOsamples/README.md) | [MLO](./mlo) |
| **PROCAR/MgO_PROCAR** | ファットバンドの重み。MgO のバンドへの O 2p の射影 | [MgO_PROCAR/](https://github.com/tkotani/ecalj/tree/main/Samples/PROCAR/MgO_PROCAR)（README なし。test.py を参照） | [その他の説明 § PROCAR](./UsageDetailed#procar-mode) |
| **PROCAR/Ni2MnGa_L21_PROCAR** | 原子ごとのファットバンド。強磁性の Heusler。収束した `rst.ni2mnga` を含む | [Ni2MnGa_L21_PROCAR/](https://github.com/tkotani/ecalj/tree/main/Samples/PROCAR/Ni2MnGa_L21_PROCAR)（README なし。test.py を参照） | [その他の説明 § PROCAR](./UsageDetailed#procar-mode) |
| **TestInstall/**（`testecalj --all` で 26 ターゲット、66 件の確認） | インストールの確認: 基底状態、GW（`gwsc`、`gw_lmfh`）、有限温度（`fe_kbt`）、eps（`job_eps`）、ChiPM、MLO のモデルの cRPA。重い GW のターゲット `cugase2_gwsc222`、`nio_gwsc444`、`pdo_gwsc443`、`gas_gwsc666` は一つずつ回す | （ディレクトリごと。`testecalj --all` で回す） | [gwsc](./gwsc)、[誘電関数](./optical) |
| **MLOQSGW/**（`GaAs`、`NiO`） | 自己エネルギーを MLO の表現で内挿する QSGW（`gwsc --mlo`） | [MLOQSGW/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/MLOQSGW/README.md) | [MLO-gwsc](./mlo_gwsc) |
| **BenchmarkTest/**（`inas2gasb2`、`inas4gasb4`） | InAs/GaSb の超格子の QSGW の 1 反復を、保存した `QPU.1run` と比べる。32 GB の GPU 1 枚では `-np2 1`。ISSP の計算機用のジョブのスクリプトは `BenchmarkTest/ISSP/`（[ISSP](../install/installISSP)） | （README なし。test.py を参照） | [GPU 版](./ecaljgpu) |
| **AFsymmetry/**（`NiO`、`NiSe`、`NiO_gwsc`） | サイトに `af = 1 / -1` を付けた反強磁性の対称性 `symgrpaf`。LDA と QSGW（GW のプログラムは両方のスピンを計算する） | [AFsymmetry/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/AFsymmetry/README.md) | [lmf](./lmf#symgrp) |
| **AFfixMMOM/**（`NiO_afsym`、`NiO_noafsym`） | 反強磁性の対のモーメントを決まった値に保つ（`mmtarget.aftest`）。NiO を `symgrpaf` あり・なしで | [AFfixMMOM/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/AFfixMMOM/README.md) | [その他の説明](./UsageDetailed#antiferro-symmetry-without-soc) |
| **Magnon/**（`Fe_mlo_magnon`） | MLO のモデルによるマグノンのスペクトル（`job_mlo_magnon`、bcc Fe）。2026-10-02 より前の版（タグ `last-wannier`）の Wannier 関数と比べる | [Magnon/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/Magnon/README.md) | [その他の説明](./UsageDetailed) |
| **kBT/scanT/**（`Si`） | 同じ入力を、いくつかの電子温度で（`scan_T.sh`、`collect_T.py`）: Si のギャップ、bcc Fe のバンドとモーメント。試験のターゲットは `Si` だけ（`Cu`、`Fe`、`GaAs` のディレクトリもある） | [kBT/scanT/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/kBT/scanT/README.md) | [kBT](./kBT) §8 |
| **FermiSurface/Cu**、**Doping/Si**、**HomoGas/es**、**BoltzTraP/Si**、**SLAB/Cu001**、**SOC/**（`FePt_MAE`、`MnGa_MAE`）、**LDAU/ReN**、**EffectiveMass/**（`GaAs`、`CdS`、`GaN`）、**Relax/LaGaO3**、**DOS/**（`ZnS`、`Fe`）、**IIR/C**、**AtomDimer/N2** | 短い例（どれも数分以内）: xcrysden のためのフェルミ面、半端な Z か背景電荷によるドーピング、Lindhard 関数と解析的な形との比較、BoltzTraP の入力、仕事関数を見るスラブ、磁気異方性エネルギー、LDA+U、SOC を入れた有効質量、原子位置の緩和、全・部分状態密度、衝突電離の割合（自己エネルギーの虚部）、箱の中の分子 | 各ディレクトリの README.md | [その他の説明](./UsageDetailed) |

## 試験にしないもの

| ディレクトリ | 物理・役割 | ディレクトリの README | 説明のページ |
|---|---|---|---|
| **MATERIALS/** | 65 物質の LDA と MLO のサンプル（簡単な 62 物質、La2CuO4、InAs/GaSb、BaTiO3。MLO のモデルの選び方は [MLO](./mlo) §9）。大きな系の QSGW: La2CuO4（20 反復、7 原子）、16 原子と 40 原子の InAs/GaSb の超格子、BaTiO3（三段階）の入力と結果を残してある。時間のかかる計算 | [MATERIALS/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/MATERIALS/README.md) | [gwsc](./gwsc) |
| **kBT/LiTi2O4** と kBT/ の残り | LiTi₂O₄（金属のスピネル）: MLO-QSGW の計算の入力（`LiTi2O4/input/qmlo`: `t_tetrakbt = -992.4`、`t_sigmaw = 1000`）とスクリプトと図、有限温度の QSGW の計算（`t_tetrakbt` > 0）、bcc Fe での対照実験、GPU の精度の報告、研究の記録。**入力と結果だけ**で、回帰試験には重すぎる（有限温度の試験のターゲットは `TestInstall/fe_kbt` と `kBT/scanT/Si`） | [kBT/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/kBT/README.md) | [kBT](./kBT) |

### `Samples/MLOsamples/` の早見（testecalj の 25 件）

[MLOsamples/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/MLOsamples/README.md) の表に、各例の役割と参照のデータがあります。主なもの:

- `GaAsSoc`、`FeSoc`、`FeMgOSoc` — 摂動としての SOC の試験（`job_mlo_soc`。それが描く SOC 入りの DFT のバンドと比べる、[MLO](./mlo) §4）
- `C`、`C.sp`、`Cu`、`SrTiO3`、`GaAs`、`Si666gwsc`（QSGW） — 非磁性の半導体・金属
- `Ag`、`Al`、`NaCl`、`SiC`、`CdTe`、`ZnO`、`TiO2` — Materials Project の構造で、既定の設定
- `Fe`、`FeCo`、`FeMgO` — 磁性金属（SOC なし）
- `Al2O3_Cr`、`GdCo5`、`GdION`、`SmP`、`NiO666lda`、`RuO2` — LDA+U・4f・反強磁性
- `Fe` の test.py は `job_mloW` を回し、オンサイトの V、W−V を中に書いた参照値と比べる。
  [MLOsamples/README.md § Regression test (testecalj Fe)](https://github.com/tkotani/ecalj/blob/main/Samples/MLOsamples/README.md#regression-test-testecalj-fe) を参照。

### `Samples/TestInstall/` が何を試すか

ディレクトリの中で `testecalj --all -np 8` として回します。ディレクトリごとの役割:

| ディレクトリ | 試す道具 |
|---|---|
| `c`、`te`、`zrt`、`co`、`cr3si6`、`felz`、`gasls`、`eras`、`crn`、`cu`、`na`、`fe`、`gdn` | LDA / lmf のインストールの確認（`fe`: スピン分極の状態密度、`gdn`: LDA+U） |
| `copt` | LDA + 力の計算 |
| `gas_eps_lmfh`、`gas_epsPP_lmfh` | GaAs の誘電関数 ε(q,ω)（局所場補正あり・なし） |
| `fe_epsPP_lmfh_chipm` | bcc Fe のスピン感受率 χ⁻⁻ |
| `si_gw_lmfh`、`gas_pw_gw_lmfh`、`si_gwsc`、`gas_gwsc`、`nio_gwsc`、`fe_gwsc` | one-shot の GW と、全体の QSGW（`gwsc`） |
| `fe_kbt` | bcc Fe の、有限温度の四面体法（`t_tetrakbt` > 0）による QSGW |
| `ni_crpa`、`srvo3_crpa` | constrained RPA |

### EPS のサンプル — 三つの違い

詳しい説明は [誘電関数](./optical) にあります。見分けるための特徴:

| ディレクトリ | sname | プログラム | BZ のメッシュ | QforEPS | 特徴 |
|---|---|---|---|---|---|
| `EPS_Cu` | `cu` | `job_eps` | 12³（lmf）/ 10³（gw） | 小さな q を 3 つ（5e-4 / 1e-3 / 2e-3 a.u.） | fcc の金属。バンドを描くスクリプトを含む |
| `EPS_GaAs` | `gaas` | `job_eps` | 4³ / 8³ | Cu と同じ 3 つの q | 閃亜鉛鉱の半導体 |
| `EPS_Ag` | `ag` | `job_eps` | 8³ / 8³ | Cu と同じ | fcc の金属、局所場補正なし。参照のデータは 2026-05-08 に作成 |

三つとも、参照の EPS のファイル（`EPS000{1,2,3}.nlfc.dat.{inter,intra}bandonly`）をコミットしてあり、test.py は、バンド内の 1/(ω+iη) の発散を skipcond で除いて比べます。
`testecalj EPS_Cu / EPS_GaAs / EPS_Ag -np 4` は、どれも 6/6 PASSED になります。

## Legacy/ — なくした

`Samples/Legacy/` の古い例は、上のディレクトリとして作り直し（2026-09-29、30）、ディレクトリは 2026-09-30 にリポジトリから外しました（ecalj のコミット `932c59a6d` より前の履歴にあります）。
最後の `TestHomoDimerAtom`（箱の中の二量体と原子、2012）は `AtomDimer/N2` になりました。

## サンプルに関わる説明の置き場

重複を減らすため、話題ごとに一つのページが正本で、この Samples のページは指すだけです。

| 話題 | 正本のページ |
|---|---|
| 入力ファイル、`Legacy2toml.py` の手順、古いキーと TOML のパスの対応 | [入力ファイル](./toml_migration) |
| `ctrlg.<sname>.toml` の構成 | [lmf](./lmf#file-structure-sections) |
| `[gw]` のキー（以前の `GWinput`）、積基底 | [GW の設定](./gwinput) |
| 有限温度、キー `t_tetrakbt` / `t_sigmaw` | [kBT](./kBT) |
| MLO-QSGW（`gwsc --mlo`） | [MLO-gwsc](./mlo_gwsc) |
| `gwsc` による QSGW の全体 | [gwsc](./gwsc) |
| 誘電関数 ε(q,ω)（`job_eps`） | [誘電関数](./optical) |
| ファットバンドの重み・`PROCAR`（`MgO_PROCAR`、`Ni2MnGa_L21_PROCAR`） | [その他の説明 § PROCAR](./UsageDetailed#procar-mode) |
| MLO のバンドとオンサイトの W（`job_mloW`） | [MLOsamples/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/MLOsamples/README.md) |
| `ecalj_auto` のバッチ | [ecalj_auto/README.md](https://github.com/tkotani/ecalj/blob/main/ecalj_auto/README.md) |
| POSCAR から QSGW までのチュートリアル | [メインドキュメント](./README_tutorial) |

新しい Samples のディレクトリを足すときの約束（`testecalj`、`comp.py`、名前の付け方）は [試験の仕組み](./dev_testecalj) にあります。
