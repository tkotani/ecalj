# Samples/kBT — 電子温度と均し方のサンプル、MLO-QSGW の結果、研究の記録

`[gw]` の `t_tetrakbt`（χ0 側）と `t_sigmaw`（Σ 側）の使い方と、それで何が変わるかを見るためのサンプル。
キーの意味と式は ecaljdoc の [kBT](https://ecalj.github.io/ecaljdoc/manual/kBT)（表 1 と §2・§3）。

**表 1**. `t_tetrakbt` と `t_sigmaw`（単位は K）

| キー | 値 | 意味 |
| --- | --- | --- |
| `t_tetrakbt`（必須） | 0 | χ0 は T = 0 のテトラヘドロン法。均さない |
| | T > 0 | χ0 は温度 T の有限温度のテトラヘドロン法。Fermi 準位は `EFERMI_kbt` |
| | −T | χ0 は T = 0 のテトラヘドロン法で、Im χ0 を Gaussian で均す。幅は温度 T の Fermi-Dirac と同じ |
| `t_sigmaw` | T | Σ の中間準位を均す Fermi-Dirac 核の温度（既定 1000） |

入力を作る `ctrlgenToml.py` は `t_tetrakbt = 300`、`t_sigmaw = 300` を書く。

## サンプル

**表 2**. このディレクトリのサンプルと記録

| ディレクトリ・ファイル | 中身 | 計算の時間 |
| --- | --- | --- |
| [`scanT/`](scanT/README.md) | **温度のスキャン**。同じ入力を温度だけ変えて QSGW で回す。半導体 Si・GaAs（ギャップが電子温度でどう変わるか）と、金属 bcc Fe・fcc Cu（バンドと磁気モーメント） | 1 つの設定が 2〜10 分 |
| [`LiTi2O4/`](LiTi2O4/README.md) | **金属スピネル LiTi₂O₄ の MLO-QSGW**（6³・9³、40 反復、tf32・fp32・fp64 の比較）のまとめと、入力、投入・描画のスクリプト | 9³ の 1 反復が 21 分（RTX 5090 × 2） |
| [`Fe/`](Fe/README.md) | bcc Fe 3000 K で Σ 側の温度だけを変えた対照の計算（2026-06、結果は 2026-09-27 に作り直し） | 40 秒 |
| `NiO/`、`Si/` | MLO-QSGW の図と、それを描くスクリプト | — |
| [`contour_test/`](contour_test/README.md) | Σc の積分路の試験 | — |
| [`bench_hgw/`](bench_hgw/README.md) | `hgw` の時間を測るスクリプト | — |

`testecalj` の試験は `scanT/Si`（Si 3000 K、1 反復）と `Samples/TestInstall/fe_kbt`（Fe 3000 K）。

## 記録

| ファイル | 中身 |
| --- | --- |
| [`MD/research_log.md`](../../MD/research_log.md) | 研究ログ（2026-09-16〜。有限温度、MLO-QSGW、GPU 高速化）。新しいものが上、時刻付き |
| [`gpu_fp32_report.md`](gpu_fp32_report.md)、[`gpu_fp32_plan.md`](gpu_fp32_plan.md) | GW の GPU 高速化の報告と計画（精度、行列積の方法の表、QSGW 1 反復の内訳） |
| [`sigma_mlo_design.md`](sigma_mlo_design.md) | MLO-QSGW（Σ を MLO で持って内挿する）の設計 |
| [`README_202606_finiteT.md`](README_202606_finiteT.md) | 2026-06 の有限温度 QSGW の計算の記録（キーの名前は当時のもの） |
| [`fp16acc.cu`](fp16acc.cu) | テンソルコアの和の精度を調べる試験 |

研究のテーマごとの要約は ecaljdoc の ForDevelopers §13 と ForDevelopers_research。
