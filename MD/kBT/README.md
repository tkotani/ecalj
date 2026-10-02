# MD/kBT — 有限温度、MLO-QSGW、GW の GPU 高速化の開発の記録（案内）

ここにあるもの（開発の記録。公開しない）と、サンプルのそばにある記録へのリンク。2026-10-02 22:03。

| 文書 | 中身 |
| --- | --- |
| [sigma_mlo_design.md](sigma_mlo_design.md) | MLO-QSGW（自己エネルギーを MLO 表現で内挿する QSGW）の設計書 |
| [gpu_fp32_plan.md](gpu_fp32_plan.md)、[gpu_fp32_report.md](gpu_fp32_report.md) | GW（hgw）の GPU 高速化の計画と報告（2026-09-27） |
| [finiteT_202606.md](finiteT_202606.md) | 2026-06 の有限温度 QSGW の計算の記録（当時のキー） |
| [LiTi2O4_finiteT_202606.md](LiTi2O4_finiteT_202606.md) | LiTi₂O₄ の 2026-06 の有限温度の記録（当時のキー。一部の解釈は後で覆った） |

サンプルのそばにある記録:

| 文書 | 中身 |
| --- | --- |
| [Samples/kBT/LiTi2O4/README.md](../../Samples/kBT/LiTi2O4/README.md) | LiTi₂O₄ の MLO-QSGW の 40 反復のまとめ（2026-09-30。tf32・fp32・fp64、表 1〜5、図 1・2） |
| [Samples/kBT/LiTi2O4/six_patterns.md](../../Samples/kBT/LiTi2O4/six_patterns.md) | {6³, 9³} × {tf32, fp32, fp64} のすべての反復の図と表（図 1〜7、表 5・6） |
| [Samples/kBT/LiTi2O4/input/qmlo/README.md](../../Samples/kBT/LiTi2O4/input/qmlo/README.md) | 本番の入力、環境の作り方、反復ごとの正常値 |
| [Samples/kBT/README.md](../../Samples/kBT/README.md) | Samples/kBT の案内（Si・GaAs・Fe・Cu の温度のスキャンなど） |

キーと方法の説明は ecaljdoc の [kBT](../../ecaljdoc/manual/kBT.md)・[mlo_gwsc](../../ecaljdoc/manual/mlo_gwsc.md)・[ecaljgpu](../../ecaljdoc/manual/ecaljgpu.md)。
