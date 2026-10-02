# MD/gpu — GW の GPU 高速化のメニュー

2026-10-02 22:15 に MD/kBT/ から分けた。

## いまの要点

- 利用者が選ぶのは精度だけ（`gwsc --prec=tf32|fp32|fp64`）。行列積の方法は機械ごとの表で自動に選ぶ
- 使い方と数値は ecaljdoc [ecaljgpu.md](../../ecaljdoc/manual/ecaljgpu.md)（正本）

## メニュー

| 行き先 | 中身 |
| --- | --- |
| [gpu_fp32_report.md](gpu_fp32_report.md) | 報告（2026-09-27）: 仕組み、段ごとの結果、ベンチマーク、tf32 の FP16 経路、Σc の非同期化、hgw 以外の段 |
| [gpu_fp32_plan.md](gpu_fp32_plan.md) | 計画書（2026-09-27）: 振り分け表、backend、方針ファイル。段 0〜6 は報告で実装済み |
| [Samples/kBT/bench_hgw/README.md](../../Samples/kBT/bench_hgw/README.md) | `hgw` の時間を測るスクリプト |
| [Samples/mptf32problem/README.md](../../Samples/mptf32problem/README.md) | 旧 `--mp`（TF32）で GW1500 に出た NaN の記録 |
| [GW1500_failures.md](../GW1500_failures.md) | GW1500 の失敗と fp32 での回し直し |
| ecaljclaude.md の [GPU 開発の教訓](../ecaljclaude.md) | 非同期の落とし穴（WB.4 の件、OpenACC の注意） |
