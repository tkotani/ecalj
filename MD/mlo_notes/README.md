# MLO — 過去の作業記録(backup)

**これは記録であって手引きではない。** 現在の方法と推奨値は
[ecaljdoc manual/mlo](https://ecalj.github.io/ecaljdoc/manual/mlo) を見よ。
経緯(method 0〜4 の関係、否定された案、途中で見つかったバグ)は
[mlo_backup.md](../mlo_backup.md) にまとめてある。

ここにあるのはさらにその前段階の生の作業記録で、**古い設定や古い評価器で
測った数値を含む**。現行の実測値と食い違うときは manual/mlo が正しい。

| ファイル | 時期 | 中身 |
|---|---|---|
| [`MLO_theory.md`](MLO_theory.md) | 2026-09-13 | method 0/1/2/3 が同じ一本の式だと気づいた回。§4 の結論はその後 Al₂O₃:Cr で覆された |
| [`MLO_optimization_log.md`](MLO_optimization_log.md) | 2026-09 | パラメータ走査の生ログ |
| [`MLO_v6_report.md`](MLO_v6_report.md) (+ `MLO_v6_figs/`) | 2026-09 | v6(自動窓)段階の報告。method 3 を採っていた頃のもので、method 3 は最終的に不採用 |
| [`mlo_nskip_cu_problem.md`](mlo_nskip_cu_problem.md) (+ `mlo_nskip_*.png`) | 2026-09-17/18 | Cu の d 模型の折れの正体（`nskip` の k 依存）。`nskip` を k の最小にした理由 |
| [`mp_20260918/README.md`](mp_20260918/README.md) (+ 図) | 2026-09-18 | Materials Project の 10 結晶で MLO の既定を試した記録 |
| [`sigma_mlo_design.md`](sigma_mlo_design.md) | 2026-09-24〜 | MLO-QSGW（Σ を MLO で持って内挿する）の設計書。現行は §9〜§13。2026-10-02 に MD/kBT/ から移した |
| [`mlo_soc_memo.md`](mlo_soc_memo.md) | 2026-09 | スピン軌道を摂動で入れる MLO の覚え書き（もとは `Samples/MLOsamples/README_SOC.md`） |

## サンプルのそばにある記録（2026-10-02）

いまの結果と参照はサンプルの側にある。ここの古い数値と比べるときはそちらを見る。

- [`Samples/MLOsamples/README.md`](../../Samples/MLOsamples/README.md) — MLO の模型のサンプル（`job_mlo`、検査 `mlo_bandcheck.py`）
- [`Samples/MLOQSGW/README.md`](../../Samples/MLOQSGW/README.md) — MLO-QSGW（`gwsc --mlo`）
- [`Samples/Magnon/Fe_mlo_magnon/README.md`](../../Samples/Magnon/Fe_mlo_magnon/README.md) — MLO のマグノン（bcc Fe、Wannier 関数との比較）
- ecaljdoc: [manual/mlo.md](../../ecaljdoc/manual/mlo.md)、[manual/mlo_gwsc.md](../../ecaljdoc/manual/mlo_gwsc.md)
