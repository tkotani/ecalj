# MD/ — 開発の記録（ここから読む）

ecalj の開発者と Claude のための文書。公開サイト（ecaljdoc）には出さない。使い方と理論の本体は [`ecaljdoc/`](../ecaljdoc/index.md)（公開サイト https://ecalj.github.io/ecaljdoc/）。

## 読む順（2026-10-02）

1. [ecaljclaude.md](ecaljclaude.md) — 決まり（記録の方針、コーディング規約）と**文書の地図**。どの文書に何があるかはここで引く
2. [handover.md](handover.md) — 引き継ぎ: 計算機、user との取り決め、ビルドと実行の落とし穴
3. [TODOandQuestion.md](TODOandQuestion.md) — やりかけ・未着手・判断待ち、やったこと
4. [research_kotani_log.md](research_kotani_log.md) — 経緯と計測。ローカルなメモを含む（最新が上、`## 日付` の下に `### 時刻`）
5. 変更の正本: 最上位の [`Changes.md`](../Changes.md)（新しい順）。その要約は ecaljdoc の [manual/whatsnew.md](../ecaljdoc/manual/whatsnew.md)

## 話題ごとのジャンプ

各表は「仕様（ecaljdoc、正本）→ サンプル → 開発の記録」の順。

### 有限温度（電子温度と均し方、`t_tetrakbt`・`t_sigmaw`・`wcsmear`）

メニューは [kBT/README.md](kBT/README.md)。

| 行き先 | 中身 |
| --- | --- |
| ecaljdoc [manual/kBT.md](../ecaljdoc/manual/kBT.md) | キーと式。§0 一覧、§2 χ0 側（§2.5 キーの変遷）、§3 Σ 側、§3.5 金属の荒れと `wcsmear`、§7 検証と限界 |
| [Samples/kBT/README.md](../Samples/kBT/README.md) | サンプルの一覧。温度のスキャン [scanT/](../Samples/kBT/scanT/README.md)（Si・GaAs・Fe・Cu）、Fe の対照 [Fe/](../Samples/kBT/Fe/README.md) |
| [Samples/kBT/LiTi2O4/README.md](../Samples/kBT/LiTi2O4/README.md) | LiTi₂O₄ の MLO-QSGW を 40 反復（6³・9³、tf32・fp32・fp64）。全反復の図は [six_patterns.md](../Samples/kBT/LiTi2O4/six_patterns.md) |
| [kBT/history.md](kBT/history.md) | 2026-06〜09 の経緯と直した誤り（古いキー `tetrakbt`・`t_sigmakbt`・`esmr`） |

### MLO（模型、SOC、U・J・cRPA）

メニューは [mlo_notes/README.md](mlo_notes/README.md)。

| 行き先 | 中身 |
| --- | --- |
| ecaljdoc [manual/mlo.md](../ecaljdoc/manual/mlo.md) | §1 最終形と入力キー、§3 部分バンドの模型、§4 SOC（`job_mlo_soc`）、§5 FeMgO（空格子球）、§6 `job_mloW`（U・J、cRPA、Löwdin の標準、表 M8）、§9 模型の選び方と検査（`mlo_bandcheck.py`、式 (12) 帯のとげ） |
| [Samples/MLOsamples/README.md](../Samples/MLOsamples/README.md) | 模型のサンプル（`job_mlo`）と参照 |
| [wannier_vs_mlo.md](wannier_vs_mlo.md) | 最大局在 Wannier 関数（2026-10-02 に外した）と MLO の比較。cRPA の U の差の理由 |
| [mlo_backup.md](mlo_backup.md) | MLO の経緯（method 0〜4、否定された案、途中のバグ） |
| [TOOLS/gadget/README.md](../TOOLS/gadget/README.md) | 標準に入れない試作: MLO の最大局在化（`mlo_maxloc.py`） |

### MLO-QSGW（Σ を MLO で持って内挿する QSGW、`gwsc --mlo`）

| 行き先 | 中身 |
| --- | --- |
| ecaljdoc [manual/mlo_gwsc.md](../ecaljdoc/manual/mlo_gwsc.md) | 方法と設定 |
| [Samples/MLOQSGW/README.md](../Samples/MLOQSGW/README.md) | GaAs・NiO の小さな試験 |
| [Samples/kBT/LiTi2O4/README.md](../Samples/kBT/LiTi2O4/README.md) | 本番の計算（上の有限温度の表と同じ） |
| [mlo_notes/sigma_mlo_design.md](mlo_notes/sigma_mlo_design.md) | 設計書（式 (1)–(17)、現行は §9〜§13） |

### マグノン（`job_mlo_magnon`）

| 行き先 | 中身 |
| --- | --- |
| ecaljdoc [manual/mlo.md](../ecaljdoc/manual/mlo.md) §6 | 相互作用 W の作り方（マグノンもこの W を使う）。Löwdin の基底でマグノンが窓によらなくなったこと、Goldstone の条件の W の倍率 η（Fe 1.24、FeCo 1.28、Ni 1.75） |
| [Samples/Magnon/README.md](../Samples/Magnon/README.md) | サンプルの一覧。Wannier 版（`job_magnon`）は git のタグ `last-wannier` |
| [Samples/Magnon/Fe_mlo_magnon/README.md](../Samples/Magnon/Fe_mlo_magnon/README.md) | bcc Fe。回し方、Wannier 版との比較の表 2・図 1（q ≤ 0.3 で MLO が 2 割高い、q ≥ 0.4 はほぼ重なる） |
| [TODOandQuestion.md](TODOandQuestion.md) §1 | 未着手「Fe のマグノン: 小さい q で 2 割高い件と、実験との比較」と、見る所 |
| [research_kotani_log.md](research_kotani_log.md) 2026-10-02 | 04:17・07:10・09:42（MLO の窓）、11:12（Löwdin の基底で窓によらなくなった）、11:35（Löwdin をメインにする判断） |

### そのほか

| 話題 | 行き先 |
| --- | --- |
| GW の GPU 高速化（`--prec=tf32\|fp32\|fp64`） | メニュー [gpu/README.md](gpu/README.md)、ecaljdoc [manual/ecaljgpu.md](../ecaljdoc/manual/ecaljgpu.md) |
| 対称性（spglib） | [symmetry_spglib.md](symmetry_spglib.md)、ecaljdoc [manual/lmf.md](../ecaljdoc/manual/lmf.md) の SYMGRP |
| GW の実装の覚え書き（χ0・W・Σ） | メニュー [implementation/README.md](implementation/README.md)、[MemoCode.md](MemoCode.md)、コードの大局 [module_map.md](module_map.md) |
| 試験（`testecalj`、Samples の組） | [developer.md](developer.md)、[testecalj_2025.md](testecalj_2025.md)、[`TOOLS/samples_tests.sh`](../TOOLS/samples_tests.sh) |
| GW1500（量産と失敗の分類） | まとめの正本 [`ecalj_auto/GW1500_status.md`](../ecalj_auto/GW1500_status.md)（経過、最終の状態の表 2、失敗の分類と回し直し、構造の不正な 16 物質）、分類ごとの一覧 [`ecalj_auto/GW1500_by_category.md`](../ecalj_auto/GW1500_by_category.md)、5 月の失敗の生の記録 [GW1500_failures.md](GW1500_failures.md)、自動の流れ [auto.md](auto.md) |
| 開発の手引き（push、ビルド、ジョブの投入） | [ForDevelopers.md](ForDevelopers.md)、研究の要約 [ForDevelopers_research.md](ForDevelopers_research.md) |
| ecaljdoc の公開 | [ecaljdoc_publish.md](ecaljdoc_publish.md) |
| 片付けたもの（trash の表、ノウハウ） | [past_log.md](past_log.md)、履歴の書き換えのハッシュの対応 [commit_map_20261002.txt](commit_map_20261002.txt) |
