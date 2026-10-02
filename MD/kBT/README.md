# MD/kBT — 有限温度（電子温度と均し方）のメニュー

2026-10-02 22:15 に整理した。ここに置くのは開発の記録 [history.md](history.md) だけで、ほかはジャンプ先。

## いまの要点

- **キー**（`[gw]`、単位は K）: χ0 側 `t_tetrakbt`（0 = T=0 のテトラヘドロン法、T > 0 = 有限温度のテトラヘドロン法、−T = Im χ0 を Gaussian で均す）、
  Σ 側 `t_sigmaw`（Σ の中間準位を均す Fermi–Dirac 核、既定 1000）。金属の荒れ（小さい q の W のプラズモン極を Σc の実軸極項が拾う）は `wcsmear`（既定）で均す。
  式と表は ecaljdoc [kBT.md](../../ecaljdoc/manual/kBT.md)（§0 一覧、§2 χ0 側、§2.5 キーの変遷、§3 Σ 側、§3.5 `wcsmear`、§7 検証と限界）
- **結果**: LiTi₂O₄ の MLO-QSGW を LDA から 40 反復（6³・9³、`gwsc --prec=tf32|fp32|fp64`）。精度の差は数 meV 以下で、反復やメッシュの差より 1 桁小さい
  → [Samples/kBT/LiTi2O4/README.md](../../Samples/kBT/LiTi2O4/README.md)
- 2026-06 の計算（`tetrakbt`・`t_sigmakbt`・`esmr` の当時のキー）の説明と、後で覆った解釈は [history.md](history.md)

## メニュー

**表 1**. 使い方とサンプル

| 行き先 | 中身 |
| --- | --- |
| ecaljdoc [kBT.md](../../ecaljdoc/manual/kBT.md) | キーの意味と式（正本） |
| [Samples/kBT/README.md](../../Samples/kBT/README.md) | サンプルの一覧（温度のスキャン `scanT/`、Fe の対照、LiTi₂O₄、試験のターゲット） |
| [Samples/kBT/scanT/README.md](../../Samples/kBT/scanT/README.md) | 温度だけを変えて回す Si・GaAs・Fe・Cu |
| [Samples/kBT/Fe/README.md](../../Samples/kBT/Fe/README.md) | Σ 側の温度だけを変えた bcc Fe の対照 |

**表 2**. LiTi₂O₄（MLO-QSGW）

| 行き先 | 中身 |
| --- | --- |
| [Samples/kBT/LiTi2O4/README.md](../../Samples/kBT/LiTi2O4/README.md) | 40 反復のまとめ（2026-09-30。表 1〜5、図 1・2） |
| [Samples/kBT/LiTi2O4/six_patterns.md](../../Samples/kBT/LiTi2O4/six_patterns.md) | {6³, 9³} × {tf32, fp32, fp64} のすべての反復の図と表 |
| [Samples/kBT/LiTi2O4/input/qmlo/README.md](../../Samples/kBT/LiTi2O4/input/qmlo/README.md) | 本番の入力、環境の作り方、反復ごとの正常値 |
| ecaljdoc [mlo_gwsc.md](../../ecaljdoc/manual/mlo_gwsc.md)、設計書 [sigma_mlo_design.md](../mlo_notes/sigma_mlo_design.md) | MLO-QSGW の方法 |

**表 3**. 記録（古いキーで書いてある）

| 行き先 | 中身 |
| --- | --- |
| [history.md](history.md) | 2026-06〜09 の経緯、当時の計算の説明、直した誤り（もとは ecaljdoc の kBT の頁の一部） |
| [Samples/kBT/finiteT_202606.md](../../Samples/kBT/finiteT_202606.md) | 2026-06 の計算の注意と、kt1 に残っているもの |
| [Samples/kBT/LiTi2O4/finiteT_202606.md](../../Samples/kBT/LiTi2O4/finiteT_202606.md) | LiTi₂O₄ の 2026-06 の run（`n666_T*`、`n999_T1000`）の説明と結果 |
| [research_log.md](../research_log.md) | 時刻入りの経緯（2026-09-16〜） |

GPU の高速化の記録は [gpu/](../gpu/README.md)。

history.md の §3 と Samples/kBT/LiTi2O4/finiteT_202606.md の §1・§3 は、同じ 2026-06 の run の説明で一部が重なる（別々に注記が足されていて、行で 24 % が同じ）。
