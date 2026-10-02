# Samples/AFfixMMOM — 反強磁性のモーメントを決めた値に保つ計算（コードの中の名前は AFTEST）

反強磁性の対の 2 つのサイトに、互いに逆向きの場をかけ、毎反復その強さを直して、モーメント (m₁ − m₂)/2 を指定した値に保つ。
モーメントが不安定な AF の金属（2020〜2021 年の NiSe の QSGW）や、モーメントを変えたときの全エネルギー E(m) を見るのに使う。

| ディレクトリ | 中身 |
|---|---|
| `NiO_afsym` | NiO（AF II）、目標 1.6 μB、`symgrpaf` あり（スピン 2 はスピン 1 から作る） |
| `NiO_noafsym` | 同じで `symgrpaf` なし（両スピンを解く）。結果は `NiO_afsym` と 10⁻⁶ eV で一致する |

## 使い方

ecaljdoc `manual/UsageDetailed.md` の「Holding the AF moment at a given value (`mmtarget.aftest`)」。作業ディレクトリの `mmtarget.aftest`（目標のモーメント）で
このモードになる。対のサイトに LDA+U のブロックが要る（`idu = 1`、`uh = jh = 0` でよい）。回し直すときは `mmagfield.aftest`・`mixmag.aftest` を消す。

## 注意

- 2026-10-01 に直した（`aa4c24129`）: 場をスピン 1 にしか入れていなかった。いまは ehk が制約の下の全エネルギー
- 残り（対の決め打ち、場の更新の利得、lmf の起動ごとの更新）は ecalj の `ecaljdoc/MD/TODOandQuestion.md`
- 調べた記録: ecalj の `ecaljdoc/MD/research_log.md` 2026-10-01 朝 06:46（表 06:46-1）
