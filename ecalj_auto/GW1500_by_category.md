# GW1500 の分類ごとの一覧

`gw1500_bycat.py` が [`gw1500_notes_20261001.tsv`](gw1500_notes_20261001.tsv) から作る（手で直さない）。分類の意味と経過は [GW1500_status.md](GW1500_status.md) §5 の表 6。採用のギャップ: 回し直し（fp32、2026-09-30〜）が収束していればその値、金属として収束したものは 0、でなければ 5 月が GOOD なら 5 月の値（TF32）。QSGW80（`scaledsigma = 0.8`）、eV。mpid は Materials Project へのリンク。

**表 1**. 分類の数

| 分類 | 数 | 意味 |
| --- | --- | --- |
| [`TOO_LARGE`](#too-large) | 2 | 原子が大きく間隔を置いて並んだだけの大きい胞（user「原子を間隔を置いていったようなもの」）。Rb8 と Sr mp-1056418（Sr は 2026-10-07 18:51 から。それまで TOO_THEORETICAL、1 原子 462 Å³、QSGW80 の値は残し MLO は作らない）。Rb8: 胞が大きすぎる（一辺 19.9 Å、4410 Å³、平面波 3.2 万、lmf --jobgw=1 が 1 ランク 35 GB）。構造も Rb-IV の高圧相を MP が圧力ゼロで緩和したもので、Rb8 の環が離れて並ぶ。2026-10-06 19:28 までの名前は INVALID_STRUCTURE（user「大きすぎる、というべきかな」）。その前（10-01〜10-06）は別の化合物の副格子 11 件も入っていた |
| [`SUSPECT_STRUCTURE`](#suspect-structure) | 4 | ICSD の備考には出ないが、非整合層状化合物の副格子と同じ形（1 原子の体積が岩塩型の約 2〜4 倍、配位 4〜5） |
| [`MAY_WRONG`](#may-wrong) | 1 | 5 月のギャップが誤り（LDA より小さい）。回し直しの値を使う |
| [`FAILED_MAY`](#failed-may) | 145 | 5 月に使える結果が無かった。回し直しで全部収束 |
| [`UNKNOWN_MAY`](#unknown-may) | 10 | LDA でもギャップが無いか 0.1 eV 程度（金属・半金属） |
| [`SUSPECT_GOOD`](#suspect-good) | 18 | 5 月は GOOD だが、QSGW80 < LDA または 2 反復目以降に 0.5 eV を超えて振動。回し直しで全部収束 |
| [`NOTCONV_MAY`](#notconv-may) | 179 | 5 月は 10 反復で収束しなかった。回し直しで全部収束 |
| [`DRIFT_GOOD`](#drift-good) | 62 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに 0.08 eV 以上動いた。回し直していない |
| [`GOOD`](#good) | 1125 | 5 月に収束、注記なし。回し直していない |
| 計 | 1546 | |

## TOO_LARGE

原子が大きく間隔を置いて並んだだけの大きい胞（user「原子を間隔を置いていったようなもの」）。Rb8 と Sr mp-1056418（Sr は 2026-10-07 18:51 から。それまで TOO_THEORETICAL、1 原子 462 Å³、QSGW80 の値は残し MLO は作らない）。Rb8: 胞が大きすぎる（一辺 19.9 Å、4410 Å³、平面波 3.2 万、lmf --jobgw=1 が 1 ランク 35 GB）。構造も Rb-IV の高圧相を MP が圧力ゼロで緩和したもので、Rb8 の環が離れて並ぶ。2026-10-06 19:28 までの名前は INVALID_STRUCTURE（user「大きすぎる、というべきかな」）。その前（10-01〜10-06）は別の化合物の副格子 11 件も入っていた。2 物質。

**表 2**. 1 原子の体積（Å³）・凸包からの距離（eV/原子）は MP の値。ICSD の備考は MP の項目のもの

| mpid | 組成 | 原子 | 1 原子の体積 | 凸包から | 5 月 | 回し直し | ICSD の備考 | 特徴 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| [mp-1056418](https://next-gen.materialsproject.org/materials/mp-1056418) | Sr | 1 | 462.2 | 1.400 | FAILED | CONVERGED 1.671738 | Strontium; Strontium cobalt oxide (1.26/0.98/2.91) - (Sr) part | 構造: Sr–Co–O の「(Sr) part」を抜き出したもの。10.99×10.99×4.42 Å の胞に Sr 1 個（Sr の鎖、隣 2 個、4.42 Å）。1 原子 462 Å³（fcc Sr 56）、凸包から 1.40 eV/原子 / 胞が大きい（1 原子 462 Å³）。QSGW80 の値は残し、MLO は作らない（skipped。2026-10-07 18:51 まで TOO_THEORETICAL） / 回し直し(run2m): CONVERGED 10 反復、ギャップ 1.671738 eV（LDA 0.436487）。構造が不正なので値に物理的な意味は無い |
| [mp-1179832](https://next-gen.materialsproject.org/materials/mp-1179832) | Rb8 | 8 | 551.3 | 0.453 | FAILED | FAIL(143) none | High pressure experimental phase; Rubidium - IV, HP | Rb-IV（高圧相、ICSD 109016）を MP が圧力ゼロで緩和。一辺 19.9 Å の胞に Rb8 の正八角形の環（Rb–Rb 4.63 Å、隣 2 個）。1 原子 551 Å³（bcc Rb 93）、凸包から 0.45 eV/原子。ギャップは環の HOMO–LUMO。平面波 3.2 万で lmf --jobgw=1 が 1 ランク 35 GB / 回し直し(run2m): FAIL(143) 0 反復、ギャップ none eV（LDA 0.104062）。構造が不正なので値に物理的な意味は無い |

## SUSPECT_STRUCTURE

ICSD の備考には出ないが、非整合層状化合物の副格子と同じ形（1 原子の体積が岩塩型の約 2〜4 倍、配位 4〜5）。4 物質。

**表 3**. 1 原子の体積（Å³）・凸包からの距離（eV/原子）は MP の値。ICSD の備考は MP の項目のもの

| mpid | 組成 | 原子 | 1 原子の体積 | 凸包から | 5 月 | 回し直し | ICSD の備考 | 特徴 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| [mp-20526](https://next-gen.materialsproject.org/materials/mp-20526) | Pb2S2 | 4 | 50.2 | 0.058 | GOOD 3.238 |  | Lead sulfide | PbS。a=4.22、c=11.3 Å、1 原子 50.2 Å³、配位 4。非整合層状化合物の PbS 層と同じ形 |
| [mp-561320](https://next-gen.materialsproject.org/materials/mp-561320) | Pb2S2 | 4 | 98.9 | 0.058 | GOOD 4.120 |  | Lead sulfide | PbS。a=4.22、c=22.3 Å、1 原子 98.9 Å³（岩塩型 26.7）、配位 5。非整合層状化合物の PbS 層と同じ形 |
| [mp-22009](https://next-gen.materialsproject.org/materials/mp-22009) | Pb2Se2 | 4 | 58.5 | 0.073 | GOOD 2.926 |  | Lead selenide | PbSe。1 原子 58.5 Å³、配位 5（岩塩型 PbSe は約 29.5、配位 6）。層状化合物の副格子と見られる |
| [mp-8936](https://next-gen.materialsproject.org/materials/mp-8936) | Sn2Se2 | 4 | 50.5 | 0.088 | GOOD 1.952 |  | Tin selenide | SnSe。a=4.29、c=11.0 Å、1 原子 50.5 Å³、配位 5。非整合層状化合物の SnSe 層と同じ形 |

## MAY_WRONG

5 月のギャップが誤り（LDA より小さい）。回し直しの値を使う。1 物質。

**表 4**

| mpid | 組成 | 原子 | 採用のギャップ | どちらの値 | LDA | 注記 |
| --- | --- | --- | --- | --- | --- | --- |
| [mp-546711](https://next-gen.materialsproject.org/materials/mp-546711) | CsClO4 | 6 | 8.953423 | 回し直し | 5.412798 | 5 月は GOOD でギャップ 0.26 eV（LDA 5.41 より小さい）。回し直しで 8.95 eV。5 月の値は誤り / 回し直し(run1): CONVERGED 5 反復、ギャップ 8.953423 eV（LDA 5.412798）、5 月との差 +8.69 eV |

## FAILED_MAY

5 月に使える結果が無かった。回し直しで全部収束。145 物質。

**表 5**

| mpid | 組成 | 原子 | 採用のギャップ | どちらの値 | LDA | 注記 |
| --- | --- | --- | --- | --- | --- | --- |
| [mp-1342](https://next-gen.materialsproject.org/materials/mp-1342) | BaO | 2 | 4.188265 | 回し直し | 1.992176 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 4.188265 eV（LDA 1.992176）、5 月との差 -0.00 eV |
| [mp-1500](https://next-gen.materialsproject.org/materials/mp-1500) | BaS | 2 | 3.920177 | 回し直し | 2.052343 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 3.920177 eV（LDA 2.052343）、5 月との差 +0.10 eV |
| [mp-1253](https://next-gen.materialsproject.org/materials/mp-1253) | BaSe | 2 | 3.573376 | 回し直し | 1.872121 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 3.573376 eV（LDA 1.872121）、5 月との差 +0.04 eV |
| [mp-2600](https://next-gen.materialsproject.org/materials/mp-2600) | BaTe | 2 | 1.947573 | 回し直し | 0.766467 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 1.947573 eV（LDA 0.766467）、5 月との差 +0.01 eV |
| [mp-1415](https://next-gen.materialsproject.org/materials/mp-1415) | CaSe | 2 | 3.848115 | 回し直し | 1.916106 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 3.848115 eV（LDA 1.916106）、5 月との差 -0.04 eV |
| [mp-1519](https://next-gen.materialsproject.org/materials/mp-1519) | CaTe | 2 | 2.963353 | 回し直し | 1.451984 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 2.963353 eV（LDA 1.451984）、5 月との差 -0.00 eV |
| [mp-2667](https://next-gen.materialsproject.org/materials/mp-2667) | CsAu | 2 | 2.366361 | 回し直し | 0.928034 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 2.366361 eV（LDA 0.928034）、5 月との差 -0.01 eV |
| [mp-22906](https://next-gen.materialsproject.org/materials/mp-22906) | CsBr | 2 | 6.971116 | 回し直し | 4.264978 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 6.971116 eV（LDA 4.264978）、5 月との差 +0.01 eV |
| [mp-571222](https://next-gen.materialsproject.org/materials/mp-571222) | CsBr | 2 | 7.049812 | 回し直し | 4.104366 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 7.049812 eV（LDA 4.104366）、5 月との差 +0.24 eV |
| [mp-22865](https://next-gen.materialsproject.org/materials/mp-22865) | CsCl | 2 | 7.977487 | 回し直し | 5.041066 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 7.977487 eV（LDA 5.041066）、5 月との差 +0.06 eV |
| [mp-573697](https://next-gen.materialsproject.org/materials/mp-573697) | CsCl | 2 | 7.812620 | 回し直し | 4.630841 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 7.812620 eV（LDA 4.630841）、5 月との差 +0.07 eV |
| [mp-1784](https://next-gen.materialsproject.org/materials/mp-1784) | CsF | 2 | 8.871848 | 回し直し | 5.245336 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 6 反復、ギャップ 8.871848 eV（LDA 5.245336）、5 月との差 +0.15 eV |
| [mp-8455](https://next-gen.materialsproject.org/materials/mp-8455) | CsF | 2 | 9.636459 | 回し直し | 6.006345 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 6 反復、ギャップ 9.636459 eV（LDA 6.006345）、5 月との差 -0.00 eV |
| [mp-632319](https://next-gen.materialsproject.org/materials/mp-632319) | CsH | 2 | 4.701782 | 回し直し | 2.302430 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 6 反復、ギャップ 4.701782 eV（LDA 2.302430）、5 月との差 +0.11 eV |
| [mp-614603](https://next-gen.materialsproject.org/materials/mp-614603) | CsI | 2 | 6.307015 | 回し直し | 3.595345 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 6.307015 eV（LDA 3.595345）、5 月との差 +0.10 eV |
| [mp-23251](https://next-gen.materialsproject.org/materials/mp-23251) | KBr | 2 | 7.277473 | 回し直し | 4.267427 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 7.277473 eV（LDA 4.267427）、5 月との差 +0.09 eV |
| [mp-570891](https://next-gen.materialsproject.org/materials/mp-570891) | KBr | 2 | 6.517218 | 回し直し | 3.742652 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 6.517218 eV（LDA 3.742652）、5 月との差 +0.01 eV |
| [mp-22898](https://next-gen.materialsproject.org/materials/mp-22898) | KI | 2 | 6.336210 | 回し直し | 3.673408 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 6.336210 eV（LDA 3.673408）、5 月との差 -0.18 eV |
| [mp-22901](https://next-gen.materialsproject.org/materials/mp-22901) | KI | 2 | 5.433731 | 回し直し | 3.081217 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 5.433731 eV（LDA 3.081217）、5 月との差 -0.04 eV |
| [mp-30373](https://next-gen.materialsproject.org/materials/mp-30373) | RbAu | 2 | 1.478791 | 回し直し | 0.167588 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 1.478791 eV（LDA 0.167588）、5 月との差 -0.01 eV |
| [mp-22867](https://next-gen.materialsproject.org/materials/mp-22867) | RbBr | 2 | 7.103416 | 回し直し | 4.083585 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 7.103416 eV（LDA 4.083585）、5 月との差 +0.10 eV |
| [mp-23303](https://next-gen.materialsproject.org/materials/mp-23303) | RbBr | 2 | 6.712522 | 回し直し | 3.940828 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 6.712522 eV（LDA 3.940828）、5 月との差 -0.01 eV |
| [mp-23295](https://next-gen.materialsproject.org/materials/mp-23295) | RbCl | 2 | 7.982697 | 回し直し | 4.679159 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 7.982697 eV（LDA 4.679159）、5 月との差 +0.04 eV |
| [mp-23299](https://next-gen.materialsproject.org/materials/mp-23299) | RbCl | 2 | 7.885537 | 回し直し | 4.800322 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 7.885537 eV（LDA 4.800322）、5 月との差 -0.02 eV |
| [mp-22903](https://next-gen.materialsproject.org/materials/mp-22903) | RbI | 2 | 6.276975 | 回し直し | 3.563558 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 6.276975 eV（LDA 3.563558）、5 月との差 +0.11 eV |
| [mp-23302](https://next-gen.materialsproject.org/materials/mp-23302) | RbI | 2 | 5.562999 | 回し直し | 3.136800 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 5.562999 eV（LDA 3.136800）、5 月との差 +0.14 eV |
| [mp-1087](https://next-gen.materialsproject.org/materials/mp-1087) | SrS | 2 | 4.172491 | 回し直し | 2.303630 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(try): CONVERGED 5 反復、ギャップ 4.172491 eV（LDA 2.303630）、5 月との差 +0.05 eV |
| [mp-2758](https://next-gen.materialsproject.org/materials/mp-2758) | SrSe | 2 | 3.767455 | 回し直し | 2.077045 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 3.767455 eV（LDA 2.077045）、5 月との差 +0.03 eV |
| [mp-1958](https://next-gen.materialsproject.org/materials/mp-1958) | SrTe | 2 | 2.989890 | 回し直し | 1.654501 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 2.989890 eV（LDA 1.654501）、5 月との差 +0.03 eV |
| [mp-23197](https://next-gen.materialsproject.org/materials/mp-23197) | TlI | 2 | 2.721309 | 回し直し | 1.607907 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 1.61 eV) / 回し直し(run1): CONVERGED 5 反復、ギャップ 2.721309 eV（LDA 1.607907） |
| [mp-568662](https://next-gen.materialsproject.org/materials/mp-568662) | BaCl2 | 3 | 8.002961 | 回し直し | 5.203188 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 8.002961 eV（LDA 5.203188）、5 月との差 +0.09 eV |
| [mp-30031](https://next-gen.materialsproject.org/materials/mp-30031) | CaI2 | 3 | 6.084609 | 回し直し | 3.470842 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 6.084609 eV（LDA 3.470842）、5 月との差 -0.01 eV |
| [mp-7988](https://next-gen.materialsproject.org/materials/mp-7988) | Cs2O | 3 | 2.279375 | 回し直し | 0.831782 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 6 反復、ギャップ 2.279375 eV（LDA 0.831782）、5 月との差 +0.25 eV |
| [mp-8426](https://next-gen.materialsproject.org/materials/mp-8426) | K2Se | 3 | 4.112822 | 回し直し | 1.952488 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 4.112822 eV（LDA 1.952488）、5 月との差 +0.26 eV |
| [mp-1747](https://next-gen.materialsproject.org/materials/mp-1747) | K2Te | 3 | 3.930307 | 回し直し | 1.953498 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 6 反復、ギャップ 3.930307 eV（LDA 1.953498）、5 月との差 +0.13 eV |
| [mp-15687](https://next-gen.materialsproject.org/materials/mp-15687) | KZnAs | 3 | 1.303605 | 回し直し | 0.139934 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 1.303605 eV（LDA 0.139934）、5 月との差 +0.02 eV |
| [mp-505297](https://next-gen.materialsproject.org/materials/mp-505297) | NbSbRu | 3 | 0.675472 | 回し直し | 0.310956 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 4 反復、ギャップ 0.675472 eV（LDA 0.310956）、5 月との差 +0.00 eV |
| [mp-1394](https://next-gen.materialsproject.org/materials/mp-1394) | Rb2O | 3 | 3.328882 | 回し直し | 1.301782 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 6 反復、ギャップ 3.328882 eV（LDA 1.301782）、5 月との差 +0.20 eV |
| [mp-8041](https://next-gen.materialsproject.org/materials/mp-8041) | Rb2S | 3 | 3.922739 | 回し直し | 1.751015 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 3.922739 eV（LDA 1.751015）、5 月との差 +0.08 eV |
| [mp-11327](https://next-gen.materialsproject.org/materials/mp-11327) | Rb2Se | 3 | 3.699034 | 回し直し | 1.592205 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 3.699034 eV（LDA 1.592205）、5 月との差 +0.25 eV |
| [mp-441](https://next-gen.materialsproject.org/materials/mp-441) | Rb2Te | 3 | 3.704224 | 回し直し | 1.660044 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 6 反復、ギャップ 3.704224 eV（LDA 1.660044）、5 月との差 +0.08 eV |
| [mp-23209](https://next-gen.materialsproject.org/materials/mp-23209) | SrCl2 | 3 | 7.816688 | 回し直し | 5.003842 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 7.816688 eV（LDA 5.003842）、5 月との差 +0.07 eV |
| [mp-2697](https://next-gen.materialsproject.org/materials/mp-2697) | SrO2 | 3 | 5.758807 | 回し直し | 2.758681 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 2.76 eV, previous QSGW iteration 5.76 eV) / 回し直し(run1): CONVERGED 5 反復、ギャップ 5.758807 eV（LDA 2.758681） |
| [mp-1207108](https://next-gen.materialsproject.org/materials/mp-1207108) | BaAlGeH | 4 | 1.364712 | 回し直し | 0.773235 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 4 反復、ギャップ 1.364712 eV（LDA 0.773235）、5 月との差 +0.00 eV |
| [mp-10969](https://next-gen.materialsproject.org/materials/mp-10969) | CdCN2 | 4 | 3.816450 | 回し直し | 2.083752 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 2.08 eV, previous QSGW iteration 3.82 eV) / 回し直し(run1): CONVERGED 5 反復、ギャップ 3.816450 eV（LDA 2.083752） |
| [mp-16763](https://next-gen.materialsproject.org/materials/mp-16763) | KYTe2 | 4 | 2.439019 | 回し直し | 1.206675 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 2.439019 eV（LDA 1.206675）、5 月との差 -0.04 eV |
| [mp-18921](https://next-gen.materialsproject.org/materials/mp-18921) | NaCoO2 | 4 | 3.375416 | 回し直し | 0.791976 | 5 月の失敗: AC_lmf-SCF-diverged@llmf_start / 回し直し(run2m(続き)): CONVERGED 14 反復、ギャップ 3.375416 eV（LDA 0.791976）、5 月との差 +0.91 eV |
| [mp-10226](https://next-gen.materialsproject.org/materials/mp-10226) | NaYS2 | 4 | 3.870547 | 回し直し | 2.117783 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 2.12 eV, previous QSGW iteration 3.62 eV) / 回し直し(run1): CONVERGED 5 反復、ギャップ 3.870547 eV（LDA 2.117783） |
| [mp-581833](https://next-gen.materialsproject.org/materials/mp-581833) | RbN3 | 4 | 6.417903 | 回し直し | 3.597908 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 5 反復、ギャップ 6.417903 eV（LDA 3.597908） |
| [mp-14437](https://next-gen.materialsproject.org/materials/mp-14437) | RbYO2 | 4 | 5.951582 | 回し直し | 3.521503 | 5 月の失敗: AC_rst-rmt-mismatch@llmf_start(lmf) / 回し直し(run1): CONVERGED 5 反復、ギャップ 5.951582 eV（LDA 3.521503）、5 月との差 +0.02 eV |
| [mp-9383](https://next-gen.materialsproject.org/materials/mp-9383) | SrHfN2 | 4 | 1.198651 | 回し直し | 0.170895 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 0.17 eV, previous QSGW iteration 1.20 eV) / 回し直し(run1): CONVERGED 5 反復、ギャップ 1.198651 eV（LDA 0.170895） |
| [mp-558189](https://next-gen.materialsproject.org/materials/mp-558189) | Ag3SI | 5 | 1.812281 | 回し直し | 0.479781 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 0.48 eV) / 回し直し(run1): CONVERGED 5 反復、ギャップ 1.812281 eV（LDA 0.479781） |
| [mp-8196](https://next-gen.materialsproject.org/materials/mp-8196) | AgNO3 | 5 | 5.055941 | 回し直し | 1.525729 | 5 月の失敗: PROD_bmix_min / 回し直し(run1): CONVERGED 6 反復、ギャップ 5.055941 eV（LDA 1.525729） |
| [mp-4559](https://next-gen.materialsproject.org/materials/mp-4559) | BaCO3 | 5 | 7.465043 | 回し直し | 4.354890 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 4.36 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 7.465043 eV（LDA 4.354890） |
| [mp-552098](https://next-gen.materialsproject.org/materials/mp-552098) | Bi2SeO2 | 5 | 0.965917 | 回し直し | 0.376638 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 0.38 eV) / 回し直し(run1): CONVERGED 4 反復、ギャップ 0.965917 eV（LDA 0.376638） |
| [mp-625548](https://next-gen.materialsproject.org/materials/mp-625548) | CdH2O2 | 5 | 4.164273 | 回し直し | 1.513245 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 6 反復、ギャップ 4.164273 eV（LDA 1.513245） |
| [mp-22535](https://next-gen.materialsproject.org/materials/mp-22535) | HfPbO3 | 5 | 3.583439 | 回し直し | 2.306512 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 5 反復、ギャップ 3.583439 eV（LDA 2.306512） |
| [mp-12725](https://next-gen.materialsproject.org/materials/mp-12725) | NbAgO3 | 5 | 2.442575 | 回し直し | 1.094585 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 1.09 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 2.442575 eV（LDA 1.094585） |
| [mp-23011](https://next-gen.materialsproject.org/materials/mp-23011) | RbClO3 | 5 | 8.810958 | 回し直し | 5.161027 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 5.16 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 8.810958 eV（LDA 5.161027） |
| [mp-29798](https://next-gen.materialsproject.org/materials/mp-29798) | TlBrO3 | 5 | 5.683209 | 回し直し | 3.387881 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 5 反復、ギャップ 5.683209 eV（LDA 3.387881） |
| [mp-29738](https://next-gen.materialsproject.org/materials/mp-29738) | TlClO3 | 5 | 6.817851 | 回し直し | 4.127036 | 5 月の失敗: PROD_bmix_min / 回し直し(run1): CONVERGED 6 反復、ギャップ 6.817851 eV（LDA 4.127036） |
| [mp-22981](https://next-gen.materialsproject.org/materials/mp-22981) | TlIO3 | 5 | 4.499910 | 回し直し | 2.703091 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 5 反復、ギャップ 4.499910 eV（LDA 2.703091） |
| [mp-14099](https://next-gen.materialsproject.org/materials/mp-14099) | ZnAgF3 | 5 | 5.471937 | 回し直し | 1.448655 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 1.45 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 5.471937 eV（LDA 1.448655） |
| [mp-1068577](https://next-gen.materialsproject.org/materials/mp-1068577) | ZrPbO3 | 5 | 3.434653 | 回し直し | 2.178316 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 2.18 eV, previous QSGW iteration 0.98 eV) / 回し直し(run1): CONVERGED 5 反復、ギャップ 3.434653 eV（LDA 2.178316） |
| [mp-22993](https://next-gen.materialsproject.org/materials/mp-22993) | AgClO4 | 6 | 6.097852 | 回し直し | 2.558213 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 6 反復、ギャップ 6.097852 eV（LDA 2.558213） |
| [mp-551758](https://next-gen.materialsproject.org/materials/mp-551758) | AgClO4 | 6 | 5.576670 | 回し直し | 2.023994 | 5 月の失敗: PROD_bmix_min / 回し直し(run1): CONVERGED 6 反復、ギャップ 5.576670 eV（LDA 2.023994） |
| [mp-551798](https://next-gen.materialsproject.org/materials/mp-551798) | AlPO4 | 6 | 1.339309 | 回し直し | 0.002741 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 7 反復、ギャップ 1.339309 eV（LDA 0.002741） |
| [mp-5606](https://next-gen.materialsproject.org/materials/mp-5606) | AlTlF4 | 6 | 7.227966 | 回し直し | 4.104807 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 6 反復、ギャップ 7.227966 eV（LDA 4.104807） |
| [mp-552934](https://next-gen.materialsproject.org/materials/mp-552934) | Ba2CuBrO2 | 6 | 4.393846 | 回し直し | 2.226497 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 6 反復、ギャップ 4.393846 eV（LDA 2.226497） |
| [mp-551456](https://next-gen.materialsproject.org/materials/mp-551456) | Ba2CuClO2 | 6 | 4.464731 | 回し直し | 2.255917 | 5 月の失敗: AC_NaN-in-Sigma(lcombined/lqpe) / 回し直し(run1): CONVERGED 6 反復、ギャップ 4.464731 eV（LDA 2.255917） |
| [mp-3307](https://next-gen.materialsproject.org/materials/mp-3307) | BaSO4 | 6 | 7.115837 | 回し直し | 4.125531 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 6 反復、ギャップ 7.115837 eV（LDA 4.125531） |
| [mp-1018656](https://next-gen.materialsproject.org/materials/mp-1018656) | Ca2H3Br | 6 | 3.825459 | 回し直し | 1.396659 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 1.40 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 3.825459 eV（LDA 1.396659） |
| [mp-36248](https://next-gen.materialsproject.org/materials/mp-36248) | H4BrN | 6 | 6.846275 | 回し直し | 4.230708 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 5 反復、ギャップ 6.846275 eV（LDA 4.230708） |
| [mp-30302](https://next-gen.materialsproject.org/materials/mp-30302) | Hf2Br2N2 | 6 | 3.453893 | 回し直し | 1.851679 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 5 反復、ギャップ 3.453893 eV（LDA 1.851679） |
| [mp-567441](https://next-gen.materialsproject.org/materials/mp-567441) | Hf2I2N2 | 6 | 2.304422 | 回し直し | 0.866233 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 0.86 eV) / 回し直し(run1): CONVERGED 5 反復、ギャップ 2.304422 eV（LDA 0.866233） |
| [mp-567299](https://next-gen.materialsproject.org/materials/mp-567299) | KC2IN2 | 6 | 7.592745 | 回し直し | 4.402887 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 5 反復、ギャップ 7.592745 eV（LDA 4.402887） |
| [mp-643043](https://next-gen.materialsproject.org/materials/mp-643043) | Rb2H2O2 | 6 | 6.427641 | 回し直し | 2.998063 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 3.00 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 6.427641 eV（LDA 2.998063） |
| [mp-1077553](https://next-gen.materialsproject.org/materials/mp-1077553) | Sr2BBrN2 | 6 | 4.096167 | 回し直し | 2.255133 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 5 反復、ギャップ 4.096167 eV（LDA 2.255133） |
| [mp-10497](https://next-gen.materialsproject.org/materials/mp-10497) | Sr2C4 | 6 | 4.424389 | 回し直し | 2.117045 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 2.12 eV, previous QSGW iteration 1.71 eV) / 回し直し(run1): CONVERGED 5 反復、ギャップ 4.424389 eV（LDA 2.117045） |
| [mp-1019269](https://next-gen.materialsproject.org/materials/mp-1019269) | Sr2H3I | 6 | 3.964839 | 回し直し | 1.793758 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 1.79 eV, previous QSGW iteration 1.71 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 3.964839 eV（LDA 1.793758） |
| [mp-7811](https://next-gen.materialsproject.org/materials/mp-7811) | SrSO4 | 6 | 6.584294 | 回し直し | 3.554214 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 3.56 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 6.584294 eV（LDA 3.554214） |
| [mp-30530](https://next-gen.materialsproject.org/materials/mp-30530) | TlClO4 | 6 | 7.632298 | 回し直し | 4.537533 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 5 反復、ギャップ 7.632298 eV（LDA 4.537533） |
| [mp-27985](https://next-gen.materialsproject.org/materials/mp-27985) | Y2Cl2O2 | 6 | 6.979708 | 回し直し | 4.890571 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 4.89 eV) / 回し直し(run1): CONVERGED 5 反復、ギャップ 6.979708 eV（LDA 4.890571） |
| [mp-23052](https://next-gen.materialsproject.org/materials/mp-23052) | Zr2I2N2 | 6 | 2.876096 | 回し直し | 1.183343 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 1.18 eV) / 回し直し(run1): CONVERGED 5 反復、ギャップ 2.876096 eV（LDA 1.183343） |
| [mp-546862](https://next-gen.materialsproject.org/materials/mp-546862) | HgPb2Cl2O2 | 7 | 3.636553 | 回し直し | 2.227586 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 2.23 eV) / 回し直し(run1): CONVERGED 5 反復、ギャップ 3.636553 eV（LDA 2.227586） |
| [mp-1077901](https://next-gen.materialsproject.org/materials/mp-1077901) | InCo3SnS2 | 7 | 0.272827 | 回し直し | 0.206730 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 4 反復、ギャップ 0.272827 eV（LDA 0.206730） |
| [mp-9583](https://next-gen.materialsproject.org/materials/mp-9583) | K2ZnF4 | 7 | 8.862448 | 回し直し | 4.233235 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 4.23 eV, previous QSGW iteration 8.84 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 8.862448 eV（LDA 4.233235） |
| [mp-571461](https://next-gen.materialsproject.org/materials/mp-571461) | K4ZnAs2 | 7 | 1.781024 | 回し直し | 0.591703 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 0.59 eV, previous QSGW iteration 2.99 eV) / 回し直し(run1): CONVERGED 5 反復、ギャップ 1.781024 eV（LDA 0.591703） |
| [mp-7983](https://next-gen.materialsproject.org/materials/mp-7983) | Mg2SiO4 | 7 | 6.349910 | 回し直し | 3.410704 | 5 月の失敗: AC_NaN-in-Sigma(lcombined/lqpe) / 回し直し(run1): CONVERGED 6 反復、ギャップ 6.349910 eV（LDA 3.410704） |
| [mp-551826](https://next-gen.materialsproject.org/materials/mp-551826) | NbTlBr4O | 7 | 2.666349 | 回し直し | 1.109565 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 1.11 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 2.666349 eV（LDA 1.109565） |
| [mp-7700](https://next-gen.materialsproject.org/materials/mp-7700) | SiB6 | 7 | 0.265274 | 回し直し | 0.044190 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 4 反復、ギャップ 0.265274 eV（LDA 0.044190） |
| [mp-545730](https://next-gen.materialsproject.org/materials/mp-545730) | Ag2C2N2O2 | 8 | 3.436113 | 回し直し | 1.134596 | 5 月の失敗: PROD_bmix_min / 回し直し(run1): CONVERGED 6 反復、ギャップ 3.436113 eV（LDA 1.134596） |
| [mp-6781](https://next-gen.materialsproject.org/materials/mp-6781) | Ag2C2N2O2 | 8 | 5.598843 | 回し直し | 2.975232 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 6 反復、ギャップ 5.598843 eV（LDA 2.975232）、5 月との差 -0.13 eV |
| [mp-552169](https://next-gen.materialsproject.org/materials/mp-552169) | Ag2Cl2O4 | 8 | 3.798930 | 回し直し | 0.651816 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 6 反復、ギャップ 3.798930 eV（LDA 0.651816） |
| [mp-2247](https://next-gen.materialsproject.org/materials/mp-2247) | Ag2N6 | 8 | 3.640062 | 回し直し | 1.384721 | 5 月の失敗: PROD_NaN_lsc / 回し直し(run1): CONVERGED 6 反復、ギャップ 3.640062 eV（LDA 1.384721） |
| [mp-571297](https://next-gen.materialsproject.org/materials/mp-571297) | Ag2N6 | 8 | 3.330038 | 回し直し | 1.470049 | 5 月の失敗: PROD_NaN_lsc / 回し直し(run1): CONVERGED 6 反復、ギャップ 3.330038 eV（LDA 1.470049） |
| [mp-468](https://next-gen.materialsproject.org/materials/mp-468) | Al2F6 | 8 | 13.197701 | 回し直し | 7.890584 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 7.89 eV, previous QSGW iteration 12.87 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 13.197701 eV（LDA 7.890584） |
| [mp-570140](https://next-gen.materialsproject.org/materials/mp-570140) | Au4Br4 | 8 | 3.474921 | 回し直し | 1.881793 | 5 月の失敗: PROD_NaN_lsc / 回し直し(run1): CONVERGED 5 反復、ギャップ 3.474921 eV（LDA 1.881793） |
| [mp-14006](https://next-gen.materialsproject.org/materials/mp-14006) | BaGeF6 | 8 | 10.568682 | 回し直し | 5.740562 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 5.74 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 10.568682 eV（LDA 5.740562） |
| [mp-1080511](https://next-gen.materialsproject.org/materials/mp-1080511) | BaMgC2O4 | 8 | 3.842877 | 回し直し | 1.007577 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 1.01 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 3.842877 eV（LDA 1.007577） |
| [mp-568388](https://next-gen.materialsproject.org/materials/mp-568388) | Bi4I4 | 8 | 1.441402 | 回し直し | 0.640065 | 5 月の失敗: PROD_KILLED / 回し直し(run1): CONVERGED 5 反復、ギャップ 1.441402 eV（LDA 0.640065）、5 月との差 +0.14 eV |
| [mp-23435](https://next-gen.materialsproject.org/materials/mp-23435) | Bi4Te2I2 | 8 | 0.331196 | 回し直し | 0.014058 | 5 月の失敗: PROD_KILLED / 回し直し(run1): CONVERGED 5 反復、ギャップ 0.331196 eV（LDA 0.014058）、5 月との差 +0.09 eV |
| [mp-3078](https://next-gen.materialsproject.org/materials/mp-3078) | Cd2Si2As4 | 8 | 1.378474 | 回し直し | 0.392957 | 5 月の失敗: PROD_NaN_lsc / 回し直し(run1): CONVERGED 5 反復、ギャップ 1.378474 eV（LDA 0.392957） |
| [mp-571100](https://next-gen.materialsproject.org/materials/mp-571100) | Cs2Ag2Br4 | 8 | 3.945622 | 回し直し | 1.581668 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 6 反復、ギャップ 3.945622 eV（LDA 1.581668）、5 月との差 -0.19 eV |
| [mp-510557](https://next-gen.materialsproject.org/materials/mp-510557) | Cs2N6 | 8 | 6.648357 | 回し直し | 4.079606 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 5 反復、ギャップ 6.648357 eV（LDA 4.079606） |
| [mp-8361](https://next-gen.materialsproject.org/materials/mp-8361) | Cs4Te4 | 8 | 2.528535 | 回し直し | 0.930500 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 5 反復、ギャップ 2.528535 eV（LDA 0.930500）、5 月との差 +0.09 eV |
| [mp-571266](https://next-gen.materialsproject.org/materials/mp-571266) | Ga2Cl6 | 8 | 7.681157 | 回し直し | 3.919182 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 6 反復、ギャップ 7.681157 eV（LDA 3.919182）、5 月との差 +0.23 eV |
| [mp-730101](https://next-gen.materialsproject.org/materials/mp-730101) | H8 | 8 | 15.230640 | 回し直し | 9.710196 | 構造: NH4D2PO4 の H（D）だけを抜き出したもの。H2 分子 4 個（H–H 0.74 Å）、密度 0.04 g/cm³ / 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 4 反復、ギャップ 15.230640 eV（LDA 9.710196）、5 月との差 +0.01 eV |
| [mp-9922](https://next-gen.materialsproject.org/materials/mp-9922) | Hf2S6 | 8 | 2.295576 | 回し直し | 0.875094 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 5 反復、ギャップ 2.295576 eV（LDA 0.875094）、5 月との差 +0.04 eV |
| [mp-570219](https://next-gen.materialsproject.org/materials/mp-570219) | In2Br6 | 8 | 4.612794 | 回し直し | 2.334327 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 5 反復、ギャップ 4.612794 eV（LDA 2.334327）、5 月との差 +0.04 eV |
| [mp-22604](https://next-gen.materialsproject.org/materials/mp-22604) | InAsF6 | 8 | 6.950402 | 回し直し | 3.347954 | 5 月の失敗: PROD_NaN_lsc / 回し直し(run1): CONVERGED 6 反復、ギャップ 6.950402 eV（LDA 3.347954） |
| [mp-9273](https://next-gen.materialsproject.org/materials/mp-9273) | K3Sb2Au3 | 8 | 1.624304 | 回し直し | 0.985786 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 4 反復、ギャップ 1.624304 eV（LDA 0.985786）、5 月との差 +0.15 eV |
| [mp-9778](https://next-gen.materialsproject.org/materials/mp-9778) | K4Ag2P2 | 8 | 2.187619 | 回し直し | 1.100432 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 5 反復、ギャップ 2.187619 eV（LDA 1.100432）、5 月との差 -0.02 eV |
| [mp-8446](https://next-gen.materialsproject.org/materials/mp-8446) | K4Cu2P2 | 8 | 2.160216 | 回し直し | 1.073362 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 5 反復、ギャップ 2.160216 eV（LDA 1.073362）、5 月との差 +0.02 eV |
| [mp-9687](https://next-gen.materialsproject.org/materials/mp-9687) | K4P2Au2 | 8 | 1.812703 | 回し直し | 0.806899 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 5 反復、ギャップ 1.812703 eV（LDA 0.806899）、5 月との差 -0.00 eV |
| [mp-863757](https://next-gen.materialsproject.org/materials/mp-863757) | K4Sb2Au2 | 8 | 1.779096 | 回し直し | 0.871825 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 5 反復、ギャップ 1.779096 eV（LDA 0.871825）、5 月との差 -0.01 eV |
| [mp-29025](https://next-gen.materialsproject.org/materials/mp-29025) | Li5Br2N | 8 | 4.487604 | 回し直し | 2.447423 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 2.45 eV) / 回し直し(run1): CONVERGED 5 反復、ギャップ 4.487604 eV（LDA 2.447423） |
| [mp-29151](https://next-gen.materialsproject.org/materials/mp-29151) | Li5NCl2 | 8 | 4.079634 | 回し直し | 2.059977 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 2.06 eV, previous QSGW iteration 3.69 eV) / 回し直し(run1): CONVERGED 5 反復、ギャップ 4.079634 eV（LDA 2.059977） |
| [mp-27419](https://next-gen.materialsproject.org/materials/mp-27419) | LiBiF6 | 8 | 6.731355 | 回し直し | 2.661866 | 5 月の失敗: PROD_bmix_min / 回し直し(run1): CONVERGED 7 反復、ギャップ 6.731355 eV（LDA 2.661866） |
| [mp-20458](https://next-gen.materialsproject.org/materials/mp-20458) | PtPbF6 | 8 | 5.109806 | 回し直し | 1.705999 | 5 月の失敗: PROD_NaN_lsc / 回し直し(run1): CONVERGED 6 反復、ギャップ 5.109806 eV（LDA 1.705999） |
| [mp-14070](https://next-gen.materialsproject.org/materials/mp-14070) | Rb2Al2O4 | 8 | 6.352957 | 回し直し | 3.277843 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 3.28 eV, previous QSGW iteration 0.42 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 6.352957 eV（LDA 3.277843） |
| [mp-7467](https://next-gen.materialsproject.org/materials/mp-7467) | Rb2Cu2O4 | 8 | 3.534965 | 回し直し | 0.800926 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 0.80 eV, previous QSGW iteration 4.91 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 3.534965 eV（LDA 0.800926） |
| [mp-29764](https://next-gen.materialsproject.org/materials/mp-29764) | Rb2H2F4 | 8 | 11.314168 | 回し直し | 6.586519 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 6.58 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 11.314168 eV（LDA 6.586519） |
| [mp-697144](https://next-gen.materialsproject.org/materials/mp-697144) | Rb2H4N2 | 8 | 4.777354 | 回し直し | 1.805169 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 1.80 eV, previous QSGW iteration 0.13 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 4.777354 eV（LDA 1.805169） |
| [mp-9274](https://next-gen.materialsproject.org/materials/mp-9274) | Rb3Sb2Au3 | 8 | 1.827489 | 回し直し | 1.118788 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 5 反復、ギャップ 1.827489 eV（LDA 1.118788）、5 月との差 -0.27 eV |
| [mp-8360](https://next-gen.materialsproject.org/materials/mp-8360) | Rb4Te4 | 8 | 2.270933 | 回し直し | 0.742770 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 5 反復、ギャップ 2.270933 eV（LDA 0.742770）、5 月との差 +0.15 eV |
| [mp-696736](https://next-gen.materialsproject.org/materials/mp-696736) | RbBe2F5 | 8 | 11.951191 | 回し直し | 6.975834 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 6 反復、ギャップ 11.951191 eV（LDA 6.975834）、5 月との差 +0.35 eV |
| [mp-570589](https://next-gen.materialsproject.org/materials/mp-570589) | Se4Br4 | 8 | 3.172827 | 回し直し | 1.394233 | 5 月の失敗: PROD_NaN_lsc / 回し直し(run1): CONVERGED 6 反復、ギャップ 3.172827 eV（LDA 1.394233） |
| [mp-691](https://next-gen.materialsproject.org/materials/mp-691) | Sn4Se4 | 8 | 1.243348 | 回し直し | 0.792248 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 4 反復、ギャップ 1.243348 eV（LDA 0.792248）、5 月との差 +0.01 eV |
| [mp-1080045](https://next-gen.materialsproject.org/materials/mp-1080045) | Sr2H2Br2O2 | 8 | 6.977302 | 回し直し | 4.137954 | 5 月の失敗: PROD_NaN_lsc / 回し直し(run1): CONVERGED 5 反復、ギャップ 6.977302 eV（LDA 4.137954） |
| [mp-555718](https://next-gen.materialsproject.org/materials/mp-555718) | SrNiF6 | 8 | 5.221480 | 回し直し | 1.192299 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 1.19 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 5.221480 eV（LDA 1.192299） |
| [mp-29027](https://next-gen.materialsproject.org/materials/mp-29027) | Ta2I4O2 | 8 | 1.876540 | 回し直し | 0.740647 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 0.74 eV, previous QSGW iteration 0.99 eV) / 回し直し(run1): CONVERGED 5 反復、ギャップ 1.876540 eV（LDA 0.740647） |
| [mp-573321](https://next-gen.materialsproject.org/materials/mp-573321) | Te2Pd2I4 | 8 | 1.137733 | 回し直し | 0.465572 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 5 反復、ギャップ 1.137733 eV（LDA 0.465572）、5 月との差 +0.01 eV |
| [mp-685022](https://next-gen.materialsproject.org/materials/mp-685022) | Te4Pb4 | 8 | 1.062651 | 回し直し | 0.666132 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 4 反復、ギャップ 1.062651 eV（LDA 0.666132）、5 月との差 -0.01 eV |
| [mp-570076](https://next-gen.materialsproject.org/materials/mp-570076) | Tl2CdSnTe4 | 8 | 0.946430 | 回し直し | 0.284114 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 5 反復、ギャップ 0.946430 eV（LDA 0.284114）、5 月との差 +0.03 eV |
| [mp-870](https://next-gen.materialsproject.org/materials/mp-870) | Tl2N6 | 8 | 2.958804 | 回し直し | 1.479471 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 5 反復、ギャップ 2.958804 eV（LDA 1.479471） |
| [mp-570051](https://next-gen.materialsproject.org/materials/mp-570051) | Tl2SnHgTe4 | 8 | 0.645070 | 回し直し | 0.076166 | 5 月の失敗: PROD_KILLED / 回し直し(run1): CONVERGED 4 反復、ギャップ 0.645070 eV（LDA 0.076166）、5 月との差 +0.02 eV |
| [mp-29047](https://next-gen.materialsproject.org/materials/mp-29047) | Tl3VO4 | 8 | 3.588622 | 回し直し | 2.331384 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 5 反復、ギャップ 3.588622 eV（LDA 2.331384） |
| [mp-720](https://next-gen.materialsproject.org/materials/mp-720) | Tl4F4 | 8 | 4.731647 | 回し直し | 2.903315 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 5 反復、ギャップ 4.731647 eV（LDA 2.903315）、5 月との差 +0.03 eV |
| [mp-555911](https://next-gen.materialsproject.org/materials/mp-555911) | TlAsF6 | 8 | 8.464341 | 回し直し | 4.593981 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 4.60 eV) / 回し直し(run1): CONVERGED 6 反復、ギャップ 8.464341 eV（LDA 4.593981） |
| [mp-8348](https://next-gen.materialsproject.org/materials/mp-8348) | TlSbF6 | 8 | 7.804920 | 回し直し | 4.143709 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 6 反復、ギャップ 7.804920 eV（LDA 4.143709）、5 月との差 +6.19 eV |
| [mp-626151](https://next-gen.materialsproject.org/materials/mp-626151) | Y2H2O4 | 8 | 6.107897 | 回し直し | 4.052821 | 5 月の失敗: PROD_bmix_min / 回し直し(run1): CONVERGED 6 反復、ギャップ 6.107897 eV（LDA 4.052821） |
| [mp-571442](https://next-gen.materialsproject.org/materials/mp-571442) | Y2I6 | 8 | 4.836823 | 回し直し | 2.798008 | 5 月の失敗: PROD_TIMEOUT_8h / 回し直し(run1): CONVERGED 5 反復、ギャップ 4.836823 eV（LDA 2.798008）、5 月との差 +0.13 eV |
| [mp-552604](https://next-gen.materialsproject.org/materials/mp-552604) | YBi2ClO4 | 8 | 2.628610 | 回し直し | 1.309736 | 5 月の失敗: PROD_NaN_lqpe / 回し直し(run1): CONVERGED 5 反復、ギャップ 2.628610 eV（LDA 1.309736） |
| [mp-551816](https://next-gen.materialsproject.org/materials/mp-551816) | YBi2IO4 | 8 | 2.460687 | 回し直し | 1.248617 | 5 月の失敗: NOGAP_GAP_COLLAPSED(LDA gap 1.25 eV) / 回し直し(run1): CONVERGED 5 反復、ギャップ 2.460687 eV（LDA 1.248617） |

## UNKNOWN_MAY

LDA でもギャップが無いか 0.1 eV 程度（金属・半金属）。10 物質。

**表 6**

| mpid | 組成 | 原子 | 採用のギャップ | どちらの値 | LDA | 注記 |
| --- | --- | --- | --- | --- | --- | --- |
| [mp-569360](https://next-gen.materialsproject.org/materials/mp-569360) | Hg3 | 3 | 0（金属） | 回し直し | none | 5 月: NOGAP_METAL_IN_LDA_TOO（金属・半金属。gwscconv がギャップを読めず止まった） / 回し直し(run2m): CONVERGED_METAL 5 反復、ギャップ none eV（LDA none） |
| [mp-9339](https://next-gen.materialsproject.org/materials/mp-9339) | Nb2P2 | 4 | 0（金属） | 回し直し | 0.004431 | 5 月: NOGAP_METAL_IN_LDA_TOO（金属・半金属。gwscconv がギャップを読めず止まった） / 回し直し(run2m): CONVERGED_METAL 5 反復、ギャップ none eV（LDA 0.004431） |
| [mp-4573](https://next-gen.materialsproject.org/materials/mp-4573) | TlSbTe2 | 4 | 0.294478 | 回し直し | 0.104106 | 5 月: NOGAP_SMALL_GAP(LDA gap 0.10 eV; possibly semimetal)（金属・半金属。gwscconv がギャップを読めず止まった） / 回し直し(run2m): CONVERGED 4 反復、ギャップ 0.294478 eV（LDA 0.104106） |
| [mp-27710](https://next-gen.materialsproject.org/materials/mp-27710) | CrB4 | 5 | 0（金属） | 回し直し | none | 5 月: NOGAP_METAL_IN_LDA_TOO（金属・半金属。gwscconv がギャップを読めず止まった） / 回し直し(run2m): CONVERGED_METAL 5 反復、ギャップ none eV（LDA none） |
| [mp-9379](https://next-gen.materialsproject.org/materials/mp-9379) | SrSn2As2 | 5 | 0.040954 | 回し直し | 0.139015 | 5 月: NOGAP_SMALL_GAP(LDA gap 0.14 eV; possibly semimetal)（金属・半金属。gwscconv がギャップを読めず止まった） / 回し直し(run2m): CONVERGED 4 反復、ギャップ 0.040954 eV（LDA 0.139015） |
| [mp-1182867](https://next-gen.materialsproject.org/materials/mp-1182867) | AlWO4 | 6 | 0（金属） | 回し直し | none | 5 月: NOGAP_METAL_IN_LDA_TOO（金属・半金属。gwscconv がギャップを読めず止まった） / 回し直し(run2m): CONVERGED_METAL 8 反復、ギャップ none eV（LDA none） |
| [mp-634](https://next-gen.materialsproject.org/materials/mp-634) | Hg3S3 | 6 | 0（金属） | 回し直し | none | 5 月: NOGAP_METAL_IN_LDA_TOO（金属・半金属。gwscconv がギャップを読めず止まった） / 回し直し(run2m): CONVERGED_METAL 5 反復、ギャップ none eV（LDA none） |
| [mp-1969](https://next-gen.materialsproject.org/materials/mp-1969) | Nb2Sb4 | 6 | 0（金属） | 回し直し | none | 5 月: NOGAP_METAL_IN_LDA_TOO（金属・半金属。gwscconv がギャップを読めず止まった） / 回し直し(run2m): CONVERGED_METAL 5 反復、ギャップ none eV（LDA none） |
| [mp-12561](https://next-gen.materialsproject.org/materials/mp-12561) | Ta2As4 | 6 | 0（金属） | 回し直し | none | 5 月: NOGAP_METAL_IN_LDA_TOO（金属・半金属。gwscconv がギャップを読めず止まった） / 回し直し(run2m): CONVERGED_METAL 5 反復、ギャップ none eV（LDA none） |
| [mp-1072660](https://next-gen.materialsproject.org/materials/mp-1072660) | V2As4 | 6 | 0（金属） | 回し直し | none | 5 月: NOGAP_METAL_IN_LDA_TOO（金属・半金属。gwscconv がギャップを読めず止まった） / 回し直し(run2m): CONVERGED_METAL 5 反復、ギャップ none eV（LDA none） |

## SUSPECT_GOOD

5 月は GOOD だが、QSGW80 < LDA または 2 反復目以降に 0.5 eV を超えて振動。回し直しで全部収束。18 物質。

**表 7**

| mpid | 組成 | 原子 | 採用のギャップ | どちらの値 | LDA | 注記 |
| --- | --- | --- | --- | --- | --- | --- |
| [mp-2894](https://next-gen.materialsproject.org/materials/mp-2894) | ScSnAu | 3 | 0.139484 | 回し直し | 0.150711 | 5 月は GOOD だが QSGW80 のギャップが LDA より小さい（履歴 0.1468, 0.1519, 0.1462） / 回し直し(run3): CONVERGED 4 反復、ギャップ 0.139484 eV（LDA 0.150711）、5 月との差 -0.01 eV |
| [mp-10423](https://next-gen.materialsproject.org/materials/mp-10423) | KAuC2 | 4 | 3.931921 | 回し直し | 1.932075 | 5 月は GOOD だが 2 反復目以降に 0.5 eV を超えて振動（履歴 3.5416, 3.7206, 3.8567, 3.4574, 2.4731, 4.4847, 4.4746, 3.8584, 3.9167, 3.9502） / 回し直し(run3): CONVERGED 5 反復、ギャップ 3.931921 eV（LDA 1.932075）、5 月との差 -0.02 eV |
| [mp-1006888](https://next-gen.materialsproject.org/materials/mp-1006888) | KYS2 | 4 | 4.019836 | 回し直し | 2.150904 | 5 月は GOOD だが 2 反復目以降に 0.5 eV を超えて振動（履歴 2.6872, 3.9244, 3.2604, 3.3863, 3.4772, 3.4245） / 回し直し(run3): CONVERGED 5 反復、ギャップ 4.019836 eV（LDA 2.150904）、5 月との差 +0.60 eV |
| [mp-1001786](https://next-gen.materialsproject.org/materials/mp-1001786) | LiScS2 | 4 | 3.043529 | 回し直し | 1.331540 | 5 月は GOOD だが 2 反復目以降に 0.5 eV を超えて振動（履歴 2.7887, 3.1778, 2.5089, 2.7130, 2.7803, 2.7531） / 回し直し(run3): CONVERGED 5 反復、ギャップ 3.043529 eV（LDA 1.331540）、5 月との差 +0.29 eV |
| [mp-22003](https://next-gen.materialsproject.org/materials/mp-22003) | NaN3 | 4 | 7.308540 | 回し直し | 3.909943 | 5 月は GOOD だが 2 反復目以降に 0.5 eV を超えて振動（履歴 6.9434, 5.3852, 6.2185, 5.8605, 5.9227, 5.8322） / 回し直し(run3): CONVERGED 5 反復、ギャップ 7.308540 eV（LDA 3.909943）、5 月との差 +1.48 eV |
| [mp-3056](https://next-gen.materialsproject.org/materials/mp-3056) | NaTlO2 | 4 | 1.981450 | 回し直し | 0.640949 | 5 月は GOOD だが 2 反復目以降に 0.5 eV を超えて振動（履歴 2.0278, 0.5805, 1.4607, 1.4096, 1.4647） / 回し直し(run3): CONVERGED 6 反復、ギャップ 1.981450 eV（LDA 0.640949）、5 月との差 +0.52 eV |
| [mp-4636](https://next-gen.materialsproject.org/materials/mp-4636) | ScCuO2 | 4 | 4.015221 | 回し直し | 2.260494 | 5 月は GOOD だが 2 反復目以降に 0.5 eV を超えて振動（履歴 3.3604, 4.7190, 4.0905, 4.1828, 4.1155） / 回し直し(run3): CONVERGED 6 反復、ギャップ 4.015221 eV（LDA 2.260494）、5 月との差 -0.10 eV |
| [mp-1066131](https://next-gen.materialsproject.org/materials/mp-1066131) | YTlS2 | 4 | 2.285675 | 回し直し | 1.374162 | 5 月は GOOD だが 2 反復目以降に 0.5 eV を超えて振動（履歴 1.8727, 2.4057, 3.0360, 2.6476, 2.5802, 2.5215, 2.5781） / 回し直し(run3): CONVERGED 5 反復、ギャップ 2.285675 eV（LDA 1.374162）、5 月との差 -0.29 eV |
| [mp-29585](https://next-gen.materialsproject.org/materials/mp-29585) | K4CdAs2 | 7 | 1.617516 | 回し直し | 0.499491 | 5 月は GOOD だが 2 反復目以降に 0.5 eV を超えて振動（履歴 1.4420, 0.4743, 1.0255, 1.1139, 1.1070） / 回し直し(run3): CONVERGED 5 反復、ギャップ 1.617516 eV（LDA 0.499491）、5 月との差 +0.51 eV |
| [mp-8753](https://next-gen.materialsproject.org/materials/mp-8753) | K4HgP2 | 7 | 1.931914 | 回し直し | 0.802632 | 5 月は GOOD だが 2 反復目以降に 0.5 eV を超えて振動（履歴 1.8568, 0.7993, 1.6384, 1.6377, 1.6410） / 回し直し(run3): CONVERGED 5 反復、ギャップ 1.931914 eV（LDA 0.802632）、5 月との差 +0.29 eV |
| [mp-1079707](https://next-gen.materialsproject.org/materials/mp-1079707) | Ca4O4 | 8 | 4.848429 | 回し直し | 2.128723 | 構造: Ca–Co 水酸化物の「(Ca(OH))-part」から H を除いたもの（Ca4O4）。c=16.75 Å の層、密度 2.00（岩塩型 CaO 3.34）。5 月は 9.97 → 3.31 → 4.78 eV と荒れた / 5 月は GOOD だが 2 反復目以降に 0.5 eV を超えて振動（履歴 9.9672, 3.3138, 4.7837, 4.7766, 4.7229） / 回し直し(run3): CONVERGED 6 反復、ギャップ 4.848429 eV（LDA 2.128723）、5 月との差 +0.13 eV |
| [mp-3829](https://next-gen.materialsproject.org/materials/mp-3829) | Cd2Sn2As4 | 8 | 0.051696 | 回し直し | 0.111058 | 5 月は GOOD だが QSGW80 のギャップが LDA より小さい（履歴 0.0595, 0.0628, 0.0628） / 回し直し(run3): CONVERGED 4 反復、ギャップ 0.051696 eV（LDA 0.111058）、5 月との差 -0.01 eV |
| [mp-9636](https://next-gen.materialsproject.org/materials/mp-9636) | CsSbF6 | 8 | 10.198779 | 回し直し | 5.227938 | 5 月は GOOD だが 2 反復目以降に 0.5 eV を超えて振動（履歴 9.5758, 10.1857, 9.4742, 9.5269, 9.5218） / 回し直し(run3): CONVERGED 6 反復、ギャップ 10.198779 eV（LDA 5.227938）、5 月との差 +0.68 eV |
| [mp-27999](https://next-gen.materialsproject.org/materials/mp-27999) | K5CuSb2 | 8 | 1.646245 | 回し直し | 0.635402 | 5 月は GOOD だが 2 反復目以降に 0.5 eV を超えて振動（履歴 1.6216, 2.5315, 1.9968, 1.9250, 1.8993） / 回し直し(run3): CONVERGED 5 反復、ギャップ 1.646245 eV（LDA 0.635402）、5 月との差 -0.25 eV |
| [mp-796276](https://next-gen.materialsproject.org/materials/mp-796276) | Mo2O6 | 8 | 0.981008 | 回し直し | 1.137048 | 5 月は GOOD だが QSGW80 のギャップが LDA より小さい（履歴 1.7149, 0.7988, 1.0927, 1.0956, 1.0694） / 回し直し(run3): CONVERGED 7 反復、ギャップ 0.981008 eV（LDA 1.137048）、5 月との差 -0.09 eV |
| [mp-1079088](https://next-gen.materialsproject.org/materials/mp-1079088) | NaAsF6 | 8 | 11.073740 | 回し直し | 5.472253 | 5 月は GOOD だが 2 反復目以降に 0.5 eV を超えて振動（履歴 11.7973, 10.4477, 11.1668, 11.1176, 11.1072） / 回し直し(run3): CONVERGED 6 反復、ギャップ 11.073740 eV（LDA 5.472253）、5 月との差 -0.03 eV |
| [mp-10474](https://next-gen.materialsproject.org/materials/mp-10474) | NaPF6 | 8 | 12.780387 | 回し直し | 7.204046 | 5 月は GOOD だが 2 反復目以降に 0.5 eV を超えて振動（履歴 12.0019, 14.5841, 13.3865, 12.9309, 12.7966, 12.7464, 12.6808, 12.7484） / 回し直し(run3): CONVERGED 6 反復、ギャップ 12.780387 eV（LDA 7.204046）、5 月との差 +0.03 eV |
| [mp-1078195](https://next-gen.materialsproject.org/materials/mp-1078195) | Si2I6 | 8 | 4.940172 | 回し直し | 2.779810 | 5 月は GOOD だが 2 反復目以降に 0.5 eV を超えて振動（履歴 4.9379, 4.1919, 5.1101, 5.0901, 5.1581） / 回し直し(run3): CONVERGED 5 反復、ギャップ 4.940172 eV（LDA 2.779810）、5 月との差 -0.22 eV |

## NOTCONV_MAY

5 月は 10 反復で収束しなかった。回し直しで全部収束。179 物質。

**表 8**

| mpid | 組成 | 原子 | 採用のギャップ | どちらの値 | LDA | 注記 |
| --- | --- | --- | --- | --- | --- | --- |
| [mp-169](https://next-gen.materialsproject.org/materials/mp-169) | C2 | 2 | 2.519014 | 回し直し | 1.663932 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1482 eV、履歴 2.3996, 3.5521, 2.9277, 2.7925, 2.9407） / 回し直し(run3): CONVERGED 4 反復、ギャップ 2.519014 eV（LDA 1.663932）、5 月との差 -0.42 eV |
| [mp-2605](https://next-gen.materialsproject.org/materials/mp-2605) | CaO | 2 | 6.601392 | 回し直し | 3.518064 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1450 eV、履歴 6.4469, 6.7180, 6.6771, 6.6035, 6.5321） / 回し直し(run3): CONVERGED 6 反復、ギャップ 6.601392 eV（LDA 3.518064）、5 月との差 +0.07 eV |
| [mp-682](https://next-gen.materialsproject.org/materials/mp-682) | NaF | 2 | 11.484239 | 回し直し | 6.428080 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1311 eV、履歴 11.2997, 11.6480, 11.5931, 11.4706, 11.4620） / 回し直し(run3): CONVERGED 6 反復、ギャップ 11.484239 eV（LDA 6.428080）、5 月との差 +0.02 eV |
| [mp-569639](https://next-gen.materialsproject.org/materials/mp-569639) | TlCl | 2 | 3.942578 | 回し直し | 2.261022 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.5346 eV、履歴 3.4055, 2.3861, 1.7880, 1.4565, 1.2534） / 回し直し(run3): CONVERGED 5 反復、ギャップ 3.942578 eV（LDA 2.261022）、5 月との差 +2.69 eV |
| [mp-1986](https://next-gen.materialsproject.org/materials/mp-1986) | ZnO | 2 | 3.037948 | 回し直し | 0.648771 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1008 eV、履歴 2.9589, 2.7614, 2.8962, 2.8394, 2.9402） / 回し直し(run3): CONVERGED 6 反復、ギャップ 3.037948 eV（LDA 0.648771）、5 月との差 +0.10 eV |
| [mp-1105](https://next-gen.materialsproject.org/materials/mp-1105) | BaO2 | 3 | 4.926768 | 回し直し | 2.129751 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1096 eV、履歴 4.7360, 5.1292, 5.1844, 5.0748, 5.1054） / 回し直し(run3): CONVERGED 6 反復、ギャップ 4.926768 eV（LDA 2.129751）、5 月との差 -0.18 eV |
| [mp-568690](https://next-gen.materialsproject.org/materials/mp-568690) | CdBr2 | 3 | 4.863640 | 回し直し | 2.601499 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.3571 eV、履歴 4.6782, 4.0737, 4.7155, 4.7898, 4.4326） / 回し直し(run3): CONVERGED 6 反復、ギャップ 4.863640 eV（LDA 2.601499）、5 月との差 +0.43 eV |
| [mp-22881](https://next-gen.materialsproject.org/materials/mp-22881) | CdCl2 | 3 | 5.994382 | 回し直し | 3.298712 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.3546 eV、履歴 4.5709, 5.6222, 5.7990, 5.4444, 5.5695） / 回し直し(run3): CONVERGED 6 反復、ギャップ 5.994382 eV（LDA 3.298712）、5 月との差 +0.42 eV |
| [mp-2278](https://next-gen.materialsproject.org/materials/mp-2278) | HgO2 | 3 | 1.765741 | 回し直し | 0.022309 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1711 eV、履歴 1.7291, 1.5795, 1.6899, 1.6899, 1.8610） / 回し直し(run3): CONVERGED 6 反復、ギャップ 1.765741 eV（LDA 0.022309）、5 月との差 -0.10 eV |
| [mp-25210](https://next-gen.materialsproject.org/materials/mp-25210) | NiO2 | 3 | 2.520079 | 回し直し | 0.916763 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1825 eV、履歴 4.4727, 2.2590, 2.1474, 2.3299, 2.1758） / 回し直し(run3): CONVERGED 6 反復、ギャップ 2.520079 eV（LDA 0.916763）、5 月との差 +0.34 eV |
| [mp-35925](https://next-gen.materialsproject.org/materials/mp-35925) | NiO2 | 3 | 3.818564 | 回し直し | 1.184401 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2662 eV、履歴 3.3566, 4.3529, 4.0668, 3.8651, 3.8006） / 回し直し(run3): CONVERGED 7 反復、ギャップ 3.818564 eV（LDA 1.184401）、5 月との差 +0.02 eV |
| [mp-22883](https://next-gen.materialsproject.org/materials/mp-22883) | PbI2 | 3 | 3.545461 | 回し直し | 2.234387 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2315 eV、履歴 3.5873, 4.4735, 3.9389, 3.8785, 3.7074） / 回し直し(run3): CONVERGED 5 反復、ギャップ 3.545461 eV（LDA 2.234387）、5 月との差 -0.16 eV |
| [mp-9813](https://next-gen.materialsproject.org/materials/mp-9813) | WS2 | 3 | 1.999899 | 回し直し | 1.314100 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1242 eV、履歴 1.3930, 1.3211, 1.4459, 1.3523, 1.3217） / 回し直し(run3): CONVERGED 4 反復、ギャップ 1.999899 eV（LDA 1.314100）、5 月との差 +0.68 eV |
| [mp-569960](https://next-gen.materialsproject.org/materials/mp-569960) | ZnBr2 | 3 | 5.545472 | 回し直し | 3.048474 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.4128 eV、履歴 5.1023, 5.7063, 5.3313, 5.6170, 5.7441） / 回し直し(run3): CONVERGED 6 反復、ギャップ 5.545472 eV（LDA 3.048474）、5 月との差 -0.20 eV |
| [mp-5770](https://next-gen.materialsproject.org/materials/mp-5770) | AgNO2 | 4 | 4.969543 | 回し直し | 1.701558 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1126 eV、履歴 1.6156, 1.0370, 1.1698, 1.0572, 1.0628） / 回し直し(run3): CONVERGED 6 反復、ギャップ 4.969543 eV（LDA 1.701558）、5 月との差 +3.91 eV |
| [mp-8039](https://next-gen.materialsproject.org/materials/mp-8039) | AlF3 | 4 | 12.837351 | 回し直し | 7.533142 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2283 eV、履歴 12.3337, 13.8761, 13.2156, 13.2866, 13.0583） / 回し直し(run3): CONVERGED 6 反復、ギャップ 12.837351 eV（LDA 7.533142）、5 月との差 -0.22 eV |
| [mp-1018098](https://next-gen.materialsproject.org/materials/mp-1018098) | Ba2BrN | 4 | 2.317010 | 回し直し | 1.013612 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2522 eV、履歴 2.9441, 2.3783, 2.2058, 2.6218, 2.9607, 3.1497, 2.4695, 2.4604, 2.2082, 2.4238） / 回し直し(run3): CONVERGED 5 反復、ギャップ 2.317010 eV（LDA 1.013612）、5 月との差 -0.11 eV |
| [mp-3915](https://next-gen.materialsproject.org/materials/mp-3915) | BaHgO2 | 4 | 4.623286 | 回し直し | 2.244048 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2914 eV、履歴 4.3522, 3.9422, 3.7591, 4.0073, 4.0505） / 回し直し(run3): CONVERGED 5 反復、ギャップ 4.623286 eV（LDA 2.244048）、5 月との差 +0.57 eV |
| [mp-569304](https://next-gen.materialsproject.org/materials/mp-569304) | C4 | 4 | 3.677225 | 回し直し | 1.848254 | 構造: 硝酸を挿入した黒鉛（Graphite, nitrated）から N・O を除いたもの。密度 1.38（黒鉛 2.26）、層間が開いたまま / 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1754 eV、履歴 3.9831, 3.6407, 3.6924, 3.8422, 3.8678） / 回し直し(run3): CONVERGED 5 反復、ギャップ 3.677225 eV（LDA 1.848254）、5 月との差 -0.19 eV |
| [mp-22936](https://next-gen.materialsproject.org/materials/mp-22936) | Ca2NCl | 4 | 4.152573 | 回し直し | 2.222162 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1008 eV、履歴 2.6339, 3.7728, 3.1456, 3.2463, 3.2232） / 回し直し(run3): CONVERGED 5 反復、ギャップ 4.152573 eV（LDA 2.222162）、5 月との差 +0.93 eV |
| [mp-23040](https://next-gen.materialsproject.org/materials/mp-23040) | Ca2PI | 4 | 2.937550 | 回し直し | 1.547321 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2151 eV、履歴 2.8001, 2.7499, 2.4916, 2.4856, 2.2765） / 回し直し(run3): CONVERGED 5 反復、ギャップ 2.937550 eV（LDA 1.547321）、5 月との差 +0.66 eV |
| [mp-4124](https://next-gen.materialsproject.org/materials/mp-4124) | CaCN2 | 4 | 5.475552 | 回し直し | 3.255472 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1802 eV、履歴 5.1758, 6.4096, 5.7632, 5.8667, 5.6865） / 回し直し(run3): CONVERGED 5 反復、ギャップ 5.475552 eV（LDA 3.255472）、5 月との差 -0.21 eV |
| [mp-7041](https://next-gen.materialsproject.org/materials/mp-7041) | CaHgO2 | 4 | 4.435741 | 回し直し | 2.134180 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.3892 eV、履歴 7.4098, 4.1163, 4.1785, 4.5441, 4.5677） / 回し直し(run3): CONVERGED 5 反復、ギャップ 4.435741 eV（LDA 2.134180）、5 月との差 -0.13 eV |
| [mp-561203](https://next-gen.materialsproject.org/materials/mp-561203) | F4 | 4 | 9.928535 | 回し直し | 2.428449 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2644 eV、履歴 9.5504, 10.3164, 10.1199, 9.9295, 9.8555） / 回し直し(run3): CONVERGED 6 反復、ギャップ 9.928535 eV（LDA 2.428449）、5 月との差 +0.07 eV |
| [mp-9889](https://next-gen.materialsproject.org/materials/mp-9889) | Ga2S2 | 4 | 2.980352 | 回し直し | 1.439901 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1315 eV、履歴 2.8969, 2.8148, 2.6709, 2.8024, 2.7951） / 回し直し(run3): CONVERGED 5 反復、ギャップ 2.980352 eV（LDA 1.439901）、5 月との差 +0.19 eV |
| [mp-632296](https://next-gen.materialsproject.org/materials/mp-632296) | H2F2 | 4 | 12.253660 | 回し直し | 6.392350 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2143 eV、履歴 13.2818, 11.9113, 12.1742, 12.1383, 12.3526） / 回し直し(run3): CONVERGED 6 反復、ギャップ 12.253660 eV（LDA 6.392350）、5 月との差 -0.10 eV |
| [mp-706](https://next-gen.materialsproject.org/materials/mp-706) | Hg2F2 | 4 | 2.891362 | 回し直し | 1.397021 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2635 eV、履歴 3.2244, 2.7793, 2.7450, 2.9804, 3.0086） / 回し直し(run3): CONVERGED 5 反復、ギャップ 2.891362 eV（LDA 1.397021）、5 月との差 -0.12 eV |
| [mp-22660](https://next-gen.materialsproject.org/materials/mp-22660) | InAgO2 | 4 | 1.975147 | 回し直し | 0.277576 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1362 eV、履歴 1.4049, 1.6665, 1.7515, 1.6154, 1.6161） / 回し直し(run3): CONVERGED 6 反復、ギャップ 1.975147 eV（LDA 0.277576）、5 月との差 +0.36 eV |
| [mp-20930](https://next-gen.materialsproject.org/materials/mp-20930) | InCuO2 | 4 | 1.955047 | 回し直し | 0.394554 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1314 eV、履歴 2.8165, 1.9590, 1.7192, 1.8505, 1.8051） / 回し直し(run3): CONVERGED 6 反復、ギャップ 1.955047 eV（LDA 0.394554）、5 月との差 +0.15 eV |
| [mp-34857](https://next-gen.materialsproject.org/materials/mp-34857) | KNO2 | 4 | 7.170460 | 回し直し | 2.372951 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 1.0537 eV、履歴 7.0565, 12.3139, 9.6802, 9.0982, 8.6264） / 回し直し(run3): CONVERGED 6 反復、ギャップ 7.170460 eV（LDA 2.372951）、5 月との差 -1.46 eV |
| [mp-546787](https://next-gen.materialsproject.org/materials/mp-546787) | KNO2 | 4 | 7.054692 | 回し直し | 2.353343 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2912 eV、履歴 6.9290, 7.2673, 6.9773, 6.8072, 6.6861） / 回し直し(run3): CONVERGED 6 反復、ギャップ 7.054692 eV（LDA 2.353343）、5 月との差 +0.37 eV |
| [mp-8175](https://next-gen.materialsproject.org/materials/mp-8175) | KTlO2 | 4 | 2.228306 | 回し直し | 0.893942 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2057 eV、履歴 2.1931, 1.8986, 2.1137, 2.2517, 2.3194） / 回し直し(run3): CONVERGED 5 反復、ギャップ 2.228306 eV（LDA 0.893942）、5 月との差 -0.09 eV |
| [mp-8001](https://next-gen.materialsproject.org/materials/mp-8001) | LiAlO2 | 4 | 9.837572 | 回し直し | 6.632788 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.8754 eV、履歴 9.0744, 9.6297, 9.2447, 9.1817, 8.3693） / 回し直し(run3): CONVERGED 5 反復、ギャップ 9.837572 eV（LDA 6.632788）、5 月との差 +1.47 eV |
| [mp-22526](https://next-gen.materialsproject.org/materials/mp-22526) | LiCoO2 | 4 | 3.194087 | 回し直し | 0.922027 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1918 eV、履歴 2.5707, 2.8438, 2.9462, 2.7544, 2.8797） / 回し直し(run3): CONVERGED 9 反復、ギャップ 3.194087 eV（LDA 0.922027）、5 月との差 +0.31 eV |
| [mp-15794](https://next-gen.materialsproject.org/materials/mp-15794) | LiYSe2 | 4 | 2.699236 | 回し直し | 1.403927 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2229 eV、履歴 2.9652, 2.5832, 2.5383, 2.7612, 2.7103） / 回し直し(run3): CONVERGED 5 反復、ギャップ 2.699236 eV（LDA 1.403927）、5 月との差 -0.01 eV |
| [mp-549706](https://next-gen.materialsproject.org/materials/mp-549706) | Mg2O2 | 4 | 6.544901 | 回し直し | 3.401128 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2320 eV、履歴 6.8317, 6.9330, 6.8118, 6.6700, 6.5797） / 回し直し(run3): CONVERGED 6 反復、ギャップ 6.544901 eV（LDA 3.401128）、5 月との差 -0.04 eV |
| [mp-9166](https://next-gen.materialsproject.org/materials/mp-9166) | MgCN2 | 4 | 6.284258 | 回し直し | 3.600962 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.6720 eV、履歴 6.9193, 6.5422, 7.2905, 6.6185, 6.6353） / 回し直し(run3): CONVERGED 5 反復、ギャップ 6.284258 eV（LDA 3.600962）、5 月との差 -0.35 eV |
| [mp-1017514](https://next-gen.materialsproject.org/materials/mp-1017514) | MgSiN2 | 4 | 6.065007 | 回し直し | 4.078558 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1270 eV、履歴 5.1661, 5.7968, 5.7694, 5.6424, 5.7593） / 回し直し(run3): CONVERGED 5 反復、ギャップ 6.065007 eV（LDA 4.078558）、5 月との差 +0.31 eV |
| [mp-7748](https://next-gen.materialsproject.org/materials/mp-7748) | NaAlO2 | 4 | 7.847294 | 回し直し | 4.661912 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.4383 eV、履歴 8.8062, 6.8088, 7.2344, 6.8435, 7.2818） / 回し直し(run3): CONVERGED 6 反復、ギャップ 7.847294 eV（LDA 4.661912）、5 月との差 +0.57 eV |
| [mp-546500](https://next-gen.materialsproject.org/materials/mp-546500) | NaCNO | 4 | 7.667180 | 回し直し | 4.202975 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 1.0322 eV、履歴 9.6311, 7.2993, 7.3406, 8.3729, 8.3547） / 回し直し(run3): CONVERGED 5 反復、ギャップ 7.667180 eV（LDA 4.202975）、5 月との差 -0.69 eV |
| [mp-27837](https://next-gen.materialsproject.org/materials/mp-27837) | NaHF2 | 4 | 12.057696 | 回し直し | 6.758239 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.4175 eV、履歴 11.9130, 11.0026, 10.8636, 11.1995, 11.2812） / 回し直し(run3): CONVERGED 6 反復、ギャップ 12.057696 eV（LDA 6.758239）、5 月との差 +0.78 eV |
| [mp-5175](https://next-gen.materialsproject.org/materials/mp-5175) | NaInO2 | 4 | 4.110343 | 回し直し | 2.011135 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.7292 eV、履歴 4.1603, 7.8133, 5.7953, 5.2157, 5.0661） / 回し直し(run3): CONVERGED 6 反復、ギャップ 4.110343 eV（LDA 2.011135）、5 月との差 -0.96 eV |
| [mp-22473](https://next-gen.materialsproject.org/materials/mp-22473) | NaInSe2 | 4 | 2.080008 | 回し直し | 0.884731 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1888 eV、履歴 2.2969, 2.1373, 2.3964, 2.2149, 2.2076） / 回し直し(run3): CONVERGED 5 反復、ギャップ 2.080008 eV（LDA 0.884731）、5 月との差 -0.13 eV |
| [mp-8830](https://next-gen.materialsproject.org/materials/mp-8830) | NaRhO2 | 4 | 3.321614 | 回し直し | 1.193512 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 5.6012 eV、履歴 6.4744, 4.0586, 3.2892, 5.7309, 0.1297） / 回し直し(run3): CONVERGED 6 反復、ギャップ 3.321614 eV（LDA 1.193512）、5 月との差 +3.19 eV |
| [mp-999460](https://next-gen.materialsproject.org/materials/mp-999460) | NaScS2 | 4 | 3.412194 | 回し直し | 1.508017 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2816 eV、履歴 3.6025, 2.3776, 2.6084, 2.6808, 2.8900） / 回し直し(run3): CONVERGED 6 反復、ギャップ 3.412194 eV（LDA 1.508017）、5 月との差 +0.52 eV |
| [mp-30041](https://next-gen.materialsproject.org/materials/mp-30041) | RbBiS2 | 4 | 2.213998 | 回し直し | 1.241437 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.7592 eV、履歴 2.5212, 6.2073, 4.0482, 3.5214, 3.2890） / 回し直し(run3): CONVERGED 5 反復、ギャップ 2.213998 eV（LDA 1.241437）、5 月との差 -1.08 eV |
| [mp-559092](https://next-gen.materialsproject.org/materials/mp-559092) | ScF3 | 4 | 11.882829 | 回し直し | 5.908148 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.8750 eV、履歴 11.8695, 17.7027, 14.7861, 14.2379, 13.9110） / 回し直し(run3): CONVERGED 6 反復、ギャップ 11.882829 eV（LDA 5.908148）、5 月との差 -2.03 eV |
| [mp-554134](https://next-gen.materialsproject.org/materials/mp-554134) | Sn2S2 | 4 | 3.997399 | 回し直し | 1.510139 | 構造: Gd–Sn–Nb–S の非整合層状化合物の SnS 層。a=4.13、c=22.4 Å、1 原子 95.7 Å³（SnS 約 24）、配位 4 / 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1797 eV、履歴 1.3981, 4.3999, 3.8971, 3.7174, 3.7793） / 回し直し(run3): CONVERGED 6 反復、ギャップ 3.997399 eV（LDA 1.510139）、5 月との差 +0.22 eV |
| [mp-23056](https://next-gen.materialsproject.org/materials/mp-23056) | Sr2BrN | 4 | 3.280951 | 回し直し | 1.851198 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.6348 eV、履歴 3.1119, 2.7926, 2.7474, 3.1918, 2.5570） / 回し直し(run3): CONVERGED 5 反復、ギャップ 3.280951 eV（LDA 1.851198）、5 月との差 +0.72 eV |
| [mp-6972](https://next-gen.materialsproject.org/materials/mp-6972) | YCuO2 | 4 | 4.370518 | 回し直し | 2.701063 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1619 eV、履歴 3.4850, 4.4042, 4.1165, 3.9546, 3.9658） / 回し直し(run3): CONVERGED 6 反復、ギャップ 4.370518 eV（LDA 2.701063）、5 月との差 +0.40 eV |
| [mp-9307](https://next-gen.materialsproject.org/materials/mp-9307) | Ba2ZnN2 | 5 | 1.535281 | 回し直し | 0.462578 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1027 eV、履歴 1.3642, 2.7823, 1.8641, 1.7912, 1.7614） / 回し直し(run3): CONVERGED 5 反復、ギャップ 1.535281 eV（LDA 0.462578）、5 月との差 -0.23 eV |
| [mp-10250](https://next-gen.materialsproject.org/materials/mp-10250) | BaLiF3 | 5 | 10.862510 | 回し直し | 6.619731 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.3028 eV、履歴 11.5070, 11.1200, 10.8117, 11.1145, 11.0945） / 回し直し(run3): CONVERGED 6 反復、ギャップ 10.862510 eV（LDA 6.619731）、5 月との差 -0.23 eV |
| [mp-644497](https://next-gen.materialsproject.org/materials/mp-644497) | BaTiO3 | 5 | 2.705998 | 回し直し | 0.322921 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.4244 eV、履歴 2.2493, 2.6584, 3.0826, 2.7961, 2.6582） / 回し直し(run3): CONVERGED 6 反復、ギャップ 2.705998 eV（LDA 0.322921）、5 月との差 +0.05 eV |
| [mp-27962](https://next-gen.materialsproject.org/materials/mp-27962) | CI4 | 5 | 3.277974 | 回し直し | 1.022817 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1763 eV、履歴 3.2155, 3.4407, 3.5269, 3.3775, 3.3506） / 回し直し(run3): CONVERGED 6 反復、ギャップ 3.277974 eV（LDA 1.022817）、5 月との差 -0.07 eV |
| [mp-8818](https://next-gen.materialsproject.org/materials/mp-8818) | Ca2ZnN2 | 5 | 1.971196 | 回し直し | 0.704174 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2865 eV、履歴 2.5443, 2.0687, 1.7297, 2.0154, 2.0162） / 回し直し(run3): CONVERGED 5 反復、ギャップ 1.971196 eV（LDA 0.704174）、5 月との差 -0.04 eV |
| [mp-7232](https://next-gen.materialsproject.org/materials/mp-7232) | Cs2HgO2 | 5 | 4.354080 | 回し直し | 1.993102 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1134 eV、履歴 4.2227, 4.0707, 4.1772, 4.2906, 4.2863） / 回し直し(run3): CONVERGED 6 反復、ギャップ 4.354080 eV（LDA 1.993102）、5 月との差 +0.07 eV |
| [mp-28873](https://next-gen.materialsproject.org/materials/mp-28873) | CsBrO3 | 5 | 7.058765 | 回し直し | 4.023388 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 1.1443 eV、履歴 3.4566, 0.3840, 5.7697, 6.1218, 4.9775） / 回し直し(run3): CONVERGED 6 反復、ギャップ 7.058765 eV（LDA 4.023388）、5 月との差 +2.08 eV |
| [mp-9816](https://next-gen.materialsproject.org/materials/mp-9816) | GeF4 | 5 | 10.539022 | 回し直し | 5.105805 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1756 eV、履歴 11.1024, 10.3504, 10.4979, 10.6734, 10.6386） / 回し直し(run3): CONVERGED 6 反復、ギャップ 10.539022 eV（LDA 5.105805）、5 月との差 -0.10 eV |
| [mp-5860](https://next-gen.materialsproject.org/materials/mp-5860) | K2HgO2 | 5 | 4.494736 | 回し直し | 2.048907 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.3196 eV、履歴 4.0425, 4.1540, 4.2867, 3.9781, 3.9671） / 回し直し(run3): CONVERGED 5 反復、ギャップ 4.494736 eV（LDA 2.048907）、5 月との差 +0.53 eV |
| [mp-4950](https://next-gen.materialsproject.org/materials/mp-4950) | KCaF3 | 5 | 10.414511 | 回し直し | 5.874176 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.4472 eV、履歴 17.2325, 11.1734, 14.4786, 14.1452, 14.0314） / 回し直し(run3): CONVERGED 6 反復、ギャップ 10.414511 eV（LDA 5.874176）、5 月との差 -3.62 eV |
| [mp-550506](https://next-gen.materialsproject.org/materials/mp-550506) | KClO3 | 5 | 8.762410 | 回し直し | 5.054930 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1575 eV、履歴 8.9647, 10.0113, 9.4169, 9.3030, 9.2594） / 回し直し(run3): CONVERGED 6 反復、ギャップ 8.762410 eV（LDA 5.054930）、5 月との差 -0.50 eV |
| [mp-3448](https://next-gen.materialsproject.org/materials/mp-3448) | KMgF3 | 5 | 12.004626 | 回し直し | 7.119649 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1896 eV、履歴 11.9008, 11.6303, 11.8311, 11.9592, 12.0207） / 回し直し(run3): CONVERGED 6 反復、ギャップ 12.004626 eV（LDA 7.119649）、5 月との差 -0.02 eV |
| [mp-19183](https://next-gen.materialsproject.org/materials/mp-19183) | Li2NiO2 | 5 | 2.783714 | 回し直し | 0.247265 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1617 eV、履歴 3.0354, 4.2034, 3.8744, 3.9458, 3.7841） / 回し直し(run3): CONVERGED 6 反復、ギャップ 2.783714 eV（LDA 0.247265）、5 月との差 -1.00 eV |
| [mp-28593](https://next-gen.materialsproject.org/materials/mp-28593) | Li3BrO | 5 | 7.160383 | 回し直し | 4.518465 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1752 eV、履歴 7.8578, 7.2093, 7.1666, 6.9913, 7.0322） / 回し直し(run3): CONVERGED 5 反復、ギャップ 7.160383 eV（LDA 4.518465）、5 月との差 +0.13 eV |
| [mp-30247](https://next-gen.materialsproject.org/materials/mp-30247) | MgH2O2 | 5 | 8.151312 | 回し直し | 4.319961 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2228 eV、履歴 8.4924, 8.1004, 7.8063, 8.0291, 8.0207） / 回し直し(run3): CONVERGED 6 反復、ギャップ 8.151312 eV（LDA 4.319961）、5 月との差 +0.13 eV |
| [mp-541989](https://next-gen.materialsproject.org/materials/mp-541989) | Na2CN2 | 5 | 5.571464 | 回し直し | 2.918104 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.3925 eV、履歴 6.0655, 5.6459, 5.0588, 5.4039, 5.4514） / 回し直し(run3): CONVERGED 5 反復、ギャップ 5.571464 eV（LDA 2.918104）、5 月との差 +0.12 eV |
| [mp-4961](https://next-gen.materialsproject.org/materials/mp-4961) | Na2HgO2 | 5 | 3.794084 | 回し直し | 1.650900 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.7717 eV、履歴 3.6948, 7.9118, 5.7284, 5.0361, 4.9567） / 回し直し(run3): CONVERGED 6 反復、ギャップ 3.794084 eV（LDA 1.650900）、5 月との差 -1.16 eV |
| [mp-560399](https://next-gen.materialsproject.org/materials/mp-560399) | NaMgF3 | 5 | 10.718812 | 回し直し | 5.696566 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.4144 eV、履歴 10.2969, 12.7744, 11.5390, 11.4363, 11.1247） / 回し直し(run3): CONVERGED 6 反復、ギャップ 10.718812 eV（LDA 5.696566）、5 月との差 -0.41 eV |
| [mp-551905](https://next-gen.materialsproject.org/materials/mp-551905) | OsO4 | 5 | 6.327797 | 回し直し | 3.168967 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2341 eV、履歴 12.1155, 6.5198, 8.7246, 8.7129, 8.9470） / 回し直し(run3): CONVERGED 6 反復、ギャップ 6.327797 eV（LDA 3.168967）、5 月との差 -2.62 eV |
| [mp-341](https://next-gen.materialsproject.org/materials/mp-341) | PbF4 | 5 | 5.132519 | 回し直し | 1.843484 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2081 eV、履歴 5.5884, 5.1507, 5.1526, 5.3607, 5.3602） / 回し直し(run3): CONVERGED 6 反復、ギャップ 5.132519 eV（LDA 1.843484）、5 月との差 -0.23 eV |
| [mp-5072](https://next-gen.materialsproject.org/materials/mp-5072) | Rb2HgO2 | 5 | 4.462417 | 回し直し | 1.970998 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1798 eV、履歴 3.0524, 3.8463, 4.0817, 3.9019, 3.9800） / 回し直し(run3): CONVERGED 5 反復、ギャップ 4.462417 eV（LDA 1.970998）、5 月との差 +0.48 eV |
| [mp-6951](https://next-gen.materialsproject.org/materials/mp-6951) | RbCdF3 | 5 | 7.039320 | 回し直し | 2.997867 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.3508 eV、履歴 6.8906, 7.4553, 7.7575, 7.5648, 7.4067） / 回し直し(run3): CONVERGED 6 反復、ギャップ 7.039320 eV（LDA 2.997867）、5 月との差 -0.37 eV |
| [mp-7849](https://next-gen.materialsproject.org/materials/mp-7849) | AlAsO4 | 6 | 8.284399 | 回し直し | 4.437228 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.3531 eV、履歴 8.0778, 8.2783, 8.3679, 8.0973, 8.0148） / 回し直し(run3): CONVERGED 6 反復、ギャップ 8.284399 eV（LDA 4.437228）、5 月との差 +0.27 eV |
| [mp-27743](https://next-gen.materialsproject.org/materials/mp-27743) | BiF5 | 6 | 5.546530 | 回し直し | 1.788071 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1233 eV、履歴 6.0730, 5.6084, 5.5242, 5.4010, 5.4394） / 回し直し(run3): CONVERGED 6 反復、ギャップ 5.546530 eV（LDA 1.788071）、5 月との差 +0.11 eV |
| [mp-1077316](https://next-gen.materialsproject.org/materials/mp-1077316) | C2O4 | 6 | 11.202313 | 回し直し | 7.454761 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1030 eV、履歴 11.2827, 11.4816, 11.3461, 11.3586, 11.4490） / 回し直し(run3): CONVERGED 6 反復、ギャップ 11.202313 eV（LDA 7.454761）、5 月との差 -0.25 eV |
| [mp-644607](https://next-gen.materialsproject.org/materials/mp-644607) | C2O4 | 6 | 10.276498 | 回し直し | 6.223344 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1693 eV、履歴 10.5616, 10.1943, 10.0932, 10.2625, 10.2588） / 回し直し(run3): CONVERGED 5 反復、ギャップ 10.276498 eV（LDA 6.223344）、5 月との差 +0.02 eV |
| [mp-571166](https://next-gen.materialsproject.org/materials/mp-571166) | Ca2Br4 | 6 | 7.075557 | 回し直し | 4.211722 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1656 eV、履歴 7.6093, 7.0396, 6.8220, 6.9876, 6.9826） / 回し直し(run3): CONVERGED 6 反復、ギャップ 7.075557 eV（LDA 4.211722）、5 月との差 +0.09 eV |
| [mp-22904](https://next-gen.materialsproject.org/materials/mp-22904) | Ca2Cl4 | 6 | 8.374065 | 回し直し | 5.055033 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2608 eV、履歴 8.9133, 8.2358, 8.2011, 8.4619, 8.4455） / 回し直し(run3): CONVERGED 6 反復、ギャップ 8.374065 eV（LDA 5.055033）、5 月との差 -0.07 eV |
| [mp-569272](https://next-gen.materialsproject.org/materials/mp-569272) | Cs4Se2 | 6 | 3.494830 | 回し直し | 1.502293 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1902 eV、履歴 3.9961, 3.6518, 3.6229, 3.7969, 3.8131） / 回し直し(run3): CONVERGED 5 反復、ギャップ 3.494830 eV（LDA 1.502293）、5 月との差 -0.32 eV |
| [mp-27716](https://next-gen.materialsproject.org/materials/mp-27716) | K2Tl2O2 | 6 | 2.327173 | 回し直し | 1.120718 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.8465 eV、履歴 9.1049, 2.5409, 3.4137, 2.9419, 2.5672） / 回し直し(run3): CONVERGED 5 反復、ギャップ 2.327173 eV（LDA 1.120718）、5 月との差 -0.24 eV |
| [mp-23856](https://next-gen.materialsproject.org/materials/mp-23856) | Li2H2O2 | 6 | 7.880123 | 回し直し | 4.018490 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1714 eV、履歴 8.1984, 7.9139, 7.8655, 7.7302, 7.9016） / 回し直し(run3): CONVERGED 6 反復、ギャップ 7.880123 eV（LDA 4.018490）、5 月との差 -0.02 eV |
| [mp-1249](https://next-gen.materialsproject.org/materials/mp-1249) | Mg2F4 | 6 | 12.150571 | 回し直し | 6.926135 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1067 eV、履歴 11.6906, 11.9279, 11.9256, 11.8189, 11.8356） / 回し直し(run3): CONVERGED 6 反復、ギャップ 12.150571 eV（LDA 6.926135）、5 月との差 +0.31 eV |
| [mp-30053](https://next-gen.materialsproject.org/materials/mp-30053) | Na2C2N2 | 6 | 9.099201 | 回し直し | 4.878076 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.4198 eV、履歴 9.1843, 9.0948, 9.0496, 9.2476, 8.8278） / 回し直し(run3): CONVERGED 6 反復、ギャップ 9.099201 eV（LDA 4.878076）、5 月との差 +0.27 eV |
| [mp-23891](https://next-gen.materialsproject.org/materials/mp-23891) | Na2H2O2 | 6 | 6.787306 | 回し直し | 2.892917 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1490 eV、履歴 6.9672, 6.8267, 6.7069, 6.8560, 6.8458） / 回し直し(run3): CONVERGED 6 反復、ギャップ 6.787306 eV（LDA 2.892917）、5 月との差 -0.06 eV |
| [mp-23940](https://next-gen.materialsproject.org/materials/mp-23940) | Na2H2O2 | 6 | 6.813937 | 回し直し | 2.924744 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1741 eV、履歴 6.4612, 6.5393, 6.4005, 6.5407, 6.5747） / 回し直し(run3): CONVERGED 6 反復、ギャップ 6.813937 eV（LDA 2.924744）、5 月との差 +0.24 eV |
| [mp-7723](https://next-gen.materialsproject.org/materials/mp-7723) | NaAlF4 | 6 | 10.604024 | 回し直し | 5.535568 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1967 eV、履歴 10.1825, 11.5519, 10.9620, 11.0450, 10.8483） / 回し直し(run3): CONVERGED 6 反復、ギャップ 10.604024 eV（LDA 5.535568）、5 月との差 -0.24 eV |
| [mp-546118](https://next-gen.materialsproject.org/materials/mp-546118) | NaClO4 | 6 | 7.898071 | 回し直し | 4.265400 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2023 eV、履歴 7.6997, 8.0070, 7.9248, 7.7866, 7.7225） / 回し直し(run3): CONVERGED 5 反復、ギャップ 7.898071 eV（LDA 4.265400）、5 月との差 +0.18 eV |
| [mp-626620](https://next-gen.materialsproject.org/materials/mp-626620) | Rb2H2O2 | 6 | 6.663139 | 回し直し | 3.409712 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 1.3001 eV、履歴 6.7731, 6.5975, 6.9117, 5.6116, 6.3739） / 回し直し(run3): CONVERGED 6 反復、ギャップ 6.663139 eV（LDA 3.409712）、5 月との差 +0.29 eV |
| [mp-5479](https://next-gen.materialsproject.org/materials/mp-5479) | RbAlF4 | 6 | 11.667692 | 回し直し | 6.889897 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2945 eV、履歴 11.2627, 12.5530, 12.3678, 12.1039, 12.0733） / 回し直し(run3): CONVERGED 6 反復、ギャップ 11.667692 eV（LDA 6.889897）、5 月との差 -0.41 eV |
| [mp-27726](https://next-gen.materialsproject.org/materials/mp-27726) | S2O4 | 6 | 6.659674 | 回し直し | 2.598452 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1368 eV、履歴 7.1482, 6.5806, 6.8169, 6.8646, 6.9537） / 回し直し(run3): CONVERGED 6 反復、ギャップ 6.659674 eV（LDA 2.598452）、5 月との差 -0.29 eV |
| [mp-7](https://next-gen.materialsproject.org/materials/mp-7) | S6 | 6 | 4.011625 | 回し直し | 1.971537 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 1.1641 eV、履歴 3.7542, 5.0870, 5.7420, 4.8792, 4.5779） / 回し直し(run3): CONVERGED 6 反復、ギャップ 4.011625 eV（LDA 1.971537）、5 月との差 -0.57 eV |
| [mp-1071820](https://next-gen.materialsproject.org/materials/mp-1071820) | Si2O4 | 6 | 8.110389 | 回し直し | 4.447308 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.4987 eV、履歴 7.2799, 7.7638, 7.0240, 6.6353, 5.7265, 5.9130, 7.1338, 7.2966, 7.6517, 7.7954） / 回し直し(run3): CONVERGED 6 反復、ギャップ 8.110389 eV（LDA 4.447308）、5 月との差 +0.32 eV |
| [mp-546794](https://next-gen.materialsproject.org/materials/mp-546794) | Si2O4 | 6 | 9.435870 | 回し直し | 5.695930 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2164 eV、履歴 9.5056, 9.3843, 9.3260, 9.2195, 9.4359） / 回し直し(run3): CONVERGED 5 反復、ギャップ 9.435870 eV（LDA 5.695930）、5 月との差 -0.00 eV |
| [mp-1873](https://next-gen.materialsproject.org/materials/mp-1873) | Zn2F4 | 6 | 8.052167 | 回し直し | 3.431755 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2458 eV、履歴 8.2590, 9.1949, 8.7264, 8.5968, 8.4806） / 回し直し(run3): CONVERGED 6 反復、ギャップ 8.052167 eV（LDA 3.431755）、5 月との差 -0.43 eV |
| [mp-8335](https://next-gen.materialsproject.org/materials/mp-8335) | Ba2ZrO4 | 7 | 5.395929 | 回し直し | 2.941367 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1285 eV、履歴 5.2045, 5.7406, 5.4663, 5.3378, 5.3527） / 回し直し(run3): CONVERGED 6 反復、ギャップ 5.395929 eV（LDA 2.941367）、5 月との差 +0.04 eV |
| [mp-570572](https://next-gen.materialsproject.org/materials/mp-570572) | C3N4 | 7 | 3.748565 | 回し直し | 1.300089 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.3132 eV、履歴 3.3068, 3.6286, 3.0740, 2.0808, 5.4038, 4.8415, 3.6905, 3.6390, 3.9521, 3.9165） / 回し直し(run3): CONVERGED 5 反復、ギャップ 3.748565 eV（LDA 1.300089）、5 月との差 -0.17 eV |
| [mp-8682](https://next-gen.materialsproject.org/materials/mp-8682) | Ca2SiO4 | 7 | 7.190660 | 回し直し | 4.220324 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.3155 eV、履歴 7.6446, 7.4103, 7.1680, 7.0446, 7.3601） / 回し直し(run3): CONVERGED 6 反復、ギャップ 7.190660 eV（LDA 4.220324）、5 月との差 -0.17 eV |
| [mp-1077904](https://next-gen.materialsproject.org/materials/mp-1077904) | CdPb2Cl2O2 | 7 | 3.954747 | 回し直し | 2.407234 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2395 eV、履歴 3.6826, 3.7988, 3.9826, 4.2221, 4.0409） / 回し直し(run3): CONVERGED 5 反復、ギャップ 3.954747 eV（LDA 2.407234）、5 月との差 -0.09 eV |
| [mp-15157](https://next-gen.materialsproject.org/materials/mp-15157) | Cs2CaF4 | 7 | 9.989961 | 回し直し | 6.004742 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1832 eV、履歴 8.9583, 6.3363, 9.7501, 9.6822, 9.8654） / 回し直し(run3): CONVERGED 6 反復、ギャップ 9.989961 eV（LDA 6.004742）、5 月との差 +0.12 eV |
| [mp-697133](https://next-gen.materialsproject.org/materials/mp-697133) | Cs2CaH4 | 7 | 5.644610 | 回し直し | 2.665492 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.8791 eV、履歴 5.4038, 5.9691, 6.0047, 5.2414, 5.1255） / 回し直し(run3): CONVERGED 6 反復、ギャップ 5.644610 eV（LDA 2.665492）、5 月との差 +0.52 eV |
| [mp-8963](https://next-gen.materialsproject.org/materials/mp-8963) | Cs2HgF4 | 7 | 5.109075 | 回し直し | 1.903564 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1353 eV、履歴 4.5373, 4.6299, 4.1729, 4.2905, 4.3082） / 回し直し(run3): CONVERGED 6 反復、ギャップ 5.109075 eV（LDA 1.903564）、5 月との差 +0.80 eV |
| [mp-31212](https://next-gen.materialsproject.org/materials/mp-31212) | K2MgF4 | 7 | 11.148002 | 回し直し | 6.454241 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 1.2301 eV、履歴 5.2535, 8.5796, 9.2650, 10.4952, 9.6980） / 回し直し(run3): CONVERGED 6 反復、ギャップ 11.148002 eV（LDA 6.454241）、5 月との差 +1.45 eV |
| [mp-22934](https://next-gen.materialsproject.org/materials/mp-22934) | K2PtCl4 | 7 | 4.619134 | 回し直し | 1.553973 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2468 eV、履歴 5.1806, 4.6352, 4.5730, 4.7387, 4.8197） / 回し直し(run3): CONVERGED 6 反復、ギャップ 4.619134 eV（LDA 1.553973）、5 月との差 -0.20 eV |
| [mp-9872](https://next-gen.materialsproject.org/materials/mp-9872) | K4BeP2 | 7 | 2.252272 | 回し直し | 0.952354 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1773 eV、履歴 1.6527, 1.3225, 1.4751, 1.6056, 1.6524） / 回し直し(run3): CONVERGED 5 反復、ギャップ 2.252272 eV（LDA 0.952354）、5 月との差 +0.60 eV |
| [mp-28627](https://next-gen.materialsproject.org/materials/mp-28627) | K4Br2O | 7 | 3.403743 | 回し直し | 0.951647 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1285 eV、履歴 3.6526, 3.2468, 3.2966, 3.4251, 3.3871） / 回し直し(run3): CONVERGED 6 反復、ギャップ 3.403743 eV（LDA 0.951647）、5 月との差 +0.02 eV |
| [mp-549308](https://next-gen.materialsproject.org/materials/mp-549308) | LiGeBO4 | 7 | 5.747476 | 回し直し | 2.586501 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1523 eV、履歴 6.7541, 5.4573, 5.5118, 5.5653, 5.6640） / 回し直し(run3): CONVERGED 6 反復、ギャップ 5.747476 eV（LDA 2.586501）、5 月との差 +0.08 eV |
| [mp-8873](https://next-gen.materialsproject.org/materials/mp-8873) | LiGeBO4 | 7 | 7.632831 | 回し直し | 4.346138 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1408 eV、履歴 7.4058, 7.6471, 7.7296, 7.8428, 7.8703） / 回し直し(run3): CONVERGED 6 反復、ギャップ 7.632831 eV（LDA 4.346138）、5 月との差 -0.24 eV |
| [mp-28591](https://next-gen.materialsproject.org/materials/mp-28591) | Na4HgP2 | 7 | 1.959142 | 回し直し | 0.798162 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2538 eV、履歴 1.7988, 2.9500, 2.3592, 2.2073, 2.1054） / 回し直し(run3): CONVERGED 5 反復、ギャップ 1.959142 eV（LDA 0.798162）、5 月との差 -0.15 eV |
| [mp-22937](https://next-gen.materialsproject.org/materials/mp-22937) | Na4I2O | 7 | 4.761284 | 回し直し | 2.336738 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2496 eV、履歴 4.1418, 5.5157, 4.8725, 4.7452, 4.6229） / 回し直し(run3): CONVERGED 6 反復、ギャップ 4.761284 eV（LDA 2.336738）、5 月との差 +0.14 eV |
| [mp-8962](https://next-gen.materialsproject.org/materials/mp-8962) | Rb2HgF4 | 7 | 5.085645 | 回し直し | 1.723822 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2100 eV、履歴 4.6337, 3.5013, 4.1973, 4.3807, 4.4073） / 回し直し(run3): CONVERGED 6 反復、ギャップ 5.085645 eV（LDA 1.723822）、5 月との差 +0.68 eV |
| [mp-7863](https://next-gen.materialsproject.org/materials/mp-7863) | Rb2Sn2O3 | 7 | 2.405706 | 回し直し | 1.030971 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1579 eV、履歴 2.3649, 2.5022, 2.6140, 2.4681, 2.4561） / 回し直し(run3): CONVERGED 5 反復、ギャップ 2.405706 eV（LDA 1.030971）、5 月との差 -0.05 eV |
| [mp-30004](https://next-gen.materialsproject.org/materials/mp-30004) | Rb4Br2O | 7 | 2.823486 | 回し直し | 0.502586 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 5.2682 eV、履歴 2.6275, 0.4777, 1.0943, 2.3450, 6.3625） / 回し直し(run3): CONVERGED 6 反復、ギャップ 2.823486 eV（LDA 0.502586）、5 月との差 -3.54 eV |
| [mp-3376](https://next-gen.materialsproject.org/materials/mp-3376) | Sr2SnO4 | 7 | 5.039359 | 回し直し | 2.712101 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1500 eV、履歴 5.4255, 5.1584, 5.0668, 5.1947, 5.2168） / 回し直し(run3): CONVERGED 5 反復、ギャップ 5.039359 eV（LDA 2.712101）、5 月との差 -0.18 eV |
| [mp-5532](https://next-gen.materialsproject.org/materials/mp-5532) | Sr2TiO4 | 7 | 4.063720 | 回し直し | 1.885675 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.6391 eV、履歴 6.8600, 2.0083, 3.6171, 3.7774, 4.2562） / 回し直し(run3): CONVERGED 7 反復、ギャップ 4.063720 eV（LDA 1.885675）、5 月との差 -0.19 eV |
| [mp-570097](https://next-gen.materialsproject.org/materials/mp-570097) | SrLi4P2 | 7 | 1.883971 | 回し直し | 0.868214 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.4129 eV、履歴 1.2886, 3.9361, 2.6857, 2.4557, 2.2728） / 回し直し(run3): CONVERGED 5 反復、ギャップ 1.883971 eV（LDA 0.868214）、5 月との差 -0.39 eV |
| [mp-571518](https://next-gen.materialsproject.org/materials/mp-571518) | WCl6 | 7 | 3.495973 | 回し直し | 1.586120 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1044 eV、履歴 4.7599, 3.4644, 4.0578, 3.9534, 3.9553） / 回し直し(run3): CONVERGED 6 反復、ギャップ 3.495973 eV（LDA 1.586120）、5 月との差 -0.46 eV |
| [mp-6070](https://next-gen.materialsproject.org/materials/mp-6070) | Ag2C2N2O2 | 8 | 4.640277 | 回し直し | 2.236520 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1995 eV、履歴 5.1089, 4.8196, 4.5193, 4.5728, 4.7188） / 回し直し(run3): CONVERGED 6 反復、ギャップ 4.640277 eV（LDA 2.236520）、5 月との差 -0.08 eV |
| [mp-1079156](https://next-gen.materialsproject.org/materials/mp-1079156) | Ag2Cl2O4 | 8 | 4.212064 | 回し直し | 0.656506 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.5221 eV、履歴 5.5393, 4.7417, 3.6603, 4.0895, 4.1823） / 回し直し(run3): CONVERGED 6 反復、ギャップ 4.212064 eV（LDA 0.656506）、5 月との差 +0.03 eV |
| [mp-3098](https://next-gen.materialsproject.org/materials/mp-3098) | Al2Cu2O4 | 8 | 3.862598 | 回し直し | 2.041132 | 5 月は 10 反復で収束せず（最後の 3 反復の幅  eV、履歴 3.8069） / 回し直し(run3): CONVERGED 6 反復、ギャップ 3.862598 eV（LDA 2.041132）、5 月との差 +0.06 eV |
| [mp-23807](https://next-gen.materialsproject.org/materials/mp-23807) | Al2H2O4 | 8 | 8.630608 | 回し直し | 5.596687 | 5 月は 10 反復で収束せず（最後の 3 反復の幅  eV、履歴 8.5507, 8.4200） / 回し直し(run3): CONVERGED 5 反復、ギャップ 8.630608 eV（LDA 5.596687）、5 月との差 +0.21 eV |
| [mp-23218](https://next-gen.materialsproject.org/materials/mp-23218) | As2I6 | 8 | 3.412576 | 回し直し | 2.100004 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.3962 eV、履歴 3.4245, 1.1146, 2.2955, 2.6231, 2.6918） / 回し直し(run3): CONVERGED 5 反復、ギャップ 3.412576 eV（LDA 2.100004）、5 月との差 +0.72 eV |
| [mp-23184](https://next-gen.materialsproject.org/materials/mp-23184) | B2Cl6 | 8 | 8.336585 | 回し直し | 4.394165 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2489 eV、履歴 8.1727, 9.7465, 9.0686, 8.8689, 8.8197） / 回し直し(run3): CONVERGED 6 反復、ギャップ 8.336585 eV（LDA 4.394165）、5 月との差 -0.48 eV |
| [mp-3104](https://next-gen.materialsproject.org/materials/mp-3104) | Ba2Zr2N4 | 8 | 2.465224 | 回し直し | 0.940354 | 5 月は 10 反復で収束せず（最後の 3 反復の幅  eV、履歴 2.4338, 2.4486） / 回し直し(run3): CONVERGED 5 反復、ギャップ 2.465224 eV（LDA 0.940354）、5 月との差 +0.02 eV |
| [mp-1078191](https://next-gen.materialsproject.org/materials/mp-1078191) | Ba4Te2O2 | 8 | 3.230588 | 回し直し | 1.741518 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1506 eV、履歴 3.1869, 3.3247, 3.1487, 3.0967, 2.9981） / 回し直し(run3): CONVERGED 5 反復、ギャップ 3.230588 eV（LDA 1.741518）、5 月との差 +0.23 eV |
| [mp-8290](https://next-gen.materialsproject.org/materials/mp-8290) | BaSnF6 | 8 | 9.693062 | 回し直し | 5.180288 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1250 eV、履歴 11.7033, 9.5621, 9.6313, 9.5063, 9.6195） / 回し直し(run3): CONVERGED 6 反復、ギャップ 9.693062 eV（LDA 5.180288）、5 月との差 +0.07 eV |
| [mp-8291](https://next-gen.materialsproject.org/materials/mp-8291) | BaTiF6 | 8 | 8.728307 | 回し直し | 4.357562 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1119 eV、履歴 8.3112, 8.5828, 8.5874, 8.6993, 8.5992） / 回し直し(run3): CONVERGED 7 反復、ギャップ 8.728307 eV（LDA 4.357562）、5 月との差 +0.13 eV |
| [mp-20463](https://next-gen.materialsproject.org/materials/mp-20463) | CaPbF6 | 8 | 8.320289 | 回し直し | 3.551001 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 1.0012 eV、履歴 9.4870, 13.7927, 11.2680, 10.5853, 10.2668） / 回し直し(run3): CONVERGED 6 反復、ギャップ 8.320289 eV（LDA 3.551001）、5 月との差 -1.95 eV |
| [mp-7922](https://next-gen.materialsproject.org/materials/mp-7922) | CaPdF6 | 8 | 6.314140 | 回し直し | 2.023140 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1691 eV、履歴 5.6216, 5.9572, 6.1395, 6.2826, 6.1135） / 回し直し(run3): CONVERGED 6 反復、ギャップ 6.314140 eV（LDA 2.023140）、5 月との差 +0.20 eV |
| [mp-4782](https://next-gen.materialsproject.org/materials/mp-4782) | CaSnF6 | 8 | 10.613974 | 回し直し | 5.496651 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1657 eV、履歴 20.3634, 8.6418, 9.8452, 10.0109, 9.9299） / 回し直し(run3): CONVERGED 6 反復、ギャップ 10.613974 eV（LDA 5.496651）、5 月との差 +0.68 eV |
| [mp-8224](https://next-gen.materialsproject.org/materials/mp-8224) | CaSnF6 | 8 | 10.662416 | 回し直し | 5.460173 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1068 eV、履歴 10.5961, 10.8020, 10.5809, 10.6080, 10.6877） / 回し直し(run3): CONVERGED 6 反復、ギャップ 10.662416 eV（LDA 5.460173）、5 月との差 -0.03 eV |
| [mp-13984](https://next-gen.materialsproject.org/materials/mp-13984) | CdPdF6 | 8 | 5.471617 | 回し直し | 1.758786 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.8548 eV、履歴 10.2481, 5.4021, 5.4225, 5.2364, 4.5677） / 回し直し(run3): CONVERGED 6 反復、ギャップ 5.471617 eV（LDA 1.758786）、5 月との差 +0.90 eV |
| [mp-13907](https://next-gen.materialsproject.org/materials/mp-13907) | CdSnF6 | 8 | 8.215922 | 回し直し | 3.772128 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 1.3383 eV、履歴 4.7790, 7.4448, 6.7781, 6.2019, 5.4398） / 回し直し(run3): CONVERGED 6 反復、ギャップ 8.215922 eV（LDA 3.772128）、5 月との差 +2.78 eV |
| [mp-29504](https://next-gen.materialsproject.org/materials/mp-29504) | Cl4F4 | 8 | 5.947637 | 回し直し | 1.846601 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2194 eV、履歴 5.6740, 6.1135, 5.5236, 5.7430, 5.7043） / 回し直し(run3): CONVERGED 6 反復、ギャップ 5.947637 eV（LDA 1.846601）、5 月との差 +0.24 eV |
| [mp-510421](https://next-gen.materialsproject.org/materials/mp-510421) | Cr2O6 | 8 | 4.237647 | 回し直し | 1.907167 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1388 eV、履歴 4.5760, 4.1494, 4.0739, 4.0106） / 回し直し(run3): CONVERGED 6 反復、ギャップ 4.237647 eV（LDA 1.907167）、5 月との差 +0.23 eV |
| [mp-542772](https://next-gen.materialsproject.org/materials/mp-542772) | Cs2Ag2Cl4 | 8 | 4.734911 | 回し直し | 2.093722 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.7394 eV、履歴 3.9103, 4.6497, 4.5523） / 回し直し(run3): CONVERGED 6 反復、ギャップ 4.734911 eV（LDA 2.093722）、5 月との差 +0.18 eV |
| [mp-24668](https://next-gen.materialsproject.org/materials/mp-24668) | Cs2H2F4 | 8 | 10.660157 | 回し直し | 6.282067 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.5788 eV、履歴 9.8755, 10.2717, 10.4543） / 回し直し(run3): CONVERGED 6 反復、ギャップ 10.660157 eV（LDA 6.282067）、5 月との差 +0.21 eV |
| [mp-510273](https://next-gen.materialsproject.org/materials/mp-510273) | Cs2Sb2O4 | 8 | 4.311969 | 回し直し | 2.346599 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1260 eV、履歴 4.4465, 4.4980, 4.3872, 4.3720） / 回し直し(run3): CONVERGED 5 反復、ギャップ 4.311969 eV（LDA 2.346599）、5 月との差 -0.06 eV |
| [mp-27422](https://next-gen.materialsproject.org/materials/mp-27422) | CsBiF6 | 8 | 7.077400 | 回し直し | 2.912728 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.3361 eV、履歴 0.1745, 1.9772, 7.1362, 7.1906, 7.4723） / 回し直し(run3): CONVERGED 6 反復、ギャップ 7.077400 eV（LDA 2.912728）、5 月との差 -0.39 eV |
| [mp-541226](https://next-gen.materialsproject.org/materials/mp-541226) | CsBrF6 | 8 | 7.192772 | 回し直し | 4.040081 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 6.1319 eV、履歴 5.9947, 0.0740, 2.0667, 6.0012, 8.1986） / 回し直し(run3): CONVERGED 6 反復、ギャップ 7.192772 eV（LDA 4.040081）、5 月との差 -1.01 eV |
| [mp-11980](https://next-gen.materialsproject.org/materials/mp-11980) | CsNbF6 | 8 | 9.226972 | 回し直し | 5.071790 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1429 eV、履歴 9.4029, 8.9795, 8.9974, 9.1387, 9.1403） / 回し直し(run3): CONVERGED 6 反復、ギャップ 9.226972 eV（LDA 5.071790）、5 月との差 +0.09 eV |
| [mp-30952](https://next-gen.materialsproject.org/materials/mp-30952) | Ga2Cl6 | 8 | 7.352972 | 回し直し | 3.880771 | 5 月は 10 反復で収束せず（最後の 3 反復の幅  eV、履歴 7.7362, 7.6817） / 回し直し(run3): CONVERGED 6 反復、ギャップ 7.352972 eV（LDA 3.880771）、5 月との差 -0.33 eV |
| [mp-588](https://next-gen.materialsproject.org/materials/mp-588) | Ga2F6 | 8 | 9.774229 | 回し直し | 4.981541 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2728 eV、履歴 9.0317, 6.1768, 7.6054, 7.8782, 7.8433） / 回し直し(run3): CONVERGED 6 反復、ギャップ 9.774229 eV（LDA 4.981541）、5 月との差 +1.93 eV |
| [mp-6949](https://next-gen.materialsproject.org/materials/mp-6949) | In2F6 | 8 | 8.137896 | 回し直し | 4.027866 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.4717 eV、履歴 7.7522, 6.7446, 6.5919, 6.9400, 7.0636） / 回し直し(run3): CONVERGED 6 反復、ギャップ 8.137896 eV（LDA 4.027866）、5 月との差 +1.07 eV |
| [mp-2437](https://next-gen.materialsproject.org/materials/mp-2437) | Ir2F6 | 8 | 4.371508 | 回し直し | 0.998863 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1392 eV、履歴 4.0828, 5.4634, 4.6499, 4.5462, 4.5107） / 回し直し(run3): CONVERGED 7 反復、ギャップ 4.371508 eV（LDA 0.998863）、5 月との差 -0.14 eV |
| [mp-3982](https://next-gen.materialsproject.org/materials/mp-3982) | K2Cu2O4 | 8 | 3.483811 | 回し直し | 0.730233 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.8992 eV、履歴 4.5659, 5.3839, 4.7090, 4.4027, 3.8098） / 回し直し(run3): CONVERGED 6 反復、ギャップ 3.483811 eV（LDA 0.730233）、5 月との差 -0.33 eV |
| [mp-24428](https://next-gen.materialsproject.org/materials/mp-24428) | K2H4N2 | 8 | 5.024534 | 回し直し | 1.934914 | 5 月は 10 反復で収束せず（最後の 3 反復の幅  eV、履歴 5.0490, 4.9640） / 回し直し(run3): CONVERGED 6 反復、ギャップ 5.024534 eV（LDA 1.934914）、5 月との差 +0.06 eV |
| [mp-19052](https://next-gen.materialsproject.org/materials/mp-19052) | K3VO4 | 8 | 6.175355 | 回し直し | 3.512554 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1571 eV、履歴 5.7678, 5.9652, 6.0167, 6.0850, 6.1738） / 回し直し(run3): CONVERGED 6 反復、ギャップ 6.175355 eV（LDA 3.512554）、5 月との差 +0.00 eV |
| [mp-7569](https://next-gen.materialsproject.org/materials/mp-7569) | KAsF6 | 8 | 10.543702 | 回し直し | 5.127676 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.8032 eV、履歴 61.0170, 12.9776, 10.2523, 10.7809, 11.0555） / 回し直し(run3): CONVERGED 6 反復、ギャップ 10.543702 eV（LDA 5.127676）、5 月との差 -0.51 eV |
| [mp-14232](https://next-gen.materialsproject.org/materials/mp-14232) | Li2B2O4 | 8 | 10.382166 | 回し直し | 7.063167 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1561 eV、履歴 10.3964, 11.2790, 10.8699, 10.8492, 10.7137） / 回し直し(run3): CONVERGED 5 反復、ギャップ 10.382166 eV（LDA 7.063167）、5 月との差 -0.33 eV |
| [mp-5840](https://next-gen.materialsproject.org/materials/mp-5840) | Li2Sc2O4 | 8 | 6.817251 | 回し直し | 3.789047 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 3.1342 eV、履歴 6.7239, 24.9336, 16.1958, 14.1226, 13.0616） / 回し直し(run3): CONVERGED 6 反復、ギャップ 6.817251 eV（LDA 3.789047）、5 月との差 -6.24 eV |
| [mp-1079483](https://next-gen.materialsproject.org/materials/mp-1079483) | LiAuF6 | 8 | 6.050029 | 回し直し | 1.904653 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2295 eV、履歴 5.9146, 6.1389, 6.1993, 5.9972, 5.9698） / 回し直し(run3): CONVERGED 6 反復、ギャップ 6.050029 eV（LDA 1.904653）、5 月との差 +0.08 eV |
| [mp-9143](https://next-gen.materialsproject.org/materials/mp-9143) | LiPF6 | 8 | 13.516106 | 回し直し | 8.046911 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 2.6132 eV、履歴 13.2709, 25.4050, 19.6095, 17.4852, 16.9963） / 回し直し(run3): CONVERGED 6 反復、ギャップ 13.516106 eV（LDA 8.046911）、5 月との差 -3.48 eV |
| [mp-3980](https://next-gen.materialsproject.org/materials/mp-3980) | LiSbF6 | 8 | 10.161204 | 回し直し | 5.111021 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.4743 eV、履歴 10.9972, 9.5229, 9.7837, 9.7319, 10.2062） / 回し直し(run3): CONVERGED 6 反復、ギャップ 10.161204 eV（LDA 5.111021）、5 月との差 -0.04 eV |
| [mp-1080679](https://next-gen.materialsproject.org/materials/mp-1080679) | LiTaF6 | 8 | 10.971758 | 回し直し | 6.245945 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.6910 eV、履歴 10.6205, 6.1081, 8.3128, 8.6016, 9.0038） / 回し直し(run3): CONVERGED 6 反復、ギャップ 10.971758 eV（LDA 6.245945）、5 月との差 +1.97 eV |
| [mp-7589](https://next-gen.materialsproject.org/materials/mp-7589) | Mg4N2F2 | 8 | 4.700299 | 回し直し | 2.340880 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.4142 eV、履歴 4.4640, 5.2561, 5.4018, 5.0570, 4.9876） / 回し直し(run3): CONVERGED 6 反復、ギャップ 4.700299 eV（LDA 2.340880）、5 月との差 -0.29 eV |
| [mp-19734](https://next-gen.materialsproject.org/materials/mp-19734) | MgPbF6 | 8 | 7.310837 | 回し直し | 3.147291 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 9.3996 eV、履歴 30.5251, 8.6004, 10.5076, 1.1080, 1.7637） / 回し直し(run3): CONVERGED 6 反復、ギャップ 7.310837 eV（LDA 3.147291）、5 月との差 +5.55 eV |
| [mp-7921](https://next-gen.materialsproject.org/materials/mp-7921) | MgPdF6 | 8 | 5.911959 | 回し直し | 1.820979 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2301 eV、履歴 6.4376, 4.3717, 5.1342, 5.3558, 5.3642） / 回し直し(run3): CONVERGED 6 反復、ギャップ 5.911959 eV（LDA 1.820979）、5 月との差 +0.55 eV |
| [mp-1095181](https://next-gen.materialsproject.org/materials/mp-1095181) | N6Cl2 | 8 | 6.507109 | 回し直し | 2.228610 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1841 eV、履歴 7.0367, 8.7138, 7.6739, 7.5383, 7.7224） / 回し直し(run3): CONVERGED 5 反復、ギャップ 6.507109 eV（LDA 2.228610）、5 月との差 -1.21 eV |
| [mp-25](https://next-gen.materialsproject.org/materials/mp-25) | N8 | 8 | 13.392300 | 回し直し | 7.118669 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1188 eV、履歴 13.7797, 13.1523, 13.2239, 13.3416, 13.3426） / 回し直し(run3): CONVERGED 5 反復、ギャップ 13.392300 eV（LDA 7.118669）、5 月との差 +0.05 eV |
| [mp-20343](https://next-gen.materialsproject.org/materials/mp-20343) | Na2Au2O4 | 8 | 3.282495 | 回し直し | 0.834700 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 1.7112 eV、履歴 3.2272, 12.6936, 8.3476, 6.8515, 6.6364） / 回し直し(run3): CONVERGED 5 反復、ギャップ 3.282495 eV（LDA 0.834700）、5 月との差 -3.35 eV |
| [mp-867515](https://next-gen.materialsproject.org/materials/mp-867515) | Na2Co2O4 | 8 | 3.678100 | 回し直し | 0.723683 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1091 eV、履歴 2.8570, 3.4392, 3.5383, 3.6123, 3.6475） / 回し直し(run3x(続き)): CONVERGED 19 反復、ギャップ 3.678100 eV（LDA 0.723683）、5 月との差 +0.03 eV |
| [mp-1087215](https://next-gen.materialsproject.org/materials/mp-1087215) | NaAsF6 | 8 | 11.187406 | 回し直し | 5.316576 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 1.8287 eV、履歴 10.9224, 20.3132, 15.7750, 14.5913, 13.9463） / 回し直し(run3): CONVERGED 6 反復、ギャップ 11.187406 eV（LDA 5.316576）、5 月との差 -2.76 eV |
| [mp-27420](https://next-gen.materialsproject.org/materials/mp-27420) | NaBiF6 | 8 | 7.241208 | 回し直し | 2.856984 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.4375 eV、履歴 8.1073, 6.7955, 7.3506, 7.1206, 6.9131） / 回し直し(run3): CONVERGED 6 反復、ギャップ 7.241208 eV（LDA 2.856984）、5 月との差 +0.33 eV |
| [mp-5955](https://next-gen.materialsproject.org/materials/mp-5955) | NaSbF6 | 8 | 10.623816 | 回し直し | 5.054327 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 1.0887 eV、履歴 11.1649, 3.6204, 7.0808, 7.9371, 8.1694） / 回し直し(run3): CONVERGED 6 反復、ギャップ 10.623816 eV（LDA 5.054327）、5 月との差 +2.45 eV |
| [mp-560008](https://next-gen.materialsproject.org/materials/mp-560008) | P2N2F4 | 8 | 9.133497 | 回し直し | 5.630316 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.4334 eV、履歴 9.4656, 9.1338, 8.7762, 9.2062, 9.2096） / 回し直し(run3): CONVERGED 6 反復、ギャップ 9.133497 eV（LDA 5.630316）、5 月との差 -0.08 eV |
| [mp-7472](https://next-gen.materialsproject.org/materials/mp-7472) | Rb2Au2O4 | 8 | 4.076355 | 回し直し | 1.390762 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 1.8813 eV、履歴 3.8953, 13.4869, 8.7181, 7.3384, 6.8368） / 回し直し(run3): CONVERGED 5 反復、ギャップ 4.076355 eV（LDA 1.390762）、5 月との差 -2.76 eV |
| [mp-9821](https://next-gen.materialsproject.org/materials/mp-9821) | RbSbF6 | 8 | 10.107447 | 回し直し | 5.106297 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1522 eV、履歴 10.0983, 10.1302, 10.2112, 10.0641, 10.0591） / 回し直し(run3): CONVERGED 6 反復、ギャップ 10.107447 eV（LDA 5.106297）、5 月との差 +0.05 eV |
| [mp-1880](https://next-gen.materialsproject.org/materials/mp-1880) | Sb2F6 | 8 | 7.305903 | 回し直し | 4.479529 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2103 eV、履歴 7.1597, 7.3202, 7.4007, 7.2090, 7.1903） / 回し直し(run3): CONVERGED 6 反復、ギャップ 7.305903 eV（LDA 4.479529）、5 月との差 +0.12 eV |
| [mp-1078300](https://next-gen.materialsproject.org/materials/mp-1078300) | Sc2F6 | 8 | 11.939590 | 回し直し | 6.051292 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 2.3103 eV、履歴 11.6013, 8.9393, 11.2231, 13.5334, 11.9499） / 回し直し(run3): CONVERGED 6 反復、ギャップ 11.939590 eV（LDA 6.051292）、5 月との差 -0.01 eV |
| [mp-24066](https://next-gen.materialsproject.org/materials/mp-24066) | Sr2H2Cl2O2 | 8 | 7.864377 | 回し直し | 4.769609 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.1902 eV、履歴 5.6615, 4.8322, 5.3526, 5.4451, 5.5428） / 回し直し(run3): CONVERGED 6 反復、ギャップ 7.864377 eV（LDA 4.769609）、5 月との差 +2.32 eV |
| [mp-623024](https://next-gen.materialsproject.org/materials/mp-623024) | TiCdF6 | 8 | 8.826431 | 回し直し | 4.394530 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.3688 eV、履歴 6.8584, 8.7834, 8.0502, 7.7964, 7.6814） / 回し直し(run3): CONVERGED 7 反復、ギャップ 8.826431 eV（LDA 4.394530）、5 月との差 +1.15 eV |
| [mp-27455](https://next-gen.materialsproject.org/materials/mp-27455) | Y2Cl6 | 8 | 7.878928 | 回し直し | 4.544462 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.3007 eV、履歴 8.7785, 7.7956, 8.0739, 7.9432, 7.7732） / 回し直し(run3): CONVERGED 6 反復、ギャップ 7.878928 eV（LDA 4.544462）、5 月との差 +0.11 eV |
| [mp-2918](https://next-gen.materialsproject.org/materials/mp-2918) | Y2Cu2O4 | 8 | 4.463069 | 回し直し | 2.727645 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.3048 eV、履歴 4.5307, 3.1156, 3.8782, 4.0487, 3.7439） / 回し直し(run3): CONVERGED 6 反復、ギャップ 4.463069 eV（LDA 2.727645）、5 月との差 +0.72 eV |
| [mp-3595](https://next-gen.materialsproject.org/materials/mp-3595) | Zn2Si2As4 | 8 | 1.713093 | 回し直し | 0.951475 | 5 月は 10 反復で収束せず（最後の 3 反復の幅  eV、履歴 1.7188, 1.7334） / 回し直し(run3): CONVERGED 4 反復、ギャップ 1.713093 eV（LDA 0.951475）、5 月との差 -0.02 eV |
| [mp-4175](https://next-gen.materialsproject.org/materials/mp-4175) | Zn2Sn2P4 | 8 | 1.577509 | 回し直し | 0.656902 | 5 月は 10 反復で収束せず（最後の 3 反復の幅  eV、履歴 1.5747, 1.5791） / 回し直し(run3): CONVERGED 5 反復、ギャップ 1.577509 eV（LDA 0.656902）、5 月との差 -0.00 eV |
| [mp-13610](https://next-gen.materialsproject.org/materials/mp-13610) | ZnPbF6 | 8 | 5.701529 | 回し直し | 1.933201 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.4479 eV、履歴 4.4228, 3.4626, 4.0115, 4.3384, 4.4593） / 回し直し(run3): CONVERGED 6 反復、ギャップ 5.701529 eV（LDA 1.933201）、5 月との差 +1.24 eV |
| [mp-13983](https://next-gen.materialsproject.org/materials/mp-13983) | ZnPdF6 | 8 | 5.102630 | 回し直し | 1.568340 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.2447 eV、履歴 3.9438, 4.7235, 4.7664, 4.5217, 4.5273） / 回し直し(run3): CONVERGED 6 反復、ギャップ 5.102630 eV（LDA 1.568340）、5 月との差 +0.58 eV |
| [mp-8256](https://next-gen.materialsproject.org/materials/mp-8256) | ZnPtF6 | 8 | 6.002285 | 回し直し | 2.007324 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 14.0367 eV、履歴 4.2993, 3.7165, 18.5010, 4.4642, 5.5725） / 回し直し(run3): CONVERGED 6 反復、ギャップ 6.002285 eV（LDA 2.007324）、5 月との差 +0.43 eV |
| [mp-13903](https://next-gen.materialsproject.org/materials/mp-13903) | ZnSnF6 | 8 | 8.365340 | 回し直し | 3.696367 | 5 月は 10 反復で収束せず（最後の 3 反復の幅 0.3282 eV、履歴 6.8517, 9.9907, 7.7057, 7.3775, 7.4443） / 回し直し(run3): CONVERGED 6 反復、ギャップ 8.365340 eV（LDA 3.696367）、5 月との差 +0.92 eV |

## DRIFT_GOOD

5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに 0.08 eV 以上動いた。回し直していない。62 物質。

**表 9**

| mpid | 組成 | 原子 | 採用のギャップ | どちらの値 | LDA | 注記 |
| --- | --- | --- | --- | --- | --- | --- |
| [mp-22913](https://next-gen.materialsproject.org/materials/mp-22913) | CuBr | 2 | 2.410 | 5 月 | 0.439 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 2.4606, 2.3185, 2.3865, 2.4095）。回し直していない |
| [mp-22905](https://next-gen.materialsproject.org/materials/mp-22905) | LiCl | 2 | 9.384 | 5 月 | 6.415 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 9.6912, 9.4714, 9.4389, 9.3841）。回し直していない |
| [mp-23870](https://next-gen.materialsproject.org/materials/mp-23870) | NaH | 2 | 6.732 | 5 月 | 3.471 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 6.6548, 6.8240, 6.7452, 6.7317）。回し直していない |
| [mp-1018141](https://next-gen.materialsproject.org/materials/mp-1018141) | BaC2 | 3 | 3.254 | 5 月 | 1.330 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 2.8517, 3.1627, 3.2186, 3.2544）。回し直していない |
| [mp-8177](https://next-gen.materialsproject.org/materials/mp-8177) | HgF2 | 3 | 3.676 | 5 月 | 0.794 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.10 eV 動いていた（履歴 3.8424, 3.3402, 3.5769, 3.6094, 3.6759）。回し直していない |
| [mp-570259](https://next-gen.materialsproject.org/materials/mp-570259) | MgCl2 | 3 | 8.381 | 5 月 | 5.256 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.10 eV 動いていた（履歴 8.4774, 8.4257, 8.3814）。回し直していない |
| [mp-2630](https://next-gen.materialsproject.org/materials/mp-2630) | SrC2 | 3 | 3.833 | 5 月 | 1.554 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 3.3838, 3.6153, 3.7708, 3.9214, 3.8500, 3.8334）。回し直していない |
| [mp-23162](https://next-gen.materialsproject.org/materials/mp-23162) | ZrCl2 | 3 | 1.882 | 5 月 | 0.971 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 1.7770, 1.9192, 1.7892, 1.8618, 1.8820）。回し直していない |
| [mp-570858](https://next-gen.materialsproject.org/materials/mp-570858) | Ag2Cl2 | 4 | 2.279 | 5 月 | 0.296 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.10 eV 動いていた（履歴 2.2404, 2.0995, 2.1834, 2.1998, 2.2794）。回し直していない |
| [mp-22894](https://next-gen.materialsproject.org/materials/mp-22894) | Ag2I2 | 4 | 3.087 | 5 月 | 1.160 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 2.9934, 3.0658, 3.0875）。回し直していない |
| [mp-1018033](https://next-gen.materialsproject.org/materials/mp-1018033) | AgRhO2 | 4 | 1.557 | 5 月 | 0.571 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 1.5928, 1.3423, 1.4710, 1.5570, 1.5574）。回し直していない |
| [mp-23227](https://next-gen.materialsproject.org/materials/mp-23227) | Cu2Br2 | 4 | 2.336 | 5 月 | 0.386 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 2.4196, 2.2432, 2.3327, 2.3358）。回し直していない |
| [mp-14116](https://next-gen.materialsproject.org/materials/mp-14116) | CuRhO2 | 4 | 1.936 | 5 月 | 0.806 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 2.0217, 1.9768, 1.9365）。回し直していない |
| [mp-23153](https://next-gen.materialsproject.org/materials/mp-23153) | I4 | 4 | 2.147 | 5 月 | 0.886 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 2.2324, 2.1523, 2.1467）。回し直していない |
| [mp-571636](https://next-gen.materialsproject.org/materials/mp-571636) | In2Cl2 | 4 | 2.467 | 5 月 | 1.382 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 2.3875, 2.5549, 2.5326, 2.4672）。回し直していない |
| [mp-571454](https://next-gen.materialsproject.org/materials/mp-571454) | LiAgC2 | 4 | 3.287 | 5 月 | 1.457 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 3.3755, 3.3072, 3.2874）。回し直していない |
| [mp-10422](https://next-gen.materialsproject.org/materials/mp-10422) | NaAuC2 | 4 | 2.969 | 5 月 | 1.246 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 2.6374, 2.7863, 2.9021, 3.0603, 2.9862, 2.9692）。回し直していない |
| [mp-20289](https://next-gen.materialsproject.org/materials/mp-20289) | NaInS2 | 4 | 3.027 | 5 月 | 1.621 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 3.4261, 2.5968, 2.9416, 2.9916, 3.0271）。回し直していない |
| [mp-999448](https://next-gen.materialsproject.org/materials/mp-999448) | NaYSe2 | 4 | 3.107 | 5 月 | 1.706 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.10 eV 動いていた（履歴 3.1934, 3.0115, 3.0501, 3.1075）。回し直していない |
| [mp-1063670](https://next-gen.materialsproject.org/materials/mp-1063670) | Pb2Se2 | 4 | 1.364 | 5 月 | 0.841 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 1.2742, 1.3277, 1.3641）。回し直していない |
| [mp-10694](https://next-gen.materialsproject.org/materials/mp-10694) | ScF3 | 4 | 11.653 | 5 月 | 6.032 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 11.5653, 11.6507, 11.6530）。回し直していない |
| [mp-665922](https://next-gen.materialsproject.org/materials/mp-665922) | SrHgO2 | 4 | 3.470 | 5 月 | 2.176 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.08 eV 動いていた（履歴 3.5140, 3.2745, 3.5540, 3.4719, 3.4703）。回し直していない |
| [mp-4255](https://next-gen.materialsproject.org/materials/mp-4255) | BaCu2S2 | 5 | 1.523 | 5 月 | 0.467 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 1.4751, 1.6084, 1.5264, 1.5231）。回し直していない |
| [mp-541837](https://next-gen.materialsproject.org/materials/mp-541837) | Bi2Se3 | 5 | 0.715 | 5 月 | 0.318 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.10 eV 動いていた（履歴 0.8131, 0.7999, 0.7149）。回し直していない |
| [mp-570231](https://next-gen.materialsproject.org/materials/mp-570231) | CsCdBr3 | 5 | 2.782 | 5 月 | 0.506 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 2.6923, 2.7642, 2.7817）。回し直していない |
| [mp-561947](https://next-gen.materialsproject.org/materials/mp-561947) | CsHgF3 | 5 | 3.431 | 5 月 | 0.507 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.10 eV 動いていた（履歴 3.5636, 3.3327, 3.4177, 3.4314）。回し直していない |
| [mp-23037](https://next-gen.materialsproject.org/materials/mp-23037) | CsPbCl3 | 5 | 3.769 | 5 月 | 1.937 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.08 eV 動いていた（履歴 3.7059, 3.8536, 3.8032, 3.7692）。回し直していない |
| [mp-552729](https://next-gen.materialsproject.org/materials/mp-552729) | KIO3 | 5 | 4.549 | 5 月 | 2.437 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.10 eV 動いていた（履歴 4.4663, 4.6479, 4.6189, 4.5490）。回し直していない |
| [mp-3614](https://next-gen.materialsproject.org/materials/mp-3614) | KTaO3 | 5 | 3.996 | 5 月 | 1.963 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 3.9014, 3.9726, 3.9964）。回し直していない |
| [mp-9610](https://next-gen.materialsproject.org/materials/mp-9610) | Li2CN2 | 5 | 5.858 | 5 月 | 3.376 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 5.5625, 5.7640, 5.7754, 5.8581）。回し直していない |
| [mp-2706](https://next-gen.materialsproject.org/materials/mp-2706) | SnF4 | 5 | 7.305 | 5 月 | 3.042 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 7.1396, 7.5076, 7.3936, 7.3069, 7.3051）。回し直していない |
| [mp-3323](https://next-gen.materialsproject.org/materials/mp-3323) | SrZrO3 | 5 | 5.498 | 5 月 | 3.048 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.08 eV 動いていた（履歴 5.4170, 5.4906, 5.4983）。回し直していない |
| [mp-1071955](https://next-gen.materialsproject.org/materials/mp-1071955) | AlPS4 | 6 | 4.804 | 5 月 | 2.420 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.08 eV 動いていた（履歴 5.4465, 4.6563, 4.7228, 4.7973, 4.8042）。回し直していない |
| [mp-23861](https://next-gen.materialsproject.org/materials/mp-23861) | Ba2H2Cl2 | 6 | 5.909 | 5 月 | 3.158 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 5.8169, 5.8998, 5.9093）。回し直していない |
| [mp-730189](https://next-gen.materialsproject.org/materials/mp-730189) | C2Br2N2 | 6 | 8.291 | 5 月 | 4.771 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 8.4710, 8.3846, 8.3677, 8.2908）。回し直していない |
| [mp-23214](https://next-gen.materialsproject.org/materials/mp-23214) | Ca2Cl4 | 6 | 8.383 | 5 月 | 5.117 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 8.4768, 8.4049, 8.3830）。回し直していない |
| [mp-1071269](https://next-gen.materialsproject.org/materials/mp-1071269) | Hg3Te3 | 6 | 1.164 | 5 月 | 0.593 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.08 eV 動いていた（履歴 1.0835, 1.1595, 1.1641）。回し直していない |
| [mp-1018809](https://next-gen.materialsproject.org/materials/mp-1018809) | Mo2S4 | 6 | 1.990 | 5 月 | 1.229 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 2.0763, 2.0088, 1.9896）。回し直していない |
| [mp-613652](https://next-gen.materialsproject.org/materials/mp-613652) | Pb2Cl2F2 | 6 | 5.840 | 5 月 | 3.562 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.08 eV 動いていた（履歴 5.7546, 5.8186, 5.8396）。回し直していない |
| [mp-569017](https://next-gen.materialsproject.org/materials/mp-569017) | Pd2I4 | 6 | 1.892 | 5 月 | 0.616 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.08 eV 動いていた（履歴 1.8198, 1.9745, 1.9016, 1.8924）。回し直していない |
| [mp-567484](https://next-gen.materialsproject.org/materials/mp-567484) | Pt2Cl4 | 6 | 3.226 | 5 月 | 0.657 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.10 eV 動いていた（履歴 3.3247, 3.2373, 3.2259）。回し直していない |
| [mp-9845](https://next-gen.materialsproject.org/materials/mp-9845) | Rb2Ca2As2 | 6 | 2.437 | 5 月 | 1.137 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.10 eV 動いていた（履歴 2.3798, 2.5360, 2.4425, 2.4372）。回し直していない |
| [mp-8453](https://next-gen.materialsproject.org/materials/mp-8453) | Rb2Na2O2 | 6 | 4.741 | 5 月 | 2.145 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 4.6459, 4.7090, 4.7409）。回し直していない |
| [mp-546279](https://next-gen.materialsproject.org/materials/mp-546279) | Sc2Br2O2 | 6 | 6.344 | 5 月 | 3.097 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.08 eV 動いていた（履歴 6.4353, 6.2626, 6.3418, 6.3435）。回し直していない |
| [mp-4359](https://next-gen.materialsproject.org/materials/mp-4359) | Sr2PdO3 | 6 | 2.610 | 5 月 | 0.062 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 2.1326, 3.1425, 2.7010, 2.6114, 2.6101）。回し直していない |
| [mp-552120](https://next-gen.materialsproject.org/materials/mp-552120) | Y2Cl2O2 | 6 | 6.501 | 5 月 | 4.256 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 4.7872, 6.5909, 6.5649, 6.5010）。回し直していない |
| [mp-20098](https://next-gen.materialsproject.org/materials/mp-20098) | Ba2PbO4 | 7 | 2.572 | 5 月 | 1.350 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.10 eV 動いていた（履歴 2.6568, 2.4744, 2.5707, 2.5722）。回し直していない |
| [mp-29342](https://next-gen.materialsproject.org/materials/mp-29342) | Ca3PCl3 | 7 | 4.026 | 5 月 | 1.625 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 4.1410, 3.9330, 4.0189, 4.0261）。回し直していない |
| [mp-27207](https://next-gen.materialsproject.org/materials/mp-27207) | K2MgCl4 | 7 | 7.965 | 5 月 | 4.600 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 7.6556, 8.0462, 7.8796, 7.9155, 7.9650）。回し直していない |
| [mp-27243](https://next-gen.materialsproject.org/materials/mp-27243) | K2PtBr4 | 7 | 3.968 | 5 月 | 1.226 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 3.9042, 4.0612, 3.9815, 3.9676）。回し直していない |
| [mp-8861](https://next-gen.materialsproject.org/materials/mp-8861) | Rb2MgF4 | 7 | 10.584 | 5 月 | 6.249 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 10.2663, 10.8874, 10.6697, 10.6195, 10.5837）。回し直していない |
| [mp-567655](https://next-gen.materialsproject.org/materials/mp-567655) | Sr2CN2Cl2 | 7 | 6.286 | 5 月 | 4.010 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 6.5630, 6.1031, 6.1981, 6.2642, 6.2857）。回し直していない |
| [mp-23297](https://next-gen.materialsproject.org/materials/mp-23297) | Br2F6 | 8 | 5.760 | 5 月 | 2.032 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.10 eV 動いていた（履歴 5.8202, 5.6641, 5.7471, 5.7597）。回し直していない |
| [mp-569416](https://next-gen.materialsproject.org/materials/mp-569416) | C8 | 8 | 4.968 | 5 月 | 1.797 | 構造: 硝酸を挿入した黒鉛から N・O を除いたもの。密度 1.67（黒鉛 2.26） / 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.08 eV 動いていた（履歴 5.0504, 4.9883, 4.9684）。回し直していない |
| [mp-23364](https://next-gen.materialsproject.org/materials/mp-23364) | Cs2Li2Cl4 | 8 | 7.838 | 5 月 | 4.882 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 7.8795, 7.7507, 7.8048, 7.8384）。回し直していない |
| [mp-5518](https://next-gen.materialsproject.org/materials/mp-5518) | Ga2Ag2Se4 | 8 | 1.554 | 5 月 | 0.162 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.09 eV 動いていた（履歴 1.4928, 1.6482, 1.5575, 1.5543）。回し直していない |
| [mp-23803](https://next-gen.materialsproject.org/materials/mp-23803) | Ga2H2O4 | 8 | 5.849 | 5 月 | 3.334 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.10 eV 動いていた（履歴 5.8770, 5.7502, 5.8429, 5.8492）。回し直していない |
| [mp-3924](https://next-gen.materialsproject.org/materials/mp-3924) | Li2Nb2O4 | 8 | 2.797 | 5 月 | 1.584 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 2.8267, 2.7048, 2.7822, 2.7972）。回し直していない |
| [mp-558312](https://next-gen.materialsproject.org/materials/mp-558312) | Li4O4 | 8 | 5.613 | 5 月 | 1.518 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに -0.08 eV 動いていた（履歴 5.5089, 5.6961, 5.6438, 5.6127）。回し直していない |
| [mp-23309](https://next-gen.materialsproject.org/materials/mp-23309) | Sc2Cl6 | 8 | 6.296 | 5 月 | 3.592 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 5.9855, 6.6573, 6.2043, 6.2940, 6.2962）。回し直していない |
| [mp-1025549](https://next-gen.materialsproject.org/materials/mp-1025549) | Tl3VSe4 | 8 | 2.330 | 5 月 | 1.523 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.10 eV 動いていた（履歴 2.2324, 2.2922, 2.3302）。回し直していない |
| [mp-9921](https://next-gen.materialsproject.org/materials/mp-9921) | Zr2S6 | 8 | 2.410 | 5 月 | 1.030 | 5 月は収束の基準を満たしたが、最後の 3 反復が同じ向きに +0.09 eV 動いていた（履歴 2.3225, 2.3966, 2.4104）。回し直していない |

## GOOD

5 月に収束、注記なし。回し直していない。1125 物質。

**表 10**

| mpid | 組成 | 原子 | 採用のギャップ | どちらの値 | LDA | 注記 |
| --- | --- | --- | --- | --- | --- | --- |
| [mp-23231](https://next-gen.materialsproject.org/materials/mp-23231) | AgBr | 2 | 2.439 | 5 月 | 0.529 |  |
| [mp-22922](https://next-gen.materialsproject.org/materials/mp-22922) | AgCl | 2 | 2.916 | 5 月 | 0.752 |  |
| [mp-22919](https://next-gen.materialsproject.org/materials/mp-22919) | AgI | 2 | 1.957 | 5 月 | 0.548 |  |
| [mp-22925](https://next-gen.materialsproject.org/materials/mp-22925) | AgI | 2 | 3.013 | 5 月 | 1.129 |  |
| [mp-2172](https://next-gen.materialsproject.org/materials/mp-2172) | AlAs | 2 | 2.216 | 5 月 | 1.358 |  |
| [mp-1330](https://next-gen.materialsproject.org/materials/mp-1330) | AlN | 2 | 6.387 | 5 月 | 4.308 |  |
| [mp-1550](https://next-gen.materialsproject.org/materials/mp-1550) | AlP | 2 | 2.490 | 5 月 | 1.462 |  |
| [mp-2624](https://next-gen.materialsproject.org/materials/mp-2624) | AlSb | 2 | 1.888 | 5 月 | 1.194 |  |
| [mp-10044](https://next-gen.materialsproject.org/materials/mp-10044) | BAs | 2 | 1.867 | 5 月 | 1.173 |  |
| [mp-1639](https://next-gen.materialsproject.org/materials/mp-1639) | BN | 2 | 6.250 | 5 月 | 4.389 |  |
| [mp-685145](https://next-gen.materialsproject.org/materials/mp-685145) | BN | 2 | 5.871 | 5 月 | 3.594 |  |
| [mp-1479](https://next-gen.materialsproject.org/materials/mp-1479) | BP | 2 | 2.063 | 5 月 | 1.225 |  |
| [mp-10680](https://next-gen.materialsproject.org/materials/mp-10680) | BaSe | 2 | 2.751 | 5 月 | 1.186 |  |
| [mp-1000](https://next-gen.materialsproject.org/materials/mp-1000) | BaTe | 2 | 2.933 | 5 月 | 1.547 |  |
| [mp-1794](https://next-gen.materialsproject.org/materials/mp-1794) | BeO | 2 | 11.369 | 5 月 | 8.192 |  |
| [mp-422](https://next-gen.materialsproject.org/materials/mp-422) | BeS | 2 | 4.580 | 5 月 | 2.934 |  |
| [mp-1541](https://next-gen.materialsproject.org/materials/mp-1541) | BeSe | 2 | 3.935 | 5 月 | 2.513 |  |
| [mp-252](https://next-gen.materialsproject.org/materials/mp-252) | BeTe | 2 | 2.881 | 5 月 | 1.925 |  |
| [mp-66](https://next-gen.materialsproject.org/materials/mp-66) | C2 | 2 | 5.808 | 5 月 | 4.172 |  |
| [mp-1672](https://next-gen.materialsproject.org/materials/mp-1672) | CaS | 2 | 4.385 | 5 月 | 2.176 |  |
| [mp-2469](https://next-gen.materialsproject.org/materials/mp-2469) | CdS | 2 | 2.336 | 5 月 | 0.822 |  |
| [mp-370](https://next-gen.materialsproject.org/materials/mp-370) | CdS | 2 | 1.329 | 5 月 |  |  |
| [mp-2691](https://next-gen.materialsproject.org/materials/mp-2691) | CdSe | 2 | 1.736 | 5 月 | 0.287 |  |
| [mp-406](https://next-gen.materialsproject.org/materials/mp-406) | CdTe | 2 | 1.664 | 5 月 | 0.404 |  |
| [mp-1057286](https://next-gen.materialsproject.org/materials/mp-1057286) | CsH | 2 | 4.596 | 5 月 | 1.959 |  |
| [mp-1056920](https://next-gen.materialsproject.org/materials/mp-1056920) | CsI | 2 | 5.748 | 5 月 | 3.405 |  |
| [mp-570354](https://next-gen.materialsproject.org/materials/mp-570354) | CuBr | 2 | 1.925 | 5 月 | 0.156 |  |
| [mp-22914](https://next-gen.materialsproject.org/materials/mp-22914) | CuCl | 2 | 2.727 | 5 月 | 0.575 |  |
| [mp-571386](https://next-gen.materialsproject.org/materials/mp-571386) | CuCl | 2 | 2.324 | 5 月 | 0.358 |  |
| [mp-22895](https://next-gen.materialsproject.org/materials/mp-22895) | CuI | 2 | 2.857 | 5 月 | 1.110 |  |
| [mp-2534](https://next-gen.materialsproject.org/materials/mp-2534) | GaAs | 2 | 1.044 | 5 月 | 0.000 |  |
| [mp-2853](https://next-gen.materialsproject.org/materials/mp-2853) | GaN | 2 | 1.507 | 5 月 | 0.219 |  |
| [mp-2490](https://next-gen.materialsproject.org/materials/mp-2490) | GaP | 2 | 2.291 | 5 月 | 1.496 |  |
| [mp-10759](https://next-gen.materialsproject.org/materials/mp-10759) | GeSe | 2 | 0.195 | 5 月 | 0.131 |  |
| [mp-2612](https://next-gen.materialsproject.org/materials/mp-2612) | GeTe | 2 | 0.441 | 5 月 | 0.411 |  |
| [mp-938](https://next-gen.materialsproject.org/materials/mp-938) | GeTe | 2 | 0.662 | 5 月 | 0.545 |  |
| [mp-632291](https://next-gen.materialsproject.org/materials/mp-632291) | H2 | 2 | 14.632 | 5 月 | 8.312 |  |
| [mp-20351](https://next-gen.materialsproject.org/materials/mp-20351) | InP | 2 | 1.230 | 5 月 | 0.366 |  |
| [mp-23193](https://next-gen.materialsproject.org/materials/mp-23193) | KCl | 2 | 8.191 | 5 月 | 4.960 |  |
| [mp-23289](https://next-gen.materialsproject.org/materials/mp-23289) | KCl | 2 | 7.763 | 5 月 | 4.660 |  |
| [mp-463](https://next-gen.materialsproject.org/materials/mp-463) | KF | 2 | 10.516 | 5 月 | 6.110 |  |
| [mp-8454](https://next-gen.materialsproject.org/materials/mp-8454) | KF | 2 | 10.544 | 5 月 | 6.313 |  |
| [mp-24084](https://next-gen.materialsproject.org/materials/mp-24084) | KH | 2 | 6.325 | 5 月 | 3.068 |  |
| [mp-23259](https://next-gen.materialsproject.org/materials/mp-23259) | LiBr | 2 | 7.595 | 5 月 | 5.070 |  |
| [mp-1138](https://next-gen.materialsproject.org/materials/mp-1138) | LiF | 2 | 13.483 | 5 月 | 8.789 |  |
| [mp-23703](https://next-gen.materialsproject.org/materials/mp-23703) | LiH | 2 | 5.546 | 5 月 | 2.626 |  |
| [mp-22899](https://next-gen.materialsproject.org/materials/mp-22899) | LiI | 2 | 6.218 | 5 月 | 4.296 |  |
| [mp-1009129](https://next-gen.materialsproject.org/materials/mp-1009129) | MgO | 2 | 6.122 | 5 月 | 3.039 |  |
| [mp-1265](https://next-gen.materialsproject.org/materials/mp-1265) | MgO | 2 | 8.326 | 5 月 | 4.891 |  |
| [mp-1315](https://next-gen.materialsproject.org/materials/mp-1315) | MgS | 2 | 4.516 | 5 月 | 2.568 |  |
| [mp-10760](https://next-gen.materialsproject.org/materials/mp-10760) | MgSe | 2 | 3.255 | 5 月 | 1.653 |  |
| [mp-22916](https://next-gen.materialsproject.org/materials/mp-22916) | NaBr | 2 | 7.077 | 5 月 | 4.110 |  |
| [mp-22851](https://next-gen.materialsproject.org/materials/mp-22851) | NaCl | 2 | 6.790 | 5 月 | 3.885 |  |
| [mp-22862](https://next-gen.materialsproject.org/materials/mp-22862) | NaCl | 2 | 8.346 | 5 月 | 5.023 |  |
| [mp-23268](https://next-gen.materialsproject.org/materials/mp-23268) | NaI | 2 | 5.908 | 5 月 | 3.509 |  |
| [mp-21276](https://next-gen.materialsproject.org/materials/mp-21276) | PbS | 2 | 0.752 | 5 月 | 0.334 |  |
| [mp-2201](https://next-gen.materialsproject.org/materials/mp-2201) | PbSe | 2 | 0.722 | 5 月 | 0.303 |  |
| [mp-2064](https://next-gen.materialsproject.org/materials/mp-2064) | RbF | 2 | 10.276 | 5 月 | 6.498 |  |
| [mp-24721](https://next-gen.materialsproject.org/materials/mp-24721) | RbH | 2 | 5.431 | 5 月 | 2.528 |  |
| [mp-149](https://next-gen.materialsproject.org/materials/mp-149) | Si2 | 2 | 1.199 | 5 月 | 0.532 |  |
| [mp-8062](https://next-gen.materialsproject.org/materials/mp-8062) | SiC | 2 | 2.342 | 5 月 | 1.324 |  |
| [mp-10013](https://next-gen.materialsproject.org/materials/mp-10013) | SnS | 2 | 0.549 | 5 月 | 0.059 |  |
| [mp-1876](https://next-gen.materialsproject.org/materials/mp-1876) | SnS | 2 | 0.789 | 5 月 | 0.626 |  |
| [mp-2693](https://next-gen.materialsproject.org/materials/mp-2693) | SnSe | 2 | 0.789 | 5 月 | 0.580 |  |
| [mp-16364](https://next-gen.materialsproject.org/materials/mp-16364) | SnTe | 2 | 0.304 | 5 月 | 0.012 |  |
| [mp-1883](https://next-gen.materialsproject.org/materials/mp-1883) | SnTe | 2 | 0.366 | 5 月 | 0.004 |  |
| [mp-2472](https://next-gen.materialsproject.org/materials/mp-2472) | SrO | 2 | 5.428 | 5 月 | 3.183 |  |
| [mp-10627](https://next-gen.materialsproject.org/materials/mp-10627) | SrS | 2 | 3.326 | 5 月 | 1.622 |  |
| [mp-10653](https://next-gen.materialsproject.org/materials/mp-10653) | SrTe | 2 | 1.363 | 5 月 | 0.347 |  |
| [mp-19717](https://next-gen.materialsproject.org/materials/mp-19717) | TePb | 2 | 1.101 | 5 月 | 0.740 |  |
| [mp-22875](https://next-gen.materialsproject.org/materials/mp-22875) | TlBr | 2 | 3.066 | 5 月 | 1.853 |  |
| [mp-568560](https://next-gen.materialsproject.org/materials/mp-568560) | TlBr | 2 | 3.726 | 5 月 | 2.060 |  |
| [mp-23167](https://next-gen.materialsproject.org/materials/mp-23167) | TlCl | 2 | 3.481 | 5 月 | 2.007 |  |
| [mp-571102](https://next-gen.materialsproject.org/materials/mp-571102) | TlI | 2 | 3.252 | 5 月 | 1.973 |  |
| [mp-2114](https://next-gen.materialsproject.org/materials/mp-2114) | YN | 2 | 0.997 | 5 月 | 0.115 |  |
| [mp-2229](https://next-gen.materialsproject.org/materials/mp-2229) | ZnO | 2 | 2.744 | 5 月 | 0.734 |  |
| [mp-10695](https://next-gen.materialsproject.org/materials/mp-10695) | ZnS | 2 | 3.677 | 5 月 | 1.918 |  |
| [mp-1190](https://next-gen.materialsproject.org/materials/mp-1190) | ZnSe | 2 | 2.729 | 5 月 | 1.060 |  |
| [mp-2176](https://next-gen.materialsproject.org/materials/mp-2176) | ZnTe | 2 | 2.434 | 5 月 | 1.011 |  |
| [mp-9519](https://next-gen.materialsproject.org/materials/mp-9519) | AgCN | 3 | 6.775 | 5 月 | 3.592 |  |
| [mp-7188](https://next-gen.materialsproject.org/materials/mp-7188) | Al2Os | 3 | 1.006 | 5 月 | 0.407 |  |
| [mp-29196](https://next-gen.materialsproject.org/materials/mp-29196) | AuCN | 3 | 3.877 | 5 月 | 1.772 |  |
| [mp-1735](https://next-gen.materialsproject.org/materials/mp-1735) | BaC2 | 3 | 3.872 | 5 月 | 1.522 |  |
| [mp-1029](https://next-gen.materialsproject.org/materials/mp-1029) | BaF2 | 3 | 10.621 | 5 月 | 6.525 |  |
| [mp-10616](https://next-gen.materialsproject.org/materials/mp-10616) | BaLiAs | 3 | 1.402 | 5 月 | 0.493 |  |
| [mp-10615](https://next-gen.materialsproject.org/materials/mp-10615) | BaLiP | 3 | 1.490 | 5 月 | 0.584 |  |
| [mp-1006878](https://next-gen.materialsproject.org/materials/mp-1006878) | BaO2 | 3 | 5.013 | 5 月 | 2.313 |  |
| [mp-1569](https://next-gen.materialsproject.org/materials/mp-1569) | Be2C | 3 | 2.095 | 5 月 | 1.198 |  |
| [mp-22965](https://next-gen.materialsproject.org/materials/mp-22965) | BiTeI | 3 | 1.880 | 5 月 | 1.188 |  |
| [mp-30068](https://next-gen.materialsproject.org/materials/mp-30068) | CIN | 3 | 7.710 | 5 月 | 4.321 |  |
| [mp-28240](https://next-gen.materialsproject.org/materials/mp-28240) | CSO | 3 | 8.442 | 5 月 | 4.750 |  |
| [mp-2482](https://next-gen.materialsproject.org/materials/mp-2482) | CaC2 | 3 | 3.468 | 5 月 | 1.218 |  |
| [mp-2741](https://next-gen.materialsproject.org/materials/mp-2741) | CaF2 | 3 | 11.041 | 5 月 | 6.930 |  |
| [mp-634859](https://next-gen.materialsproject.org/materials/mp-634859) | CaO2 | 3 | 6.255 | 5 月 | 2.717 |  |
| [mp-241](https://next-gen.materialsproject.org/materials/mp-241) | CdF2 | 3 | 6.418 | 5 月 | 2.900 |  |
| [mp-567259](https://next-gen.materialsproject.org/materials/mp-567259) | CdI2 | 3 | 3.986 | 5 月 | 2.184 |  |
| [mp-567677](https://next-gen.materialsproject.org/materials/mp-567677) | GeI2 | 3 | 3.115 | 5 月 | 2.016 |  |
| [mp-644272](https://next-gen.materialsproject.org/materials/mp-644272) | HCN | 3 | 11.108 | 5 月 | 6.458 |  |
| [mp-644385](https://next-gen.materialsproject.org/materials/mp-644385) | HCN | 3 | 11.304 | 5 月 | 6.682 |  |
| [mp-550893](https://next-gen.materialsproject.org/materials/mp-550893) | HfO2 | 3 | 5.679 | 5 月 | 3.563 |  |
| [mp-985829](https://next-gen.materialsproject.org/materials/mp-985829) | HfS2 | 3 | 2.471 | 5 月 | 0.982 |  |
| [mp-985831](https://next-gen.materialsproject.org/materials/mp-985831) | HfSe2 | 3 | 1.604 | 5 月 | 0.390 |  |
| [mp-11869](https://next-gen.materialsproject.org/materials/mp-11869) | HfSnPd | 3 | 0.440 | 5 月 | 0.396 |  |
| [mp-971](https://next-gen.materialsproject.org/materials/mp-971) | K2O | 3 | 3.836 | 5 月 | 1.737 |  |
| [mp-1062676](https://next-gen.materialsproject.org/materials/mp-1062676) | K2Pt | 3 | 3.046 | 5 月 | 1.321 |  |
| [mp-1022](https://next-gen.materialsproject.org/materials/mp-1022) | K2S | 3 | 4.344 | 5 月 | 2.138 |  |
| [mp-20134](https://next-gen.materialsproject.org/materials/mp-20134) | KCN | 3 | 8.686 | 5 月 | 4.813 |  |
| [mp-10161](https://next-gen.materialsproject.org/materials/mp-10161) | KZnSb | 3 | 1.396 | 5 月 | 0.281 |  |
| [mp-1960](https://next-gen.materialsproject.org/materials/mp-1960) | Li2O | 3 | 8.076 | 5 月 | 4.804 |  |
| [mp-1153](https://next-gen.materialsproject.org/materials/mp-1153) | Li2S | 3 | 5.254 | 5 月 | 3.178 |  |
| [mp-2286](https://next-gen.materialsproject.org/materials/mp-2286) | Li2Se | 3 | 4.575 | 5 月 | 2.809 |  |
| [mp-2530](https://next-gen.materialsproject.org/materials/mp-2530) | Li2Te | 3 | 3.815 | 5 月 | 2.389 |  |
| [mp-5920](https://next-gen.materialsproject.org/materials/mp-5920) | LiAlGe | 3 | 0.475 | 5 月 |  |  |
| [mp-3161](https://next-gen.materialsproject.org/materials/mp-3161) | LiAlSi | 3 | 0.481 | 5 月 | 0.026 |  |
| [mp-11390](https://next-gen.materialsproject.org/materials/mp-11390) | LiGaSi | 3 | 0.531 | 5 月 | 0.112 |  |
| [mp-12558](https://next-gen.materialsproject.org/materials/mp-12558) | LiMgAs | 3 | 2.286 | 5 月 | 1.293 |  |
| [mp-570213](https://next-gen.materialsproject.org/materials/mp-570213) | LiMgBi | 3 | 1.362 | 5 月 | 0.430 |  |
| [mp-9124](https://next-gen.materialsproject.org/materials/mp-9124) | LiZnAs | 3 | 1.685 | 5 月 | 0.509 |  |
| [mp-7575](https://next-gen.materialsproject.org/materials/mp-7575) | LiZnN | 3 | 2.008 | 5 月 | 0.749 |  |
| [mp-10182](https://next-gen.materialsproject.org/materials/mp-10182) | LiZnP | 3 | 1.931 | 5 月 | 1.183 |  |
| [mp-408](https://next-gen.materialsproject.org/materials/mp-408) | Mg2Ge | 3 | 0.648 | 5 月 | 0.140 |  |
| [mp-1367](https://next-gen.materialsproject.org/materials/mp-1367) | Mg2Si | 3 | 0.635 | 5 月 | 0.163 |  |
| [mp-30034](https://next-gen.materialsproject.org/materials/mp-30034) | MgBr2 | 3 | 6.955 | 5 月 | 4.124 |  |
| [mp-23205](https://next-gen.materialsproject.org/materials/mp-23205) | MgI2 | 3 | 5.531 | 5 月 | 3.445 |  |
| [mp-1434](https://next-gen.materialsproject.org/materials/mp-1434) | MoS2 | 3 | 1.712 | 5 月 | 1.176 |  |
| [mp-2352](https://next-gen.materialsproject.org/materials/mp-2352) | Na2O | 3 | 4.988 | 5 月 | 2.071 |  |
| [mp-648](https://next-gen.materialsproject.org/materials/mp-648) | Na2S | 3 | 4.688 | 5 月 | 2.364 |  |
| [mp-1266](https://next-gen.materialsproject.org/materials/mp-1266) | Na2Se | 3 | 4.170 | 5 月 | 1.957 |  |
| [mp-2784](https://next-gen.materialsproject.org/materials/mp-2784) | Na2Te | 3 | 3.809 | 5 月 | 1.921 |  |
| [mp-9437](https://next-gen.materialsproject.org/materials/mp-9437) | NbFeSb | 3 | 0.803 | 5 月 | 0.511 |  |
| [mp-510753](https://next-gen.materialsproject.org/materials/mp-510753) | NiO2 | 3 | 2.416 | 5 月 | 0.852 |  |
| [mp-315](https://next-gen.materialsproject.org/materials/mp-315) | PbF2 | 3 | 6.869 | 5 月 | 4.383 |  |
| [mp-22893](https://next-gen.materialsproject.org/materials/mp-22893) | PbI2 | 3 | 3.492 | 5 月 | 2.192 |  |
| [mp-617](https://next-gen.materialsproject.org/materials/mp-617) | PtO2 | 3 | 2.986 | 5 月 | 1.286 |  |
| [mp-762](https://next-gen.materialsproject.org/materials/mp-762) | PtS2 | 3 | 1.786 | 5 月 | 1.013 |  |
| [mp-1115](https://next-gen.materialsproject.org/materials/mp-1115) | PtSe2 | 3 | 0.866 | 5 月 | 0.355 |  |
| [mp-30459](https://next-gen.materialsproject.org/materials/mp-30459) | ScNiBi | 3 | 0.353 | 5 月 | 0.185 |  |
| [mp-3432](https://next-gen.materialsproject.org/materials/mp-3432) | ScNiSb | 3 | 0.496 | 5 月 | 0.270 |  |
| [mp-569779](https://next-gen.materialsproject.org/materials/mp-569779) | ScSbPd | 3 | 0.394 | 5 月 | 0.277 |  |
| [mp-7173](https://next-gen.materialsproject.org/materials/mp-7173) | ScSbPt | 3 | 0.809 | 5 月 | 0.622 |  |
| [mp-14](https://next-gen.materialsproject.org/materials/mp-14) | Se3 | 3 | 2.154 | 5 月 | 1.109 |  |
| [mp-1170](https://next-gen.materialsproject.org/materials/mp-1170) | SnS2 | 3 | 2.701 | 5 月 | 1.400 |  |
| [mp-665](https://next-gen.materialsproject.org/materials/mp-665) | SnSe2 | 3 | 1.484 | 5 月 | 0.574 |  |
| [mp-981](https://next-gen.materialsproject.org/materials/mp-981) | SrF2 | 3 | 10.638 | 5 月 | 6.727 |  |
| [mp-10614](https://next-gen.materialsproject.org/materials/mp-10614) | SrLiP | 3 | 1.583 | 5 月 | 0.713 |  |
| [mp-31454](https://next-gen.materialsproject.org/materials/mp-31454) | TaSbRu | 3 | 1.039 | 5 月 | 0.601 |  |
| [mp-19](https://next-gen.materialsproject.org/materials/mp-19) | Te3 | 3 | 0.906 | 5 月 | 0.562 |  |
| [mp-5967](https://next-gen.materialsproject.org/materials/mp-5967) | TiCoSb | 3 | 1.310 | 5 月 | 1.072 |  |
| [mp-1008680](https://next-gen.materialsproject.org/materials/mp-1008680) | TiGePt | 3 | 1.228 | 5 月 | 0.852 |  |
| [mp-924130](https://next-gen.materialsproject.org/materials/mp-924130) | TiNiSn | 3 | 0.589 | 5 月 | 0.458 |  |
| [mp-1062030](https://next-gen.materialsproject.org/materials/mp-1062030) | TiS2 | 3 | 2.321 | 5 月 | 0.010 | 構造: Pb–Ti–S の非整合層状化合物の「Ti S2-part」。c=11.25 Å（1T-TiS2 は 5.70）、1 原子 37.9 Å³（約 19）。5 月は 1 反復目に 0.65 → 2.49 eV と跳んだ |
| [mp-30847](https://next-gen.materialsproject.org/materials/mp-30847) | TiSnPt | 3 | 1.059 | 5 月 | 0.765 |  |
| [mp-567636](https://next-gen.materialsproject.org/materials/mp-567636) | VFeSb | 3 | 0.568 | 5 月 | 0.279 |  |
| [mp-31455](https://next-gen.materialsproject.org/materials/mp-31455) | VSbRu | 3 | 0.421 | 5 月 | 0.137 |  |
| [mp-28797](https://next-gen.materialsproject.org/materials/mp-28797) | YHSe | 3 | 2.510 | 5 月 | 1.394 |  |
| [mp-30460](https://next-gen.materialsproject.org/materials/mp-30460) | YNiBi | 3 | 0.449 | 5 月 | 0.200 |  |
| [mp-11520](https://next-gen.materialsproject.org/materials/mp-11520) | YNiSb | 3 | 0.602 | 5 月 | 0.289 |  |
| [mp-4964](https://next-gen.materialsproject.org/materials/mp-4964) | YSbPt | 3 | 0.617 | 5 月 | 0.062 |  |
| [mp-1008601](https://next-gen.materialsproject.org/materials/mp-1008601) | ZrAg2 | 3 | 0.045 | 5 月 | 0.006 |  |
| [mp-31451](https://next-gen.materialsproject.org/materials/mp-31451) | ZrCoBi | 3 | 1.153 | 5 月 | 1.000 |  |
| [mp-1565](https://next-gen.materialsproject.org/materials/mp-1565) | ZrO2 | 3 | 5.299 | 5 月 | 3.200 |  |
| [mp-1186](https://next-gen.materialsproject.org/materials/mp-1186) | ZrS2 | 3 | 2.171 | 5 月 | 0.924 |  |
| [mp-2076](https://next-gen.materialsproject.org/materials/mp-2076) | ZrSe2 | 3 | 1.385 | 5 月 | 0.375 |  |
| [mp-22689](https://next-gen.materialsproject.org/materials/mp-22689) | ZrSnPd | 3 | 0.256 | 5 月 | 0.067 |  |
| [mp-570301](https://next-gen.materialsproject.org/materials/mp-570301) | Ag2Br2 | 4 | 2.439 | 5 月 | 0.501 |  |
| [mp-570687](https://next-gen.materialsproject.org/materials/mp-570687) | Ag2Cl2 | 4 | 2.770 | 5 月 | 0.670 |  |
| [mp-567809](https://next-gen.materialsproject.org/materials/mp-567809) | Ag2I2 | 4 | 2.930 | 5 月 | 1.289 |  |
| [mp-568927](https://next-gen.materialsproject.org/materials/mp-568927) | Ag2I2 | 4 | 2.157 | 5 月 | 0.644 |  |
| [mp-29678](https://next-gen.materialsproject.org/materials/mp-29678) | AgBiS2 | 4 | 0.972 | 5 月 | 0.416 |  |
| [mp-661](https://next-gen.materialsproject.org/materials/mp-661) | Al2N2 | 4 | 6.067 | 5 月 | 4.087 |  |
| [mp-7885](https://next-gen.materialsproject.org/materials/mp-7885) | AlAgS2 | 4 | 2.866 | 5 月 | 1.521 |  |
| [mp-3748](https://next-gen.materialsproject.org/materials/mp-3748) | AlCuO2 | 4 | 3.336 | 5 月 | 2.030 |  |
| [mp-2793](https://next-gen.materialsproject.org/materials/mp-2793) | Au2Se2 | 4 | 0.627 | 5 月 | 0.135 |  |
| [mp-2653](https://next-gen.materialsproject.org/materials/mp-2653) | B2N2 | 4 | 7.490 | 5 月 | 5.386 |  |
| [mp-7991](https://next-gen.materialsproject.org/materials/mp-7991) | B2N2 | 4 | 6.543 | 5 月 | 3.803 |  |
| [mp-984](https://next-gen.materialsproject.org/materials/mp-984) | B2N2 | 4 | 6.543 | 5 月 | 4.328 |  |
| [mp-1008559](https://next-gen.materialsproject.org/materials/mp-1008559) | B2P2 | 4 | 1.964 | 5 月 | 1.135 |  |
| [mp-1018099](https://next-gen.materialsproject.org/materials/mp-1018099) | Ba2NCl | 4 | 2.221 | 5 月 | 0.963 |  |
| [mp-1018096](https://next-gen.materialsproject.org/materials/mp-1018096) | Ba2NF | 4 | 2.112 | 5 月 | 0.906 |  |
| [mp-1008500](https://next-gen.materialsproject.org/materials/mp-1008500) | Ba2O2 | 4 | 4.282 | 5 月 | 2.401 |  |
| [mp-7487](https://next-gen.materialsproject.org/materials/mp-7487) | Ba2O2 | 4 | 4.246 | 5 月 | 2.628 |  |
| [mp-27869](https://next-gen.materialsproject.org/materials/mp-27869) | Ba2PCl | 4 | 2.424 | 5 月 | 1.064 |  |
| [mp-1008505](https://next-gen.materialsproject.org/materials/mp-1008505) | Ba2Pd2 | 4 | 0.293 | 5 月 | 0.123 |  |
| [mp-571093](https://next-gen.materialsproject.org/materials/mp-571093) | BaAlSiH | 4 | 1.339 | 5 月 | 0.685 |  |
| [mp-1018095](https://next-gen.materialsproject.org/materials/mp-1018095) | BaGaGeH | 4 | 0.687 | 5 月 | 0.068 |  |
| [mp-1018094](https://next-gen.materialsproject.org/materials/mp-1018094) | BaGaSnH | 4 | 0.568 | 5 月 | 0.000 |  |
| [mp-2542](https://next-gen.materialsproject.org/materials/mp-2542) | Be2O2 | 4 | 10.963 | 5 月 | 7.661 |  |
| [mp-23301](https://next-gen.materialsproject.org/materials/mp-23301) | BiF3 | 4 | 6.400 | 5 月 | 3.992 |  |
| [mp-1066781](https://next-gen.materialsproject.org/materials/mp-1066781) | Br2Cl2 | 4 | 4.373 | 5 月 | 1.666 |  |
| [mp-23154](https://next-gen.materialsproject.org/materials/mp-23154) | Br4 | 4 | 3.718 | 5 月 | 1.212 |  |
| [mp-2516584](https://next-gen.materialsproject.org/materials/mp-2516584) | C4 | 4 | 2.400 | 5 月 | 1.448 |  |
| [mp-47](https://next-gen.materialsproject.org/materials/mp-47) | C4 | 4 | 5.240 | 5 月 | 3.580 |  |
| [mp-28554](https://next-gen.materialsproject.org/materials/mp-28554) | Ca2AsI | 4 | 2.843 | 5 月 | 1.453 |  |
| [mp-23009](https://next-gen.materialsproject.org/materials/mp-23009) | Ca2BrN | 4 | 3.867 | 5 月 | 2.305 |  |
| [mp-568177](https://next-gen.materialsproject.org/materials/mp-568177) | CaAlSiH | 4 | 1.031 | 5 月 | 0.264 |  |
| [mp-672](https://next-gen.materialsproject.org/materials/mp-672) | Cd2S2 | 4 | 2.424 | 5 月 | 0.895 |  |
| [mp-1070](https://next-gen.materialsproject.org/materials/mp-1070) | Cd2Se2 | 4 | 1.776 | 5 月 | 0.338 |  |
| [mp-12779](https://next-gen.materialsproject.org/materials/mp-12779) | Cd2Te2 | 4 | 1.723 | 5 月 | 0.441 |  |
| [mp-9146](https://next-gen.materialsproject.org/materials/mp-9146) | CdHgO2 | 4 | 2.023 | 5 月 | 0.315 |  |
| [mp-22848](https://next-gen.materialsproject.org/materials/mp-22848) | Cl4 | 4 | 5.794 | 5 月 | 2.093 |  |
| [mp-7896](https://next-gen.materialsproject.org/materials/mp-7896) | Cs2O2 | 4 | 3.947 | 5 月 | 1.529 |  |
| [mp-29266](https://next-gen.materialsproject.org/materials/mp-29266) | Cs2S2 | 4 | 3.900 | 5 月 | 1.496 |  |
| [mp-635413](https://next-gen.materialsproject.org/materials/mp-635413) | Cs3Bi | 4 | 1.067 | 5 月 | 0.110 |  |
| [mp-10378](https://next-gen.materialsproject.org/materials/mp-10378) | Cs3Sb | 4 | 1.656 | 5 月 | 0.480 |  |
| [mp-1065405](https://next-gen.materialsproject.org/materials/mp-1065405) | CsAuC2 | 4 | 4.037 | 5 月 | 1.985 |  |
| [mp-28650](https://next-gen.materialsproject.org/materials/mp-28650) | CsBr2F | 4 | 4.717 | 5 月 | 1.634 |  |
| [mp-22990](https://next-gen.materialsproject.org/materials/mp-22990) | CsICl2 | 4 | 4.979 | 5 月 | 2.094 |  |
| [mp-581024](https://next-gen.materialsproject.org/materials/mp-581024) | CsK2Sb | 4 | 2.023 | 5 月 | 0.821 |  |
| [mp-22917](https://next-gen.materialsproject.org/materials/mp-22917) | Cu2Br2 | 4 | 2.783 | 5 月 | 0.928 |  |
| [mp-24093](https://next-gen.materialsproject.org/materials/mp-24093) | Cu2H2 | 4 | 1.366 | 5 月 | 0.448 |  |
| [mp-1007922](https://next-gen.materialsproject.org/materials/mp-1007922) | Cu2I2 | 4 | 2.506 | 5 月 | 0.934 |  |
| [mp-22863](https://next-gen.materialsproject.org/materials/mp-22863) | Cu2I2 | 4 | 3.033 | 5 月 | 1.499 |  |
| [mp-569346](https://next-gen.materialsproject.org/materials/mp-569346) | Cu2I2 | 4 | 2.828 | 5 月 | 1.148 |  |
| [mp-1933](https://next-gen.materialsproject.org/materials/mp-1933) | Cu3N | 4 | 0.857 | 5 月 | 0.288 |  |
| [mp-29643](https://next-gen.materialsproject.org/materials/mp-29643) | CuAsSe2 | 4 | 0.947 | 5 月 | 0.662 |  |
| [mp-672285](https://next-gen.materialsproject.org/materials/mp-672285) | CuCSN | 4 | 3.291 | 5 月 | 2.069 |  |
| [mp-804](https://next-gen.materialsproject.org/materials/mp-804) | Ga2N2 | 4 | 3.291 | 5 月 | 1.914 |  |
| [mp-1018275](https://next-gen.materialsproject.org/materials/mp-1018275) | Ga2P2 | 4 | 1.417 | 5 月 | 0.718 |  |
| [mp-11342](https://next-gen.materialsproject.org/materials/mp-11342) | Ga2Se2 | 4 | 2.583 | 5 月 | 1.113 |  |
| [mp-4280](https://next-gen.materialsproject.org/materials/mp-4280) | GaCuO2 | 4 | 2.083 | 5 月 | 1.008 |  |
| [mp-632229](https://next-gen.materialsproject.org/materials/mp-632229) | H2Br2 | 4 | 7.694 | 5 月 | 4.334 |  |
| [mp-632326](https://next-gen.materialsproject.org/materials/mp-632326) | H2Cl2 | 4 | 9.265 | 5 月 | 5.365 |  |
| [mp-24504](https://next-gen.materialsproject.org/materials/mp-24504) | H4 | 4 | 15.709 | 5 月 | 8.967 |  |
| [mp-23177](https://next-gen.materialsproject.org/materials/mp-23177) | Hg2Br2 | 4 | 3.587 | 5 月 | 2.022 |  |
| [mp-22897](https://next-gen.materialsproject.org/materials/mp-22897) | Hg2Cl2 | 4 | 4.073 | 5 月 | 2.365 |  |
| [mp-22859](https://next-gen.materialsproject.org/materials/mp-22859) | Hg2I2 | 4 | 2.854 | 5 月 | 1.612 |  |
| [mp-1525634](https://next-gen.materialsproject.org/materials/mp-1525634) | I4 | 4 | 2.138 | 5 月 | 0.912 | MP の summary が無い（ID が消えたか統合された） |
| [mp-22870](https://next-gen.materialsproject.org/materials/mp-22870) | In2Br2 | 4 | 2.032 | 5 月 | 1.114 |  |
| [mp-571555](https://next-gen.materialsproject.org/materials/mp-571555) | In2Cl2 | 4 | 2.183 | 5 月 | 1.148 |  |
| [mp-23202](https://next-gen.materialsproject.org/materials/mp-23202) | In2I2 | 4 | 2.136 | 5 月 | 1.396 |  |
| [mp-21405](https://next-gen.materialsproject.org/materials/mp-21405) | In2Se2 | 4 | 1.890 | 5 月 | 0.966 |  |
| [mp-22691](https://next-gen.materialsproject.org/materials/mp-22691) | In2Se2 | 4 | 1.280 | 5 月 | 0.212 |  |
| [mp-20162](https://next-gen.materialsproject.org/materials/mp-20162) | InAgS2 | 4 | 0.728 | 5 月 | 0.052 |  |
| [mp-568516](https://next-gen.materialsproject.org/materials/mp-568516) | K3Bi | 4 | 0.948 | 5 月 | 0.095 |  |
| [mp-10159](https://next-gen.materialsproject.org/materials/mp-10159) | K3Sb | 4 | 1.724 | 5 月 | 0.608 |  |
| [mp-10102](https://next-gen.materialsproject.org/materials/mp-10102) | KAgC2 | 4 | 4.585 | 5 月 | 2.411 |  |
| [mp-27418](https://next-gen.materialsproject.org/materials/mp-27418) | KAuO2 | 4 | 3.462 | 5 月 | 0.875 |  |
| [mp-15724](https://next-gen.materialsproject.org/materials/mp-15724) | KNa2Sb | 4 | 1.814 | 5 月 | 0.626 |  |
| [mp-8188](https://next-gen.materialsproject.org/materials/mp-8188) | KScO2 | 4 | 6.217 | 5 月 | 3.517 |  |
| [mp-8409](https://next-gen.materialsproject.org/materials/mp-8409) | KYO2 | 4 | 6.051 | 5 月 | 3.897 |  |
| [mp-16238](https://next-gen.materialsproject.org/materials/mp-16238) | Li2AgSb | 4 | 0.828 | 5 月 | 0.069 |  |
| [mp-1378](https://next-gen.materialsproject.org/materials/mp-1378) | Li2C2 | 4 | 5.777 | 5 月 | 3.195 |  |
| [mp-15988](https://next-gen.materialsproject.org/materials/mp-15988) | Li2CuSb | 4 | 0.945 | 5 月 | 0.497 |  |
| [mp-568273](https://next-gen.materialsproject.org/materials/mp-568273) | Li2I2 | 4 | 6.230 | 5 月 | 4.286 |  |
| [mp-570935](https://next-gen.materialsproject.org/materials/mp-570935) | Li2I2 | 4 | 6.605 | 5 月 | 4.188 |  |
| [mp-23222](https://next-gen.materialsproject.org/materials/mp-23222) | Li3Bi | 4 | 0.997 | 5 月 | 0.444 |  |
| [mp-2251](https://next-gen.materialsproject.org/materials/mp-2251) | Li3N | 4 | 2.468 | 5 月 | 1.135 |  |
| [mp-2074](https://next-gen.materialsproject.org/materials/mp-2074) | Li3Sb | 4 | 1.488 | 5 月 | 0.744 |  |
| [mp-29868](https://next-gen.materialsproject.org/materials/mp-29868) | LiAuC2 | 4 | 2.997 | 5 月 | 1.443 |  |
| [mp-9158](https://next-gen.materialsproject.org/materials/mp-9158) | LiCuO2 | 4 | 2.319 | 5 月 | 0.317 |  |
| [mp-8002](https://next-gen.materialsproject.org/materials/mp-8002) | LiGaO2 | 4 | 4.747 | 5 月 | 3.711 |  |
| [mp-24199](https://next-gen.materialsproject.org/materials/mp-24199) | LiHF2 | 4 | 11.959 | 5 月 | 8.081 |  |
| [mp-10618](https://next-gen.materialsproject.org/materials/mp-10618) | LiInSe2 | 4 | 1.695 | 5 月 | 0.665 |  |
| [mp-2659](https://next-gen.materialsproject.org/materials/mp-2659) | LiN3 | 4 | 6.495 | 5 月 | 3.602 |  |
| [mp-14115](https://next-gen.materialsproject.org/materials/mp-14115) | LiRhO2 | 4 | 3.136 | 5 月 | 1.310 |  |
| [mp-15788](https://next-gen.materialsproject.org/materials/mp-15788) | LiYS2 | 4 | 2.991 | 5 月 | 1.755 |  |
| [mp-1018040](https://next-gen.materialsproject.org/materials/mp-1018040) | Mg2Se2 | 4 | 4.468 | 5 月 | 2.433 |  |
| [mp-1039](https://next-gen.materialsproject.org/materials/mp-1039) | Mg2Te2 | 4 | 3.921 | 5 月 | 2.213 |  |
| [mp-570747](https://next-gen.materialsproject.org/materials/mp-570747) | N4 | 4 | 13.238 | 5 月 | 6.774 |  |
| [mp-672233](https://next-gen.materialsproject.org/materials/mp-672233) | N4 | 4 | 13.718 | 5 月 | 6.931 |  |
| [mp-11180](https://next-gen.materialsproject.org/materials/mp-11180) | Na2C2 | 4 | 6.236 | 5 月 | 3.370 |  |
| [mp-723629](https://next-gen.materialsproject.org/materials/mp-723629) | Na2I2 | 4 | 5.423 | 5 月 | 2.863 |  |
| [mp-8860](https://next-gen.materialsproject.org/materials/mp-8860) | Na3As | 4 | 1.150 | 5 月 | 0.122 |  |
| [mp-4541](https://next-gen.materialsproject.org/materials/mp-4541) | NaCuO2 | 4 | 2.927 | 5 月 | 0.340 |  |
| [mp-5077](https://next-gen.materialsproject.org/materials/mp-5077) | NaLi2Sb | 4 | 1.424 | 5 月 | 0.624 |  |
| [mp-1066400](https://next-gen.materialsproject.org/materials/mp-1066400) | NaN3 | 4 | 6.712 | 5 月 | 3.986 |  |
| [mp-570538](https://next-gen.materialsproject.org/materials/mp-570538) | NaN3 | 4 | 7.222 | 5 月 | 3.826 |  |
| [mp-2964](https://next-gen.materialsproject.org/materials/mp-2964) | NaNO2 | 4 | 6.964 | 5 月 | 2.304 |  |
| [mp-7017](https://next-gen.materialsproject.org/materials/mp-7017) | NaNbN2 | 4 | 1.755 | 5 月 | 0.880 |  |
| [mp-7914](https://next-gen.materialsproject.org/materials/mp-7914) | NaScO2 | 4 | 6.530 | 5 月 | 3.901 |  |
| [mp-5475](https://next-gen.materialsproject.org/materials/mp-5475) | NaTaN2 | 4 | 2.299 | 5 月 | 1.251 |  |
| [mp-157](https://next-gen.materialsproject.org/materials/mp-157) | P4 | 4 | 0.747 | 5 月 | 0.015 |  |
| [mp-19921](https://next-gen.materialsproject.org/materials/mp-19921) | Pb2O2 | 4 | 2.414 | 5 月 | 1.440 |  |
| [mp-1018115](https://next-gen.materialsproject.org/materials/mp-1018115) | Pb2S2 | 4 | 1.619 | 5 月 | 0.989 |  |
| [mp-726184](https://next-gen.materialsproject.org/materials/mp-726184) | Pb2S2 | 4 | 3.737 | 5 月 | 1.730 | 構造: Pb–Ti–S の非整合層状化合物の PbS 層。a=4.22、c=16.6 Å、1 原子 73.6 Å³（岩塩型 PbS 26.7、最近接 2.99 Å・配位 6 に対し 2.66 Å・配位 5） |
| [mp-727323](https://next-gen.materialsproject.org/materials/mp-727323) | Pb2S2 | 4 | 3.231 | 5 月 | 1.729 | 構造: Pb–Ti–S の非整合層状化合物の「Pb S-part」。a=4.22、c=11.2 Å、1 原子 49.7 Å³（岩塩型 PbS 26.7） |
| [mp-288](https://next-gen.materialsproject.org/materials/mp-288) | Pt2S2 | 4 | 1.199 | 5 月 | 0.210 |  |
| [mp-7895](https://next-gen.materialsproject.org/materials/mp-7895) | Rb2O2 | 4 | 4.534 | 5 月 | 1.624 |  |
| [mp-558071](https://next-gen.materialsproject.org/materials/mp-558071) | Rb2S2 | 4 | 4.063 | 5 月 | 1.463 |  |
| [mp-10421](https://next-gen.materialsproject.org/materials/mp-10421) | RbAuC2 | 4 | 4.079 | 5 月 | 2.106 |  |
| [mp-5841](https://next-gen.materialsproject.org/materials/mp-5841) | RbCuC2 | 4 | 3.877 | 5 月 | 1.873 |  |
| [mp-8145](https://next-gen.materialsproject.org/materials/mp-8145) | RbScO2 | 4 | 5.414 | 5 月 | 3.192 |  |
| [mp-8176](https://next-gen.materialsproject.org/materials/mp-8176) | RbTlO2 | 4 | 1.796 | 5 月 | 0.870 |  |
| [mp-999265](https://next-gen.materialsproject.org/materials/mp-999265) | RbYS2 | 4 | 3.569 | 5 月 | 2.130 |  |
| [mp-973185](https://next-gen.materialsproject.org/materials/mp-973185) | ScAgO2 | 4 | 3.516 | 5 月 | 2.129 |  |
| [mp-12908](https://next-gen.materialsproject.org/materials/mp-12908) | ScAgSe2 | 4 | 1.714 | 5 月 | 0.543 |  |
| [mp-6980](https://next-gen.materialsproject.org/materials/mp-6980) | ScCuS2 | 4 | 1.970 | 5 月 | 0.784 |  |
| [mp-13313](https://next-gen.materialsproject.org/materials/mp-13313) | ScTlSe2 | 4 | 1.169 | 5 月 | 0.483 |  |
| [mp-13314](https://next-gen.materialsproject.org/materials/mp-13314) | ScTlTe2 | 4 | 0.683 | 5 月 | 0.108 |  |
| [mp-7140](https://next-gen.materialsproject.org/materials/mp-7140) | Si2C2 | 4 | 3.528 | 5 月 | 2.440 |  |
| [mp-29803](https://next-gen.materialsproject.org/materials/mp-29803) | Si2H2 | 4 | 3.784 | 5 月 | 2.033 |  |
| [mp-165](https://next-gen.materialsproject.org/materials/mp-165) | Si4 | 4 | 0.978 | 5 月 | 0.320 |  |
| [mp-2097](https://next-gen.materialsproject.org/materials/mp-2097) | Sn2O2 | 4 | 0.782 | 5 月 | 0.245 |  |
| [mp-545552](https://next-gen.materialsproject.org/materials/mp-545552) | Sn2O2 | 4 | 2.375 | 5 月 | 1.392 |  |
| [mp-1379](https://next-gen.materialsproject.org/materials/mp-1379) | Sn2S2 | 4 | 1.917 | 5 月 | 0.925 |  |
| [mp-559676](https://next-gen.materialsproject.org/materials/mp-559676) | Sn2S2 | 4 | 0.847 | 5 月 | 0.365 |  |
| [mp-727322](https://next-gen.materialsproject.org/materials/mp-727322) | Sn2S2 | 4 | 2.691 | 5 月 | 1.394 | 構造: Sn–Ti–S の非整合層状化合物の SnS 層。a=4.09、c=11.0 Å、1 原子 45.3 Å³ |
| [mp-8781](https://next-gen.materialsproject.org/materials/mp-8781) | Sn2S2 | 4 | 3.030 | 5 月 | 1.491 | 構造: Sn–Nb–S の非整合層状化合物の SnS 層。a=4.12、c=11.0 Å、1 原子 46.8 Å³ |
| [mp-2168](https://next-gen.materialsproject.org/materials/mp-2168) | Sn2Se2 | 4 | 0.641 | 5 月 | 0.291 |  |
| [mp-690794](https://next-gen.materialsproject.org/materials/mp-690794) | Sr2HN | 4 | 2.797 | 5 月 | 1.719 |  |
| [mp-569677](https://next-gen.materialsproject.org/materials/mp-569677) | Sr2IN | 4 | 2.796 | 5 月 | 1.764 |  |
| [mp-23033](https://next-gen.materialsproject.org/materials/mp-23033) | Sr2NCl | 4 | 3.185 | 5 月 | 1.773 |  |
| [mp-980057](https://next-gen.materialsproject.org/materials/mp-980057) | SrAlGeH | 4 | 1.440 | 5 月 | 0.640 |  |
| [mp-570485](https://next-gen.materialsproject.org/materials/mp-570485) | SrAlSiH | 4 | 1.236 | 5 月 | 0.532 |  |
| [mp-978847](https://next-gen.materialsproject.org/materials/mp-978847) | SrGaGeH | 4 | 1.128 | 5 月 | 0.290 |  |
| [mp-9382](https://next-gen.materialsproject.org/materials/mp-9382) | SrZrN2 | 4 | 1.149 | 5 月 | 0.215 |  |
| [mp-8927](https://next-gen.materialsproject.org/materials/mp-8927) | TaCuN2 | 4 | 1.203 | 5 月 | 0.464 |  |
| [mp-999135](https://next-gen.materialsproject.org/materials/mp-999135) | TaCuN2 | 4 | 1.020 | 5 月 | 0.409 |  |
| [mp-19963](https://next-gen.materialsproject.org/materials/mp-19963) | TiFe2Sn | 4 | 0.058 | 5 月 |  |  |
| [mp-568949](https://next-gen.materialsproject.org/materials/mp-568949) | Tl2Br2 | 4 | 3.510 | 5 月 | 2.053 |  |
| [mp-571079](https://next-gen.materialsproject.org/materials/mp-571079) | Tl2Cl2 | 4 | 3.889 | 5 月 | 2.345 |  |
| [mp-557835](https://next-gen.materialsproject.org/materials/mp-557835) | Tl2F2 | 4 | 4.725 | 5 月 | 2.960 |  |
| [mp-22858](https://next-gen.materialsproject.org/materials/mp-22858) | Tl2I2 | 4 | 3.182 | 5 月 | 2.026 |  |
| [mp-554310](https://next-gen.materialsproject.org/materials/mp-554310) | TlBiS2 | 4 | 0.810 | 5 月 | 0.350 |  |
| [mp-29662](https://next-gen.materialsproject.org/materials/mp-29662) | TlBiSe2 | 4 | 0.454 | 5 月 | 0.107 |  |
| [mp-27438](https://next-gen.materialsproject.org/materials/mp-27438) | TlBiTe2 | 4 | 0.642 | 5 月 | 0.358 |  |
| [mp-22566](https://next-gen.materialsproject.org/materials/mp-22566) | TlInS2 | 4 | 1.451 | 5 月 | 0.479 |  |
| [mp-19390](https://next-gen.materialsproject.org/materials/mp-19390) | WO3 | 4 | 1.533 | 5 月 | 0.388 |  |
| [mp-1067744](https://next-gen.materialsproject.org/materials/mp-1067744) | YTlSe2 | 4 | 1.879 | 5 月 | 1.157 |  |
| [mp-1065514](https://next-gen.materialsproject.org/materials/mp-1065514) | YTlTe2 | 4 | 1.237 | 5 月 | 0.571 |  |
| [mp-2133](https://next-gen.materialsproject.org/materials/mp-2133) | Zn2O2 | 4 | 3.176 | 5 月 | 0.750 |  |
| [mp-560588](https://next-gen.materialsproject.org/materials/mp-560588) | Zn2S2 | 4 | 3.777 | 5 月 | 1.986 |  |
| [mp-380](https://next-gen.materialsproject.org/materials/mp-380) | Zn2Se2 | 4 | 2.771 | 5 月 | 1.100 |  |
| [mp-22995](https://next-gen.materialsproject.org/materials/mp-22995) | Ag3SI | 5 | 1.534 | 5 月 | 0.396 |  |
| [mp-8579](https://next-gen.materialsproject.org/materials/mp-8579) | BaAg2S2 | 5 | 1.976 | 5 月 | 0.732 |  |
| [mp-8281](https://next-gen.materialsproject.org/materials/mp-8281) | BaCd2As2 | 5 | 0.849 | 5 月 | 0.037 |  |
| [mp-8279](https://next-gen.materialsproject.org/materials/mp-8279) | BaCd2P2 | 5 | 1.307 | 5 月 | 0.531 |  |
| [mp-10437](https://next-gen.materialsproject.org/materials/mp-10437) | BaCu2Se2 | 5 | 1.322 | 5 月 | 0.394 |  |
| [mp-23818](https://next-gen.materialsproject.org/materials/mp-23818) | BaLiH3 | 5 | 4.329 | 5 月 | 2.012 |  |
| [mp-632808](https://next-gen.materialsproject.org/materials/mp-632808) | BaLiH3 | 5 | 3.459 | 5 月 | 1.071 |  |
| [mp-8280](https://next-gen.materialsproject.org/materials/mp-8280) | BaMg2As2 | 5 | 2.005 | 5 月 | 1.012 |  |
| [mp-29209](https://next-gen.materialsproject.org/materials/mp-29209) | BaMg2Bi2 | 5 | 1.252 | 5 月 | 0.409 |  |
| [mp-8278](https://next-gen.materialsproject.org/materials/mp-8278) | BaMg2P2 | 5 | 2.095 | 5 月 | 1.052 |  |
| [mp-9567](https://next-gen.materialsproject.org/materials/mp-9567) | BaMg2Sb2 | 5 | 1.672 | 5 月 | 0.831 |  |
| [mp-3163](https://next-gen.materialsproject.org/materials/mp-3163) | BaSnO3 | 5 | 2.348 | 5 月 | 0.509 |  |
| [mp-5020](https://next-gen.materialsproject.org/materials/mp-5020) | BaTiO3 | 5 | 4.204 | 5 月 | 2.258 |  |
| [mp-5777](https://next-gen.materialsproject.org/materials/mp-5777) | BaTiO3 | 5 | 3.973 | 5 月 | 2.052 |  |
| [mp-5986](https://next-gen.materialsproject.org/materials/mp-5986) | BaTiO3 | 5 | 3.387 | 5 月 | 1.677 |  |
| [mp-3834](https://next-gen.materialsproject.org/materials/mp-3834) | BaZrO3 | 5 | 5.330 | 5 月 | 2.984 |  |
| [mp-1017552](https://next-gen.materialsproject.org/materials/mp-1017552) | Bi2O3 | 5 | 1.866 | 5 月 | 1.113 |  |
| [mp-27910](https://next-gen.materialsproject.org/materials/mp-27910) | Bi2Te2S | 5 | 0.819 | 5 月 | 0.422 |  |
| [mp-29666](https://next-gen.materialsproject.org/materials/mp-29666) | Bi2Te2Se | 5 | 0.786 | 5 月 | 0.384 |  |
| [mp-34202](https://next-gen.materialsproject.org/materials/mp-34202) | Bi2Te3 | 5 | 0.868 | 5 月 | 0.500 |  |
| [mp-31406](https://next-gen.materialsproject.org/materials/mp-31406) | Bi2TeSe2 | 5 | 1.403 | 5 月 | 0.812 |  |
| [mp-7223](https://next-gen.materialsproject.org/materials/mp-7223) | Ca3AsN | 5 | 1.809 | 5 月 | 0.650 |  |
| [mp-31149](https://next-gen.materialsproject.org/materials/mp-31149) | Ca3BiN | 5 | 1.295 | 5 月 | 0.356 |  |
| [mp-13146](https://next-gen.materialsproject.org/materials/mp-13146) | Ca3N2 | 5 | 2.400 | 5 月 | 1.380 |  |
| [mp-11824](https://next-gen.materialsproject.org/materials/mp-11824) | Ca3PN | 5 | 1.902 | 5 月 | 0.718 |  |
| [mp-1018630](https://next-gen.materialsproject.org/materials/mp-1018630) | CaBe2As2 | 5 | 1.488 | 5 月 | 0.606 |  |
| [mp-1017548](https://next-gen.materialsproject.org/materials/mp-1017548) | CaBe2P2 | 5 | 1.663 | 5 月 | 0.651 |  |
| [mp-7067](https://next-gen.materialsproject.org/materials/mp-7067) | CaCd2As2 | 5 | 0.938 | 5 月 | 0.070 |  |
| [mp-9570](https://next-gen.materialsproject.org/materials/mp-9570) | CaCd2P2 | 5 | 1.504 | 5 月 | 0.664 |  |
| [mp-7430](https://next-gen.materialsproject.org/materials/mp-7430) | CaCd2Sb2 | 5 | 0.801 | 5 月 | 0.177 |  |
| [mp-1070099](https://next-gen.materialsproject.org/materials/mp-1070099) | CaGa2As2 | 5 | 0.375 | 5 月 | 0.113 |  |
| [mp-23879](https://next-gen.materialsproject.org/materials/mp-23879) | CaH2O2 | 5 | 7.330 | 5 月 | 4.022 |  |
| [mp-9564](https://next-gen.materialsproject.org/materials/mp-9564) | CaMg2As2 | 5 | 2.251 | 5 月 | 1.162 |  |
| [mp-29208](https://next-gen.materialsproject.org/materials/mp-29208) | CaMg2Bi2 | 5 | 1.288 | 5 月 | 0.435 |  |
| [mp-5795](https://next-gen.materialsproject.org/materials/mp-5795) | CaMg2N2 | 5 | 3.593 | 5 月 | 2.181 |  |
| [mp-9565](https://next-gen.materialsproject.org/materials/mp-9565) | CaMg2Sb2 | 5 | 1.625 | 5 月 | 0.758 |  |
| [mp-5893](https://next-gen.materialsproject.org/materials/mp-5893) | CaSiO3 | 5 | 6.187 | 5 月 | 3.356 |  |
| [mp-7986](https://next-gen.materialsproject.org/materials/mp-7986) | CaSnO3 | 5 | 3.923 | 5 月 | 1.128 |  |
| [mp-5827](https://next-gen.materialsproject.org/materials/mp-5827) | CaTiO3 | 5 | 3.820 | 5 月 | 1.718 |  |
| [mp-9571](https://next-gen.materialsproject.org/materials/mp-9571) | CaZn2As2 | 5 | 1.396 | 5 月 | 0.292 |  |
| [mp-9569](https://next-gen.materialsproject.org/materials/mp-9569) | CaZn2P2 | 5 | 1.693 | 5 月 | 0.707 |  |
| [mp-542112](https://next-gen.materialsproject.org/materials/mp-542112) | CaZrO3 | 5 | 5.436 | 5 月 | 3.042 |  |
| [mp-505824](https://next-gen.materialsproject.org/materials/mp-505824) | Cs2PdC2 | 5 | 3.120 | 5 月 | 1.703 |  |
| [mp-505825](https://next-gen.materialsproject.org/materials/mp-505825) | Cs2PtC2 | 5 | 2.434 | 5 月 | 1.196 |  |
| [mp-30056](https://next-gen.materialsproject.org/materials/mp-30056) | CsCaBr3 | 5 | 6.979 | 5 月 | 4.157 |  |
| [mp-7104](https://next-gen.materialsproject.org/materials/mp-7104) | CsCaF3 | 5 | 11.348 | 5 月 | 6.911 |  |
| [mp-644203](https://next-gen.materialsproject.org/materials/mp-644203) | CsCaH3 | 5 | 5.928 | 5 月 | 2.814 |  |
| [mp-568544](https://next-gen.materialsproject.org/materials/mp-568544) | CsCdCl3 | 5 | 4.212 | 5 月 | 1.458 |  |
| [mp-8399](https://next-gen.materialsproject.org/materials/mp-8399) | CsCdF3 | 5 | 6.983 | 5 月 | 3.155 |  |
| [mp-1068340](https://next-gen.materialsproject.org/materials/mp-1068340) | CsGeBr3 | 5 | 2.559 | 5 月 | 1.133 |  |
| [mp-570223](https://next-gen.materialsproject.org/materials/mp-570223) | CsGeBr3 | 5 | 1.445 | 5 月 | 0.563 |  |
| [mp-22988](https://next-gen.materialsproject.org/materials/mp-22988) | CsGeCl3 | 5 | 3.844 | 5 月 | 1.977 |  |
| [mp-28377](https://next-gen.materialsproject.org/materials/mp-28377) | CsGeI3 | 5 | 1.963 | 5 月 | 1.058 |  |
| [mp-8401](https://next-gen.materialsproject.org/materials/mp-8401) | CsMgF3 | 5 | 10.821 | 5 月 | 6.473 |  |
| [mp-600089](https://next-gen.materialsproject.org/materials/mp-600089) | CsPbBr3 | 5 | 3.103 | 5 月 | 1.574 |  |
| [mp-5811](https://next-gen.materialsproject.org/materials/mp-5811) | CsPbF3 | 5 | 5.077 | 5 月 | 2.693 |  |
| [mp-1069538](https://next-gen.materialsproject.org/materials/mp-1069538) | CsPbI3 | 5 | 2.287 | 5 月 | 1.235 |  |
| [mp-27214](https://next-gen.materialsproject.org/materials/mp-27214) | CsSnBr3 | 5 | 1.473 | 5 月 | 0.474 |  |
| [mp-1070375](https://next-gen.materialsproject.org/materials/mp-1070375) | CsSnCl3 | 5 | 1.994 | 5 月 | 0.738 |  |
| [mp-614013](https://next-gen.materialsproject.org/materials/mp-614013) | CsSnI3 | 5 | 1.240 | 5 月 | 0.475 |  |
| [mp-8397](https://next-gen.materialsproject.org/materials/mp-8397) | CsSrF3 | 5 | 10.306 | 5 月 | 6.140 |  |
| [mp-1017567](https://next-gen.materialsproject.org/materials/mp-1017567) | Hf2SN2 | 5 | 2.079 | 5 月 | 0.940 |  |
| [mp-1069533](https://next-gen.materialsproject.org/materials/mp-1069533) | Hf2SN2 | 5 | 3.246 | 5 月 | 1.913 |  |
| [mp-1017565](https://next-gen.materialsproject.org/materials/mp-1017565) | In2Se3 | 5 | 0.935 | 5 月 | 0.077 |  |
| [mp-10408](https://next-gen.materialsproject.org/materials/mp-10408) | K2CN2 | 5 | 5.608 | 5 月 | 3.124 |  |
| [mp-1068977](https://next-gen.materialsproject.org/materials/mp-1068977) | K2PdC2 | 5 | 2.859 | 5 月 | 1.402 |  |
| [mp-540584](https://next-gen.materialsproject.org/materials/mp-540584) | K2PdO2 | 5 | 3.782 | 5 月 | 1.303 |  |
| [mp-1068941](https://next-gen.materialsproject.org/materials/mp-1068941) | K2PdS2 | 5 | 2.621 | 5 月 | 0.897 |  |
| [mp-1069706](https://next-gen.materialsproject.org/materials/mp-1069706) | K2PdSe2 | 5 | 2.048 | 5 月 | 0.639 |  |
| [mp-976876](https://next-gen.materialsproject.org/materials/mp-976876) | K2PtC2 | 5 | 2.120 | 5 月 | 0.893 |  |
| [mp-7928](https://next-gen.materialsproject.org/materials/mp-7928) | K2PtS2 | 5 | 2.877 | 5 月 | 1.018 |  |
| [mp-8621](https://next-gen.materialsproject.org/materials/mp-8621) | K2PtSe2 | 5 | 2.414 | 5 月 | 0.740 |  |
| [mp-1068813](https://next-gen.materialsproject.org/materials/mp-1068813) | K2Te2Pd | 5 | 1.589 | 5 月 | 0.605 |  |
| [mp-8623](https://next-gen.materialsproject.org/materials/mp-8623) | K2Te2Pt | 5 | 1.855 | 5 月 | 0.603 |  |
| [mp-9200](https://next-gen.materialsproject.org/materials/mp-9200) | K3AuO | 5 | 2.167 | 5 月 | 0.445 |  |
| [mp-28166](https://next-gen.materialsproject.org/materials/mp-28166) | K3BrO | 5 | 3.219 | 5 月 | 0.876 |  |
| [mp-28171](https://next-gen.materialsproject.org/materials/mp-28171) | K3IO | 5 | 3.167 | 5 月 | 1.005 |  |
| [mp-22958](https://next-gen.materialsproject.org/materials/mp-22958) | KBrO3 | 5 | 6.975 | 5 月 | 3.624 |  |
| [mp-10175](https://next-gen.materialsproject.org/materials/mp-10175) | KCdF3 | 5 | 6.899 | 5 月 | 2.911 |  |
| [mp-23737](https://next-gen.materialsproject.org/materials/mp-23737) | KMgH3 | 5 | 5.292 | 5 月 | 2.160 |  |
| [mp-6920](https://next-gen.materialsproject.org/materials/mp-6920) | KNO3 | 5 | 7.369 | 5 月 | 2.850 |  |
| [mp-4342](https://next-gen.materialsproject.org/materials/mp-4342) | KNbO3 | 5 | 2.884 | 5 月 | 1.170 |  |
| [mp-5246](https://next-gen.materialsproject.org/materials/mp-5246) | KNbO3 | 5 | 3.460 | 5 月 | 1.615 |  |
| [mp-7375](https://next-gen.materialsproject.org/materials/mp-7375) | KNbO3 | 5 | 3.669 | 5 月 | 1.879 |  |
| [mp-935811](https://next-gen.materialsproject.org/materials/mp-935811) | KNbO3 | 5 | 2.662 | 5 月 | 1.228 |  |
| [mp-5878](https://next-gen.materialsproject.org/materials/mp-5878) | KZnF3 | 5 | 8.192 | 5 月 | 3.706 |  |
| [mp-7608](https://next-gen.materialsproject.org/materials/mp-7608) | Li2PdO2 | 5 | 3.538 | 5 月 | 1.087 |  |
| [mp-3216](https://next-gen.materialsproject.org/materials/mp-3216) | Li2ZrN2 | 5 | 2.990 | 5 月 | 1.680 |  |
| [mp-2646](https://next-gen.materialsproject.org/materials/mp-2646) | Mg3Sb2 | 5 | 0.855 | 5 月 | 0.135 |  |
| [mp-9514](https://next-gen.materialsproject.org/materials/mp-9514) | MgAl2C2 | 5 | 2.573 | 5 月 | 1.767 |  |
| [mp-865185](https://next-gen.materialsproject.org/materials/mp-865185) | MgBe2As2 | 5 | 1.376 | 5 月 | 0.572 |  |
| [mp-11917](https://next-gen.materialsproject.org/materials/mp-11917) | MgBe2N2 | 5 | 5.874 | 5 月 | 4.101 |  |
| [mp-1017628](https://next-gen.materialsproject.org/materials/mp-1017628) | MgBe2P2 | 5 | 1.736 | 5 月 | 0.757 |  |
| [mp-4823](https://next-gen.materialsproject.org/materials/mp-4823) | Na2PdC2 | 5 | 2.115 | 5 月 | 0.848 |  |
| [mp-4366](https://next-gen.materialsproject.org/materials/mp-4366) | Na2PtC2 | 5 | 1.367 | 5 月 | 0.307 |  |
| [mp-22313](https://next-gen.materialsproject.org/materials/mp-22313) | Na2PtO2 | 5 | 3.806 | 5 月 | 1.368 |  |
| [mp-28602](https://next-gen.materialsproject.org/materials/mp-28602) | Na3ClO | 5 | 5.014 | 5 月 | 2.320 |  |
| [mp-3136](https://next-gen.materialsproject.org/materials/mp-3136) | NaNbO3 | 5 | 2.836 | 5 月 | 1.282 |  |
| [mp-4170](https://next-gen.materialsproject.org/materials/mp-4170) | NaTaO3 | 5 | 4.181 | 5 月 | 2.117 |  |
| [mp-1069193](https://next-gen.materialsproject.org/materials/mp-1069193) | PdN2Cl2 | 5 | 1.229 | 5 月 | 0.128 |  |
| [mp-10918](https://next-gen.materialsproject.org/materials/mp-10918) | Rb2PdC2 | 5 | 2.901 | 5 月 | 1.551 |  |
| [mp-1069771](https://next-gen.materialsproject.org/materials/mp-1069771) | Rb2PdSe2 | 5 | 2.188 | 5 月 | 0.730 |  |
| [mp-10919](https://next-gen.materialsproject.org/materials/mp-10919) | Rb2PtC2 | 5 | 2.359 | 5 月 | 1.049 |  |
| [mp-7929](https://next-gen.materialsproject.org/materials/mp-7929) | Rb2PtS2 | 5 | 2.756 | 5 月 | 1.028 |  |
| [mp-8622](https://next-gen.materialsproject.org/materials/mp-8622) | Rb2PtSe2 | 5 | 2.405 | 5 月 | 0.794 |  |
| [mp-1070356](https://next-gen.materialsproject.org/materials/mp-1070356) | Rb2Te2Pd | 5 | 1.408 | 5 月 | 0.549 |  |
| [mp-1070498](https://next-gen.materialsproject.org/materials/mp-1070498) | Rb2Te2Pt | 5 | 2.093 | 5 月 | 0.761 |  |
| [mp-4405](https://next-gen.materialsproject.org/materials/mp-4405) | Rb3AuO | 5 | 1.847 | 5 月 | 0.070 |  |
| [mp-30055](https://next-gen.materialsproject.org/materials/mp-30055) | Rb3BrO | 5 | 2.580 | 5 月 | 0.346 |  |
| [mp-28872](https://next-gen.materialsproject.org/materials/mp-28872) | RbBrO3 | 5 | 6.920 | 5 月 | 3.805 |  |
| [mp-3654](https://next-gen.materialsproject.org/materials/mp-3654) | RbCaF3 | 5 | 10.590 | 5 月 | 6.323 |  |
| [mp-23949](https://next-gen.materialsproject.org/materials/mp-23949) | RbCaH3 | 5 | 6.243 | 5 月 | 2.923 |  |
| [mp-571458](https://next-gen.materialsproject.org/materials/mp-571458) | RbGeI3 | 5 | 1.112 | 5 月 | 0.500 |  |
| [mp-7482](https://next-gen.materialsproject.org/materials/mp-7482) | RbHgF3 | 5 | 3.425 | 5 月 | 0.414 |  |
| [mp-27193](https://next-gen.materialsproject.org/materials/mp-27193) | RbIO3 | 5 | 4.955 | 5 月 | 2.758 |  |
| [mp-8402](https://next-gen.materialsproject.org/materials/mp-8402) | RbMgF3 | 5 | 12.007 | 5 月 | 7.112 |  |
| [mp-21043](https://next-gen.materialsproject.org/materials/mp-21043) | RbPbF3 | 5 | 5.717 | 5 月 | 2.979 |  |
| [mp-3525](https://next-gen.materialsproject.org/materials/mp-3525) | Sb2Te2Se | 5 | 0.355 | 5 月 | 0.093 |  |
| [mp-1201](https://next-gen.materialsproject.org/materials/mp-1201) | Sb2Te3 | 5 | 0.502 | 5 月 | 0.236 |  |
| [mp-571550](https://next-gen.materialsproject.org/materials/mp-571550) | Sb2TeSe2 | 5 | 0.810 | 5 月 | 0.501 |  |
| [mp-8612](https://next-gen.materialsproject.org/materials/mp-8612) | Sb2TeSe2 | 5 | 0.733 | 5 月 | 0.434 |  |
| [mp-1818](https://next-gen.materialsproject.org/materials/mp-1818) | SiF4 | 5 | 13.021 | 5 月 | 7.589 |  |
| [mp-571436](https://next-gen.materialsproject.org/materials/mp-571436) | SnI4 | 5 | 3.359 | 5 月 | 1.096 |  |
| [mp-1069819](https://next-gen.materialsproject.org/materials/mp-1069819) | SnN2F2 | 5 | 3.609 | 5 月 | 1.000 |  |
| [mp-9306](https://next-gen.materialsproject.org/materials/mp-9306) | Sr2ZnN2 | 5 | 1.634 | 5 月 | 0.714 |  |
| [mp-570008](https://next-gen.materialsproject.org/materials/mp-570008) | Sr3BiN | 5 | 1.014 | 5 月 | 0.263 |  |
| [mp-7752](https://next-gen.materialsproject.org/materials/mp-7752) | Sr3SbN | 5 | 1.064 | 5 月 | 0.254 |  |
| [mp-1068011](https://next-gen.materialsproject.org/materials/mp-1068011) | SrAlH3 | 5 | 1.197 | 5 月 | 0.563 |  |
| [mp-7771](https://next-gen.materialsproject.org/materials/mp-7771) | SrCd2As2 | 5 | 0.875 | 5 月 | 0.000 |  |
| [mp-8277](https://next-gen.materialsproject.org/materials/mp-8277) | SrCd2P2 | 5 | 1.402 | 5 月 | 0.550 |  |
| [mp-4551](https://next-gen.materialsproject.org/materials/mp-4551) | SrHfO3 | 5 | 6.141 | 5 月 | 3.563 |  |
| [mp-24419](https://next-gen.materialsproject.org/materials/mp-24419) | SrLiH3 | 5 | 4.354 | 5 月 | 1.637 |  |
| [mp-863260](https://next-gen.materialsproject.org/materials/mp-863260) | SrMg2As2 | 5 | 2.229 | 5 月 | 1.260 |  |
| [mp-1069695](https://next-gen.materialsproject.org/materials/mp-1069695) | SrMg2Bi2 | 5 | 1.340 | 5 月 | 0.401 |  |
| [mp-10550](https://next-gen.materialsproject.org/materials/mp-10550) | SrMg2N2 | 5 | 2.957 | 5 月 | 1.692 |  |
| [mp-9566](https://next-gen.materialsproject.org/materials/mp-9566) | SrMg2Sb2 | 5 | 1.773 | 5 月 | 0.933 |  |
| [mp-1017439](https://next-gen.materialsproject.org/materials/mp-1017439) | SrSiO3 | 5 | 4.766 | 5 月 | 2.321 |  |
| [mp-546973](https://next-gen.materialsproject.org/materials/mp-546973) | SrSnO3 | 5 | 3.091 | 5 月 | 0.708 |  |
| [mp-5229](https://next-gen.materialsproject.org/materials/mp-5229) | SrTiO3 | 5 | 3.597 | 5 月 | 1.720 |  |
| [mp-7770](https://next-gen.materialsproject.org/materials/mp-7770) | SrZn2As2 | 5 | 1.243 | 5 月 | 0.163 |  |
| [mp-8276](https://next-gen.materialsproject.org/materials/mp-8276) | SrZn2P2 | 5 | 1.558 | 5 月 | 0.680 |  |
| [mp-19845](https://next-gen.materialsproject.org/materials/mp-19845) | TiPbO3 | 5 | 2.722 | 5 月 | 1.557 |  |
| [mp-20459](https://next-gen.materialsproject.org/materials/mp-20459) | TiPbO3 | 5 | 3.221 | 5 月 | 1.623 |  |
| [mp-553875](https://next-gen.materialsproject.org/materials/mp-553875) | Zr2SN2 | 5 | 1.724 | 5 月 | 0.703 |  |
| [mp-20337](https://next-gen.materialsproject.org/materials/mp-20337) | ZrPbO3 | 5 | 3.964 | 5 月 | 2.632 |  |
| [mp-27863](https://next-gen.materialsproject.org/materials/mp-27863) | Al2Cl2O2 | 6 | 8.935 | 5 月 | 5.674 |  |
| [mp-10910](https://next-gen.materialsproject.org/materials/mp-10910) | Al4Ru2 | 6 | 0.462 | 5 月 | 0.090 |  |
| [mp-7848](https://next-gen.materialsproject.org/materials/mp-7848) | AlPO4 | 6 | 9.618 | 5 月 | 5.623 |  |
| [mp-2455](https://next-gen.materialsproject.org/materials/mp-2455) | As4Os2 | 6 | 1.002 | 5 月 | 0.631 |  |
| [mp-766](https://next-gen.materialsproject.org/materials/mp-766) | As4Ru2 | 6 | 0.757 | 5 月 | 0.457 |  |
| [mp-947](https://next-gen.materialsproject.org/materials/mp-947) | Au4S2 | 6 | 3.200 | 5 月 | 2.001 |  |
| [mp-3277](https://next-gen.materialsproject.org/materials/mp-3277) | BAsO4 | 6 | 8.260 | 5 月 | 4.535 |  |
| [mp-3589](https://next-gen.materialsproject.org/materials/mp-3589) | BPO4 | 6 | 11.012 | 5 月 | 7.372 |  |
| [mp-27724](https://next-gen.materialsproject.org/materials/mp-27724) | BPS4 | 6 | 4.128 | 5 月 | 1.955 |  |
| [mp-9899](https://next-gen.materialsproject.org/materials/mp-9899) | Ba2Ag2P2 | 6 | 0.804 | 5 月 | 0.057 |  |
| [mp-1205316](https://next-gen.materialsproject.org/materials/mp-1205316) | Ba2Ag2Sb2 | 6 | 0.691 | 5 月 | 0.050 |  |
| [mp-23070](https://next-gen.materialsproject.org/materials/mp-23070) | Ba2Br2F2 | 6 | 7.430 | 5 月 | 4.665 |  |
| [mp-10293](https://next-gen.materialsproject.org/materials/mp-10293) | Ba2C4 | 6 | 4.393 | 5 月 | 2.017 |  |
| [mp-23432](https://next-gen.materialsproject.org/materials/mp-23432) | Ba2Cl2F2 | 6 | 8.347 | 5 月 | 5.177 |  |
| [mp-8644](https://next-gen.materialsproject.org/materials/mp-8644) | Ba2F4 | 6 | 9.773 | 5 月 | 5.774 |  |
| [mp-24424](https://next-gen.materialsproject.org/materials/mp-24424) | Ba2H2Br2 | 6 | 5.356 | 5 月 | 2.953 |  |
| [mp-23862](https://next-gen.materialsproject.org/materials/mp-23862) | Ba2H2I2 | 6 | 4.654 | 5 月 | 2.604 |  |
| [mp-1018651](https://next-gen.materialsproject.org/materials/mp-1018651) | Ba2H3I | 6 | 4.108 | 5 月 | 2.023 |  |
| [mp-22951](https://next-gen.materialsproject.org/materials/mp-22951) | Ba2I2F2 | 6 | 6.112 | 5 月 | 3.848 |  |
| [mp-13277](https://next-gen.materialsproject.org/materials/mp-13277) | Ba2Li2P2 | 6 | 1.890 | 5 月 | 0.910 |  |
| [mp-10485](https://next-gen.materialsproject.org/materials/mp-10485) | Ba2Li2Sb2 | 6 | 1.618 | 5 月 | 0.700 |  |
| [mp-684](https://next-gen.materialsproject.org/materials/mp-684) | Ba2S4 | 6 | 3.719 | 5 月 | 1.520 |  |
| [mp-7547](https://next-gen.materialsproject.org/materials/mp-7547) | Ba2Se4 | 6 | 2.413 | 5 月 | 0.635 |  |
| [mp-2150](https://next-gen.materialsproject.org/materials/mp-2150) | Ba2Te4 | 6 | 1.480 | 5 月 | 0.358 |  |
| [mp-30139](https://next-gen.materialsproject.org/materials/mp-30139) | Be2Br4 | 6 | 7.891 | 5 月 | 4.935 |  |
| [mp-23267](https://next-gen.materialsproject.org/materials/mp-23267) | Be2Cl4 | 6 | 9.431 | 5 月 | 6.095 |  |
| [mp-570886](https://next-gen.materialsproject.org/materials/mp-570886) | Be2I4 | 6 | 6.137 | 5 月 | 3.872 |  |
| [mp-1072399](https://next-gen.materialsproject.org/materials/mp-1072399) | Be5Pt | 6 | 0.355 | 5 月 | 0.139 |  |
| [mp-5046](https://next-gen.materialsproject.org/materials/mp-5046) | BeSO4 | 6 | 10.600 | 5 月 | 7.007 |  |
| [mp-23072](https://next-gen.materialsproject.org/materials/mp-23072) | Bi2Br2O2 | 6 | 3.264 | 5 月 | 2.134 |  |
| [mp-22939](https://next-gen.materialsproject.org/materials/mp-22939) | Bi2Cl2O2 | 6 | 3.728 | 5 月 | 2.510 |  |
| [mp-22987](https://next-gen.materialsproject.org/materials/mp-22987) | Bi2I2O2 | 6 | 2.557 | 5 月 | 1.415 |  |
| [mp-23074](https://next-gen.materialsproject.org/materials/mp-23074) | Bi2O2F2 | 6 | 4.999 | 5 月 | 3.416 |  |
| [mp-28944](https://next-gen.materialsproject.org/materials/mp-28944) | Bi2Te2Cl2 | 6 | 2.277 | 5 月 | 1.293 |  |
| [mp-27989](https://next-gen.materialsproject.org/materials/mp-27989) | C2Br2N2 | 6 | 8.516 | 5 月 | 4.829 |  |
| [mp-27502](https://next-gen.materialsproject.org/materials/mp-27502) | C2N2Cl2 | 6 | 9.170 | 5 月 | 5.569 |  |
| [mp-1077906](https://next-gen.materialsproject.org/materials/mp-1077906) | C2O4 | 6 | 10.089 | 5 月 | 6.274 |  |
| [mp-2232](https://next-gen.materialsproject.org/materials/mp-2232) | C2S4 | 6 | 6.623 | 5 月 | 3.312 |  |
| [mp-1071625](https://next-gen.materialsproject.org/materials/mp-1071625) | C2Se4 | 6 | 4.579 | 5 月 | 2.278 |  |
| [mp-22888](https://next-gen.materialsproject.org/materials/mp-22888) | Ca2Br4 | 6 | 7.146 | 5 月 | 4.318 |  |
| [mp-1575](https://next-gen.materialsproject.org/materials/mp-1575) | Ca2C4 | 6 | 4.290 | 5 月 | 1.901 |  |
| [mp-917](https://next-gen.materialsproject.org/materials/mp-917) | Ca2C4 | 6 | 4.531 | 5 月 | 2.152 |  |
| [mp-27546](https://next-gen.materialsproject.org/materials/mp-27546) | Ca2Cl2F2 | 6 | 9.075 | 5 月 | 5.442 |  |
| [mp-24422](https://next-gen.materialsproject.org/materials/mp-24422) | Ca2H2Br2 | 6 | 6.632 | 5 月 | 3.882 |  |
| [mp-23859](https://next-gen.materialsproject.org/materials/mp-23859) | Ca2H2Cl2 | 6 | 6.441 | 5 月 | 3.524 |  |
| [mp-24204](https://next-gen.materialsproject.org/materials/mp-24204) | Ca2H2I2 | 6 | 5.406 | 5 月 | 3.278 |  |
| [mp-1205322](https://next-gen.materialsproject.org/materials/mp-1205322) | Ca2H4 | 6 | 2.792 | 5 月 | 0.413 |  |
| [mp-24809](https://next-gen.materialsproject.org/materials/mp-24809) | Ca2H4 | 6 | 2.765 | 5 月 | 0.393 |  |
| [mp-471](https://next-gen.materialsproject.org/materials/mp-471) | Cd2As4 | 6 | 0.768 | 5 月 | 0.102 |  |
| [mp-27934](https://next-gen.materialsproject.org/materials/mp-27934) | Cd2Br4 | 6 | 4.957 | 5 月 | 2.616 |  |
| [mp-28248](https://next-gen.materialsproject.org/materials/mp-28248) | Cd2I4 | 6 | 3.942 | 5 月 | 2.116 |  |
| [mp-1492](https://next-gen.materialsproject.org/materials/mp-1492) | Cd3Te3 | 6 | 0.728 | 5 月 | 0.003 |  |
| [mp-1077383](https://next-gen.materialsproject.org/materials/mp-1077383) | Cs2Au2S2 | 6 | 3.509 | 5 月 | 1.699 |  |
| [mp-574599](https://next-gen.materialsproject.org/materials/mp-574599) | Cs2Au2Se2 | 6 | 3.263 | 5 月 | 1.561 |  |
| [mp-22935](https://next-gen.materialsproject.org/materials/mp-22935) | Cs2Br2F2 | 6 | 5.400 | 5 月 | 2.400 |  |
| [mp-541037](https://next-gen.materialsproject.org/materials/mp-541037) | Cs2Cu2O2 | 6 | 2.861 | 5 月 | 1.058 |  |
| [mp-6973](https://next-gen.materialsproject.org/materials/mp-6973) | Cs2Na2S2 | 6 | 4.320 | 5 月 | 2.261 |  |
| [mp-8658](https://next-gen.materialsproject.org/materials/mp-8658) | Cs2Na2Se2 | 6 | 3.851 | 5 月 | 1.882 |  |
| [mp-5339](https://next-gen.materialsproject.org/materials/mp-5339) | Cs2Na2Te2 | 6 | 3.574 | 5 月 | 1.789 |  |
| [mp-1541909](https://next-gen.materialsproject.org/materials/mp-1541909) | Cs2Te2Au2 | 6 | 2.722 | 5 月 | 1.236 |  |
| [mp-573755](https://next-gen.materialsproject.org/materials/mp-573755) | Cs2Te2Au2 | 6 | 2.344 | 5 月 | 1.038 |  |
| [mp-13548](https://next-gen.materialsproject.org/materials/mp-13548) | Cs4Pt2 | 6 | 3.103 | 5 月 | 1.681 |  |
| [mp-1077814](https://next-gen.materialsproject.org/materials/mp-1077814) | Cs4S2 | 6 | 3.783 | 5 月 | 1.843 |  |
| [mp-9384](https://next-gen.materialsproject.org/materials/mp-9384) | CsAu3S2 | 6 | 3.441 | 5 月 | 2.030 |  |
| [mp-9386](https://next-gen.materialsproject.org/materials/mp-9386) | CsAu3Se2 | 6 | 3.318 | 5 月 | 1.940 |  |
| [mp-553303](https://next-gen.materialsproject.org/materials/mp-553303) | CsCu3O2 | 6 | 3.388 | 5 月 | 1.463 |  |
| [mp-7786](https://next-gen.materialsproject.org/materials/mp-7786) | CsCu3S2 | 6 | 3.070 | 5 月 | 1.627 |  |
| [mp-1077811](https://next-gen.materialsproject.org/materials/mp-1077811) | Cu2Ag2S2 | 6 | 1.378 | 5 月 | 0.607 |  |
| [mp-8911](https://next-gen.materialsproject.org/materials/mp-8911) | Cu2Ag2S2 | 6 | 1.001 | 5 月 | 0.339 |  |
| [mp-1072589](https://next-gen.materialsproject.org/materials/mp-1072589) | Cu2GeS3 | 6 | 0.477 | 5 月 | 0.229 |  |
| [mp-361](https://next-gen.materialsproject.org/materials/mp-361) | Cu4O2 | 6 | 1.922 | 5 月 | 0.571 |  |
| [mp-2008](https://next-gen.materialsproject.org/materials/mp-2008) | Fe2As4 | 6 | 0.420 | 5 月 | 0.222 |  |
| [mp-20027](https://next-gen.materialsproject.org/materials/mp-20027) | Fe2P4 | 6 | 0.636 | 5 月 | 0.358 |  |
| [mp-1522](https://next-gen.materialsproject.org/materials/mp-1522) | Fe2S4 | 6 | 1.294 | 5 月 | 0.745 |  |
| [mp-760](https://next-gen.materialsproject.org/materials/mp-760) | Fe2Se4 | 6 | 0.679 | 5 月 | 0.203 |  |
| [mp-19880](https://next-gen.materialsproject.org/materials/mp-19880) | Fe2Te4 | 6 | 0.387 | 5 月 | 0.049 |  |
| [mp-8211](https://next-gen.materialsproject.org/materials/mp-8211) | Ga2Ge2Te2 | 6 | 1.061 | 5 月 | 0.108 |  |
| [mp-570875](https://next-gen.materialsproject.org/materials/mp-570875) | Ga4Os2 | 6 | 1.061 | 5 月 | 0.652 |  |
| [mp-1072429](https://next-gen.materialsproject.org/materials/mp-1072429) | Ga4Ru2 | 6 | 0.443 | 5 月 | 0.133 |  |
| [mp-29403](https://next-gen.materialsproject.org/materials/mp-29403) | GaCuI4 | 6 | 3.506 | 5 月 | 1.727 |  |
| [mp-1072104](https://next-gen.materialsproject.org/materials/mp-1072104) | Ge2O4 | 6 | 3.849 | 5 月 | 1.220 |  |
| [mp-470](https://next-gen.materialsproject.org/materials/mp-470) | Ge2O4 | 6 | 4.246 | 5 月 | 1.571 |  |
| [mp-1071032](https://next-gen.materialsproject.org/materials/mp-1071032) | Ge2S4 | 6 | 2.817 | 5 月 | 0.894 |  |
| [mp-7582](https://next-gen.materialsproject.org/materials/mp-7582) | Ge2S4 | 6 | 4.166 | 5 月 | 2.448 |  |
| [mp-10074](https://next-gen.materialsproject.org/materials/mp-10074) | Ge2Se4 | 6 | 2.875 | 5 月 | 1.517 |  |
| [mp-568346](https://next-gen.materialsproject.org/materials/mp-568346) | Hf2Br2N2 | 6 | 3.830 | 5 月 | 1.936 |  |
| [mp-541911](https://next-gen.materialsproject.org/materials/mp-541911) | Hf2N2Cl2 | 6 | 3.956 | 5 月 | 2.206 |  |
| [mp-1018721](https://next-gen.materialsproject.org/materials/mp-1018721) | Hf2O4 | 6 | 6.865 | 5 月 | 4.517 |  |
| [mp-23292](https://next-gen.materialsproject.org/materials/mp-23292) | Hg2Br4 | 6 | 4.490 | 5 月 | 2.199 |  |
| [mp-570172](https://next-gen.materialsproject.org/materials/mp-570172) | Hg2I2Br2 | 6 | 4.109 | 5 月 | 2.054 |  |
| [mp-23173](https://next-gen.materialsproject.org/materials/mp-23173) | Hg2I4 | 6 | 3.849 | 5 月 | 2.008 |  |
| [mp-23192](https://next-gen.materialsproject.org/materials/mp-23192) | Hg2I4 | 6 | 2.599 | 5 月 | 0.996 |  |
| [mp-1077107](https://next-gen.materialsproject.org/materials/mp-1077107) | Hg3O3 | 6 | 2.599 | 5 月 | 1.199 |  |
| [mp-7826](https://next-gen.materialsproject.org/materials/mp-7826) | Hg3O3 | 6 | 2.585 | 5 月 | 1.203 |  |
| [mp-9252](https://next-gen.materialsproject.org/materials/mp-9252) | Hg3S3 | 6 | 2.803 | 5 月 | 1.533 |  |
| [mp-1018722](https://next-gen.materialsproject.org/materials/mp-1018722) | Hg3Se3 | 6 | 1.808 | 5 月 | 0.943 |  |
| [mp-358](https://next-gen.materialsproject.org/materials/mp-358) | Hg3Te3 | 6 | 1.117 | 5 月 | 0.551 |  |
| [mp-1077098](https://next-gen.materialsproject.org/materials/mp-1077098) | Hg6 | 6 | 1.794 | 5 月 | 1.094 |  |
| [mp-27703](https://next-gen.materialsproject.org/materials/mp-27703) | In2Br2O2 | 6 | 3.969 | 5 月 | 1.958 |  |
| [mp-27702](https://next-gen.materialsproject.org/materials/mp-27702) | In2Cl2O2 | 6 | 4.507 | 5 月 | 2.427 |  |
| [mp-20790](https://next-gen.materialsproject.org/materials/mp-20790) | InPS4 | 6 | 3.850 | 5 月 | 2.189 |  |
| [mp-16236](https://next-gen.materialsproject.org/materials/mp-16236) | K2Ag2Se2 | 6 | 1.654 | 5 月 | 0.305 |  |
| [mp-7077](https://next-gen.materialsproject.org/materials/mp-7077) | K2Au2S2 | 6 | 3.669 | 5 月 | 1.948 |  |
| [mp-9881](https://next-gen.materialsproject.org/materials/mp-9881) | K2Au2Se2 | 6 | 3.383 | 5 月 | 1.706 |  |
| [mp-1077078](https://next-gen.materialsproject.org/materials/mp-1077078) | K2C2N2 | 6 | 8.964 | 5 月 | 4.868 |  |
| [mp-1018736](https://next-gen.materialsproject.org/materials/mp-1018736) | K2Cd2As2 | 6 | 0.898 | 5 月 | 0.023 |  |
| [mp-7435](https://next-gen.materialsproject.org/materials/mp-7435) | K2Cu2Se2 | 6 | 1.773 | 5 月 | 0.143 |  |
| [mp-7436](https://next-gen.materialsproject.org/materials/mp-7436) | K2Cu2Te2 | 6 | 2.016 | 5 月 | 0.544 |  |
| [mp-634676](https://next-gen.materialsproject.org/materials/mp-634676) | K2H2S2 | 6 | 5.672 | 5 月 | 2.950 |  |
| [mp-722181](https://next-gen.materialsproject.org/materials/mp-722181) | K2Li2S2 | 6 | 4.813 | 5 月 | 2.631 |  |
| [mp-8430](https://next-gen.materialsproject.org/materials/mp-8430) | K2Li2S2 | 6 | 5.079 | 5 月 | 2.907 |  |
| [mp-8756](https://next-gen.materialsproject.org/materials/mp-8756) | K2Li2Se2 | 6 | 4.422 | 5 月 | 2.411 |  |
| [mp-4495](https://next-gen.materialsproject.org/materials/mp-4495) | K2Li2Te2 | 6 | 4.020 | 5 月 | 2.289 |  |
| [mp-1019089](https://next-gen.materialsproject.org/materials/mp-1019089) | K2Mg2As2 | 6 | 2.343 | 5 月 | 1.050 |  |
| [mp-1019105](https://next-gen.materialsproject.org/materials/mp-1019105) | K2Mg2Bi2 | 6 | 1.290 | 5 月 | 0.367 |  |
| [mp-1018737](https://next-gen.materialsproject.org/materials/mp-1018737) | K2Mg2P2 | 6 | 2.851 | 5 月 | 1.551 |  |
| [mp-7089](https://next-gen.materialsproject.org/materials/mp-7089) | K2Mg2Sb2 | 6 | 2.343 | 5 月 | 1.177 |  |
| [mp-6948](https://next-gen.materialsproject.org/materials/mp-6948) | K2Na2O2 | 6 | 4.957 | 5 月 | 2.317 |  |
| [mp-29055](https://next-gen.materialsproject.org/materials/mp-29055) | K2Sb4 | 6 | 0.721 | 5 月 | 0.127 |  |
| [mp-3481](https://next-gen.materialsproject.org/materials/mp-3481) | K2Sn2As2 | 6 | 1.021 | 5 月 | 0.342 |  |
| [mp-3486](https://next-gen.materialsproject.org/materials/mp-3486) | K2Sn2Sb2 | 6 | 0.851 | 5 月 | 0.305 |  |
| [mp-7421](https://next-gen.materialsproject.org/materials/mp-7421) | K2Zn2As2 | 6 | 1.536 | 5 月 | 0.230 |  |
| [mp-7437](https://next-gen.materialsproject.org/materials/mp-7437) | K2Zn2P2 | 6 | 1.975 | 5 月 | 0.757 |  |
| [mp-7438](https://next-gen.materialsproject.org/materials/mp-7438) | K2Zn2Sb2 | 6 | 1.666 | 5 月 | 0.394 |  |
| [mp-559278](https://next-gen.materialsproject.org/materials/mp-559278) | K4S2 | 6 | 3.838 | 5 月 | 1.778 |  |
| [mp-5347](https://next-gen.materialsproject.org/materials/mp-5347) | KAlF4 | 6 | 11.813 | 5 月 | 6.783 |  |
| [mp-546633](https://next-gen.materialsproject.org/materials/mp-546633) | KClO4 | 6 | 8.928 | 5 月 | 5.218 |  |
| [mp-1019888](https://next-gen.materialsproject.org/materials/mp-1019888) | KNa2BN2 | 6 | 3.865 | 5 月 | 1.641 |  |
| [mp-545359](https://next-gen.materialsproject.org/materials/mp-545359) | KNa2CuO2 | 6 | 3.764 | 5 月 | 1.382 |  |
| [mp-9244](https://next-gen.materialsproject.org/materials/mp-9244) | Li2B2C2 | 6 | 1.980 | 5 月 | 0.967 |  |
| [mp-9562](https://next-gen.materialsproject.org/materials/mp-9562) | Li2Be2As2 | 6 | 1.770 | 5 月 | 0.949 |  |
| [mp-9915](https://next-gen.materialsproject.org/materials/mp-9915) | Li2Be2P2 | 6 | 2.004 | 5 月 | 1.070 |  |
| [mp-9575](https://next-gen.materialsproject.org/materials/mp-9575) | Li2Be2Sb2 | 6 | 1.417 | 5 月 | 0.951 |  |
| [mp-9919](https://next-gen.materialsproject.org/materials/mp-9919) | Li2Zn2Sb2 | 6 | 1.028 | 5 月 | 0.352 |  |
| [mp-29149](https://next-gen.materialsproject.org/materials/mp-29149) | Li4NCl | 6 | 3.362 | 5 月 | 1.768 |  |
| [mp-29771](https://next-gen.materialsproject.org/materials/mp-29771) | Mg2C4 | 6 | 4.663 | 5 月 | 2.518 |  |
| [mp-1072956](https://next-gen.materialsproject.org/materials/mp-1072956) | Mg2F4 | 6 | 12.419 | 5 月 | 6.633 |  |
| [mp-23710](https://next-gen.materialsproject.org/materials/mp-23710) | Mg2H4 | 6 | 6.542 | 5 月 | 3.446 |  |
| [mp-2815](https://next-gen.materialsproject.org/materials/mp-2815) | Mo2S4 | 6 | 2.189 | 5 月 | 1.276 |  |
| [mp-1018807](https://next-gen.materialsproject.org/materials/mp-1018807) | Mo2Se4 | 6 | 1.965 | 5 月 | 1.177 |  |
| [mp-1634](https://next-gen.materialsproject.org/materials/mp-1634) | Mo2Se4 | 6 | 1.748 | 5 月 | 1.065 |  |
| [mp-9573](https://next-gen.materialsproject.org/materials/mp-9573) | Na2Be2As2 | 6 | 2.288 | 5 月 | 1.243 |  |
| [mp-9574](https://next-gen.materialsproject.org/materials/mp-9574) | Na2Be2Sb2 | 6 | 1.753 | 5 月 | 0.907 |  |
| [mp-7433](https://next-gen.materialsproject.org/materials/mp-7433) | Na2Cu2Se2 | 6 | 1.315 | 5 月 | 0.008 |  |
| [mp-7434](https://next-gen.materialsproject.org/materials/mp-7434) | Na2Cu2Te2 | 6 | 1.752 | 5 月 | 0.579 |  |
| [mp-634657](https://next-gen.materialsproject.org/materials/mp-634657) | Na2H2S2 | 6 | 5.720 | 5 月 | 2.971 |  |
| [mp-8452](https://next-gen.materialsproject.org/materials/mp-8452) | Na2Li2S2 | 6 | 5.165 | 5 月 | 2.999 |  |
| [mp-5962](https://next-gen.materialsproject.org/materials/mp-5962) | Na2Mg2As2 | 6 | 2.192 | 5 月 | 0.927 |  |
| [mp-7090](https://next-gen.materialsproject.org/materials/mp-7090) | Na2Mg2Sb2 | 6 | 1.945 | 5 月 | 1.073 |  |
| [mp-12433](https://next-gen.materialsproject.org/materials/mp-12433) | Na2Sn2N2 | 6 | 1.814 | 5 月 | 0.905 |  |
| [mp-29529](https://next-gen.materialsproject.org/materials/mp-29529) | Na2Sn2P2 | 6 | 1.251 | 5 月 | 0.446 |  |
| [mp-13097](https://next-gen.materialsproject.org/materials/mp-13097) | Na2Zn2As2 | 6 | 1.376 | 5 月 | 0.313 |  |
| [mp-4824](https://next-gen.materialsproject.org/materials/mp-4824) | Na2Zn2P2 | 6 | 1.974 | 5 月 | 0.717 |  |
| [mp-8374](https://next-gen.materialsproject.org/materials/mp-8374) | Na4S2 | 6 | 3.526 | 5 月 | 1.464 |  |
| [mp-486](https://next-gen.materialsproject.org/materials/mp-486) | Ni2P4 | 6 | 0.494 | 5 月 | 0.336 |  |
| [mp-29443](https://next-gen.materialsproject.org/materials/mp-29443) | P2I4 | 6 | 3.363 | 5 月 | 1.585 |  |
| [mp-2319](https://next-gen.materialsproject.org/materials/mp-2319) | P4Os2 | 6 | 1.092 | 5 月 | 0.772 |  |
| [mp-28266](https://next-gen.materialsproject.org/materials/mp-28266) | P4Pd2 | 6 | 0.759 | 5 月 | 0.346 |  |
| [mp-1413](https://next-gen.materialsproject.org/materials/mp-1413) | P4Ru2 | 6 | 0.822 | 5 月 | 0.464 |  |
| [mp-23008](https://next-gen.materialsproject.org/materials/mp-23008) | Pb2Br2F2 | 6 | 4.016 | 5 月 | 2.369 |  |
| [mp-22964](https://next-gen.materialsproject.org/materials/mp-22964) | Pb2Cl2F2 | 6 | 5.188 | 5 月 | 3.251 |  |
| [mp-22969](https://next-gen.materialsproject.org/materials/mp-22969) | Pb2I2F2 | 6 | 3.374 | 5 月 | 1.911 |  |
| [mp-540789](https://next-gen.materialsproject.org/materials/mp-540789) | Pb2I4 | 6 | 3.568 | 5 月 | 2.242 |  |
| [mp-567503](https://next-gen.materialsproject.org/materials/mp-567503) | Pb2I4 | 6 | 3.582 | 5 月 | 2.204 |  |
| [mp-569595](https://next-gen.materialsproject.org/materials/mp-569595) | Pb2I4 | 6 | 3.645 | 5 月 | 2.250 |  |
| [mp-1018888](https://next-gen.materialsproject.org/materials/mp-1018888) | Pd2Cl4 | 6 | 3.604 | 5 月 | 0.603 |  |
| [mp-1018891](https://next-gen.materialsproject.org/materials/mp-1018891) | Pd2Cl4 | 6 | 2.416 | 5 月 | 0.403 |  |
| [mp-569008](https://next-gen.materialsproject.org/materials/mp-569008) | Pd2Cl4 | 6 | 3.211 | 5 月 | 0.531 |  |
| [mp-1285](https://next-gen.materialsproject.org/materials/mp-1285) | Pt2O4 | 6 | 1.679 | 5 月 | 0.487 |  |
| [mp-7868](https://next-gen.materialsproject.org/materials/mp-7868) | Pt2O4 | 6 | 3.401 | 5 月 | 1.312 |  |
| [mp-9010](https://next-gen.materialsproject.org/materials/mp-9010) | Rb2Au2S2 | 6 | 3.637 | 5 月 | 1.816 |  |
| [mp-9731](https://next-gen.materialsproject.org/materials/mp-9731) | Rb2Au2Se2 | 6 | 3.327 | 5 月 | 1.607 |  |
| [mp-9846](https://next-gen.materialsproject.org/materials/mp-9846) | Rb2Ca2Sb2 | 6 | 2.576 | 5 月 | 1.333 |  |
| [mp-696804](https://next-gen.materialsproject.org/materials/mp-696804) | Rb2H2S2 | 6 | 5.610 | 5 月 | 2.767 |  |
| [mp-8751](https://next-gen.materialsproject.org/materials/mp-8751) | Rb2Li2S2 | 6 | 4.671 | 5 月 | 2.565 |  |
| [mp-9250](https://next-gen.materialsproject.org/materials/mp-9250) | Rb2Li2Se2 | 6 | 4.138 | 5 月 | 2.139 |  |
| [mp-8799](https://next-gen.materialsproject.org/materials/mp-8799) | Rb2Na2S2 | 6 | 4.252 | 5 月 | 2.131 |  |
| [mp-568018](https://next-gen.materialsproject.org/materials/mp-568018) | Rb2Sb4 | 6 | 1.042 | 5 月 | 0.196 |  |
| [mp-9008](https://next-gen.materialsproject.org/materials/mp-9008) | Rb2Te2Au2 | 6 | 1.947 | 5 月 | 0.814 |  |
| [mp-721469](https://next-gen.materialsproject.org/materials/mp-721469) | Rb4S2 | 6 | 3.616 | 5 月 | 1.601 |  |
| [mp-9385](https://next-gen.materialsproject.org/materials/mp-9385) | RbAu3Se2 | 6 | 2.675 | 5 月 | 1.476 |  |
| [mp-550759](https://next-gen.materialsproject.org/materials/mp-550759) | RbClO4 | 6 | 8.967 | 5 月 | 5.325 |  |
| [mp-568100](https://next-gen.materialsproject.org/materials/mp-568100) | ReNCl4 | 6 | 2.903 | 5 月 | 0.964 |  |
| [mp-28051](https://next-gen.materialsproject.org/materials/mp-28051) | Sb2Te2I2 | 6 | 1.411 | 5 月 | 0.767 |  |
| [mp-2695](https://next-gen.materialsproject.org/materials/mp-2695) | Sb4Os2 | 6 | 0.638 | 5 月 | 0.312 |  |
| [mp-20928](https://next-gen.materialsproject.org/materials/mp-20928) | Sb4Ru2 | 6 | 0.251 | 5 月 | 0.003 |  |
| [mp-147](https://next-gen.materialsproject.org/materials/mp-147) | Se6 | 6 | 2.555 | 5 月 | 1.286 |  |
| [mp-6947](https://next-gen.materialsproject.org/materials/mp-6947) | Si2O4 | 6 | 8.099 | 5 月 | 4.958 |  |
| [mp-7905](https://next-gen.materialsproject.org/materials/mp-7905) | Si2O4 | 6 | 8.237 | 5 月 | 4.554 |  |
| [mp-8352](https://next-gen.materialsproject.org/materials/mp-8352) | Si2O4 | 6 | 8.911 | 5 月 | 5.299 |  |
| [mp-1602](https://next-gen.materialsproject.org/materials/mp-1602) | Si2S4 | 6 | 4.901 | 5 月 | 2.691 |  |
| [mp-7583](https://next-gen.materialsproject.org/materials/mp-7583) | Si2S4 | 6 | 4.829 | 5 月 | 3.026 |  |
| [mp-1179443](https://next-gen.materialsproject.org/materials/mp-1179443) | Si2Se4 | 6 | 3.711 | 5 月 | 1.792 |  |
| [mp-568264](https://next-gen.materialsproject.org/materials/mp-568264) | Si2Se4 | 6 | 3.692 | 5 月 | 1.796 |  |
| [mp-550172](https://next-gen.materialsproject.org/materials/mp-550172) | Sn2O4 | 6 | 2.785 | 5 月 | 0.637 |  |
| [mp-856](https://next-gen.materialsproject.org/materials/mp-856) | Sn2O4 | 6 | 3.152 | 5 月 | 0.846 |  |
| [mp-9984](https://next-gen.materialsproject.org/materials/mp-9984) | Sn2S4 | 6 | 2.649 | 5 月 | 1.219 |  |
| [mp-10667](https://next-gen.materialsproject.org/materials/mp-10667) | Sr2Ag2P2 | 6 | 0.837 | 5 月 | 0.021 |  |
| [mp-23024](https://next-gen.materialsproject.org/materials/mp-23024) | Sr2Br2F2 | 6 | 7.638 | 5 月 | 5.040 |  |
| [mp-22957](https://next-gen.materialsproject.org/materials/mp-22957) | Sr2Cl2F2 | 6 | 8.661 | 5 月 | 5.666 |  |
| [mp-12557](https://next-gen.materialsproject.org/materials/mp-12557) | Sr2Cu2As2 | 6 | 0.528 | 5 月 | 0.381 |  |
| [mp-552537](https://next-gen.materialsproject.org/materials/mp-552537) | Sr2CuBrO2 | 6 | 4.510 | 5 月 | 2.401 |  |
| [mp-1019258](https://next-gen.materialsproject.org/materials/mp-1019258) | Sr2F4 | 6 | 9.640 | 5 月 | 5.738 |  |
| [mp-24423](https://next-gen.materialsproject.org/materials/mp-24423) | Sr2H2Br2 | 6 | 5.852 | 5 月 | 3.508 |  |
| [mp-23860](https://next-gen.materialsproject.org/materials/mp-23860) | Sr2H2Cl2 | 6 | 6.385 | 5 月 | 3.667 |  |
| [mp-24205](https://next-gen.materialsproject.org/materials/mp-24205) | Sr2H2I2 | 6 | 5.388 | 5 月 | 3.199 |  |
| [mp-1179094](https://next-gen.materialsproject.org/materials/mp-1179094) | Sr2H4 | 6 | 3.284 | 5 月 | 1.048 |  |
| [mp-23759](https://next-gen.materialsproject.org/materials/mp-23759) | Sr2H4 | 6 | 3.338 | 5 月 | 1.078 |  |
| [mp-29754](https://next-gen.materialsproject.org/materials/mp-29754) | Sr2Li2N2 | 6 | 1.168 | 5 月 | 0.248 |  |
| [mp-13276](https://next-gen.materialsproject.org/materials/mp-13276) | Sr2Li2P2 | 6 | 2.350 | 5 月 | 1.246 |  |
| [mp-1950](https://next-gen.materialsproject.org/materials/mp-1950) | Sr2S4 | 6 | 3.001 | 5 月 | 1.146 |  |
| [mp-22945](https://next-gen.materialsproject.org/materials/mp-22945) | Te2Rh2Cl2 | 6 | 1.586 | 5 月 | 0.640 |  |
| [mp-602](https://next-gen.materialsproject.org/materials/mp-602) | Te4Mo2 | 6 | 1.431 | 5 月 | 0.948 |  |
| [mp-267](https://next-gen.materialsproject.org/materials/mp-267) | Te4Ru2 | 6 | 0.685 | 5 月 | 0.209 |  |
| [mp-1019322](https://next-gen.materialsproject.org/materials/mp-1019322) | Te4W2 | 6 | 1.679 | 5 月 | 1.131 |  |
| [mp-27849](https://next-gen.materialsproject.org/materials/mp-27849) | Ti2Br2N2 | 6 | 1.791 | 5 月 | 0.536 |  |
| [mp-27850](https://next-gen.materialsproject.org/materials/mp-27850) | Ti2N2Cl2 | 6 | 1.948 | 5 月 | 0.529 |  |
| [mp-390](https://next-gen.materialsproject.org/materials/mp-390) | Ti2O4 | 6 | 3.906 | 5 月 | 2.019 |  |
| [mp-27484](https://next-gen.materialsproject.org/materials/mp-27484) | Tl4O2 | 6 | 1.157 | 5 月 | 0.630 |  |
| [mp-551470](https://next-gen.materialsproject.org/materials/mp-551470) | Tl4O2 | 6 | 1.135 | 5 月 | 0.629 |  |
| [mp-224](https://next-gen.materialsproject.org/materials/mp-224) | W2S4 | 6 | 2.219 | 5 月 | 1.287 |  |
| [mp-1821](https://next-gen.materialsproject.org/materials/mp-1821) | W2Se4 | 6 | 2.042 | 5 月 | 1.265 |  |
| [mp-28260](https://next-gen.materialsproject.org/materials/mp-28260) | WBr4O | 6 | 2.756 | 5 月 | 0.998 |  |
| [mp-27754](https://next-gen.materialsproject.org/materials/mp-27754) | WCl4O | 6 | 3.883 | 5 月 | 1.675 |  |
| [mp-10219](https://next-gen.materialsproject.org/materials/mp-10219) | Y2O2F2 | 6 | 7.585 | 5 月 | 5.032 |  |
| [mp-3637](https://next-gen.materialsproject.org/materials/mp-3637) | Y2O2F2 | 6 | 7.272 | 5 月 | 4.973 |  |
| [mp-10086](https://next-gen.materialsproject.org/materials/mp-10086) | Y2S2F2 | 6 | 2.758 | 5 月 | 1.174 |  |
| [mp-22909](https://next-gen.materialsproject.org/materials/mp-22909) | Zn2Cl4 | 6 | 6.935 | 5 月 | 3.831 |  |
| [mp-567279](https://next-gen.materialsproject.org/materials/mp-567279) | Zn2Cl4 | 6 | 6.799 | 5 月 | 3.579 |  |
| [mp-1071319](https://next-gen.materialsproject.org/materials/mp-1071319) | Zn3Te3 | 6 | 1.104 | 5 月 | 0.106 |  |
| [mp-571195](https://next-gen.materialsproject.org/materials/mp-571195) | Zn3Te3 | 6 | 2.175 | 5 月 | 0.836 |  |
| [mp-545756](https://next-gen.materialsproject.org/materials/mp-545756) | ZnSO4 | 6 | 8.766 | 5 月 | 4.571 |  |
| [mp-541912](https://next-gen.materialsproject.org/materials/mp-541912) | Zr2Br2N2 | 6 | 2.909 | 5 月 | 1.551 |  |
| [mp-570157](https://next-gen.materialsproject.org/materials/mp-570157) | Zr2Br2N2 | 6 | 3.414 | 5 月 | 1.733 |  |
| [mp-580886](https://next-gen.materialsproject.org/materials/mp-580886) | Zr2I2N2 | 6 | 2.004 | 5 月 | 0.806 |  |
| [mp-542791](https://next-gen.materialsproject.org/materials/mp-542791) | Zr2N2Cl2 | 6 | 3.251 | 5 月 | 1.820 |  |
| [mp-2574](https://next-gen.materialsproject.org/materials/mp-2574) | Zr2O4 | 6 | 6.232 | 5 月 | 3.828 |  |
| [mp-8231](https://next-gen.materialsproject.org/materials/mp-8231) | Zr2S2O2 | 6 | 1.737 | 5 月 | 0.681 |  |
| [mp-23485](https://next-gen.materialsproject.org/materials/mp-23485) | Ag2HgI4 | 7 | 2.628 | 5 月 | 1.039 |  |
| [mp-570256](https://next-gen.materialsproject.org/materials/mp-570256) | Ag2HgI4 | 7 | 2.540 | 5 月 | 0.948 |  |
| [mp-5928](https://next-gen.materialsproject.org/materials/mp-5928) | Al2CdS4 | 7 | 4.468 | 5 月 | 2.572 |  |
| [mp-3159](https://next-gen.materialsproject.org/materials/mp-3159) | Al2CdSe4 | 7 | 3.434 | 5 月 | 1.884 |  |
| [mp-7909](https://next-gen.materialsproject.org/materials/mp-7909) | Al2CdTe4 | 7 | 2.714 | 5 月 | 1.517 |  |
| [mp-7906](https://next-gen.materialsproject.org/materials/mp-7906) | Al2HgS4 | 7 | 3.486 | 5 月 | 1.902 |  |
| [mp-3038](https://next-gen.materialsproject.org/materials/mp-3038) | Al2HgSe4 | 7 | 2.671 | 5 月 | 1.326 |  |
| [mp-7910](https://next-gen.materialsproject.org/materials/mp-7910) | Al2HgTe4 | 7 | 2.050 | 5 月 | 1.034 |  |
| [mp-9254](https://next-gen.materialsproject.org/materials/mp-9254) | Al2Te5 | 7 | 1.742 | 5 月 | 0.911 |  |
| [mp-7907](https://next-gen.materialsproject.org/materials/mp-7907) | Al2ZnSe4 | 7 | 3.788 | 5 月 | 2.066 |  |
| [mp-7908](https://next-gen.materialsproject.org/materials/mp-7908) | Al2ZnTe4 | 7 | 2.722 | 5 月 | 1.557 |  |
| [mp-1591](https://next-gen.materialsproject.org/materials/mp-1591) | Al4C3 | 7 | 2.240 | 5 月 | 1.316 |  |
| [mp-9321](https://next-gen.materialsproject.org/materials/mp-9321) | Ba2HfS4 | 7 | 2.213 | 5 月 | 0.710 |  |
| [mp-3359](https://next-gen.materialsproject.org/materials/mp-3359) | Ba2SnO4 | 7 | 4.871 | 5 月 | 2.585 |  |
| [mp-3813](https://next-gen.materialsproject.org/materials/mp-3813) | Ba2ZrS4 | 7 | 1.866 | 5 月 | 0.509 |  |
| [mp-545788](https://next-gen.materialsproject.org/materials/mp-545788) | Ba3ZnN2O | 7 | 1.571 | 5 月 | 0.445 |  |
| [mp-8300](https://next-gen.materialsproject.org/materials/mp-8300) | Ba4As2O | 7 | 1.628 | 5 月 | 0.628 |  |
| [mp-9774](https://next-gen.materialsproject.org/materials/mp-9774) | Ba4Sb2O | 7 | 1.708 | 5 月 | 0.617 |  |
| [mp-550751](https://next-gen.materialsproject.org/materials/mp-550751) | BaBeSiO4 | 7 | 7.050 | 5 月 | 4.262 |  |
| [mp-1077866](https://next-gen.materialsproject.org/materials/mp-1077866) | BaZnBO3F | 7 | 6.347 | 5 月 | 3.453 |  |
| [mp-676250](https://next-gen.materialsproject.org/materials/mp-676250) | Bi2Te4Pb | 7 | 1.041 | 5 月 | 0.660 |  |
| [mp-571653](https://next-gen.materialsproject.org/materials/mp-571653) | C3N4 | 7 | 4.485 | 5 月 | 2.363 |  |
| [mp-13650](https://next-gen.materialsproject.org/materials/mp-13650) | Ca2GeO4 | 7 | 5.323 | 5 月 | 2.730 |  |
| [mp-27294](https://next-gen.materialsproject.org/materials/mp-27294) | Ca3AsBr3 | 7 | 3.546 | 5 月 | 1.490 |  |
| [mp-28069](https://next-gen.materialsproject.org/materials/mp-28069) | Ca3AsCl3 | 7 | 3.973 | 5 月 | 1.639 |  |
| [mp-30315](https://next-gen.materialsproject.org/materials/mp-30315) | Ca3BN3 | 7 | 1.412 | 5 月 | 0.434 |  |
| [mp-8789](https://next-gen.materialsproject.org/materials/mp-8789) | Ca4As2O | 7 | 2.500 | 5 月 | 1.111 |  |
| [mp-551873](https://next-gen.materialsproject.org/materials/mp-551873) | Ca4Bi2O | 7 | 1.931 | 5 月 | 0.780 |  |
| [mp-5380](https://next-gen.materialsproject.org/materials/mp-5380) | Ca4P2O | 7 | 2.504 | 5 月 | 1.143 |  |
| [mp-13660](https://next-gen.materialsproject.org/materials/mp-13660) | Ca4Sb2O | 7 | 2.014 | 5 月 | 0.865 |  |
| [mp-1025377](https://next-gen.materialsproject.org/materials/mp-1025377) | CdAg2I4 | 7 | 3.489 | 5 月 | 1.568 |  |
| [mp-4452](https://next-gen.materialsproject.org/materials/mp-4452) | CdGa2S4 | 7 | 3.568 | 5 月 | 1.982 |  |
| [mp-3772](https://next-gen.materialsproject.org/materials/mp-3772) | CdGa2Se4 | 7 | 2.575 | 5 月 | 1.225 |  |
| [mp-13949](https://next-gen.materialsproject.org/materials/mp-13949) | CdGa2Te4 | 7 | 1.999 | 5 月 | 0.946 |  |
| [mp-22304](https://next-gen.materialsproject.org/materials/mp-22304) | CdIn2Se4 | 7 | 2.020 | 5 月 | 0.801 |  |
| [mp-568032](https://next-gen.materialsproject.org/materials/mp-568032) | CdIn2Se4 | 7 | 1.967 | 5 月 | 0.719 |  |
| [mp-568661](https://next-gen.materialsproject.org/materials/mp-568661) | CdIn2Se4 | 7 | 1.795 | 5 月 | 0.605 |  |
| [mp-21374](https://next-gen.materialsproject.org/materials/mp-21374) | CdIn2Te4 | 7 | 1.729 | 5 月 | 0.715 |  |
| [mp-23752](https://next-gen.materialsproject.org/materials/mp-23752) | Cs2MgH4 | 7 | 5.345 | 5 月 | 2.672 |  |
| [mp-1078085](https://next-gen.materialsproject.org/materials/mp-1078085) | Cs2PdCl4 | 7 | 4.590 | 5 月 | 1.567 |  |
| [mp-551629](https://next-gen.materialsproject.org/materials/mp-551629) | CsLiMoO4 | 7 | 6.911 | 5 月 | 4.257 |  |
| [mp-23353](https://next-gen.materialsproject.org/materials/mp-23353) | Cu2HgI4 | 7 | 1.850 | 5 月 | 0.563 |  |
| [mp-568598](https://next-gen.materialsproject.org/materials/mp-568598) | Cu2HgI4 | 7 | 1.554 | 5 月 | 0.273 |  |
| [mp-557373](https://next-gen.materialsproject.org/materials/mp-557373) | Cu2WS4 | 7 | 2.703 | 5 月 | 1.476 |  |
| [mp-8976](https://next-gen.materialsproject.org/materials/mp-8976) | Cu2WS4 | 7 | 2.615 | 5 月 | 1.393 |  |
| [mp-4809](https://next-gen.materialsproject.org/materials/mp-4809) | Ga2HgS4 | 7 | 2.831 | 5 月 | 1.444 |  |
| [mp-4730](https://next-gen.materialsproject.org/materials/mp-4730) | Ga2HgSe4 | 7 | 1.965 | 5 月 | 0.801 |  |
| [mp-16337](https://next-gen.materialsproject.org/materials/mp-16337) | Ga2HgTe4 | 7 | 1.457 | 5 月 | 0.577 |  |
| [mp-2371](https://next-gen.materialsproject.org/materials/mp-2371) | Ga2Te5 | 7 | 1.550 | 5 月 | 0.896 |  |
| [mp-27948](https://next-gen.materialsproject.org/materials/mp-27948) | GeBi2Te4 | 7 | 1.053 | 5 月 | 0.734 |  |
| [mp-14790](https://next-gen.materialsproject.org/materials/mp-14790) | GeTe4As2 | 7 | 0.737 | 5 月 | 0.601 |  |
| [mp-3167](https://next-gen.materialsproject.org/materials/mp-3167) | Hg2GeSe4 | 7 | 1.195 | 5 月 | 0.294 |  |
| [mp-655275](https://next-gen.materialsproject.org/materials/mp-655275) | HgC2S2N2 | 7 | 4.251 | 5 月 | 1.914 |  |
| [mp-20731](https://next-gen.materialsproject.org/materials/mp-20731) | In2HgSe4 | 7 | 1.564 | 5 月 | 0.473 |  |
| [mp-19765](https://next-gen.materialsproject.org/materials/mp-19765) | In2HgTe4 | 7 | 1.277 | 5 月 | 0.412 |  |
| [mp-643257](https://next-gen.materialsproject.org/materials/mp-643257) | K2H4Pd | 7 | 5.475 | 5 月 | 2.504 |  |
| [mp-23956](https://next-gen.materialsproject.org/materials/mp-23956) | K2MgH4 | 7 | 6.201 | 5 月 | 3.156 |  |
| [mp-27138](https://next-gen.materialsproject.org/materials/mp-27138) | K2PdBr4 | 7 | 3.758 | 5 月 | 0.994 |  |
| [mp-22956](https://next-gen.materialsproject.org/materials/mp-22956) | K2PdCl4 | 7 | 4.315 | 5 月 | 1.244 |  |
| [mp-28163](https://next-gen.materialsproject.org/materials/mp-28163) | K2PdF4 | 7 | 6.603 | 5 月 | 1.805 |  |
| [mp-7502](https://next-gen.materialsproject.org/materials/mp-7502) | K2Sn2O3 | 7 | 2.471 | 5 月 | 1.027 |  |
| [mp-29584](https://next-gen.materialsproject.org/materials/mp-29584) | K4BeAs2 | 7 | 2.042 | 5 月 | 0.726 |  |
| [mp-28302](https://next-gen.materialsproject.org/materials/mp-28302) | K4CdP2 | 7 | 1.701 | 5 月 | 0.697 |  |
| [mp-29484](https://next-gen.materialsproject.org/materials/mp-29484) | K4HgAs2 | 7 | 1.546 | 5 月 | 0.591 |  |
| [mp-11719](https://next-gen.materialsproject.org/materials/mp-11719) | K4ZnP2 | 7 | 1.885 | 5 月 | 0.803 |  |
| [mp-6867](https://next-gen.materialsproject.org/materials/mp-6867) | KCaCO3F | 7 | 7.400 | 5 月 | 4.050 |  |
| [mp-865427](https://next-gen.materialsproject.org/materials/mp-865427) | KSrCO3F | 7 | 7.086 | 5 月 | 3.831 |  |
| [mp-29009](https://next-gen.materialsproject.org/materials/mp-29009) | Li2MgBr4 | 7 | 7.129 | 5 月 | 4.361 |  |
| [mp-8874](https://next-gen.materialsproject.org/materials/mp-8874) | LiSiBO4 | 7 | 9.483 | 5 月 | 6.223 |  |
| [mp-11175](https://next-gen.materialsproject.org/materials/mp-11175) | LiZnPS4 | 7 | 4.605 | 5 月 | 2.560 |  |
| [mp-7604](https://next-gen.materialsproject.org/materials/mp-7604) | Mg3NF3 | 7 | 7.026 | 5 月 | 3.849 |  |
| [mp-16755](https://next-gen.materialsproject.org/materials/mp-16755) | MgAl2S4 | 7 | 3.936 | 5 月 | 1.795 |  |
| [mp-9479](https://next-gen.materialsproject.org/materials/mp-9479) | MgAl2Se4 | 7 | 2.765 | 5 月 | 0.865 |  |
| [mp-558836](https://next-gen.materialsproject.org/materials/mp-558836) | MoF6 | 7 | 7.839 | 5 月 | 4.039 |  |
| [mp-643260](https://next-gen.materialsproject.org/materials/mp-643260) | Na2H4Pd | 7 | 4.531 | 5 月 | 2.094 |  |
| [mp-697096](https://next-gen.materialsproject.org/materials/mp-697096) | Na2H4Pt | 7 | 4.427 | 5 月 | 1.939 |  |
| [mp-28599](https://next-gen.materialsproject.org/materials/mp-28599) | Na4Br2O | 7 | 4.713 | 5 月 | 2.213 |  |
| [mp-8752](https://next-gen.materialsproject.org/materials/mp-8752) | Na4CdP2 | 7 | 1.845 | 5 月 | 0.673 |  |
| [mp-10426](https://next-gen.materialsproject.org/materials/mp-10426) | Nb2O5 | 7 | 3.090 | 5 月 | 1.657 |  |
| [mp-1595](https://next-gen.materialsproject.org/materials/mp-1595) | Nb2O5 | 7 | 3.881 | 5 月 | 1.761 |  |
| [mp-24517](https://next-gen.materialsproject.org/materials/mp-24517) | Rb2CaH4 | 7 | 6.229 | 5 月 | 2.965 |  |
| [mp-29873](https://next-gen.materialsproject.org/materials/mp-29873) | Rb2CdCl4 | 7 | 5.686 | 5 月 | 2.910 |  |
| [mp-31002](https://next-gen.materialsproject.org/materials/mp-31002) | Rb2Te5 | 7 | 1.005 | 5 月 | 0.226 |  |
| [mp-1077942](https://next-gen.materialsproject.org/materials/mp-1077942) | Rb2TeSe4 | 7 | 2.035 | 5 月 | 0.824 |  |
| [mp-863745](https://next-gen.materialsproject.org/materials/mp-863745) | RbSrCO3F | 7 | 5.332 | 5 月 | 3.851 |  |
| [mp-8560](https://next-gen.materialsproject.org/materials/mp-8560) | SF6 | 7 | 12.546 | 5 月 | 6.096 |  |
| [mp-31507](https://next-gen.materialsproject.org/materials/mp-31507) | Sb2Te4Pb | 7 | 0.678 | 5 月 | 0.398 |  |
| [mp-10955](https://next-gen.materialsproject.org/materials/mp-10955) | SnHg2Se4 | 7 | 0.997 | 5 月 | 0.114 |  |
| [mp-27947](https://next-gen.materialsproject.org/materials/mp-27947) | SnSb2Te4 | 7 | 0.519 | 5 月 | 0.307 |  |
| [mp-8299](https://next-gen.materialsproject.org/materials/mp-8299) | Sr4As2O | 7 | 2.254 | 5 月 | 1.099 |  |
| [mp-8298](https://next-gen.materialsproject.org/materials/mp-8298) | Sr4P2O | 7 | 2.316 | 5 月 | 1.107 |  |
| [mp-554867](https://next-gen.materialsproject.org/materials/mp-554867) | Ta2O5 | 7 | 1.850 | 5 月 | 0.039 |  |
| [mp-624688](https://next-gen.materialsproject.org/materials/mp-624688) | Ta2O5 | 7 | 4.468 | 5 月 | 2.316 |  |
| [mp-1025299](https://next-gen.materialsproject.org/materials/mp-1025299) | Tl2CdTe4 | 7 | 0.526 | 5 月 | 0.000 |  |
| [mp-29889](https://next-gen.materialsproject.org/materials/mp-29889) | Tl2PdCl4 | 7 | 2.985 | 5 月 | 0.778 |  |
| [mp-9791](https://next-gen.materialsproject.org/materials/mp-9791) | Tl3AsS3 | 7 | 1.792 | 5 月 | 0.959 |  |
| [mp-7684](https://next-gen.materialsproject.org/materials/mp-7684) | Tl3AsSe3 | 7 | 1.319 | 5 月 | 0.654 |  |
| [mp-8393](https://next-gen.materialsproject.org/materials/mp-8393) | Tl3SbS3 | 7 | 2.134 | 5 月 | 1.377 |  |
| [mp-720166](https://next-gen.materialsproject.org/materials/mp-720166) | V2O5 | 7 | 2.419 | 5 月 | 1.180 |  |
| [mp-28483](https://next-gen.materialsproject.org/materials/mp-28483) | WBr6 | 7 | 2.148 | 5 月 | 0.876 |  |
| [mp-546864](https://next-gen.materialsproject.org/materials/mp-546864) | Y2CN2O2 | 7 | 5.727 | 5 月 | 3.695 |  |
| [mp-1078049](https://next-gen.materialsproject.org/materials/mp-1078049) | ZnGa2S4 | 7 | 3.501 | 5 月 | 1.831 |  |
| [mp-5350](https://next-gen.materialsproject.org/materials/mp-5350) | ZnGa2S4 | 7 | 3.824 | 5 月 | 2.158 |  |
| [mp-15776](https://next-gen.materialsproject.org/materials/mp-15776) | ZnGa2Se4 | 7 | 2.746 | 5 月 | 1.313 |  |
| [mp-15777](https://next-gen.materialsproject.org/materials/mp-15777) | ZnGa2Te4 | 7 | 1.943 | 5 月 | 0.944 |  |
| [mp-22253](https://next-gen.materialsproject.org/materials/mp-22253) | ZnIn2S4 | 7 | 1.140 | 5 月 | 0.016 |  |
| [mp-541937](https://next-gen.materialsproject.org/materials/mp-541937) | ZnIn2S4 | 7 | 1.147 | 5 月 | 0.036 |  |
| [mp-22607](https://next-gen.materialsproject.org/materials/mp-22607) | ZnIn2Se4 | 7 | 2.132 | 5 月 | 0.843 |  |
| [mp-20832](https://next-gen.materialsproject.org/materials/mp-20832) | ZnIn2Te4 | 7 | 1.769 | 5 月 | 0.735 |  |
| [mp-559470](https://next-gen.materialsproject.org/materials/mp-559470) | Ag2Pb2Br2O2 | 8 | 1.724 | 5 月 | 0.082 |  |
| [mp-580941](https://next-gen.materialsproject.org/materials/mp-580941) | Ag4I4 | 8 | 3.091 | 5 月 | 1.381 |  |
| [mp-1079720](https://next-gen.materialsproject.org/materials/mp-1079720) | Ag4O4 | 8 | 0.624 | 5 月 | 0.022 |  |
| [mp-9631](https://next-gen.materialsproject.org/materials/mp-9631) | Al2Ag2O4 | 8 | 3.216 | 5 月 | 1.549 |  |
| [mp-5782](https://next-gen.materialsproject.org/materials/mp-5782) | Al2Ag2S4 | 8 | 3.494 | 5 月 | 1.799 |  |
| [mp-14091](https://next-gen.materialsproject.org/materials/mp-14091) | Al2Ag2Se4 | 8 | 2.622 | 5 月 | 1.036 |  |
| [mp-14092](https://next-gen.materialsproject.org/materials/mp-14092) | Al2Ag2Te4 | 8 | 2.215 | 5 月 | 0.983 |  |
| [mp-25469](https://next-gen.materialsproject.org/materials/mp-25469) | Al2Cl6 | 8 | 8.011 | 5 月 | 5.033 |  |
| [mp-4979](https://next-gen.materialsproject.org/materials/mp-4979) | Al2Cu2S4 | 8 | 3.090 | 5 月 | 1.767 |  |
| [mp-8016](https://next-gen.materialsproject.org/materials/mp-8016) | Al2Cu2Se4 | 8 | 2.320 | 5 月 | 0.932 |  |
| [mp-8017](https://next-gen.materialsproject.org/materials/mp-8017) | Al2Cu2Te4 | 8 | 2.242 | 5 月 | 1.050 |  |
| [mp-23933](https://next-gen.materialsproject.org/materials/mp-23933) | Al2H6 | 8 | 5.125 | 5 月 | 2.146 |  |
| [mp-9579](https://next-gen.materialsproject.org/materials/mp-9579) | Al2Tl2Se4 | 8 | 1.404 | 5 月 | 0.430 |  |
| [mp-561703](https://next-gen.materialsproject.org/materials/mp-561703) | Au2C2Cl2O2 | 8 | 5.596 | 5 月 | 2.996 |  |
| [mp-505366](https://next-gen.materialsproject.org/materials/mp-505366) | Au4Br4 | 8 | 3.096 | 5 月 | 1.457 |  |
| [mp-32780](https://next-gen.materialsproject.org/materials/mp-32780) | Au4Cl4 | 8 | 3.239 | 5 月 | 1.702 |  |
| [mp-27725](https://next-gen.materialsproject.org/materials/mp-27725) | Au4I4 | 8 | 3.074 | 5 月 | 1.609 |  |
| [mp-23225](https://next-gen.materialsproject.org/materials/mp-23225) | B2Br6 | 8 | 6.723 | 5 月 | 3.563 |  |
| [mp-23189](https://next-gen.materialsproject.org/materials/mp-23189) | B2I6 | 8 | 4.870 | 5 月 | 2.530 |  |
| [mp-1086652](https://next-gen.materialsproject.org/materials/mp-1086652) | Ba2Ag2S2F2 | 8 | 2.770 | 5 月 | 1.257 |  |
| [mp-1079647](https://next-gen.materialsproject.org/materials/mp-1079647) | Ba2Ag2Se2F2 | 8 | 2.573 | 5 月 | 1.086 |  |
| [mp-16742](https://next-gen.materialsproject.org/materials/mp-16742) | Ba2Ag2Te2F2 | 8 | 2.667 | 5 月 | 1.321 |  |
| [mp-1095151](https://next-gen.materialsproject.org/materials/mp-1095151) | Ba2Cd2As2F2 | 8 | 1.950 | 5 月 | 0.773 |  |
| [mp-1079603](https://next-gen.materialsproject.org/materials/mp-1079603) | Ba2Cd2P2F2 | 8 | 2.235 | 5 月 | 0.928 |  |
| [mp-1078693](https://next-gen.materialsproject.org/materials/mp-1078693) | Ba2Cd2Sb2F2 | 8 | 1.176 | 5 月 | 0.120 |  |
| [mp-1080650](https://next-gen.materialsproject.org/materials/mp-1080650) | Ba2Cu2S2F2 | 8 | 3.079 | 5 月 | 1.535 |  |
| [mp-9195](https://next-gen.materialsproject.org/materials/mp-9195) | Ba2Cu2Se2F2 | 8 | 2.754 | 5 月 | 1.313 |  |
| [mp-13287](https://next-gen.materialsproject.org/materials/mp-13287) | Ba2Cu2Te2F2 | 8 | 2.205 | 5 月 | 0.924 |  |
| [mp-10322](https://next-gen.materialsproject.org/materials/mp-10322) | Ba2Hf2N4 | 8 | 2.592 | 5 月 | 1.110 |  |
| [mp-18943](https://next-gen.materialsproject.org/materials/mp-18943) | Ba2Ni2O4 | 8 | 2.181 | 5 月 | 0.151 |  |
| [mp-7808](https://next-gen.materialsproject.org/materials/mp-7808) | Ba2P6 | 8 | 1.409 | 5 月 | 0.507 |  |
| [mp-4009](https://next-gen.materialsproject.org/materials/mp-4009) | Ba2Pd2S4 | 8 | 1.450 | 5 月 | 0.227 |  |
| [mp-239](https://next-gen.materialsproject.org/materials/mp-239) | Ba2S6 | 8 | 3.043 | 5 月 | 1.272 |  |
| [mp-7548](https://next-gen.materialsproject.org/materials/mp-7548) | Ba2Se6 | 8 | 2.249 | 5 月 | 0.874 |  |
| [mp-8234](https://next-gen.materialsproject.org/materials/mp-8234) | Ba2Te6 | 8 | 1.652 | 5 月 | 0.810 |  |
| [mp-1079320](https://next-gen.materialsproject.org/materials/mp-1079320) | Ba2Zn2As2F2 | 8 | 1.336 | 5 月 | 0.090 |  |
| [mp-548469](https://next-gen.materialsproject.org/materials/mp-548469) | Ba2Zn2S2O2 | 8 | 4.053 | 5 月 | 2.044 |  |
| [mp-7394](https://next-gen.materialsproject.org/materials/mp-7394) | BaAg2GeS4 | 8 | 1.730 | 5 月 | 0.365 |  |
| [mp-569790](https://next-gen.materialsproject.org/materials/mp-569790) | BaAg2GeSe4 | 8 | 1.050 | 5 月 |  |  |
| [mp-555166](https://next-gen.materialsproject.org/materials/mp-555166) | BaAg2SnS4 | 8 | 1.626 | 5 月 | 0.333 |  |
| [mp-569114](https://next-gen.materialsproject.org/materials/mp-569114) | BaAg2SnSe4 | 8 | 1.003 | 5 月 |  |  |
| [mp-613910](https://next-gen.materialsproject.org/materials/mp-613910) | BaNiF6 | 8 | 5.373 | 5 月 | 1.283 |  |
| [mp-19799](https://next-gen.materialsproject.org/materials/mp-19799) | BaPbF6 | 8 | 6.644 | 5 月 | 2.934 |  |
| [mp-5588](https://next-gen.materialsproject.org/materials/mp-5588) | BaSiF6 | 8 | 12.111 | 5 月 | 7.452 |  |
| [mp-7599](https://next-gen.materialsproject.org/materials/mp-7599) | Be4O4 | 8 | 10.489 | 5 月 | 7.231 |  |
| [mp-568301](https://next-gen.materialsproject.org/materials/mp-568301) | Bi2Br6 | 8 | 4.208 | 5 月 | 2.628 |  |
| [mp-30200](https://next-gen.materialsproject.org/materials/mp-30200) | Bi2CO5 | 8 | 2.592 | 5 月 | 0.945 |  |
| [mp-22849](https://next-gen.materialsproject.org/materials/mp-22849) | Bi2I6 | 8 | 3.490 | 5 月 | 2.360 |  |
| [mp-569157](https://next-gen.materialsproject.org/materials/mp-569157) | Bi2I6 | 8 | 3.333 | 5 月 | 2.170 |  |
| [mp-11875](https://next-gen.materialsproject.org/materials/mp-11875) | C4O4 | 8 | 12.293 | 5 月 | 5.919 |  |
| [mp-729933](https://next-gen.materialsproject.org/materials/mp-729933) | C4O4 | 8 | 11.752 | 5 月 | 5.711 |  |
| [mp-570002](https://next-gen.materialsproject.org/materials/mp-570002) | C8 | 8 | 4.587 | 5 月 | 3.064 |  |
| [mp-7801](https://next-gen.materialsproject.org/materials/mp-7801) | Ca2Ge2N4 | 8 | 4.105 | 5 月 | 2.738 |  |
| [mp-642725](https://next-gen.materialsproject.org/materials/mp-642725) | Ca2H2Cl2O2 | 8 | 8.270 | 5 月 | 5.239 |  |
| [mp-9122](https://next-gen.materialsproject.org/materials/mp-9122) | Ca2P6 | 8 | 0.437 | 5 月 | 0.076 |  |
| [mp-7204](https://next-gen.materialsproject.org/materials/mp-7204) | Ca2Zn2S2O2 | 8 | 4.416 | 5 月 | 2.487 |  |
| [mp-28553](https://next-gen.materialsproject.org/materials/mp-28553) | Ca4I2N2 | 8 | 3.545 | 5 月 | 2.261 |  |
| [mp-8255](https://next-gen.materialsproject.org/materials/mp-8255) | CaPtF6 | 8 | 6.611 | 5 月 | 2.482 |  |
| [mp-4953](https://next-gen.materialsproject.org/materials/mp-4953) | Cd2Ge2As4 | 8 | 0.380 | 5 月 | 0.132 |  |
| [mp-3668](https://next-gen.materialsproject.org/materials/mp-3668) | Cd2Ge2P4 | 8 | 1.572 | 5 月 | 0.621 |  |
| [mp-644222](https://next-gen.materialsproject.org/materials/mp-644222) | Cd2H2Cl2O2 | 8 | 4.749 | 5 月 | 2.195 |  |
| [mp-4666](https://next-gen.materialsproject.org/materials/mp-4666) | Cd2Si2P4 | 8 | 1.982 | 5 月 | 1.266 |  |
| [mp-5213](https://next-gen.materialsproject.org/materials/mp-5213) | Cd2Sn2P4 | 8 | 1.003 | 5 月 | 0.173 |  |
| [mp-1087540](https://next-gen.materialsproject.org/materials/mp-1087540) | CdCu2SiTe4 | 8 | 1.078 | 5 月 | 0.171 |  |
| [mp-5858](https://next-gen.materialsproject.org/materials/mp-5858) | CdPtF6 | 8 | 5.733 | 5 月 | 2.113 |  |
| [mp-19178](https://next-gen.materialsproject.org/materials/mp-19178) | Co2Ag2O4 | 8 | 1.837 | 5 月 | 0.415 |  |
| [mp-715566](https://next-gen.materialsproject.org/materials/mp-715566) | Cr2O6 | 8 | 4.365 | 5 月 | 1.849 |  |
| [mp-23454](https://next-gen.materialsproject.org/materials/mp-23454) | Cs2Ag2Br4 | 8 | 4.068 | 5 月 | 1.788 |  |
| [mp-10101](https://next-gen.materialsproject.org/materials/mp-10101) | Cs2Ag2C4 | 8 | 4.833 | 5 月 | 2.530 |  |
| [mp-567558](https://next-gen.materialsproject.org/materials/mp-567558) | Cs2Ag2Cl4 | 8 | 4.663 | 5 月 | 2.018 |  |
| [mp-14069](https://next-gen.materialsproject.org/materials/mp-14069) | Cs2Al2O4 | 8 | 7.060 | 5 月 | 4.324 |  |
| [mp-29506](https://next-gen.materialsproject.org/materials/mp-29506) | Cs2Bi2O4 | 8 | 4.126 | 5 月 | 2.282 |  |
| [mp-552131](https://next-gen.materialsproject.org/materials/mp-552131) | Cs2Cl2O4 | 8 | 5.559 | 5 月 | 2.122 |  |
| [mp-553310](https://next-gen.materialsproject.org/materials/mp-553310) | Cs2Cu2O4 | 8 | 3.155 | 5 月 | 0.865 |  |
| [mp-5038](https://next-gen.materialsproject.org/materials/mp-5038) | Cs2Ga2S4 | 8 | 4.977 | 5 月 | 2.882 |  |
| [mp-23057](https://next-gen.materialsproject.org/materials/mp-23057) | Cs2Li2Br4 | 8 | 7.023 | 5 月 | 4.006 |  |
| [mp-8978](https://next-gen.materialsproject.org/materials/mp-8978) | Cs2Nb2N4 | 8 | 2.498 | 5 月 | 1.415 |  |
| [mp-1079694](https://next-gen.materialsproject.org/materials/mp-1079694) | Cs4I4 | 8 | 5.409 | 5 月 | 2.839 |  |
| [mp-9820](https://next-gen.materialsproject.org/materials/mp-9820) | CsAsF6 | 8 | 10.661 | 5 月 | 5.544 |  |
| [mp-1079080](https://next-gen.materialsproject.org/materials/mp-1079080) | CsCd4As3 | 8 | 0.798 | 5 月 | 0.069 |  |
| [mp-1101093](https://next-gen.materialsproject.org/materials/mp-1101093) | CsZn4P3 | 8 | 1.103 | 5 月 | 0.290 |  |
| [mp-12954](https://next-gen.materialsproject.org/materials/mp-12954) | Cu2B2S4 | 8 | 2.682 | 5 月 | 1.659 |  |
| [mp-983565](https://next-gen.materialsproject.org/materials/mp-983565) | Cu2B2Se4 | 8 | 2.271 | 5 月 | 1.409 |  |
| [mp-23116](https://next-gen.materialsproject.org/materials/mp-23116) | Cu2Bi2Se2O2 | 8 | 0.852 | 5 月 | 0.172 |  |
| [mp-562090](https://next-gen.materialsproject.org/materials/mp-562090) | Cu2C2Cl2O2 | 8 | 5.445 | 5 月 | 2.312 |  |
| [mp-559044](https://next-gen.materialsproject.org/materials/mp-559044) | Cu2C2S2N2 | 8 | 3.556 | 5 月 | 2.035 |  |
| [mp-1078469](https://next-gen.materialsproject.org/materials/mp-1078469) | Cu2P2N4 | 8 | 3.147 | 5 月 | 1.896 |  |
| [mp-557510](https://next-gen.materialsproject.org/materials/mp-557510) | Cu3TeS3Cl | 8 | 1.649 | 5 月 | 0.731 |  |
| [mp-871](https://next-gen.materialsproject.org/materials/mp-871) | Fe4Si4 | 8 | 0.282 | 5 月 | 0.124 |  |
| [mp-5342](https://next-gen.materialsproject.org/materials/mp-5342) | Ga2Ag2S4 | 8 | 2.371 | 5 月 | 0.882 |  |
| [mp-556916](https://next-gen.materialsproject.org/materials/mp-556916) | Ga2Ag2S4 | 8 | 2.238 | 5 月 | 0.797 |  |
| [mp-4899](https://next-gen.materialsproject.org/materials/mp-4899) | Ga2Ag2Te4 | 8 | 1.238 | 5 月 | 0.135 |  |
| [mp-5238](https://next-gen.materialsproject.org/materials/mp-5238) | Ga2Cu2S4 | 8 | 1.932 | 5 月 | 0.784 |  |
| [mp-4840](https://next-gen.materialsproject.org/materials/mp-4840) | Ga2Cu2Se4 | 8 | 1.231 | 5 月 | 0.077 |  |
| [mp-3839](https://next-gen.materialsproject.org/materials/mp-3839) | Ga2Cu2Te4 | 8 | 1.219 | 5 月 | 0.223 |  |
| [mp-2507](https://next-gen.materialsproject.org/materials/mp-2507) | Ga4S4 | 8 | 3.085 | 5 月 | 1.583 |  |
| [mp-556742](https://next-gen.materialsproject.org/materials/mp-556742) | Ga4S4 | 8 | 3.181 | 5 月 | 1.582 |  |
| [mp-1572](https://next-gen.materialsproject.org/materials/mp-1572) | Ga4Se4 | 8 | 2.310 | 5 月 | 0.805 |  |
| [mp-1943](https://next-gen.materialsproject.org/materials/mp-1943) | Ga4Se4 | 8 | 2.294 | 5 月 | 0.815 |  |
| [mp-10009](https://next-gen.materialsproject.org/materials/mp-10009) | Ga4Te4 | 8 | 1.388 | 5 月 | 0.373 |  |
| [mp-1025397](https://next-gen.materialsproject.org/materials/mp-1025397) | Ge4Ru4 | 8 | 0.265 | 5 月 | 0.141 |  |
| [mp-2242](https://next-gen.materialsproject.org/materials/mp-2242) | Ge4S4 | 8 | 1.700 | 5 月 | 1.010 |  |
| [mp-700](https://next-gen.materialsproject.org/materials/mp-700) | Ge4Se4 | 8 | 1.395 | 5 月 | 0.796 |  |
| [mp-1080459](https://next-gen.materialsproject.org/materials/mp-1080459) | Ge4Te4 | 8 | 0.728 | 5 月 | 0.390 |  |
| [mp-643364](https://next-gen.materialsproject.org/materials/mp-643364) | H2Pb2Cl2O2 | 8 | 4.365 | 5 月 | 2.786 |  |
| [mp-15622](https://next-gen.materialsproject.org/materials/mp-15622) | Hf2Se6 | 8 | 1.377 | 5 月 | 0.293 |  |
| [mp-1224](https://next-gen.materialsproject.org/materials/mp-1224) | Hg4O4 | 8 | 2.487 | 5 月 | 1.103 |  |
| [mp-1546049](https://next-gen.materialsproject.org/materials/mp-1546049) | I2Cl6 | 8 | 3.831 | 5 月 | 1.475 |  |
| [mp-27729](https://next-gen.materialsproject.org/materials/mp-27729) | I2Cl6 | 8 | 4.348 | 5 月 | 1.621 |  |
| [mp-19833](https://next-gen.materialsproject.org/materials/mp-19833) | In2Ag2S4 | 8 | 1.587 | 5 月 | 0.229 |  |
| [mp-22386](https://next-gen.materialsproject.org/materials/mp-22386) | In2Ag2Te4 | 8 | 0.964 | 5 月 | 0.042 |  |
| [mp-504535](https://next-gen.materialsproject.org/materials/mp-504535) | In2H2O4 | 8 | 3.888 | 5 月 | 1.880 |  |
| [mp-632711](https://next-gen.materialsproject.org/materials/mp-632711) | In2H2O4 | 8 | 4.064 | 5 月 | 1.940 |  |
| [mp-19795](https://next-gen.materialsproject.org/materials/mp-19795) | In4S4 | 8 | 2.197 | 5 月 | 1.268 |  |
| [mp-630528](https://next-gen.materialsproject.org/materials/mp-630528) | In4S4 | 8 | 0.865 | 5 月 | 0.214 |  |
| [mp-1079260](https://next-gen.materialsproject.org/materials/mp-1079260) | In4Se4 | 8 | 1.806 | 5 月 | 0.504 |  |
| [mp-20485](https://next-gen.materialsproject.org/materials/mp-20485) | In4Se4 | 8 | 1.446 | 5 月 | 0.305 |  |
| [mp-559521](https://next-gen.materialsproject.org/materials/mp-559521) | InBi2S4Cl | 8 | 2.274 | 5 月 | 1.412 |  |
| [mp-27397](https://next-gen.materialsproject.org/materials/mp-27397) | Ir2Br6 | 8 | 3.061 | 5 月 | 1.265 |  |
| [mp-27666](https://next-gen.materialsproject.org/materials/mp-27666) | Ir2Cl6 | 8 | 3.962 | 5 月 | 1.623 |  |
| [mp-1084785](https://next-gen.materialsproject.org/materials/mp-1084785) | K2Al2O4 | 8 | 5.634 | 5 月 | 2.695 |  |
| [mp-10165](https://next-gen.materialsproject.org/materials/mp-10165) | K2Al2Te4 | 8 | 2.915 | 5 月 | 1.635 |  |
| [mp-9486](https://next-gen.materialsproject.org/materials/mp-9486) | K2AlF5 | 8 | 11.087 | 5 月 | 6.419 |  |
| [mp-30988](https://next-gen.materialsproject.org/materials/mp-30988) | K2Bi2O4 | 8 | 3.210 | 5 月 | 1.700 |  |
| [mp-1079645](https://next-gen.materialsproject.org/materials/mp-1079645) | K2CdSnSe4 | 8 | 2.906 | 5 月 | 1.439 |  |
| [mp-31368](https://next-gen.materialsproject.org/materials/mp-31368) | K2Cl2O4 | 8 | 5.956 | 5 月 | 2.097 |  |
| [mp-5256](https://next-gen.materialsproject.org/materials/mp-5256) | K2Cu2C4 | 8 | 3.999 | 5 月 | 1.942 |  |
| [mp-23846](https://next-gen.materialsproject.org/materials/mp-23846) | K2H2F4 | 8 | 11.748 | 5 月 | 6.887 |  |
| [mp-19851](https://next-gen.materialsproject.org/materials/mp-19851) | K2In2Te4 | 8 | 2.033 | 5 月 | 0.891 |  |
| [mp-28994](https://next-gen.materialsproject.org/materials/mp-28994) | K2Li4As2 | 8 | 1.813 | 5 月 | 0.644 |  |
| [mp-827](https://next-gen.materialsproject.org/materials/mp-827) | K2N6 | 8 | 7.315 | 5 月 | 4.031 |  |
| [mp-7938](https://next-gen.materialsproject.org/materials/mp-7938) | K2Nb2S4 | 8 | 1.517 | 5 月 | 0.713 |  |
| [mp-10417](https://next-gen.materialsproject.org/materials/mp-10417) | K2Sb2O4 | 8 | 3.620 | 5 月 | 1.880 |  |
| [mp-11703](https://next-gen.materialsproject.org/materials/mp-11703) | K2Sb2S4 | 8 | 2.802 | 5 月 | 1.443 |  |
| [mp-568968](https://next-gen.materialsproject.org/materials/mp-568968) | K2SnHgSe4 | 8 | 2.616 | 5 月 | 1.267 |  |
| [mp-14206](https://next-gen.materialsproject.org/materials/mp-14206) | K3Ag3As2 | 8 | 1.935 | 5 月 | 1.107 |  |
| [mp-14205](https://next-gen.materialsproject.org/materials/mp-14205) | K3Cu3As2 | 8 | 2.131 | 5 月 | 1.280 |  |
| [mp-7439](https://next-gen.materialsproject.org/materials/mp-7439) | K3Cu3P2 | 8 | 1.986 | 5 月 | 1.275 |  |
| [mp-9911](https://next-gen.materialsproject.org/materials/mp-9911) | K3SbS4 | 8 | 4.171 | 5 月 | 1.997 |  |
| [mp-7642](https://next-gen.materialsproject.org/materials/mp-7642) | K4Ag2As2 | 8 | 2.083 | 5 月 | 1.020 |  |
| [mp-1084770](https://next-gen.materialsproject.org/materials/mp-1084770) | K4Bi2Au2 | 8 | 1.705 | 5 月 | 0.842 |  |
| [mp-1079679](https://next-gen.materialsproject.org/materials/mp-1079679) | K4Cd2Sn2 | 8 | 0.456 | 5 月 | 0.003 |  |
| [mp-15684](https://next-gen.materialsproject.org/materials/mp-15684) | K4Cu2As2 | 8 | 2.064 | 5 月 | 0.958 |  |
| [mp-10381](https://next-gen.materialsproject.org/materials/mp-10381) | K4Cu2Sb2 | 8 | 2.039 | 5 月 | 1.068 |  |
| [mp-2672](https://next-gen.materialsproject.org/materials/mp-2672) | K4O4 | 8 | 5.436 | 5 月 | 2.318 |  |
| [mp-2072](https://next-gen.materialsproject.org/materials/mp-2072) | K4Te4 | 8 | 2.281 | 5 月 | 0.547 |  |
| [mp-14018](https://next-gen.materialsproject.org/materials/mp-14018) | K6As2 | 8 | 0.997 | 5 月 | 0.042 |  |
| [mp-569940](https://next-gen.materialsproject.org/materials/mp-569940) | K6Bi2 | 8 | 0.784 | 5 月 | 0.001 |  |
| [mp-7897](https://next-gen.materialsproject.org/materials/mp-7897) | K6P2 | 8 | 1.144 | 5 月 | 0.148 |  |
| [mp-14017](https://next-gen.materialsproject.org/materials/mp-14017) | K6Sb2 | 8 | 1.296 | 5 月 | 0.305 |  |
| [mp-12015](https://next-gen.materialsproject.org/materials/mp-12015) | KAg2AsO4 | 8 | 2.830 | 5 月 | 0.044 |  |
| [mp-12532](https://next-gen.materialsproject.org/materials/mp-12532) | KAg2PS4 | 8 | 2.639 | 5 月 | 1.004 |  |
| [mp-9490](https://next-gen.materialsproject.org/materials/mp-9490) | KAg2SbS4 | 8 | 1.822 | 5 月 | 0.382 |  |
| [mp-4266](https://next-gen.materialsproject.org/materials/mp-4266) | KAsF6 | 8 | 10.324 | 5 月 | 5.176 |  |
| [mp-4608](https://next-gen.materialsproject.org/materials/mp-4608) | KPF6 | 8 | 12.094 | 5 月 | 7.257 |  |
| [mp-643033](https://next-gen.materialsproject.org/materials/mp-643033) | KPH2SO3 | 8 | 7.140 | 5 月 | 4.174 |  |
| [mp-1087504](https://next-gen.materialsproject.org/materials/mp-1087504) | KZn4P3 | 8 | 1.422 | 5 月 | 0.361 |  |
| [mp-4586](https://next-gen.materialsproject.org/materials/mp-4586) | Li2Al2Te4 | 8 | 3.485 | 5 月 | 2.306 |  |
| [mp-555874](https://next-gen.materialsproject.org/materials/mp-555874) | Li2As2S4 | 8 | 1.978 | 5 月 | 1.024 |  |
| [mp-1078724](https://next-gen.materialsproject.org/materials/mp-1078724) | Li2As2Se4 | 8 | 1.409 | 5 月 | 0.754 |  |
| [mp-1205315](https://next-gen.materialsproject.org/materials/mp-1205315) | Li2Bi2O4 | 8 | 1.907 | 5 月 | 0.807 |  |
| [mp-7611](https://next-gen.materialsproject.org/materials/mp-7611) | Li2CaGeO4 | 8 | 6.693 | 5 月 | 3.982 |  |
| [mp-7610](https://next-gen.materialsproject.org/materials/mp-7610) | Li2CaSiO4 | 8 | 7.864 | 5 月 | 5.113 |  |
| [mp-5048](https://next-gen.materialsproject.org/materials/mp-5048) | Li2Ga2Te4 | 8 | 2.810 | 5 月 | 1.569 |  |
| [mp-19896](https://next-gen.materialsproject.org/materials/mp-19896) | Li2GePbS4 | 8 | 3.350 | 5 月 | 2.093 |  |
| [mp-5488](https://next-gen.materialsproject.org/materials/mp-5488) | Li2In2O4 | 8 | 3.855 | 5 月 | 1.954 |  |
| [mp-20187](https://next-gen.materialsproject.org/materials/mp-20187) | Li2In2Se4 | 8 | 2.950 | 5 月 | 1.463 |  |
| [mp-20782](https://next-gen.materialsproject.org/materials/mp-20782) | Li2In2Te4 | 8 | 2.477 | 5 月 | 1.304 |  |
| [mp-7936](https://next-gen.materialsproject.org/materials/mp-7936) | Li2Nb2S4 | 8 | 1.402 | 5 月 | 0.792 |  |
| [mp-3524](https://next-gen.materialsproject.org/materials/mp-3524) | Li2P2N4 | 8 | 5.514 | 5 月 | 3.618 |  |
| [mp-1079885](https://next-gen.materialsproject.org/materials/mp-1079885) | Li2Sb2S4 | 8 | 0.777 | 5 月 | 0.122 |  |
| [mp-8179](https://next-gen.materialsproject.org/materials/mp-8179) | Li2Tl2O4 | 8 | 1.391 | 5 月 | 0.374 |  |
| [mp-841](https://next-gen.materialsproject.org/materials/mp-841) | Li4O4 | 8 | 5.843 | 5 月 | 1.982 |  |
| [mp-757](https://next-gen.materialsproject.org/materials/mp-757) | Li6As2 | 8 | 1.719 | 5 月 | 0.748 |  |
| [mp-2341](https://next-gen.materialsproject.org/materials/mp-2341) | Li6N2 | 8 | 2.806 | 5 月 | 1.411 |  |
| [mp-736](https://next-gen.materialsproject.org/materials/mp-736) | Li6P2 | 8 | 1.822 | 5 月 | 0.775 |  |
| [mp-7955](https://next-gen.materialsproject.org/materials/mp-7955) | Li6Sb2 | 8 | 1.449 | 5 月 | 0.625 |  |
| [mp-9144](https://next-gen.materialsproject.org/materials/mp-9144) | LiAsF6 | 8 | 10.811 | 5 月 | 5.552 |  |
| [mp-1078799](https://next-gen.materialsproject.org/materials/mp-1078799) | LiNbF6 | 8 | 8.706 | 5 月 | 4.995 |  |
| [mp-1016200](https://next-gen.materialsproject.org/materials/mp-1016200) | Mg2Ge2As4 | 8 | 1.628 | 5 月 | 0.574 |  |
| [mp-34903](https://next-gen.materialsproject.org/materials/mp-34903) | Mg2Ge2P4 | 8 | 2.193 | 5 月 | 1.335 |  |
| [mp-864954](https://next-gen.materialsproject.org/materials/mp-864954) | Mg2Mo2N4 | 8 | 1.578 | 5 月 | 0.947 |  |
| [mp-1016197](https://next-gen.materialsproject.org/materials/mp-1016197) | Mg2Si2As4 | 8 | 1.831 | 5 月 | 1.131 |  |
| [mp-2961](https://next-gen.materialsproject.org/materials/mp-2961) | Mg2Si2P4 | 8 | 2.128 | 5 月 | 1.217 |  |
| [mp-562561](https://next-gen.materialsproject.org/materials/mp-562561) | Mo2O6 | 8 | 3.639 | 5 月 | 1.749 |  |
| [mp-154](https://next-gen.materialsproject.org/materials/mp-154) | N8 | 8 | 13.930 | 5 月 | 7.319 |  |
| [mp-10166](https://next-gen.materialsproject.org/materials/mp-10166) | Na2Al2Se4 | 8 | 3.822 | 5 月 | 2.002 |  |
| [mp-10163](https://next-gen.materialsproject.org/materials/mp-10163) | Na2Al2Te4 | 8 | 2.310 | 5 月 | 1.229 |  |
| [mp-22984](https://next-gen.materialsproject.org/materials/mp-22984) | Na2Bi2O4 | 8 | 2.315 | 5 月 | 1.040 |  |
| [mp-561075](https://next-gen.materialsproject.org/materials/mp-561075) | Na2CdSnS4 | 8 | 3.474 | 5 月 | 1.658 |  |
| [mp-30965](https://next-gen.materialsproject.org/materials/mp-30965) | Na2Cl2O4 | 8 | 6.403 | 5 月 | 2.197 |  |
| [mp-10164](https://next-gen.materialsproject.org/materials/mp-10164) | Na2Ga2Te4 | 8 | 1.099 | 5 月 | 0.122 |  |
| [mp-22483](https://next-gen.materialsproject.org/materials/mp-22483) | Na2In2Te4 | 8 | 1.543 | 5 月 | 0.548 |  |
| [mp-3744](https://next-gen.materialsproject.org/materials/mp-3744) | Na2Nb2O4 | 8 | 2.710 | 5 月 | 1.431 |  |
| [mp-7937](https://next-gen.materialsproject.org/materials/mp-7937) | Na2Nb2S4 | 8 | 1.371 | 5 月 | 0.601 |  |
| [mp-7939](https://next-gen.materialsproject.org/materials/mp-7939) | Na2Nb2Se4 | 8 | 1.210 | 5 月 | 0.521 |  |
| [mp-10572](https://next-gen.materialsproject.org/materials/mp-10572) | Na2P2N4 | 8 | 6.904 | 5 月 | 4.620 |  |
| [mp-5414](https://next-gen.materialsproject.org/materials/mp-5414) | Na2Sb2S4 | 8 | 1.807 | 5 月 | 0.799 |  |
| [mp-557179](https://next-gen.materialsproject.org/materials/mp-557179) | Na2Sb2S4 | 8 | 1.439 | 5 月 | 0.832 |  |
| [mp-1078419](https://next-gen.materialsproject.org/materials/mp-1078419) | Na3PSe4 | 8 | 2.726 | 5 月 | 1.008 |  |
| [mp-8703](https://next-gen.materialsproject.org/materials/mp-8703) | Na3SbSe4 | 8 | 2.538 | 5 月 | 0.919 |  |
| [mp-7392](https://next-gen.materialsproject.org/materials/mp-7392) | Na4Ag2Sb2 | 8 | 1.223 | 5 月 | 0.549 |  |
| [mp-7773](https://next-gen.materialsproject.org/materials/mp-7773) | Na4As2Au2 | 8 | 0.877 | 5 月 | 0.261 |  |
| [mp-1078709](https://next-gen.materialsproject.org/materials/mp-1078709) | Na4Bi2Au2 | 8 | 0.718 | 5 月 | 0.230 |  |
| [mp-15685](https://next-gen.materialsproject.org/materials/mp-15685) | Na4Cu2As2 | 8 | 1.342 | 5 月 | 0.588 |  |
| [mp-7639](https://next-gen.materialsproject.org/materials/mp-7639) | Na4Cu2P2 | 8 | 1.392 | 5 月 | 0.580 |  |
| [mp-2400](https://next-gen.materialsproject.org/materials/mp-2400) | Na4S4 | 8 | 3.356 | 5 月 | 1.067 |  |
| [mp-7774](https://next-gen.materialsproject.org/materials/mp-7774) | Na4Sb2Au2 | 8 | 0.774 | 5 月 | 0.210 |  |
| [mp-15668](https://next-gen.materialsproject.org/materials/mp-15668) | Na4Se4 | 8 | 2.271 | 5 月 | 0.317 |  |
| [mp-1136](https://next-gen.materialsproject.org/materials/mp-1136) | Na6As2 | 8 | 1.134 | 5 月 | 0.084 |  |
| [mp-1598](https://next-gen.materialsproject.org/materials/mp-1598) | Na6P2 | 8 | 1.395 | 5 月 | 0.388 |  |
| [mp-7956](https://next-gen.materialsproject.org/materials/mp-7956) | Na6Sb2 | 8 | 1.396 | 5 月 | 0.358 |  |
| [mp-550070](https://next-gen.materialsproject.org/materials/mp-550070) | Nb2Br4O2 | 8 | 2.379 | 5 月 | 0.700 |  |
| [mp-549720](https://next-gen.materialsproject.org/materials/mp-549720) | Nb2I4O2 | 8 | 1.691 | 5 月 | 0.597 |  |
| [mp-5621](https://next-gen.materialsproject.org/materials/mp-5621) | NbCu3S4 | 8 | 2.684 | 5 月 | 1.608 |  |
| [mp-4043](https://next-gen.materialsproject.org/materials/mp-4043) | NbCu3Se4 | 8 | 2.293 | 5 月 | 1.319 |  |
| [mp-991676](https://next-gen.materialsproject.org/materials/mp-991676) | NbCu3Te4 | 8 | 1.575 | 5 月 | 0.924 |  |
| [mp-1091384](https://next-gen.materialsproject.org/materials/mp-1091384) | NbTl3S4 | 8 | 3.446 | 5 月 | 2.291 |  |
| [mp-1025396](https://next-gen.materialsproject.org/materials/mp-1025396) | NbTl3Se4 | 8 | 2.973 | 5 月 | 1.980 |  |
| [mp-12957](https://next-gen.materialsproject.org/materials/mp-12957) | O8 | 8 | 4.349 | 5 月 | 1.322 |  |
| [mp-27529](https://next-gen.materialsproject.org/materials/mp-27529) | P2I6 | 8 | 3.804 | 5 月 | 2.001 |  |
| [mp-1954](https://next-gen.materialsproject.org/materials/mp-1954) | P3N5 | 8 | 4.263 | 5 月 | 2.477 |  |
| [mp-12883](https://next-gen.materialsproject.org/materials/mp-12883) | P8 | 8 | 5.042 | 5 月 | 3.042 |  |
| [mp-20878](https://next-gen.materialsproject.org/materials/mp-20878) | Pb4O4 | 8 | 3.168 | 5 月 | 2.108 |  |
| [mp-550714](https://next-gen.materialsproject.org/materials/mp-550714) | Pb4O4 | 8 | 3.372 | 5 月 | 2.184 |  |
| [mp-1078500](https://next-gen.materialsproject.org/materials/mp-1078500) | Pb4S4 | 8 | 1.718 | 5 月 | 1.116 |  |
| [mp-1079172](https://next-gen.materialsproject.org/materials/mp-1079172) | Pb4Se4 | 8 | 1.750 | 5 月 | 1.090 |  |
| [mp-29521](https://next-gen.materialsproject.org/materials/mp-29521) | Rb2Bi2O4 | 8 | 3.929 | 5 月 | 1.981 |  |
| [mp-549590](https://next-gen.materialsproject.org/materials/mp-549590) | Rb2Cl2O4 | 8 | 6.040 | 5 月 | 2.149 |  |
| [mp-5425](https://next-gen.materialsproject.org/materials/mp-5425) | Rb2Cu2C4 | 8 | 4.014 | 5 月 | 1.951 |  |
| [mp-22255](https://next-gen.materialsproject.org/materials/mp-22255) | Rb2In2Te4 | 8 | 2.228 | 5 月 | 1.014 |  |
| [mp-28237](https://next-gen.materialsproject.org/materials/mp-28237) | Rb2Li2Br4 | 8 | 6.942 | 5 月 | 4.075 |  |
| [mp-28243](https://next-gen.materialsproject.org/materials/mp-28243) | Rb2Li2Cl4 | 8 | 8.025 | 5 月 | 4.887 |  |
| [mp-743](https://next-gen.materialsproject.org/materials/mp-743) | Rb2N6 | 8 | 6.893 | 5 月 | 4.051 |  |
| [mp-10418](https://next-gen.materialsproject.org/materials/mp-10418) | Rb2Sb2O4 | 8 | 4.100 | 5 月 | 2.154 |  |
| [mp-7650](https://next-gen.materialsproject.org/materials/mp-7650) | Rb2Sc2O4 | 8 | 6.210 | 5 月 | 3.349 |  |
| [mp-1080691](https://next-gen.materialsproject.org/materials/mp-1080691) | Rb2SnHgTe4 | 8 | 1.567 | 5 月 | 0.531 |  |
| [mp-16764](https://next-gen.materialsproject.org/materials/mp-16764) | Rb2Y2Te4 | 8 | 2.581 | 5 月 | 1.266 |  |
| [mp-16319](https://next-gen.materialsproject.org/materials/mp-16319) | Rb6Sb2 | 8 | 0.906 | 5 月 | 0.596 |  |
| [mp-9819](https://next-gen.materialsproject.org/materials/mp-9819) | RbAsF6 | 8 | 10.532 | 5 月 | 5.423 |  |
| [mp-27421](https://next-gen.materialsproject.org/materials/mp-27421) | RbBiF6 | 8 | 6.652 | 5 月 | 2.802 |  |
| [mp-975140](https://next-gen.materialsproject.org/materials/mp-975140) | RbZn4P3 | 8 | 1.180 | 5 月 | 0.013 |  |
| [mp-27871](https://next-gen.materialsproject.org/materials/mp-27871) | Rh2Br6 | 8 | 2.611 | 5 月 | 1.016 |  |
| [mp-27770](https://next-gen.materialsproject.org/materials/mp-27770) | Rh2Cl6 | 8 | 3.388 | 5 月 | 1.333 |  |
| [mp-2068](https://next-gen.materialsproject.org/materials/mp-2068) | Rh2F6 | 8 | 3.617 | 5 月 | 0.773 |  |
| [mp-28099](https://next-gen.materialsproject.org/materials/mp-28099) | S4Br4 | 8 | 4.097 | 5 月 | 1.964 |  |
| [mp-236](https://next-gen.materialsproject.org/materials/mp-236) | S4N4 | 8 | 5.407 | 5 月 | 2.908 |  |
| [mp-23281](https://next-gen.materialsproject.org/materials/mp-23281) | Sb2I6 | 8 | 2.899 | 5 月 | 2.024 |  |
| [mp-1080747](https://next-gen.materialsproject.org/materials/mp-1080747) | Sc2Al2C2O2 | 8 | 2.631 | 5 月 | 1.105 |  |
| [mp-972937](https://next-gen.materialsproject.org/materials/mp-972937) | Sc2Al2C2O2 | 8 | 1.538 | 5 月 | 0.280 |  |
| [mp-3642](https://next-gen.materialsproject.org/materials/mp-3642) | Sc2Cu2O4 | 8 | 4.169 | 5 月 | 2.400 |  |
| [mp-13312](https://next-gen.materialsproject.org/materials/mp-13312) | Sc2Tl2S4 | 8 | 1.609 | 5 月 | 0.818 |  |
| [mp-1078181](https://next-gen.materialsproject.org/materials/mp-1078181) | Si2Br6 | 8 | 6.483 | 5 月 | 3.807 |  |
| [mp-8311](https://next-gen.materialsproject.org/materials/mp-8311) | Si3NiP4 | 8 | 0.708 | 5 月 | 0.273 |  |
| [mp-2488](https://next-gen.materialsproject.org/materials/mp-2488) | Si4Os4 | 8 | 0.821 | 5 月 | 0.497 |  |
| [mp-189](https://next-gen.materialsproject.org/materials/mp-189) | Si4Ru4 | 8 | 0.408 | 5 月 | 0.231 |  |
| [mp-8289](https://next-gen.materialsproject.org/materials/mp-8289) | Sn2F6 | 8 | 5.880 | 5 月 | 2.786 |  |
| [mp-1078644](https://next-gen.materialsproject.org/materials/mp-1078644) | Sn4O4 | 8 | 2.233 | 5 月 | 1.326 |  |
| [mp-2231](https://next-gen.materialsproject.org/materials/mp-2231) | Sn4S4 | 8 | 1.429 | 5 月 | 1.010 |  |
| [mp-1078562](https://next-gen.materialsproject.org/materials/mp-1078562) | Sr2Ag2S2F2 | 8 | 2.648 | 5 月 | 1.069 |  |
| [mp-1080438](https://next-gen.materialsproject.org/materials/mp-1080438) | Sr2Ag2Te2F2 | 8 | 2.464 | 5 月 | 1.089 |  |
| [mp-12444](https://next-gen.materialsproject.org/materials/mp-12444) | Sr2Cu2S2F2 | 8 | 3.038 | 5 月 | 1.505 |  |
| [mp-21228](https://next-gen.materialsproject.org/materials/mp-21228) | Sr2Cu2Se2F2 | 8 | 2.582 | 5 月 | 1.077 |  |
| [mp-1091416](https://next-gen.materialsproject.org/materials/mp-1091416) | Sr2Cu2Te2F2 | 8 | 2.699 | 5 月 | 1.445 |  |
| [mp-11108](https://next-gen.materialsproject.org/materials/mp-11108) | Sr2P6 | 8 | 0.348 | 5 月 | 0.021 |  |
| [mp-1175](https://next-gen.materialsproject.org/materials/mp-1175) | Sr2S6 | 8 | 2.738 | 5 月 | 0.878 |  |
| [mp-9517](https://next-gen.materialsproject.org/materials/mp-9517) | Sr2Ti2N4 | 8 | 2.045 | 5 月 | 0.695 |  |
| [mp-1080135](https://next-gen.materialsproject.org/materials/mp-1080135) | Sr2Zn2As2F2 | 8 | 1.803 | 5 月 | 0.541 |  |
| [mp-10748](https://next-gen.materialsproject.org/materials/mp-10748) | TaCu3S4 | 8 | 3.067 | 5 月 | 1.885 |  |
| [mp-4081](https://next-gen.materialsproject.org/materials/mp-4081) | TaCu3Se4 | 8 | 2.630 | 5 月 | 1.571 |  |
| [mp-9295](https://next-gen.materialsproject.org/materials/mp-9295) | TaCu3Te4 | 8 | 1.807 | 5 月 | 1.115 |  |
| [mp-7562](https://next-gen.materialsproject.org/materials/mp-7562) | TaTl3S4 | 8 | 3.538 | 5 月 | 2.376 |  |
| [mp-10644](https://next-gen.materialsproject.org/materials/mp-10644) | TaTl3Se4 | 8 | 3.066 | 5 月 | 2.104 |  |
| [mp-2552](https://next-gen.materialsproject.org/materials/mp-2552) | Te2O6 | 8 | 3.174 | 5 月 | 1.446 |  |
| [mp-569766](https://next-gen.materialsproject.org/materials/mp-569766) | Te4I4 | 8 | 1.523 | 5 月 | 0.648 |  |
| [mp-1079574](https://next-gen.materialsproject.org/materials/mp-1079574) | Te4Pb4 | 8 | 1.244 | 5 月 | 0.796 |  |
| [mp-9920](https://next-gen.materialsproject.org/materials/mp-9920) | Ti2S6 | 8 | 1.145 | 5 月 | 0.108 |  |
| [mp-27801](https://next-gen.materialsproject.org/materials/mp-27801) | Tl2Ag2I4 | 8 | 2.736 | 5 月 | 1.438 |  |
| [mp-568890](https://next-gen.materialsproject.org/materials/mp-568890) | Tl2CdGeTe4 | 8 | 0.934 | 5 月 | 0.267 |  |
| [mp-9580](https://next-gen.materialsproject.org/materials/mp-9580) | Tl2Ga2Se4 | 8 | 1.356 | 5 月 | 0.433 |  |
| [mp-3785](https://next-gen.materialsproject.org/materials/mp-3785) | Tl2Ga2Te4 | 8 | 1.085 | 5 月 | 0.452 |  |
| [mp-569246](https://next-gen.materialsproject.org/materials/mp-569246) | Tl2HgGeTe4 | 8 | 0.573 | 5 月 | 0.015 |  |
| [mp-1078621](https://next-gen.materialsproject.org/materials/mp-1078621) | Tl2In2S4 | 8 | 1.764 | 5 月 | 0.728 |  |
| [mp-20042](https://next-gen.materialsproject.org/materials/mp-20042) | Tl2In2S4 | 8 | 1.167 | 5 月 | 0.483 |  |
| [mp-22232](https://next-gen.materialsproject.org/materials/mp-22232) | Tl2In2Se4 | 8 | 1.444 | 5 月 | 0.547 |  |
| [mp-22791](https://next-gen.materialsproject.org/materials/mp-22791) | Tl2In2Te4 | 8 | 1.177 | 5 月 | 0.490 |  |
| [mp-1079426](https://next-gen.materialsproject.org/materials/mp-1079426) | Tl2SiHgSe4 | 8 | 1.746 | 5 月 | 0.828 |  |
| [mp-1078676](https://next-gen.materialsproject.org/materials/mp-1078676) | Tl2SnHgSe4 | 8 | 1.615 | 5 月 | 0.644 |  |
| [mp-5513](https://next-gen.materialsproject.org/materials/mp-5513) | Tl3VS4 | 8 | 2.946 | 5 月 | 1.895 |  |
| [mp-322](https://next-gen.materialsproject.org/materials/mp-322) | Tl4S4 | 8 | 1.430 | 5 月 | 0.461 |  |
| [mp-1836](https://next-gen.materialsproject.org/materials/mp-1836) | Tl4Se4 | 8 | 0.840 | 5 月 | 0.066 |  |
| [mp-19412](https://next-gen.materialsproject.org/materials/mp-19412) | VAg3O4 | 8 | 2.106 | 5 月 | 0.397 |  |
| [mp-510657](https://next-gen.materialsproject.org/materials/mp-510657) | VCu3O4 | 8 | 0.717 | 5 月 | 0.219 |  |
| [mp-3762](https://next-gen.materialsproject.org/materials/mp-3762) | VCu3S4 | 8 | 1.959 | 5 月 | 0.998 |  |
| [mp-21855](https://next-gen.materialsproject.org/materials/mp-21855) | VCu3Se4 | 8 | 1.657 | 5 月 | 0.782 |  |
| [mp-991652](https://next-gen.materialsproject.org/materials/mp-991652) | VCu3Te4 | 8 | 1.135 | 5 月 | 0.526 |  |
| [mp-19443](https://next-gen.materialsproject.org/materials/mp-19443) | W2O6 | 8 | 1.769 | 5 月 | 0.276 |  |
| [mp-12903](https://next-gen.materialsproject.org/materials/mp-12903) | Y2Ag2Te4 | 8 | 1.580 | 5 月 | 0.777 |  |
| [mp-546011](https://next-gen.materialsproject.org/materials/mp-546011) | Y2Zn2As2O2 | 8 | 1.962 | 5 月 | 1.107 |  |
| [mp-12509](https://next-gen.materialsproject.org/materials/mp-12509) | Y2Zn2P2O2 | 8 | 2.040 | 5 月 | 1.199 |  |
| [mp-553243](https://next-gen.materialsproject.org/materials/mp-553243) | YBi2BrO4 | 8 | 2.530 | 5 月 | 1.293 |  |
| [mp-4008](https://next-gen.materialsproject.org/materials/mp-4008) | Zn2Ge2As4 | 8 | 1.028 | 5 月 | 0.008 |  |
| [mp-4524](https://next-gen.materialsproject.org/materials/mp-4524) | Zn2Ge2P4 | 8 | 1.963 | 5 月 | 1.117 |  |
| [mp-4763](https://next-gen.materialsproject.org/materials/mp-4763) | Zn2Si2P4 | 8 | 1.981 | 5 月 | 1.199 |  |
| [mp-10281](https://next-gen.materialsproject.org/materials/mp-10281) | Zn4S4 | 8 | 3.684 | 5 月 | 1.953 |  |
| [mp-1079889](https://next-gen.materialsproject.org/materials/mp-1079889) | ZnAg2SnS4 | 8 | 1.299 | 5 月 | 0.007 |  |
| [mp-6408](https://next-gen.materialsproject.org/materials/mp-6408) | ZnCu2GeS4 | 8 | 1.471 | 5 月 | 0.403 |  |
| [mp-1078498](https://next-gen.materialsproject.org/materials/mp-1078498) | ZnCu2SiTe4 | 8 | 0.959 | 5 月 | 0.166 |  |
| [mp-1079541](https://next-gen.materialsproject.org/materials/mp-1079541) | ZnCu2SnS4 | 8 | 1.125 | 5 月 | 0.125 |  |
| [mp-1683](https://next-gen.materialsproject.org/materials/mp-1683) | Zr2Se6 | 8 | 1.519 | 5 月 | 0.440 |  |
