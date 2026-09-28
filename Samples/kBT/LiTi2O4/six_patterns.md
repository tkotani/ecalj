# LiTi₂O₄ MLO-QSGW: {9³, 6³} × {tf32, fp32, fp64} と空球入りの 9³ tf32 の MLO バンド（2026-09-28 から）

`update_six_rows.sh` が kt1 の run から MLO バンドを取ってきて、下の 4 枚と数値（npz）を作り直す（同じファイル名に上書き）。
経過と判断は `../kBT_research.md`（2026-09-28 の 14:38 から）。

- 枠の題の黒: その列の run（いまのコード）。橙: 代わりのもの（`[09-26 code]`・`[09-27 code]` は前のコードの仮、`[job_band]` は LDA のバンド）。
  括弧の日時はその反復が終わった時刻（`llmf.<N>run`）。灰色の「no band」はまだ無いところ
- 赤い × は Σ のメッシュ点（内挿がそこで正確）。9³ は x = 2/9, 4/9, 6/9, 8/9、6³ は 1/3, 2/3
- 数値: `liti2o4_six_rows12_data.npz`（12 列の図の各枠: x、描いたバンド、元のファイルの kt1 での場所、印、日時）、`liti2o4_six_rows_data.npz`（t2g）、
  `liti2o4_six_conv_data.npz`（反復ごとの動き）

## データの在りか（図を描き直すためのもの）

| 所 | 中身 |
| --- | --- |
| この directory（git） | 描いた数値の npz: `liti2o4_six_rows12_data.npz`（t2g と b45〜58）、`liti2o4_six_rows_data.npz`、`liti2o4_six_conv_data.npz`。各枠の x・バンド・元のファイル・日時 |
| 手元 `~/data/liti2o4_six/six_cache/<run>/` | kt1 から取ってきた MLO バンドのファイル（全バンド）と各 run の反復の終わりの時刻（`<run>.times`）。`update_six_rows.sh` の既定の置き場 |
| kt1 `/mnt/data1/LiTi2O4_kbt_runs/<run>/` | 元の MLO バンドのファイル `bnd_mlo_iter<N>.dat`（09-26 の鎖は `mloband/bnd_iter<N>.dat`）、描き直しに要る `QSGW.<N>run/`・`snap/iter<N>/` |
| kt1 `archive/mlo_bands_20260928.tar.gz` | 2026-09-28 21:44 時点のバンドのファイル 62 本の控え |

## 図 1: 14 列（6 通りと空球入りの 9³ tf32、それぞれに t2g と 2.6〜8 eV）

右の枠: eg の 8 本と 54〜58（黒の細線）、バンド 53（赤紫、4〜7 eV で分散）。計算ごとに 2 枚組で、黒い縦線が計算の境、1 つおきに背景を薄く塗った

![12 columns](liti2o4_six_rows12.png)

## 図 2: 占有の t2g の 2 本（Γ から x = 0.5）

緑が一番下、黒が 2 本目。各枠に沈み込み（x < 0.3）とガタつき（x < 0.5 の 2 階差分の平均）

![occupied pair](liti2o4_six_occ.png)

## 図 3: MLO バンドの反復ごとの動き

前の反復からの差の最大（meV、対数）。左から t2g、占有の 2 本、eg（2.5〜4.8 eV）。点線は前のコードの代わり

![change per iteration](liti2o4_six_conv.png)

## 図 4: t2g だけの 7 列（空球入りの 9³ tf32 を 9³ tf32 の隣に）

![six columns](liti2o4_six_rows.png)
