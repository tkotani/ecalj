# kbt/ — 有限温度テトラヘドロン法（method B'）作業ディレクトリ

## 命名規則
run ディレクトリ: **`n<mesh>_dq<deltaq>_T<TTTT>K[_<note>]`**
- `n666` / `n444` … k メッシュ（6³ / 4³）
- `dq0.3` / `dq0.1` … `deltaq_scale`（offset-Gamma q0 head のシフト量）
- `T<TTTT>K` … tetrakbt 温度（4桁ゼロ詰め）。`T0000K` は tetrakbt off（=T=0 厳密 lindtet6 基準）
- `_efTfix` … 有限温度フェルミ準位（heftet/leftet 離散和）を使った検証 run

## 構成
- `runs/`   本番 gwsc run（上記命名）
- `plots/`  バンド重ね描き（*.png）と gnuplot スクリプト（*.glt）
- `notes/`  実装ノート（`finite_T_tetrahedron_note.md`）
- `scripts/` run 起動スクリプト
- `debug/`  NaN調査アーカイブ（stale GPU binary 由来の偽 NaN 追跡記録・scratch cases）

## runs 一覧
| dir | mesh | dq | T | 備考 |
|---|---|---|---|---|
| n666_dq0.3_T0000K | 6³ | 0.3 | off | T=0 基準（tetrakbt off）|
| n666_dq0.3_T0500K | 6³ | 0.3 | 500 | |
| n666_dq0.3_T1000K | 6³ | 0.3 | 1000 | 実用域 |
| n666_dq0.3_T1000K_efTfix | 6³ | 0.3 | 1000 | 有限温度 E_F 検証 |
| n666_dq0.3_T2000K | 6³ | 0.3 | 2000 | |
| n666_dq0.3_T3000K | 6³ | 0.3 | 3000 | **K-Γ中央に異常スパイク（offset-Gamma artifact）** |
| n666_dq0.1_T2000K | 6³ | 0.1 | 2000 | |
| n666_dq0.1_T3000K | 6³ | 0.1 | 3000 | dq小→**異常消滅**（dq0.3 と対比）|
| n444_dq0.3_T2000K | 4³ | 0.3 | 2000 | バグ調査用 高速メッシュ |
| n444_dq0.3_T3000K | 4³ | 0.3 | 3000 | バグ調査用 高速メッシュ |

## 主要知見
- **B' は T→0 で lindtet6 に厳密収束**、≤1000K で滑らかに安定。
- **3000K の K-Γ中央バンド異常は offset-Gamma (deltaq) アーティファクト**。
  dq0.3→−1.24eV スパイク、dq0.1→消滅（`plots/overlay_3000K_dq03_vs_dq01.png`）。
  tetrakbt/B' 本体の問題ではなく q→0 head の高温×大 deltaq での破綻。

注: `plots/*.glt` 内のデータパスは旧 `plots/data/<name>` を指している可能性あり。
再描画時は `runs/n666_dq0.3_T<...>K` へ要修正。
