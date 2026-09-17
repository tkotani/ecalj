# Materials Project から取った 10 結晶で MLO の既定を試す（DFT レベル、2026-09-18 未明）

目的: `mlo_method = 4` の既定（`mlo_delta = mlo_w = 2.0 eV`、`nskip` は k に依らない自動決定）が
サンプル以外の結晶でもそのまま動くかを、QSGW 抜き（LDA）で見る。

手順（全部 `gwinit` の既定のまま、手で触ったのは `mlo_lm` のコメントを外しただけ）:
```
POSCAR (MP の primitive) → vasp2ctrl → ctrls.<name> → ctrlgenToml.py <name>
→ mlo_lm の "! 1 X  1 2 3 4 5 6 7 8 9" 行を全原子で有効化（s,p,d の全チャネル）
→ lmfa → lmf (nkabc 8³) → getsyml → job_band → job_mlo
```
評価: MLO の各点から最近接の DFT バンドまでの差。窓は金属 [E_F−8, E_F+2]、絶縁体
[E_F−8, CBM+3] eV（`mlo_delta=2` の設計どおり CBM+2 までが目標なので +3 は少し厳しめ）。
「edge」は金属 [−2,+1]、絶縁体 [VBM−2, CBM+1]。

| 結晶 | mp-id | 型 | nskip(自動) | 窓 rms / max (meV) | edge rms / max |
|---|---|---|---|---|---|
| Ag | mp-124 | 金属 | 0 | 15 / 126 | 57 / 123 |
| Al | mp-134 | 金属（自由電子的） | 0 | 74 / 237 | 91 / 237 |
| MgO | mp-1265 | gap 4.9 | 0 | 2.1 / 15 | 1.0 / 3.6 |
| NaCl | mp-22862 | gap 5.0 | 3 (Na 2p) | 2.9 / 20 | 1.5 / 12 |
| SiC (3C) | mp-8062 | gap 1.3 | 0 | 3.3 / 18 | 2.4 / 7.7 |
| GaAs | mp-2534 | gap 0.35 | 5 (Ga 3d LO) | 14 / 130 | 5.8 / 30 |
| CdTe | mp-406 | gap 0.40 | 0 | 28 / 257 | 8.7 / 34 |
| **ZnO** (wurtzite) | mp-2133 | gap 0.75 | 0 | **476 / 1391** | 512 / 1007 |
| ZnO, d 模型に LO を使う (`--mlo_lod`) | | | 0 | **0.8 / 16** | 0.4 / 2.4 |
| TiO2 (rutile) | mp-2657 | gap 1.7 | 6 (Ti 3p LO) | 0.8 / 7.3 | 0.8 / 5.5 |
| Ni | mp-23 | 金属 | 0 | 21 / 146 | 17 / 146 |

図: `fig_<name>.png`（灰 DFT、赤 MLO、青破線 = 窓の上端）。

## 見えたこと

1. **既定で動く。** 絶縁体は数 meV、遷移金属も 15〜30 meV。Al は自由電子帯なので
   E_F+2 eV 以上は外れる（74 meV は上端付近の寄与）。
2. **ZnO だけ壊れた（476 meV）。原因は Zn 3d の扱い。** `ctrlgenToml` は Zn 3d を局所軌道
   （`pz = 3.9`）にするが、`mlo_lm` の d チャネルは常に EH 関数（4d 的）を取り、LO は
   使わない（`m_HamPMT` の選択ループは `k_table==3` を捨てる）。Ga 3d（−15 eV、`nskip=5`
   で射影子から外れる）ならそれでよいが、Zn 3d は −6 eV で O 2p と混成し窓の中にある
   （`nskip=0`）ので、EH の d では表せない。試しに d チャネルに LO を使う実験フラグ
   `--mlo_lod` を付けて回すと **0.8 meV**。GaAs に同じことをすると 14 → 20 meV と少し悪化
   （深い 3d を模型に入れ、伝導帯の 4d 的成分を失う）。
   → **半芯 LO が浅い（窓に入る）ときは LO を模型関数に、深いときは EH を使って LO 状態は
   nskip で外す**、という切り替えが要る。自動化するならエネルギー基準（LO 帯が
   E_F − 10 eV より上なら LO を使う、など）。`--mlo_lod` は実験用に残してある
   （既定 off、マニュアル未記載）。
3. **Ni は非磁性に落ちた**（`mmom = 0.6` を `[[spec]]` に置いて `nspin = 2` にしても
   0.002 μB に収束）。MLO 自体は問題ないが、`ctrlgenToml` の既定（`readp/pnufix` など）と
   磁性の相性は別途確認が要る。
4. 自動 `nskip`（k での最小値）は NaCl の Na 2p、GaAs の Ga 3d、TiO2 の Ti 3p を
   正しく数え、折れは出ていない。

作業場所: `~/trash/mp_mlo_20260918/`（入力・ログ・`eval.py`）。
