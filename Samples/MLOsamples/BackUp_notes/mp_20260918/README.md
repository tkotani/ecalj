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
| ZnO, 浅い LO の自動判定後（既定） | | | 0 | **0.8 / 16** | 0.4 / 2.4 |
| TiO2 (rutile) | mp-2657 | gap 1.7 | 6 (Ti 3p LO) | 0.8 / 7.3 | 0.8 / 5.5 |
| Ni（非磁性に落ちた run） | mp-23 | 金属 | 0 | 21 / 146 | 17 / 146 |
| Ni, `mix="A3"` で 0.66 μB | mp-23 | 強磁性 | 0 | 20 / 148 (↑), 20 / 134 (↓) | 19 / 136, 16 / 134 |
| Fe, `mix="A3"` で 2.24 μB | mp-13 | 強磁性 | 3 (Fe 3p LO) | 17 / 218 (↑), 18 / 147 (↓) | 19 / 218, 20 / 147 |

図: `fig_<name>.png`（灰 DFT、赤 MLO、青破線 = 窓の上端）。磁性は `fig_fe_mag.png`, `fig_ni_mag.png`。

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
   → **自動化した（同日夜、commit 参照）。** `m_HamPMT` が模型選択の前に PMT ハミルトニアンを
   全 k で対角化し、半芯 LO（0 < pz < 10 かつ int(pz) < int(pnu)）ごとに「LO 部分空間への
   射影重み > 1/2 の占有状態」の最高エネルギーを取る。それが **E_F − 10 eV より上なら
   浅い LO として LO を模型関数に**（EH は捨てる）、下なら従来どおり EH を使い LO 状態は
   nskip で外す。ZnO の Zn 3d は −3.8 eV → LO を使って 0.8 meV、GaAs Ga 3d −14.7 eV、
   NaCl Na 2p −20.7、TiO2 Ti 3p −32.8、Fe 3p −50.9 は深い → 変化なし。
   拡張 LO（pz が価電子殻より上、例 RuO2 の Ru pz=5.5）は候補にしない（試すと 14 → 254 meV）。
   18 サンプルで参照が動いたのは NiO だけ（Ni 3d LO が浅い判定、16.6 → 16.4 meV、
   最大 52 meV の差、参照更新）。Gd/Sm の 4f は pz と pnu が同じ殻なので候補外だが、
   試しに LO を使うと GdCo5 12.6 → 7.4、SmP 76 → 40 meV と良くなる — 同殻 LO の扱いは
   要検討。`--mlo_lod` は廃止（見え消し）。
3. **Ni・Fe が非磁性に落ちた — 原因は `ctrlgenToml` 既定の混合 `mix = "B3", b = 0.2`。**
   Fe (mp-13) で 1 反復目は 2.13 μB なのに 2 反復目で 0.02 μB に潰れる。Anderson `A3`
   （b=0.3）または `B3` でも `b=0.5` なら 2.24 μB を保つ（TestInstall/fe は `A6, b=0.5`）。
   `ctrlgenToml.py --nspin=2` は `mix = "A3", b = 0.3` を書くようにし、非磁性の既定
   `B3` の行と [iter] の説明に「磁性なら A3」を明記した（commit 参照）。修正後の
   Fe / Ni の MLO は上の表のとおり両スピン 17〜20 meV。
4. 自動 `nskip`（k での最小値）は NaCl の Na 2p、GaAs の Ga 3d、TiO2 の Ti 3p を
   正しく数え、折れは出ていない。

作業場所: `~/trash/mp_mlo_20260918/`（入力・ログ・`eval.py`）。
