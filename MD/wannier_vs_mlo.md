# MLO と Wannier 関数のずれ（2026-10-02、Wannier の経路を外す前の最後の比較）

ecalj は 2026-10-02 に Wannier 関数（最大局在化、`genMLWFx`・`hmaxloc`・`hpsig_MPI`・`wanplot`・`hmagnon`）の経路を外し、
cRPA（`job_mloW --crpa`）・広がり（`mlo_spread.py`）・マグノン（`job_mlo_magnon`）を MLO だけで行うことにした（user の判断）。
その前に両者を同じ入力で比べた。比較の一式は git のタグ `last-wannier`（`dbcd6e51d`）の `Samples/WannierVsMLO`（`run.sh`、`ujk.py`、`wan_spread.py`）。
図つきのページ: https://claude.ai/artifact/EKfcmW1oGuWfcZygjVeqV1 （private）。利用者向けの要約は ecaljdoc mlo §6 の表 M7。
細かい経緯は研究ログ 2026-10-02 00:04・00:32。

## 1. 比べ方

- 同じ DFT（LDA）、同じ `[mlo] mlo_lm`（Ni d 5 本、SrVO₃ の V t₂g 3 本: lm = 5 6 8）。入力は当時の `TestInstall/ni_crpa`・`srvo3_crpa`
  に `mlo_nkabc`（Ni 8³、SrVO₃ 6³）を足したもの。Wannier の窓は Ni E_F ± 10 eV、SrVO₃ −1.9〜+3.0 eV（内窓 −1.9〜0）
- cRPA は両方とも重みの方式（PRB 83, 121101）: χ₀ の遷移 kn → k+q n′ に 1 − p_kn p_k+q,n′。
  Wannier: p_kn = Σ_m |⟨ψ_kn|w_m⟩|²（`m_wan_wfs` が `pkm4crpa` を書いた）。MLO: p_kn = [C (C†C)⁻¹ C†]_nn（`m_mlo_wfs` の `write_pkm4crpa_mlo`）
- MLO の規格化: 実空間の 2 乗積分 N_i = O_ii(R=0) で、k によらない定数で割る（ecaljdoc mlo 式 (7a)、`HamRsMLO` の末尾の記録）
- 広がり: ⟨r²⟩ = (2/N_k) Σ_k Σ_b w_b [O_ii(k) − Re M_ii(k,b)]、⟨r⟩ = −(1/N_k) Σ w_b b Im M_ii（ecaljdoc mlo 式 (7d)(7e)）。
  Wannier は `huumat --dwnb=wan` が `readgeig_mlw: ngpmx<ngp(iq)` で止まった（使われていなかった経路）ので、`UUU`（バンド間の ⟨u_k|u_k+b⟩）と
  `MLWU`（ゲージ行列 dnk）から M^W = dnk† UU dnk を組んだ。Marzari–Vanderbilt の形で `hmaxloc` の値（Ni 1.4357、1.5396）を再現

## 2. 結果

*表 1*. ω = 0、R = 0 のオンサイトの平均（eV）。U = (ii|ii)、U′ = (ii|jj)、J = (ij|ji)。Ω は bohr²（括弧は MV の形）

| 系（GW の k メッシュ） | 方法 | Ω | v の U | RPA の U | cRPA の U | cRPA の U′ | cRPA の J | (U_cRPA − U_RPA)/v |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SrVO₃ t₂g（4³） | MLO | 7.51 | 15.774 | 0.715 | 3.125 | 2.202 | 0.437 | 0.153 |
| SrVO₃ t₂g（4³） | Wannier | 6.37（6.13） | 15.992 | 0.730 | 3.149 | 2.210 | 0.444 | 0.151 |
| SrVO₃ t₂g（2³） | MLO | 5.92 | 17.750 | 1.026 | 3.512 | 2.459 | 0.495 | 0.140 |
| SrVO₃ t₂g（2³） | Wannier | 4.51（4.09） | 16.662 | 0.983 | 3.204 | 2.201 | 0.473 | 0.133 |
| Ni d（4³） | MLO | t₂g 2.02、e_g 1.90 | 23.801 | 1.405 | 2.838 | 1.543 | 0.645 | 0.060 |
| Ni d（4³） | Wannier | t₂g 1.47、e_g 1.58（1.44、1.54） | 26.002 | 1.575 | 3.779 | 2.239 | 0.768 | 0.085 |

*表 2*. 規格化の仕方と cRPA の U（eV）。W_r は p_kn だけで決まり規格化によらない。U は 1/N_i² で変わる

| 系 | なし（生） | k ごと（廃止） | GW メッシュで実空間 | 模型のメッシュで実空間（採用） | Wannier |
| --- | --- | --- | --- | --- | --- |
| Ni d（4³） | 0.407（N 0.38） | 2.875 | 2.841 | 2.838 | 3.779 |
| SrVO₃ t₂g（2³） | 0.106（N 0.18） | 3.349 | 3.124 | 3.512 | 3.204 |

*表 3*. bcc Fe のマグノン、Γ→H の Im R の極大（eV）。[`Samples/Magnon/Fe_mlo_magnon/magnon_peaks.npz`](../Samples/Magnon/Fe_mlo_magnon/magnon_peaks.npz)。MLO の値は規格化の前後で同じ

| q (2π/a) | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Wannier（`Fe_magnon`） | 0.068 | 0.133 | 0.197 | 0.197 | 0.209 | 0.319 | 0.532 | 0.740 | 0.697 | 0.637 |
| MLO、`mlo_nkabc` 8³ | 0.075 | 0.133 | 0.216 | 0.431 | 0.431 | 0.810 | 0.913 | 0.999 | 0.969 | 0.969 |
| MLO、4³ | 0.122 | 0.203 | 0.251 | 0.486 | 0.431 | 0.860 | 0.913 | 1.029 | 0.999 | 0.969 |
| MLO、4³、`mlo_w` 11 eV | 0.093 | 0.175 | 0.236 | 0.371 | 0.300 | 0.676 | 0.740 | 0.913 | 0.860 | 0.834 |

## 3. どうずれたか

1. **SrVO₃（t₂g が孤立）は合う**: 4³ で cRPA の U が 0.8 %、J が 1.6 %。p_kn がどちらもほぼ 0 か 1 で、除く遷移が同じ。
   2³ は粗すぎる（模型のメッシュで決めた N_i が GW の箱の上で 1.06 になり U が 1 割ずれる。広がりも小さく出る）
2. **Ni（d が s と絡む）は MLO の cRPA の U が 25 % 低い**:
   - RPA では同じ（U/v: MLO 0.059、Wannier 0.061）。遮蔽の強さに差は無い
   - cRPA で取り除く遮蔽が MLO は Wannier の 7 割（(U_cRPA − U_RPA)/v: 0.060 対 0.085）。差の大部分はここ
   - MLO の d が広い（Ω: t₂g +38 %、e_g +21 %）。v の U が 8 % 小さい
   - 部分空間の乗り方: 各 k で重みの大きい 5 本以外の重みが MLO 0.12、Wannier 0.28（大きい 5 本の和 4.88 対 4.72、
     5 本のうち p < 0.9 の割合 0.084 対 0.219）。Wannier の d は広い窓から作られ s・p 的なバンドにも乗るので、d ↔ s,p の遷移も一部除かれる
   - `mlo_w` 2 → 11 eV で 3.28 eV（−13 %）まで上がる（2026-10-01 23:21 の、Löwdin で規格化していた時の測定）
3. **Fe のマグノン**: q ≤ 0.3 は MLO 8³ が 20 meV 以内。q ≥ 0.4 は MLO が 1.5〜2.5 倍高い。`mlo_w` 11 eV で下がる。
   MLO の d が広く W が小さいこと（2 と同じ向き）と関係がありそうだが、確かめていない

**判断**: ずれは、d が他のバンドと絡む系で模型の部分空間をどう切り出すかの違い。どちらが正しいかは決めていない。
次の手は、実験（Fe のマグノン分散）との比較と、MLO の窓（`mlo_delta`、`mlo_w`）への依存の測定（[`MD/TODOandQuestion.md`](TODOandQuestion.md)）。
（2026-10-02 の追記: 窓への依存は Löwdin にすると消えた。§3a）

## 3a. Löwdin を標準にした後（2026-10-02 20:08 追記）

2026-10-02 から ecalj の MLO は Löwdin で直交化した関数（MLO の部分空間の射影 Wannier 関数、ecaljdoc mlo §6）。上の §2・§3 は生の MLO の値。

*表 3a*. Löwdin の MLO（窓は既定 (2, 2)、Fe は `mlo_nkabc` 8³、模型を Löwdin の H̃(R) にした版）と Wannier

| 量 | Löwdin の MLO | Wannier | 生の MLO（上の表） |
| --- | --- | --- | --- |
| Ni d、4³: cRPA の U (eV) | 2.90 | 3.78 | 2.84 |
| Ni d、4³: RPA の U (eV) | 1.43 | 1.58 | — |
| Fe のマグノン q = 0.1 / 0.3 / 0.4 / 0.6 (eV) | 0.085 / 0.229 / 0.159 / 0.339 | 0.068 / 0.197 / 0.197 / 0.319 | 0.075 / 0.216 / 0.431 / 0.810 |

- **cRPA の U**: 差は部分空間の取り方（§3 の 2）で、Löwdin でもほぼ同じ。窓を広げると cRPA の U は (4, 2) 3.13、(6, 2) 3.30 eV と動くが、RPA の U は 1.43〜1.52 eV
  とあまり動かない。cRPA は「部分空間の中の遮蔽を除く」量なので曖昧さがあり、遮蔽を全部入れた RPA の U は基底にほとんど依らない（user 2026-10-02）。
  直すものではなく、部分空間の定義の違いとして扱う
- **Fe のマグノン**: q ≥ 0.4 は Löwdin で Wannier 版にほぼ合う（生の MLO の 1.5〜2.5 倍のずれが消えた。生の MLO の窓への依存は MLO が直交していないことから来ていた）。
  q ≤ 0.3 は Löwdin が 2 割高い。q ≈ 0.4 は Stoner の連続体に入る所でピークの位置は定義しにくい。Goldstone の倍率 η は Fe 1.24（η ≠ 1 の原因は TODO）
- **残り**: 実験（Fe のマグノン分散）との比較（[`MD/TODOandQuestion.md`](TODOandQuestion.md)）

## 4. 再現

```bash
git checkout last-wannier && python3 InstallAll.py --fc gfortran --bindir <bindir>
cd Samples/WannierVsMLO && ./run.sh -np 8
```
