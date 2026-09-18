# kBT — 有限温度 QSGW (`tetrakbt` + `t_sigmakbt`) のサンプル

電子温度を χ₀ 側 (`tetrakbt`) と Σ 側 (`t_sigmakbt`) の**両方**に入れた QSGW の
実例。手法とその限界は
[ecaljdoc: kBT — 有限温度の自己エネルギー計算](https://ecalj.github.io/manual/kBT)。

> **これは入力と結果を置いてあるだけである。**
> GW を何十回も反復するので計算が重く、`testecalj` のターゲットにはしていない。
> 手法が何を変えるかを、収束した結果そのもので確認するためのもの。

ブランチ: `t_sigmakbt` (Σ 側。χ₀ 側の `tetrakbt` は `main` にある)

| | 中身 | 大きさ |
|---|---|---|
| [`LiTi2O4/`](LiTi2O4/) | 金属スピネル。$6^3$/$9^3$ のメッシュ収束、反復収束、温度依存、`deltaq` の罠 | 66 MB |
| [`Fe/`](Fe/) | **`t_sigmakbt` だけを切り替えた対照実験。** Σ 側が効くことの最短の証拠 | 0.8 MB |

---

## 1. まずここだけ読む — 使うときの注意

### `t_sigmakbt` は `tetrakbt = true` と必ず併用する

Σ 側が読む `EFERMI_kbt` を書くのは `tetrakbt` を有効にした `heftet` である。
単独で指定すると警告が一行出て、**有限温度が適用されないまま静かに走る**。

```
sigmakbt_setup: WARNING t_sigmakbt>0 but EFERMI_kbt missing (need tetrakbt/heftet).
                Sigma-side finite-T NOT applied.
```

### `t_sigmakbt == t_tetrakbt` にする

そうして初めて χ₀ ($W$) と Σ ($G$) が一つの温度を共有する。
片方だけ温めると何を計算しているのか分からなくなる。Fe のサンプルは
**それを承知で片方だけ切った対照実験**であって、実用の設定ではない。

### `esmr` は既定 (0.01 Ry) のままにする

`t_sigmakbt > 0` にしても Σ 側の**状態範囲**は依然 `esmr` で決まっている
(窓 = ±10·esmr)。ここのサンプルの設定では窓が ±7.9 $k_BT$ なので切り捨ては
3.7e-4 で無害だが、`esmr` を下げたり $T$ を上げたりすると**警告なしに**
占有数が切り捨てられる (5000 K で 4%、esmr=0.002 で 17%)。
→ [doc §7.3](https://ecalj.github.io/manual/kBT)

### 高温では `deltaq_scale` も小さくする

3000 K・`deltaq_scale = 0.3` では K–Γ 中央に −1.24 eV の偽のスパイクが出る。
`deltaq_scale = 0.1` で消える。offset-Gamma の $q\to0$ head が高温 × 大きい
`deltaq` で破綻するもので、`tetrakbt` のバグではない
([`LiTi2O4/README.md` §3.5](LiTi2O4/README.md))。

### $T$ について収束を確認する

`t_tetrakbt` は物理的な温度というより **QSGW の反復を安定化する正則化**である。
求めたい量が $T$ について収束していること、できれば同じ $T$ で 2 つの k メッシュが
一致することを確認すること (それが `LiTi2O4/` の $6^3$ vs $9^3$)。

---

## 2. どのバイナリで走らせたか (結果の信頼性)

**ここの結果はすべて 2026-06-15〜06-20 の実行である。**

### 有限温度に固有 — `tetrakbt` のペア選別 (コミット済みの窓が実在する)

**`6b83b86e7` (2026-06-09、method B′ の導入) から `37e6fbc23` (2026-06-13) までの
コミットで `tetrakbt = true` を使った run は、高温で χ₀ をなめらかに過小評価
している。** `tetwt5` の上流のペア選別が sharp θ のままだったためで、落ちる殻の
重みの 1000 K で ~1%、3000 K で ~18%。現在は `wocc = 12*kbt` で窓を広げ、
`fbound`/`tolpair` の厳密上界で刈っている。`usetetrakbt` の分岐の中の話なので、
T=0 の経路 (sharp θ が正しい) には影響しない。

### 有限温度とは無関係 — 実軸極項の OOB ガード (コミットはされていない)

`m_sxcf_sc.f90` の実軸極項のビン詰めには 2026-06-13 の `2298e75e9` まで
**範囲外ガードが無かった**。`findloc` が 0 を返す (`we` が `freq_r(nw)` を超える) と
`nttp(-1)` / `wgtiw(:,-1)` へ書き込む。温度によらず全ての gwsc が通る経路である。

この修正を作る途中、`ixs < 2` という**下限を 1 つ間違えた**版が一時的に存在した。
配列は `nttp(0:nw)` / `wgtiw(:,0:nw)` / `freq_r(0:nw)` (`freq_r(0)=0`) と 0 始まりなので
正しい下限は `ixs >= 1` であり、`ixs >= 2` にすると補間三つ組の下端 bin 0 =
**静的 W のスロット**が除外される。そこが on-shell (ω_ε≈0) の極項 ≈ −Wc(0)/2 の
置き場所なので、それが丸ごと落ちて si_gwsc の QPU が 3.54 eV ずれた。
**ただしこの版はコミットされていない** — リポジトリの歴史は「ガード無し →
正しいガード」で、`ixs < 2` は kt1 の作業ツリーにだけ 2026-06-11〜06-13 存在した。
したがってこれは **kt1 の `runs/` にある当時の run についての注意**であって
(`*_oobfix` が付いているものは再計算版)、ecalj の利用者には関係しない。

ここに置いた run はどちらの窓にも入っていない。

---

## 3. kt1 に残っているもの

計算はすべて kt1 の `~/LiTi2O4/kbt/` と `~/sigmakbt_test/` で回した。
ここには入力と結果の抜粋しか置いていない。

run ディレクトリの命名規則: **`n<mesh>_dq<deltaq>_T<TTTT>K[_<note>]`**

- `n666` / `n999` / `n444` … k メッシュ
- `dq0.1` / `dq0.3` … `deltaq_scale` (offset-Gamma の $q_0$ head のシフト量)
- `T<TTTT>K` … `t_tetrakbt` (4 桁ゼロ詰め)。`T0000K` は `tetrakbt` off
- `_sigmakbtNNNN` … Σ 側も同じ温度にしたもの (ここのサンプルはこれ)
- `_scf10`, `_baseline`, `_oobfix`, `_efTfix` … 検証用の枝

kt1 側にあってここに持ってきていないもの:

| kt1 の場所 | 中身 |
|---|---|
| `LiTi2O4/kbt/runs/n*_T1000K_sigmakbt1000/band/iter*/` | 1000 K の全反復のバンド生データ (6³ 17 MB, 9³ 26 MB)。ここには最終反復だけと、全反復をまとめた PDF を置いてある |
| `LiTi2O4/kbt/notes/finite_T_tetrahedron_note.md` | method B′ の式番号付き実装ノート |
| `LiTi2O4/kbt/notes/offset_gamma_deltaq_anomaly.md` | `deltaq` 異常の調査 (3 層の根拠、結論確定) |
| `LiTi2O4/kbt/runs/REPORT_20260613.md` | 上の OOB バグ発見の記録 |
| `LiTi2O4/kbt/runs/GK_*`, `orbital_GK_*` | Γ–K の軌道分解 (本件とは別テーマ) |
| `LiTi2O4/kbt/runs/bisect_*`, `diagA_*` | バグ追跡用の枝 |

`LiTi2O4/kbt/plots/*.glt` のデータパスは古い場所を指していることがある。

### 9³ の 3000 K は無い (探さなくてよい)

9³ で反復を追える run は **1000 K の `n999_dq0.1_T1000K_sigmakbt1000` (45 反復) だけ**である。
3000 K の 9³ は 3 つとも途中で終わっている:

| kt1 の run | 状態 |
|---|---|
| `n999_dq0.1_T3000K_sigmakbt3000` | **一度も走っていない。** 7 ファイル 44 KB で、`chain.log` に「2026-06-16 17:27:20 chain(9^3) started; waiting for n666_dq0.1_T1000K_sigmakbt1000/done.txt」とあるだけ。連鎖が発火しなかった |
| `n999_dq0.1_T3000K_baseline` | 反復 1 の途中 (`hgw --jobgw=1`) で止まっている。QPU もバンドも無い |
| `n999_dq0.1_T3000K_scf10` | 反復 2 まで。バンドは `band_iter1.png` の 1 枚だけ |
| `bisect_999_{dq015,fp64,sm05}` | 3000 K の 1 反復。バグ追跡用の診断であって結果ではない |

したがって **3000 K のメッシュ比較 (6³ vs 9³) はできない。** メッシュ依存の話は
1000 K の対 (§3.3)、Σ 側の有無の話は 3000 K の 6³ の対 (§3.4) と、別々に見ること。

---

## 4. これから詰めるべき課題

研究ログ（日付順、試したこと・数字・仮説）は [kBT_research.md](kBT_research.md)。

一覧は [ecaljdoc: kBT §9](https://ecalj.github.io/ecaljdoc/manual/kBT#9-これから詰めるべき課題-2026-09-17)。
要点だけ:

- `t_sigmakbt` は **300 K でも金属で $\Sigma - v_{xc}$ が 0.1 eV 動く** (Fe $5^3$, 1 反復:
  χ₀ だけなら ≤1 meV、Σ も温めると max 0.12 eV)。$E_F$ 近傍の少数の状態の分数占有が
  原因で、メッシュ収束と `esmr` との関係が未整理。gwinit のテンプレは 2026-09-17 に
  有限温度キーを見え消しに戻した (新規 ctrlg は T=0)。
- 絶縁体の `EFERMI_kbt` は 2026-09-17 まで書かれていなかった (heftet 修正済み)。
  ギャップ内の $E_F$ の決め方は暫定。
- `testecalj` に有限温度の回帰ターゲットが無い。
- $9^3$ 3000 K (kt1) は 2 反復で停止中。再開前に `ctrlg_absorb.py liti2o4`。

