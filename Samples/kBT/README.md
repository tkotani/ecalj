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

**ここの結果はすべて 2026-06-15〜06-20 の実行で、下の 2 件の修正が入った後のもの。**

- **2026-06-13 `m_sxcf_sc.f90` の OOB ガード修正。** `ixs < 2` が正当な
  $\omega_\epsilon\approx0$ の実軸極項 (静的 $W$ のビン) を全部捨てていた。
  これ以前のビルドで走らせた run は **QP エネルギーの絶対値と $E_F$ 近傍の形が
  信用できない** (si_gwsc で 3.54 eV ずれた)。kt1 の `runs/` には修正前の run も
  残っているので、日付で判別すること (`*_oobfix` が付いているものは再計算版)。
- **2026-06-13 `tetwt5.f90` のペア選別の有限温度化** (`37e6fbc23`)。
  それ以前は上流のペア選別が sharp θ のままで、高温で χ₀ をなめらかに過小評価
  していた (1000 K で ~1%、3000 K で落ちる殻の重みの ~18%)。
  現在は `wocc = 12·kBT` で窓を広げ、`fbound`/`tolpair` で厳密上界を使って
  刈っている。

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
