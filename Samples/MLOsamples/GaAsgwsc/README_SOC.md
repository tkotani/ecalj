# MLO with SOC-as-perturbation (temp, WIP)

## 使い方

```
job_mlo_soc gaas -np 8
```

内部で `lmf gaas --writeham --mkprocar --noinv --mlo -vnspin=2 -vso=0 --socmatrix --phispinsym`
と `mlo gaas --socmatrix` を実行。

`-vnspin=2 -vso=0` は ctrl 変数のオーバーライド。`--socmatrix` でSOC行列を別ファイル
`__HamiltonianPMTsoc` に出力。`--phispinsym` で up/dn の radial function を揃える
（SOC を摂動として扱う標準設定）。

## 出力

- `band_MLO_spin1.dat`: 2N spinor 固有値（Ry）を xdat, eigval, 1, band index で出力
- `band_MLO_spin2.dat`: 空（SOC でスピンは混合するため）

## GaAs Γ点での確認値

```
band 1-2:   -0.82928   (1s 由来, 2重縮退)
band 3-4:    0.11952   (j=1/2 split-off, 2重縮退)
band 5-8:    0.14419   (j=3/2, 4重縮退)
```

Δ_SO = 0.14419 − 0.11952 = **0.02467 Ry ≈ 0.336 eV** （実験値 0.34 eV と良一致）

## パイプライン

1. **lmf (lmf)**: nspin=2 lso=0 でスカラー相対論DFT → `__HamiltonianPMT` (H, S)
   - 同時に aughsoc で SOC 行列を計算 → `__HamiltonianPMTsoc` (hammhso[i,j,io=1..3])
   - io=1,2,3 は up-up, dn-dn, up-dn 成分（L·S 行列要素）
2. **mlo (mlo)**: PMT → MLO 基底変換
   - Hreduction で `cmlo` (PMT固有状態基底) と `zMLO = evecpmt * cmlo` (PMT基底関数基底) を取得
   - SOC を MLO に射影: `hammhso^MLO_ij = sum(dconjg(zMLO[:,i]) * hammhsop * zMLO[:,j])`
   - 実空間 `hammhsor` → 最終的に MLO バンド生成時に 2N スピノル Hamiltonian を組み立てて対角化

## 既知の限界 / WIP

1. **hammhso の symmetrization をスキップ**
   [m_HamPMT.f90](../../../SRC/subroutines/m_HamPMT.f90) で、spatial symmetry operation の一部は
   スピン軸を回すため、hammhso 3 成分 (up-up, dn-dn, up-dn) を独立にスカラー rotation
   すると symmetry 平均で SOC が相殺されてしまう。現状は identity (igg=1) のみ使用。
   正式には spin rotation matrix `D^(1/2)(R)` を含めた 2N×2N スピノル行列として
   symmetrize する必要あり。

2. **`--skiphammsoc` flag 追加済み**（[m_bandcal.f90](../../../SRC/subroutines/m_bandcal.f90)）
   lso=1 で hamm への SOC 加算をスキップする用。まだ lso=1 経路の完全動作確認は未実施
   （mrechsoc サイズや hammhsop の spinor 次元対応が必要）。

3. **テスト範囲**
   - GaAs Γ点 Δ_SO が実験値と一致することを確認
   - 他の k 点、他材料での妥当性検証はまだ
