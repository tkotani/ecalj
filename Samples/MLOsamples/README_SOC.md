# MLO + SOC-as-perturbation (memo)

## 概要

MLO (MTO-based Localized Orbital) ハミルトニアンに、フル DFT/QSGW 計算とは別途に
SOC を摂動として加える方式。lmf の aughsoc を `--socmatrix` で別ファイル
(`__HamiltonianPMTsoc`) に書き出し、mlo 側で PMT→MLO 縮約と k-空間再構成を行う。

## パイプライン

`job_mlo_soc <target> -np N` で自動実行される手順（始めに前の run の `bnd*.spin*` と `band_MLO_spin2.dat` を消す）:

1. **SOC Ef 決定** — `lmf --quit=band --efermi=efermi_soc --ctrlg:ham.nspin=2 --ctrlg:ham.phispinsym=true --ctrlg:ham.so=1`
   - フル LS のスピノルハミルトニアンを全 BZ メッシュで解き、Fermi 準位を `efermi_soc` に書く
     （`--efermi=efermi_soc`。`efermi.lmf` は変えない）

1b. **SOC の DFT のバンド**（2026-10-01 から）— 同じフラグで `lmf --band --efermi=efermi_soc`
   - 対称線の上のスピン軌道ありの DFT のバンドを `bnd*.spin1`・`bandplot.isp1.glt` に書く。MLO と比べる相手。
     `qplist.dat` もこの計算が書き、MLO の窓の基準のバンド端になる
2. **H + V_SO 書き出し** — `lmf --writeham --mkprocar --noinv --mlo --ctrlg:ham.nspin=2 --ctrlg:ham.phispinsym=true --ctrlg:ham.so=0 --socmatrix`
   - スカラー相対論のハミルトニアンを `__HamiltonianPMT` に、SOC 行列を
     `__HamiltonianPMTsoc` に書く
3. **MLO バンド計算** — `mlo --ctrlg:ham.nspin=2 --ctrlg:ham.phispinsym=true --socmatrix --efermi=efermi_soc`
   - Hreduction で MLO 基底、`zMLO = evecpmt * cmlo` で PMT-basis 係数を計算
   - SOC 行列を MLO に射影 → 実空間 `hammhsor` → `HamRsMLO` に追記
   - symmetry line 上で 2N×2N スピノルを対角化 → `band_MLO_spin1.dat` (2N 固有値)
4. **プロットファイル書き換え** — `bandplot_MLO.isp1.glt` の `ef=` を `efermi_soc` の値に

## 実装上の重要点

- **m_HamPMT.f90 symmetrization**: SOC は L·S なので、little-group 回転および
  iqibz→iqbz 回転のとき、軌道部分 `rotmatt` と spin 部分 `D^(1/2)(R)` (SU(2)) を
  テンソル積した full spinor 変換 `(rotmatt ⊗ D) V (rotmatt ⊗ D)†` として作用。
  `so3_to_su2`, `spinor_rotate` を追加済み。
- **PMT→MLO projection**: `cmlo` は PMT 固有状態基底、`hammhsop` は PMT 基底関数
  基底。`zMLO = evecpmt * cmlo` で後者へ変換して sandwich:
  `hammhso^MLO = zMLO^† hammhsop zMLO`。
- **m_bandcal.f90 `--skiphammsoc`**: lso=1 で hamm に SOC を加算しないフラグ
  （摂動後段適用のため）。nsp=1 収束 rst を nsp=2 で使う際、rst サイズが
  倍化する既知の挙動。初期 rst は nsp=1 で良い。

## テストディレクトリ

| ディレクトリ | 系 | yrange (eV) | SOC Ef (Ry) | 特徴 |
|------------|-----|-------------|-------------|------|
| `GaAsSoc` | GaAs QSGW | [-2:2] | 0.1442 | 半導体, As p Δ_SO=0.336 eV |
| `FeSoc` | bcc Fe DFT | [-10:15] | 0.0201 | 磁性金属, Fe 3d |
| `FeMgOSoc` | FeMgO 多層 DFT | [-10:15] | 0.0720 | 磁性, Fe 3d + MgO |

共通ファイル構成:
- 入力: `ctrlg.<sname>.toml`, `syml.*`, `qplist.dat`, `rst.*`, `atmpnu.*.*`
- QSGW の場合は `sigm.*` も (GaAsSoc のみ)
- `bnd00N.spin1`: 試料に同梱のものは **SOC なし**の DFT（Test 1 の非 SOC の MLO の比べ相手）。`job_mlo_soc` を回すと段 1b が SOC ありで描き直す
- `bandplot.isp1.glt`: プロットテンプレート
- `band_MLO_spin{1,2}.dat`: 非SOC MLO 参照（nspin=2 なら per-spin 2 ファイル）
- `band_MLO_spin1.soc.dat`: SOC MLO 参照 (2N スピノル)
- `test.py`: job_mlo + job_mlo_soc 両方実行、許容誤差 7.4e-5 Ry (≈0.001 eV)

## 実行例

```bash
cd /home/takao/ecalj/Samples/MLOsamples
testecalj -np 8 GaAsSoc
testecalj -np 8 FeSoc
testecalj -np 8 FeMgOSoc
```

各テスト末尾に SOC プロット表示コマンドが案内される。手動で実行する場合:
```bash
cd <target>_work
gnuplot -p bandplot_MLO.isp1.glt
```
- **赤点**: MLO-SOC bands (`job_mlo_soc` 出力)
- **黒線**: `job_mlo_soc` の段 1b が描いたスピン軌道ありの DFT のバンド（`bnd00N.spin1`）

誤差は `mlo_bandcheck.py .`（SOC の MLO と SOC なしの DFT を比べようとすると WARNING を出す）。GaAsSoc は Δ_SO が DFT 0.337 / MLO 0.335 eV、
ギャップの誤差 +0.010 eV、rms 0.010 eV（ecaljdoc mlo §4 の表 M2）。

注意: 段 1b はディレクトリの `bnd*`・`qplist.dat` を SOC のものにする。同じディレクトリで後から SOC なしの `job_mlo` を回すときは、
先に `job_band` を回し直す（MLO の窓の基準のバンド端は、最後に回したバンドの計算の `qplist.dat` から取る。2026-10-01、TODO）

## 関連コミット

- `5daa001f` MLO SOC-as-perturbation pipeline (job_mlo_soc) + fixes
- `91e21e28` proper spinor symmetrization (spatial + spin rotation)
- `b9d3b2bf` Hreduction: cmlo optional, zMLO output
- `8ab9f905` job_mlo_soc: SOC-corrected Ef via lmf --quit=band
- `53b0330a`, `f6cb274d` GaAs minimal test → GaAsSoc
- `521610f8`, `f22a812d` FeMgO minimal test → FeMgOSoc
- `dec748dc` FeSoc minimal test
- `c8c9d75e`, `69ff2d49`, `708baf88`, `328acbdb` SOC DFT band overlay + yrange
