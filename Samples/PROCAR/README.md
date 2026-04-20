# PROCAR: band weight decomposition samples

ecalj で band structure と同時に、各原子/軌道からの **weight**（fat band）を
書き出すサンプル。`lmf --mkprocar --band:fn=syml` で `PROCAR.UP[.DN]` が生成され、
後処理で atom/spin 別に重みを可視化する。

## サンプル

| ディレクトリ | 系 | 可視化対象 | リファレンス |
|------------|-----|----------|-----------|
| [`MgO_PROCAR`](MgO_PROCAR/) | MgO | O 2p weight | `bw.dat` |
| [`Ni2MnGa_L21_PROCAR`](Ni2MnGa_L21_PROCAR/) | Ni₂MnGa Heusler (磁性) | atom/spin 別 4 本 | `bandweight_atomN.spinM` |

## 実行方法 (testecalj)

```bash
cd /home/takao/ecalj/Samples/PROCAR
testecalj -np 8 MgO_PROCAR
testecalj -np 8 Ni2MnGa_L21_PROCAR
```

各テストは:
1. 入力ディレクトリを `{name}_work/` にコピー
2. `lmfa` → `lmf` でセルフコンシステント計算
3. `job_band` でバンドを生成（`--fatband` 付きで PROCAR も同時生成）
4. 後処理（gnuplot / BandWeight.py）で PDF 出力
5. リファレンスファイル（`bw.dat` や `bandweight_*.spin*`）と数値比較

テスト終了後、プロットファイル PDF のパスが表示される:
```bash
evince <workdir>/mgoWeight.pdf       # MgO
evince <workdir>/fatband.pdf          # Ni2MnGa
```

## パイプラインの仕組み

### MgO_PROCAR（単一軌道 weight）

1. `lmfa mgo`（原子計算）
2. `lmf mgo`（SC 計算）
3. `job_band mgo -np N` で PROCAR.UP.* 生成
4. `lmf --mkprocar --band:fn=syml mgo` で fat band 計算
5. `PROCAR.UP.*` を連結 → `BandWeight.py` で O(2p) 重みを `bw.dat` へ抽出
6. `gnuplot bnds.gnu.mgoW` で `mgoWeight.pdf` 作成

### Ni2MnGa_L21_PROCAR（4チャンネル fat band）

1. `lmfa` → `lmf`（SC）
2. `job_band ni2mnga --fatband --emin=-5 --emax=5` で atom/spin 別 weight を自動生成
   → `bandweight_atom1.spin1`, `bandweight_atom2.spin1`, `*.spin2`
3. `gnuplot fatband.glt` で `fatband.pdf` 作成

## 手動実行（testecalj を使わない場合）

例：MgO_PROCAR を手動で回す
```bash
cd MgO_PROCAR
lmfa mgo
mpirun -np 8 lmf mgo
job_band mgo -np 8
rm -rf PROCAR*
mpirun -np 8 lmf --mkprocar --band:fn=syml mgo
cat PROCAR.UP.* >> PROCAR.UP
./BandWeight.py > bw.dat
gnuplot bnds.gnu.mgoW
evince mgoWeight.pdf
```

他の原子/軌道の weight を見たい場合は `BandWeight.py` を編集して条件を変更する。
`fatband.glt` (Ni2MnGa 側) も plot 内容をカスタマイズ可能。

## 関連フラグ

- `--mkprocar`: PROCAR.UP[.DN] ファイルを書き出し（per-k, per-band, per-orbital weight）
- `--band:fn=syml`: symmetry line の k点で band 計算
- `--fatband --emin=X --emax=Y`: 自動で energy range [emin, emax] で fat band 計算＋atom別分離
