# BoltzTraP の入力ファイル（Si、LDA）

> 説明の本体: ecaljdoc の [manual/UsageDetailed.md](../../../ecaljdoc/manual/UsageDetailed.md)「BoltTrap」（サイト https://ecalj.github.io/ecaljdoc/manual/UsageDetailed）。サンプルの一覧は [manual/samples.md](../../../ecaljdoc/manual/samples.md)。

`lmf --boltztrap` は、輸送係数を計算するプログラム BoltzTraP の一般形式（GENE）の入力を三つ書いて止まります（表 1）。
自己無撞着計算のあと、細かい k メッシュを指定して一度だけ走らせます。ここでは Si を LDA で計算します。

## 実行

```bash
lmfa si > llmfa
mpirun -np 4 lmf si > llmf
mpirun -np 4 lmf si --boltztrap --ctrlg:bz.nkabc=[12,12,12] > llmf_boltztrap
```

三つ目は `Exit 0 ... Done boltztrap: boltztrap.* are generated` で終わります。`rst.si` と `efermi.lmf` は変わりません。

表 1: `lmf --boltztrap` が書くファイル（書いているのは `SRC/subroutines/boltztrap.f90`）

| ファイル | 中身 | BoltzTraP での名前 |
| --- | --- | --- |
| `si.intrans_template.boltztrap` | 計算の指定。3 行目にフェルミエネルギー（Ry）、エネルギーの刻み、幅、電子数 | `si.intrans` |
| `si.struct.boltztrap` | 基本並進ベクトル（bohr）、対称操作の数、対称操作（格子ベクトルを基底にした整数の 3×3 行列） | `si.struct` |
| `si.energy.isp1.boltztrap` | k 点の数、k 点ごとに「k（逆格子ベクトルを基底にした成分）とバンドの数」、続けて固有値（Ry） | `si.energy` |

スピン分極のときは `isp2` のファイルもできます。k 点はメッシュの既約な点だけです。
`si.intrans_template.boltztrap` は雛形で、温度の範囲（`Tmax`、刻み）などを書き換えて使います。

## 見るところ

表 2: 得られた値（2026-09-30、gfortran）

| 量 | 値 |
| --- | --- |
| 全エネルギー（`save.si` の `c` 行、ehf） | −15730.734181 eV |
| `si.intrans_template.boltztrap` のフェルミエネルギー | 0.2243688571 Ry（自己無撞着計算の値 = 価電子帯の頂上） |
| 同 電子数 | 8.0000 |
| `si.energy.isp1.boltztrap` の k 点 | 72 点（12×12×12 の既約点）、各点 14 バンド |
| Γ 点の固有値（1〜5 番） | −0.654858, 0.224387, 0.224387, 0.224387, 0.411057 Ry |
| 12×12×12 の点での価電子帯の頂上 / 伝導帯の底 | 0.224387 / 0.259552 Ry（差 0.478 eV） |

## 注意

- `si.struct.boltztrap` の対称操作は結晶の空間群のもの（Si は 48 個）です。
- 2026-09-29 までのコードでは、2 プロセス以上の `lmf --boltztrap` は segmentation fault で終わり、LDA の計算では対称操作が 0 個と書かれました。
- フェルミエネルギーは `rst.si` に入っている自己無撞着計算の値で、`--boltztrap` の k メッシュで決め直したものではありません。
- BoltzTraP そのものはここでは走らせていません。

## 時間

`testecalj Si -np 4` は 9〜12 秒（16 コアの計算機、4 並列）。

## テスト

`testecalj Si -np 4`（一つ上のディレクトリで）。`save.si` の全エネルギー、`si.intrans_template.boltztrap` の数値、
`si.energy.isp1.boltztrap` の全部の固有値（許容 2×10⁻⁴ Ry）、`si.struct.boltztrap` を比べます。

もとは `Samples/Legacy/BOLZTRAP`（説明一行だけで入力は無かった）です。入力は `Samples/TestInstall/si_gwsc/ctrlg.si.toml` の LDA の部分です。
