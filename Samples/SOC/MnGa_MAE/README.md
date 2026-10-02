# MnGa_MAE: L1₀ MnGa の磁気異方性エネルギー（力の定理）

> 説明の本体: ecaljdoc の [manual/UsageDetailed.md](../../../ecaljdoc/manual/UsageDetailed.md)「Spin-orbit coupling」（サイト https://ecalj.github.io/ecaljdoc/manual/UsageDetailed）。サンプルの一覧は [manual/samples.md](../../../ecaljdoc/manual/samples.md)。

`Samples/SOC/FePt_MAE` と同じ手順を MnGa に適用する。スピン軌道相互作用（SOC）を入れない自己無撞着計算（SCF）の
ポテンシャルの上で、SOC を入れたバンドエネルギー（`sev(eV)=`、占有状態の固有値の和）をスピン軸 001 と 110 について
一回ずつ計算し、差から磁気異方性エネルギー（MAE）を求める（Liqin Ke, Phys. Rev. B 99, 054418 (2019) の付録 A）。

$$ E_{\rm MAE} = E_{\rm band}^{\rm SOC}(110) - E_{\rm band}^{\rm SOC}(001) \tag{1}$$

正なら 001（c 軸）が磁化容易軸である。FePt に比べて MAE が一桁小さいので、k メッシュと SCF の収束がより効く。

## 入力と実行

入力の要点（`phispinsym = true`、SCF は `so = 0`、一回計算でだけ `--ctrlg:ham.so=1 --ctrlg:ham.socaxis=[..]`）、実行の手順、
試験の回し方は [../FePt_MAE/README.md](../FePt_MAE/README.md) と同じ（`fept` を `mnga` に読み替える）。MnGa に固有なのは次の 2 点:

- `[bz] nkabc = [6, 6, 6]`: 試験用の粗いメッシュで、**MAE は収束していない**（旧サンプルは 12×12×12。表 2）
- `lmxa = 6`, `pwmode = 11`, `xcfun = 103`（GGA-PBE）

## 結果の見方

`mae.py` は 4 本のログの最後の `sev(eV)=` の行を読み、`energies.txt`（表 1）と `mae.txt` を書く。

表 1: 得られた値（6×6×6、gfortran、MPI 4 並列、2026-09-30）。mmom は全磁気モーメント（μB）。

| 計算 | mmom | ehf (eV) | sev (eV) |
|---|---|---|---|
| SCF, so=0 | 2.4875 | -84429.159506 | -440.879464 |
| 一回計算 so=0（対称性なし） | 2.4875 | -84429.159469 | -440.881396 |
| 一回計算 so=1, 軸 001 | 2.4862 | -84429.179338 | -440.901265 |
| 一回計算 so=1, 軸 110 | 2.4864 | -84429.178947 | -440.900874 |

- 式 (1) から $E_{\rm MAE} = -440.900874 - (-440.901265) = 0.391$ meV。001 が容易軸。
- 表 1 の 1 行目と 2 行目の sev の差（1.9 meV）は、対称性を使うかどうかによる数値誤差の目安である。
  これは MAE より大きいが、式 (1) は同じ条件（対称性なし、同じポテンシャル）の二つの計算の差なので打ち消し合う。
- SOC によるバンドエネルギーの利得（so=1 と so=0 の sev の差）は 001 で -19.869 meV、110 で -19.478 meV。

表 2: k メッシュと MAE。

| nkabc | MAE (meV) | 出所 |
|---|---|---|
| 6×6×6 | 0.391 | このサンプル |
| 12×12×12 | 0.399 | 旧サンプルの `save.mnga` |

軌道モーメントの表（`IORBTM`）の読み方と、2026-09-30 までの版の誤りの注意は [../FePt_MAE/README.md](../FePt_MAE/README.md) と同じ。

## 計算時間

MPI 4 並列（16 コアの計算機）で、SCF が約 2.4 分、一回計算が 1 本 30〜50 秒、全体で約 4.5 分。
計算機が他のジョブで混んでいるとき（負荷が 16 コアに対して 27）は 8〜11 分かかった。

## 試験の内容（`test.py`）

- `save.mnga`: 最初の反復と収束後の ehf, ehk, mmom（`test1_check`）
- `energies.txt`: 表 1 の値（許容 2e-3）
- `mae.txt`: MAE と SOC による利得（許容 0.02 meV）

## 由来

`Samples/Legacy/SOCAXIS/MnGa`（説明は `Samples/Legacy/SOCAXIS/READMEsoaxis.md`）を今の入力形式と命令に直したもの。
