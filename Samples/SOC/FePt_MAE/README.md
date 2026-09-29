# FePt_MAE: L1₀ FePt の磁気異方性エネルギー（力の定理）

スピン軌道相互作用（SOC）を入れない自己無撞着計算（SCF）のポテンシャルの上で、SOC を入れたバンドエネルギーを
スピン軸 001 と 110 について一回ずつ計算し、その差から磁気異方性エネルギー（MAE）を求める
（力の定理。Liqin Ke, Phys. Rev. B 99, 054418 (2019) の付録 A）。
バンドエネルギーは占有状態の固有値の和で、`save.fept` の `sev(eV)=` に出る。MAE は

$$ E_{\rm MAE} = E_{\rm band}^{\rm SOC}(110) - E_{\rm band}^{\rm SOC}(001) \tag{1}$$

で定義する。正なら 001（c 軸）が磁化容易軸である。

## 入力の要点（`ctrlg.fept.toml`）

- `[ham] phispinsym = true`: 動径関数を上向きと下向きのスピンで共通にする。`so=1` でスピン軸を 001 以外にするときに必要で、
  SCF も同じ設定で行う（`rst.fept` と整合させるため）。
- `[ham] so = 0`: SCF は SOC なし。一回計算のときだけコマンド行で `--ctrlg:ham.so=1` と `--ctrlg:ham.socaxis=[..]` を与える。
  スピン軸は 001, 110, 100, 010 が使える。
- `[bz] nkabc = [6, 6, 4]`: 試験用の粗いメッシュで、**MAE は収束していない**（旧サンプルは 8×8×6。表 2）。
  実際の計算ではメッシュを増やして収束を確かめる（旧サンプルの説明は 16×16×12 を例に挙げている）。
- `lmxa = 6`, `pwmode = 11`, `xcfun = 103`（GGA-PBE）。

## 実行

```bash
lmfa fept > llmfa
# 1. SCF (SOC なし)
mpirun -np 4 lmf fept > llmf_scf
# 2. 一回計算。--nosym で空間群の対称性を使わず、--quit=band でバンドエネルギーまでで止める（rst.fept は書き換えない）
mpirun -np 4 lmf fept --nosym --quit=band --ctrlg:ham.so=0 > llmf_so0
mpirun -np 4 lmf fept --nosym --quit=band --ctrlg:ham.so=1 --ctrlg:ham.socaxis=[0,0,1] > llmf_so001
mpirun -np 4 lmf fept --nosym --quit=band --ctrlg:ham.so=1 --ctrlg:ham.socaxis=[1,1,0] > llmf_so110
# 3. sev を拾って式 (1) を計算する
python3 mae.py
```

試験は `Samples/SOC` で `testecalj FePt_MAE -np 4`。手で実行するときは、このディレクトリを別の場所に写してその中で行う
（`testecalj` はこのディレクトリをまるごと写して使うので、`rst.*` などをここに残さない）。

## 結果の見方

`mae.py` は 4 本のログの最後の `sev(eV)=` の行を読み、`energies.txt`（表 1）と `mae.txt` を書く。

表 1: 得られた値（6×6×4、gfortran、MPI 4 並列、2026-09-30）。mmom は全磁気モーメント（μB）。

| 計算 | mmom | ehf (eV) | sev (eV) |
|---|---|---|---|
| SCF, so=0 | 3.2573 | -535698.907637 | -367.524187 |
| 一回計算 so=0（対称性なし） | 3.2573 | -535698.907632 | -367.525380 |
| 一回計算 so=1, 軸 001 | 3.2153 | -535699.220996 | -367.838744 |
| 一回計算 so=1, 軸 110 | 3.2312 | -535699.218095 | -367.835843 |

- 式 (1) から $E_{\rm MAE} = -367.835843 - (-367.838744) = 2.901$ meV。001 が容易軸。
- 表 1 の 1 行目と 2 行目の sev の差（1.2 meV）は、対称性を使うかどうかによる数値誤差の目安である。
  ehf と ehk が一致していること（SCF で 1e-5 eV 以内）も収束の確認になる。
- SOC によるバンドエネルギーの利得（so=1 と so=0 の sev の差）は 001 で -313.364 meV、110 で -310.463 meV。
- 同じ手順でスピン軸を 100 と 010 にすると sev = -367.835554 eV（二つは対称性で等しい）で、001 との差は 3.190 meV。

表 2: k メッシュと MAE。

| nkabc | MAE (meV) | 出所 |
|---|---|---|
| 6×6×4 | 2.901 | このサンプル |
| 8×8×6 | 2.555 | 旧サンプルの `save.fept` |

`llmf_so001` などの `IORBTM: orbital moments` の表は、MT 球の中の軌道モーメントの z 成分（μB）を l ごとに示す。
スピン軸 001 では Fe 0.070、Pt 0.048。2026-03-30 から 2026-09-30 までの版は、`pwmode = 11` でこの表に誤った値を出す
（sev と mmom はどの版でも同じ）。

## 計算時間

MPI 4 並列（16 コアの計算機）で、SCF が約 2 分、一回計算が 1 本 20〜30 秒、全体で 3.6〜4.5 分（3 回の実測。他のジョブと同時）。

## 試験の内容（`test.py`）

- `save.fept`: 最初の反復と収束後の ehf, ehk, mmom（`test1_check`）
- `energies.txt`: 表 1 の値（許容 2e-3）
- `mae.txt`: MAE と SOC による利得（許容 0.02 meV）

## 由来

`Samples/Legacy/SOCAXIS/FePt`（説明は `Samples/Legacy/SOCAXIS/READMEsoaxis.md`）を今の入力形式と命令に直したもの。
