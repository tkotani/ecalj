# MLO-QSGW のサンプル

`gwsc --mlo` で回す QSGW。自己エネルギー $\Sigma$ を MLO 表現で持ち、k 空間で内挿する（MLO Sigma インターポレーション）。
理論と設定の説明は ecaljdoc の [MLO-gwsc](https://ecalj.github.io/ecaljdoc/manual/mlo_gwsc) にある。

| ディレクトリ | 系 | k 点 | MLO | 所要時間（-np 4） |
|---|---|---|---|---|
| `GaAs/` | 閃亜鉛鉱 GaAs、非磁性 | 2×2×2 | Ga, As の s+p+d | 約 1 分 |
| `NiO/`  | NiO、AF (Ni↑, Ni↓, O) | 2×2×2 | Ni の s+p+d, O の s+p | 約 3 分 |

どちらも TestInstall の `gas_gwsc` / `nio_gwsc` の入力に `[mlo]` を足しただけ。k 点は粗く、物理量の精度は見ていない。コードが通ることの確認用。

## 実行

```bash
cd ecalj/Samples/MLOQSGW
testecalj -np 4 GaAs NiO
```

`GaAs_work/`、`NiO_work/` に計算結果ができ、`QPU`（NiO は `QPD` も）と `log.<sname>` の `fp evl` を参照と比べる。
中身は次と同じ。

```bash
gwsc 2 -np 4 --mlo gaas
```

LDA から 2 反復。2 反復目の lmf は、1 反復目の $\Sigma^{\rm MLO}$ を内挿した $\Sigma$ で解く。

## 入力で要るもの

`ctrlg.<sname>.toml` の末尾の `[mlo]`:

- `mlo_nkabc` は `[bz] nkabc` と同じにする（`gwsc` が開始前に調べて、違えば止まる）。
- `mlo_lm` に MLO にする軌道を書く。
- `[gw] mixbeta = 0.5`。$\Sigma^{\rm MLO}$ も `sigm` と同じく混合する（2026-09-30 から既定。それより前のコードでは `ECALJ_MLO_MIX=1` が要った）。

## 残るファイル

| ファイル | 中身 |
|---|---|
| `QMLO_SigRs` | 実空間の $\Sigma^{\rm MLO}$。次の反復・再開はこれを読む |
| `QMLO_z` | $\Sigma^{\rm MLO}$ を書いたときの $z^{\rm MLO}$（全 k 点・スピンで 1 ファイル） |
| `HamRsMLO` | MLO の索引と窓 |
| `sigm` | QSGW の $\Sigma$（PMT 基底）。これがあると再開になる |

`__QMLO_*` は作業ファイルで、消してよい（`__QMLO_mixsig` は混合の履歴。消すと次の反復は線形混合から始まる）。

## 参照値

`-np 4`、gfortran-14 で作った。各反復の終わりの lmf のバンドギャップ（`llmf.1run`、`llmf.2run` の `gap`）:

| | 1 反復目 | 2 反復目 |
|---|---|---|
| GaAs | 0.723 eV | 1.030 eV |
| NiO  | 1.587 eV | 2.066 eV |
