# Samples/AFsymmetry — 反強磁性の対称性（`symgrpaf`）を使う計算

反強磁性体では、「空間の操作 g とスピンの反転」の組が対称操作になる。これを `symgrpaf` に書くと、`lmf` はスピン 1 だけを解き、
スピン 2 の固有関数と密度を対称操作で作る（SCF の計算が約半分になる）。

## 入力

`ctrlg.<sname>.toml` に次の 2 つを書く。

```toml
symgrp   = "find"
symgrpaf = "find"             # AF の操作は spglib が af の印から見つける（2026-10-02）

[[site]]
atom = "Niup"
pos  = [0.0, 0.0, 0.0]
af   = 1                      # 反強磁性の対。スピンを反転すると互いに移るサイトに af = 1 と af = -1 を書く
[[site]]
atom = "Nidn"
pos  = [1.0, 1.0, 1.0]
af   = -1
```

- `symgrpaf = "find"`: spglib が af の印から磁気空間群を求め、時間反転つきの操作（空間の操作 g とスピンの反転の組）で AF のモードになる。
  生成元（`symgrp` と同じ書き方。NiO は `i:( 1 1 1 )`、NiAs 型の NiSe は `r2z:( 0 0 0.5 )`）は従来の探し方（`ECALJ_SYMFIND=ecalj`）のときだけ使い、
  そこで `"find"` と書くと止まる（対を同じ種にまとめると磁気の並びが消えるため）
- `pwmode = 11` にする（スピン 2 の関数をスピン 1 の関数の回転で作るので、q ごとの |q+G| の組が要る。`pwmode = 1` は NiO の 4×4×4 で `rotwave: q+G rotation error`）
- 対になる 2 つのサイトは別の種（`[[spec]]` の `Niup` と `Nidn`）にし、`mmom` の初期値を逆符号にする
- スピン軌道相互作用を完全に入れる計算（`so = 1`）との併用は確かめていない
- LDA+U の `idu` が 10 以上の設定（`idu = [0, 0, 12]` など）と併用すると `lmf` が止まる

## サンプル

**表 1**. サンプルと試験

| ディレクトリ | 中身 | 試験が比べるもの | 時間（8 コア） |
| --- | --- | --- | --- |
| `NiO/` | NiO の LDA。保存してある `rst.nio` から 3 回の SCF | `out.lmf.nio` の全エネルギーなど | 10 秒 |
| `NiSe/` | NiAs 型の NiSe の LDA。同じく 3 回の SCF | `out.lmf.nise` | 10 秒 |
| `NiO_gwsc/` | NiO の QSGW（2×2×2、LDA から 1 反復）。`Samples/TestInstall/nio_gwsc` の入力に `symgrpaf` と `af` を足したもの | `QPU`、`log.nio` | 40 秒 |

```bash
cd Samples/AFsymmetry
testecalj NiO NiSe NiO_gwsc -np 8
```

## QSGW での扱い

`symgrpaf` を書いた QSGW では、SCF（`lmf`）は対称性を使ってスピン 1 だけを解き、GW のプログラム（`hgw`、`hsfp0_sc`、`hqpe_sc`）は
両方のスピンを計算する。結果は `symgrpaf` を書かない計算と一致する（表 2）。

**表 2**. NiO、QSGW 2 反復（2×2×2）。`symgrpaf` あり（`NiO_gwsc` の入力）と、なし（`TestInstall/nio_gwsc` の入力）

| | LDA のギャップ（eV） | 1 反復後 | 2 反復後 |
| --- | --- | --- | --- |
| `symgrpaf` あり | 1.0253 | 2.0018 | 2.9862 |
| `symgrpaf` なし | 1.0253 | 2.0016 | 2.9845 |

2 反復後の `QPU` の差は、SEx・SEc で最大 0.003 eV、QP エネルギーのずれ（`dSEnoZ`）で最大 0.001 eV。`symgrpaf` ありでは `QPU` と `QPD`（2 つのスピン）の
QP エネルギーが一致する（gfortran、4 コア、2026-09-30）。

2026-09-29 までのコードでは、`symgrpaf` を書いた QSGW は 1 反復目の `hqpe_sc` で止まった（2026-09-30 に直した）。
