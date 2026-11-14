# 力と原子位置の緩和 {#forces-and-atomic-position-relaxation}

LDA/GGA で原子の位置を緩和できます（単位胞の最適化はしません）。サンプル: [`Samples/Relax/LaGaO3`](https://github.com/tkotani/ecalj/tree/main/Samples/Relax/LaGaO3)（直方晶のペロブスカイト、20 原子）。

`ctrlg.<sname>.toml` の `[dyn]` の節で設定します（表 1）。`[dyn].mode` が 0 でなければ、力も計算します。
緩和させる方向は、各 `[[site]]` に `relax = [0, 0, 1]` のように書きます（1 の方向だけ動かす）。

表 1. `[dyn]` のキー（既定値は `m_lmfinit.f90` のもの）

| キー | 既定 | 意味 |
|---|---|---|
| `mode` | 0 | 0: 緩和しない、4: 共役勾配法、5: Fletcher-Powell、6: Broyden |
| `nit` | 1 | 緩和の段の数の上限 |
| `hess` | true | ヘッセ行列を読む |
| `xtol` | 1e-3 | 変位の収束の条件。>0: 長さを使う、<0: 最大値を使う、=0: 使わない |
| `gtol` | 0 | 勾配の収束の条件。>0: 長さを使う、<0: 最大値を使う、=0: 使わない |
| `step` | 0.015 | 最初の（そして最大の）一歩の長さ |
| `nkill` | 0 | この段の数の後でヘッセ行列を捨てる |

LaGaO3 のサンプルの値:
```toml
[dyn]
mode = 5
nit  = 20
hess = true
xtol = 0.0001
gtol = 0.0
step = 0.015
```

* 注意: 緩和のモードでは、原子の位置のファイル `AtomPos.<sname>` を使います。このファイルがあれば、`ctrlg.<sname>.toml` の位置の代わりに、そこに書かれた位置を使います
  （`SRC/subroutines/main_lmf.f90` の `call ReadAtomPos(irpos)`）。新しく始めるときは、古い `AtomPos.<sname>` を消してください。
* ペロブスカイトの LaGaO3 の緩和がうまくいくことは確かめています。ただ、緩和のモードはあまり使っていないので、注意して使ってください。
