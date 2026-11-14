# G のスペクトル関数

::: warning 工事中
プログラムは動きますが、このページ（マニュアル）はまだ整っていません。
:::

$\langle {\bf q} n|\Sigma(\omega)|{\bf q} n\rangle$ の計算のしかた

対角要素 $\langle \psi({\bf q},n)|\Sigma_{\rm c}(\omega) |\psi({\bf q},n)\rangle$ を、`gw_lmfh` で one-shot の GW の計算をした後で、
プログラム `hsfp0` のモード `--job=4` で計算します
（`gw_lmfh` は $W$ をファイル `__WVR.*` と `__WVI.*` に残し、`hsfp0 --job=4` がそれを読みます）。
`sigm.*` があれば、`gwsc` と同じように読みます。

計算する $\bf q$ とバンドの番号 $n$ を決める必要があります。これらは `ctrlg.<sname>.toml` の `[gw]` のキーで、`gw_lmfh` と同じものです
（キーの意味は [GW の設定](./gwinput) の `QforGW`、`EMAXforGW`、`--readQforGW`）。

(A) q 点

`QforGW` の q 点を使います（一行に一つ、`2pi/alat` の単位）。
```toml
[gw]
QforGW = """
 0.0 0.0 0.0
 0.5 0.0 0.0
 1.0 0.0 0.0
"""
```
`QforGW` が無いか `QforGWIBZ = true` なら、メッシュ `n1n2n3` のすべての既約な k 点で自己エネルギーを計算します。A(ω) にはこれが要るかもしれません。
シフトしたメッシュの方法のため、メッシュの上に無い q 点も `QforGW` に書けます。

(B) バンドの番号

バンドはエネルギーで選びます。`EMINforGW` と `EMAXforGW`（eV、フェルミエネルギーから）が、自己エネルギーを計算するバンドの範囲です。
```toml
[gw]
EMINforGW = -2.0
EMAXforGW =  3.0
```

(C) エネルギーのメッシュ

自己エネルギーのメッシュは $\omega - E_{\rm F} = 0.01\,{\rm Ry} \times i$（$i=0,\pm1,\pm2,...$）で、実軸の $W$ のメッシュの最高のエネルギーの 2 倍までです。

Σ の虚部は、Im W(ω) と Green 関数の極の畳み込みで与えられることに注意してください
（`[gw].t_sigmaw`（K）が、極の Fermi-Dirac のなましの幅を与えます）。Im W（ω=0 の近く）の分解能は `[gw].HistBin_dw` で決まります。

実行:
```bash
gw_lmfh si -np 24
mpirun -np 24 hsfp0 si --job=4 > lsc_spec
```
SEComg.UP（DN）のファイルができます。ファイルの一行は、`main_hsfp0.f90` で
```
    write(ifoutsec,"(4i5,3f10.6,3x,f10.6,2x,2f16.8,x,3f16.8)")
         iw,itq(i),ip,is, q(1:3,ip), wibz(ip), eqx(i,ip,is),
         (omega(i,iw)-ef)*rydberg(),  hartree*zsec(iw,i,ip)
```
と書かれます。エネルギーは eV です。

**表 1**. `SEComg.UP`（`SEComg.DN`）の列

| 列 | 変数 | 意味 |
|---|---|---|
| 1 | `iw` | ω の番号 |
| 2 | `itq(i)` | バンドの番号 |
| 3 | `ip` | q 点の番号（`QforGW` の順） |
| 4 | `is` | スピンの番号 |
| 5-7 | `q` | q ベクトル（直交座標、2pi/alat の単位） |
| 8 | `wibz` | 既約 BZ での q 点の重み |
| 9 | `eqx` | 固有値（eV、フェルミエネルギーから） |
| 10 | `(omega(i,iw)-ef)*rydberg()` | フェルミエネルギーから測った ω |
| 11, 12 | `hartree*zsec(iw,i,ip)` | 自己エネルギーの実部と虚部 |

一つの q 点の一つのバンドのデータが一つの塊で、塊は 2 行の空行で区切られます。
`QforGW` を変えたら、`gw_lmfh` を走らせ直してください（`QforGW` の q 点での固有関数を作ります）。

メッシュの上に無い `QforGW` の q 点についての注意:

- `EMAXforGW` を設定してください。そのような q では、固有値を解くプログラムは `nband` より少ない状態しか返さず、残りを 1d20 で埋めます。
  上限が無いと、この埋めた状態が $\Sigma$ に入り、`hsfp0` は `sxcf 222: |w-e| out of range` で止まります。
- `EMINforGW` や `EMAXforGW` を変えたら、交換の部分もやり直してください（`gw_lmfh` をはじめから）。交換のファイル（`SEXU`、`XCU`）と相関のファイル（`SECU`）は、同じバンドを持っている必要があります。そうでないと `hqpe` が `SEXU` の終わりで止まります。
- `hgw` と `hx0fp0` の標準出力には、Γ のまわりの小さな q について `epsWVR: iq iw omg eps(wFC) eps(woLFC)` の行があります。
  そこの $\varepsilon(\omega)$ で、別の計算をしなくてもプラズモン（Re $\varepsilon = 0$）が分かります。鋭いプラズモンの構造に重なる準位は実軸のメッシュに敏感なので、
  `HistBin_ratio` を細かくして（例えば 1.03 → 1.015）確かめてください。

* 例:
```
  plot 'SEComg.UP' u ($10):($11) w l,'' u ($10):($12) w l
```
  で、Re Σc(ω) と Im Σc(ω) の図が描けます（表 1 の列）。

## 積分したスペクトル関数（状態密度）
積分したスペクトル関数（状態密度）を得るには、すべてのスペクトル関数（すべての q 点とすべてのバンド）を重ね合わせる必要があります。各 q 点の縮重度（多重度）に注意してください。
SEComg のファイルから組み立てる必要があります。q 点の重みは 8 列目にあります（表 1）。qg4gw の画面の出力（lqg4gw）のキーワード `Multiplicity` の行も見てください。

ともあれ、「それをする価値があるか」を考えてください。
結果を確かめるには、和則（スペクトルの重みの和）を使います。実部と虚部の関係（ヒルベルト変換）にも注意してください。
