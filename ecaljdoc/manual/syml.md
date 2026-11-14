# `syml.<sname>`

バンドプロットの対称ラインのファイルです。
`syml.<sname>` は `getsyml` で作れます。

## バンドプロットの対称ラインとブリルアンゾーンの図を得る
`syml.<sname>` は `ctrlg.<sname>.toml`（古い形式では `ctrl.<sname>`）から作ります。
バンドプロットには `syml.<sname>` が要ります。
作った後で、`job_band` のために `syml.<sname>` を手で直すのも簡単です。

ecalj/GetSyml に getsyml.py があります。これは
[seekpath](https://github.com/giovannipizzi/seekpath/) と [spglib](https://anaconda.org/conda-forge/spglib) にもとづいています。

### 使い方
`InstallAll.py` でインストールするときに、`ecalj/GetSyml/getsyml.py` への symlink `getsyml` を作ります。
次のように実行します。
```
getsyml nio
(or)
getsyml ctrls.nio
```
バンドプロットの対称ラインと一緒に、3 次元のブリルアンゾーンを表示します。
例は [BZsamples](https://ecalj.sakura.ne.jp/BZgetsyml/) にあります。
対称ラインは、ecalj 用の `syml.<sname>` ファイルに書かれます。
`syml` の分割数は大まかな方法で決めているので、必要なら直してください。

### getsyml を使ったときの引用
  論文にするときは、ecalj の引用に加えて、次の引用をお願いします（[ライセンス](./license) を参照）。

   1.Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, 
     Band structure diagram paths based on crystallography,
     Comp. Mat. Sci. 128, 140 (2017)
     
   2.実装に欠かせないライブラリである spglib も引用してください。
     A. Togo, K. Shinohara, I. Tanaka, Sci. Technol. Adv. Mater., Meth. 4, 2384822 (2024).

* spglib と seekpath のライセンス: [ライセンス](./license)。

### （開発者向けのメモ）
   a. lmchk を直して、妥当な値を与えるのに必要な情報を書かせる。
     例えば ndiv（ラインに沿ったメッシュの数）。
   b. 計算の数値精度。
     lmchk の出力を読んでいるので、np.set_printoptions(precision=16) は意味がない。

 対称ラインのファイルは、選んだ対称ラインに沿ってエネルギーバンドを描くとき、
    あるいはフェルミ面のような等エネルギー面を作るときの入力です。
    このファイル（その名前は、"Command-line switches" の節で説明している
    コマンドライン引数 --band の修飾子で指定する）は、いくつかの形をとれます。

### syml の書式
    決まった対称ラインに沿ってバンドを作ります。
    次の例は、単純立方格子の X→Γ と Γ→M のラインの入力です。
       21  .5 0  0    0  0  0         X     Gamma
       21  0  0  0    .5 .5 0         Gamma  M
       0   0  0  0    0  0  0
    最初の数は、各ラインに沿った点の数です。
    次の 6 つは、始点と終点の q 点です。
    最後の行はすべて 0 にします。
