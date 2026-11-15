# ドーピングのしかた

第一原理計算は、電荷が中性の条件で行います。
単位胞の中の電子の数を整数でない値にするには、二つの方法（または両方の組み合わせ）があります。

 1. 原子番号 `[[spec]].z` をフラクショナルな値（整数でない値）にする。
 2. 背景電荷 `[bz].zbak` で価電子の数を決める。

以下の例は Si（Z=14、単位胞に 2 原子）のドーピングです。

## フラクショナルな Z

例えば Z=14.2 とします。Z と一緒に、原子の `[[spec]].q` も書く必要があります。
Z=14.2 なら `q = [2, 2.2]` などとします（`q` は s、p、d、f の順の価電子の初期の数です）。
つまり `q = [2, 2.2]` は、球対称な原子の密度を計算するときに原子が中性になるよう、p に 0.2 個の電子を足すという意味です。

`lmfa` を `grep conf` と一緒に走らせると、初期の電子の配置が分かります。
`lmfa` の出力に Q=0（電荷の中性の確認）が出ます。`lmfa` が正しく終わったことを確かめてください。

`lmf` を走らせ、画面の出力を `llmf` に残します。
```bash
 grep 'z=' llmf
```
で原子の Z が分かります。
```bash
 grep -A5 'Charges:' llmf
```
で
```
 Charges:  valence     8.40000   cores    20.00000   nucleii   -28.40000
    hom background     0.00000   deviation from neutrality:     -0.00000
```
が出ます。原子核は $-2 \times 14.2 = -28.4$（Z から）です。

## 背景電荷 zbak

`[bz].zbak` を設定します。`zbak = 0.4` なら、背景電荷は −0.4|e| で、電子の数は 0.4 少なくなります
（コア + 価電子 + zbak = |原子核の全電荷|）。
```
Charges:  valence     7.60000   cores    20.00000   nucleii   -28.00000
    hom background     0.40000   deviation from neutrality:     -0.00000
```
つまり 7.6 + 0.4 + 20.0 = 28.0 です。`zbak = 0.2` なら valence は 7.8、cores は 20.0 になります。

## フラクショナルな Z と zbak の組み合わせ

Z=14.2、`q = [2, 2.2]` と、`zbak = 0.3` を一緒に使うと、
```
   Charges:  valence     8.10000   cores    20.00000   nucleii   -28.40000
    hom background     0.30000   deviation from neutrality:     -0.00000
```
となります。単位胞に 2 原子なので、
- 原子核 = $-2 \times 14.2 = -28.4$（Z から）
- コア = $2 \times 10 = 20$（Z の整数部分の原子から）
- 背景 = 0.3（zbak から）

で、価電子 = 28.4 − 20 − 0.3 = 8.1 です。

## 最初の反復での注意 {#background-charge-and-fractional-z}

最初の反復では、`Charges:` は例えば
```
  Charges:  valence     8.00000   cores    20.00000   nucleii   -28.00000
  hom background     0.12300   deviation from neutrality: 0.12300
```
のように出ます。初期の密度を原子の密度の重ね合わせで作るためで、中性からのずれが 0 でないように見えます。
次の反復からは `deviation from neutrality: 0` になります（電荷は保存されます）。
