# Effective Screening Medium（ESM） {#esm-effective-screening-medium}

真空層のあるスラブは周期的な結晶ではありませんが、平面波（なめらかな密度）の計算は、Poisson 方程式を周期的なものとして解きます。
ESM（Otani & Sugino, [PRB 73, 115407 (2006)](https://doi.org/10.1103/PhysRevB.73.115407)。ecalj への実装は金沢大学の小幡正雄による）は、
静電の周期的な Green 関数を、スラブの両側に本当に欲しい境界条件（真空、または金属の電極）を持つものに置き換えます。
スラブに電場を掛けることもできます。ESM と QSGW を組み合わせた例は
[PRB 101, 205120 (2020)](https://doi.org/10.1103/PhysRevB.101.205120) にあります。

**真空層のあるスラブでは、必ず ESM を使ってください。** 使わないと、クーロンポテンシャルの $G=0$ 成分が真空の準位でなく周期的な約束で決まり、
固有値の全体が別のゼロ点の上に乗ります。Fe/MgO のサンプル（`Samples/MLOsamples/FeMgO`、真空 30.5 a.u.）では **4.4 eV** ずれますが、
全エネルギー、`Vesav`、磁気モーメントは 8〜13 桁一致するので、フェルミエネルギーを比べるまで何もおかしく見えません。

`ctrlg.<sname>.toml` の `[esm]` の節で設定します（下に、すべてのキー、`boundary` の値、電場の掛け方、
やめた `esm_input.dat` の `ctrlg_absorb.py` による変換）。サンプル: `Samples/MLOsamples/FeMgO`。

## `[esm]` の節 {#the-esm-section}

```toml
[esm]
boundary  = "vac/slab/vac"  # vacuum(-z)/slab/vacuum(+z)
origin    = -8.63717        # (a.u.) z-translation of the density
shiftmode = 0               # 0: origin is absolute, 1: in units of the cell length
zb        = [23.907104, -23.907104]   # (a.u.) boundaries z1esm, z2esm
potential = [0.0, 0.0]      # (Ry) on the +z / -z sides
field     = [0.0, 0.0]      # (Ry/a.u.) on the +z / -z sides
```

**`[esm]` の節が無いことは、ESM を使わないことを意味します。** バルクの結晶ではそれが正しく、スラブでは誤りです。
そこで `lmf` は、`c` が `a`、`b` の長い方の 3 倍より大きく、`[esm]` が無いときに警告を出します。

`ctrlgenToml.py` は、新しいファイルのすべてに、コメントにした `[esm]` の雛形（`# === ESM` の見出しと、キーごとに一行）を `[ham]` の直後に書きます。
スラブでは、コメントを外して数値を注意して設定してください。バルクの結晶ではコメントのままにします（`[esm]` が無ければ ESM は使わない）。

表 1. `[esm]` のキー

| キー | 意味 | 古い `esm_input.dat` |
|---|---|---|
| `boundary` | 両側の境界条件（表 2） | 1 行目、`jesm` |
| `origin` | 解く前に密度に掛ける z 方向の平行移動。コードの中では `-origin` を使うので、古いファイルの **2 行目と同じ数**です | 2 行目、`tresm` |
| `shiftmode` | `0`: `origin` は a.u.、`1`: 単位胞の長さを単位とする | 2 行目、`jtresm` |
| `zb` | 二つの境界面 `z1esm`、`z2esm`（a.u.）。省くと `±alat·plat[2][2]/2` | 3 行目 |
| `potential` | +z / −z の側に固定するポテンシャル（Ry）。`metal/...` と `:v-e` / `:e-v` の形で使う | 4 行目 |
| `field` | +z / −z の側に固定する電場（Ry/a.u.）。`:field` の形で使う | 5 行目 |

### `boundary` の値

表 2. `boundary` の値

| 値 | 系 | `jesm` |
|---|---|---|
| `"off"`（または節が無い） | ESM を使わない。ふつうの周期境界 | 0 |
| `"vac/slab/vac"` | 両側が真空 | 1 |
| `"metal/slab/metal"` | 両側が金属の電極 | 2 |
| `"vac/slab/metal"` | −z が真空、+z が金属 | 3 |
| `"metal/slab/vac"` | −z が金属、+z が真空 | 4 |
| `"vac/slab/vac:field"` | 両側が真空で、`field` を掛ける | 5 |
| `"metal/slab/metal:v-e"` | −z に `potential`、+z に `field` | 6 |
| `"metal/slab/metal:e-v"` | −z に `field`、+z に `potential` | 7 |
| `"periodic:esm"` | `off` と同じ結果を、ESM の経路で計算する | 10 |

整数そのものも受け付けるので、機械的に変換したファイルもそのまま動きます。ただし `1` と `3` は見ただけでは区別できないので、文字列の方を使ってください。

**スラブに電場を掛けるには**、`"vac/slab/vac:field"` で `field = [E, -E]` とするか、`"metal/slab/metal"` で `potential = [V1, V2]` として二つの電極の電位を固定します。

## `esm_input.dat` からの移行 {#migration-from-esm-input-dat}

2026-09-16 より前は、ESM を別の位置引数のファイル `esm_input.dat`（キーの無い数値の 6 行）で設定していました。
これには二つの問題がありました。`esmsmves.f90` を読まないと数値の意味が分からないこと、そして**ファイルが無いと ESM が黙って切れる**こと（`lmf` は一行書いて続けていた）です。
そのため、その一つのファイルを忘れてサンプルのディレクトリを写すと、エラー無しに物理が変わっていました。

いまの Fortran は `ctrlg.<sname>.toml` だけを読みます（2026-09-17 から）。`esm_input.dat` が残っていると、
`lmf` / `lmfa` / `lmchk` は変換の道具の名前を出して**止まります**。無視せずに止めるのは、スラブがうっかり ESM 無しで走らないためです。

```bash
ctrlg_absorb.py <sname>      # esm_input.dat -> [esm] in ctrlg.<sname>.toml
                             # (and PB.<sname>.toml -> [product_basis], if present)
```

表 3. `ctrlg_absorb.py` のすること

| 状況 | `ctrlg_absorb.py` のすること |
|---|---|
| `ctrlg.<sname>.toml` にすでに `[esm]` がある | TOML の方を使う。`esm_input.dat` は `esm_input.dat.bk` に移し、その冒頭に「この設定は使われなかった」と記録する |
| TOML に `[esm]` が無い | ファイルを変換して `[esm]` の部分を（コメント付きで、GW の節の前に）入れ、元のファイルを `esm_input.dat.bk` に移して、冒頭に移した先を書く |

`Legacy2toml.py` も、古いディレクトリを移すときに同じ変換をします。

## 例: Fe/MgO のスラブ

`Samples/MLOsamples/FeMgO`（`c/a = 8.50`、真空 30.5 a.u.、それを埋める空の球 7 個）:

```toml
[esm]
boundary  = "vac/slab/vac"
origin    = -8.63717
shiftmode = 0
zb        = [23.907104, -23.907104]   # alat*plat[2][2]/2 = 23.903861; set by hand here
potential = [0.0, 0.0]
field     = [0.0, 0.0]
```

同じ `rst.<sname>` で、この節があるときと無いときを比べると表 4 のとおりです。

表 4. `[esm]` の有無（FeMgO）

| | E_F (Ry) | Vesav (eV) | `val*vef` (Ry) |
|---|---|---|---|
| `[esm]` あり | −0.38286150 | −2.12295183 | −257.64 |
| `[esm]` なし | −0.05833040 | −2.12295183 | −253.66 |

`Vesav` は 8 桁一致するので、密度と静電は同じです。動くのはエネルギーのゼロ点だけで、0.3245 Ry = 4.4155 eV です。
どちらだったかは `llmf_ef`（または `lmf` のどのログでも）で確かめられます。

```
effective screening medium method jesm=  1        <- ESM on
system : vaccum(-z)/slab/vaccum(+z)

esmsmves: ESM is off (no [esm] section in ctrlg.<sname>.toml)   <- ESM off
```

## 実装の場所 {#where-it-is-implemented}

`SRC/subroutines/esmsmves.f90`（解く部分。`smves` から呼ぶ）と、`SRC/subroutines/m_lmfinit.f90` の `readesm`（`[esm]` を読む。`esm_input.dat` が残っていれば止める）。
