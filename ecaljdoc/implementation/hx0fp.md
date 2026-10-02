# Implementation of $W$ with RPA response function

> ⚠️ **TOML migration (2026-05)** — `hx0fp0` / `hrcxq` now read `ctrlg.<sname>.toml` only. See [TOML migration](../manual/toml_migration).

## 実行プログラム
- `hrcxq` : QSGW で使用
- `hx0fp0` : 誘電関数等で使用

## 並列化
以下のMPI並列化が実装されている。
 - $q$点: 第一還元ゾーン中のk点数。 offset Γ 処理を含む。
 - $k$点: 第一ゾーン中のk点数 `[gw].n1n2n3` で指定した量 (legacy: `n1n2n3` in `GWinput`)
 - MPB: 積基底

> [!important]
> デフォルトでは `hrcxq` では q点並列とMPB並列が利用され、`hx0fp0` では q点並列とk点並列が使用されるが、3軸並列を同時に使うことも可能である。その場合は、以下の実行時引数で並列数を指定する必要がある。
> - MPB並列数: `--nb=XX`
> - k点並列数: `--nk=YY`

## $𝛱(ω)$の平滑化

$\Im 𝛱(ω)$ を計算した後に, 振動数について平滑化を行うオプションがある。
平滑化には以下のガウシアンフィルターを使用する． `GaussianFilterX0` を検索。
$$
\begin{align}
\Im \overline{𝛱}(ω_i) &= \sum_j G(𝜔_i, 𝜔_j){\Im 𝛱}(ω_j) \\
G(𝜔_i, 𝜔_j) &= e^{-\frac{(𝜔_i - 𝜔_j)^2}{2𝜎^2}}/\sum_{k} e^{-\frac{(𝜔_i - 𝜔_k)^2}{2𝜎^2}}
\end{align}
$$

> [!NOTE]
> - $𝛱(ω)$は通常は$ω>0$のみ計算している。
> - 範囲外の寄与はゼロパディングを用いてる。

<!-- ## $W^c(𝒒,𝜔)$の平滑化 -->
<!-- 誘電関数$𝜀_{𝜇𝜈}(𝒒,𝜔) ≃ 0$となる場合はプラズマ振動が発生し$W^c = (1/𝜀 - 1)v$に異常が表れるが数値計算上の有限の$𝜔$メッシュでは,  -->
<!-- $$ -->
<!-- \begin{align} -->
<!-- 𝜀_{𝜇𝜈}(𝒒,𝜔) = 𝛿_{𝜇𝜈}(1 - iη) - \sqrt{v_𝜇(𝒒)} {𝛱}_{𝜇𝜈}(𝒒,ω)\sqrt{v_𝜈(𝒒)} -->
<!-- \end{align} -->
<!-- $$ -->
<!---->
