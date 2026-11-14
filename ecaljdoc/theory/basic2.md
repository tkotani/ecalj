# 基本2

::: warning 工事中
このページは書きかけです。
:::

- 時間秩序, 先進, 遅延の関係について
- ヒルベルト変換/クラマース・クロニッヒ関係について

Green関数やそれに付随する物理量は時間秩序化されており、複素関数としての解析性は悪い。
Wickの定理が使用でき摂動計算では有用であるため、時間秩序化されたものを用いるが、解析性が良い関数にした方が良いこともある。
例えば、クラマース・クロニッヒ関係は, 実部(虚部)の情報から虚部(実部)を構築でき、数値計算上も計算を効率的に行う技術として使えるが、上半面も下半面でも正則でない時間秩序Green関数はそのままではKK関係が使用できない。
また、KK関係は数値計算精度を確認するバロメーターとしても使えるので、解析性の良い遅延・先進関数との関係を調べておく。

## ヒルベルトの公式
$$
\lim_{𝜂 → 0 }\int_{-∞}^{∞} \frac{𝜒(ω')}{ω' - ω ± i 𝜂} dω' = 𝒫 \int_{-∞}^{∞} \frac{𝜒(ω')}{ω' - ω } dω' ∓ iπ𝜒(ω)
\tag{K1}
$$

## 先進, 遅延, 時間秩序の関係
$\mathrm{R}$は遅延, $\mathrm{A}$は先進を示す。

### グリーン関数の場合
時間秩序化されていない場合, 場の演算子の反交換関係のプロパゲータが遅延(先進)グリーン関数。

$$ \begin{align}
iG(𝒓_1t_1, 𝒓_2t_2) &= \braket{\Psi_0|T[\hat 𝜓_\mathrm{H}(𝒓_1,t_1)\hat 𝜓^†_\mathrm{H}(𝒓_2,t_2)]|\Psi_0}  \\
iG^\mathrm{R}(𝒓_1t_1, 𝒓_2t_2) &=  \braket{\Psi_0|\{\hat 𝜓_\mathrm{H}(𝒓_1,t_1),\hat 𝜓^†_\mathrm{H}(𝒓_2,t_2)\}|\Psi_0} 𝜃(t_1-t_2) \\
iG^\mathrm{A}(𝒓_1t_1, 𝒓_2t_2) &= -\braket{\Psi_0|\{\hat 𝜓_\mathrm{H}(𝒓_1,t_1),\hat 𝜓^†_\mathrm{H}(𝒓_2,t_2)\}|\Psi_0} 𝜃(t_2-t_1) \\
\tag{K2}
\end{align} $$
フーリエ表示は、
$$ \begin{align}
G(𝒓_1, 𝒓_2;𝜔)  &=  \sum_{n'}^{\mathrm{unocc}} \frac{𝜙_{n'}(𝒓_1)𝜙^*_{n'}(𝒓_2)}{𝜔-𝜀_{n'}/\hbar + i𝛿} +
    \sum_{n}^{\mathrm{occ}} \frac{𝜙_{n}(𝒓_1)𝜙^*_{n}(𝒓_2)}{𝜔 - 𝜀_{n}/\hbar -i𝛿} \\ &= G_{\mathrm{ele}} + G_{\mathrm{hole}} \\
G^\mathrm{R}(𝒓_1, 𝒓_2;𝜔)  &=  \sum_{n'}^{\mathrm{unocc}} \frac{𝜙_{n'}(𝒓_1)𝜙^*_{n'}(𝒓_2)}{𝜔-𝜀_{n'}/\hbar + i𝛿} +
    \sum_{n}^{\mathrm{occ}} \frac{𝜙_{n}(𝒓_1)𝜙^*_{n}(𝒓_2)}{𝜔 - 𝜀_{n}/\hbar +i𝛿} \\ &= G_{\mathrm{ele}} + G^†_{\mathrm{hole}} \\
G^\mathrm{A}(𝒓_1, 𝒓_2;𝜔)  &=  \sum_{n'}^{\mathrm{unocc}} \frac{𝜙_{n'}(𝒓_1)𝜙^*_{n'}(𝒓_2)}{𝜔-𝜀_{n'}/\hbar - i𝛿} +
    \sum_{n}^{\mathrm{occ}} \frac{𝜙_{n}(𝒓_1)𝜙^*_{n}(𝒓_2)}{𝜔 - 𝜀_{n}/\hbar -i𝛿} \\ &= G^†_{\mathrm{ele}} + G_{\mathrm{hole}}
\tag{K3}
\end{align} $$
となる。$†$ は複素共役をとって$𝒓_1, 𝒓_2$をいれかえるという意味。
$G^\mathrm{R}$, $G^\mathrm{A}$は$ω$面での極の位置が、下半面(遅延)、上半面(先進)と分離されていてクラマース・クロニッヒ(KK)関係を満たす。
$$
\Re G^\mathrm{R,A} (𝒓_1, 𝒓_2;𝜔) = ± 𝒫 \int_{-∞}^{∞} \frac{dω'}{π} \frac{ \Im G^{\mathrm{R,A}}(𝒓_1, 𝒓_2;𝜔') }{ω' - ω}
\tag{K4}
$$
である。

### 密度応答の場合
密度応答関数は、時間反転対称性がある。

$$ \begin{align}
iD(𝒓_1t_1, 𝒓_2t_2) &= \braket{T[\hat {𝛿𝜌}(𝒓_1,t_1)\hat {𝛿𝜌}(𝒓_2,t_2)]} \\ &= G(𝒓_1t_1,𝒓_2t_2)G(𝒓_2t_2,𝒓_1t_1) \\
iD^\mathrm{R}(𝒓_1t_1, 𝒓_2t_2) &= \braket{[\hat {𝛿𝜌}(𝒓_1,t_1),\hat {𝛿𝜌}(𝒓_2,t_2)]} 𝜃(t_1-t_2) \\
iD^\mathrm{A}(𝒓_1t_1, 𝒓_2t_2) &= -\braket{[\hat {𝛿𝜌}(𝒓_1,t_1),\hat {𝛿𝜌}(𝒓_2,t_2)]} 𝜃(t_2-t_1)
\tag{K5}
\end{align} $$

先進関数には $G^\mathrm{A}$ と同じく負号が付く．

> [!INFO]
> フェッター ワレッカの章末問題
>  $\Re D(ω) = \Re D^\mathrm{R}(ω)$ と $\Im D(ω) \mathrm{sgn}(ω) = \Im D^\mathrm{R}(ω)$

### $GW$自己エネルギーの場合

$$ \begin{align}
𝛴(𝒓_1t_1, 𝒓_2t_2) &=iG(𝒓_1t_1, 𝒓_2t_2)W(𝒓_1t_1, 𝒓_2t_2) \\
                  &= iG_\mathrm{ele}(𝒓_1t_1, 𝒓_2t_2)W(𝒓_1t_1, 𝒓_2t_2) \\
                  &\quad + iG_\mathrm{hole}(𝒓_1t_1, 𝒓_2t_2)W(𝒓_1t_1, 𝒓_2t_2)
\tag{K6}
\end{align} $$
$$ \begin{align}
𝛴_{nn}(𝜔) &= (G_\mathrm{ele}W)(𝜔)_{nn} + (G_\mathrm{hole}W)(𝜔)_{nn} = 𝛴_\mathrm{ele}^{nn}(𝜔)+𝛴_\mathrm{hole}^{nn}(𝜔) \\
𝛴^\mathrm{R}_{nn}(𝜔) &= 𝛴_\mathrm{ele}^{nn}(𝜔)+(𝛴_\mathrm{hole}^{†})^{nn}(𝜔) \\
\tag{K7}
\end{align} $$

$G_\mathrm{ele}$ はすでに $𝜃(t_1-t_2)$ を, $G_\mathrm{hole}$ は $𝜃(t_2-t_1)$ を含むので, $𝛴_\mathrm{ele}=iG_\mathrm{ele}W$ は $t_1>t_2$ だけ, $𝛴_\mathrm{hole}=iG_\mathrm{hole}W$ は $t_1\lt t_2$ だけに値を持つ（$𝛴_\mathrm{ele}(𝜔)$ の極は下半面, $𝛴_\mathrm{hole}(𝜔)$ の極は上半面）．$G^\mathrm{R}$ と同じく正孔の部分を $†$ にしたものが遅延自己エネルギーで, 絶対零度で $𝜔$ をフェルミ準位から測ると $\Re𝛴^\mathrm{R}_{nn}(𝜔)=\Re𝛴_{nn}(𝜔)$, $\Im𝛴^\mathrm{R}_{nn}(𝜔)=\Im𝛴_{nn}(𝜔)\,\mathrm{sgn}(𝜔)$ である．

## ヒルベルト変換

### 数値計算
主値積分があるが、基本は通常の数値積分方を使う。ただし発散点を数値的に扱わないようにする。
データ点は$ω_{C_i}$にあるとして, 近接のデータ点までの中間地点を $ω_{L_i}$と$ω_{R_i}$とする。
- 第ゼロ近似

$$ \begin{align}
F(ω) &= 𝒫 \int_{-∞}^{∞} dω' \frac{f(ω')}{ω' -ω} \\ &= \sum_{i} 𝒫 \int_{ω_{L_i}}^{ω_{R_i}} dω' \frac{f(ω')}{ ω' -ω} \\
& \simeq \sum_{i} f(ω_{C_i})𝒫 \int_{ω_{L_i}}^{ω_{R_i}} dω' \frac{1}{ω'-ω} \\
&= \sum_{i} f(ω_{C_i}) \log{\frac{|ω_{R_i} -ω|} {| ω_{L_i}-ω|}}
\tag{K8}
\end{align} $$

ここで $ω_{L_i}\lt ω_{C_i}\lt ω_{R_i}$（$ω_{L_i}$ が左, $ω_{R_i}$ が右の中間点）．

 通常は虚部(実部)と同じメッシュ点に対して実部(虚部)を計算する。
