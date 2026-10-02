#準粒子方程式とHedin方程式の導出について書く


## Self-energy

$$
G(x_1, x_2) = G_0(x_1, x_2) + \int dx_3 dx_4  G_0(x_1, x_3) \varSigma(x_3, x_4) G(x_4, x_2)
$$

## Hedin's Equation
$$
\begin{aligned}
P(x_1, x_2) &= -i \int dx_3 dx_4  G(x_4, x_2) G(x_2, x_3) \Gamma(x_3, x_4, x_1) \\
W(x_1, x_2) &= v(x_1, x_2) + \int dx_3 dx_4  v(x_1, x_3) P(x_3, x_4) W(x_4, x_2) \\
\varSigma(x_1, x_2) &= i \int dx_3 dx_4 \, G(x_1, x_4) \Gamma(x_4, x_2, x_3) W(x_1^+, x_3) \\
\Gamma(x_1, x_2, x_3) &= \delta(x_1, x_2) \delta(x_1, x_3) + \int dx_4 dx_5 dx_6 dx_7 
\frac{\delta \varSigma(x_1, x_2)}{\delta G(x_4, x_5)} G(x_4, x_6) G(x_7, x_5) \Gamma(x_6, x_7, x_3)
\end{aligned}
$$
