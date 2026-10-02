# Optical Properties

## dielectric function in micro-scale

$$
\begin{align}
𝛿𝜌 = 𝜒𝛿𝜙_{ext} = 𝑃𝛿𝜙 \\
𝜖𝛿𝜙 = 𝛿𝜙_\mathrm{ext} \\
𝛿𝜙 = 𝜖^{-1}𝛿𝜙_\mathrm{ext}
\end{align}
$$
where:
- $𝛿𝜌$: Change in charge density
- $𝜒$: Electric susceptibility
- $𝛿𝜙_{ext}$: External potential change
- $𝑃$: Polarization function
- $𝜖$: dielectric function
- $𝛿𝜙 = 𝛿𝜙_\mathrm{ext} + 𝛿𝜙_\mathrm{ind}$: total potential change, which is sum of changes in the external potential and induced one
::: info
the short formula is used above equations 
$$
(𝜒𝛿𝜙_{ext})(𝒓,𝒓',𝜔)= \int d𝒓'' 𝜒(𝒓,𝒓'',𝜔)𝛿𝜙_{ext}(𝒓'',𝒓',𝜔)
$$
:::
From these equations, we obtain the following equations straightforwardly
$$
\begin{align}
𝜖 &= 1 - 𝑣𝑃  \\
𝜖^{-1} &= 1 + 𝑣 𝜒 \\
𝜒 &= 𝑃 + 𝑃 𝑣 𝜒
\end{align}
$$
where:
- $𝑣$: Coulomb potential

## RPA Approxiamation:
$$
𝑃 = 𝛱
$$
where
 - $𝛱$: Kohn-Sham density response function

$$
 \epsilon(\mathbf{r}, \mathbf{r}', \omega) = \delta(\mathbf{r} - \mathbf{r}') - \int d\mathbf{r}'' \, v(\mathbf{r} - \mathbf{r}'') 𝛱(\mathbf{r}'', \mathbf{r}', \omega)
$$

## macroscale dielectric function
$$
\hat{𝒒}⋅𝜖_\mathrm{mac}(𝒒,𝜔)⋅\hat{𝒒}= 1/[𝜖^{-1}(𝒒,𝜔)]_{0,0}
$$

## without local field correction
$$
𝜖_\mathrm{mac}(𝒒,𝜔) ≃ 𝜖_{0,0}(𝒒,𝜔)
$$
in this case, the macroscoptic dielectric function is obtained by following way
$$
𝜖_{0,0}(𝒒,𝜔) = 1 - 𝑣_0(𝒒)𝑃_{0,0}(𝒒,𝜔)
$$
only, $0,0$ component of $𝑃$ is needed, being the low computational cost.

## local field correction
$$
\hat{𝒒}⋅𝜖_\mathrm{mac}(𝒒,𝜔)⋅\hat{𝒒}= 1/[𝜖^{-1}(𝒒,𝜔)]_{0,0}
$$
in this case, all matrix element of $𝜖_{μν}(𝒒,𝜔)$ is needed to calculate the inverse of matrix.
See the [Adler](https://dx.doi.org/10.1103/physrev.126.413) for the details.

## relation with optical functions
- $𝑛$: refractive index
$$
𝑛 = \sqrt{ \frac{ 𝜖_\mathrm{mac,1} + \sqrt{𝜖^2_\mathrm{mac,1} + 𝜖^2_\mathrm{mac,2}}}{2}} 
$$
- $𝜅$: extinction coeffcient
$$
𝜅 = \sqrt{ \frac{-𝜖_\mathrm{mac,1} + \sqrt{𝜖^2_\mathrm{mac,1} + 𝜖^2_\mathrm{mac,2}}}{2}} \\
$$
- $𝑅$: reflectance
$$
𝑅 = \frac{( 𝑛 - 1)^2 + 𝜅^2} {( 𝑛 + 1)^2 + 𝜅^2}
$$
