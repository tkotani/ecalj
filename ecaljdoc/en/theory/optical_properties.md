# Optical properties

::: warning Unfinished
This page is unfinished.
:::

> English version of the page [theory/optical_properties](../../theory/optical_properties.md). The pages at the top of ecaljdoc are the reference; this one was made from theory/optical_properties.md of commit c173e3b69 (source blob 9c39003eda09).

## Dielectric function on the microscopic scale

$$
\begin{align}
𝛿𝜌 = 𝜒𝛿𝜙_{ext} = 𝑃𝛿𝜙 \\
𝜖𝛿𝜙 = 𝛿𝜙_\mathrm{ext} \\
𝛿𝜙 = 𝜖^{-1}𝛿𝜙_\mathrm{ext}
\tag{O1}
\end{align}
$$
where:
- $𝛿𝜌$: Change in charge density
- $𝜒$: density response function to the external potential (the full response)
- $𝛿𝜙_{ext}$: External potential change
- $𝑃$: polarization function (the response to the total potential)
- $𝜖$: dielectric function
- $𝛿𝜙 = 𝛿𝜙_\mathrm{ext} + 𝛿𝜙_\mathrm{ind}$: total potential change, which is sum of changes in the external potential and induced one
::: info
Short notation is used in the equations above, for example
$$
(𝜒𝛿𝜙_{ext})(𝒓,𝜔)= \int d𝒓'' 𝜒(𝒓,𝒓'',𝜔)𝛿𝜙_{ext}(𝒓'',𝜔)
\tag{O2}
$$
:::
From these equations, we obtain the following equations straightforwardly
$$
\begin{align}
𝜖 &= 1 - 𝑣𝑃  \\
𝜖^{-1} &= 1 + 𝑣 𝜒 \\
𝜒 &= 𝑃 + 𝑃 𝑣 𝜒
\tag{O3}
\end{align}
$$
where:
- $𝑣$: Coulomb potential

## Random phase approximation (RPA)
$$
𝑃 = 𝛱
\tag{O4}
$$
where
 - $𝛱$: Kohn-Sham density response function

$$
 \epsilon(\mathbf{r}, \mathbf{r}', \omega) = \delta(\mathbf{r} - \mathbf{r}') - \int d\mathbf{r}'' \, v(\mathbf{r} - \mathbf{r}'') 𝛱(\mathbf{r}'', \mathbf{r}', \omega)
\tag{O5}
$$

## Macroscopic dielectric function
$$
\hat{𝒒}⋅𝜖_\mathrm{mac}(𝒒,𝜔)⋅\hat{𝒒}= 1/[𝜖^{-1}(𝒒,𝜔)]_{0,0}
\tag{O6}
$$
Here $0,0$ means the component $𝑮=𝑮'=0$.

## Without the local-field correction
$$
𝜖_\mathrm{mac}(𝒒,𝜔) ≃ 𝜖_{0,0}(𝒒,𝜔)
\tag{O7}
$$
In this case the macroscopic dielectric function is obtained as
$$
𝜖_{0,0}(𝒒,𝜔) = 1 - 𝑣_0(𝒒)𝑃_{0,0}(𝒒,𝜔)
\tag{O8}
$$
Only the $0,0$ component of $𝑃$ is needed, so the computational cost is low.

## With the local-field correction
Eq. (O6) is used as it is. In this case all the matrix elements of $𝜖_{μν}(𝒒,𝜔)$ are needed to calculate the inverse matrix.
See the [Adler](https://dx.doi.org/10.1103/physrev.126.413) for the details.

## Relation to optical constants
- $𝑛$: refractive index
$$
𝑛 = \sqrt{ \frac{ 𝜖_\mathrm{mac,1} + \sqrt{𝜖^2_\mathrm{mac,1} + 𝜖^2_\mathrm{mac,2}}}{2}} 
\tag{O9}
$$
- $𝜅$: extinction coefficient
$$
𝜅 = \sqrt{ \frac{-𝜖_\mathrm{mac,1} + \sqrt{𝜖^2_\mathrm{mac,1} + 𝜖^2_\mathrm{mac,2}}}{2}}
\tag{O10}
$$
- $𝑅$: reflectance
$$
𝑅 = \frac{( 𝑛 - 1)^2 + 𝜅^2} {( 𝑛 + 1)^2 + 𝜅^2}
\tag{O11}
$$

$𝜖_\mathrm{mac,1}$ and $𝜖_\mathrm{mac,2}$ are the real and imaginary parts of $𝜖_\mathrm{mac}$.
