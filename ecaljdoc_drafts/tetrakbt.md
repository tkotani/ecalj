# Finite-temperature polarization: `tetrakbt`

*(draft for ecaljdoc; target location: manual/gw section)*

## Purpose

For metals, the tetrahedron integration of the polarization function
χ⁰(q,ω) resolves the Fermi surface sharply. When the Fermi surface has
strong nesting (quasi-1D bands, flat sheets), the small-q response becomes
nearly singular and the QSGW self-consistency cycle may oscillate: each
iteration shifts the Fermi surface, W changes sharply, and Σ overreacts.

`tetrakbt` evaluates χ⁰ with Fermi–Dirac occupations at a finite electronic
temperature T. This broadens the Fermi surface physically (Landau-damping
of the sharp response) and is the recommended regularization for such
instabilities. The ω-smearing option `SmearX0` only smooths Im χ⁰ along
frequency and does not remove the underlying k-space sharpness.

## Input keys (`ctrlg.<sname>.toml`, section `[gw]`)

| key | type | default | meaning |
|---|---|---|---|
| `tetrakbt` | logical | `false` | use the finite-T tetrahedron for χ⁰ |
| `t_tetrakbt` | real | `300.0` | electronic temperature in Kelvin |

kBT[eV] ≈ T[K]/11604. Example:

```toml
[gw]
tetrakbt   = true
t_tetrakbt = 2000.0    # kBT = 0.172 eV
```

## What it does

- χ⁰ occupation factors become f(ε−E_F;T) inside the tetrahedron weights
  (energy-convolution method B′: the exact T=0 routine is convolved with the
  thermal kernel −f′(E); exact T→0 limit, sum rules preserved).
- `heftet` (in the GW chain) determines the finite-T Fermi level by bisection
  on the discrete k-sum of f(ε−E_F;T) and writes it to the file `EFERMI_kbt`.
  The χ⁰ tetrahedron reads it, so occupations and E_F refer to the same T.
- The self-energy side (poles of G, `esmr` smearing, exchange) remains at
  T = 0. (A finite-T Σ is under development and not part of releases.)

## Choosing the temperature

kBT should cover the energy spacing of your k-mesh near the Fermi surface,
roughly kBT ≳ v_F·|b|/N for an N×N×N mesh. Then results become insensitive
to the mesh (two meshes agree at fixed T), and the QSGW iteration stabilizes.
Typical useful range for transition-metal compounds at N = 6–10: T ≈
1000–3000 K. Converge your observable in T; the T→0 extrapolation recovers
the sharp-FS limit.

## Output check

The W-builder (hx0fp0 / combined hgw) prints at startup:

```
tetrakbt_init: T[K], kbt[Ry], kbt[eV]  2000  0.12667E-01  0.17234E+00
tetwt5: job efermi usetetrakbt:= 1 0.244900D+00 T
```

Verify the temperature and that `usetetrakbt` is `T`.
