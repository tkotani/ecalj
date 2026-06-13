# GW diagnostics: --dumpW, --WVR2ptRaxis, epsWVR

*(draft for ecaljdoc; target location: manual/gw section)*

## `--dumpW` (gwsc / hgw)

The combined `hgw` builds W per q-point in MPI shared memory and normally
leaves no W on disk (only the header `__WV.d`). With

```bash
gwsc 1 -np N --gpu --mp --fp32 --dumpW <sname>
```

each q-group root writes `__WVR.<iq>` (real axis) and `__WVI.<iq>`
(imaginary axis) after the Γ-cell correction (`W0w0i`), in the same
direct-access record format as the split `hx0fp0` flow:
`nblochpmx × nblochpmx` complex records, one per frequency, in the
product (E) basis; Wc = W − v. Use for offline analysis of Wc(q,ω).

## ε(q0,ω) head without extra cost: `epsWVR` lines

The W builder always prints, for the offset-Γ points,

```
epsWVR: iq iw_R omg(iw) eps(wFC) eps(woLFC)   <iq> <iw>  <omg>  <Re,Im eps_LFC>  <Re,Im eps_noLFC>
```

(ω in Ry). Grepping these lines gives the small-q dielectric head with and
without local fields; `−Im(1/ε)` locates plasmon poles (zeros of Re ε).

## `--WVR2ptRaxis` (gwsc / gw_lmfh)

The correlation pole term evaluates Wc at ω_ε = |ω−ε|/2 by 3-point Lagrange
interpolation in ω² on the real-axis histogram (`alagr3zz`). Lagrange end
weights can become negative, so a strongly curved Wc could in principle be
overshot. `--WVR2ptRaxis` switches both Σ implementations (m_sxcf_sc for
hgw; sxcf_fal2 for hsfp0) to 2-point linear interpolation on the bracketing
pair: convex weights in [0,1], no overshoot by construction.

Intended as an error-bound diagnostic: run with and without the flag; the
difference bounds the ω-interpolation error of the pole term. If the
difference is large, the histogram (`HistBin_dw`, `HistBin_ratio`) is too
coarse where Wc has sharp structure (e.g. plasmon poles); note that pointwise
sampling of a quasi-discrete pole does not converge by mesh refinement
alone — broaden it physically first (`tetrakbt`, or `SmearX0` of at least a
few bin widths).
