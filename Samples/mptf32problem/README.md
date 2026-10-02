# mptf32problem — TF32 breaks QSGW for ill-conditioned dielectrics (mp-8196, AgNO3)

> The TF32 of this page is the `--mp` of June 2026, which ran every GPU product in TF32, also those that build W.
> In the present `gwsc` (`--prec=tf32|fp32|fp64`) `--mp` alone means `--prec=tf32`: W (chi0, the inverse dielectric
> matrix) is built in FP32 and only the final products of Sigma_c take 10-bit inputs; `--mp --fp32` is `--prec=fp32`
> (`MD/kBT/gpu_fp32_report.md`, ecaljdoc [manual/ecaljgpu](https://ecalj.github.io/ecaljdoc/manual/ecaljgpu)).
> AgNO3 has not been rerun with the present tf32: the `tf32` numbers below, and `reference/`, are those of the old mode.

This sample documents a precision problem of the GPU mixed-precision GW path with TF32 in every
product (`--mp` of June 2026), and the `--fp32` fix.

## Symptom

Running QSGW on **mp-8196 (AgNO3, silver nitrate)** with the GPU
mixed-precision build and TF32 in every product (`gwsc --gpu --mp` of June 2026) **diverged**: the static
self-energy blew up and `lmf` failed to converge ("lmf failed even after
reducing bmix to minimum"). The same calculation on CPU, or with
`--gpu --mp --fp32`, converges normally.

| case | flags | sigma_m (lqpe) | sev (eV) | result |
|------|-------|----------------|----------|--------|
| cpu  | (none)              | 529.45  | -204.67 | converged ✅ |
| tf32 | `--gpu --mp`        | **55488** | **-5059** | diverges ❌ |
| fp32 | `--gpu --mp --fp32` | 529.46  | -204.67 | converged ✅ |

`--fp32` matches the CPU reference to 4–5 digits; TF32 is off by ~100×.

## Cause

That `--mp` ran all complex GEMMs (`cmm_d` in `m_blas.f90`) on the tensor cores
with **TF32** (`CUBLAS_COMPUTE_32F_FAST_TF32`, 10-bit mantissa, ~1e-3 relative
error). AgNO3's molecular NO3- anion produces a very strong local field:
the head dielectric constant is ~1536 while the local-field-corrected value is
~3.5, so the dielectric matrix `(1 - v·X0)` has condition number ~400.
The matrix inversion (done in FP64) faithfully **amplifies** TF32's ~1e-3
input error into ~40% error on `W = eps^-1·v`, which destroys the
correlation self-energy `SEc`. (The corruption seeds in the X0 accumulation,
upstream of the inversion — see the SEc-vs-SEx evidence below.)

The fix `--fp32` (now `--prec=fp32`) switches the GEMM compute type to true **FP32**
(`CUBLAS_COMPUTE_32F`, 23-bit mantissa, ~1e-7). Storage stays single
precision; only the matmul precision changes, so it was only ~7% slower than
that TF32 (and ~3–4× faster than the FP64 path). FP64 would be overkill.

## Evidence — `reference/QPU.{cpu,tf32,fp32}.head`

The quasiparticle table columns are: `q(3)  state  SEx  SExcore  SEc  vxc ...`.
Look at the correlation self-energy `SEc` of the deep valence states 4 and 5:

```
state   SEx        SEc(cpu)   SEc(tf32)
  4    -26.744     +6.725     -64.955     <- TF32 destroys SEc
  5    -23.271     +5.448     -93.980
```

`SEx` (exchange, bare Coulomb) is nearly the same; only `SEc` (from the
screened W) explodes under TF32. `--fp32` reproduces the CPU values exactly.

## Reproduce

```sh
./reproduce.sh [NP] [NP2]        # default NP=8 MPI ranks, NP2=2 GPU ranks
```

It regenerates the inputs from `ctrls.mp-8196` with the modern flow
(`ctrlgenToml.py --ssig=0.8 mp-8196`) and runs one QSGW iteration in each of
the three precision modes, printing sigma_m / sev / SEc. Compare against
`reference/summary.txt` and `reference/QPU.*.head`. Its `tf32` case calls `gwsc --gpu --mp`,
which is `--prec=tf32` with the present `gwsc`: it is not the mode that gave the `tf32` reference.

Requires the GPU build (`hgw_mp_gpu`, `hsfp0_sc_mp_gpu`, ...) and
`ctrlgenToml.py`, `gwsc` on PATH.

## Practical guidance

Use `--prec=fp32` (= `--mp --fp32`) for systems with strong local fields / ill-conditioned
dielectric matrices — heavy element + small molecular anion
(NO3, N3, ClO, BrO, IO, PO4, ...), oxides, d/f electrons. `--prec=tf32` (= `--mp` alone)
builds W in FP32 as well and uses 10-bit inputs only in the products of Sigma_c; it has not
been tested on AgNO3. When in doubt, `--prec=fp32`. The GW1500 batch drivers
(`ecalj_auto/`) pass `--mp --fp32` for this reason.

See `Changes.txt` (2026-06-04) for the implementation.
