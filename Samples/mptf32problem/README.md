# mptf32problem — TF32 breaks QSGW for ill-conditioned dielectrics (mp-8196, AgNO3)

This sample documents and reproduces a precision problem in the GPU
mixed-precision (`--mp`) GW path, and the `--fp32` fix.

## Symptom

Running QSGW on **mp-8196 (AgNO3, silver nitrate)** with the GPU
mixed-precision build (`gwsc --gpu --mp`) **diverges**: the static
self-energy blows up and `lmf` fails to converge ("lmf failed even after
reducing bmix to minimum"). The same calculation on CPU, or with
`--gpu --mp --fp32`, converges normally.

| case | flags | sigma_m (lqpe) | sev (eV) | result |
|------|-------|----------------|----------|--------|
| cpu  | (none)              | 529.45  | -204.67 | converged ✅ |
| tf32 | `--gpu --mp`        | **55488** | **-5059** | diverges ❌ |
| fp32 | `--gpu --mp --fp32` | 529.46  | -204.67 | converged ✅ |

`--fp32` matches the CPU reference to 4–5 digits; TF32 is off by ~100×.

## Cause

`--mp` runs the complex GEMMs (`cmm_d` in `m_blas.f90`) on the tensor cores
with **TF32** (`CUBLAS_COMPUTE_32F_FAST_TF32`, 10-bit mantissa, ~1e-3 relative
error). AgNO3's molecular NO3- anion produces a very strong local field:
the head dielectric constant is ~1536 while the local-field-corrected value is
~3.5, so the dielectric matrix `(1 - v·X0)` has condition number ~400.
The matrix inversion (done in FP64) faithfully **amplifies** TF32's ~1e-3
input error into ~40% error on `W = eps^-1·v`, which destroys the
correlation self-energy `SEc`. (The corruption seeds in the X0 accumulation,
upstream of the inversion — see the SEc-vs-SEx evidence below.)

The fix `--fp32` switches the GEMM compute type to true **FP32**
(`CUBLAS_COMPUTE_32F`, 23-bit mantissa, ~1e-7). Storage stays single
precision; only the matmul precision changes, so it is only ~7% slower than
TF32 (and ~3–4× faster than the FP64 path). FP64 would be overkill.

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
`reference/summary.txt` and `reference/QPU.*.head`.

Requires the GPU build (`hgw_mp_gpu`, `hsfp0_sc_mp_gpu`, ...) and
`ctrlgenToml.py`, `gwsc` on PATH.

## Practical guidance

Use `--fp32` for systems with strong local fields / ill-conditioned
dielectric matrices — heavy element + small molecular anion
(NO3, N3, ClO, BrO, IO, PO4, ...), oxides, d/f electrons. The default `--mp`
(TF32) is fine for homogeneous, well-conditioned systems (e.g. simple
semiconductors), where the per-element rounding averages out. When in doubt,
`--fp32`: it is only ~7% slower and always correct. The GW1500 batch drivers
(`ecalj_auto/`) pass `--fp32` by default for this reason.

See `Changes.txt` (2026-06-04) for the implementation.
