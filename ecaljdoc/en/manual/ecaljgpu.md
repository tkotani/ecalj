# How to use the ecalj GPU version

> English version of the page [manual/ecaljgpu](../../manual/ecaljgpu.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/ecaljgpu.md of commit 4f74d9132 (source blob 0d87eb4d5d30).

ecalj can perform QSGW calculations on GPUs.
* The GPU is used for calculating the self-energy and the screened Coulomb interaction.
* At present `lmf` (the DFT calculation) runs on the CPU.
> [!IMPORTANT]
> You need to specify the number of parallel processes separately for the part running on the CPU and the part running on the GPU.

* To use GPUs, you need to compile the GPU version of ecalj with the NVIDIA compiler.
Requirements: an NVIDIA GPU and NVHPC (including cuBLAS and cuSOLVER). Compile with `python3 InstallAll.py --fc nvfortran --gpu` (the CPU version is built at the same time). We recommend running GPU calculations on a single node. When using multiple nodes, check that the MPI processes are distributed to the nodes correctly.
For installation, initial setup, and the GPU test calculation on Kugui, see [installISSP](../install/installISSP.md).

## Recommended settings on ISSP Kugui
* Use one ACC node (a node with GPUs). Each node has 4 GPUs.
* Set the number of MPI processes for the part running on the CPU (given by `-np`) to 64.
* Set the number of MPI processes for the part running on the GPU (given by `-np2`) to 4, the same as the number of GPUs.

An example qsub file follows.

```bash job.sh
#!/bin/sh
#PBS -q i1accs
#PBS -l select=1:ncpus=64:mpiprocs=64:ompthreads=1

ulimit -s unlimited
gwsc 2 -np 64 -np2 4 --gpu inas2gasb2 > lgwsc
```

`--gpu` takes no argument. The first number `2` is the number of iterations (the first argument of `gwsc`).

`i1accs` is a test queue (maximum time 30 minutes); use `F1accs` or others for production runs.

Excerpt of the `[gw]` section of `ctrlg.<sname>.toml`
```toml
[gw]
t_tetrakbt   = 0            # (K) required: 0 = no smearing of chi0; -T = Im chi0 smeared as wide as Fermi-Dirac at T; T > 0 = finite-T tetrahedron (manual/kBT)
KeepEigen = false
zmel_batch_gb = 4
```

* `[gw].KeepEigen = false`
* `[gw].zmel_batch_gb = 4`

For the meaning of each, see [gwinput](../manual/gwinput.md).

### Mixed precision: `--mp` and `--fp32` {#gpu-mp-fp32}

`--mp` and `--mp --fp32` are the old spellings of `--prec=tf32` and `--prec=fp32`, and can still be used with the same meaning
(the correspondence is in the table of the section "Choosing the precision and automatic selection of the matrix operations" below).

With `tf32`, only the last product of Σc has its input mantissa reduced to 10 bits; the side that builds $W$ (the sum for $\chi_0$, the inverse of the dielectric matrix) is computed in FP32.
If the side that builds $W$ is also done with a 10-bit mantissa, then **in systems with heavy elements + molecular anions (NO3 / N3 / ClO / BrO /
PO4 etc.), where the local field is strong and the condition number of ε(q,ω) is large**, errors of ~1e-3 are
amplified by $(1-v\chi_0)^{-1}$, $W$ / $\Sigma^c$ are destroyed, and QSGW diverges
(`lmf failed`) or gives NaN (`lgw` / `lqpe`). Examples and data to reproduce this are in
[`Samples/mptf32problem/`](https://github.com/tkotani/ecalj/tree/main/Samples/mptf32problem).
The storage is single precision for both `tf32` and `fp32`, so the memory usage is the same.

```bash
gwsc 5 -np 64 -np2 4 --gpu --prec=fp32 inas2gasb2 > lgwsc
```

### Choosing the precision and automatic selection of the matrix operations (2026-09) {#gpu-precision}

`gwsc --prec=tf32|fp32|fp64` selects the precision of the GW part. These are the same as the old `--mp` and `--mp --fp32`; either spelling works.

| `--prec` | Executable | Single-precision products | Old spelling |
| --- | --- | --- | --- |
| `fp64` | double-precision version (`hgw_gpu` etc.) | — | `--gpu` |
| `fp32` | mixed-precision version (`hgw_mp_gpu` etc.) | FP32 | `--gpu --mp --fp32` |
| `tf32` | mixed-precision version | only the Σc products with a 10-bit input mantissa (TF32 or FP16, chosen by the table), others FP32 | `--gpu --mp` |

With `tf32`, the precision is lowered only for the last product of Σc (the product of W and zmel, zsec). If the precision of the side that builds W (the sum for χ0, the inverse of the dielectric matrix) is lowered,
small errors are magnified by $(1-v\chi_0)^{-1}$. For one `hgw` run of LiTi2O4 6³ (kt1, RTX 5090 ×2, the version of the night of 2026-09-27):

| | Time | Difference from double precision (Re Σc within E_F ±1 eV, maximum) |
| --- | --- | --- |
| fp32 | 367 s | 0.02 meV |
| tf32, Σc products in FP16 (chosen by the RTX 5090 table) | 173 s | 1.3 meV |
| tf32, Σc products in TF32 (`--linalg=tf32.cgemm.large=cublas,tf32.cgemm.small=cublas`) | — | 0.8 meV |
| all products in TF32 (`--mp` until 2026-09-27) | — | 5 meV (the NiO gap shifts by 0.1 eV) |

One QSGW iteration (`gwsc`, LiTi2O4 6³, tf32, 2 GPUs and 60 CPU cores) takes about 212 s, 80% of which is `hgw`.

`--mp` alone has also been made the same as `tf32`. Most of the tf32 error is a uniform shrinking of Σc (the relative change ΔΣ/Σ of Re Σc is −5×10⁻⁵ to −1×10⁻⁴),
which comes from the way the tensor cores accumulate sums (truncation of the intermediate sums). The input rounding is round-to-nearest; it cannot be reduced further cheaply.
It is assumed that tf32 is used for the intermediate iterations and that the last iterations are done in fp32 with `--prec-final=fp32:N` (since the error is not random but a uniform shrinking, it should vanish in the last fp32 iterations.
**As of 2026-09-28 this usage has not actually been run and checked**. For 6³, even 10 iterations of tf32 alone gave bands agreeing with fp32 within rms 0.3 meV).

Writing e.g. `--prec-final=fp32:2` runs only the last 2 iterations of that run in a different precision (tf32 in between and fp32 at the end, etc.).

**Which method is used** on the GPU for the matrix products and the inverse of the dielectric matrix is determined by a table for each precision. The user chooses only the precision.

- The table is measured on that GPU by `linalgtune_gpu` at the end of `InstallAll.py --gpu`, and written to `<bindir>/ecalj_linalg_policy.toml`
  (when the GPU is used by other jobs, the measurement is skipped; you can run `linalgtune_gpu` by hand later)
- When there is no table, or the GPU name differs, everything uses cuBLAS (the same behavior as before)
- The table used is printed once at the beginning of the standard output of `hgw` etc. (`linalg policy (fp32, from …)`)
- Candidate methods: for complex products, `cublas`, `realsgemm` (recasts the complex product into one real SGEMM), `realhgemm` (its FP16 version: input in FP16,
  sums in FP32. The mantissa is 10 bits, the same as TF32, and on the RTX 5090 it is 1.8 times faster than TF32. By its precision it is selected only in the `tf32` row), `gemmul8` (Ozaki scheme, INT8).
  For the inverse, `lu64` (FP64 LU), `mixed1` / `mixed2` (FP32 LU + 1 / 2 Newton steps). The measuring tool selects only methods that satisfy the error bound of that precision
  (FP32: relative 1e-5, FP64: 1e-12, etc.) and are at least 10% faster than cuBLAS
- To override by hand, list the rows like `--linalg=fp32.cgemm.large=realsgemm,fp32.epsinv.large=mixed1` (written as a `gwsc` option, it is passed to everything)
- The Σc products of `tf32` are determined by the `tf32` row of the table. `--linalg=tf32.cgemm.large=gemmul8:6` gives the Ozaki scheme with 6 splits,
  and the precision of the Σc products rises from 3e-4 of TF32 to 1e-5 (correspondingly slower). With Ozaki, one split changes the precision by about one digit
- Example (RTX 5090 ×2, one `hgw` run of LiTi2O4): the fp32 table is "single-precision complex products recast, double-precision products by GEMMul8 14, inverse by mixed1",
  and the tf32 row (the Σc products) becomes the FP16 recast (`realhgemm`). 6³ went from 833 s before to 367 s with fp32 and 173 s with tf32 (night of 2026-09-27);
  9³ went from 5004 s to 2990 s with fp32 and 1975 s with tf32 (the version of the morning of the same day, with the tetrahedron-weight helper). For one fp32 `hgw` run of 6³, the difference from the double-precision calculation (maximum of Re Σc over all states) also decreases from 2.4e-4 eV with the versions up to 09-26 (cuBLAS complex products) to 4.2e-5 eV with the recast (`realsgemm`) (morning of 2026-09-27; restricted to E_F ±1 eV it is the 0.02 meV of the table above)
- GEMMul8 is not used for small products (any side less than 64, or m·n·k less than 1e8) regardless of the table (the fixed cost of the splitting is larger than the product)

### Running several jobs on the same machine {#gpu-multijob}

Before running an executable that uses GPUs, `gwsc` takes per-GPU locks (`/tmp/ecalj_res/gpu<N>.lock`, `flock`) for as many GPUs as needed,
and passes only the GPUs it obtained via `CUDA_VISIBLE_DEVICES`. If not enough are free, it waits (and shows what it is waiting for).

| Environment variable | Meaning |
| --- | --- |
| `CUDA_VISIBLE_DEVICES` | if set, GPUs are taken only from among those GPUs |
| `ECALJ_GPU_WAIT=<seconds>` | stop waiting after that many seconds and abort (default: wait forever) |
| `ECALJ_GPU_LOCK=0` | do not lock |

A lock is released when the process ends, so locks of crashed jobs do not remain.

### Screen output and logs (2026-09-27)

The GW programs (`hvccfp0`, `hsfp0_sc`, `hgw`, etc.) write the screen output of rank 0 directly to the `gwsc` logs (`lvcc`, `lsxC`, `lgw`, etc.).
The progress (`Time:` for each q point, the progress for each (q, k)) and the memory usage (`Memused (CPU) … (GPU) …`) can be followed there (`tail -f lgw`).
The screen output of the other ranks is discarded. To see it per rank, `--fullstdo` (written for `gwsc`, it is passed to everything) produces `stdout.<rank>.<prog>`, and
`gwsc` moves them to `STDOUT/`. Errors of any rank are printed with the rank number to standard error (the output of `gwsc`).
(Until 2026-09-27, all ranks (including rank 0) wrote `stdout.<rank>.<prog>`.)

### Computing the tetrahedron weights first on the CPU {#gpu-tetwt-cpu}

While the GPU is running `hvccfp0` and `hgw`, `gwsc --gpu` runs the same `hgw` without GPU on the free CPU cores (the `-np` number) (`--tetwt_write`),
and writes the tetrahedron weights of each q to `__TETWT.<iq>.<isp>`. `hgw` reads them and skips the weight calculation (about 6 s per q for LiTi2O4 6³).
The results are bit-identical to those computed by itself. When the files are missing or do not match, it computes them itself. To disable this, use `--no-tetwt-helper`.

### MPS settings

Write the following in `~/.bashrc`.
```bash ~/.bashrc
if which nvidia-cuda-mps-control > /dev/null 2>&1 ; then
  export CUDA_MPS_PIPE_DIRECTORY=$(pwd)/nvidia-mps-$(hostname)
  export CUDA_MPS_LOG_DIRECTORY=$(pwd)/nvidia-log-$(hostname)
  echo "start nvidia-cuda-mps-control at" $(hostname)
  nvidia-cuda-mps-control -d
fi
```
