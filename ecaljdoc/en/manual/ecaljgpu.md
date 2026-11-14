# How to use the GPU version of ecalj

> English version of the page [manual/ecaljgpu](../../manual/ecaljgpu.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/ecaljgpu.md of commit 2efea4c73 (source blob c5e85a8d7201).

<div class="clfix-old">

> ⚠️ The GPU version of `gwsc` also reads only `ctrlg.<sname>.toml`. Convert a legacy `ctrl.<sname>` / `GWinput` directory with `Legacy2toml.py <sname>` before running. See [Input file](./toml_migration).

</div>

ecalj can perform QSGW calculations on GPUs.
* The GPU is used for the calculation of the self-energy and the screened Coulomb interaction.
* At present `lmf` (the DFT calculation) runs on the CPU.
> [!IMPORTANT]
> You must specify the number of parallel processes separately for the part run on the CPU and the part run on the GPU.

* To use the GPU, the GPU version of ecalj must be compiled with the NVIDIA compiler.
<span class="clfix">Requirements: an NVIDIA GPU and NVHPC (with cuBLAS and cuSOLVER). Compile with `python3 InstallAll.py --fc nvfortran --gpu` (the CPU version is built as well). Running the GPU calculation on one node is recommended; on several nodes, make sure that the MPI processes are distributed over the nodes correctly. (Moved from the old English section at the end.)</span>
For installation, initial setup and GPU test calculations on Kugui, see [installISSP](../install/installISSP.md).

## Recommended settings on ISSP Kugui
* Use one ACC node (a node with GPUs). One node has 4 GPUs.
* Set the number of MPI processes for the part run on the CPU (given by `-np`) to 64.
* Set the number of MPI processes for the part run on the GPU (given by `-np2`) to 4, equal to the number of GPUs.

An example qsub file:

<div class="clfix-old">

```bash job.sh
#!/bin/sh
#PBS -q i1accs
#PBS -l select=1:ncpus=64:mpiprocs=64:ompthreads=1

ulimit -s unlimited
gwsc -np 64 -np2 4 --gpu 2 inas2gasb2 > lgwsc
```

</div>

<div class="clfix">

```bash job.sh
#!/bin/sh
#PBS -q i1accs
#PBS -l select=1:ncpus=64:mpiprocs=64:ompthreads=1

ulimit -s unlimited
gwsc 2 -np 64 -np2 4 --gpu inas2gasb2 > lgwsc
```

`--gpu` takes no argument; the number `2` is the number of iterations (the first argument of `gwsc`).

</div>

`i1accs` is a test queue (maximum 30 minutes); use `F1accs` or similar for production runs.

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

### Mixed precision: `--mp` and `--fp32`

`--mp` and `--mp --fp32` are old spellings of `--prec=tf32` and `--prec=fp32`, and can still be used with the same meaning
(the correspondence is in the table of the section "Choosing the precision and automatic selection of matrix operations" below).

With `tf32`, only the last product of Σc has its input mantissa reduced to 10 bits; the side that builds $W$ (the sum for $\chi_0$, the inverse of the dielectric matrix) is computed in FP32.
If the side that builds $W$ also uses a 10-bit mantissa, then **in systems of heavy elements + molecular anions (NO3 / N3 / ClO / BrO /
PO4 etc.), where local fields are strong and the condition number of ε(q,ω) is large**, an error of ~1e-3 is
amplified by $(1-v\chi_0)^{-1}$, $W$ / $\Sigma^c$ are destroyed, and QSGW diverges
(`lmf failed`) or gives NaN (`lgw` / `lqpe`). Examples and data to reproduce them are in
[`Samples/mptf32problem/`](https://github.com/tkotani/ecalj/tree/main/Samples/mptf32problem).
Storage is single precision for both `tf32` and `fp32`, so the memory usage is the same.

```bash
gwsc 5 -np 64 -np2 4 --gpu --prec=fp32 inas2gasb2 > lgwsc
```

### Choosing the precision and automatic selection of matrix operations (2026-09)

`gwsc --prec=tf32|fp32|fp64` selects the precision of the GW part. It is the same as the old `--mp`, `--mp --fp32`; either spelling works.

| `--prec` | Executable | Single-precision products | Old spelling |
| --- | --- | --- | --- |
| `fp64` | double-precision version (`hgw_gpu` etc.) | — | `--gpu` |
| `fp32` | mixed-precision version (`hgw_mp_gpu` etc.) | FP32 | `--gpu --mp --fp32` |
| `tf32` | mixed-precision version | only the Σc product has a 10-bit input mantissa (TF32 or FP16, chosen by the table); the rest FP32 | `--gpu --mp` |

`tf32` lowers the precision only of the last product of Σc (the product of W and zmel, zsec). Lowering the precision of the side that builds W (the sum for χ0, the inverse of the dielectric matrix)
lets small errors be magnified by $(1-v\chi_0)^{-1}$. For one `hgw` run of LiTi2O4 6³ (kt1, RTX 5090 ×2, version of the evening of 2026-09-27):

| | Time | Difference from double precision (Re Σc within E_F ±1 eV, max) |
| --- | --- | --- |
| fp32 | 367 s | 0.02 meV |
| tf32, Σc product in FP16 (chosen by the RTX 5090 table) | 173 s | 1.3 meV |
| tf32, Σc product in TF32 (`--linalg=tf32.cgemm.large=cublas,tf32.cgemm.small=cublas`) | — | 0.8 meV |
| all products in TF32 (`--mp` until 2026-09-27) | — | 5 meV (the NiO gap shifts by 0.1 eV) |

One QSGW iteration (`gwsc`, LiTi2O4 6³, tf32, 2 GPUs, 60 CPU cores) takes about 212 s, 80% of which is `hgw`.

Specifying only `--mp` is now also the same as `tf32`. Most of the tf32 error is a uniform shrinking of Σc (<del class="clfix">by a factor of −5×10⁻⁵ to −1×10⁻⁴ of Re Σc</del><span class="clfix">a relative change ΔΣ/Σ of Re Σc of −5×10⁻⁵ to −1×10⁻⁴</span>),
which comes from the way tensor cores accumulate sums (truncation of intermediate sums). Input rounding is round-to-nearest, and this cannot be reduced further cheaply.
The assumption is that tf32 is used for the intermediate iterations and the last iterations are run in fp32 with `--prec-final=fp32:N` (the error is not random but a uniform shrinking, so it should vanish in the final fp32 iterations.
**As of 2026-09-28 this usage has not actually been run and checked**. For 6³, even 10 iterations with tf32 only gave bands agreeing with fp32 to rms 0.3 meV).

Writing e.g. `--prec-final=fp32:2` runs only the last 2 iterations of that run at another precision (tf32 in between, fp32 at the end, etc.).

**Which method is used** on the GPU for matrix products and the inverse of the dielectric matrix is decided by a table for each precision. The user chooses only the precision.

- The table is measured on that GPU by `linalgtune_gpu` at the end of `InstallAll.py --gpu` and written to `<bindir>/ecalj_linalg_policy.toml`
  (when the GPU is being used by other jobs, the measurement is skipped; run `linalgtune_gpu` by hand later)
- When there is no table, or the GPU name differs, everything uses cuBLAS (the same behavior as before)
- The table used is printed once at the beginning of the standard output of `hgw` etc. (`linalg policy (fp32, from …)`)
- Candidate methods: for complex products `cublas`, `realsgemm` (recasts a complex product as one real SGEMM), `realhgemm` (its FP16 version; input in FP16,
  sums in FP32; the mantissa is 10 bits as in TF32, and on RTX 5090 it is 1.8 times faster than TF32; by precision it is chosen only in the `tf32` row), `gemmul8` (Ozaki scheme, INT8).
  For the inverse, `lu64` (LU in FP64), `mixed1` / `mixed2` (LU in FP32 + 1 / 2 Newton steps). The measuring tool selects, among those satisfying the error promise of that precision
  (FP32: relative 1e-5, FP64: 1e-12, etc.), only those at least 10% faster than cuBLAS
- To override by hand, list rows as in `--linalg=fp32.cgemm.large=realsgemm,fp32.epsinv.large=mixed1` (written as a `gwsc` option, it is passed to everything)
- The Σc product for `tf32` is decided by the `tf32` row of the table. `--linalg=tf32.cgemm.large=gemmul8:6` gives the Ozaki scheme with 6 splits,
  and raises the precision of the Σc product from 3e-4 of TF32 to 1e-5 (correspondingly slower). In the Ozaki scheme one split changes the precision by about one digit
- Example (RTX 5090 ×2, one `hgw` run of LiTi2O4): the table in fp32 is "single-precision complex products recast, double-precision products GEMMul8 14, inverse mixed1",
  and the `tf32` row (the Σc product) becomes the FP16 recast (`realhgemm`). For 6³, the former 833 s became 367 s with fp32 and 173 s with tf32 (evening of 2026-09-27);
  for 9³, 5004 s became 2990 s with fp32 and 1975 s with tf32 (version of the morning of the same day, with the tetrahedron-weight helper). <del class="clfix">The difference from the double-precision calculation in Re Σc also drops from 2.4e-4 to 4.2e-5 eV</del><span class="clfix">In one fp32 `hgw` run of 6³, the difference from the double-precision calculation (maximum of Re Σc over all states) also drops from 2.4e-4 eV with the version up to 09-26 (cuBLAS complex products) to 4.2e-5 eV with the recast products (`realsgemm`) (morning of 2026-09-27; within E_F ±1 eV it is the 0.02 meV of the table above)</span>
- GEMMul8 is not used for small products (any dimension below 64, or m·n·k < 1e8) regardless of the table (the fixed cost of the splitting exceeds the product)

### Running several jobs on the same machine

Before running an executable that uses the GPU, `gwsc` takes per-GPU locks (`/tmp/ecalj_res/gpu<N>.lock`, `flock`) for as many GPUs as needed,
and passes only the GPUs it obtained via `CUDA_VISIBLE_DEVICES`. If not enough are free, it waits (and shows what it is waiting for).

| Environment variable | Meaning |
| --- | --- |
| `CUDA_VISIBLE_DEVICES` | if set, GPUs are taken only from among these |
| `ECALJ_GPU_WAIT=<seconds>` | stop waiting after that many seconds and abort (default: wait forever) |
| `ECALJ_GPU_LOCK=0` | do not lock |

Locks are released when the process ends, so locks of crashed jobs do not remain.

### Screen output and logs (2026-09-27)

The GW programs (`hvccfp0`, `hsfp0_sc`, `hgw`, etc.) write the screen output of rank 0 directly into the `gwsc` logs (`lvcc`, `lsxC`, `lgw`, etc.).
Progress (`Time:` per q point, progress per (q, k)) and memory usage (`Memused (CPU) … (GPU) …`) can be followed there (`tail -f lgw`).
The screen output of other ranks is discarded. To see it per rank, `--fullstdo` (written for `gwsc`, it is passed to everything) produces `stdout.<rank>.<prog>`,
and `gwsc` moves them to `STDOUT/`. Errors of any rank go to standard error (the output of `gwsc`) with the rank number.
(Until 2026-09-27, all ranks (including rank 0) wrote `stdout.<rank>.<prog>`.)

### Computing the tetrahedron weights first on the CPU

While the GPU is running `hvccfp0` and `hgw`, `gwsc --gpu` runs the same `hgw` without GPU on the free CPU cores (the `-np` count) (`--tetwt_write`),
and writes the tetrahedron weights of each q to `__TETWT.<iq>.<isp>`. `hgw` reads them and skips the weight calculation (about 6 s per q for LiTi2O4 6³).
The result is bit-identical to computing them itself. When the files are missing or do not match, it computes them itself. To disable, use `--no-tetwt-helper`.

### MPS settings

Add the following to `~/.bashrc`.
```bash ~/.bashrc
if which nvidia-cuda-mps-control > /dev/null 2>&1 ; then
  export CUDA_MPS_PIPE_DIRECTORY=$(pwd)/nvidia-mps-$(hostname)
  export CUDA_MPS_LOG_DIRECTORY=$(pwd)/nvidia-log-$(hostname)
  echo "start nvidia-cuda-mps-control at" $(hostname)
  nvidia-cuda-mps-control -d
fi
```

<div class="clfix">

The old English section at the end ("Using the GPU Version") repeats the text above and is removed. What was not in the text above (requirements, the compile command, the advice to use one node) has been moved to the beginning. The note that GPU runs use 2 GB when `zmel_batch_gb` is not set (or is below 2) is left to the description in [gwinput](../manual/gwinput.md).

</div>

<div class="clfix-old">

--------------

**Using the GPU Version**

The GPU can be utilized in QSGW calculations:
- The GPU is used for the calculation of self-energy and screened Coulomb interaction.
- `lmf` is executed on the CPU.

**Requirements:**
1. NVIDIA GPU
2. Compiler: NVHPC with GPU-accelerated libraries: cuBLAS, cuSolver

**Compilation**
```bash
python3 InstallAll.py --fc nvfortran --gpu
```

**Execution**
- The GPU version can be used by specifying the `--gpu` option with `gwsc`.
- The number of MPI processes for the GPU code can be specified with `-np2`. This should typically match the number of GPUs available.
Example:
<div class="clfix-old">

```
gwsc -np 64 -np2 4 --gpu 1 $id > lgwsc
```

</div>

<div class="clfix">

```
gwsc 1 -np 64 -np2 4 --gpu $id > lgwsc
```

</div>

**Tips**

1. The GPU version uses fewer MPI processes compared to the CPU version, so the available memory per MPI process is larger
   This allows for larger batch sizes in the calculation, potentially improving computation speed.
   The variable controlling the batch size is `zmel_batch_gb` in the `[gw]` section of `ctrlg.<sname>.toml`.
   For GPU calculations, set the value to around 4 (representing 4GB).
   ```toml
   [gw]
   zmel_batch_gb = 4
   ```
   > [!TIP]
   > If not set (or set below 2), 2 GB is used in GPU calculations.

2. When handling large systems, add the following to `[gw]` to prevent memory exhaustion:
   ```toml
   [gw]
   KeepEigen = false
   ```

**Notes**
- It is recommended to run GPU calculations on a single node.
- When using multiple nodes, ensure that MPI processes are correctly distributed across nodes.

</div>
