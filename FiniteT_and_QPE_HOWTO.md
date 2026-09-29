# HOWTO: finite-T chi0 (t_tetrakbt) and one-shot QP energies on arbitrary k-lines (gw_lmfh + QforGW)

User guide for the features added 2026-06-12/13 (see Changes.txt for the changelog,
README for the one-paragraph summary). Sections 2 to 4 were exercised on
LiTi2O4 (metallic spinel) on kt1 (2x RTX5090); concrete numbers quoted come
from those runs.

## 1. Finite-temperature chi0: `t_tetrakbt`

In `ctrlg.<sname>.toml`, section `[gw]`, both in kelvin:

    t_tetrakbt = 2000    # required. 0: T=0 tetrahedron, no smearing of chi0;
                         #  T > 0: finite-temperature tetrahedron at T (heftet writes the Fermi level EFERMI_kbt);
                         # -T: T=0 tetrahedron, Im chi0 smeared by a Gaussian as wide as Fermi-Dirac at T
    t_sigmaw   = 1000    # Fermi-Dirac width of the levels in Sigma (default 1000; the gwinit template writes 300)

Older keys: `SmearX0` stops the programs (`ctrlg_update.py ctrlg.<sname>.toml` turns it into `t_tetrakbt = -T` and
adds a missing `t_tetrakbt`); `esmr` and `t_sigmakbt` are read only when `t_sigmaw` is absent; `tetrakbt = true` is ignored.
What the three cases compute, how to choose the temperatures, and what makes metallic QSGW oscillate (a sharp plasmon
pole of W; `wcsmear`, `t_sigmaw`, `t_tetrakbt < 0`):
ecaljdoc [manual/kBT](https://ecalj.github.io/ecaljdoc/manual/kBT), sections 0, 2 and 3.

Check in the output: with `t_tetrakbt > 0` the stdout of hx0fp0 / hgw has

    tetrakbt_init: T[K], kbt[Ry], kbt[eV]  2000  0.12667E-01 0.17234E+00

## 2. gw_lmfh: one-shot GW (G0W0-style QPE) — modern flow

`gw_lmfh` now uses the same TOML gate as gwsc (needs `ctrlg.<sname>.toml`
in cwd) and accepts GPU/MP options:

    gw_lmfh <sname> -np 60                          # CPU
    gw_lmfh <sname> -np 60 -np2 2 --gpu             # FP64 GPU W (hvccfp0/hx0fp0_gpu)
    gw_lmfh <sname> -np 60 -np2 2 --gpu --mp --fp32 # fast path: MP W (FP32 GEMM),
                                                    # SP WVR/WVI promoted on read

Notes:

  - `-np2` = ranks for the GPU steps (typically 1 rank per GPU).
  - `--mp` makes hvccfp0/hx0fp0 write single-precision WVR/WVI; the CPU hsfp0
    detects this (`hsfp0: single-precision WVR/WVI detected`) and promotes on
    read. Always combine with `--fp32` (true-FP32 GEMMs; plain TF32 corrupts
    ill-conditioned dielectric matrices — see the 2026-06-04 mp-8196 entry).
  - With `--gpu`, the correlation step `--job=12` runs on `hsfp0_gpu`
    (CUDA-Fortran offload, validated to 3e-13 eV against CPU at production
    scale; ~13x per-rank). Exchange steps stay on the CPU binary.
  - W files: `__WV.d`, `__WVR.<iq>`, `__WVI.<iq>` (the unprefixed `WV.d` is a
    legacy name no longer used anywhere).

## 3. QP energies on an arbitrary k-line (e.g. Gamma-K)

The QPU file from gwsc contains only the regular mesh q-points; band-path
features (kinks at band crossings) live BETWEEN mesh points and are created by
off-diagonal Sigma + rediagonalization, so they are invisible in mesh-QPU.
`gw_lmfh` + `QforGW` evaluates the DIAGONAL Sigma directly at any q, with no
interpolation.

### Recipe

1) Add the k-list and a target-band window to `ctrlg.<sname>.toml`:

    [gw]
    ...
    EMINforGW = -6.0   # eV relative to EF: lowest  Sigma target band
    EMAXforGW =  6.0   # eV relative to EF: highest Sigma target band
    QforGW = """
     0.000  0.000  0.000
     0.075 -0.075  0.000
     ...                      # 3 reals per line, Cartesian units of 2pi/alat
     0.750 -0.750  0.000      # (same convention as syml / QPU q-columns)
    """

2) Run (from a converged LDA `rst.<sname>`, or on top of a QSGW rst+sigm):

    gw_lmfh <sname> -np 60 -np2 2 --gpu --mp --fp32

3) Read `QPU`: rows with your q-vectors; eQP(noZ) = eLDA + dSEnoZ where
   dSEnoZ = SEx + SExcore + SEc - vxc (column 10 in the table part).

### Pitfalls (each one cost us a crashed or wasted run)

  - **1d20-padded states**: at off-mesh q the eigensolver stores fewer states
    than nband; the rest are padded with 1d20. Without EMAXforGW the Sigma
    targets include padded states and hsfp0 aborts with
    `sxcf 222: |w-e| out of range`. ALWAYS set EMAXforGW for QforGW runs.
  - **Window changes need an exchange rerun**: SEXU/XCU (from --job=3/--job=11)
    and SECU (--job=12) must have the same target-state count. If you edit
    EMIN/EMAXforGW after the exchange stage, rerun
    `hsfp0 --job=3` and `hsfp0 --job=11` before `hqpe`, otherwise hqpe dies
    with an EOF on SEXU.
  - **Memory (CPU correlation)**: each rank holds the full W(omega) stack in
    QPE mode (~ngb^2 * nw * 16B; 4.6 GB/rank for LiTi2O4 6^3). 60 ranks OOM'd
    a 256 GB node; 20 ranks fit. The GPU binary keeps this on the device.
  - **Diagonal only**: gw_lmfh gives Sigma_nn. Band kinks at (near-)crossings
    are driven by off-diagonal Sigma within the degenerate subspace — compare
    with the full QSGW band plot, or extract |Sigma_12| from the 2x2 inverse
    problem (E+ - E-)^2 = (d1-d2)^2 + 4|Sigma12|^2 using LDA splitting +
    diagonal-QPE splitting + QSGW band splitting (all alignment-free
    differences).

## 4. Diagnostics for W and the real-axis integration

  - `gwsc ... --dumpW` : the combined hgw normally keeps W in shared memory
    only. This writes `__WVR.<iq>/__WVI.<iq>` (split-flow record format) for
    offline analysis of Wc(q,omega).
  - eps head without extra cost: hgw/hx0fp0 stdout already contains lines
    `epsWVR: iq iw omg eps(wFC) eps(woLFC)` for the offset-Gamma q0 points —
    grep them to inspect the small-q dielectric function / loss function
    (plasmon poles at Re eps = 0).
  - `--WVR2ptRaxis` (gwsc and gw_lmfh): switch the Sc real-axis pole-term
    interpolation from 3-point Lagrange (end weights can be negative) to
    2-point linear (convex, overshoot-free). Run with and without to bound
    the omega-interpolation error; on LiTi2O4 6^3 the schemes agree to
    rms 6e-4 eV (interpolation NOT the source of the metallic instability).
  - Real-axis mesh: `HistBin_dw` / `HistBin_ratio` set the exponential
    omega-histogram. Beware: states whose we=|omega-e|/2 falls on sharp
    plasmon structure are sensitive to the local bin width; refine ratio
    (1.03 -> 1.015 -> ...) to test. Point-sampling of a quasi-discrete
    plasmon pole does NOT converge by refinement alone — give the pole a
    width first (`t_tetrakbt > 0`, or `t_tetrakbt < 0` with a Gaussian wider than the bins).

## 5. Branch status

  - `main`: everything above, all validated (testecalj si_gwsc / si_gw_lmfh /
    gas_pw_gw_lmfh PASS).
  - `gwkbt-dev`: finite-T Sigma (gwkbt Stage A/B) WIP. Known issue fixed on
    the branch (Stage B entry2 static-bin registration was dead); requires
    design review + re-validation against gwkbt_test_plasmonpole.py before
    any production use.
