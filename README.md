---
project: ecalj README at https://github.com/tkotani/ecalj
author: takao kotani
email: takaokotani@gmail.com
---
ecalj documents is at [ecaljdoc](https://ecalj.github.io/ecaljdoc/)

Developers, and AI sessions that start without memory: read ecaljdoc `manual/ForDevelopers.md` first (reading order,
build and test, machines, how jobs are submitted, a digest of the research log), and `ecaljclaude.md` (how code comments,
commits and the research log are written; Claude Code reads it through `CLAUDE.md`).

## 2026-09-28  `[gw] t_tetrakbt` is required; `t_tetrakbt < 0` replaces `SmearX0`

- `t_tetrakbt` (K) must be written in `[gw]`. Its sign selects how chi0 is smeared: `0` none, `T > 0` the finite-T
  tetrahedron, `-T` Im chi0 smeared along omega by a Gaussian as wide as the Fermi-Dirac distribution at T (what
  `SmearX0` did; `SmearX0` now stops the program). The gwinit template writes `t_tetrakbt = 300` and `t_sigmaw = 300`.
- `ctrlg_update.py ctrlg.<sname>.toml` (the former `ctrlg_drop_counts.py`) converts an older file without changing the
  calculation (`SmearX0 = 0.0057` Ha becomes `t_tetrakbt = -992.4`; a missing `t_tetrakbt` becomes 0). The Samples are converted.
- The formulas and the history of the keys: ecaljdoc `manual/kBT.md` §2; the key: `manual/gwinput.md`.

## 2026-09-28  ctrlg: no `nbas` / `nspec` in `[struc]` (lmf stops on them)

- The numbers of sites and species are the numbers of `[[site]]` / `[[spec]]` tables. A count written in `[struc]`
  used to win over the tables, so sites added to the tables were silently left out. Now `lmf` stops with a message.
  Fix an older file with `ctrlg_update.py ctrlg.<sname>.toml` (it removes the two lines, and if `nbas` was
  smaller than the tables it comments out the `[[site]]` blocks beyond it). To leave a site out, comment out its
  whole `[[site]]` block. `Legacy2toml.py` no longer writes them. ecaljdoc `manual/lmf.md`.
- `gwsc` keeps `efermi.lmf`, `QMLO_SigRs` and `QMLO_z` in `QSGW.<N>run/` as well, so that the bands of every
  iteration can be drawn afterwards (MLO-QSGW; 600 MB per iteration for LiTi2O4 on 9^3).

## 2026-09-28  Finite temperature (t_tetrakbt > 0): core exchange and the thermal quadrature

- The core exchange (`hsfp0_sc --job=3`) keeps its Fermi level below the valence. With `t_tetrakbt > 0` it was set
  to `EFERMI_kbt`, which let occupied valence states into SExcore; a second bug (fixed the same day) had hidden this
  except where the core states are fewer than the occupied bands. SExcore is now independent of T.
- The thermal kernel of the finite-T tetrahedron weights and of `EFERMI_kbt` is integrated on 4 panels
  [-6,-1,0,1,6] x 5-point Gauss-Legendre (20 points as before; one 20-point rule had no node within 0.92 kBT of E_F).
  bcc Fe at 3000 K: the error of `EFERMI_kbt` 9.4 -> 1.2 meV. Details: ecaljdoc `manual/kBT.md` 7.2 and 7.5.

## 2026-09-27  GPU GW: one precision switch; the method of each product is measured at installation

    gwsc 5 -np 32 -np2 2 --gpu --prec=fp32 <sname>    # or tf32 / fp64
    gwsc 10 -np 32 -np2 2 --gpu --prec=tf32 --prec-final=fp32:2 <sname>   # 8 iterations tf32, last 2 fp32

- `--prec`: fp64 = double precision; fp32 = mixed precision with FP32 products (= `--mp --fp32`);
  tf32 = fp32 but the products of Sigma_c with 10-bit-mantissa inputs (TF32, or FP16 where the table
  chooses it; on RTX 5090 the LiTi2O4 6^3 `hgw` takes 173 s against 367 s in fp32, and Re Sigma_c near E_F is within
  1.3 meV of double precision; numbers as of the night of 2026-09-27).
  `--mp` alone = tf32 (until 2026-09-27 it was TF32 in every product: 5 meV there, the NiO gap off by
  0.1 eV, for only 10% less time; W suffers through (1 - v chi0)^-1).
- Sigma_c runs asynchronously on the GPU (kernels and products on one stream; the pole weights are
  computed on the CPU meanwhile). LiTi2O4 6^3 `hgw` on 2 RTX 5090: fp32 424 -> 412 s, tf32 293 -> 242 s.

- **Which method runs each GPU matrix product and the epstilde inverse** (cuBLAS, `realsgemm` = a complex
  product as one real SGEMM, GEMMul8 = Ozaki on INT8, mixed-precision inverse) comes from one table,
  `<bindir>/ecalj_linalg_policy.toml`. `InstallAll.py --gpu` writes it at the end by measuring every method on
  this GPU (`linalgtune_gpu`, about a minute; skipped when the GPU is busy, `--notune` to skip). Without the
  file everything is cuBLAS, as before. `--linalg=fp32.cgemm.large=realsgemm,...` overrides rows.
- **GPU programs of two jobs on one machine do not collide**: `gwsc` takes a lock per GPU (`/tmp/ecalj_res/gpu<N>.lock`) for
  its GPU programs and waits while they are in use (`ECALJ_GPU_WAIT=<s>` to give up, `ECALJ_GPU_LOCK=0` to switch off).
  The lock is taken per step: the CPU steps (lmf with many ranks) of two jobs still compete, and two `hgw` on one node
  can run out of /dev/shm, so run one QSGW job per node.
- **Tetrahedron weights on the idle CPU cores**: with `--gpu`, `hgw --tetwt_write` runs next to `hgw` on the
  `-np` cores and hands the weights over in `__TETWT.*` files (bit-identical; `--no-tetwt-helper` to switch off).
- **The rest of a QSGW iteration** (evening and night): hvccfp0 on the GPU, lmf --jobgw=1 (H and H without xc in
  one pass), mlo --mlofreeze (Sigma^MLO only), hqpe_sc, the core exchange. One iteration of LiTi2O4 6^3 (tf32,
  2 RTX 5090, 60 CPU cores): 342 -> 212 s, 80% of it in `hgw` (LiTi2O4 6^3 `hgw`: fp32 367 s, tf32 173 s).
- **Console output of the GW programs**: rank 0 writes to the log of the step (`lgw`, `lsxC`, ...: progress and
  `Memused`), the other ranks write nothing; `--fullstdo` gives the per-rank `stdout.<rank>.<prog>` files as before.
  Errors of any rank reach stderr with the rank number.
- Report and numbers: `Samples/kBT/gpu_fp32_report.md`; user guide: ecaljdoc `manual/ecaljgpu.md`; for developers:
  ecaljdoc `manual/ForDevelopers.md` §11 and `ecaljclaude.md` (read by Claude Code through `CLAUDE.md`).

## 2026-06-13  Changelog: finite-T chi0, hsfp0 GPU, gw_lmfh GPU/MP flow, fixes

New in the GW chain (commits 37e6fbc2..2a04e767; user guide: FiniteT_and_QPE_HOWTO.md):

- **tetrakbt / t_tetrakbt** (`[gw]` in ctrlg toml): finite-temperature tetrahedron
  for chi0 with a consistent finite-T Fermi level (heftet writes `EFERMI_kbt`;
  m_tetwt consumes it). Physically broadens the Fermi surface — the recommended
  regularization for metallic QSGW instabilities (sharp nesting response).
  (Note 2026-09-28: the instability turned out to come from a sharp plasmon pole of W at the first-shell q
  hit by the real-axis pole term of Sigma_c, not from nesting; the remedy is `wcsmear` (default) and
  `t_sigmaw`, with `t_tetrakbt < 0` (formerly `SmearX0`) for W if needed. ecaljdoc manual/kBT.md 3.5.)
- **hsfp0_gpu** (new binary, gpu variant): CUDA-Fortran offload of the one-shot
  correlation W contractions. Validated against CPU at production scale
  (LiTi2O4 6^3: max |CPU-GPU| 3e-13 eV); ~13x per-rank speedup.
- **gw_lmfh**: TOML-mode gate (ctrlg+PB, same as gwsc), new options
  `--gpu --mp --fp32 -np2 N`. With `--gpu`, hvccfp0/hx0fp0 run as GPU variants
  and `--job=12` runs on hsfp0_gpu; with `--mp`, the single-precision WVR/WVI
  written by hx0fp0_mp* are promoted on read by the double-precision hsfp0.
  NOTE: if you change EMAXforGW/EMINforGW between stages, rerun hsfp0
  --job=3/--job=11 before hqpe (SEXU/SECU state counts must match).
- **Arbitrary-q QP energies**: `[blocks] QforGW` (3 reals per line, Cartesian
  2pi/alat) + the gw_lmfh flow evaluates diagonal Sigma at any q (e.g. band
  lines). Combine with EMINforGW/EMAXforGW (eV vs EF) to limit target bands
  (avoids 1d20-padded states at off-mesh q).
- **Diagnostics**: `--dumpW` (gwsc/hgw) persists the streaming-SHM W to
  `__WVR.<iq>/__WVI.<iq>` for offline Wc(q,omega) analysis;
  `--WVR2ptRaxis` switches the Sc real-axis pole interpolation to 2-point
  linear (overshoot-free) — bounds the omega-interpolation error.
- **Fixes**: real-axis pole binning OOB guard (findloc-miss wrote nttp(-1));
  hsfp0 imag-axis/zwz0 hand loops replaced by zgemm (~8x); gw_lmfh stale-W
  cleanup glob matched legacy names and never ran.
- **Branch `gwkbt-dev`**: finite-T GxW Stage A/B (gwkbt / gwkbt_boson keys) is
  isolated there as WIP — the Stage B entry2 static-bin fix on that branch
  needs design review and re-validation before merging.

## What's new (2026-06 .. 09)

See `HIGHLIGHTS_2026-06_09.md`: one input file `ctrlg.<sname>.toml`
(PB / esm_input.dat / GWinput.toml retired; `ctrlg_absorb.py` converts old
directories), MLO in its final form (`mlo_method = 4`), finite-T QSGW samples,
`gw_lmfh` on GPU, and the bug fixes.

## 2026-05  Quick start (the new TOML flow)

Fortran binaries (lmf, lmfa, lmchk, gwsc, hsfp0, ...) read one file only:

  - `ctrlg.<sname>.toml`  -- ctrl + GW driver sections ([gw] [mlo] [blocks]);
                             [product_basis] closes the file with the cut-offs
                             and the per-atom tables nlx / valence / core
                             (GW path only, written by gwinit, not hand-edited)

  Leftover side files from older layouts (`PB.<sname>.toml`, `esm_input.dat`)
  are not read: the binaries abort and name `ctrlg_absorb.py <sname>`, which
  folds them into ctrlg. `GWinput.toml` is dead and can be removed.

### Starting from scratch (POSCAR or hand-written ctrls)

    # 1. prepare ctrls.<sname>  (basic structure: atoms, lattice, ...)
    # 2. generate the TOML pair:
    ctrlgenToml.py <sname>            # writes ctrlg.<sname>.toml
    #    add --skipgw if you do not need GW (saves ~0.5 s)
    # 3. run as usual:
    lmfa <sname>
    lmf  <sname>
    gwsc 5 -np N <sname>              # GW (when needed)
    gwsc 5 -np N --gpu --mp --fp32 <sname>   # GPU mixed precision; --fp32 uses
                                             # true FP32 (not TF32) in the GEMMs.
                                             # Needed for ill-conditioned dielectrics
                                             # (heavy element + molecular anion, e.g.
                                             # NO3/N3/ClO): without it TF32 corrupts
                                             # W/SEc and QSGW diverges or yields NaN.
                                             # (2026-09-28: since 2026-09-27 --mp alone is
                                             # --prec=tf32, which keeps W in FP32.)

`ctrlg.<sname>.toml` contains every ctrl/GWinput key with inline
comments (units, role, defaults).  Edit it directly; no re-conversion
step is required.

### Migrating an old (pre 2026-05) directory

    cd <your-old-dir>                 # has ctrl.<sname> and GWinput
    Legacy2toml.py <sname>            # writes ctrlg.<sname>.toml
    # legacy ctrl.<sname> / GWinput remain on disk but are no longer read.
    lmf <sname> ... # usual workflow

`Legacy2toml.py --help` documents every step.

### Old `-vfoo=bar` overrides

`%const` was removed.  Run-time overrides now use TOML-path syntax:

    OLD:  lmf si -vnk=8 -vmetal=3
    NEW:  lmf si --ctrlg:bz.nkabc=[8,8,8] --ctrlg:bz.metal=3

The `--ctrlg:<path>=val` form is text-substituted into the TOML in memory
before parsing; the file on disk is never modified.

When `Legacy2toml.py` sees `-vNAME=VAL` it prints a 3-level diagnostic:

  - **WARN**  `NAME` is not in `%const`            -> the override is a no-op.
  - **INFO**  `NAME` maps to a TOML path           -> use `--ctrlg:<path>=val`
                                                     at run time instead;
                                                     no reconversion needed.
  - **ERROR** `NAME` changes topology              -> save the result as a
                                                     variant file
                                                     `ctrlg.<sname>.<tag>.toml`
                                                     and switch via `cp`.
                                                     (Example:
                                                     `Samples/TestInstall/te`.)

### Tab completion

`InstallAll.py` appends a guarded source line to `~/.bashrc`
(skip with `--no-bashrc`) so new shells pick up tab-completion
for the ecalj toolchain.

What completes:

  - `lmf <TAB>` / `lmfa <TAB>` / `lmchk <TAB>` -> the full
    `ctrlg.<sname>.toml` filenames in cwd.  The binary strips
    `ctrlg.` and `.toml` at startup, so `lmf ctrlg.nio.toml` and
    `lmf nio` are equivalent.  `lmf ctrlg.nio` (no `.toml` suffix)
    and `lmf ctrl.nio` (legacy text format) abort with a hint.
  - `lmf nio --<TAB>` -> every registered cmdopt0 / cmdopt2 flag
    (`--writeham`, `--jobgw=`, ...).  The flag list is dumped at
    install time by `lmf --listcmdopt`, so a registry edit is
    picked up on the next `./InstallAll.py` run.
  - `lmf nio --ctrlg:bz.<TAB>` -> every dotted-path key that
    actually exists in the cwd's `ctrlg.<sname>.toml`, parsed
    live by `python3 tomllib` on each TAB so an edit shows up
    immediately.  `[[spec]]` / `[[site]]` arrays expand to
    `spec.1.r`, `spec.2.r`, ...
  - `mpirun -np 8 lmf <TAB>` works the same -- the completion
    walks `mpirun`'s argument list, recognises the ecalj binary
    that follows, and delegates.
  - `Legacy2toml.py <TAB>` -> `ctrl.*`, `ctrlgenToml.py <TAB>`
    -> `ctrls.*`.

Because the completion list comes from the same source of truth
the binary uses at startup (the cmdopt registry + the on-disk
TOML), TAB-completed input never trips the strict typo /
key-not-found checks in `m_cmdopt_registry::validate_arglist`.
