---
project: ecalj README at https://github.com/tkotani/ecalj
author: takao kotani
email: takaokotani@gmail.com
---
ecalj documents is at [ecaljdoc](https://ecalj.github.io/ecaljdoc/)

## 2026-06-13  Changelog: finite-T chi0, hsfp0 GPU, gw_lmfh GPU/MP flow, fixes

New in the GW chain (commits 37e6fbc2..2a04e767; user guide: FiniteT_and_QPE_HOWTO.md):

- **tetrakbt / t_tetrakbt** (`[gw]` in ctrlg toml): finite-temperature tetrahedron
  for chi0 with a consistent finite-T Fermi level (heftet writes `EFERMI_kbt`;
  m_tetwt consumes it). Physically broadens the Fermi surface — the recommended
  regularization for metallic QSGW instabilities (sharp nesting response).
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
