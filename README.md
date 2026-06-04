---
project: ecalj README at https://github.com/tkotani/ecalj
author: takao kotani
email: takaokotani@gmail.com
---
ecalj documents is at [ecaljdoc](https://ecalj.github.io/ecaljdoc/)

## 2026-05  Quick start (the new TOML flow)

Fortran binaries (lmf, lmfa, lmchk, gwsc, hsfp0, ...) read only:

  - `ctrlg.<sname>.toml`  -- merged ctrl + GW driver sections + PB cut-offs
  - `PB.<sname>.toml`     -- per-atom product-basis tables (GW path only)

### Starting from scratch (POSCAR or hand-written ctrls)

    # 1. prepare ctrls.<sname>  (basic structure: atoms, lattice, ...)
    # 2. generate the TOML pair:
    ctrlgenToml.py <sname>            # writes ctrlg.<sname>.toml + PB.<sname>.toml
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
    Legacy2toml.py <sname>            # writes ctrlg.<sname>.toml + PB.<sname>.toml
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
