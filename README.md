---
project: ecalj README at https://github.com/tkotani/ecalj
author: takao kotani
email: takaokotani@gmail.com
---
ecalj documents is at [ecaljdoc](https://ecalj.github.io/ecaljdoc/)

## 2026-05  Quick start (the new TOML flow)

Fortran binaries (lmf, lmfa, lmchk, gwsc, hsfp0, ...) read only:

  - `ctrlG.<sname>.toml`  -- merged ctrl + GW driver sections + PB cut-offs
  - `PB.toml`             -- per-atom product-basis tables (sname-free)

### Starting from scratch (POSCAR or hand-written ctrls)

    # 1. prepare ctrls.<sname>  (basic structure: atoms, lattice, ...)
    # 2. generate the TOML pair:
    ctrlgenToml.py <sname>            # writes ctrlG.<sname>.toml + PB.toml
    #    add --skipgw if you do not need GW (saves ~0.5 s)
    # 3. run as usual:
    lmfa <sname>
    lmf  <sname>
    gwsc 5 -np N <sname>              # GW (when needed)

`ctrlG.<sname>.toml` contains every ctrl/GWinput key with inline
comments (units, role, defaults).  Edit it directly; no re-conversion
step is required.

### Migrating an old (pre 2026-05) directory

    cd <your-old-dir>                 # has ctrl.<sname> and GWinput
    Legacy2toml.py <sname>            # writes ctrlG.<sname>.toml + PB.toml
    # legacy ctrl.<sname> / GWinput remain on disk but are no longer read.
    lmf <sname> ... # usual workflow

`Legacy2toml.py --help` documents every step.

### Old `-vfoo=bar` overrides

`%const` was removed.  Run-time overrides now use TOML-path syntax:

    OLD:  lmf si -vnk=8 -vmetal=3
    NEW:  lmf si -v[bz.nkabc]=[8,8,8] -v[bz.metal]=3

The `-v[<path>]=val` form is text-substituted into the TOML in memory
before parsing; the file on disk is never modified.

When `Legacy2toml.py` sees `-vNAME=VAL` it prints a 3-level diagnostic:

  - **WARN**  `NAME` is not in `%const`            -> the override is a no-op.
  - **INFO**  `NAME` maps to a TOML path           -> use `-v[<path>]=val`
                                                     at run time instead;
                                                     no reconversion needed.
  - **ERROR** `NAME` changes topology              -> save the result as a
                                                     variant file
                                                     `ctrlG.<sname>.<tag>.toml`
                                                     and switch via `cp`.
                                                     (Example:
                                                     `Samples/TestInstall/te`.)

### Tab completion

`InstallAll.py` appends a guarded source line to `~/.bashrc`
(skip with `--no-bashrc`) so new shells pick up `<sname>`
tab-completion for `Legacy2toml.py`, `ctrlgenToml.py`,
`lmf`, `lmfa`, `lmchk`, `gwsc`, `gw_lmfh`, `eps_lmfh`,
`epsPP_lmfh`, `genMLWF`, `genMLWFx`.  Candidates are taken
from `ctrl.*` / `ctrls.*` / `ctrlG.*.toml` in cwd, depending
on which input each script reads.
