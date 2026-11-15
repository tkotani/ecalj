# Input files {#toml-migration-2026-05}

> English version of the page [manual/toml_migration](../../manual/toml_migration.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/toml_migration.md of commit 296542c3a (source blob 90d0c155ab0e).

ecalj consists of the one-body calculation part (the DFT part) and the GW part.

The current ecalj reads only the TOML file `ctrlg.<sname>.toml` (hereafter ctrlg) as its basic input file. That is, every calculation starts from this ctrlg alone.
The old format (before 2026-05) required `ctrl.<sname>` and `GWinput`. For backward compatibility, the converter
`Legacy2toml.py`, which converts them to ctrlg, is kept. To reproduce past calculations, convert them beforehand to make ctrlg ([section below](#migration-of-an-existing-legacy-directory)).

## The basic input file is a single TOML file {#the-one-toml-file}

| File | Role |
|---|---|
| `ctrlg.<sname>.toml` | The basic input file of ecalj. Every calculation starts from this file. It consists of sections such as `[struc]`, `[[site]]`, `[[spec]]`, `[bz]`, `[iter]`, `[ham]`, `[gw]`, `[mlo]`, `[blocks]`, `[product_basis]` (and also `[options]`, `[str]`, `[ewald]`, `[esm]`, `[dyn]`, etc.) |

The basic input file is the single ctrlg. The basic way is to generate ctrlg automatically from the crystal structure (`ctrlgenToml.py`) and edit it for use.
Alternatively, generate it with options at generation time.
For automatic generation, convert a POSCAR (VASP format) into the structure-only file `ctrls.<sname>` with `vasp2ctrl`, and then make ctrlg with `ctrlgenToml.py <sname>`.
If you leave it to an AI, you can regard POSCAR as the input.
Detailed descriptions of the keys are in [lmf](./lmf) for the ctrl side (`[struc]`, `[[site]]`, `[[spec]]`, `[bz]`, `[ham]`, etc.) and in [GW settings](./gwinput) for the GW side (`[gw]` etc.).

**The per-atom tables at the end**

- Used in GW calculations. Usually not edited by hand.
- They are passed to the programs that build the mixed product basis (`hbasfp0` / `hvccfp0` etc.) used by `gwsc`, the eps tools, and the W side of `mlo`.
- Written by `gwinit` (through `ctrlgenToml.py` or `Legacy2toml.py`). Leave them as generated.
- Until 2026-09 they were in a separate file `PB.<sname>.toml`. This file is no longer read.
  If it remains, the binaries stop and tell you to use `ctrlg_absorb.py <sname>`. This script moves the tables into `ctrlg`.

**When the GW sections are not needed (`--skipgw`)**

- For DFT-only use (`lmf` / `lmfa` / band plots only), pass `--skipgw` to `ctrlgenToml.py`.
- It skips the steps `lmfa → lmf --jobgw=0 → gwinit` and does not write `[gw]` / `[mlo]` / `[blocks]` / `[product_basis]` to the output.

**Adding the GW sections later (`--addgw`, keeping hand edits)**

- Use `ctrlgenToml.py <sname> --addgw`.
- It appends `[gw]` / `[mlo]` / `[blocks]` / `[product_basis]` without regenerating the keys of the ctrl side.
  Changes you made in `[bz]`, `[ham]`, `[[spec]]`, etc. are kept as they are.
- It does not run when `[gw]` already exists (to prevent unnoticed duplication).
- You must **not** simply rerun `ctrlgenToml.py <sname>` for this purpose. The usual flow rewrites the whole
  `ctrlg.<sname>.toml` from `ctrls.<sname>` (the previous file is moved to `ctrlg.<sname>.toml.bakup`, but only one generation of backup is kept).

`<sname>` is the material extension chosen by the user (`si`, `cu`, `gaas`, `fe`, etc.). A working directory usually has
exactly one `ctrlg.<sname>.toml`. When run inside that directory, `lmf` finds that file even without arguments.

### Layout of the generated files {#file-layout}

`ctrlgenToml.py`, `gwinit`, `Legacy2toml.py`, and the formatting tool `pylib/toml_tidy.py` all produce files with the same layout.
Only the major headings (the run of `# === X ===` comments and a lone `[section]` heading) are separated by blank lines.
Within a section (between `[[site]]` or `[[spec]]` tables, a multi-line block and what precedes and follows it, the three product-basis tables), `# ----` rules are used as separators.
The `[blocks]` section appears only when there are old `<...>` tags with no other place to go (`QPNT`, `QforEPSL`, `hrotr`); it holds them as multi-line strings, which the GW programs read as they are.
The order of the sections is `[gw]` `[mlo]` `[blocks]` `[product_basis]`, and `[product_basis]` is always last (`toml_tidy.py` arranges them in this order).

## Migration of an old-format directory {#migration-of-an-existing-legacy-directory}

```bash
# inside an old directory containing ctrl.<sname> [+ GWinput]
Legacy2toml.py <sname>
# produces ctrlg.<sname>.toml
# resume normal workflow (lmf, gwsc, job_eps, ...)
```

`Legacy2toml.py` is at `ecalj/SRC/exec/Legacy2toml.py` (`InstallAll.py` puts it into `~/bin/` together with the binaries).
Running it any number of times gives the same result. If there are command-line `-v` overrides that do not survive the conversion, it prints `[INFO] / [WARN] / [ERROR]` diagnostics.

### Correspondence between old keys and TOML paths {#legacy-toml-map}

Many of the old `Category_Token` keys correspond directly to `[section.key]` when lowercased. Those whose name or form changed are as follows (Table 1).
The meanings of the keys can be looked up in [the old ctrl document](./lmf_legacy_ctrl).

Table 1. Old keys and TOML paths

| Old `Category_Token` | Location in TOML | Note |
|---|---|---|
| `STRUC_ALAT` / `STRUC_PLAT` | `[struc].alat` / `.plat` | |
| `SITE_ATOM=`, `POS=`, `XPOS=` | `[[site]] atom = ...; pos = [...]` (or `xpos`) | one `[[site]]` table per atom |
| `SPEC_ATOM=`, `Z=`, `R=`, `RSMH=`, `EH=`, ... | `[[spec]] atom = ...; z = ...; r = ...; rsmh = [...]` | one `[[spec]]` table per species. `R/W=` and `R/A=` are the quoted key `"r/a"` |
| `BZ_NKABC` / `BZ_METAL` / `BZ_TETRA` | `[bz].nkabc` / `.metal` / `.tetra` | |
| `HAM_NSPIN` / `HAM_SO` / `HAM_XCFUN` / `HAM_GMAX` / `HAM_PWMODE` / `HAM_PWEMAX` | `[ham].nspin` etc. | |
| `HAM_FORCES` / `HAM_OVEPS` / `HAM_FRZWF` / `HAM_REL` | `[ham].forces` etc. | |
| `HAM_ScaledSigma` | `[ham].scaledsigma` | lowercase |
| `HAM_READP` / `HAM_PNUFIX` | `[ham].readp` / `.pnufix` | `HAM_V0FIX` is the command-line option `--v0fix` |
| `esm_input.dat` (a separate file) | the **`[esm]`** section | dropped on 2026-09-16. It is not read, so convert it with `ctrlg_absorb.py <sname>` ([ESM](./esm)) |
| `ITER_NIT` / `ITER_MIX` / `ITER_CONV` / `ITER_CONVC` / `ITER_UMIX` / `ITER_TOLU` | `[iter].nit` / `.mix` / ... | |
| `SYMGRP` / `SYMGRPAF` | top-level `symgrp = "..."` (and `symgrpaf`) | not put inside a section |
| `EWALD_TOL` | `[ewald].tol` | rarely touched |
| `DYN_MODE` / `DYN_NIT` / `DYN_HESS` / ... | `[dyn].mode` etc. | |
| `IO_VERBOS` / `IO_TIM` | top-level `verbose` / `time` | moved out of `[io]` (the section was removed) |
| `<Worb>...</Worb>` block | `[mlo].mlo_lm = """ ... """` | multi-line string |
| `<QforEPS>...</QforEPS>` block | `[gw].QforEPS = """ ... """` (same for `QforGW`) | |
| `n1n2n3` (GWinput) | `[gw].n1n2n3 = [k1,k2,k3]` | integer vector |
| `HistBin_dw` / `HistBin_ratio` / `niw` (GWinput) | `[gw].HistBin_dw` etc. | names are the same as in GWinput (including case) |
| `esmr`, `SmearX0` / `GaussianFilterX0`, `GaussSmear`, `delta`, `dw`, `omg_c`, `WgtQ0P` (GWinput) | `[gw].t_sigmaw` and `[gw].t_tetrakbt` (K; `t_tetrakbt` is required) | old keys and how they are handled: Table 1 of [GW settings](./gwinput), Table 2 in §3.5 of [kBT](./kBT) |
| upper limits of the product basis (`tolerance`, `lcutmx` of `<PRODUCT_BASIS>`) | `[product_basis].pb_tolerance` / `.pb_lcutmx` | |
| per-atom tables of the product basis (`nlx`, `valence`, `core`) | `[product_basis].nlx` / `.valence` / `.core` (at the end of the file; before 2026-09, the separate `PB.<sname>.toml`) | |

The full schema is in [`SRC/exec/pylib/ctrl_schema.py`](https://github.com/tkotani/ecalj/blob/main/SRC/exec/pylib/ctrl_schema.py) (the reference for the conversion).
Annotations and end-of-line comments are added by [`SRC/exec/pylib/toml_comments.py`](https://github.com/tkotani/ecalj/blob/main/SRC/exec/pylib/toml_comments.py).

### `esm_input.dat` (slabs) {#esm-input-dat-slabs}

The separate file `esm_input.dat` was dropped on 2026-09-16 and became the `[esm]` section of `ctrlg.<sname>.toml`. Fortran does not read this file.
If it remains, `lmf` stops and names the conversion script `ctrlg_absorb.py <sname>`. This script adds `[esm]` to ctrlg if it is not there yet,
moves the original file to `esm_input.dat.bk`, and writes at its head where the settings were moved. `Legacy2toml.py` does the same when converting the old format.
The same script also takes a remaining `PB.<sname>.toml` (the per-atom product-basis tables; a separate file until 2026-09) into `[product_basis]`.
The GW binaries stop in the same way if this file exists. See [ESM](./esm).

### An older `ctrlg.<sname>.toml`: `ctrlg_update.py` {#an-older-ctrlg-sname-toml-ctrlg-update-py}

A `ctrlg.<sname>.toml` written by a previous version of the conversion scripts may contain things that make the programs stop.

- `[struc]` must not have `nbas` / `nspec` (the numbers of sites and species are determined by the numbers of `[[site]]` / `[[spec]]` tables).
- `[gw]` needs `t_tetrakbt` (K) and must not have `SmearX0` / `SmearX0q0`.
- `[gw]` needs `t_sigmaw`, `HistBin_dw`, `HistBin_ratio`, `QforEPSau`, `EMAXforGW` (since 2026-11-14).

```bash
ctrlg_update.py ctrlg.<sname>.toml [...]    # in place; the original is kept as <file>.bak_update
```

This removes `nbas` / `nspec` (if `nbas` was smaller than the number of `[[site]]` tables, the later tables are commented out),
rewrites `SmearX0 = s` as `t_tetrakbt = -T` of the same width, adds `t_tetrakbt = 0` where the key is missing, and also adds the other required keys with the defaults of the previous code, so that the calculation does not change ([GW settings](./gwinput)).
With `--no-backup`, no backup of the original file is kept. The rules are in [lmf.md](./lmf#file-structure-sections) and §3 of [kBT](./kBT), and
the old keys of `[gw]` are in Table 1 of [GW settings](./gwinput).

## Run-time overrides with `--ctrlg:` {#run-time-ctrlg-overrides}

To change a key for a single run only, add `--ctrlg:<section.key>=<value>` (e.g. `lmf si --ctrlg:bz.nkabc=[8,8,8]`). The file is not rewritten. Details of the syntax, the handling of keys not in the file, and the dropped old syntax (`-v<NAME>=<VAL>` etc.; using it stops the program) are in [Command-line options](./cmdopts#run-time-ctrlg-overrides).

## Samples {#migrated-samples-use-these-as-templates}

The list of sample `ctrlg.<sname>.toml` files to use as templates is in [Samples](./samples). The examples of the old `Samples/Legacy/` have also been rebuilt as the directories listed there. Convert your own old-format input as described above.

## Summaries of results (worked-out reference values) {#result-summaries-worked-out-reference-values}

- [MLOsamples/README.md § Reference value: bcc Fe Fe-3d on-site W](https://github.com/tkotani/ecalj/blob/main/Samples/MLOsamples/README.md#reference-value-bcc-fe-fe-3d-on-site-w-with---mlo_feb4) — the on-site screened W of Fe-3d in bcc Fe is about 1.5 eV (reproduced with job_mloW, BZ 2×2×2 and 4×4×4)
- [MLOsamples/README.md § Regression test (testecalj Fe)](https://github.com/tkotani/ecalj/blob/main/Samples/MLOsamples/README.md#regression-test-testecalj-fe) — a test that automatically compares the diagonal V/W-V with the reference values written in it
- MLOsamples/README_SOC.md (`MD/mlo_notes/mlo_soc_memo.md`) — the version that treats SOC as a perturbation (`job_mlo_soc`)

## Tools for running jobs automatically (`ecalj_auto`) {#auto-job-runner-ecalj-auto}

To run QSGW / GW jobs for many materials in bulk, see
[ecalj_auto/README.md](https://github.com/tkotani/ecalj/blob/main/ecalj_auto/README.md) and
[ecalj_auto/README_slot_scheduler.md](https://github.com/tkotani/ecalj/blob/main/ecalj_auto/README_slot_scheduler.md).
Which driver scripts work with the current programs is written in `MD/auto.md`. `jobtemplate.{kugui,ohtaka,ucgw,...}` are templates for SLURM/PBS clusters.

## Common pitfalls {#common-gotchas}

- **Old `Worb` blocks** — `gwinput2toml.py` (the part of `Legacy2toml.py` that handles GWinput) used to take the last one when `<Worb>` blocks
  were duplicated. It now takes the first one.
- **`mlo_emax`** — read only by the old `mlo_method = 0` and `1`. The default `mlo_method = 4` does not use it.
- **Missing required keys such as `t_tetrakbt`, `SmearX0`, `nbas` / `nspec`** — the programs stop. Bring the file to the current rules with `ctrlg_update.py` (above).
- **Code generation of gfortran 13.3 / 14.2** — a known code-generation bug is worked around in `m_HamPMT.f90` (a one-line bait
  `aaa = trim(aaa) // ' '`). See master commits `35ee9f225` and `f693848ce`.

Specific to nvfortran 26.1: the run-time crash of `findloc` on logical arrays was fixed by adding
`use m_nvfortran, only: findloc` to the 7 affected source files (commit `f693848ce`).
