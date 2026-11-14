# Input file {#toml-migration-2026-05}

> English version of the page [manual/toml_migration](../../manual/toml_migration.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/toml_migration.md of commit ca2563cd3 (source blob 3a11b235a9e7).

ecalj consists of the one-body part (the DFT part) and the GW part.

The present ecalj reads only the TOML-format file `ctrlg.<sname>.toml` (ctrlg below) as its basic input file; that is, every calculation starts from this ctrlg alone.
The old format (before 2026-05) needed `ctrl.<sname>` and `GWinput`. For backward compatibility the converter
`Legacy2toml.py`, which turns them into a ctrlg, is kept. To reproduce a past calculation, convert the input into a ctrlg
first ([below](#migration-of-an-existing-legacy-directory)).

## The basic input file is one TOML file {#the-one-toml-file}

| file | role |
|---|---|
| `ctrlg.<sname>.toml` | The basic input file of ecalj. Every calculation starts from this file. It consists of sections such as `[struc]`, `[[site]]`, `[[spec]]`, `[bz]`, `[iter]`, `[ham]`, `[gw]`, `[mlo]`, `[blocks]`, `[product_basis]` (and also `[options]`, `[str]`, `[ewald]`, `[esm]`, `[dyn]`, etc.) |

The basic input file is one ctrlg. The usual way is to generate ctrlg automatically from the crystal structure (`ctrlgenToml.py`) and edit it;
or options are given when it is generated.
The automatic generation converts a POSCAR (VASP format) with `vasp2ctrl` into the structure-only file `ctrls.<sname>`, and then makes the ctrlg
with `ctrlgenToml.py <sname>`. If you leave it to an AI, you may think of the POSCAR as the input.
The keys are described in detail in [lmf](./lmf) for the ctrl side (`[struc]`, `[[site]]`, `[[spec]]`, `[bz]`, `[ham]`, ...) and in [GW settings](./gwinput) for the GW side (`[gw]`, ...).

**The per-atom tables at the end**

- They are used in GW calculations; normally they are not edited by hand.
- They feed the mixed-product-basis generator (`hbasfp0` / `hvccfp0` etc.) used by `gwsc`, the eps tools and `mlo`'s W side.
- They are written by `gwinit` (through `ctrlgenToml.py` or `Legacy2toml.py`); leave them as generated.
- Until 2026-09 they were a separate `PB.<sname>.toml`. Such a file is not read any more:
  the binaries abort and point at `ctrlg_absorb.py <sname>`, which moves the tables into `ctrlg`.

**When the GW sections are not needed (`--skipgw`)**

- For a pure DFT workflow (just `lmf` / `lmfa` / band plots), pass `--skipgw` to `ctrlgenToml.py`.
- It skips the `lmfa → lmf --jobgw=0 → gwinit` sub-step and omits `[gw]` / `[mlo]` / `[blocks]` / `[product_basis]` from the output.

**Adding the GW sections later (`--addgw`, keeping your hand edits)**

- Use `ctrlgenToml.py <sname> --addgw`.
- It appends `[gw]` / `[mlo]` / `[blocks]` / `[product_basis]` without regenerating the ctrl-side keys;
  your edits to `[bz]`, `[ham]`, `[[spec]]` etc. are kept as they are.
- It refuses to run if `[gw]` is already present (to avoid silent duplication).
- Do **not** simply re-run `ctrlgenToml.py <sname>` for this purpose: the default flow overwrites the whole
  `ctrlg.<sname>.toml` from `ctrls.<sname>` (the previous file is moved to `ctrlg.<sname>.toml.bakup`, a one-generation backup only).

`<sname>` is the user-defined material extension (e.g. `si`, `cu`, `gaas`,
`fe`, ...). One working directory normally has exactly one
`ctrlg.<sname>.toml`; `lmf` auto-detects it without a positional argument
when run inside the directory.

<div class="clfix">

### Layout of the generated files {#file-layout}

`ctrlgenToml.py`, `gwinit`, `Legacy2toml.py` and the formatter `pylib/toml_tidy.py` all produce the same layout.
A blank line separates only the big headings (the `# === X ===` comment runs and single `[section]` headers).
Everything inside a section (one `[[site]]` or `[[spec]]` table from the next, a multi-line block from its neighbours, the three product-basis tables) is divided by a `# ----` rule.
A `[blocks]` section only appears for the few legacy `<...>` tags that have no better home (`QPNT`, `QforEPSL`, `hrotr`), kept as multi-line strings and parsed opaquely by the GW driver.
The sections come in the order `[gw]` `[mlo]` `[blocks]` `[product_basis]`; `[product_basis]` always closes the file (`toml_tidy.py` enforces the order).

</div>

## Migration of an existing legacy directory

```bash
# inside an old directory containing ctrl.<sname> [+ GWinput]
Legacy2toml.py <sname>
# produces ctrlg.<sname>.toml
# resume normal workflow (lmf, gwsc, job_eps, ...)
```

`Legacy2toml.py` lives at `ecalj/SRC/exec/Legacy2toml.py` (installed in
`~/bin/` along with the binaries by `InstallAll.py`). It is idempotent
and prints `[INFO] / [WARN] / [ERROR]` diagnostics for any
command-line `-v` overrides that would not survive the conversion.

<div class="clfix">

### Legacy keys and their TOML paths {#legacy-toml-map}

Most legacy `Category_Token` keys map 1:1 onto `[section.key]` after lower-casing. Those renamed or restructured are in Table 1.
The meaning of a key can be looked up in the [legacy ctrl reference](./lmf_legacy_ctrl).

Table 1. Legacy keys and TOML paths

| legacy `Category_Token` | TOML location | note |
|---|---|---|
| `STRUC_ALAT` / `STRUC_PLAT` | `[struc].alat` / `.plat` | |
| `SITE_ATOM=`, `POS=`, `XPOS=` | `[[site]] atom = ...; pos = [...]` (or `xpos`) | one `[[site]]` table per atom |
| `SPEC_ATOM=`, `Z=`, `R=`, `RSMH=`, `EH=`, ... | `[[spec]] atom = ...; z = ...; r = ...; rsmh = [...]` | one `[[spec]]` table per species; `R/W=` and `R/A=` kept as `"r/a"` quoted-key |
| `BZ_NKABC` / `BZ_METAL` / `BZ_TETRA` | `[bz].nkabc` / `.metal` / `.tetra` | |
| `HAM_NSPIN` / `HAM_SO` / `HAM_XCFUN` / `HAM_GMAX` / `HAM_PWMODE` / `HAM_PWEMAX` | `[ham].nspin` etc. | |
| `HAM_FORCES` / `HAM_OVEPS` / `HAM_FRZWF` / `HAM_REL` | `[ham].forces` etc. | |
| `HAM_ScaledSigma` | `[ham].scaledsigma` | lowercase |
| `HAM_READP` / `HAM_PNUFIX` | `[ham].readp` / `.pnufix` | `HAM_V0FIX` is the command-line option `--v0fix` |
| `esm_input.dat` (separate file) | **`[esm]`** section | retired 2026-09-16; not read — `ctrlg_absorb.py <sname>` converts it ([ESM](./esm)) |
| `ITER_NIT` / `ITER_MIX` / `ITER_CONV` / `ITER_CONVC` / `ITER_UMIX` / `ITER_TOLU` | `[iter].nit` / `.mix` / ... | |
| `SYMGRP` / `SYMGRPAF` | top-level `symgrp = "..."` (and `symgrpaf`) | not nested in a section |
| `EWALD_TOL` | `[ewald].tol` | rarely touched |
| `DYN_MODE` / `DYN_NIT` / `DYN_HESS` / ... | `[dyn].mode` etc. | |
| `IO_VERBOS` / `IO_TIM` | top-level `verbose` / `time` | promoted out of `[io]` (section dropped) |
| `<Worb>...</Worb>` block | `[mlo].mlo_lm = """ ... """` | multi-line string |
| `<QforEPS>...</QforEPS>` block | `[gw].QforEPS = """ ... """` (likewise `QforGW`) | |
| `n1n2n3` (GWinput) | `[gw].n1n2n3 = [k1,k2,k3]` | int vector |
| `HistBin_dw` / `HistBin_ratio` / `niw` (GWinput) | `[gw].HistBin_dw` etc. | unchanged names (upper and lower case as in GWinput) |
| `esmr`, `SmearX0` / `GaussianFilterX0`, `GaussSmear`, `delta`, `dw`, `omg_c`, `WgtQ0P` (GWinput) | `[gw].t_sigmaw` and `[gw].t_tetrakbt` (K; `t_tetrakbt` is required) | the old keys and what happens to them: [GW settings](./gwinput) Table 1, [kBT](./kBT) §2.5 Table 2 |
| product-basis cut-offs (`tolerance`, `lcutmx` from `<PRODUCT_BASIS>`) | `[product_basis].pb_tolerance` / `.pb_lcutmx` | |
| product-basis per-atom tables (`nlx`, `valence`, `core`) | `[product_basis].nlx` / `.valence` / `.core` (end of the file; a separate `PB.<sname>.toml` before 2026-09) | |

The full schema lives in [`SRC/exec/pylib/ctrl_schema.py`](https://github.com/tkotani/ecalj/blob/main/SRC/exec/pylib/ctrl_schema.py) (the conversion source of truth).
Annotations and inline comments are applied by [`SRC/exec/pylib/toml_comments.py`](https://github.com/tkotani/ecalj/blob/main/SRC/exec/pylib/toml_comments.py).

</div>

### `esm_input.dat` (slabs)

The separate positional `esm_input.dat` was retired on 2026-09-16 and became
the `[esm]` section of `ctrlg.<sname>.toml`. The Fortran does not read
the file: a leftover one makes `lmf` abort and name the converter,
`ctrlg_absorb.py <sname>`, which appends `[esm]` to the ctrlg if it is not
there yet and moves the original to `esm_input.dat.bk` with a header saying
where the settings went. `Legacy2toml.py` does the same during a legacy
conversion. The same script folds a leftover `PB.<sname>.toml` (the
per-atom product-basis tables, a separate file until 2026-09) into
`[product_basis]`; the GW binaries abort on that file likewise.
See [ESM](./esm).

### An older `ctrlg.<sname>.toml`: `ctrlg_update.py`

The programs stop on things that a `ctrlg.<sname>.toml` written by an
earlier version of the converters may carry:

- `[struc]` must not have `nbas` / `nspec` (the numbers of sites and
  species are those of the `[[site]]` / `[[spec]]` tables);
- `[gw]` must have `t_tetrakbt` (K), and must not have `SmearX0` /
  `SmearX0q0`;
- `[gw]` must have `t_sigmaw`, `HistBin_dw`, `HistBin_ratio`, `QforEPSau` and `EMAXforGW` (since 2026-11-14).

```bash
ctrlg_update.py ctrlg.<sname>.toml [...]    # in place; the original is kept as <file>.bak_update
```

removes `nbas` / `nspec` (if `nbas` was smaller than the number of
`[[site]]` tables, it comments out the tables beyond it), rewrites
`SmearX0 = s` as `t_tetrakbt = -T` of the same width, and adds
`t_tetrakbt = 0` where the key is missing, and adds the other required keys with the former code defaults ([GW settings](./gwinput)), so that the calculation stays
the same. `--no-backup` leaves no copy of the original. The rules are
described in [lmf.md](./lmf#file-structure-sections) and
[kBT](./kBT) §2; the old keys of `[gw]` are in Table 1 of
[GW settings](./gwinput).

## Run-time `--ctrlg:` overrides

<div class="clfix-old">

`%const` constants embedded in the legacy `ctrl.<sname>` are baked into
`ctrlg.<sname>.toml` at conversion time. Run-time overrides use
`--ctrlg:<dotted.path>=<value>`, processed in-memory by
`m_toml_override.f90` (no on-disk rewrite):

```bash
# OLD (the programs stop with a message that names the new form)
lmf si -vnk=8 -vmetal=3
lmf si -v[bz.nkabc]=[8,8,8] -v[bz.metal]=3
lmf si --toml.bz.nkabc=[8,8,8] --toml.bz.metal=3

# NEW
lmf si --ctrlg:bz.nkabc=[8,8,8] --ctrlg:bz.metal=3
```

The text after `--ctrlg:` is the dotted TOML path (`section.key`); the
value parses as TOML (`[8,8,8]` for an integer vector, `true`/`false`
lowercase for bools, quoted strings, etc.). The override layer
auto-quotes a bare identifier value (`find` → `"find"`) so that shell
quote-stripping is forgiving.

> ⚠️ **Retired override forms (all abort with a hint).**
> Earlier iterations of the migration accepted several override
> spellings (`-v<NAME>=<VAL>`, `-v[<path>]=<VAL>`, `--[<path>]=<VAL>`,
> `--<a.b.c>=<VAL>`, `--toml.<path>=<VAL>`, `--pr=N`, `--time=N,M`,
> `--phispinsym`). They are all gone; any of them on the cmdline
> triggers a one-line abort pointing at the canonical
> `--ctrlg:<path>=<value>` form.
>
> - `--ctrlg:bz.nkabc=[8,8,8]` — directly set the value at TOML path
>   `[bz].nkabc` in `ctrlg.<sname>.toml` for this run, in memory.
>   No `%const` indirection; no on-disk rewrite. <del class="clfix">The path must match
>   a real key in the schema.</del><span class="clfix">A key that the file does not carry is added (since 2026-09-30); a misspelled key is added too and read by no program, so check the log line ([cmdopts](./cmdopts)).</span> Each applied override is logged on
>   rank 0 so it appears in `llmf` etc.
>
> Practical consequences:
>
> - Old habits like `-vmetal=3`, `-vnspin=2 -vso=0` abort with a
>   migration hint. Use `--ctrlg:bz.metal=3`,
>   `--ctrlg:ham.nspin=2 --ctrlg:ham.so=0`.
> - `Legacy2toml.py` prints `[INFO] / [WARN] / [ERROR]` diagnostics
>   for any `%const` overrides it encounters during conversion to
>   help spot silent breakages.

</div>

<div class="clfix">

To change a key for one run, add `--ctrlg:<section.key>=<value>` (e.g. `lmf si --ctrlg:bz.nkabc=[8,8,8]`); the file is not rewritten. The details of the syntax, what happens to a key the file does not carry, and the retired forms (`-v<NAME>=<VAL>` and others, which make the programs stop) are in [Command-line options](./cmdopts#run-time-ctrlg-overrides).

</div>

## Samples {#migrated-samples-use-these-as-templates}

The list of sample `ctrlg.<sname>.toml` inputs to use as templates is on [Samples](./samples). The old examples of `Samples/Legacy/` were rebuilt as the directories listed there. Convert your own legacy input as described above.

## Result summaries (worked-out reference values)

- [MLOsamples/README.md § Reference value: bcc Fe Fe-3d on-site W](https://github.com/tkotani/ecalj/blob/main/Samples/MLOsamples/README.md#reference-value-bcc-fe-fe-3d-on-site-w-with---mlo_feb4) — bcc Fe Fe-3d on-site screened W ~1.5 eV (job_mloW reproduction, 2×2×2 vs 4×4×4 BZ)
- [MLOsamples/README.md § Regression test (testecalj Fe)](https://github.com/tkotani/ecalj/blob/main/Samples/MLOsamples/README.md#regression-test-testecalj-fe) — automated diagonal V/W-V check against inline reference values
- MLOsamples/README_SOC.md (`MD/mlo_notes/mlo_soc_memo.md`) — SOC-as-perturbation variant (`job_mlo_soc`)

## Auto-job runner (`ecalj_auto`)

For batch QSGW / GW jobs across many materials, see
[ecalj_auto/README.md](https://github.com/tkotani/ecalj/blob/main/ecalj_auto/README.md)
and [ecalj_auto/README_slot_scheduler.md](https://github.com/tkotani/ecalj/blob/main/ecalj_auto/README_slot_scheduler.md).
`MD/auto.md` says which of its drivers run with the present
programs; `jobtemplate.{kugui,ohtaka,ucgw,...}` are the templates for
SLURM/PBS clusters.

## Common gotchas

- <del class="clfix">**Retired override forms** — `lmf` aborts on `-vnspin=2 -vso=0`, `-v[ham.so]=1`, `--toml.ham.so=1`, `--pr=50`, `--time=5,5`, `--phispinsym`. Each abort prints the canonical `--ctrlg:<path>=<value>` replacement.</del>
- **Stale `Worb` blocks** — `gwinput2toml.py` (the GWinput leg of
  `Legacy2toml.py`) used to last-write-wins on duplicate `<Worb>` blocks;
  now keeps the first.
- **`mlo_emax`** — read only by the older `mlo_method = 0` and `1`; the
  default `mlo_method = 4` does not use it.
- **A required key such as `t_tetrakbt` missing, `SmearX0`, `nbas` / `nspec`** — the programs
  stop; `ctrlg_update.py` (above) brings the file up to the present rules.
- **gfortran 13.3 / 14.2 codegen** — a documented codegen bug was
  worked around in `m_HamPMT.f90` (one-line bait
  `aaa = trim(aaa) // ' '`); see commits `35ee9f225` and `f693848ce`
  on master.

For nvfortran 26.1 specifically: `findloc` runtime crashes on logical
arrays were fixed by adding `use m_nvfortran, only: findloc` to the
seven affected source files (commit `f693848ce`).
