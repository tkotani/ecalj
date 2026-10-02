# Command-line options for `lmf` / `lmfa` / `lmchk` (and friends)

## Cmdline argument structure (the minimum)

The Fortran binaries accept four token forms. Anything else aborts at
startup with a one-line hint. All parsing lives in
[`SRC/subroutines/m_cmdopt_registry.f90`](https://github.com/tkotani/ecalj/blob/main/SRC/subroutines/m_cmdopt_registry.f90)
(option table + typo / retired-syntax check) and
[`SRC/subroutines/m_toml_override.f90`](https://github.com/tkotani/ecalj/blob/main/SRC/subroutines/m_toml_override.f90)
(TOML override application).

| Form | Meaning | Where it goes |
|---|---|---|
| `<sname>` | positional, the ctrlg file extension | `m_ext::sname` |
| `ctrlg.<sname>.toml` | positional, the full TOML filename (e.g. after TAB completion) | `m_ext::sname` |
| `--ctrlg:<dotted.path>=<value>` | runtime TOML override (applied in memory; file untouched) | `m_toml_override` |
| `--<flag>` | registered cmdopt0 — presence bool | `c0_<flag>` in `m_cmdopt_registry` |
| `--<flag>=<value>` | registered cmdopt2 — typed value | `c2_<flag>` in `m_cmdopt_registry` |

The two positional forms are equivalent: `lmf nio` and
`lmf ctrlg.nio.toml` both end up with `sname='nio'`. m_ext_init
rejects anything in between -- `ctrlg.nio` (no `.toml` suffix) and
`ctrl.nio` (the legacy text format) abort at startup with a hint
pointing at the two canonical forms (and at `Legacy2toml.py` for
the latter).

A single argument-validation pass runs at startup
(`m_cmdopt_registry::validate_arglist`, called from `m_args::m_setargs`).
It aborts on:
- any `--<word>` / `-<word>` that is not in the cmdopt0 table;
- any `--<word>=<value>` whose bare prefix is not in the cmdopt2 table;
- any of the retired forms below (with a per-form migration hint).

Adding a new option is one line each in `m_cmdopt_registry.f90`:
declare `c0_<name>` (or `c2_<name>`) and add a matching
`call set0('--<name>', c0_<name>, narg, arglist)`
(resp. `get2('--<name>', outs, narg, arglist)`) in the relevant
`load_*_registry` subroutine. The typo-detection table is built as a
side effect — there is no parallel "known-flag list" to keep in sync.

This page is the catalogue of registered options. It is written by hand
and checked against
`grep -E "call set0\('--|get2\('--" SRC/subroutines/m_cmdopt_registry.f90`
(or `lmf --listcmdopt`).

## How to read the table

- **Flags** are cmdopt0 (no value): write them verbatim, e.g. `--writeham`.
- **Key=value** are cmdopt2: always include `=`, e.g. `--jobgw=1`.
  Unrecognised `--<word>=<value>` aborts with a typo hint.
- **TOML override** lives in a different namespace: `--ctrlg:<dotted.path>=<value>`
  (handled in `m_toml_override.f90`). It does not appear in this table.
- Most flags work in `lmf` (and `lmf --jobgw=…` which acts as the
  GW driver). `lmfa` / `lmchk` honour the few that apply to the
  atomic-program / structure-check phases — they are explicitly
  marked below.

## Retired (now abort with migration hint)

These used to be `cmdopt0` / `cmdopt2` shortcuts that overrode a TOML
key directly. They are gone and any of these on the cmdline triggers
an abort + one-line hint pointing at the canonical
`--ctrlg:<path>=<value>` form.

| Retired form | Replacement |
|---|---|
| `--pr=N` / `--pr N` / `-pr N` | `--ctrlg:verbose=N` |
| `--time=N[,M]` | `--ctrlg:time=[N,M]` |
| `--phispinsym` (cmdline flag) | `--ctrlg:ham.phispinsym=true` |
| `-v<name>=<value>` (legacy `%const`) | `--ctrlg:<dotted.path>=<value>` |
| `-v[<path>]=<value>` (older bracketed) | `--ctrlg:<dotted.path>=<value>` |
| `--[<path>]=<value>` (intermediate) | `--ctrlg:<dotted.path>=<value>` |
| `--<a.b.c>=<value>` (intermediate dotted) | `--ctrlg:<a.b.c>=<value>` |
| `--toml.<path>=<value>` (intermediate) | `--ctrlg:<path>=<value>` |
| `-emin=` / `-emax=` / `-ndos=` (single dash, legacy) | `--emin=` / `--emax=` / `--ndos=` |
| `-EfermiShifteV=` (single dash, legacy) | `--EfermiShifteV=` |
| `--terse` / `-terse` (lmchk overlap-table) | dropped — use `--ctrlg:verbose=N` (N≤40) |
| `--no-iactiv` (legacy lmfa flag) | dropped — silently ignored before, now aborts |

## Common (lmf / lmfa / lmchk)

| Option | Type | Effect | Used in |
|---|---|---|---|
| `--help` | flag | Print Usage banner + link to this page, then exit | `m_ext.f90` |
| `--debug` | flag | Global debug printout | many |
| `--show_time` | flag | Print `hambl` timing | `hambl.f90` |
| `--fullstdo` | flag | Every rank prints (not only rank 0). In the GW programs each rank then writes its own `stdout.<rank>.<prog>` (gwsc moves them to `STDOUT/`); without it rank 0 writes to the log of the run and the other ranks to /dev/null (2026-09-27) | `m_mpi.f90` |

## SCF / electronic structure (lmf core)

| Option | Type | Effect | Used in |
|---|---|---|---|
| `--etot` | flag | Total-energy diagnostics on | `m_lmfinit.f90` |
| `--vbmonly` | flag | Print approximate VBM/CBM relative to vacuum and exit | `main_lmf.f90` |
| `--getq` | flag | Output channel-resolved occupation Q | `main_lmf.f90` |
| `--v0fix` | flag | Freeze the spherical potential V₀ during SCF | `m_lmfinit.f90` |
| `--keepmixm` | flag | Keep the mixing history `__mixm.<sname>` of a previous `lmf` run (without it `lmf` discards the file at start) | `m_lmfp.f90` |
| `--noinv` | flag | Disable inversion symmetry in BZ sampling | `m_lmfinit.f90` |
| `--nosym` | flag | Disable spatial symmetry (forces no-symmetry SCF) | `m_lmfinit.f90` |
| `--nosymdm` | flag | No symmetry on the LDA+U density matrix | various |
| `--shorten` | flag | Shorten lattice vectors (`lmchk` action) | `lmaux.f90` |
| `--slat` | flag | Print Ewald shorted lattice (`lmchk`) | `lmaux.f90` |
| `--afsym` | flag | Use AF (Shubnikov) symmetry; halves the BZ sampling | `m_qplist.f90` |
| `--onesp` | flag | Force `nspxx=1` even when `nsp=2` | `m_qplist.f90` |
| `--avoidgamma` | flag | Skip Γ-point in q sampling | various |

## Quit-points (debugging / staged runs)

| Option | Effect | Used in |
|---|---|---|
| `--quit=show` | Stop after `--show` echo | `m_lmfinit.f90` |
| `--quit=ham` | Stop after `m_hamindex_init` | `main_lmf.f90` |
| `--quit=mkpot` | Stop after `m_mkpot` (potential build) | `m_bndfp.f90` |
| `--quit=dmat` | Stop after first density-matrix construction | `main_lmf.f90` |
| `--quit=band` | Stop after the band loop (no density update) | `m_lmfinit.f90` |
| `--quitecore` | Stop after core treatment | `m_bndfp.f90` |

## Bands / DOS / Fermi surface

| Option | Type | Effect | Used in |
|---|---|---|---|
| `--band` | flag | Compute eigenvalues along `syml.<sname>` (band plot) | `m_lmfinit.f90` |
| `--fullmesh` | flag | Full BZ mesh; used by Fermi surface / PROCAR / boltztrap | `m_qplist.f90` |
| `--fermisurface` | flag | Fermi surface mode (Xcrysden `.xsf`) | `m_bndfp.f90` |
| `--mkprocar` | flag | Emit PROCAR (atomic/orbital projected weights) | `m_bndfp.f90` |
| `--writepdos` | flag | Write `pdos.*` files | `main_lmf.f90` |
| `--dos` / `--pdos` / `--tdos` | flag | DOS / partial DOS / total DOS | `m_lmfinit.f90` |
| `--tdostetf` | flag | Disable tetrahedron for `--tdos` (use Methfessel–Paxton instead) | `m_bndfp.f90` |
| `--wdsawada` | flag | Sawada-style simple DOS output, then exit | `main_lmf.f90` |
| `--eigen-at-k` | flag | Write eigenvalues at user-supplied k | `m_qplist.f90` |
| `--writeeigen` | flag | Same as `--eigen-at-k` (alias) | `m_bndfp.f90` |
| `--allband` | flag | Output all bands at every k | `m_bandcal.f90` |
| `--boltztrap` | flag | Write `BoltzTraP` input on a full mesh | `m_bndfp.f90` |
| `--cls` | flag | Core-level spectroscopy weights (see `m_clsmode.f90`) | `m_bndfp.f90` |
| `--corehole` | flag | Core-hole calculation | various |
| `--efermi=<file>` | kv | Fermi-level file that `lmf` writes (and reads in band-plot mode) and `mlo` reads, instead of `efermi.lmf`. `job_band`, `job_dos`, `job_fermisurface` pass `efermi.lmf.job_band` etc. and `job_mlo_soc` passes `efermi_soc`, so that `efermi.lmf` keeps the value of the last self-consistent `lmf` (2026-09-26) | `m_cmdopt_registry.f90` (`efermi_file`) |

## GW driver and friends (`lmf --jobgw=…`, `gwsc` internals)

| Option | Type | Effect | Used in |
|---|---|---|---|
| `--jobgw=0` | kv | `lmf` writes the index file `__HAMindex0` (for `qg4gw` and `gwinit`) and stops | `main_lmf.f90` |
| `--jobgw=1` | kv | `lmf` produces the GW input (replaces legacy `lmfgw-MPIK`). `hgw --jobgw=1` is the GW step of `gwsc` | `main_lmf.f90`, `main_hgw.f90` |
| `--writeham` | flag | Write `HamiltonianPMT.*` for downstream tools (mlo / writeham) | `main_lmf.f90` |
| `--writev0` | flag | Write V₀ (spherical part) on the radial mesh | `locpot.f90` |
| `--writesene` | flag | Write self-energy file `__sene.*` (debug) | `m_bndfp.f90` |
| `--wsig_fbz` | flag | Write `sigm` reshaped to the full BZ | `m_bndfp.f90` |
| `--use_sigm_fbz` | flag | Read `sigm_fbz` instead of the regular `sigm.*` | `m_bndfp.f90` |
| `--readQforGW` | flag | `qg4gw`: take the q-vectors of `[gw] QforGW` (one-shot GW, `gw_lmfh` passes it) | `m_q0p.f90` |
| `--n1n2n3eps` | flag | Honour `n1n2n3eps` mesh in `qg4gw`/`mkqg` for ε/χ⁺⁻ | `mkqg.f90` |
| `--job=N` | kv | Stage selector in `hbasfp0`, `hvccfp0`, `hsfp0_sc`, `hsfp0`, `hx0fp0`, `hqpe`, `qg4gw`, `heftet`. Range is per-binary; the values that `gwsc` uses are in the console output shown in [gwsc](./gwsc) | many |
| `--ntqxx` | flag | Limit number of QP states (debug) | various |
| `--nk=` | kv | k-parallel group size in `mlo_magnon` (overrides the default from the number of ranks) | `main_mlo_magnon.f90` |
| `--cutuu=`, `--nww=`, `--nb=`, `--Wtype=` | — | 廃止（2026-10-02、AHC・Wannier のマグノンを外した）。与えると止まる | `m_cmdopt_registry.f90` |
| `--dumpW` | flag | `hgw`: write W−v, which is kept in memory, to `__WVR.<iq>` / `__WVI.<iq>` for analysis (`gwsc --dumpW`) | `main_hgw.f90` |
| `--WVR2ptRaxis` | flag | Pole term of Σc on the real axis: W_c at the pole by linear interpolation of 2 points instead of 3 points (`gwsc --WVR2ptRaxis`) | `m_sxcf_sc.f90` |
| `--wcsmear` | flag | Same as `[gw] wcsmear = true` (the key is preferred; it is true by default) | `m_sxcf_sc.f90` |
| `--sigma_tf32` | flag | MP GPU programs at `--use_fp32`: run only the products of Σc (W·zmel on both axes, zsec) with 10-bit-mantissa inputs (TF32, or FP16 when the table chooses `realhgemm`), with the rows of level `tf32` of the linear-algebra table (`gwsc --prec=tf32` passes it) | `m_sxcf_sc.f90`, `m_linalg_policy.f90` |
| `--tetwt_write` | flag | `hgw`: only write the tetrahedron weights of every q to `__TETWT.<iq>.<isp>` (no GPU, no W, no Σ) and stop. `gwsc --gpu` runs this on the CPU cores next to `hgw`, which then reads the files (bit-identical; a missing or non-matching file is computed in `hgw`) | `main_hgw.f90`, `x0kf_v4h.f90` |
| `--qibzonly` | flag | Compute only IBZ q-points (skip auxiliary q₀) | various |
| `--geteta` | flag | Get η damping diagnostics for ε(q,ω) | various |
| `--interbandonly` / `--intrabandonly` | flag | Limit interband / intraband transitions in ε | various |
| `--testso` | flag | SO sanity test path | various |
| `--tetwtk` | flag | χ0: compute the tetrahedron weights k by k in the k loop instead of for all k at once (less memory) | `x0kf_v4h.f90` |
| `--tetraw` | flag | PDOS mode of `lmf`: write the tetrahedron data `tetradata.dat` and the k list `qplistf.dat` | `m_procar.f90` |

## MLO / Wannier-replacement (`job_mlo`, `job_mloW`, `job_mlo_soc`)

| Option | Type | Effect | Used in |
|---|---|---|---|
| `--mlo` | flag | MLO mode. `lmf --jobgw=1 --mlo` writes `__cmlo.data`/`.info` (step a') when `HamRsMLO` exists; `mlo --mlo` builds `HamRsMLO` (the MLO index) and, when `__HamiltonianGW.info` is there (job_mloW, job_mlo_magnon), `__cmlo` from it | `sugw.f90`, `m_HamPMT.f90` |
| `--mlog` | flag | MLO + GW combined output | various |
| `--mlo_diagnorm`, `--mlo_feb4`, `--mlo_ortho`, `--mlo_orthonorm`, `--gs` | — | 廃止（2026-10-02）。与えると止まる。MLO は実空間で 2 乗積分が 1 になるよう、軌道ごとに k によらない定数で規格化する（[mlo](./mlo) §6）。k ごとの規格化と直交化は、軌道の形と内挿したバンドを変えるので使わない | `m_cmdopt_registry.f90` |
| `--mlofreeze` | flag | **MLO-gwsc（開発中）**: 既存の `HamRsMLO` の MLO 索引（種にする MTO チャネルと窓の基準）を使い、`QMLO_SigRs` だけ書き直す（Σ^MLO(q) → 実空間。H・O は扱わず、`__HamiltonianPMT` も読まない、2026-09-27）。$\tilde\chi$ の係数は毎反復作り直す。`gwsc --mlo` が反復内の `mlo` に渡す。→ [MLO-gwsc](./mlo_gwsc) | `m_HamPMT.f90` |
| `--dwnb=mlo` | kv | `huumat_MPI`: the basis of ⟨u_k\|u_k+b⟩ is the MLOs (`job_mlo_magnon`, `mlo_spread.py`). `--dwnb=wan` (the Wannier functions) was removed 2026-10-02 | `main_huumat.f90` |
| `--ahc`, `--mloahc`, `--AHCMAT`, `--UUMAT`, `--cmlo`, `--q2q1test` | — | 廃止（2026-10-02、AHC・lmfham2・試験用）。与えると止まる | `m_cmdopt_registry.f90` |
| `--writedw` | flag | Default-on: write `dw` weights in `m_procar` | `m_procar.f90` |
| `--nowritedw` | flag | Suppress `dw` writing | `m_procar.f90` |
| `--socmatrix` | flag | Write SOC matrix (`__HamiltonianPMTsoc`) | `m_HamPMT.f90` |

## Density / potential I/O

| Option | Type | Effect | Used in |
|---|---|---|---|
| `--density` | flag | Write smoothed density to `xsf` / file (see `m_bndfp.f90`) | `m_bndfp.f90` |
| `--wrhomt` | flag | Write true MT density `rhoMT.<ib>` (debug) | `locpot.f90` |
| `--wpotmt` | flag | Write true MT potential `vtrue.<ib>` | `locpot.f90` |
| `--vesatom` | flag | Output `lmfa` electrostatic potential | `freeat.f90` |
| `--vesdat` | flag | Write the per-rank electrostatic data | `locpot.f90` |
| `--espot` | flag | Electrostatic-potential calculation | `m_bndfp.f90` |
| `--estaticall` | flag | Output full electrostatic potential | `m_bndfp.f90` |
| `--eszero` | flag | Set Madelung correction to zero (diagnostic) | various |
| `--novxc` | flag | Skip XC contribution (debug) | `locpot.f90` |

## Structure / symmetry (`lmchk` actions)

| Option | Type | Effect | Used in |
|---|---|---|---|
| `--getwsr` | flag | Compute touching-MT radii; used by `ctrlgenToml.py` | `lmaux.f90` |
| `--kchk` | flag | Check whether the k-mesh is compatible with the space group | `bzmesh.f90` |
| `--ylmc` | flag | Print `Ylm` analysis | various |

## SOC

| Option | Type | Effect | Used in |
|---|---|---|---|
| `--socmatrix` | flag | (see MLO section) | `m_HamPMT.f90` |
| `--skiphammsoc` | flag | Don't add the SOC term to `hamm` (perturbation path) | `m_bandcal.f90` |
| `--ctrlg:ham.phispinsym=true` | TOML | Spin-averaged radial wfns (formerly `--phispinsym`); set via TOML override or in `[ham]` of `ctrlg.<sname>.toml`. | `m_toml_override.f90` |

## GPU / linear-algebra knobs

| Option | Type | Effect | Used in |
|---|---|---|---|
| `--gpu` | flag | Take the GPU code path (`InstallAll.py` grabs `/tmp/gpu.lock`; `gwsc` takes per-GPU locks in `/tmp/ecalj_res/`, see [GPU 版 § 複数のジョブ](./ecaljgpu#同じ機械で複数のジョブを走らせるとき)) | `InstallAll.py`, `m_*.f90` |
| `--linalg=<row>,<row>,…` | kv | Override rows of the GPU linear-algebra table (which method runs each product and the epstilde inverse), e.g. `--linalg=fp32.cgemm.large=realsgemm,fp32.epsinv.large=mixed1`. Row: `<tf32\|fp32\|fp64>.<cgemm\|zgemm\|dgemm\|epsinv>.<small\|large\|minm\|minn\|mink\|minmnk>`; methods `cublas`, `realsgemm`, `realhgemm` (FP16 inputs, FP32 sums; tf32 rows only), `gemmul8[:moduli]`, `lu64`, `mixed1`, `mixed2`. Normally the table comes from `<bindir>/ecalj_linalg_policy.toml` written by `linalgtune_gpu` at installation (`ECALJ_LINALG_POLICY` overrides the path). See [GPU 版 § 自動選択](./ecaljgpu#精度の選び方と行列演算の自動選択-2026-09) | `m_linalg_policy.f90` |
| `--use_gemmul8` | flag | Emulate the large GPU matrix products with GEMMul8 (Ozaki-II on INT8 tensor cores). Single-precision complex: 7 moduli, only when all sizes ≥ 1000 (faster and more accurate than cuBLAS FP32 there). Double precision: 14 moduli from m·n·k ≥ 1e8 (FP64 accuracy, 12× native FP64 on RTX 5090). Needs `InstallAll.py --gemmul8`; `ECALJ_GEMMUL8_MODULI_C`/`_Z`/`_D` and `ECALJ_GEMMUL8_FAST` override the defaults | `m_blas.f90`, `m_gemmul8.f90` |
| `--modifiedGS` | flag | Modified Gram–Schmidt orthogonalisation | various |
| `--diag=default` / `--diag=tridiag` / `--diag=chefsi` | kv | Diagonaliser selector (LAPACK default, tridiag, ChebFSI) | various |

## Post-processing helpers (`hmlo_ovlppair`, eps tools)

| Option | Type | Effect | Used in |
|---|---|---|---|
| `--sp1=` / `--sp2=` | kv | spin selectors in `hmlo_ovlppair` | `hmlo_ovlppair.f90` |
| `--emin=` / `--emax=` / `--ndos=` | kv | DOS energy window / resolution | various |
| `--EfermiShifteV=` | kv | Shift `EF` (eV) before DOS / spectrum | various |

## Diagnostic / debug-only flags

These exist for development and are unlikely to be useful in production
runs. Behaviour is not guaranteed to be stable.

`--debugbndfp`, `--debugpwmat`, `--debugsugw`, `--debugzmel`, `--zmel0`,
`--normcheck`, `--showdmat`, `--x0test`, `--q2q1test`, `--cvK:`,
`--skip_qvalcheck`, `--skipbstruxinit`, `--skipCPHI`, `--skipGS`,
`--skiphammsoc`, `--skip1d` / `--skip2nd` / `--skip2ndd` / `--skip2ndp`
/ `--skip2nds`, `--skipd` / `--skipf` / `--skiplo` (all currently dead
code in `m_HamPMT.f90`), `--wanatom`.

## Run-time TOML overrides (`--ctrlg:<path>=<value>`)

The only accepted override syntax is `--ctrlg:<dotted.path>=<value>`.
Values are TOML-typed (bool lowercase `true`/`false`, strings quoted,
arrays in `[...]`); bare identifiers are auto-quoted so shell
quote-stripping is forgiving.

```bash
lmf si --ctrlg:bz.nkabc=[8,8,8] --ctrlg:ham.scaledsigma=0.8
lmf si --ctrlg:ham.so=1 --ctrlg:ham.nspin=2          # SOC as perturbation
lmf si --ctrlg:ham.phispinsym=true                   # bool lowercase
lmf si --ctrlg:verbose=50 --ctrlg:time=[5,5]         # top-level keys
lmf si --ctrlg:spec.1.r=2.5                          # override SPEC[1].r
lmf si --ctrlg:symgrp=find                           # bare identifier auto-quoted
```

The path mirrors the TOML structure (`section.key` /
`section.<idx>.key` for arrays-of-tables, 1-indexed). Applied
overrides are logged on rank 0. A key that the file does not carry is
added (below the header of its section; a top-level key at the top of the
file; with a new `[section]` when the section is absent), and the log says
so — until 2026-09-30 such an override was skipped and the run went on with
the default value. A misspelled key is added too and is read by no program:
check the log line. See
[`toml_migration`](./toml_migration#run-time-ctrlg-overrides) for
the full spec.

## When this page goes stale

If the source has a flag this page doesn't list (or vice versa), trust
the source. Refresh with:

```bash
lmf --listcmdopt | sort        # the registered names (cmdopt2 with a trailing =)
```

and reconcile by hand. We deliberately do not auto-generate this page
because the categorical grouping and the prose notes are the value;
a flat alphabetical dump is what `grep` already gives you.

## Tab completion (bash)

`InstallAll.py` ships a bash-completion script that knows about both
the cmdopt registry and the TOML override syntax. After installation:

- `lmf <TAB>` (or `lmfa`, `lmchk`) -> offers every `ctrlg.<sname>.toml`
  in the current directory as the positional `<sname>`.
- `lmf nio --w<TAB>` -> the cmdopt0/cmdopt2 names that start with
  `--w` (`--writedw`, `--writeham`, `--writepdos`, ...).
- `lmf nio --ctrlg:<TAB>` (or any prefix like `--ctrlg:bz.`) -> the
  dotted-path keys that **actually exist** in the cwd's
  `ctrlg.<sname>.toml`. `[[spec]]` / `[[site]]` arrays expand into
  `spec.1.r`, `spec.2.r`, ... automatically.
- Trailing space is suppressed for `--ctrlg:` and any cmdopt2 (so
  `--ct<TAB>` becomes `--ctrlg:` with the cursor still on the colon,
  ready for `bz.nkabc=...`).

The completion has two data sources, each with a different freshness:

| Completion target | Source | When it updates |
|---|---|---|
| `--ctrlg:<dotted.path>` | the cwd's `ctrlg.<sname>.toml` parsed live by `python3` (tomllib) at every TAB | always reflects the on-disk TOML |
| `--<flag>` / `--<flag>=` | `BINDIR/ecalj_cmdopts.list` (dumped by `lmf --listcmdopt` during `InstallAll.py`) | refreshed on every `./InstallAll.py` run |

The mapping mirrors the parser:

- A `ctrlg.*.toml` key only appears as a completion if it really exists
  in that file — so a TAB-completed `--ctrlg:` path never triggers the
  "key not found in target section" abort from `m_toml_override`.
- A `--<flag>` only appears if it is registered in
  `m_cmdopt_registry::load_cmdopt0_registry` /
  `load_cmdopt2_registry` — so a TAB-completed flag never triggers
  the "unknown option" abort from `validate_arglist`.

Both `~/bin/ecalj_complete.bash` and `~/bin/ecalj_cmdopts.list` are
written by `InstallAll.py`; a guarded `[ -f .../ecalj_complete.bash ]
&& source ...` line is appended once to `~/.bashrc`. To opt out, run
`InstallAll.py --no-bashrc`; to disable later, delete the
`# >>> ecalj bash completion ... <<<` block from `~/.bashrc`.

The `lmf --listcmdopt` entry point (added together with the
completion plumbing) dumps every registered cmdopt0 / cmdopt2 name
on stdout, one per line, with a trailing `=` on cmdopt2 entries so
the completion script can distinguish them from bare flags:

```text
$ lmf --listcmdopt | head -3
--afsym
--allband
--avoidgamma
$ lmf --listcmdopt | grep '=$' | head -3
--jobgw=
--job=
--nk=
```

Use the same `--listcmdopt` to generate completion data for any
external scripting that needs to know the registry contents at
install time.
