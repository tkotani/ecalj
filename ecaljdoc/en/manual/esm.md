# Effective Screening Medium (ESM) {#esm-effective-screening-medium}

> English version of the page [manual/esm](../../manual/esm.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/esm.md of commit 8ef3c4868 (source blob 0aeeb03f7022).

A slab with a vacuum layer is not a periodic crystal, but a plane-wave (smooth-density) code solves Poisson's equation as if it were. ESM
(Otani & Sugino, [PRB 73, 115407 (2006)](https://doi.org/10.1103/PhysRevB.73.115407); implemented in ecalj by M. Obata, Kanazawa Univ.)
replaces the periodic electrostatics with the boundary condition you actually want (vacuum, or a metal electrode) on
each side. It also lets you apply an electric field to a slab. ESM combined
with QSGW is used in
[PRB 101, 205120 (2020)](https://doi.org/10.1103/PhysRevB.101.205120).

**Turn it on for every slab with a vacuum layer.** Without it the $G=0$ component of the Coulomb potential is fixed by the periodic
convention instead of by the vacuum level, and the whole eigenvalue spectrum sits on a different zero — 4.4 eV on the Fe/MgO sample (`Samples/MLOsamples/FeMgO`, 30.5 a.u. of vacuum) —
while the total energy, `Vesav` and the magnetic moment agree to 8–13 digits, so nothing looks wrong until you compare Fermi energies.

Configured by the `[esm]` section of `ctrlg.<sname>.toml` (below: all keys, the `boundary` values, how to apply a field,
and the conversion of the retired `esm_input.dat` by `ctrlg_absorb.py`).
Sample: `Samples/MLOsamples/FeMgO`.

## The `[esm]` section {#the-esm-section}

```toml
[esm]
boundary  = "vac/slab/vac"  # vacuum(-z)/slab/vacuum(+z)
origin    = -8.63717        # (a.u.) z-translation of the density
shiftmode = 0               # 0: origin is absolute, 1: in units of the cell length
zb        = [23.907104, -23.907104]   # (a.u.) boundaries z1esm, z2esm
potential = [0.0, 0.0]      # (Ry) on the +z / -z sides
field     = [0.0, 0.0]      # (Ry/a.u.) on the +z / -z sides
```

**No `[esm]` section at all means ESM is off.** That is the right thing for
a bulk crystal and the wrong thing for a slab, so `lmf` prints a warning when
`c` is more than 3× the longer of `a`, `b` and no `[esm]` is present.

`ctrlgenToml.py` writes a commented-out `[esm]` template (with the
`# === ESM` heading and one line per key) into every new file, right after
`[ham]`; for a slab, uncomment it and set the numbers with care. A bulk
crystal leaves it commented — no `[esm]` section means ESM off.

Table 1. The keys of `[esm]`

| key | meaning | legacy `esm_input.dat` |
|---|---|---|
| `boundary` | which boundary condition on each side (table 2) | line 1, `jesm` |
| `origin` | z-translation applied to the density before solving. The code uses `-origin` internally, so this is **the same number as line 2** of the old file | line 2, `tresm` |
| `shiftmode` | `0`: `origin` is in a.u.; `1`: in units of the cell length | line 2, `jtresm` |
| `zb` | the two boundary planes `z1esm`, `z2esm` in a.u. Omit to get `±alat·plat[2][2]/2` | line 3 |
| `potential` | fixed potential on the +z / -z side (Ry). Used by `metal/...` and the `:v-e` / `:e-v` forms | line 4 |
| `field` | fixed electric field on the +z / -z side (Ry/a.u.). Used by the `:field` form | line 5 |

### `boundary` values

Table 2. The values of `boundary`

| value | system | `jesm` |
|---|---|---|
| `"off"` (or the section absent) | ESM disabled — plain periodic | 0 |
| `"vac/slab/vac"` | vacuum on both sides | 1 |
| `"metal/slab/metal"` | metal electrode on both sides | 2 |
| `"vac/slab/metal"` | vacuum on −z, metal on +z | 3 |
| `"metal/slab/vac"` | metal on −z, vacuum on +z | 4 |
| `"vac/slab/vac:field"` | vacuum both sides, `field` imposed | 5 |
| `"metal/slab/metal:v-e"` | `potential` on −z, `field` on +z | 6 |
| `"metal/slab/metal:e-v"` | `field` on −z, `potential` on +z | 7 |
| `"periodic:esm"` | same result as `off`, computed through the ESM path | 10 |

A bare integer is also accepted, so mechanically converted files keep
working; the strings are preferred because `1` and `3` are not
distinguishable by eye.

**To apply an electric field to a slab**, use `"vac/slab/vac:field"` with
`field = [E, -E]`, or hold the two electrodes at fixed potentials with
`"metal/slab/metal"` and `potential = [V1, V2]`.

## Migration from `esm_input.dat` {#migration-from-esm-input-dat}

Before 2026-09-16, ESM was configured by a separate positional file
`esm_input.dat` — six lines of bare numbers with no keys. It had two
problems: you could not tell what the numbers meant without reading
`esmsmves.f90`, and **a missing file disabled ESM silently** (`lmf` printed
one line and carried on). Copying a sample directory without that one file
therefore changed the physics without any error.

The Fortran reads `ctrlg.<sname>.toml` and nothing else (since 2026-09-17). A leftover
`esm_input.dat` makes `lmf` / `lmfa` / `lmchk` **abort** with a message
naming the converter — abort rather than ignore, so a slab never runs
without ESM by accident.

```bash
ctrlg_absorb.py <sname>      # esm_input.dat -> [esm] in ctrlg.<sname>.toml
                             # (and PB.<sname>.toml -> [product_basis], if present)
```

Table 3. What `ctrlg_absorb.py` does

| situation | what `ctrlg_absorb.py` does |
|---|---|
| `ctrlg.<sname>.toml` already has `[esm]` | the TOML wins. `esm_input.dat` is moved to `esm_input.dat.bk`, whose header records that its settings were **not** used |
| no `[esm]` in the TOML | the file is converted, the `[esm]` block is inserted (with comments, before the GW sections), and the original is moved to `esm_input.dat.bk` with a header saying where it went |

`Legacy2toml.py` performs the same conversion when it migrates a legacy
directory.

## Worked example — Fe/MgO slab

`Samples/MLOsamples/FeMgO` (`c/a = 8.50`, 30.5 a.u. of vacuum, 7 empty
spheres filling it):

```toml
[esm]
boundary  = "vac/slab/vac"
origin    = -8.63717
shiftmode = 0
zb        = [23.907104, -23.907104]   # alat*plat[2][2]/2 = 23.903861; set by hand here
potential = [0.0, 0.0]
field     = [0.0, 0.0]
```

With and without this section, on the same `rst.<sname>`, table 4 compares the two.

Table 4. With and without `[esm]` (FeMgO)

| | E_F (Ry) | Vesav (eV) | `val*vef` (Ry) |
|---|---|---|---|
| `[esm]` present | −0.38286150 | −2.12295183 | −257.64 |
| `[esm]` absent | −0.05833040 | −2.12295183 | −253.66 |

`Vesav` agrees to 8 digits, so the density and the electrostatics are the
same — only the zero of energy moves, by 0.3245 Ry = 4.4155 eV. Check
`llmf_ef` (or any `lmf` log) for which one you got:

```
effective screening medium method jesm=  1        <- ESM on
system : vaccum(-z)/slab/vaccum(+z)

esmsmves: ESM is off (no [esm] section in ctrlg.<sname>.toml)   <- ESM off
```

## Where it is implemented {#where-it-is-implemented}

`SRC/subroutines/esmsmves.f90` (the solver; called from `smves`) and
`SRC/subroutines/m_lmfinit.f90::readesm` (reads `[esm]`; stops if `esm_input.dat` is left).
