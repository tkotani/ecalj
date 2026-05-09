# GetStarted/GaAs — minimal worked example for the tutorial

This directory carries the GaAs sample referenced from the
[GetStarted](https://ecalj.github.io/ecaljdoc/manual/README_tutorial#getstarted)
walk-through in `ecaljdoc`. It is the smallest end-to-end input
seed: just enough to follow Steps 1 → 6 of the tutorial without
hunting for a POSCAR.

## What's in here

| file | role |
|---|---|
| `ctrls.gaas` | Step 1 output (lightweight structure-only seed; lattice + 2 sites: Ga, As). |
| `ctrlG.gaas.toml` | Step 2 output (full TOML input read by `lmf` / `lmfa` / `gwsc`). Generated from `ctrls.gaas` by `ctrlgenToml.py gaas`. |
| `PB.toml` | Step 2 by-product (per-atom product-basis tables, GW path only — normally not edited by hand). |

## How to reproduce

In a clean copy of this directory (say `cp -r GaAs gaas_work && cd gaas_work`):

```bash
# Step 2: regenerate the TOML pair from the lightweight seed.
ctrlgenToml.py gaas

# Step 3: LDA self-consistency.
lmfa gaas
mpirun -np 4 lmf gaas

# Step 4-5: band path + plot.
getsyml gaas
job_band gaas -np 4

# Step 6: QSGW (small loop count).
gwsc -np 8 1 gaas
```

The tutorial explains each step in detail and shows what to expect
in the console output.

## Source

`ctrls.gaas` is the same content the tutorial shows at
[Step 1](https://ecalj.github.io/ecaljdoc/manual/README_tutorial#step-1-convert-poscar-to-ctrls):

```text
#id = GaAs
%const bohr=0.529177 a=5.65325/bohr
STRUC
     ALAT={a}
     PLAT=0 0.5 0.5  0.5 0 0.5  0.5 0.5 0
SITE
     ATOM=Ga POS=0.0 0.0 0.0
     ATOM=As POS=0.25 0.25 0.25
```

`alat` evaluates to `5.65325 Å / 0.529177 Å·bohr⁻¹ = 10.6831 a.u.`
when `ctrlgenToml.py` bakes the `%const` math in.
