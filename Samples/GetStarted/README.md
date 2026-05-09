# GetStarted — minimal samples for the ecaljdoc tutorial

The GetStarted walk-through in
[ecaljdoc](https://ecalj.github.io/ecaljdoc/manual/README_tutorial#getstarted)
takes you from a structure file to a QSGW band plot. To save users
hunting for a starting POSCAR / `ctrls`, each subdirectory here
ships the minimal input pair (`ctrls.<sname>` + `ctrlG.<sname>.toml`
+ `PB.toml`) that the tutorial assumes.

## Samples

| dir | system | relevant tutorial steps |
|---|---|---|
| [GaAs/](./GaAs/) | GaAs (zinc-blende, 2 atoms / cell, non-magnetic) | Steps 1-6 (LDA → band plot → QSGW) |

Each sub-directory has its own `README.md` explaining what's
inside and how to reproduce.
