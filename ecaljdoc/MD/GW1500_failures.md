# GW1500 production failures — record and --fp32 recovery

Failed calculations from the GW1500 QSGW80 production run on kt1
(`~/DATA/gw1500/`, snapshot 2026-05-10). These crashes are valuable evidence:
most trace to the TF32 mixed-precision bug documented in this sample
(`README.md`). Raw signatures are in `GW1500_failed.log` (verbatim copy of the
production `failed.log`).

## Categories (62 total)

| cat | signature | count | --fp32 |
|---|---|---|---|
| A | `lmf failed even after reducing bmix to minimum` | 6 | fixed |
| B | `NaN detected in lsc` (correlation) | 8 | fixed |
| C | `NaN detected in lqpe` | 22 | fixed |
| D | `TIMEOUT after 8h` | 22 | not addressed (speed, separate axis) |
| E | KILLED (interrupted) ×3 + `Command failed lmf` ×1 | 4 | rerun / legacy path |

[A]+[B]+[C] = 36 precision-caused failures. All 36 were re-run with the modern
flow (`ctrls -> ctrlgenToml.py --ssig=0.8 -> gwsc --gpu --mp --fp32`, 1 QSGW
iteration from LDA): **36/36 completed with no NaN and no divergence**
(sigma_m back to O(100-1000) vs the broken O(50000)). mp-8196 and mp-581833
were further verified to match the CPU reference; mp-8196 converges to 0.1 eV
in the full gwscconv flow (gap 5.05 eV).

## [A]+[B]+[C] — precision failures and --fp32 1-iter result

| mpid | nat | 組成 | cat | --fp32 | sigma_m |
|---|---|---|---|---|---|
| mp-27419 | 8 | LiBiF | A | OK | 548.1 |
| mp-545730 | 8 | AgCNO | A | OK | 798.9 |
| mp-29738 | 5 | TlClO | A | OK | 704.7 |
| mp-551758 | 6 | AgClO | A | OK | 568.9 |
| mp-626151 | 8 | YHO | A | OK | 1215.1 |
| mp-8196 | 5 | AgNO3 | A | OK | 539.6 |
| mp-1080045 | 8 | SrHBrO | B | OK | 799.9 |
| mp-20458 | 8 | PtPbF | B | OK | 358.2 |
| mp-2247 | 8 | AgN | B | OK | 723.3 |
| mp-22604 | 8 | InAsF | B | OK | 621.6 |
| mp-3078 | 8 | CdSiAs | B | OK | 208.3 |
| mp-570140 | 8 | AuBr | B | OK | 199.0 |
| mp-570589 | 8 | SeBr | B | OK | 484.3 |
| mp-571297 | 8 | AgN | B | OK | 498.8 |
| mp-581833 | 4 | RbN3 | C | OK | 339.8 |
| mp-22981 | 5 | TlIO | C | OK | 366.0 |
| mp-22535 | 5 | HfPbO | C | OK | 206.4 |
| mp-29798 | 5 | TlBrO | C | OK | 407.6 |
| mp-625548 | 5 | CdHO | C | OK | 381.9 |
| mp-1077553 | 6 | SrBBrN | C | OK | 177.5 |
| mp-22993 | 6 | AgClO | C | OK | 960.5 |
| mp-30530 | 6 | TlClO | C | OK | 628.5 |
| mp-36248 | 6 | HBrN | C | OK | 597.9 |
| mp-3307 | 6 | BaSO | C | OK | 577.5 |
| mp-30302 | 6 | HfBrN | C | OK | 332.5 |
| mp-551798 | 6 | AlPO | C | OK | 243.9 |
| mp-5606 | 6 | AlTlF | C | OK | 374.5 |
| mp-567299 | 6 | KCIN | C | OK | 683.6 |
| mp-552934 | 6 | BaCuBrO | C | OK | 216.3 |
| mp-1077901 | 7 | InCoSnS | C | OK | 275.1 |
| mp-7700 | 7 | SiB | C | OK | 79.5 |
| mp-29047 | 8 | TlVO | C | OK | 661.3 |
| mp-510557 | 8 | CsN | C | OK | 640.0 |
| mp-552169 | 8 | AgClO | C | OK | 811.2 |
| mp-552604 | 8 | YBiClO | C | OK | 206.7 |
| mp-870 | 8 | TlN | C | OK | 527.3 |

## [D] TIMEOUT (22) — NOT fixed by --fp32 (speed / convergence difficulty)

- mp-6781
- mp-685022
- mp-691
- mp-696736
- mp-720
- mp-730101
- mp-8348
- mp-8360
- mp-8361
- mp-8446
- mp-863757
- mp-571266
- mp-9778
- mp-9922
- mp-9687
- mp-570219
- mp-571100
- mp-9273
- mp-571442
- mp-9274
- mp-573321
- mp-570076

## [E] other (4)

- KILLED (interrupted when GW1500 was stopped; rerun resolves): mp-570051 mp-568388 mp-23435 
- Command failed lmf --jobgw=1 (old bin2 path): mp-1179832 

## Where the raw data lives (kt1)

- `~/DATA/gw1500/failed.log` — master signatures (copied here as `GW1500_failed.log`).
- `~/DATA/gw1500/*_FAILBACKUP/` — per-material backups of failed runs (~69 GB;
  many empty or heavy binaries; at risk under disk pressure — the durable record
  is this document + GW1500_failed.log).
- `~/DATA/gw1500fp32/*_v2/` — the --fp32 re-runs (summary_v2.tsv).
- Detailed 3-way evidence (CPU / TF32 / FP32) for mp-8196: `reference/` in this sample.

## Not yet done
- Full gwscconv convergence (not just 1 iter) for the 35 beyond mp-8196.
- CPU ground-truth spot checks beyond mp-8196 / mp-581833.
- [D] TIMEOUT class needs a separate fix (max-iter, heavy-element convergence).
