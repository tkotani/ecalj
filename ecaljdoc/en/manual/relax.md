# Forces and atomic position relaxation {#forces-and-atomic-position-relaxation}

> English version of the page [manual/relax](../../manual/relax.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/relax.md of commit 1fe9c7640 (source blob edb69065eb51).

Atomic positions can be relaxed in LDA/GGA (no cell optimization). Sample: [`Samples/Relax/LaGaO3`](https://github.com/tkotani/ecalj/tree/main/Samples/Relax/LaGaO3) (orthorhombic perovskite, 20 atoms).

Set the `[dyn]` section of `ctrlg.<sname>.toml` (Table 1). When `[dyn].mode` is not 0, the forces are computed too.
The directions to relax are given in each `[[site]]` as `relax = [0, 0, 1]` or so (only the directions with 1 move).

Table 1. Keys of `[dyn]` (defaults from `m_lmfinit.f90`)

| key | default | meaning |
|---|---|---|
| `mode` | 0 | 0: no relaxation, 4: conjugate gradients, 5: Fletcher-Powell, 6: Broyden |
| `nit` | 1 | maximum number of relaxation steps |
| `hess` | true | read the Hessian matrix |
| `xtol` | 1e-3 | convergence criterion in displacements; >0: use length, <0: use max value, =0: do not use |
| `gtol` | 0 | convergence criterion in gradients; >0: use length, <0: use max value, =0: do not use |
| `step` | 0.015 | initial (and maximum) step length |
| `nkill` | 0 | remove the Hessian after this many steps |

The values of the LaGaO3 sample:
```toml
[dyn]
mode = 5
nit  = 20
hess = true
xtol = 0.0001
gtol = 0.0
step = 0.015
```

* Caution: in the relaxation mode the atomic-position file `AtomPos.<sname>` is used. If it exists, the positions written there replace those of `ctrlg.<sname>.toml`
  (`call ReadAtomPos(irpos)` in `SRC/subroutines/main_lmf.f90`). Delete an old `AtomPos.<sname>` when you start anew.
* The relaxation of the perovskite LaGaO3 was checked to work well. The relaxation mode is not used much, so use it with care.
