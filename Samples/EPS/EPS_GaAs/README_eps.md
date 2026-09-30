# ε(q,ω) in the RPA: GaAs

The manual is ecaljdoc `manual/optical.md`; this note is about the settings of this sample.

## Run

1. Self-consistent calculation first (`lmfa`, `lmf`; for QSGW, `gwsc`, which leaves `sigm.<sname>`).
2. Then
   ```
   job_eps gaas -np 8               # without the local-field correction (fast)
   job_eps gaas -np 8 --lcf         # with the local-field correction (much slower)
   job_eps gaas -np 8 --decompose   # interband and intraband parts separately (the test of this sample)
   ```

## Settings in `ctrlg.gaas.toml`

### q points: `[gw] QforEPS`

q = 0 cannot be used, so ε is computed at small q along a line and the q → 0 limit is taken from them
(or the smallest q is used). Too small a q is numerically unstable.

```toml
[gw]
QforEPSau = true   # QforEPS in a.u.; false (default): in units of 2π/alat
QforEPS = """
 0 0 0.00050
 0 0 0.00100
 0 0 0.00200
"""
```

### k mesh: `[gw] n1n2n3`

Sharp structures of Im ε need many k points (about 20×20×20 for GaAs). Without `--lcf` this is affordable;
with `--lcf` start from about 12×12×12 (ε(0,0) of MgO is already reasonable there).

### Mixed product basis (for `--lcf`)

To make `--lcf` affordable, reduce `[product_basis] pb_lcutmx` (e.g. 2 per atom) and `[gw] QpGcut_cou` (about 2.0).
Check with a small k mesh that the result does not change.

### Energy mesh: `[gw] HistBin_dw`, `HistBin_ratio`

A finer mesh along the real axis (ω_i = b(e^{a(i-1)} − 1), a = ratio − 1, dw = b·a) costs little:

```toml
HistBin_dw    = 2e-4   # (a.u.) bin width at ω = 0; 1e-5 is fine (metals)
HistBin_ratio = 1.01   # 1.03 is safer
```

## Output

`EPS0001.nlfc.dat`, `EPS0002.nlfc.dat`, … for the 1st, 2nd, … q of `QforEPS` (`nlfc`: no local field; with `--lcf`
also `EPS000n.dat`). With `--decompose`, `.interbandonly` and `.intrabandonly`. Columns: q (1–3), ω in Ry (4), Re ε and Im ε (5, 6), Re and Im of 1/ε (7, 8):

```gnuplot
set datafile fortran
set xrange [0:20]
plot "EPS0001.nlfc.dat" u ($4*13.605):($6) w lp
```

## Caution

The √(E − E0) onset of Im ε at the absorption edge is hard to reproduce with the tetrahedron method unless the k mesh is
very fine (80×80×80 for GaAs). ε(q → 0, ω = 0) is much easier to converge.
