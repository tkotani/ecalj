# m_tetrakbt (finite-temperature tetrahedron) — KNOWN BUG REPORT

**Status: BROKEN. Do NOT enable `tetrakbt = true` for production until fixed.**
Date: 2026-06-08

## TL;DR
`m_tetrakbt` (the finite-temperature tetrahedron weights for chi0, selected by
`[gw] tetrakbt = true`) activates correctly but produces a **wrong chi0**. For metals
this makes the dielectric matrix singular at some q and the screened Coulomb `W` comes
out **NaN** (then `SEc` = NaN). The bug exists even in the T->0 limit. It is an
"idea-only sketch" (cf. the in-source comment *"Temperature factor kbt will be
implemented by H.Okumura soon!"* and the leftover `if(.TRUE.)/else` dev scaffolding).

## Symptoms
- `[gw] tetrakbt = true` (any `t_tetrakbt`): the run prints `tetrakbt_init: T[K] ...`
  and `tetwt5: ... usetetrakbt:= T` (activation + kBT are correct), `SEx` is unchanged
  (exchange does not use chi0), but **`SEc` = NaN** at certain q.
- Example: Na, GW 4x4x4, the NaN first appears at q = (0.25,0.25,0.5) (all real-frequency
  W elements NaN). Verified the NaN originates in the `epstilde` inversion (`m_llw`
  `zminv`), i.e. chi0 is finite-but-wrong -> epsilon singular -> W NaN.

## How it was diagnosed (reproducible)
One-shot comparison on Na against the proven T=0 routine `lindtet6`:
1. Build a Na GW dir (`lmfa; lmf --jobgw=0; gwinit`), set `[gw] tetrakbt=true`,
   `t_tetrakbt=1.0` (~T=0).
2. In `tetwt5x_dtet4`, when `usetetrakbt`, also call `lindtet6` into a scratch array and
   compare per-tetrahedron `sum(wtthis(:,0))` (both columns are the same quantity:
   non-matrix-linear weight = `pvn*(integb-intega)`; see `intttvc6` line ~796 vs
   `tetrakbt` accumulation).
3. Result at T=1K (essentially T=0, should match lindtet6 exactly):
   - **34489 tetrahedra disagree by >10%**; **max overestimate = 4451x (445100%)**.
   - **99.8% of the bad tetrahedra are Fermi-surface-crossing** (occupied band straddles
     E_F: corner energies eocc-E_F of mixed sign).

## Root cause (algorithm)
The correct weight is
```
W(ihis) = ∫_{bin} dω ∫_tet d3k [f(ea)-f(eb)] δ(ω-(eb-ea))
```
`lindtet6` (T=0) computes this **exactly** by geometric subdivision: it splits the
tetrahedron at ea=E_F (occupied) and eb=E_F (unoccupied) and integrates the energy
delta only over the occupied×unoccupied sub-volume — correct regardless of the ω-bin mesh.

`m_tetrakbt` instead uses a **midpoint factorization**:
```
W(ihis) ≈ G(ihis) × <f(ea)-f(eb)>           (eaf_triangle evaluated at the bin CENTER ω)
   G(ihis) = pvn*(integb-intega)            (geometric δ-weight over the WHOLE tetrahedron)
```
This is only valid when the occupation is ~constant across each ω-bin. It **breaks** when:
- the FS crosses the tetrahedron (occupation varies with ω), AND/OR
- the ω-bin is wider than the tetrahedron's transition-energy spread — which is the norm
  at high ω on the exponential `frhis` mesh (`HistBin_dw`, `HistBin_ratio`).

When the occupied corner sits at a different ω than where the geometric weight is
concentrated, the bin-center (plus the `eaf_triangle` ecenter clamping at `el`/`eh`)
evaluates the occupation at the **wrong ω** -> the weight lands in the **wrong bin** with
the **wrong (huge) magnitude**.

### Smoking-gun per-bin dump (worst Na tetrahedron, T=1K)
```
eocc-E_F = [-0.0081, 0.1257, 0.1257, 0.1257]   (one corner barely occupied)
eunocc-E_F ≈ 3.8 (fully unoccupied)
  bin 137 (ω=1.840):  lindtet6 = 1.17e-6   tetrakbt = 0.0      <- correct weight, missed
  bin 138 (ω=1.934):  lindtet6 = 0.0       tetrakbt = 1.32e-4  <- wrong bin, ~113x too big
```
The occupied corner has the LOW transition energy (1.84); the bulk geometric weight is at
the HIGH transition energy (1.93). The factorization multiplies the bulk geometric weight
by a mis-evaluated occupation -> gross overestimate at the wrong frequency.

## Fix attempts so far
- **factri eps-lifting / degeneracy guards**: did NOT help (the bug is not in `factri`;
  forcing `factri = factri0` (exact T=0 occupation) still left all 34489 mismatches).
- **(A') sub-bin ω-integration of the occupation** (integrate `g(ω)<f>(ω)` over fine
  sub-bins instead of bin-center): **works** — max error 4451x -> 0.77, numerically stable
  and mass-conserving — but a naive uniform `Nsub=200` is **far too heavy** (~hundreds×
  more `integtetn`/`eaf_triangle` calls; Na hgw went from ~1 min to >12 min and unfinished).
  Residual 0.77 at T=1K is the sharp-step resolution limit of sub-binning (vanishes at
  realistic finite T, but the architecture is still approximate + expensive).

## Recommended fix (planned)
**(B') Energy convolution** — exact and minimal:
```
f(e) = ∫ dE (-f'(E)) θ(E-e)   =>   W_T(ω) = ∫ dE (-f'(E)) · lindtet6(ω; E_F = E)
```
Call the proven exact T=0 routine with the Fermi level shifted to E, weighted by the
thermal kernel `-f'(E)` (a bell curve of width kBT). Advantages: geometry is exact (no
sub-bin error), reuses `lindtet6` verbatim, only new code is the E-quadrature loop
(Gauss, NE~16), cost controllable (apply the NE× only to tetrahedra within ~10 kBT of the
FS; far tetrahedra need lindtet6 once). Abandons the broken midpoint factorization
entirely.

Sequencing (agreed plan):
1. (this report)
2. Simplify the current tetrahedron code into an equivalent, cleaner form.
3. Add finite-temperature on top via (B').
