# TODO (what is left in development) {#dev-todo}

> English version of the page [manual/dev_todo](../../manual/dev_todo.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/dev_todo.md of commit 418976791 (source blob f1613824a485).

::: info About this page
A list of what is not yet done in the development of ecalj (made 2026-11-17).
The master copy is [`MD/TODOandQuestion.md`](https://github.com/tkotani/ecalj/blob/main/MD/TODOandQuestion.md) in the repository;
this page takes from it the items that are useful to people outside. The progress of work on each computer is not copied. Where they differ, the master copy is right.
The date in parentheses in each item is the day it was found. Marks: **[not started]**, **[in progress]** started and halfway, **[undecided]** the policy is not decided, **[kept]** not urgent.
:::

## Open items written in other pages {#in-pages}

They are not copied here; these point to the sections.

- MLO: [Remaining problems](./mlo#mlo-problems) (one fixed $w$ cannot satisfy semiconductors and metals at the same time, etc.)
- MLO-gwsc: the section "Open issues and future work" of [MLO-gwsc](./mlo_gwsc) (comparison with `pwmode = 11`, use with antiferromagnetic symmetry, effect of the number of MLOs, etc.)
- kBT: the section "Remaining issues" of [kBT](./kBT) (how to choose `t_sigmaw` in metals, `EFERMI_kbt` of insulators, the Bose term, etc.; 8 items)
- Unfinished theory pages: [Optical properties](../theory/optical_properties), [Basics 2](../theory/basic2), [Representation of wave functions in GW](../theory/gw_basis),
  [Product basis](../theory/gw_product_basis), [Matrix elements used in GW](../theory/gw_zmel)

## MLO {#mlo}

- **[not started] EuO with EH2 added stops where the model is made** (2026-10-01). With both EH and EH2 as seeds, `zhev_tk4` drops directions close to linear dependence.
  Atoms with 4f are excluded from criterion 2 of [How to choose the model](./mlo#mlo-choice), so this is not urgent
- **[kept] The cause of η ≠ 1, the factor in the Goldstone condition of magnons** (2026-10-02). With Löwdin MLOs, Fe 1.24, FeCo 1.28, Ni 1.75. Why Ni is large is not known
- **[not started] Fe magnons: 20 % too high at small q** (2026-10-02). Löwdin MLOs give 20 % higher energies than the Wannier version at q ≤ 0.3 (0.085 vs 0.068 eV at q = 0.1).
  q ≥ 0.4 agrees. Comparison with the experimental dispersion (neutron scattering) is not done either (the data are not in the repository).
  Sample: [`Samples/Magnon/Fe_mlo_magnon`](https://github.com/tkotani/ecalj/tree/main/Samples/Magnon/Fe_mlo_magnon)
- **[kept] Automatic empty spheres (ES)** (2026-10-01). Whether to build [`ctrlg_addes.py`](https://github.com/tkotani/ecalj/blob/main/SRC/exec/ctrlg_addes.py) into the generation of ctrlg.
  For slabs (structures with thick vacuum), MLO works when the stacking is extended outside the surface atoms and ES (s,p) are put there (made a rule on 2026-10-08). `ctrlg_addes.py` does not yet place them this way
- **[not started] ES lower the symmetry** (2026-10-05). `ctrlg_addes.py` takes non-overlapping voids from the largest, so it can take only part of a set of voids equivalent by symmetry.
  The irreducible q points increase and GW becomes 2–6 times heavier. Take equivalent voids together
- **[not started] `gwinit` should write `mlo_lm` of an ES (z = 0) to match its basis** (2026-10-05). It now writes s,p, so for an ES with s only, a p that is not in the basis enters

## Symmetry and GW {#gw}

- **[not started] Revive EIBZ and symmetrize Σ (and χ0) with the group of operations that keep q** (2026-10-02).
  Changing the order of the symmetry operations changes the QP of the first QSGW iteration of AF NiO (R-3m, 12 operations) by up to 15 meV. This comes from truncations
  that are not compatible with the rotations, such as the upper limit of the number of bands or the cut at `emax_sigm`; it is not a bug (the test reference was set to the values with the present order).
  Averaging Σ(q) as $\Sigma\leftarrow\frac{1}{N_g}\sum_g U_g\Sigma U_g^\dagger$ removes this difference and also reduces the k sum.
  [`eibzgen.f90`](https://github.com/tkotani/ecalj/blob/main/SRC/subroutines/eibzgen.f90) for EIBZ is still there.
  First symmetrize Σ of NiO 4³ with numpy, and see whether the difference goes away and how long it takes
- **[not started] GW without time reversal (`npm = 2`)** (2026-10-02). For a magnet with spin-orbit coupling inside GW and without space inversion, χ0, W and Σc at negative frequencies are needed separately.
  Now χ0 stops at the start of `x0kf_zxq`, and Σ stops at `npm=2 need to be examined` in `m_sxcf_sc_count.f90`.
  Plan: (a) with space inversion, make $W(\mathbf q)$ at negative frequencies from $W_{\mathbf G\mathbf G'}(\mathbf q,-\omega)=W_{-\mathbf G',-\mathbf G}(-\mathbf q,\omega)$,
  (b) without it, carry the frequency integrals over to the negative side. Needed for altermagnets (RuO₂ and MnTe have space inversion, so (a))
- **[not started] Compute the core exchange without the product basis (MPB)** (2026-11-15). Now each `gwsc` iteration makes an MPB for the core before computing Σ^Cx
  (below (I15) in [Implementation notes](./dev_impl)). The core is only inside the MT spheres, so it should be computable directly as integrals inside the MTs
- **[not started] What is left of `hgw`** (2026-04). Share W within a node (`MPI_Win_allocate_shared`), and move Sx and the core exchange of `hsfp0_sc` into `hgw` (low priority)
- **[not started] Memory of `lmf --jobgw=1`** (2026-10-01). Each rank holds `ppovl(ngp,ngp)` of `GEIGpart` (32·ngp² bytes, the square of the cell volume).
  Ideas: decide the parallelism from ngp, keep only the Cholesky factor, make $O\cdot x$ by FFT and solve iteratively

## Finite temperature {#finite-t}

- **[not started] Check how the finite-temperature tetrahedron method relates to the T = 0 one** (2026-10-02). Whether the weights of `t_tetrakbt > 0` (`m_tetrakbt`) agree, as T → 0,
  with those of `t_tetrakbt = 0` (`tetwt5`). With the same k mesh, compare the imaginary part of χ0 and the QP energies, lowering T to 300, 100, 30, 10 K.
  Also the combination with the width `t_sigmaw` on the Σ side. Other issues are in "Remaining issues" of [kBT](./kBT)

## Input, k points, SCF {#input}

- **[in progress] A rule that sets the k mesh by spacing** (2026-10-05). Now the default of `ctrlgenToml.py` is 8×8×8 and GW is 4×4×4 for all materials.
  With the spacing Δk, $n_i=\max(n_{\min},\lceil|\mathbf b_i|/\Delta k\rceil)$; it was decided that lmf uses Δk = 0.26 Å⁻¹ (8 for Si) and GW 0.34 Å⁻¹ (6 for Si, mostly 4 for cells of 8 atoms).
  `kmesh.py` is on the branch `kmesh` and not yet in main. `ctrlgenToml.py` should write the defaults of `[bz] nkabc` and `[gw] n1n2n3` with it
- **[not started] A rule for the k points of GW in metals** (2026-10-05). For Na and Al (sp) and Fe, Ni and Cu (3d), vary the GW mesh and find the step where the QP levels,
  the d band width and the exchange splitting settle within 0.05 eV. Decide whether this is ordered by $k_F$ or by the d density of states, and make a rule that raises the step only for metals
- **[not started] Make the automatic adjustment of the mixing parameter b complete inside lmf** (2026-10-05). Now Python (`bmix_reduction` in [`pylib/dft.py`](https://github.com/tkotani/ecalj/blob/main/SRC/exec/pylib/dft.py))
  lowers b and restarts lmf from the beginning every time lmf fails. This is slow, and does not help people who run lmf directly.
  Idea: in the SCF loop of `m_lmfp.f90`, on signs of divergence (qdiff growing several times in a row, NaN), go back to the density at the start of the iteration, discard the mixing history, halve b and continue.
  Write it at the end of the line in the save file, like `b=0.20->0.10 at it 7`
- **[undecided] Settings for molecules** (2026-11-14). `ctrlgenToml.py --systype=molecule` writes `tetra = true`, `w = 0.01`, but
  the comments of `toml_comments.py` and [Samples/AtomDimer/N2](https://github.com/tkotani/ecalj/tree/main/Samples/AtomDimer) use `tetra = false`, `w = 0.001`.
  When molecules are computed, compare both and make one agree with the other
- **[undecided] Temperatures in the test inputs**. The tests use `t_tetrakbt = 262`, `t_sigmaw = 0`. Making them the template values (300, 300) moves the references of gas_gwsc, fe_gwsc and others
- **[not started] `rst` and `dmats` of `Samples/MLOsamples/RuO2`** were made while LDA+U with `pwmode = 11` had an error (2026-03-30 to 09-30). Remake them

## Code {#code}

- **[in progress] Comments around the electrostatics and the core in PMT** (2026-11-13). The wrong comments were fixed. Left: (a) the equation numbers in the comments are those of the paper (JPSJ 84, 034702);
  add those of [Band-structure method of ecalj](../theory/pmt) too, (b) remove the old `locpt2` turned into comments after `locpot` (about 200 lines)
- **[not started] Bring in code for the anomalous Hall conductivity (AHC)** (2026-10-02). The former AHC (`hahc`, `job_AHC`) was removed together with the Wannier path (git tag `last-wannier`).
  Put new code in
- **[not started] What is left of AFTEST (fixing antiferromagnetic moments)** (2026-10-01). (a) The update overshoots when lowering the moment, (b) the pair is fixed to sites 1, 2 and blocks 1, 2,
  (c) every start of lmf updates the field and rewrites `mmagfield.aftest` (also in `job_band`)
- **[not started] Files in `SRC/` that seem unused** (2026-10-02). `SRC/.#Memo4rotation`, and in `SRC/exec/`: `cmake_install.cmake`, `hello.py`, `platform`,
  `lmf2.py`, `lmchk.py`, `pylmfa`, `pysample`, `ohtaka`, `epsPPd`, `epsPPsaito`, `job_senefbz`, `readeps_dig2.py`, `auto_kauto.py`. Confirm they are unused and remove them
- **[not started] Marks in the code asking for a check**: `need check` on terms of the force (`frcbandsym`, `fh_rv` in `m_bndfp.f90`),
  and in `m_bandcal.f90` "refactor m_igv2x so that any rank can read the data of any k point"

## GW1500 and Materials Project {#gw1500}

GW1500 is a database of QSGW for about 1500 materials of the Materials Project ([DOSnpSupplement](https://github.com/tkotani/DOSnpSupplement)).

- **[not started] `auto_mpquery.py` does not work with the present Materials Project** (2026-10-01). The validation of `mp_api` stops on the new form of IDs. Upgrade `mp_api`, or read REST directly
- **[undecided] Materials with doubtful structures** (2026-10-01). MP entries that are a sublattice taken out were selected. Whether to remove `INVALID_STRUCTURE` (12) and `SUSPECT_STRUCTURE` (4) from the set (now only noted).
  Reject or mark them with `auto_mpquery.py` when selecting
- **[not started] Materials where the 8³ mesh misses band extrema** (51 materials). Graphite-like carbon is a semimetal but gets a gap. Determine the gap including the points on the path, or use a finer mesh
- **[not started] Check the materials where QSGW80 is smaller than LDA**: Mo₂O₆ (mp-796276), CdSnAs₂ (mp-3829, experiment 0.26 eV)
- **[undecided] QSGW80 of slab materials depends on the vacuum thickness** (2026-10-08). For the same PbS slab, LDA 1.728 / 1.738 eV and QSGW80 3.725 / 4.152 eV (vacuum 13.7 / 19.4 Å; screening from periodic images).
  Note it in the database, or rerun with the vacuum made equal
- **[not started] Review of materials where MLO does not pass** (2026-10-06). Look one by one at the materials of the "Later list" of the database: C₃N₄ mp-570572, Rb₂Sc₂O₄ mp-7650, H₂ mp-632291,
  SnS mp-554134, layered BN mp-685145, H8 mp-730101, etc. To try: p in the ES basis, ES added in the gaps, the H basis (s,p of EH2).
  Use the ES threshold of 2.0 a.u. for layered crystals too (tell molecular from layered by the dimension of the bond connectivity)
- **[not started] With `wcsmear = false` and tf32, `hgw` falls in the middle of Σc** (2026-10-02, MgO, no message). `wcsmear = false` is not the default, but is used for comparisons

## Tests {#tests}

- **[not started] Make the test of ES on slabs a sample** (2026-10-08). For a PbS slab (mp-561320), the way ES are placed made a clear difference between MLO passing and failing.
  QSGW80 takes 5 hours, so look for a form that shows the same difference with LDA or a small slab

## Documents (ecaljdoc) {#docs}

- **[undecided] Red-pen marks kept**: 2 places in [MLO](./mlo) ("0.25 Ry" in the loss history, the note on "12 systems"). Not fully checked against the code or data
- **[not started] Equations running out of the text** (2026-11-14). At a width of 1300px kBT, pmt, mlo and gw; at 960px also basic, basic2 and mlo_gwsc. Break them with `aligned`
- **[not started] The menu heading "GPU version of QSGW" in Japanese**
- **[not started] What is left of the duplicated explanations** (2026-11-14). The main document (README_tutorial) and other pages (the `save` lines, lmchk, the SOC steps, etc.), duplication within kBT and mlo,
  old statements in the READMEs of Samples
