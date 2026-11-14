# Removed code (past code) {#past-code}

> English version of the page [manual/past_code](../../manual/past_code.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/past_code.md of the tree of 31d8516dc (source blob c53f7d78e0c5).

A list of code that was used in the past and is not used any more. It has been removed from the repository but stays in the
git history: take it out from the "commit before removal" in the table below. Look here also for code that may be used again
(prototypes and the like).

```bash
git show <commit>:<path>                 # show one file
git checkout <commit> -- <path>          # bring one file back into the working tree
git show <commit> --stat                 # what changed in that commit
```

Do this in a local clone of ecalj. The more detailed developer record (what was decided and why) is table 1 of `MD/past_log.md`
in the repository (not on this site).

*Table 1*. Removed code (newest first)

| removed | what | where it was | commit before removal |
| --- | --- | --- | --- |
| 2026-11-14 | Five modules used by nothing: `m_bz_integ` (BZ integral by the average over the tetrahedron corners), `m_pbindex` (index table of the product basis), `m_readqgcou` (reads QGcou), `m_HamRsMLO` (reads the old `HamRsMLO`), `m_qdwh_eigensolver` (prototype GPU eigensolver: generalized eigenvalue problem by the QDWH polar decomposition, 2026-03; may be used again) | `SRC/subroutines/m_bz_intg.f90`, `m_pbindex.f90`, `m_readQGcou.f90`, `m_qdwh_eigensolver.f90`, `module m_HamRsMLO` at the end of `m_HamPMT.f90` | `31d8516dc` |
| 2026-11-14 | The input routine `m_keyvalue` (`getkeyvalue`) that read GWinput; unused since the TOML input of 2026-05 | `SRC/subroutines/keyvalue.f90` | `31d8516dc` |
| 2026-10-07 | Three scripts of the input flow before 2026-05 | `SRC/exec/ctrlgenM1.py`, `auto_creplot.py`, `auto_job_mp.py` | `72967dec5` |
| 2026-10-02 | The Wannier-function (maximal localization) path: `hmaxloc`, `hpsig`, `wanplot`, `hsocmat`, the Wannier magnons `hmagnon`, the scripts `genMLWF`, `job_magnon`, etc. Replaced by MLO | `SRC/main/`, `SRC/subroutines/` (`main_wannier_*`, `m_wannier_maxloc`, ...), `SRC/exec/` | tag `last-wannier` (`dbcd6e51d`); there `Samples/WannierVsMLO/run.sh` reproduces the comparison with MLO |
| 2026-10-02 | AHC (anomalous Hall conductivity): `hahc`, `x0kf_ahc`, `job_AHC`, `hx0ahc.py` | `SRC/main/`, `SRC/subroutines/`, `SRC/exec/` | tag `last-wannier` |
| 2026-10-02 | lmfham2 (the old iterative MLO): `lmfham2`, `genMLO` | `SRC/main/`, `SRC/subroutines/`, `SRC/exec/` | tag `last-wannier` |
| 2026-10-02 | The symmetry routines `m_mksym_util.f90` (contents moved to `m_symop_util`, `m_symderive`, `m_symfind`, `m_mksym`) | `SRC/subroutines/m_mksym_util.f90` | `b889a00b8` |
| 2026-10-01 | Retired scripts | `SRC/exec_legacy/` | `a1c30d1a1` |
| 2026-10-01 | Prototype of automatic job runs | `jobauto/` | `a1c30d1a1` |
| 2026-10-01 | Old tools of TOOLS (about 60 items) | `TOOLS/` | `c6a1b6b7f` |
