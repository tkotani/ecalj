ecaljdoc/manual/UsageDetailed.md（サイトのページ「その他の説明」）から 2026-11-14 に開発者向けのメモとして外した節を、そのまま写したもの（利用者に要る事実はページの各節へ移した）。

### ecalj/Samples/
For the canonical layout and per-tree role table see the
[Samples overview page](../ecaljdoc/manual/samples.md) (or the upstream
[Samples/README.md](https://github.com/tkotani/ecalj/blob/main/Samples/README.md)).
Highlights touched in this page:

* `MLOsamples/` — new MTO Localized Orbital basis (Wannier replacement),
  with on-site W via `job_mloW`. **The method itself is written up in
  [MLO](../ecaljdoc/manual/mlo.md)** (what theta is, why `mlo_method = 4`, measured accuracy of the
  samples with band plots). See also the per-dir
  [README](https://github.com/tkotani/ecalj/blob/main/Samples/MLOsamples/README.md).
* the older Wannier90-style examples (CuMLWFs, La2CuO4, NiOMLWF, SrVO3MLWF; cRPA via the
  Juelich group's formulation) are not in the repository any more; `MLOsamples/` replaces them.
* effective mass calculation — `Samples/EffectiveMass/GaAs` (QSGW bands with SOC, the mass
  mode of `syml.<sname>`).

## How to do version up? git minimum

Be careful to do version up ecalj. It may cause another problem.
If you are not good at git, make another clone. 
Do not mix up with previous version (e.g. where you install)

>cd ecalj  
>git log  
   This shows what version you use now.

>git diff > gitdiff_backup
This is to save your changes added to the original (to a file git_diff_backup ) for safe.
I recommend you do take git diff >foobar as backup.
>git stash also move your changes to stash.

>git checkout -f
     CAUTION!!!: this delete your changes in ecalj/.
     This recover files controlled by git to the original which was just downloaded.

>git pull
    This takes all new changes. It is safer to use `git fetch` and `gitk --all` (`git checkout FETCH_HEAD -b youbranch`) to checkout your local branch.


I think it is recommended to use 
>gitk --all 

and read this document. Difference can be easily taken,
e.g. by 
>git diff d2281:README 81d27:README 
(here d2281 and 81d27 are several digits of the begining of its version id). 

>git show 81d27:README 

is also useful.  

## Bandplot for FSMOMMETHOD/=0
  Even when you use `[bz].fsmommethod` /= 0 in `ctrlg.<sname>.toml` for `gwsc`,
  you need to set `fsmommethod = 0` (or remove the key) when you run `job_band`.
  [If you run `job_band` with `fsmommethod` /= 0, it make a shift (adding bias magnetic field).]
  Anyway confusing. If necessary, I have to make it straight.

## We use only CORE1 treatment only (exchange only core)
See 10.1103/PhysRevB.76.165106 (Eq.35 and after).
Now I usually not use CORE2 (CORE1 only).

## a little unstable when metal GGA, especially when we have large empty regions. (still problematic?)
(negative density points appears: console output). I am not so sure about this is serious or not.

## Do we use VWN or GGA for QSGW? ===
In principle, QSGW results should not depend on VWN or GGA.  But there is minor dependence, because
1. frozen core density.
2. core eigenfunctions.
3. radial basis functions
4. Slight numerical reason
(This is probably because Sigma-interpolation procedure But not exactly figured out yet →
affect about 0.02eV as for band gap for GaAs.). In anyway, use VWN (`[ham] xcfun = 1`)
as standard. And such technical things affects, 0.05 eV level of error for band gap.

## numerical accuaray of QSGW band energies
As for the band gaps 0.1 eV is a kind of limiting accuaray for usual semiconductors with 0~3 eV band gap materials.
For Si or something simple, you may be able to have better accuracy. However, QSGW is on the so many parameters.
So, just say ~0.1eV level of accuracy. Convergence check is expected for the technical developments in future.

## Anisotropic Q0P 
   In cases, we are worrying about the anisotropy. The offset Gamma method uses the 1/10 of gamma cell (`[gw] Q0Pchoice = 1`; see [gwinput](../ecaljdoc/manual/gwinput.md)).
   It seems not so bad, however, we may need to examine the anisotropy problem again in future.  
   It is problematic to use unbalanced k points for anisotropic cell. See Copmuter Physics Comm. 176(2007)1-13.
