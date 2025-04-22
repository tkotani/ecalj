---
project: ecalj README at https://github.com/tkotani/ecalj
author: takao kotani
email: takaokotani@gmail.com
---
LICENCE: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.html)
# ecalj: a suite of first-principles elecronic structure calculations

## What can we do in ecalj package?

1. **All electron full-potential PMT method**
   
   The PMT method means; a mixed basis method of two kinds of augmented waves, that is, APW+MTO.
   In other words, the PMT method= the linearized (APW+MTO) method, which is unique except the [Questaal](https://www.questaal.org/) having the same origin with ecalj. Our recent research shows that very localized MTOs (damping factor exp(-kappa r) where kappa \sim 1 a.u), together with APW (cutoff is \approx 3 Ry) works well to get reasonable convergences. We can perform atomic-position rexalxation at GGA/LDA level. Because of including APWs, we can describe the scattering states very well.
   
   The current PMT formulation is given in

   [1][KotaniKinoAkai2015, PMT formalism](Document/PAPERandPRESENTATION/KotaniKinoAkai2015FormulationPMT.pdf)   
   [2][KotaniKino2013, PMT applied to diatomic molecules](Document/PAPERandPRESENTATION/KotaniKino2013PMTMolecule.pdf).

   Since we have automatic settings for basis set paremeters, 
   we don't need to be bothered with the parameter settings. Just crystal structure (cif,POSCAR) are needed for calculation. 

   In principle, it is possible to perform reasonable calculations just from crystal structures and
   very minimum setting. We have documents included about how to confirm calculations to be publicaiton quality. (but someghing still in progress..)

2. **the PMT-QSGW method** 
   
   The PMT-QSGW means 
   *the Quasiparticle self-consistent GW method (QSGW) based on the PMT method*.
   We can perform QSGW with litte settings by hand (for nonmag case, virtually no settings by hand).
   After converged, we can easily make band plots of the QSGW calculaiton without the Wanneir interpolation. This is because an interpolation scheme of self-energy is internally built in.
   We can handle even metals, Fermi surface as well. This may be an unique feature differenciate ecalj GW from those implemeneted in other packages.

   Paralellization is minimum and we now working on it on GPU. Roughly speaking, our current PMT-QSGW is applicable upto ~16 atoms a cell (light atoms are easier to handle) for research purpose.

   [3][Kotani2014,Formulation of PMT-QSGW method](Document/PAPERandPRESENTATION/Kotani2014QSGWinPMT.pdf)

   [4][PMT-QSGW applied to a variety of insulators](/Document/PAPERandPRESENTATION/deguchi2016.pdf)

3. **The Model Hamiltonian with Wannier functions** 
   
   We can generate the effective model (Maxloc Wannier and effective interaction between Wannier funcitons). 
   This is originally from codes by Dr.Miyake, Dr.Sakuma, and Dr.Kino. The cRPA given by Juelich group is implemented. We are now replacing this with a new version MLO (Muffin-Tin-orbail-based localized orbital).
   
4. **Dielectric functions and magnetic susceptibilities**
   
   We can calculate GW-related quantities such as **dielectric functions, spectrum
    function** of the Green's functions and so on (some of functions are obsolate. Need fixing).**Magnetic fluctuation** and so on are implemented but not documented well (ask us).
    See [papers by TK](https://scholar.google.co.jp/citations?user=xJV_H4YAAAAJ&hl=ja).

5. **Utilitys are included**
    
   For example, a converter between POSCAR(VASP) and our crystal strucrue file 'ctrls.*' are included.
   A command `job_mp [num of cores] [mp id]` access to material project to get crystal structure, and perform LDA/GGA calculaitons, showing DOS and BAND automatically.
   As it is combined with seekpath and spglib, 
   Brillowin Zone and symmetry lines are automatically drawn with matplotlib (we need minor fixing...).
   
## [ecaljdoc](https://ecalj.github.io/ecaljdoc/)
This is out main document sites for ecalj.
Here is [Install and InstallTest](https://ecalj.github.io/ecaljdoc/install/install).


