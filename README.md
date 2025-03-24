---
project: ecalj README at https://github.com/tkotani/ecalj
author: takao kotani
email: takaokotani@gmail.com
---
LICENCE: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.html)
# ecalj: a suite of first-principles elecronic structure calculations

## What can we do in ecalj package? (on github)

1. **All electron full-potential PMT method**
   
   Here the PMT method means; a mixed basis method of two kinds of augmented waves, that is, APW+MTO.
   In other words, the PMT method= the linearized (APW+MTO) method, which is unique except the [Questaal](https://www.questaal.org/) having the same origin with ecalj. Our recent research shows that very localized MTOs (damping factor exp(-kappa r) where kappa \sim 1 a.u), together with APW (cutoff is \approx 3 Ry) works well to get reasonable convergences. We can perform atomic-position rexalxation at GGA/LDA level. Because of including APWs, we can describe the scattering states very well.
   
   The current PMT formulation is given in

   [1][KotaniKinoAkai2015, PMT formalism](Document/PAPERandPRESENTATION/KotaniKinoAkai2015FormulationPMT.pdf)   
   [2][KotaniKino2013, PMT applied to diatomic molecules](Document/PAPERandPRESENTATION/KotaniKino2013PMTMolecule.pdf).

   Since we have automatic settings for basis set paremeters, 
   we don't need to be bothered with the parameter settings. Just crystal structure (cif,POSCAR) are needed for calculation. (In 2022, we implemented a command, job_mp, to perform automatic calculations only from the id of material project. job_mp seems working well, although we need some further developments.)

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
   
## Install and automatic test
  After we get package from github, together with some standard tools, we do execute an automatic installer and test system. This is just by a command. Tests are for Linux with gfortran and ifort. See [README_Install.org](README_Install.org).

## Manuals, Documents, Samples for ecalj 
We have manuals and presentations of ecalj at [.Document/](./Document). Especially, [ecaljmanual.pdf](./Document/Manual/ecaljmanual.pdf) is the main manual. Read "LDA/GGA calculations and Plots" at first (but some part become obsolate---I have to revise)
QSGW calculations can be performed easily after you learend the procedure. Band plot and so on are in the same manner of the LDA calculations, althogh QSGW is computationally expensive (roughly 100~1000 times expensive. We ara going to use GPU now.).

It is instractive to reproduce samples in [./Document/PAPERandPRESENTATION/deguchi2016.pdf](). We can set up templates for your calculations. Ask us.

* How to perform paper quality calculations with minimum costs? See [./README_hints.org]()
* Usage minimum. QSGW for Si
See Section.4. of [./Document/Manual/ecaljmanual.pdf]().


## Very minimum to run ecalj

In Japanese, see a page by Dr.Gomi at http://gomisai.blog75.fc2.com/blog-entry-675.html (and others. Use search engine.)

1. write structure file ctrls.si by hand 
    (you can generate ctrls from POSCAR(VASP) with vasp2ctrl in
    ecalj/StructureTool/, thus cif --> POSCAR ---> ctrls is also possible.)

2. conver ctrls.si to ctrl.si by ctrlgenM1.py si 
   (without argument, it gives a help). 
   Then you have default ctrl.si (rename ctrlgenM1.ctr.si to ctrl.si). 
   Edit number of k points, spin (nsp=0 or 1) and so on if necessary.

3. Run "lmfa si" to prepare atoms (very quick).Then run 'mpirun -np 4 lmf-MPIK si'.
    This generates rst.si, which contains self-consistent density in LDA.
    Postprocessing for energy bands are job_band si, job_tdos, job_pdos are also available.
    For job_band, you need symmetry line file syml.si, which can be generated at the method implemented in GetSyml/

>>> NOTE: If you like to skip steps (1)-(3),  run ./job_materials.py Si at ./MATERIALS/. Then
>>>```bash
>>> >cd Si
>>> >cp ../syml.si
>>> >job_band si
>>>```
>>> This shows energy bands in LDA in gnuplot. To generate syml.si, we can use
>>> ecalj/GetSyml/getsyml.py. When it is correctly installed (see below),  
>>> ```bash
>>> >getsyml si
>>> ```
>>> should generate a syml.si from ctrl.si. You can edit it and run job_band.

4. For PMT-QSGW, make GWinput.tmp by mkGWIN_v2 si. Copy GWinput.tmp as GWinput.

5. Then run a script gwsc, e.g. "gwsc 2 si -np 3" 
    (2+1 iteration with 3 nodes).

6. To continue calculation do "gwsc 5 si -np 3" again.
   To start, you need ctrl.si rst.si QGpsi sigm.si. Latest gwsc is given in python3 (adopted from that by Dr.suzuki)

7. For band, dos, and pdos plot, 
    we have scripts which almost automatically makes these plot in
    gnuplot. Thus easy to modify these plots at your desposal.
    For example, job_band is for band plot. But symmetry line path file syml.si is required.
    The syml can be generated by getsyml.py, which also visualise the pathes in the BZ.


## Memo for StructureTool/ and Getsyml/
In any calculations, we first need to supply crystal structure correctly.
In the case of ecalj, we write it ctrls.*. 
All calculaitons can be performed from it.

For this purpose, we have converters between POSCAR
(VASP's crystal structure file, Cartesian setting is needed; 'conversion bug for Fractional aug2019') 
and ctrls.*(that for ecalj). In addition, we have a simple script to invoke crystal strucrure
viewer, usually VESTA. It is in [[file:StructureTool/README.txt][StructureTool/]].

Furthermore, we have a tool to generate BZ and symmetry lines on it for
band plot in [./GetSyml/]()
The symmetry line is written into syml.* and used for the
band plot mode, job_band. The BZ and the lines are visualized.

#### Install the viever at StructureTool/
-----
Here we use VESTA at http://jp-minerals.org/vesta/.
Download it, and expand it to a directory. 
VESTA can handle kinds of format of crystal structure.

Then make a softlike by
>  ln -s ~/ecalj/StructureTool/viewvesta.py ~/bin/viewvesta  
>  ln -s ~/ecalj/StructureTool/ctrl2vasp.py ~/bin/ctrl2vasp  
>  ln -s ~/ecalj/StructureTool/vasp2ctrl.py ~/bin/vasp2ctrl  
 
With this procedure we can run command viewvesta, ctrl2vasp,
vasp2ctrl from console as long as you have ~/bin/ in the command
search path. In my case, .bashrc have a line
  export PATH=$HOME/bin:$HOME/VESTA-x86_64:$PATH  

It depends on your machine. (after editing .bashrc, you have to do
"source ~/.bashrc" to reflect changes).

Set the variable of VESTA=, at the begining of 
~/ecalj/StructureTool/viewvesta.py to let it know where is VESTA.

#### Symmetry line finder at GetSyml/
-----
This is to generate symmetry lines. syml.* from ctrl.* in ecalj/GetSyml/
In the directory, we have getsyml.py, which is based on the
[seekpath](https://github.com/giovannipizzi/seekpath/)
on top of [spglib](https://github.com/spglib/spglib).
See Lincence.txt in it. Folllowing citations are required.
 1.Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, 
    Band structure diagram paths based on crystallography, Comp. Mat. Sci. 128, 140 (2017) 
 2.You should also cite spglib that is an essential library used in the implementation.

## How to do version up? minimum for git
Be careful to do version up. It may cause another problem.
But it is not so difficult to move it back to original version if you have knowledge of git.
An important things is keeping your changes by yourself.

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
    This takes all new changes.


I think it is recommended to use 
>gitk --all 

and read this document. Difference can be easily taken,
e.g. by 
>git diff d2281:README 81d27:README 
(here d2281 and 81d27 are several digits of the begining of its version id). 

>git show 81d27:README 

is also useful.  

## Licence 
- [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.html)
- For publications, we hope to make a citation cleary to this homepage such as;
  
  [foobar] ecalj package available from https://github.com/tkotani/ecalj/.




<!--- In addition, we are going to organize [[file:Document/Manual/!CategoryAndToken.org][Document/

Manual/  
: # CategoryAndToken.org]] which explain switches in ctrl file
The ctrl file have grammer of tree-type of ordering, "Category->token->subtoken".
Varieties of samples are at [[file:MATERIALS/][MATERIALS/]] (some of them are obsolate), 
see [[file:./MATERIALS/README_MATERIALS.org][./MATERIALS/README_MATERIALS.org]]
[[./Document/Manual/CategoryAndToken.org]]
(link to titus/ is broken now).
*** MEMO
** For 4f, we need modification to GWinput.tmp
   Old memo is at [[./Document/Manual/GdQSGW4.pdf][./Document/Manual/GdQSGW4.pdf]]
   Latest version automatically set default for 4f systems.
** (for previous users): known bug(or not) for spin susceptibility mode
(This mode is now obsolate because we are switching to a new method
with localized basis for spin susceptibility.)
T.Kotani thinks epsPP\_lmfh\_chipm branch may/(or may not) have a bug
(because of symmetrization). The bug may be near
#+begin_src f90
          if (is==nspinmx) then 
            symmetrize=.true.
            call x0kf_v4hz(npm,ncc,... 
#+end_src
in SRC/main/hx0fp0.m.F
(This bug may be from a few years ago, after I implemented EIBZ mode).
I think  "if (is==nspinmx.or.chipm) then" may be necessary
especially for cases with more than two atoms in the cell
(thus fe\_epsPP\_lmfh test may not work for this case...)
A possible test is by removing symmetrization -> use eibzsym=F. 

We have old documents at [[./Document/LMF@2009/]]
These are back up files at year2009. We still have some meaningful information in it.
But this is very detailed and mainly for developers.

See [[file:Document/Manual/GdQSGW4.pdf][Document/Mamual/GdQSGW4.pdf]]

* 4f system
Default setting might be not good enough. We are preparing a paper. 

* GaussianFilterX0.
This switch in GWinput is ver'y useful and promising (probably) 
to stabilize the convergence of metallic cases
(when many bands are located at the Fermi level).


--->

