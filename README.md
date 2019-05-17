ecalj package
=============================
This is read me at https://github.com/tkotani/ecalj. 
A first-principle electronic structure calculation package in f90 and
python, especially for the PMT-QSGW. Query to takaokotani@gmail.com.
 
Overview
--------------------------
1.  All electron full-potential PMT method: PMT= a mixed basis method of two
   kinds of augmented waves, that is, L(APW+MTO). 
   Relaxiation of atomic positions is possible in GGA/LDA.
   Our recent development shows that very localized MTO (damping
   factor is \sim 1 a.u), together with APW
   (cutoff is \sim 2 to 3 Ry) works well to get reasonable convergences.
   In principle, it is possible to perform default calculations just
   from atomic structures.
   http://journals.jps.jp/doi/abs/10.7566/JPSJ.83.094711
   
2. The PMT-QSGW method, that is,
   the Quasiparticle self-consistent GW method (QSGW) based on the PMT method. 
   In addion, we can calculate dielectric functions, 
   spectrum function of the Green's functions and so on.
   GW-related codes are in ~/ecalj/fpgw/.
   For paralellized calculations, 
   we can use lmf-MPIK in LDA level (although we still have so much room to improve it).
   The PMT allows us to perform
   the QSGW calculations virtually automatically.
   http://journals.jps.jp/doi/abs/10.7566/JPSJ.83.094711

3.  Wannier function generator and effective model generator
   (Maxloc Wannier and effective interaction between Wannier funcitons). 
   This is adopted from codes by Dr.Miyake,Dr.Sakuma, and Dr.Kino.
   See fpgw/Wannier/README.

Utilities such as a converter between POSCAR(VASP) and our crystal strucrue file
'ctrls.*' are included (slightly buggy).


####(for previous users): known bug(or not) for spin susceptibility mode
(This mode is now obsolate because we are switching to a new method
with localized basis for spin susceptibility.)
T.Kotani thinks epsPP\_lmfh\_chipm branch may/(or may not) have a bug
(because of symmetrization). The bug may be near

          if (is==nspinmx) then 
            symmetrize=.true.
            call x0kf_v4hz(npm,ncc,... 

in fpgw/main/hx0fp0.m.F
(This bug may be from a few years ago, after I implemented EIBZ mode).
I think  "if (is==nspinmx.or.chipm) then" may be necessary
especially for cases with more than two atoms in the cell
(thus fe\_epsPP\_lmfh test may not work for this case...)
A possible test is by removing symmetrization---> use eibzsym=F. 


Request for your publication
--------------------------------
For your publications, make a citation cleary to this homepage such as;

[1] ecalj package at https://github.com/tkotani/ecalj/. 
Its one-body part is developed based on the Questaal at 
https://lordcephei.github.io/. Its GW part is originally taken from ecalj.

in the references on the same footing of other papers.


Install and Test 
------------------------
See ecalj/README_Install.org


Manuals and Samples of ecalj 
--------
We have documents in ecalj/Document/
Especially,  ecalj/Document/Manual/ecaljmanual is the main document.
Samples are at MATERIALS/README.org


How to perform paper quality calculations with minimum costs?
--------
See ecalj/README_hints.org.


Structure tool
-----------------
In any calculations, we first supply crystal structure correctly.
In the case of ecalj, it is written in ctrls.*.
For this purpose, we have converters between POSCAR
(VASP's crystal structure file) and ctrls(that for ecalj). 
In addition, we have a simple script to invoke crystal strucrure
viewer, usually VESTA. It is in ecalj/Structure/tool/.
Further more, we have a tool to generate BZ and symmetry lines on it for
band plot. The symmetry line is written into syml.* and used for the
band plot mode, job_band.


Install the viever
----------------
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


Usage minimum. QSGW for si
------
See ecalj/Document/Manual/ecaaljmanual.pdf
Here I show its very minimum.
-------------------------------------------
(1) Write structure file ctrls.si by hand 
    (you can have ctrls from POSCAR(VASP) with vasp2ctrl in
    ecalj/StructureTool/.)

(2) conver ctrls.si to ctrl.si by ctrlgen2.py si --nk=6 
   (without argument, it shows help). 
   Then you have default ctrl.si (rename ctrlgen2.ctr.si to ctrl.si). 

(3) Run "lmfa si" to prepare atom.

NOTE: If you like to skip them,  run ./job_materials.py Si at /home/takao/ecalj/MATERIALS.
 >cd Si
 >cp ../syml.si
 >job_band_nspin1 si
This shows you band by LDA. To generate syml.si, we can use
ecalj/GetSyml/getsyml.py. When it is correctly installed (see below), 
$getsyml si
should generate a syml.si from ctrl.si. You can edit it and run job_band.

(4) For PMT-QSGW, make GWinput.tmp by mkGWIN_v2 si.
    Copy GWinput.tmp as GWinput. (you supply three numbers for the
    command mkGIWN_V2.)

(5) Then run a script gwsc, e.g. "gwsc 2 si -np 3" 
    (2+1 iteration with 3 nodes).

(6) To continue calculation do "gwsc 5 si -np 3" again.
    (To start, you need ctrl.si rst.si QGpsi ESEAVR sigm.si)
    When you start from these files, 0th iteration is skipped
   ---thus we have just five iteration.

(7) For band, dos, and pdos plot, 
    we have scripts which almost automatically makes these plot in
    gnuplot. Thus easy to modify these plots at your desposal.


Symmetry line finder.
----------
This is to generate symmetry lines. syml.* from ctrl.* in ecalj/GetSyml/
In the directory, we have getsyml.py, which is based on the seekpath
https://github.com/giovannipizzi/seekpath/
See Lincence.txt in it.
 Folllowing citations are required.
    1.Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, 
       Band structure diagram paths based on crystallography, Comp. Mat. Sci. 128, 140 (2017) 
    2.You should also cite spglib that is an essential library used in the implementation.


How to do version up? 
-----
Be careful to do version up. It may cause another problem.
But it is not so difficult to move it back to original version if you use git.
An important things is keeping your changes by yourself.
Especially your own Make.inc.* files (see InstalAll.ifort).

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
e.g. by >git diff d2281:README 81d27:README (here d2281 and 81d27 are
several digits of the begining of its version id). 
>git show 81d27:README is also useful.  



History (not maintained well).
---------
* March 2019: this document is cleaned up slightly
* March 2016: new histgram bin m_freq.F 
  (HistBin_ratio and HisBin_dw are used to specify new mesh.
* March 2016:  wklm(1) is only used (only f_L for l=m=0 is used. 
  See Eq.28 in JPSJ83,094711(2014).)


