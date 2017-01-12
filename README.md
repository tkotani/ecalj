ecalj package
=============================
This is read me at https://github.com/tkotani/ecalj. 
A first-principle electronic structure calculation package in
f90, especially for the PMT-QSGW. 

Tutorial course: we have course at CMD workshops held by Osaka university (every
March and Sep). http://phoenix.mp.es.osaka-u.ac.jp/CMD/index_en.html

We have another home page at http://pmt.sakura.ne.jp/wiki/, but
not well-organized yet, little in English yet. We will renew it.
 
Overview
--------------------------
1.  All electron full-potential PMT method: PMT= a mixed basis method of two
   kinds of augmented waves, that is, L(APW+MTO). 
   Relaxiation of atomic positions is possible in GGA/LDA and LDA+U.
   Our recent development shows that very localized MTO (damping
   factor is \sim 1 a.u), together with APW
   (cutoff is \sim 2 to 4 Ry) works well to get reasonable convergences.
   In principle, it is possible to perform default calculations just
   from atomic structures.
   http://journals.jps.jp/doi/abs/10.7566/JPSJ.83.094711
   
2. The PMT-QSGW method, that is,
   the Quasiparticle self-consistent GW method (QSGW) based on the PMT method. 
   In addion, we can calculate dielectric functions, 
   spectrum function of the Green's functions and so on.
   GW-related codes are in ~/ecalj/fpgw/.
   For paralellized calculations, 
   we can use lmf-MPIK and mpi version of hvccfp0,hx0fp0_sc,hsfp0_sc.
   (although we still have so much room to improve it).
   The PMT allows us to perform
   the QSGW calculations virtually automatically.
   http://journals.jps.jp/doi/abs/10.7566/JPSJ.83.094711

3.  Wannier function generator and effective model generator
   (Maxloc Wannier and effective interaction between Wannier funcitons). 
   This is adopted from codes by Dr.Miyake,Dr.Sakuma, and Dr.Kino.
   See fpgw/Wannier/README.

Utilities such as a converter between POSCAR(VASP) and our crystal strucrue file
'ctrls.*' are included (slightly buggy).


!! CAUTION for know bug(or not) for spin susceptibility mode!!! (apr2105).
T.Kotani thinks epsPP\_lmfh\_chipm branch may/(or may not) have a bug
(because of symmetrization). It may be near
------------------
          if (is==nspinmx) then 
            symmetrize=.true.
            call x0kf_v4hz(npm,ncc,... 
-----------------
in fpgw/main/hx0fp0.m.F
(This bug may be from a few years ago, after I implemented EIBZ mode).
I think  "if (is==nspinmx.or.chipm) then" may be necessary
especially for cases with more than two atoms in the cell
(thus fe_epsPP_lmfh test may not work for this case...)
A possible test is remove symmetrization---> use eibzsym=F.
If necessary, let me know...


Requirement for using ecalj
--------------------------------
For your publications, make a citation directly to this homepage as;

[1] ecalj package at https://github.com/tkotani/ecalj/. 
Its one-body part is developed based on the Questaal at 
https://lordcephei.github.io/. Its GW part is originally taken from ecalj.

in the references on the same footing of other papers.


Install and Test 
--------------------------------
Follow these steps explained below.  
However, You can run steps (1)-(2) by a command InstallAll.foobar at ecalj/.
When install procedure have finished, we have all required binaries and
shell scripts in your \verb+~/bin/+ directory).
(or somewhere else where BINDIR specified in InstallAll.*).

(0) Get ecalj package and get tools.  
(1) make for single-core LDA part, make for MPIK LDA part,   
    make for MPIK GW part,
    set crystal structure tools.
(2) Install test  

In the following explanation, we assume gfortran in ubuntu.
But we can also use ifort and others in your environment with
minimum changes in makefiles. Let me have your InstallAll.foobar; it is very helpful for us.
For for small systems such as Si and GaAs, 
we can use even Ubuntu + gfortran + note PC for test purpose to observe how QSGW works.

#### (0) Get ecalj package and get tools.
--- Let us assume you use ubuntu. ---
You need following tools and libraries to be installed.  
>sudo apt-get install git  #version control and to get source from github  
>sudo apt-get install gitk #git gui.   
>sudo apt-get install gfortran      # GFORTRAN  
>sudo apt-get install openmpi-bin libopenmpi-dev #  or openmpi-dev ?
>sudo apt-get install libfftw3-3     or something else # FFTW  
>sudo apt-get install libblas3gf     or something else # BLAS  
>sudo apt-get install liblapack3gf   or something else  # LAPACK  
>sudo apt-get install csh bash tcsh gawk  # shells  

memo: I think etags is automatically installed when you install 
      emacs in the latest ubuntu. etags is essentially needed only for developers.

Note that python 2.x is also assumed 
(usually already in ubuntu. Type \>python (ctrl+D for quit)).
Version ctrl is by git (which makes things easier, but not
necessarily required for installation).

After you have installed git (version control software), 
you can get ecalj package by  
>git clone https://github.com/tkotani/ecalj.git # Get source code  

for ecalj. or get it in the form *.zip 
from the page https://github.com/tkotani/ecalj
(push zip button). I recommend you to use git, 
to check your changes (\>git diff), know version id, and upgrade.
After you did the above git clone command, a directory ecalj/ appears
(under the directory at which you did git clone).

We can check history of ecalj code development by
"\>gik --all" at ecalj/ directory after you did git clone.

#### (2) is performed by InstallAll.foobar.
Run InstallAll.ifort (in the case of intel fortran) at ecalj/. 
We can invoke this command as;

>cd ecalj  
>./InstallAll.ifort  
(To clean all, do ./CleanAll.ifort).  

Please look into the script "InstallAll.ifort". It is a small text file.
It contains the setting of your BINDIR= directory,
to which the InstallAll.ifort will copy all binaries and scripts.
It internally uses three machine-comilar dependent files;  
  a.lm7K/MAKEINC/Make.inc.ifort (for single core version )  
  b.lm7K/MAKEINC/Make.inc.ifort_mpik (k-point paralell version )  
  c.fpgw/exec/make.inc.ifort  (this is only for mpi-omp version ).  
At the last stage of the script, it runs automatic tests.
(You can neglect failure for nio_gwsc; it may show one-failure among two checks).
The test may use ten minutes or more... Have a coffee!
  
InstallAll.ifort may not work for your environment.
The you may prepare your own InstallAll.foobar,
in which you have to set compilar, linker, compilar options.

When InstallAll.ifort works well, it will show OK! signs finally.
(one last test (nio_gwsc) may fail in cases, but usually no problem).

<B> you don't need to read follwings when InstallAll.foobar works fine </B>
##### (1) make single core LDA part (it is in ecalj/lm7K/).
Let us assume gfortran case.
Move to ecalj/lm7K/, then do "make PLATFORM=gfortran LIBMATH=xxx". 
Then make process start. (LIBMATH= specify BLAS,Lapack, and fftw.)
The main makefile is ecalj/lm7K/Makefile, which contains lines
>  PLATFORM=gfortran   #default is PLATFORM=gfortran  
>  ...  
>  include MAKEINC/Make.inc.$(PLATFORM)  

This means that this Makefile uses ecalj/lm7K/MAKEINC/Make.inc.gfortran
as a part of the Makefile. Thus we have to check settings in it 
in advance to run "make PLATFORM=...".
LIBMATH= should contain path to the math libraries, FFTW, BLAS and LAPACK.
An example is   
  LIBMATH="/usr/lib/libfftw3.so.3 /usr/lib/liblapack.so.3gf
  /usr/lib/libblas.so.3gf"  
Compilar options FFLAGS=, FFLAGS_LESS=... can be choosed by your own
manner. But usually Make.inc.gfortran works without changes
(let me know your changes; I like to include it in ecalj).

Let us think about an ifort case.
In this case, we run  
>make PLATFORM=ifort LIBMATH='-mkl'   

There are several MAKEINC/Make.inc.ifort*
(not _mpik*) with which we installed to machines. 
You can choose one of them or you can set your own Make.inc.ifort.*
(compilar, compilar options, math library).

Warning messages like ": warning: ignoring old commands for target `/vxcnls.o'" is
just because of problem of Makefile. you can neglect this. We will fix it..

Parallel make like  
>make -j24 PLATFORM=gfortran  

can speed it up for multicore machines(24 core in this case). 
But it stops because of dependency is not well-described in our current Makefile. 
In such case, repeat it a few times, or repeat it without -j24.

Finally run  
>make PLATFORM=gfortran install  

This just copy required files (binaries and scripts) to your ~/bin.
(check it in Makefile). If you like to copy them to ~/xxx instead of~/bin,
make with BINDIR=xxx.

(For CMD workshop participants: run  
>make PLATFORM=ifort.cmd LIBMATH='-mkl' BINDIR=~/bin


#### WARN! Install problems ---
* I saw that current ecalj with gfortran4.6 or 4.7 works fine with
  FFLAGS=-O2, but failed with FFLAGS=-O3. (I am not sure now).
* ifort12 may need FFLAGS=-O1 in MAKEINC/Make.inc.ifort. 
  -->Try InstallAll.ifort12.
* We may need -heap-arrays 100 (when zero, we had a problem in a version
  of ifort). In cases, -heap-arrays option did not generate working binaries.
  However, I think "ulimit -s unlimited" before QSGW calculations and
  so on works OK. So, maybe we don't need -heap-arrays option.
* mpiifort works for liker, but mpif90 did not... (but opposite case
  maybe). Need to set this in lm7K/MAKEINC/Make.inc.ifort
  lm7K/MAKEINC/Make.inc.ifort_mpik and fpgw/exec/make.inc.ifort
  (FC and LK variables).


##### (2)-i make MPI LDA part.
lmf-MPIK and lmfgw-MPIK are needed for gwsc (srcipt for QSGW). 
These are k-point parallel version of lmf, and gw driver lmfgw. To
make it, do  
"make PLATFORM=gfortran_mpik".  
For ifort, set PLATFORM=ifort_mpik.  
Then Makefile includes ecalj/lm7K/MAKEINC/Make.inc.ifort_mpik.
You may need to add -heap-arrays 1 (for large calculations. Because we
use large stacksize) to ecalj/lm7K/MAKEINC/Make.inc.ifort_mpi, but I
am not so sure about this.

(For CMD workshop participants: run  
 >make PLATFORM=ifort_mpik.cmd LIBMATH='-mkl'

which corresponds to MAKEINC/Make.inc.ifort_mpik.cmd)

*Clean up:  
If something wrong. do "make clean" or "make cleanall" and start over.
Look into Makefile if you like to know what they do.
"make cleanall" removes all *.o *.a modules, and binaries.

* Move binaries to your bin by 
>make install
at ecalj/lm7K. It just moves all requied binaries to your ~/bin.
In advance, you have to make your bin, and add ~/bin to  your path
(e.g. "export PATH=$HOME/bin:$PATH" in .bashrc . Then login again or "source .bashrc")

##### (2)-ii Installation for fpgw/
At ecalj/fpgw/exec/ directory, you have to a softlink make.inc such as
>lrwxrwxrwx 1 takao takao 17 Aug 25 13:18 make.inc -> make.inc.gfortran

For each machine you have to prepare your own make.inc.foobar 
(There are samples. Here is the case of make.inc.ifort.cmd), 
and do  
>ln -s make.inc.ifort.cmd make.inc  

to make a soft like make.inc -> make.inc.cmd


*** Q. What is soft link foo -> bar?  A. "foo" is an alias of the file "bar"  ***

Then you have to run  
>make  
>make install  
>make install2  

Before this, you have to set blas and lapack in fpge/exec/make.inc.
(for ifort, -mkl is enough. LIBMATH= should be the same as that in Make.inc.*.
"make install" copy requied files to your ~/bin.

* Caution!: we often see "Segmentation fault"due to stacksize limit 
(See the size by a command "ulimit -a"). 
It is needed to run "ulimit -s unimited" in the job-submition script 
or before running GW jobs. 


##### (2)-iii Install test
We have to check whether binaries works fine or not.
Move to ecalj/TestInstall. Then type make (with no arguments). 
It shows help about how to do test.
To test all of binaries, just do
>make all  

All tests may require ~10min or a little more.  (nio_gwsc takes ~300sec)
In cases, nio_gwsc fails, showing  
 >FAILED: nio_gwsc QPU compared by ./bin/dqpu  
 >PASSED: nio_gwsc diffnum  

However, we do not need to care its failure sign. (so nio_gwsc test
must be improved...). (numerically small differences).

Help of make (no arguments) at ecalj/TestInstall, shows
>make lmall   !tests only LDA part.  
>make gwall   !tests only GW part.  

* NOTE (nov19 2014 kino):   
In TestInstall/Makefile.define, we define  
LMF=lmf  
LMFP=lmf-MPIK  
(it is possible to use "LMFP=lmf-MPI" instead(for future development).
If we set LMFP=$(LMF), tests are done with using lmf, not with using lmf-MPIK.

* NOTE: in principle, repeat make should do nothing when all binaries
are correctly generated. However, because of some problem in makefile, 
you may see some make procedure is repeated. You can neglect it as
long as "All test are passed!" is shown in the (4)install test.

##### (2)-iv Structure tool.
This is not necessary if you don't need to need converter between PROCAR and ctrl/ctrls
(crystal structure file in ecalj).
In any calculations, we first have to supply crystal structure correctly.
To help this, we have a converter between POSCAR(VASP's crystal
structure file) and ctrls(that for ecalj). 
In addition, we have a simple tool to invoke crystal strucrure viewer.
It is in \verb+ecalj/Structure/tool/.

In advance, install a viever of crystal structure for POSCAR.
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



### How to do version up? ###

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


### Documents of ecalj ###
We have documents in ecalj/Document/
Especially, 
ecalj/Document/Manual/ecaljmanual.
is the main document.


###  Usage minimum. (e.g, PMT-QSGW(gwsc) for si) ###
Read ecalj/Document/Manual/ecaaljmanual.pdf
<pre>
Here is its very minimum.
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
This shows you band by LDA.

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
-------------------------------------------


### New: symmetry line finder (jan.2017) ###
Get symmetry line. syml.* from ctrl.* in ecalj/GetSyml/

In this directory, we have getsyml.py, which is based on the seekpath
https://github.com/giovannipizzi/seekpath/
See Lincence.txt
Folllowing citations are required.
    1.Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, 
       Band structure diagram paths based on crystallography, Comp. Mat. Sci. 128, 140 (2017) 
    2.You should also cite spglib that is an essential library used in the implementation.

</pre>

History
<pre>
* March 2016: new histgram bin m_freq.F 
  (HistBin_ratio and HisBin_dw are used to specify new mesh.

* March 2016:  wklm(1) is only used (only f_L for l=m=0 is used. 
  See Eq.28 in JPSJ83,094711(2014).)
</pre>


