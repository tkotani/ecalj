ecalj 
============
 A first-principle electronic structure calculation package in
 f90, especially for the PMT-QSGW. 

 We have a tutorial course at CMD workshops held by Osaka university (every
 March and Sep). http://phoenix.mp.es.osaka-u.ac.jp/CMD/index_en.html

 This is read me at https://github.com/tkotani/ecalj. 
 (checked at Dec2014.This README is still away from perfect...)

Overview
------------

1.  All electron full-potential PMT method: a mixed basis method of two
   kinds of augmented waves, that is, L(APW+MTO). 
   Relaxiation of atomic positions is possible in GGA/LDA and LDA+U.
   Our recent development shows that very localized MTO (damping
   factor is \sim 1 a.u), together with APW
   (cutoff is \sim 3 to 4 Ry) works well to get reasonable convergences.
   In principle, it is possible to perform default calculations just
   from atomic structures.
   
2. The PMT-QSGW method, that is,
   the Quasiparticle self-consistent GW method (QSGW) based on the PMT method. 
   In addion, we can calculate dielectric functions, spin fluctuations,
   spectrum function of the Green's functions and so on.
   GW-related codes are in ~/ecalj/fpgw/.
   For paralellized calculations, 
   we can use lmf-MPIK and mpi version of hvccfp0,hx0fp0_sc,hsfp0_sc.
   (although we still have so much room to improve it).
   The PMT allows us to perform
   the QSGW calculations virtually automatically
   (10.7566/JPSJ.83.094711). http://journals.jps.jp/doi/abs/10.7566/JPSJ.83.094711
 

3.  Wannier function generator and effective model generator
   (Maxloc Wannier and effective interaction between Wannier funcitons). 
   This is mainly from Dr.Miyake,Dr.Sakuma, and Dr.Kino.
   See fpgw/Wannier/README.

We have another home page at http://pmt.sakura.ne.jp/wiki/, but
not well-organized yet, little in English yet. We will renew it.

Utilities such as a converter between POSCAR(VASP) and our crystal strucrue file
'ctrls.*' are included. (slightly buggy; let T.Kotani know problems in
it; note that we should supply numerically accurate atomic positions to judge 
crystal symmetry automatically).

We need acknowledgment as
----------------------------------------------------------
For your publications, please make two citations directly 
to this homepage as;

[1] ecalj package at https://github.com/tkotani/ecalj/. 
Its one-body part is developed based on Ref.[2].  
[2] LMsuit package at http://www.lmsuite.org/. 
Its GW part is adopted mainly from Ref.[1].

in addition to our papers (in the references).


Install and test
-----------------------------
Follow these steps explained below.  
(0) Get ecalj package and get tools.  
 --- You can run folling (1)-(4) by a command InstallAll.foobar
at ecalj/  
(1) make for single-core LDA part,   
(2) make for MPIK LDA part,   
(3) make for MPIK GW part.  
(4) Install test  
(5) crystal structure tools  (not necessary).
In the following explanation, we assume gfortran in ubuntu. 
But we can also use ifort, and others in your environment with
minimum changes in makefiles.


##### (0) Get ecalj package and get tools.
--- Let us assume you use ubuntu. ---
You need following tools and libraries to be installed.  
>sudo apt-get install git  #version control and to get source from github  
>sudo apt-get install gitk #git gui.   
>sudo apt-get install gfortran      # GFORTRAN  
>sudo apt-get install openmpi-bin openmpi-dev     # MPI  
>sudo apt-get install libfftw3-3     or something else # FFTW  
>sudo apt-get install libblas3gf     or something else # BLAS  
>sudo apt-get install liblapack3gf   or something else  # LAPACK  
>sudo apt-get install etags csh bash tcsh gawk  # shells  

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

#### Automatic command from (1) through (4)  
Following procedures, (1)-(4), are done automatically by a script,
for example, InstallAll.ifort (in the case of intel fortran) as;

>cd ecalj  
>./InstallAll.ifort  
(To clean all, do ./CleanAll.ifort).  

Please look into the script "InstallAll.ifort". 
It contains the setting of your BINDIR= directory,
to which the InstallAll.ifort will copy all binaries and scripts.
It internally uses three machine-comilar dependent files;  
  a.lm7K/MAKEINC/Make.inc.ifort (for single core version -->(1))  
  b.lm7K/MAKEINC/Make.inc.ifort_mpik (k-point paralell version -->(2)  
  c.fpgw/exec/make.inc.ifort  (this is only for mpi-omp version -->(3)).  
At the last stage of the script, it runs automatic tests.
(You can neglect failure for nio_gwsc; it may show one-failure among two checks).
The test may use ten minutes or more... Have a coffee!
  
InstallAll.ifort may not work for your environment.
The you may prepare your own InstallAll.foobar,
in which you have to set compilar, linker, compilar options.

When InstallAll.ifort works well, it will show OK! sign finally.
(one last test (nio_gwsc) may fail in cases, but usually no problem).
The skip following (1) thru (4). (5) is also not necessary.


---
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

You may need to add -heap-arrays 100 (or 0)   
(Because we use large stacksize; but I am not so sure about this.
 In a case, -heap-arrays option did not generate working binaries.
 However, I think "ulimit -s unlimited" before QSGW calculations and
 so on works OK. So, maybe we don't need -heap-arrays option.)

Warning messages like ": warning: ignoring old commands for target `/vxcnls.o'" is
just because of problem of Makefile. you can neglect this. We will fix it..

I saw that current ecalj with gfortran4.6 or 4.7 works fine with
FFLAGS=-O2, but failed with FFLAGS=-O3. (I am not sure now).

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

which corresponds to MAKEINC/Make.inc.ifort.cmd
We did not set "-heap-array 100" option since it give not working
binaries.)  

##### (2) make MPI LDA part.
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

##### (3) Installation for fpgw/
At ecalj/fpgw/exec/ directory, you have to a softlink make.inc such as
>lrwxrwxrwx 1 takao takao 17 Aug 25 13:18 make.inc -> make.inc.gfortran

For each machine you have to prepare your own make.inc.foobar 
(There are samples. Here is the case of make.inc.ifort.cmd), 
and do  
>ln -s make.inc.ifort.cmd make.inc  

to make a soft like make.inc -> make.inc.cmd

* Q. What is soft link foo -> bar?  A. "foo" is an alias of the file "bar".  

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


##### (4) Install test
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


##### (5) Structure tool.
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


<pre>
###### Documents of ecalj #####################
We have three documents.

1. GetStarted:
ecalj/Document/Manual/EcaljGetstarted.pdf
Almost written OK.
It gives minimum instruction of ecalj.

2. Usage:
ecalj/Document/Manual/EcaljUsage.pdf
Not yet completed.

3. Details:
Two papers (but not yet submitted) are
  ecalj/Document/Manual/pmttheory17.pdf (for PMT)
  ecalj/Document/Manual/pmtqsgw13.pdf   (for PMT-QSGW)
A little extended versions in japanese are
ecalj/Document/Manual/ecalj_pmttheory17j.pdf (for PMT)
ecalj/Document/Manual/ecalj_pmqsgw13j.pdf   (for PMT-QSGW)
We will unify this...


#######  Usage minimum. (e.g, PMT-QSGW(gwsc) for si)  ###################
Read ecalj/Document/Manual/EcaljGetstarted.pdf
Here is very minimum.
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

###### PROCASR (VASP format) generator #########
PROCASR mode for lmf (not yet in lmf-MPIK)
Band weight decomposition. 
(Size of circles show the size of components. Superposed on band
plot). See ecalj/MATERIALS/MgO_PROCAR/README.
May19.2014



#####################################################
Usage problems, Q&A error message.

1.Bandplot for FSMOMMETHOD/=0 
Even when you use FSMMOMMETHOD/=0 in GWinput for gwsc, 
yuo need to set FSMOMMETHOD=0 (or comment it out) when you run job_band_nspin2.
[If you run job_band_nspin2 with FSMOMMETHOD/=0, it make a shift 
 (adding bias magnetic field).]

2.Note that ctrlgenM1.py automatically set this for --systype=molecule.
   Then we have 
       TETRA=0
       N=-1  #Negative is the Fermi distribution function W= gives temperature.
       W=0.001 #W=0.001 corresponds to T=157K as shown in console
   In addiiton, FSMOM (n_up-n_down) is needed (FSMOMMETHOD=1)if we
   have magnetic moment.

3. core>evalence message.
   Ecore is grater than Evalence.
   For save, we do not allow this.
   Complare ECORE file and valence levels, shown in log file or
   console output.

4. If you see a error message from lmf (e.g., internally called in the gwsc script).;
  Exit -1 rdsigm: Bloch sum deviates more than allowed tolerance (tol=5e-6)
You have to enlarge RSRNG so that lmf finsh normally.
See FAQ at the bottom of http://titus.phy.qub.ac.uk/packages/LMTO/gw.html
E.g. try RSRNGE=8 for Ni...

5. Back ground charge and fractional Z.
   You can use fractional numbers for ATOM_Z, and also can set
   valence charge by BZ_ZBAK (I removed BZ_VAL).
   You see console out put, e.g,
     "Charges:  valence    19.80000   cores     8.00000   nucleii   -28.00000
        hom background       .20000   deviation from neutrality:      0.00000
   . This is a case with BZ_ZBAK=.2.

   NOTE: at the first iteration, Charges: shows such as
     Charges:  valence     8.00000   cores    20.00000   nucleii   -28.00000
      hom background     0.12300   deviation from neutrality: 0.12300
      because of the initial condition by superposition of atoms. It show
      deviation seems nonzero. But charge should be conserved from the
      next iteration.

6. Not converged in metal. --->mixing may help
For example, if you try metal such as Bi2Sr2CuO6, it may fail at LDA/GGA level.
Then use ITER MIX=A2,b=.2. or something (.2 means it only mix 20% of output to give
new input for next iteration). Then I see convergence. (b is the
mixing parameter.

7. Use PZ or not.
   If spillout of core is not so small (more then 0.05 or something.),
   it is better to use PZ(lo). Treat the core as valecne.
   Bi4d is such a case. Maybe use PZ=0,0,4.9

8.Core treatment 
  See 10.1103/PhysRevB.76.165106 (Eq.35 and after).
  Now I usually not use CORE2 (CORE1 only).

9. ERROR EXIT! rgwina: 2nd wrong l valence
   This may be because you use wrong GWinput.
   Back it up. And run mkGWIN_lmf2 (any n1 n2 n3 is fine).

10. Known bug.
    Error occurs when system is anisotroic such as CuAlTe2. 
    Temporary fix is "Add token NPWPAD=100 in HAM category".
    (guess of used APW fails (more than expected)).
    CuAlTe2,CuGaTe2 cases.

11. Known bug
    a little unstable when metal GGA, especially when we have large
    empty regions.

==========================================
Wanneir Branch now under developing (imported from T.Miyake's Wannier code and H.Kino's).
   A. make at ecalj/fpgw/Wannier/ directory, and do make, and make install. 
      (need to check Makefile first). You first have to install fpgw/exec/ in advance.
   B. Samples are at these directories. 
      MATERIALS/CuMLWFs (small samples),
      MATERIALS/CuMLWF/
      MATERIALS/CuMLWFs/
      MATERIALS/FeMLWF/      
      MATERIALS/NiOMLWF/
      MATERIALS/SrVO3MLWF/
   C. With GWinput and ctrl.*, run 
      >genMLWF
      at these directories.
      In GWinput, we supply settings to generate Wannier funcitons. (Sorry,not documentet yet..)
   D. After genMLWF, do
      >grep Wan lwmatK*
      then compare these with Result.grepWanlwmatK
      These are onsite effective interactions (diagonal part only shown).
      *.xsf are for plotting the Maximally localized Wannier funcitons.
Anyway, documentaion on Wannier is on the way half.
Time consuming part (and also the advantage) is for effective interaction in RPA.
Look into the shell script genMLWF; you can skip last part if you don't need the effective interaction.


--------- this is the end for users -----------


##################################################################
######    LOG and MEMO for developers                 ############
##################################################################
=====
Real TIME and Memory usage command. (H.Kino). rev.jun23.2014
We set comment lines such as
----
!TIME0
  ... here lines of fortran codes ...
!TIME1 "foo bar"
----
in fortran files, e.g., in fpgw/main/hsfp0.sc.m.F.
This reports computational time at the end of output from each node.
(see stdout.{rankID}.* files).

For example, if we see a line (in stdout.*)
>timediff,id= 106 (foo_bar_____)    420.00 sec     4 times
, we used the CPU time 420 sec (roughly =4 x 105 sec) for
the lines between '!TIME0' and '!TIME1 "foo bar"'.

See fpgw/exec/makefile. At its bottom,
gawk replaces '!TIME*' with fortran codes (time_xxx.F are
generated automatically by this procedure).

--
NOTE: !TIME can be nested as
!TIME0
  ... code1 ...
!TIME0
  ... code2 ...
!TIME1 "foobar1"
!TIME1 "foobar2"


=====
1.branch x0mpi_merge (6b5e3, Aug2013): 
contains a try to paralellize x0kf
(imaginary part of dielectric function) more. It just show a
strategy, not efficient yet.

2. With FSMOM, Efermi is not uniquely given in job_band_nspin2*.
It is given by a bndfp-bzwtsf-bzwsf L300 block
   if ((.not. lfill) .or. (metal .and. (nkp .eq. 1))) then
   (bisection method to determine a middle of LUMO and HOMO).
It can give some energy between LUMO and HOMO. 
Small changes of computational condition can give large change. 
But no problem.

3. Line length for fortran; Add .emacs the following three lines.
(add-hook 'fortran-mode-hook
	  '(lambda ()
	     (setq fortran-line-length 132)))



======
2014 nov23:
job_band si -np 4 [options], where
options can be "-vso=1 -vnspin=2" for SO case.
--->Need to add documents to GetStarted.

===
--ssig option (ScaledSigma option).
Need explanation...

======
PDOS: sigm_fbz is required.
(when cp sigm,rst,GWinput ->LDA-like result.
 Then cp sigm_fbz ->it fails.
 Need to make new directory, and copy rst,sigm_fbz.)
And how to check it. (whether 


======
mixbeta:
takao@TT4:~/ecalj/fpgw$ grep mixbeta */*.F
main/hqpe.sc.m.F:      call getkeyvalue("GWinput","mixbeta",beta,default=1d0,status=ret)
mixing parameter on sigm file.
As the default beta is unitiy, mixsigm and mixsigma files are 


=======
Check convergecne on QSGW.
grep rms lqpe*


======
other to DO
(1)
 modify shell script to bash:
 use function. (learn git, push your new package to github).
(2) 
 fortran.
 Read dataflow and data structure.
 Look into gwsc and 

 You can insert your own note in souce code.

</pre>

end of file
--------- 


