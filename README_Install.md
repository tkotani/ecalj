# Install memo for ecalj (checked 2025Apr)

We use python and fortran in ecalj package, where we require packages and libraries for python and fortran.
In addition, we need some tools such as cmake.

## How to install 
### Step 1. Get ecalj package and get tools.
#### tools  and libraries
We may need to install following tools and libraries.We need
> git, gfortran, openmpi, bash, cmake intel-mel

+ I use mpirun (Open MPI) 4.1.6 for ubuntu24.
We need python 3 (usually already in ubuntu. Type \>python (ctrl+D for quit)).
+ git makes things easier. Especiall for version up. >git diff at ecalj/ shows orginal and your modification.

With pip, we need to install
> plotly, spglib, seekpath, numpy

If your system is old, use mise or something to use latest python. It install things locally at ./local.

#### Get ecalj package 
```bash
git clone https://github.com/tkotani/ecalj.git # Get source code  
```
After you did the above git clone command, a directory ecalj/ appears.
We can check history of ecalj by ">gitk --all" at ecalj/ directory after you got git clone.

### Steps 2 InstallAll
Run
```bash
FC=gfortran ./InstallAll
```
at ecalj/. It generates ecalj/SRC/Makefile with cmake, then compiles all the fortran code.
Then we move to  ecalj/SRC/TestInstall and run tests.
If you see "OK! All passed!", you have succeeded.  
(FC=foobar can be ifort or nvfortran)
In advance, Check InstallAll. We will copy all binaries to BINDIR.

Set command path BINDIR. For example, write
```
PATH="~/bin/:$PATH"
```
in your .bashrc when you move all ecalj binaries to your ~/bin.

### Install VEST and getsyml
It is convenient to see structures with VESTA.
(I installed VESTA-gtk3.tar.bz2 (ver. 3.5.8, built on Aug 11 2022, 23.8MB) on ubuntu 24)
At ecalj/StructureTool/, we have 'viewvesta' command. Try 
```
viewvesta ctrl.si
```
to check the structure in viewer.At /StructureTool, we have exchange converters, 
```vasp2ctrl``` and ```ctrl2vasp```. These allows convert structures with POSCAR.

In addition, we need to install getsyml.py to obtain symmetry line for band plot.
Generated syml.* is used for the band plot in ecalj. (syml is a little strange... we will fix).
As long as you have spglib and seekpath, we don't need extra things to do.
But here is a memo for install [./GetSyml/README.org](./GetSyml/README.org).


#### Additional memo
* When InstallAll have finished, we have all required binaries and shell scripts in your ~/bin/ directory, 
or somewhere else where BINDIR specified in InstallAll

* Clean up by CleanAll:  
If something wrong, run "./CleanAll" at ecalj/ and redo installation.
Look into CleanAll. You may need to do 'rm -f CMakefiles CMakeCache.txt' at ecalj/SRC/exec/.

* Compile fortran only.
To compile fortran source only, move to ecalj/SRC/exec/ and run
>FC=fortran cmake . -D CMAKE_BUILD_TYPE=Debug
, for example. You may look into CMakeLists.txt.
Remove CMakeCache.txt and CMakeFiles/ if you want to recompile.

* Compilar bug: In cases, we have troubles due to the compilar.
Usually we use -O2 in CMakeList.txt. 
But we may need to use -O1 or -O0 for some files to avoid compilar bugs.
We may set some conditional compilation settings. May need to examine CMakelists.text

* Souce codes, Test, make system are under SRC/
SRC/ 
├── TestInstall : Root of Install test 
├── exec        : CMakeLists.txt and scripts
├── main        : All main *.f90
└── subroutines : All subrouitnes. *.f90. 
All fortran codes are in main/ and subrouitnes/ 
We have a CMakeLists.txt which generates Makefile. Look into it.

* Install test system at ecalj/SRC/TestInstall.
We have a test system with make at ecalj/SRC/TestInstall. Look into test.py and testecalj.py.
These controls all the test. 

I think it is not so difficult to add your own test to testecalj.py.
You have to compute something at first. Then inputs and minimum results are stored in a directory.
Then you describe the test in testecalj.py.

To test all of binaries, just do
>./test.py
>./test.py gwall  !tests only GW part.  
>./test.py si_gwsc  nio_gwsc !test si_gwsc and nio_gwsc only.

* openmpi failed on ubuntu22.   I obseved that gfortran+openmpi failed for ubuntu22. Use mpich.
  But I don't know cucrrent status.
 
