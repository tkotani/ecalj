
#### AIX/HITACHI SR11000 ###

# Never tested.


# ... patch section
$(SLATSM)(huntx.o): 
	$(FC) $(FFLAGS_LESS3) -c huntx.f
	ar rv $(SLATSM) huntx.o
	rm -f huntx.o
 
$(SLATSM)(hunti.o): 
	$(FC) $(FFLAGS_LESS3) -c hunti.f
	ar rv $(SLATSM) hunti.o
	rm -f hunti.o
 
$(SLATSM)(polcof.o): 
	$(FC) $(FFLAGS_LESS3) -c polcof.f
	ar rv $(SLATSM) polcof.o
	rm -f polcof.o
 
$(SLATSM)(rdfiln.o): 
	$(FC) $(FFLAGS_LESS3) -c rdfiln.f
	ar rv $(SLATSM) rdfiln.o
	rm -f rdfiln.o

# ... C compiler and flags
CC = cc
CFLAGS =  -DINTEL_IFC  -g


# ... Fortran compiler and flags, and linker ... for the INTEL IA32
FC = mpif90 -64
LK = mpif90 -64
FFLAGS = -Os  -noparallel -fixed=132 -cpp $(CPP_SW)
FFLAGS_LESS = -O2 -noparallel -fixed=132 -cpp $(CPP_SW)
FFLAGS_LESS2 = -O1 -noparallel -fixed=132 -cpp $(CPP_SW)
FFLAGS_LESS3 = -O0 -noparallel -fixed=132 -cpp $(CPP_SW)
FFLAGS_NONE = -g  -noparallel -fixed=132 -cpp $(CPP_SW)

# --- Libraries LIBSLA,LIBLOC ---
ECAL = ..
LIBSLA = $(ECAL)/slatsm/slatsm.a
LIBLOC = -L/home3/kino/kit/fftw-3.1.2/.libs -lfftw3 \
-L/opt/intel/mkl/9.1/lib/em64t -lmkl_lapack -lmkl -pthread -parallel
#-lsvml


# --- Shell variables to be supplied by configure ---
SHELL = /bin/sh
RANLIB = ranlib
AR = ar -X64 
ARFLAGS = -rv

# ... cpp
#CCOMP = ./ccomp
#CCOMPDIR = .
#CCOMP_SW =  -dINTEL_IFC -dF90 -dFFTW -dIN_PLACE
CPP_SW =   -DF90 -DFFTW -DIN_PLACE -DHASGETARG -DHASIARGC -DDOL_FDATE  -DMPI -UMPE 

# ifort v10
# -DHASGETARG -DHASNARGS -DFDATE
# hitachi/SR11000 
# -DHASGETARG -DHASIARGC -DDOL_FDATE 

# MPI
# -DMPI -UMPE 



# ... for compilers such as xlf that output numbers 0.nnn as .nnn
ADD0 = 
# ... Compiler extensions
FC_IS_F90 = yes
FC_AUTOARRAY = 
FC_POINTER = 


# ... Path to, and MPI-specific arguments for, the MPI fortran compiler
F90M = 


# ... Fortran-C linkage
FC_UPPERCASE = 0
FC_UNDERSCORE = 1
NOPASSEDARGS = 1
CMAIN = MAIN__
CMFLAGS = $(CFLAGS) -DFC_UNDERSCORE=1 -DFC_UPPERCASE=0 -DNOPASSEDARGS=1 -DCMAIN=MAIN__ -DNARGFCALL=nargs_ -DADD_TO_NARGFCALL=0 -DARGFCALL=getarg_



### See SLATSM48ORG/startup/README
#
#    CMAIN           is the name of the entry point function, e.g. CMAIN=main
#
#    ... the following describe how fortran function names are mangled
#    FC_UNDERSCORE   if 0, means fortran appends no underscore to a function
#                    if 1, means fortran appends an underscore to a function
#                    if 2, means fortran appends two underscores to function names that already contain underscores
#
#    FC_UPPERCASE    if 1, means fortran converts a function name to upper case
#
#
#    ... The following are used for extracting command-line arguments from fortran
#        There are two function calls fmain.c supplies for the command-line argument:
#        nargc()        returns the number of arguments
#        gtargc(iarg,s) returns particular argument iarg in string s.
#
#
#    NOPASSEDARGS=#  if #=0 then argc, argv are passed to program entry point CMAIN
#                    Functions nargc and gtargc just return data from supplied argc, argv.
#
#                    if #=1 then argc, argv are not passed to program entry point CMAIN
#                    and in an initialization step fmain.c uses the function calls to
#                    extract the information and store it locally.
#                    fmain then functions in the same way as NOPASSEDARGS = 0.
#
#    ... the following three tokens are used when NOPASSEDARGS=1
#    NARGFCALL=fn    name of function call that returns # arguments
#                    It need not be supplied; however, 
#                    nargc() always returns 0.
#
#    ADD_TO_NARGFCALL=strn  (optionally used in conjunction with NARGFCALL=fn)
#                    number-of-arguments = fn + strn
#                    Designed for implementations that return something different from # args, e.g.
#                    # args = iargc() + 1
#
#    ARGFCALL=fn     name of function call that returns one entry in argv
#                    It need not be supplied; however, gtargc is then not defined.
