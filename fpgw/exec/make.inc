### you have to set FC,LK, FFLAGS, LIBMATH, ECAL, and BINDIR ###

# ##### for a gfortran. ########
# FC = mpif90  #compilar
# LK = mpif90  #linker
# FFLAGS=  -O3  $(CPPSWITCH_INTELLINUXIFC) -fomit-frame-pointer \
#         -funroll-loops -ffast-math -ffixed-line-length-132  ### Compilar options
# LIBMATH= /usr/lib/libfftw3.so.3 /usr/lib/liblapack.so.3gf \
#          /usr/lib/libblas.so.3gf   ### FFT+ BLAS+LAPACK
# ECAL   = $(HOME)/ecal     ### Root of ecal
# BINDIR = $(HOME)/bin      ### to which "make install (and install2) copy binaries.

# ##### for an intel. ########
FC=mpiifort -heap-arrays 100 -mcmodel=medium 
LK=mpiifort -heap-arrays 100 -mcmodel=medium 
#FC=mpiifort -C -CB -g  -mcmodel=medium 
#LK=mpiifort -C -CB -heap-arrays 100 -g  -mcmodel=medium 

#FC = mpif90
#LK = mpif90
FFLAGS = -O3 -heap-arrays 100 -cpp $(CPPSWITCH_INTELLINUXIFC) -132 
LIBMATH = -mkl
ECAL   = $(HOME)/ecal     ### Root of ecal
BINDIR = $(HOME)/bin     ### to which "make install (and install2) copy binaries.



#############################################################
#####   you usually don't need to touch followings      #####
#############################################################

LKFLAGS2 = $(LIBMATH) 

## CPP SWITCHES ###################
CPPSWITCH_INTELLINUXIFC  = \
-DEXPAND_ISWAP  -DEXPAND_VDV   -DCOMMONLL  -UDUMMY_ETIME -DEXPAND_MELPLN2         \
-DUSE_X0KBLAS   -DX0KBLAS_DIV  -UEXPAND_SUBSTITUTION     -UCOMMENTOUTfor_PARALLEL \
-DMbytes_X0KBLAS_DIV=2        -DNWORD_RECORDSIZE=1     -DEXPAND_SORTEA


### Don't change para_g = .o ... below (or modify it if you know how this work) 
#### don't need to read here #####NoteStart
# Some groups of .f sources are compiled into .c*_o files.  (* is 1 to 4).
# The compile options are in FFLAGS_c*. The others are with .o and FFLAGS. See makefile and Search para_g or so.
# ---> It cause a problem if a source file foo.f, which compiled into foo.c*_o contains USE module, 
#      because checkmodule does now just support *.o. In such a case, you have to modify checkmodule by yourself.
#      (This note is by takao. Oct.2003)
##################################NoteEnd
para_g = .o     # ppbafp.f  psi2bc1.f psi2bc.f See makefile.
sxcf_g = .o     # sxcf.f
x0kf_g = .o     # x0kf.f
hqpe_g = .o     # hqpe.f
tet5_g = .o
