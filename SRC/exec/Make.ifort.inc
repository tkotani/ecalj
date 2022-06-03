### NOTE. we use mpi compilar, but lmfa, lmdos lmchk lmf2gw lmfham1 are only for -np 1
FC = mpiifort
LK = mpiifort
#-qopenmp 
#-xHost -mcmodel=medium -heap-arrays 100

LIBMATH= -mkl  

##############################

### INTEL FORTRAN LINUX ###
module = -module $(moddir) 

FFLAGS = -O2 -132 $(module)
FFLAGS_LESS = -O0 -132 $(module)
FFLAGS_NONE = -O0 -132 $(module)

#FFLAGS = -O2 -132 -cpp $(CPP_SW) $(module)
#FFLAGS_LESS = -O0 -132 -cpp $(CPP_SW) $(module)
#FFLAGS_NONE = -O0 -132 -cpp $(CPP_SW) $(module)

#FFLAGS_OMP= -qopenmp -O3 -cpp $(CPPSWITCH_INTELLINUXIFC)
#-openmp
#FFLAGS=-O0 -check bounds -traceback -g -cpp $(CPPSWITCH_INTELLINUXIFC)
#FFLAGS_OMP=-O0 -check bounds -traceback -g -cpp $(CPPSWITCH_INTELLINUXIFC)
#

CPP_SW =
#-DMPIK 

LIBLOC = ${LIBMATH}
LKFLAGS2 = $(LIBMATH) 

