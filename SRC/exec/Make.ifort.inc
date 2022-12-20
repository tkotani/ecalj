### NOTE. we use mpi compilar, but lmfa, lmdos lmchk lmf2gw lmfham1 are only for -np 1
#-qopenmp -xHost -mcmodel=medium -heap-arrays 100

FC = mpiifort
LK = mpiifort
LIBMATH= -mkl  
LIBLOC = ${LIBMATH}
LKFLAGS2 = $(LIBMATH) 
module = -module $(moddir) 

#FFLAGS = -O2 -132 $(module)
#FFLAGS_LESS = -O0 -132 $(module)
#FFLAGS_NONE = -O0 -132 $(module)
FFLAGS = -O2 -132 $(module)
FFLAGS_LESS = -O0 -132 $(module)
FFLAGS_NONE = -O0 -132 $(module)


