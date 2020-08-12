### NOTE. we use mpi compilar, but lmfa, lmdos lmchk lmf2gw lmfham1 are only for -np 1
FC = mpif90 
LK = mpif90 
#-qopenmp 
#-xHost -mcmodel=medium -heap-arrays 100

LIBMATH= -mkl  

##############################

### INTEL FORTRAN LINUX ###
FFLAGS = -O2 -132 -cpp $(CPP_SW) $(module)
FFLAGS_LESS = -O0 -132 -cpp $(CPP_SW) $(module)

#FFLAGS_OMP= -qopenmp -O3 -cpp $(CPPSWITCH_INTELLINUXIFC)
#-openmp
#FFLAGS=-O0 -check bounds -traceback -g -cpp $(CPPSWITCH_INTELLINUXIFC)
#FFLAGS_OMP=-O0 -check bounds -traceback -g -cpp $(CPPSWITCH_INTELLINUXIFC)
#

CPP_SW = -DHASIARGC  -DHASGETARG  -DFDATE  -DHASCPUTIME -DHASGETENVIRONMENTVARIABLE  -DNOQUAD  \
 -DNOT_OVERLAP_VAR -DMPIK \
 -DEXPAND_ISWAP  -DEXPAND_VDV   -DCOMMONLL  -UDUMMY_ETIME -DEXPAND_MELPLN2         \
 -DUSE_X0KBLAS   -DX0KBLAS_DIV  -UEXPAND_SUBSTITUTION     -UCOMMENTOUTfor_PARALLEL \
 -DMbytes_X0KBLAS_DIV=2        -DNWORD_RECORDSIZE=1 -DEXPAND_SORTEA \
 -DUSE_GEMM_FOR_SUM

LIBLOC = ${LIBMATH}
LKFLAGS2 = $(LIBMATH) 

