### NOTE. we use mpi compilar, but lmfa, lmdos lmchk lmf2gw lmfham1 are only for -np 1
#-qopenmp -xHost -mcmodel=medium -heap-arrays 100

FC = mpiifort
LK = mpiifort
LIBMATH= -mkl  
LIBLOC = ${LIBMATH}
LKFLAGS2 = $(LIBMATH) 
# -check -traceback -debug extended

module = -module $(moddir) 
FFLAGS = -O2       -init:snan $(module)
FFLAGS_LESS = -O0  -init:snan $(module)
FFLAGS_NONE = -O0   -init:snan $(module)

#FFLAGS = -O2 -traceback  -init:snan $(module)
#FFLAGS_LESS = -O0  -traceback -init:snan $(module)
#FFLAGS_NONE = -O0  -traceback -init:snan $(module)


### conditional compile to avoid compilar bug for ucgw 2023feb ifort 18.0.5.274 
$(obj_path)/m_qplist.o: $(subr)m_qplist.f90
	$(FC) $(FFLAGS_NONE) -c $< -o $@
