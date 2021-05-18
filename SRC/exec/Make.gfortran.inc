### NOTE. we use mpi compilar, but lmfa, lmdos lmchk lmf2gw lmfham1 are only for -np 1
FC = mpif90 
LK = mpif90 
LIBMATH = -lmkl_rt

#FC = gfortran
#LK = gfortran
# LIBMATH =/usr/lib/x86_64-linux-gnu/libblas.so.3 \
# /usr/lib/x86_64-linux-gnu/liblapack.so.3 /usr/lib/x86_64-linux-gnu/libfftw3.so.3

module = -J$(moddir) -I$(moddir)

##############################
FFLAGS_COMMON= -g -fimplicit-none -ffixed-line-length-132 $(CPP_SW) $(module)
FFLAGS = -O2  $(FFLAGS_COMMON)
FFLAGS_LESS = -O0   $(FFLAGS_COMMON)
FFLAGS_LESS2 = -O0  $(FFLAGS_COMMON)
FFLAGS_LESS3 = -O0  $(FFLAGS_COMMON)
FFLAGS_NONE = -O0   $(FFLAGS_COMMON)
FLAGS_LESS67 = -O0  $(FFLAGS_COMMON)

#-fomit-frame-pointer -funroll-loops  -ffast-math
# safer option (debug)
#  -pg for profiling??? how to use gprof
#FFLAGS =      -g  -fbacktrace $(FFLAGS_COMMON)
#FFLAGS_LESS = -g  -fbacktrace $(FFLAGS_COMMON)
#FFLAGS_LESS2 = -g  -fbacktrace $(FFLAGS_COMMON)
#FFLAGS_LESS3 = -g  -fbacktrace $(FFLAGS_COMMON)
#FFLAGS_NONE = -g  -fbacktrace $(FFLAGS_COMMON)

CPP_SW = -DMPIK 
LIBLOC = ${LIBMATH}
LKFLAGS2 = $(LIBMATH) 

##### Patch section #############
$(obj_path)/huntx.o: $(subr)/huntx.F
	$(FC) $(FFLAGS_LESS3) -c $< -o $@
$(obj_path)/hunti.o: $(subr)/hunti.F
	$(FC) $(FFLAGS_LESS3) -c $< -o $@
$(obj_path)/polcof.o: $(subr)/polcof.F
	$(FC) $(FFLAGS_LESS3) -c $< -o $@
$(obj_path)/rdfiln.o: $(subr)/rdfiln.F
	$(FC) $(FFLAGS_LESS3) -c $< -o $@
$(obj_path)/ropbes.o: $(subr)/ropbes.F
	$(FC) $(FFLAGS_LESS3) -c $< -o $@
$(obj_path)/ropyln.o: $(subr)/ropyln.F
	$(FC) $(FFLAGS_LESS3) -c $< -o $@
$(obj_path)/spcgrp.o  : $(subr)/spcgrp.F 
	$(FC) $(FFLAGS_LESS) -c $< -o $@
$(obj_path)/m_toksw.o  : $(subr)/m_toksw.F 
	$(FC) $(FFLAGS_NONE) -c $< -o $@
$(obj_path)/m_rdctrl.o : $(subr)/m_rdctrl.F
	$(FC) $(FFLAGS_NONE) -c $< -o $@
$(obj_path)/m_rdctrlchk.o : $(subr)/m_rdctrlchk.F
	$(FC) $(FFLAGS_NONE) -c $< -o $@
$(obj_path)/m_struc_def.o : $(subr)/m_struc_def.F 
	$(FC) $(FFLAGS_NONE) -c $< -o $@
$(obj_path)/m_struc_func.o : $(subr)/m_struc_func.F 
	$(FC) $(FFLAGS_NONE) -c $< -o $@
$(obj_path)/rhogkl.o  : $(subr)/rhogkl.F 
	$(FC) $(FFLAGS_LESS) -c $< -o $@

