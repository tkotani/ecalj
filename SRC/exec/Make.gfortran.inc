### NOTE. we use mpi compilar, but lmfa, lmdos lmchk lmf2gw lmfham1 are only for -np 1
FC = mpif90 -Wsurprising -Wmaybe-uninitialized
#
LK = mpif90 
LIBMATH = -lmkl_rt

module = -J$(moddir) -I$(moddir)

# ##############################

#2023feb: For newer gfortran, we need -fallow-argument-mismatch
FFLAGS_COMMON= -g -fimplicit-none -finit-integer=NaN -finit-real=NaN $(module)
#FFLAGS_COMMON= -g -fallow-argument-mismatch -fimplicit-none -finit-integer=NaN -finit-real=NaN $(module)
FFLAGS = -O2  $(FFLAGS_COMMON)
FFLAGS_LESS = -O0   $(FFLAGS_COMMON)
FFLAGS_LESS2 = -O0  $(FFLAGS_COMMON)
FFLAGS_LESS3 = -O0  $(FFLAGS_COMMON)
FFLAGS_NONE = -O0   $(FFLAGS_COMMON)


#-fomit-frame-pointer -funroll-loops  -ffast-math
# safer option (debug)
#  -pg for profiling??? how to use gprof
# FFLAGS =      -g  -fbacktrace $(FFLAGS_COMMON)
# FFLAGS_LESS = -O0 -g  -fbacktrace $(FFLAGS_COMMON)
# FFLAGS_LESS2 = -O0 -g  -fbacktrace $(FFLAGS_COMMON)
# FFLAGS_LESS3 = -O0 -g  -fbacktrace $(FFLAGS_COMMON)
# FFLAGS_NONE = -O0 -g  -fbacktrace $(FFLAGS_COMMON)

#nfvortran???
#+FC = mpif90 
#+LIBMATH = -Mmkl
#+module = -module $(moddir) -I$(moddir)
#+FFLAGS_COMMON=  $(module)


#-DMPIK 
LIBLOC = ${LIBMATH}
LKFLAGS2 = $(LIBMATH) 

$(obj_path)/iopen.o: $(maxloc)/iopen.f90
	$(FC) $(FFLAGS_LESS3) -c $< -o $@
	
# ##### Patch section #############
# $(obj_path)/huntx.o: $(subr)/huntx.f90
# 	$(FC) $(FFLAGS_LESS3) -c $< -o $@
# $(obj_path)/hunti.o: $(subr)/hunti.f90
# 	$(FC) $(FFLAGS_LESS3) -c $< -o $@
# $(obj_path)/polcof.o: $(subr)/polcof.f90
# 	$(FC) $(FFLAGS_LESS3) -c $< -o $@
# $(obj_path)/rdfiln.o: $(subr)/rdfiln.f90
# 	$(FC) $(FFLAGS_LESS3) -c $< -o $@
# $(obj_path)/ropbes.o: $(subr)/ropbes.f90
# 	$(FC) $(FFLAGS_LESS3) -c $< -o $@
# $(obj_path)/ropyln.o: $(subr)/ropyln.f90
# 	$(FC) $(FFLAGS_LESS3) -c $< -o $@
# $(obj_path)/spcgrp.o  : $(subr)/spcgrp.f90
# 	$(FC) $(FFLAGS_LESS) -c $< -o $@
# $(obj_path)/m_toksw.o  : $(subr)/m_toksw.f90
# 	$(FC) $(FFLAGS_NONE) -c $< -o $@
# $(obj_path)/m_rdctrl.o : $(subr)/m_rdctrl.f90
# 	$(FC) $(FFLAGS_NONE) -c $< -o $@
# $(obj_path)/m_rdctrlchk.o : $(subr)/m_rdctrlchk.f90
# 	$(FC) $(FFLAGS_NONE) -c $< -o $@
# $(obj_path)/m_struc_def.o : $(subr)/m_struc_def.f90 
# 	$(FC) $(FFLAGS_NONE) -c $< -o $@
# $(obj_path)/m_struc_func.o : $(subr)/m_struc_func.f90 
# 	$(FC) $(FFLAGS_NONE) -c $< -o $@
# $(obj_path)/rhogkl.o  : $(subr)/rhogkl.f90
# 	$(FC) $(FFLAGS_LESS) -c $< -o $@

