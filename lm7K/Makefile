######################################################################################################
# How to install. 
# 1. Choose machine dependent $PLATFORM below.
#
#   ./MAKEINC/Make.inc.$PLATFORM should exist.
#    Then you have to set math library and so on in the Make.inc.$PLATFORM.
#
# 2. Need gnumake. MAKE=gmake
#
# 3. Type 
#    prompt>make
#        or 
#    prompt>make PLATFORM=gfortran_mpik
#    or so.
#
#    If you write your own MAKEINC/Make.inc.$PLATFORM, please send it to us.
#
#    There are make options .PHONY targets.
#
#  * This make includes slatsm/make.inc subs/make.inc fp/make.inc. In addition,
#  * Be careful---> you can not add any space at the end of PLATFORM.
#  * When you switch to different platform in a machine, you have to clean up FPLOTdir/*.o *.o in advance by hand, 
#    or use 'make clean'.
######################################################################################################
BINDIR=$(HOME)/bin

# It will reads MAKEINC/Make.inc.$(PLATFORM).
#PLATFORM=fujitsufortran.fx10
#PLATFORM=fujitsufortran
#PLATFORM=fujitsufortran_mpik

PLATFORM=gfortran
LIBMATH=/usr/local/lib/libblas.a /usr/local/lib/liblapack.a /usr/local/lib/libfftw3.a

#PLATFORM=gfortran.fpp
#PLATFORM=gfortran_mpi
#PLATFORM=gfortran_mpik

#PLATFORM=ifort
#PLATFORM=ifort_mpi
#PLATFORM=ifort_mpik


###########################################################################
# Don't edit lines below 
###########################################################################


obj_path=obj.$(PLATFORM)

include slatsm/make.inc
include subs/make.inc
include fp/make.inc

#################################
.PHONY:  progs clean cleanbin cleanall veryclean distclean check checkgw install \
libs slatsm-dir subs-dir fp-dir mod-dir slatsm-lib subs-lib fp-lib


#################################


fplotprogs := fplot pldos

ALL: libs
	$(MAKE) -f Makefile progs
	-[ -f slatsm/slatsm.a ] && (cd FPLOTdir; $(MAKE) PLATFORM=$(PLATFORM) LIBLOC='$(LIBLOC)'; cp -f $(fplotprogs) ../)
#	(cd slatsm; make -f Makefile SLATSM=$(SLATSM))
#	(cd subs; make -f Makefile SUBS=$(SUBS))
#	(cd fp; make -f Makefile FP=$(FP))
#	(cd FPLOTdir; make;cp $(fplotprogs) ../)
#	make -f Makefile $(programs)
#
#	(cd nc; make -f Makefile )
#	(cd tb; make -f Makefile )


# Generate main routines
libs: mod-dir slatsm-dir subs-dir fp-dir
	$(MAKE) -f Makefile slatsm-lib
	$(MAKE) -f Makefile subs-lib
	$(MAKE) -f Makefile fp-lib

#slatsm-lib: slatsm-dir
#subs-lib: subs-dir
#fp-lib: fp-dir subs-lib


# initialization to check dependency of modules
init: moduledepends.inc

moduledepends.inc:
	exec ../TOOLS/checkmodule subs/*.F slatsm/*.F fp/*.F FPLOTdir/*.F
moduledepends.$(PLATFORM).inc: moduledepends.inc
	sed -e 's/\//\/$(obj_path)\//g' moduledepends.inc >moduledepends.$(PLATFORM).inc

dep:
	@echo  'This generate a call-caller data set for lmv7.F fp/*.F and subs/\*.F subs/*.F90 slatsm/*.F .'
	@echo  'No arguments required. This drives ../TOOLS/f_calltree.py '
	@echo  'HELP --> ../TOOLS/FparserTools/f_calltree.py --help, and read CallCaller.sh'
	@echo
	@echo  '--- Now generating a file 'callcaller.dat' ... Wait!!! It takes 1 minute or so!'
	@echo '        If you like to apply this to other programs, modify this script'
	@echo  ' NOTE: T.Kotani is not sure whether this is relaiable enough or not... let me know something wrong...'
	../TOOLS/FparserTools/f_calltree.py lmv7.F fp/*.F subs/*.F subs/*.F90 slatsm/*.F >callcaller.dat 2>callcaller.err
	-egrep -e '^(ERROR|Error)' callcaller.err
	@echo
	@echo '------------------------------------------------------------------------------'
	@echo '--- If no ERROR is shown above (if ERROR is not in callcaller.err), it is succeeded. ---'
	@echo '       Note that Unsed files might be used by other mainprogram.'
	@echo '--- If ERROR is shown above, look into callcaller.err. Something wrong.'
	@echo 
	@echo ' If you want to make a callcaller-tree picture, try'
	@echo ' >GenCCtree.sh callcaller.dotdata'
	@echo ' --> Then you get ccmap.ps.; it is better to use smaller callcaller.dotdata(need to modify this script to make it).' 
	@echo ' Note that you need graphviz for GenCCtree.sh. as apt-get install graphviz'

########################################
# general rules
########################################
.SUFFIXES: .o .F90 .F .mod

%.o:%.F
	$(FC) $(FFLAGS) -c $<

%.o:%.F90
	$(FC) $(FFLAGS) -c $<

sla_obj=$(addprefix slatsm/$(obj_path)/,$(sla_list))
sla_obj_path:=slatsm/$(obj_path)
slatsm-lib: $(sla_obj)  
slatsm-dir:
	[ -d $(sla_obj_path) ] || mkdir  $(sla_obj_path) 
$(sla_obj_path)/%.o: slatsm/%.F
	$(FC) $(FFLAGS) -c $< -o $@ 

subs_obj:=$(addprefix subs/$(obj_path)/,$(subs_list))
subs_obj_path:=subs/$(obj_path)
subs-lib: $(subs_obj)  
subs-dir:
	[ -d $(subs_obj_path) ] || mkdir  $(subs_obj_path) 
$(subs_obj_path)/%.o: subs/%.F
	$(FC) $(FFLAGS) -c $< -o $@ 
$(subs_obj_path)/%.o: subs/%.F90
	$(FC) $(FFLAGS) -c $< -o $@ 

fp_obj=$(addprefix fp/$(obj_path)/,$(fp_list))
fp_obj_path=fp/$(obj_path)
fp-lib: $(fp_obj)  
fp-dir:
	[ -d $(fp_obj_path) ] || mkdir  $(fp_obj_path) 
$(fp_obj_path)/%.o: fp/%.F
	$(FC) $(FFLAGS) -c $< -o $@ 


#subs/obj.gfortran/m_rdctrl.o:subs/obj.gfortran/m_gtv.o

###############################################
# module dependency (fixed). What *.F contains what modules.
################################################
# $(subs_obj_path)/m_rdctrl.o:    $(subs_obj_path)/m_gtv.o $(subs_obj_path)/m_toksw.o $(subs_obj_path)/m_struc_def.o $(subs_obj_path)/m_struc_func.o 
# $(subs_obj_path)/m_rdctrlchk.o: $(subs_obj_path)/m_gtv.o $(subs_obj_path)/m_toksw.o
# $(subs_obj_path)/rdctrl2.o:     $(subs_obj_path)/m_gtv.o $(subs_obj_path)/m_rdctrl.o $(subs_obj_path)/m_struc_def.o $(subs_obj_path)/m_struc_func.o
# $(subs_obj_path)/rdctrlchk.o:   $(subs_obj_path)/m_gtv.o $(subs_obj_path)/m_rdctrlchk.o
# $(fp_obj_path)/bndfp.o: $(subs_obj_path)/m_rdctrl.o
# #
# $(subs_obj_path)/struc_main.o: $(susb_obj_path)/m_struc_def.o $(susb_obj_path)/m_struc_func.o
# $(subs_obj_path)/struc_sub.o: $(susb_obj_path)/m_struc_def.o $(susb_obj_path)/m_struc_func.o
# $(subs_obj_path)/m_struc_func.o: $(susb_obj_path)/m_struc_def.o 

########################################
# machie dependent rules
# I must load Make.inc here
########################################

#include Make.mod.dependency
include MAKEINC/Make.inc.$(PLATFORM)
include moduledepends.$(PLATFORM).inc
#include MAKEINC/Make.inc.$(PLATFORM).patch #this is for less optimization to avoid compilar bugs.

############################################
# version information
############################################
# This must be placed after including Make.inc.PLATFORM
result_uname_a := $(shell uname -a | head -c 120)
result_linked_at := $(shell LANG=C; date)
#result_etc_issue := $(shell cat /etc/issue|head -1 -c120)
result_etc_issue := $(shell grep -v \"^\$\" /etc/issue | head -1 )
result_git_branch := $(shell cat ../.git/HEAD| sed 's/^ref: *//')
result_git_commit := $(shell cat ../.git/`cat ../.git/HEAD| sed 's/^ref: *//' `)
result_fflags := $(shell echo ${FFLAGS} | head -1 -c 120)
result_libloc := $(shell echo ${LIBLOC} | head -1 -c 120)
#result_fc_version := is placed at Make.inc.PLATFORM

VERSIONFLAGS = -D___UNAME_A___="\"$(result_uname_a)\""   \
    -D___LINK_TIME___="\"$(result_linked_at)\"" \
    -D___ETC_ISSUE___="\"$(result_etc_issue)\"" \
    -D___GIT_BRANCH___="\"$(result_git_branch)\"" \
    -D___GIT_COMMIT___="\"$(result_git_commit)\"" \
    -D___PLATFORM___="\"$(PLATFORM)\"" \
    -D___FC_VERSION___="\"$(result_fc_version)\"" \
    -D___FFLAGS___="\"$(result_fflags)\"" \
    -D___LIBLOC___="\"$(result_libloc)\""

mod-dir:
	[ -d MODDIR ] || mkdir  MODDIR 
	[ -d $(moddir) ] || mkdir  $(moddir) 

########################################
# Generate main routines
########################################
lmchk.o: lmv7util.F
	$(FC) $(FFLAGS_NONE) -DLMCHK -DFP -c lmv7util.F -o $@
#lm.o: lmv7.F
#	$(FC) $(FFLAGS) -DLM -DNC -c lmv7.F -o $@
# non colinear version
#lm.o: lmv7.F
#	$(FC) $(FFLAGS) -DLM -DNC -c lmv7.F -o $@
lmfa.o: lmfav7.F
	$(FC) $(FFLAGS_NONE) -DLMFA -c lmfav7.F  -o $@
lmf.o: lmv7.F
	$(FC) $(FFLAGS_NONE) -DLMF  -c lmv7.F -o $@
lmdos.o: lmv7util.F
	$(FC) $(FFLAGS_NONE) -DLMDOS -c lmv7util.F -o $@
lmfgw.o: lmv7.F
	$(FC) $(FFLAGS_NONE) -DLMF -DLMFGWD -c lmv7.F -o $@
#lmscell.o: lmv7.F
#	$(FC) $(FFLAGS) -DLMSCELL -c lmv7.F -o $@
#tbe.o: lmv7.F
#	$(FC) $(FFLAGS) -DTBE -c lmv7.F -o $@
lmf2gw.o: lmf2gw.F
	$(FC) $(FFLAGS_NONE) -c lmf2gw.F -o $@
lmfham.o: lmfham.F
	$(FC) $(FFLAGS_NONE) -c lmfham.F -o $@

########################################
# programs
########################################
show_programinfo.o: make_fc_version ../show_programinfo.F
	$(FC) $(FFLAGS)  $(VERSIONFLAGS)  $(CPP_SW) -c ../show_programinfo.F -o show_programinfo.o

lmchk:	lmchk.o $(LIBFP) $(LIBSUBS) $(LIBSLA) show_programinfo.o
	$(LK) $(LKFLAGS1) $@.o show_programinfo.o $(LIBFP) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@
 
#lmscell: lmscell.o $(LIBSUBS) $(LIBSLA)
#	$(LK) $(LKFLAGS1) $@.o $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

lmdos:	lmdos.o $(LIBSUBS) $(LIBSLA) show_programinfo.o
	$(LK) $(LKFLAGS1) $@.o show_programinfo.o $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@
 
lmfa:	lmfa.o $(LIBSUBS) $(LIBFP) $(LIBSLA) show_programinfo.o
	$(LK) $(LKFLAGS1) $@.o show_programinfo.o $(LIBFP) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@
 
lmf:	lmf.o  $(LIBSUBS) $(LIBFP) $(LIBSLA) show_programinfo.o
	$(LK) $(LKFLAGS1) $@.o show_programinfo.o $(LIBFP) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@
#lm:	lm.o  $(LIBSUBS) 
#	$(LK) $(LKFLAGS1) $@.o $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@
# non-colinear verion
#lm:	lm.o  $(LIBSUBS) 
#	$(LK) $(LKFLAGS1) $@.o $(LIBNC) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@
 
lmfgw:	lmfgw.o $(LIBSUBS) $(LIBFP)  $(LIBSLA) show_programinfo.o
	$(LK) $(LKFLAGS1) $@.o  show_programinfo.o $(LIBFP) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

lmf2gw:	lmf2gw.o show_programinfo.o
	$(FC) $(LKFLAGS1) $@.o show_programinfo.o $(LKFLAGS2) -o $@

lmfham:	lmfham.o $(LIBSUBS) $(LIBFP)  $(LIBSLA) show_programinfo.o
	$(LK) $(LKFLAGS1) $@.o  show_programinfo.o $(LIBFP) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

#tbe:	tbe.o $(LIBSUBS) $(LIBTB)
#	$(LK) $(LKFLAGS1) $@.o $(LIBTB) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

#rdcmd:	rdcmd.o $(LIBSUBS) $(LIBSLA)
#	$(LK) $(LKFLAGS1) $@.o $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

#lm67:	lm67src/lm67.F
#	(cd lm67src; $(FC) $(FFLAGS_LESS67) lm67.F $(LIBLOC) -o lm67 ; mv lm67 ../)

####### for MPI
lmf-MPI.o: lmv7.F
	$(FC) $(FFLAGS_NONE) -DLMF -c lmv7.F -o $@

lmf-MPI:        lmf-MPI.o $(LIBSUBS) $(LIBFP) $(LIBSLA) show_programinfo.o
	$(LK) $(LKFLAGS1) $@.o show_programinfo.o $(LIBFP) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

lmfgw-MPI.o: lmv7.F
	$(FC) $(FFLAGS_NONE) -DLMF -DLMFGWD -c lmv7.F -o $@

lmfgw-MPI:      lmfgw-MPI.o $(LIBSUBS) $(LIBFP)  $(LIBSLA) show_programinfo.o
	$(LK) $(LKFLAGS1) $@.o show_programinfo.o $(LIBFP) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

####### for MPIK
lmf-MPIK.o: lmv7.F
	$(FC) $(FFLAGS_NONE) -DLMF -c lmv7.F -o $@

lmf-MPIK:       lmf-MPIK.o $(LIBSUBS) $(LIBFP) $(LIBSLA) show_programinfo.o
	$(LK) $(LKFLAGS1) $@.o show_programinfo.o $(LIBFP) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

lmfgw-MPIK.o: lmv7.F
	$(FC) $(FFLAGS_NONE) -DLMF -DLMFGWD -c lmv7.F -o $@

lmfgw-MPIK:     lmfgw-MPIK.o $(LIBSUBS) $(LIBFP)  $(LIBSLA) show_programinfo.o
	$(LK) $(LKFLAGS1) $@.o show_programinfo.o $(LIBFP) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@


##########################################
# .PHONY targets
##########################################
progs: $(programs)

clean:
	rm -f $(programs) $(fplotprogs) */*.a *.o 
	(cd FPLOTdir; rm *.o;$(MAKE) PLATFORM=$(PLATFORM) clean)

cleanbin:
	rm -f $(programs) $(fplotprogs)
	(cd FPLOTdir; $(MAKE) PLATFORM=$(PLATFORM) clean)


cleanall:
	rm -f $(programs) $(fplotprogs) $(LIBSUBS) $(LIBFP) $(LIBSLA) */obj.*/*.o  *.o  */*.a 
	(cd FPLOTdir; $(MAKE) PLATFORM=$(PLATFORM) cleanobj)
	rm -rf MODDIR
	mkdir MODDIR
 
veryclean:
	for p in $(programs) ; do \
	  rm -f $$p $$p.o ; \
	done
	testing/test.lm --quiet --clean
	fp/test/test.fp --all --quiet --clean
 
distclean:
	for p in $(programs) ; do \
	  rm -f $$p $$p.o ; \
	done
	testing/test.lm --quiet --clean
	fp/test/test.fp --all --quiet --clean
	rm -f v7input/subs.a subs/subs.a fp/subs.a 
	rm -f v7input/{m_gtv,m_rdctrl,m_toksw}.mod

check:
	(cd fp/test;./test.fp --all --quiet) 

checkgw:
	echo 'For GW driver, run samples in TESTinstallGW.'

#	testing/test.lm --quiet $(ADD0)


install: 
	cp $(programs) $(fplotprogs) ctrlgenM1.py job* run_arg  $(BINDIR)

### this is to avoid compilar bug of ifort bug (Version 14.0.1.106)
#$(subs_obj_path)/xcpbe.o : subs/xcpbe.F90
#	$(FC) $(FFLAGS_LESS) -c $< -o $(subs_obj_path)/xcpbe.o
