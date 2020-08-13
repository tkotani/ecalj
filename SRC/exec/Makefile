### This file is not machine-dependent. #####
####################################################
# How to install. 
# 1. Set $PLATFORM of your machine. 
#   Make.inc.$PLATFORM should exist.
# 2.You may need gnumake. MAKE=gmake
# 3. Type like
#    > make PLATFORM=gfort
####################################################
PLATFORM = gfortran
BINDIR   = $(HOME)/bin

#-------------------------------------------------------
obj_path= OBJ.$(PLATFORM)
main  = ../main/
subr  = ../subroutines/
maxloc= ../wanniergw/
moddir = $(obj_path)
module = -J$(moddir) -I$(moddir)
tags = .

progslm = lmf-MPIK lmfgw-MPIK lmfa lmdos lmchk lmf2gw  lmfham1
progsgw = hbasfp0 hvccfp0 hx0fp0 hsfp0 hef hqpe hqpe_qsgw qg4gw gwinit heftet hmergewv rdata4gw_v2 convgwin hx0fp0_sc hsfp0_sc hqpe_sc kino_input_test hecor eout eout2 hmagnon hhomogas
progsw  = huumat hpsig hmaxloc huumat_MPI hpsig_MPI wanplot hwmatK hwmatK_MPI htbplot # hx0fp0_mlw 
script  = cleargw* dqpu eps* gw_* gwg* gws* gwp*  mkG* chipm_mat genMLWF genMLWF_vw genMLWF2 genMLWFmod run_argloc RPAWannier* FLEX_interaction.py

all: mkobj progs
mkobj:
	@echo './$(obj_path)/ for .o and .mod '
	@[ -d $(obj_path) ] || mkdir  $(obj_path) 
moduledepends.$(PLATFORM).inc: 
	@echo 'Generate module-dependency moduledepends.$(PLATFORM).inc'
	@[ -d moduledepends.inc ] || exec ../../TOOLS/checkmodule $(subr)*.F $(main)*.F $(main)*.F $(maxloc)*.F
	sed -e "s/\.\.\/\w*\//$(obj_path)\//g" moduledepends.inc |uniq >moduledepends.$(PLATFORM).inc

include moduledepends.$(PLATFORM).inc #this is automatically generated
include Make.$(PLATFORM).inc
include makesrc.inc #define subslm, GWLIB_list

.PHONY:  progs clean cleanall init echotest mkdep
.SUFFIXES: .o .F90 .F .mod

ELIB = $(addprefix ./$(obj_path)/,$(ELIB_list))
MLOC = $(addprefix ./$(obj_path)/,$(MLOC_list))
MPI  = $(addprefix ./$(obj_path)/,$(MPI_list))

$(obj_path)/%.o: $(subr)%.F
	$(FC) $(FFLAGS) -c $< -o $@ 
$(obj_path)/%.o: $(subr)%.F90
	$(FC) $(FFLAGS) -c $< -o $@ 
$(obj_path)/%.o: $(maxloc)%.F
	$(FC) $(FFLAGS) -c $< -o $@ 
$(obj_path)/%.o: $(maxloc)%.F90
	$(FC) $(FFLAGS) -c $< -o $@ 
$(obj_path)/%.o: $(main)%.F
	$(FC) $(FFLAGS) -c $< -o $@ 
$(obj_path)/%.o: $(main)%.F90
	$(FC) $(FFLAGS) -c $< -o $@ 

#%.o:%.F
#	$(FC) $(FFLAGS) -c $<
#%.o:%.F90
#	$(FC) $(FFLAGS) -c $<

echotest:
	echo $(MPI)
	echo $(ELIB)
	echo $(MLOC)

### version info for show_programinfo.F
result_fc_version :=$(shell LANG=us; ${FC} -v 2>&1|tail -1)
result_uname_a := $(shell uname -a | head -c 120)
result_linked_at := $(shell LANG=us; date)
#result_etc_issue := $(shell cat /etc/issue|head -1 -c120)
result_etc_issue := $(shell grep -v \"^\$\" /etc/issue | head -1 )
result_git_branch := $(shell cat ../../.git/HEAD| sed 's/^ref: *//')
result_git_commit := $(shell cat ../../.git/`cat ../../.git/HEAD| sed 's/^ref: *//' `)
result_fflags := $(shell echo ${FFLAGS} | head -1 -c 120)
result_libloc := $(shell echo ${LIBLOC} | head -1 -c 120)
VERSIONFLAGS = -D___UNAME_A___="\"$(result_uname_a)\""   \
    -D___LINK_TIME___="\"$(result_linked_at)\"" \
    -D___ETC_ISSUE___="\"$(result_etc_issue)\"" \
    -D___GIT_BRANCH___="\"$(result_git_branch)\"" \
    -D___GIT_COMMIT___="\"$(result_git_commit)\"" \
    -D___PLATFORM___="\"$(PLATFORM)\"" \
    -D___FC_VERSION___="\"$(result_fc_version)\"" \
    -D___FFLAGS___="\"$(result_fflags)\"" \
    -D___LIBLOC___="\"$(result_libloc)\""

########################################
# Generate o of main routines for lmf 
########################################
$(obj_path)/lmfa.o: $(main)lmfav7.F 
	$(FC) $(FFLAGS_NONE) -DLMFA -c $< -o $@
$(obj_path)/lmchk.o: $(main)lmv7util.F 
	$(FC) $(FFLAGS_NONE) -DLMCHK -DFP -c $< -o $@
$(obj_path)/lmdos.o: $(main)lmv7util.F 
	$(FC) $(FFLAGS_NONE) -DLMDOS -c $< -o $@
$(obj_path)/lmf2gw.o: $(main)lmf2gw.F 
	$(FC) $(FFLAGS_NONE) -c $< -o $@
$(obj_path)/lmfham1.o: $(main)lmfham1.F
	$(FC) $(FFLAGS_NONE) -c $< -o $@
$(obj_path)/lmf-MPIK.o: $(main)lmv7.F $(ELIB)
	$(FC) $(FFLAGS_NONE) -DLMF -c $(main)lmv7.F -o $@
$(obj_path)/lmfgw-MPIK.o: $(main)lmv7.F
	$(FC) $(FFLAGS_NONE) -DLMF -DLMFGWD -c $(main)lmv7.F -o $@
$(obj_path)/show_programinfo.o: $(subr)show_programinfo.F  $(main)*.F $(subr)*.F 
	$(FC) $(FFLAGS) $(VERSIONFLAGS) $(CPP_SW) -c $< -o $@

#$(obj_path)/lmf.o: $(main)lmv7.F $(ELIB)
#	$(FC) $(FFLAGS_NONE) -DLMF -c $< -o $@
#$(obj_path)/lmfgw.o: $(main)lmv7.F 
#	$(FC) $(FFLAGS_NONE) -DLMF -DLMFGWD -c $< -o $@
#$(obj_path)/show_programinfo_mpik.o:  $(subr)show_programinfo.F $(main)*.F $(subr)*.F 
#	$(FC) $(FFLAGS)  $(VERSIONFLAGS) $(CPP_SW) -c $< -o $@

########################################
# progslm
########################################
lmfa:	$(obj_path)/lmfa.o $(ELIB)  $(MPI) 
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $(obj_path)/$@.o $(ELIB)  $(MPI) $(LIBLOC) $(LKFLAGS2) -o $@
lmchk:	$(obj_path)/lmchk.o $(ELIB) $(MPI) 
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $(obj_path)/$@.o $(ELIB)  $(MPI) $(LIBLOC) $(LKFLAGS2) -o $@
lmdos:	$(obj_path)/lmdos.o $(ELIB)  $(MPI) 
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $(obj_path)/$@.o $(ELIB)  $(MPI) $(LIBLOC) $(LKFLAGS2) -o $@
lmf2gw:	$(obj_path)/lmf2gw.o 
	@echo $(LK) $@.o 'link with *.o'
	@$(FC) $(LKFLAGS1) $(obj_path)/$@.o $(ELIB)  $(MPI) $(LIBLOC) $(LKFLAGS2) -o $@
lmfham1: $(obj_path)/lmfham1.o $(ELIB)  $(MPI) 
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $(obj_path)/$@.o $(ELIB)  $(MPI) $(LIBLOC) $(LKFLAGS2) -o $@
lmf-MPIK: $(obj_path)/lmf-MPIK.o $(ELIB) $(MPI) 
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $(obj_path)/$@.o $(ELIB) $(MPI)  $(LIBLOC) $(LKFLAGS2) -o $@
lmfgw-MPIK: $(obj_path)/lmfgw-MPIK.o $(ELIB)  $(MPI) 
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $(obj_path)/$@.o $(ELIB) $(MPI)  $(LIBLOC) $(LKFLAGS2) -o $@

# lmf:	$(obj_path)/lmf.o $(ELIB) 
# 	@echo $(LK) $@.o 'link with *.o'
# 	@$(LK) $(LKFLAGS1) $(obj_path)/$@.o $(ELIB) $(LIBLOC) $(LKFLAGS2) -o $@
# lmfgw:	$(obj_path)/lmfgw.o $(ELIB)
# 	@echo $(LK) $@.o 'link with *.o'
# 	@$(LK) $(LKFLAGS1) $(obj_path)/$@.o $(ELIB) $(LIBLOC) $(LKFLAGS2) -o $@

###################################################
progs: $(progslm) $(progsgw) $(progsw)

clean:
	rm -f  $(progslm)  $(progsgw) $(progsw)
cleanall:
	rm -rf $(progslm) $(progsgw) $(progsw) $(obj_path) moduledepends.* 
install: 
	cp $(script) $(progslm) $(progsgw) $(progsw) ctrlgenM1.py job* run_arg  $(BINDIR)
mkdep: moduledepends.$(PLATFORM).inc
dep:
	@echo 'This generate a call-caller data set for lmv7.F fp/*.F and subs/\*.F subs/*.F90 slatsm/*.F .'
	@echo 'No arguments required. This drives ../TOOLS/f_calltree.py '
	@echo 'HELP --> ../TOOLS/FparserTools/f_calltree.py --help, and read CallCaller.sh'
	@echo
	@echo '--- Now generating a file 'callcaller.dat' ... Wait!!! It takes 1 minute or so!'
	@echo '   If you like to apply this to other programs modify this script'
	@echo ' NOTE: T.Kotani is not sure whether this is reliable enough .. let me know something wrong...'
	../TOOLS/FparserTools/f_calltree.py lmv7.F ../srclm/*.F* srcgw/*.F* >callcaller.dat 2>callcaller.err
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
doxygen:
	cd $(tags);doxygen;cd ./latex;make
	echo 'fpgw/latex/refman.pdf generated'

#---------------------------------------------------
hecor:  $(main)hecor.o   $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI)  $(ELIB) $(LKFLAGS2) -o $@
eout:   $(main)eout.o   $(MPI)  $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI)  $(ELIB) $(LKFLAGS2) -o $@
eout2: $(main)eout2.o   $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI)  $(ELIB) $(LKFLAGS2) -o $@
gwinit:	$(main)gwinit.m.o $(MPI)    $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI)  $(ELIB) $(LKFLAGS2) -o $@
qg4gw:	$(main)qg4gw.m.o  $(MPI)   $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI)  $(ELIB) $(LKFLAGS2) -o $@
rdata4gw_v2: $(main)rdata4gw_v2.m.o  $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI) $(ELIB) $(LKFLAGS2) -o $@

hbasfp0: $(obj_path)/hbasfp0.m.o  $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI) $(ELIB) $(LKFLAGS2) -o $@
hvccfp0: $(obj_path)/hvccfp0.m.o  $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI) $(ELIB) $(LKFLAGS2) -o $@
hx0fp0:	 $(obj_path)/hx0fp0.m.o  $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI) $(ELIB) $(LKFLAGS2) -o $@
hmagnon: $(obj_path)/hmagnon.m.o  $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI) $(ELIB) $(LKFLAGS2) -o $@
hhomogas: $(obj_path)/hhomogas.m.o  $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI) $(ELIB) $(LKFLAGS2) -o $@
hx0fp0_sc: $(obj_path)/hx0fp0.sc.m.o  $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI)  $(ELIB) $(LKFLAGS2) -o $@
hsfp0:	$(obj_path)/hsfp0.m.o  $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI) $(ELIB) $(LKFLAGS2) -o $@
hsfp0_sc: $(obj_path)/hsfp0.sc.m.o  $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI)  $(ELIB) $(LKFLAGS2) -o $@
heftet:	$(obj_path)/heftet.m.o  $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI)  $(ELIB) $(LKFLAGS2) -o $@
hef:	$(obj_path)/hef.m.o  $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI) $(ELIB) $(LKFLAGS2) -o $@
hqpe:	$(obj_path)/hqpe.m.o   $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI) $(ELIB) $(LKFLAGS2) -o $@
hqpe_sc:$(obj_path)/hqpe.sc.m.o   $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI)  $(ELIB) $(LKFLAGS2) -o $@
hqpe_qsgw: $(obj_path)/hqpe_qsgw.m.o   $(MPI)  $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI) $(ELIB) $(LKFLAGS2) -o $@
hmergewv: $(obj_path)/hmergewv.m.o $(MPI)  $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI) $(ELIB) $(LKFLAGS2) -o $@
convgwin: $(obj_path)/convgwin.o $(MPI)  $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI) $(ELIB) $(LKFLAGS2) -o $@
kino_input_test: $(obj_path)/kino_input_test.o  $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $<  $(MPI)  $(ELIB) $(LKFLAGS2) -o $@

hmaxloc: $(obj_path)/hmaxloc.o $(MLOC)  $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(MLOC)  $(MPI)  $(ELIB) $(LKFLAGS2) -o $@
hpsig: $(obj_path)/hpsig.o $(MLOC)  $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(MLOC)  $(MPI)  $(ELIB) $(LKFLAGS2) -o $@
hpsig_MPI: $(obj_path)/hpsig_MPI.o $(MLOC)  $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'	
	@$(LK) $(LKFLAGS1) $< $(MLOC)  $(MPI)  $(ELIB) $(LKFLAGS2) -o $@
huumat:	$(obj_path)/huumat.o $(MLOC)  $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(MLOC)  $(MPI)  $(ELIB) $(LKFLAGS2)  -o $@
huumat_MPI: $(obj_path)/huumat_MPI.o $(MLOC)  $(MPI)  $(ELIB) $(MPIMLOC)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(MLOC)  $(MPI)  $(ELIB) $(MPIMLOC) $(LKFLAGS2) -o $@
hwmatK:	$(obj_path)/hwmatK.o    $(MLOC)  $(MPI) $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(MLOC)  $(MPI)  $(ELIB) $(LKFLAGS2)  -o $@ 
hwmatK_MPI: $(obj_path)/hwmatK_MPI.o $(MLOC)  $(MPI)  $(ELIB) $(MPIMLOC)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(MLOC)  $(MPI)  $(ELIB) $(MPIMLOC) $(LKFLAGS2)  -o $@ 
wanplot : $(obj_path)/wanplot.o ${MLOC} ${GWLIB} $(MPI)  $(ELIB)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(MLOC) $(MPI)  $(ELIB) $(LKFLAGS2)  -o $@ 
htbplot :  $(obj_path)/htbplot.o
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(LKFLAGS2) -o $@
#----------------------------
tags: $(progs)
	cd $(tags);etags $(main)/*.F $(subr)/*.F $(maxloc)/*.F

