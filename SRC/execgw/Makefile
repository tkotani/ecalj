### This file is not machine-dependent. #####
####################################################
# How to install. 
# 1. Set $PLATFORM of your machine. 
#   Make.inc.$PLATFORM should exist.
# 2.You may need gnumake. MAKE=gmake
# 3. Type like
#    > make PLATFORM=gfortran
####################################################
PLATFORM=gfortran  # Machine-specific compiler flags 
include make.inc.$(PLATFORM) 
BINDIR = $(HOME)/ecalj/bin

#-------------------------------------------------------
main   = ../maingw/
gwsrc  = ../srcgw/
maxloc = ../wanniergw/
moddir = MODDIR
module = -J$(moddir) -I$(moddir)
tags   = ../

progs  = hbasfp0 hvccfp0 hx0fp0 hsfp0 hef hqpe hqpe_qsgw qg4gw gwinit heftet hmergewv rdata4gw_v2 convgwin hx0fp0_sc hsfp0_sc hqpe_sc kino_input_test hecor eout eout2 hmagnon hhomogas
progsw= huumat hpsig hmaxloc huumat_MPI hpsig_MPI wanplot hwmatK hwmatK_MPI htbplot # hx0fp0_mlw 
progs2 = $(progs) $(progsw) $(tags)TAGS 
scriptw = genMLWF genMLWF_vw genMLWF2 genMLWFmod run_argloc RPAWannier* FLEX_interaction.py
script = cleargw* dqpu eps* gw* mkG* chipm_mat $(scriptw)


#### You can choose these options. all is default.

all : $(progs2) 
	[ -d $(moddir) ] || mkdir  $(moddir) 

clean:  
	 rm -f  $(progs) $(progsw)

install:  
	 cp  $(progs) $(progsw) $(BINDIR)

install2:  
	 cp  $(script) $(BINDIR)

cleanall:  
	 rm -f  $(progs2) $(main)*.o $(gwsrc)*.o  *.mod  $(moddir)/*.mod $(tote)*.o  $(maxloc)*.o

doxygen:
	cd $(tags);doxygen;cd ./latex;make
	echo 'fpgw/latex/refman.pdf generated'

dep:
	@echo  'This generate a call-caller data set for fpgw/'
	@echo  'HELP --> ../TOOLS/FparserTools/f_calltree.py --help, and read CallCaller.sh'
	@echo
	@echo  '--- Now generating a file 'callcaller.dat' ... Wait!!! It takes 1 minute or so!'
	@echo '        If you like to apply this to other programs, modify this script'
	@echo  ' NOTE: T.Kotani is not sure whether this is relaiable enough or not... let me know something wrong...'
	$(tags)/../TOOLS/FparserTools/f_calltree.py $(main)/*.F $(gwsrc)/*.F $(tote)/*.F >callcaller.dat 2>callcaller.err
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

# This is necesaly to compile *.f in right order.
# When you recompile and link, just repeat 'make' (not necessary to repeat 'make init').
# When checkmodule recompile source, you have to repeat 'make'.
init:
	rm moduledepends.*
moduledepends.inc: 
	exec ../../TOOLS/checkmodule $(gwsrc)*.F $(main)*.F $(maxloc)*.F

########################################################
MPI =  $(gwsrc)MPI_fpgw2.o

GWLIB =  \
 $(gwsrc)m_homoelectron.o \
 $(gwsrc)m_readwan.o \
 $(gwsrc)diagcv3.o \
 $(gwsrc)m_w0w0i.o \
 $(gwsrc)getwemax.o \
 $(gwsrc)genallcf_dump.o \
 $(gwsrc)wse.o \
 $(gwsrc)bzints2.o \
 $(gwsrc)wintzsg.o \
 $(gwsrc)gintxx.o \
 $(gwsrc)gwinput_v2.o \
 $(gwsrc)pplmat.o \
 $(gwsrc)rs.o \
 $(gwsrc)conv2gwinput.o \
 $(gwsrc)getbzdata1.o \
 $(gwsrc)getgv2.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt5.o \
 $(gwsrc)m_tetrakbt.o \
 $(gwsrc)m_tetwt.o \
 $(gwsrc)x0kf_v4h.o \
 $(gwsrc)cinvrx.o \
 $(gwsrc)zsvd.o \
 $(gwsrc)m_zmel.o \
 $(gwsrc)m_freq.o \
 $(gwsrc)m_hamindex.o\
 $(gwsrc)readpomat.o \
 $(gwsrc)keyvalue.o \
 $(gwsrc)rppovl.o \
 $(gwsrc)nocctotg.o \
 $(gwsrc)ppbafp.fal.o \
 $(gwsrc)psi2b_v2.o \
 $(gwsrc)psi2b_v3.o \
 $(gwsrc)wfacx.o \
 $(gwsrc)sortea.o \
 $(gwsrc)rydberg.o \
 $(gwsrc)polinta.o \
 $(gwsrc)efsimplef.o \
 $(gwsrc)extension.o \
 $(gwsrc)rangedq.o \
 $(gwsrc)nword.o \
 $(gwsrc)scg.o \
 $(gwsrc)matm.o \
 $(gwsrc)rdpp.o \
 $(gwsrc)mptauof.o \
 $(gwsrc)genallcf_mod.o \
 $(gwsrc)rgwinf_mod.o \
 $(gwsrc)rotdlmm.o \
 $(gwsrc)iopen.o \
 $(gwsrc)cputid.o \
 $(gwsrc)rw.o \
 $(gwsrc)ext.o \
 $(gwsrc)ext2.o \
 $(gwsrc)cross.o \
 $(gwsrc)mate.o \
 $(gwsrc)mate1.o \
 $(gwsrc)icopy.o \
 $(gwsrc)bib1.o \
 $(gwsrc)index.o \
 $(gwsrc)idxk.o \
 $(gwsrc)maxnn.o \
 $(gwsrc)reindx.o \
 $(gwsrc)iprint.o \
 $(gwsrc)bz.o \
 $(gwsrc)bzmesh.o \
 $(gwsrc)genqbz.o \
 $(gwsrc)switches.o \
 $(gwsrc)rwbzdata.o \
 $(gwsrc)llnew.o  \
 $(gwsrc)readeigen.o \
 $(gwsrc)readqg.o \
 $(gwsrc)iqindx.o  \
 $(gwsrc)alloclist.o \
 $(gwsrc)m_pkm4crpa.o \
 $(gwsrc)m_anf.o \
 $(gwsrc)qpe1.sc.o \
 $(gwsrc)icompvv2.o \
 $(gwsrc)iopenxx.o \
 $(gwsrc)qpe1.o \
 $(gwsrc)mopen.o \
 $(gwsrc)checksymlon.o \
 $(gwsrc)mkqg.o \
 $(gwsrc)m_q0p.o \
 $(gwsrc)q0irre.o \
 $(gwsrc)basnfp.o \
 $(gwsrc)excore.o \
 $(gwsrc)mkjp.o \
 $(gwsrc)strxq.o \
 $(gwsrc)sxcf_fal2.sc.o \
 $(gwsrc)sxcf_fal2.o \
$(gwsrc)amix.o \
$(gwsrc)dsifa.o \
$(gwsrc)dsisl.o \
$(gwsrc)dsidi.o \
$(gwsrc)diagcv2.o \
$(gwsrc)wronkj.o \
$(gwsrc)rxx.o \
$(gwsrc)hsmq.o \
$(gwsrc)u_lat_0.o \
$(gwsrc)mklegw.o \
$(gwsrc)bessl.o \
$(gwsrc)lgen.o \
$(gwsrc)hansr5.o \
$(gwsrc)hansr4.o \
$(gwsrc)lattc.o \
$(gwsrc)qdist.o \
$(gwsrc)dlmtor.o \
$(gwsrc)dpcopy.o \
$(gwsrc)dpadd.o \
$(gwsrc)dpzero.o \
$(gwsrc)ropyln.o \
$(gwsrc)ropcsm.o \
$(gwsrc)rpaq.o \
$(gwsrc)m_readeps.o \
$(gwsrc)m_readq0p.o \
$(gwsrc)m_readQGcou.o\
$(gwsrc)m_pbindex.o \
$(gwsrc)m_rotMPB.o \
$(gwsrc)m_iqindx_qtt.o \
$(gwsrc)rotwave.o \
$(gwsrc)rotcg.o \
$(gwsrc)rotgvec.o\
$(gwsrc)m_readVcoud.o\
$(gwsrc)m_itq.o\
$(gwsrc)nan.o\
$(gwsrc)m_readgwinput.o\
$(gwsrc)m_qbze.o\
$(gwsrc)m_eibz.o\
$(gwsrc)m_gaussianfilter.o\
$(gwsrc)m_llw.o\
$(gwsrc)m_readfreq_r.o\
$(gwsrc)m_readepswklm.o\
$(gwsrc)qpe1.qsgw.o


hecor:  $(main)hecor.o $(GWLIB)  $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
eout:   $(main)eout.o  $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
eout2: $(main)eout2.o  $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
gwinit:	$(main)gwinit.m.o $(MPI)  $(GWLIB) 
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
qg4gw:	$(main)qg4gw.m.o  $(MPI) $(GWLIB) 
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
rdata4gw_v2: $(main)rdata4gw_v2.m.o $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
hbasfp0: $(main)hbasfp0.m.o $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
hvccfp0: $(main)hvccfp0.m.o $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
hx0fp0:	 $(main)hx0fp0.m.o $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
hmagnon: $(main)hmagnon.m.o $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
hhomogas: $(main)hhomogas.m.o $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
hx0fp0_sc: $(main)hx0fp0.sc.m.o $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
hsfp0:	$(main)hsfp0.m.o $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
hsfp0_sc: $(main)hsfp0.sc.m.o $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
heftet:	$(main)heftet.m.o $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
hef:	$(main)hef.m.o $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
hqpe:	$(main)hqpe.m.o  $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
hqpe_sc:$(main)hqpe.sc.m.o  $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
hqpe_qsgw: $(main)hqpe_qsgw.m.o  $(GWLIB) $(MPI) 
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
hmergewv: $(main)hmergewv.m.o $(MPI) $(GWLIB) 
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
convgwin: $(main)convgwin.o $(MPI) $(GWLIB) 
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
kino_input_test: $(main)kino_input_test.o $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(GWLIB) $(MPI) $(LKFLAGS2) -o $@

##################################
MLOC = \
$(maxloc)maxloc0.o \
$(maxloc)maxloc1.o \
$(maxloc)maxloc2.o \
$(maxloc)maxloc3.o \
$(maxloc)maxloc4.o \
$(maxloc)m_read_Worb.o\
$(maxloc)m_DATA4GW.o\
$(maxloc)m_QG.o\
$(maxloc)m_FFT3D.o\
$(maxloc)mymath.o\
$(maxloc)myYlm.o\
$(maxloc)wfn2dx_2.o\
$(maxloc)m_MLWF.o\
$(maxloc)wmatK.o \
$(maxloc)wansub.o \
$(maxloc)cubeformat.o\
$(maxloc)xsfformat2.o\
$(maxloc)expand_mesh.o\
$(maxloc)wfn2dx_abc.o\
$(maxloc)wfnrho_abc.o\
$(maxloc)write_mesh2d.o\
 $(maxloc)RSMPI_mod.o \
 $(maxloc)efsimplef_RSMPI.o\
 $(maxloc)RSMPI_qkgroup_mod.o \
 $(maxloc)RSMPI_rotkindex_mod.o \
 $(maxloc)RSMPI_utils.o\
 $(maxloc)ent_MPI.o\
 $(maxloc)wmatK_mpi.o

hmaxloc: $(maxloc)hmaxloc.o $(MLOC) $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(MLOC) $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
hpsig: $(maxloc)hpsig.o $(MLOC) $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(MLOC) $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
hpsig_MPI: $(maxloc)hpsig_MPI.o $(MLOC) $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(MLOC) $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
huumat:	$(maxloc)huumat.o $(MLOC) $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(MLOC) $(GWLIB) $(MPI)  $(LKFLAGS2)  -o $@ 
huumat_MPI: $(maxloc)huumat_MPI.o $(MLOC) $(GWLIB) $(MPI) $(MPIMLOC)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(MLOC) $(GWLIB) $(MPI) $(MPIMLOC) $(LKFLAGS2) -o $@
hwmatK:	$(maxloc)hwmatK.o    $(MLOC) $(GWLIB) $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(MLOC) $(GWLIB) $(MPI) $(LKFLAGS2)  -o $@ 
hwmatK_MPI: $(maxloc)hwmatK_MPI.o $(MLOC) $(GWLIB) $(MPI) $(MPIMLOC)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(MLOC) $(GWLIB) $(MPI) $(MPIMLOC) $(LKFLAGS2)  -o $@ 
wanplot : $(maxloc)wanplot.o ${MLOC} ${GWLIB} $(MPI)
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(MLOC) $(GWLIB) $(MPI) $(LKFLAGS2)  -o $@ 
htbplot :  $(maxloc)htbplot.o
	@echo $(LK) $@.o 'link with *.o'
	@$(LK) $(LKFLAGS1) $< $(LKFLAGS2) -o $@

####################################################
$(tags)TAGS: $(progs)
	cd $(tags);etags ./*/*/*.F ./*/*.F

# --- Make rules ---
.SUFFIXES:
.SUFFIXES: .F .o
.F.o:
	$(FC) $(FFLAGS) $*.F -c -o $*.o

### fix 
$(gwsrc)wintzsg.o :  $(gwsrc)wintzsg.F
	$(FC)  $(FFLAGS) $(gwsrc)wintzsg.F  -c -o $*.o

$(maxloc)wanplot.o: $(maxloc)cubeformat.o $(maxloc)xsfformat2.o  

include moduledepends.inc


