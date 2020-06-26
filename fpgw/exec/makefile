### I think that you don't needs to modify this file. ###
### This file is not machine-dependent. #####
### Machine dependence in make.inc


# ---- Machine-specific compiler flags ---
#include make.inc.ifort_asahi_kino 
#include make.inc.thinkpad_gfortran_tkotani
PLATFORM=ifort
LIBMATH=/usr/lib/x86_64-linux-gnu/libfftw3.so.3 /usr/lib/liblapack.so.3gf /usr/lib/libblas.so.3gf

#PLATFORM=ifort
#LIBMATH=-mkl

include make.inc.$(PLATFORM)

BINDIR = $(HOME)/ecaljdev/bin

#-------------------------------------------------------
# src directories
main    = ../main/
gwsrc  = ../gwsrc/
tote = ../tote/
maxloc = ../Wannier/
tags   = ../

# tag directory
#
#progs  = hbasfp0 hvccfp0 hx0fp0 hsfp0 hef hqpe hchknw qg4gw gwinit heftet hmergewv hparainfo hbndout rdata4gw_v2 convgwin hx0fp0_sc hsfp0_sc hqpe_sc kino_input_test hecor eout eout2 h_uumatrix hsigmconv
# lmf_exec
#progs  = hbasfp0 hvccfp0 hx0fp0 hsfp0 hef hqpe hchknw qg4gw gwinit heftet hmergewv hparainfo hbndout rdata4gw_v2  hx0fp0_fal hx0fp1 

progs  = hbasfp0 hvccfp0 hx0fp0 hsfp0 hef hqpe hqpe_qsgw qg4gw gwinit heftet hmergewv rdata4gw_v2 convgwin hx0fp0_sc hsfp0_sc hqpe_sc kino_input_test hecor eout eout2 hmagnon hhomogas

progsw= huumat hpsig hmaxloc huumat_MPI hpsig_MPI wanplot hwmatK hwmatK_MPI htbplot # hx0fp0_mlw 

# progs  = hbasfp0 hvccfp0 hx0fp0 hsfp0 hef hqpe hchknw qg4gw gwinit heftet hmergewv hbndout rdata4gw_v2 convgwin hx0fp0_sc hsfp0_sc hqpe_sc kino_input_test hecor eout eout2 h_uumatrix hsigmconv hwmat hmaxloc huumat qpwf hpsig hnocc_mlw hx0fp0_mlw hphig

# hmaxloc1D
progs2 = $(progs) $(progsw) $(tags)TAGS 
#checkmod 

#script = cleargw* dqpu dtote eps* ex* gw* hqpemetal* inf* lmgw* plotg save* tote_lmfh2 xqp mkG*
scriptw = genMLWF genMLWF_vw genMLWF2 genMLWFmod run_argloc RPAWannier* FLEX_interaction.py

script = cleargw* dqpu eps* gw* mkG* chipm_mat $(scriptw)

#### You can choose these options. all is default.

all :$(progs2) 

clean:  
	 rm -f  $(progs)

install:  
	 cp  $(progs)  $(BINDIR)

install2:  
	 cp  $(script) $(BINDIR)

cleanall:  
	 rm -f  $(progs2) $(main)*.o $(gwsrc)*.o  *.mod  $(tote)*.o  $(maxloc)*.o

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
	rm -f $(main)time_hsfp0.sc.m.F   
	rm -f $(main)time_hx0fp0.sc.m.F  
	rm -f $(gwsrc)time_sxcf_fal2.sc.F
	rm -f $(gwsrc)time_rppovl.F      
	rm -f $(gwsrc)time_x0kf_v4h.F    
	rm -f $(gwsrc)time_ppbafp.fal.F  
	exec ../../TOOLS/checkmodule ../gwsrc/*.F ../main/*.F ../tote/*.F ../Wannier/*.F

checkmod:
	init

$(maxloc)wanplot.o: $(maxloc)cubeformat.o $(maxloc)xsfformat2.o  

OBJS2=\
$(maxloc)m_DATA4GW.o\
$(maxloc)m_QG.o\
$(maxloc)m_FFT3D.o\
$(maxloc)mymath.o\
$(maxloc)myYlm.o\
$(maxloc)wfn2dx_2.o\
$(maxloc)m_MLWF.o\
$(maxloc)cubeformat.o\
$(maxloc)xsfformat2.o\
$(maxloc)expand_mesh.o\
$(maxloc)wfn2dx_abc.o\
$(maxloc)wfnrho_abc.o\
$(maxloc)write_mesh2d.o 


#../../lm7K/subs/m_hamindex.F
# m_hamindex

## tete (total energy)  ###############################
## these are experimental code ###
ECOR = \
$(tote)hecor.o 

EO= \
$(tote)eout.o \

EO2= \
$(tote)eout2.o 

hecor: $(ECOR) $(GWLIB)   $(MPI)  $(GWLIB) $(COMM)
	$(LK) $(LKFLAGS1) $(ECOR) $(GWLIB) $(MPI)  $(COMM) $(LKFLAGS2) -o $@

eout: $(EO)   $(GWLIB) $(MPI)
	$(LK) $(LKFLAGS1) $(EO) $(GWLIB) $(MPI) $(LKFLAGS2) -o $@

eout2: $(EO2)  $(GWLIB) $(MPI)
	$(LK) $(LKFLAGS1) $(EO2) $(GWLIB) $(MPI) $(LKFLAGS2) -o $@
########################################################


 kino_input_test =  \
 $(main)kino_input_test.o 

 convg =  \
 $(main)convgwin.o 

 GWINIT =  \
 $(main)gwinit.m.o 

 QG =  \
 $(main)qg4gw.m.o 

 RDAT_v2 = \
 $(main)rdata4gw_v2.m.o 

 BAS = \
 $(main)hbasfp0.m.o

 VCC= \
 $(main)hvccfp0.m.o 

 SXC_SC = \
 $(main)hsfp0.sc.m.o

 SXC = \
 $(main)hsfp0.m.o 

 heftet = \
 $(main)heftet.m.o 

 # hnocc_mlw = \
 # $(maxloc)hnocc_mlw.o 

 hef = \
 $(main)hef.m.o 

 X0_SC = \
 $(main)hx0fp0.sc.m.o 

 X0 = \
 $(main)hx0fp0.m.o

 X0oku = \
 $(gwsrc)m_homoelectron.o \
 $(gwsrc)m_readwan.o \
 $(gwsrc)diagcv3.o \
 $(main)hmagnon.m.o

 X0homo = \
 $(gwsrc)m_homoelectron.o \
 $(gwsrc)m_readwan.o \
 $(gwsrc)diagcv3.o \
 $(main)hhomogas.m.o

 MPI =  $(gwsrc)MPI_fpgw2.o

 GWLIB =  \
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
 $(gwsrc)tetwt5$(tet5_g) \
 $(gwsrc)m_tetrakbt.o \
 $(gwsrc)m_tetwt.o \
 $(gwsrc)x0kf_v4h$(x0kf_g) \
 $(gwsrc)cinvrx.o \
 $(gwsrc)zsvd.o \
 $(gwsrc)m_zmel.o \
 $(gwsrc)m_freq.o \
 $(gwsrc)m_hamindex.o\
 $(gwsrc)readpomat.o \
 $(gwsrc)keyvalue.o \
 $(gwsrc)rppovl.o \
 $(gwsrc)nocctotg.o \
 $(gwsrc)ppbafp.fal$(para_g) \
 $(gwsrc)psi2b_v2$(para_g) \
 $(gwsrc)psi2b_v3$(para_g) \
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
 $(gwsrc)sxcf_fal2.sc$(sxcf_g) \
 $(gwsrc)sxcf_fal2$(sxcf_g) \
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
$(gwsrc)m_eibz.o


 QPE_QSGW = \
 $(main)hqpe_qsgw.m.o\
 $(gwsrc)qpe1.qsgw.o

 QPE_SC = \
 $(main)hqpe.sc.m$(hqpe_g)

 QPE = \
 $(main)hqpe.m$(hqpe_g)

 MERGE = \
  $(main)hmergewv.m.o


############################################

 gwinit:	$(GWINIT)  $(MPI)  $(GWLIB) 
	 $(LK) $(LKFLAGS1) $(GWINIT) $(MPI) $(GWLIB)  $(LKFLAGS2) -o $@

 # qpwf:		$(maxloc)qpwf.o $(GWLIB) $(MPI)  $(COMM) 
 # 	$(LK) $(LKFLAGS1) $(maxloc)qpwf.o $(GWLIB) $(MPI) $(COMM)  $(LKFLAGS2) -o $@

 qg4gw:		$(QG)  $(MPI) $(GWLIB) $(COMM)  
	 $(LK) $(LKFLAGS1) $(QG) $(MPI) $(GWLIB) $(COMM)    $(LKFLAGS2) -o $@

 rdata4gw_v2:	$(RDAT_v2)  $(MPI)  $(COMM) $(GWLIB)
	 $(LK) $(LKFLAGS1) $(RDAT_v2) $(MPI) $(COMM)  $(GWLIB) $(LKFLAGS2) -o $@

 hbasfp0:	$(BAS)  $(MPI)  $(COMM) $(GWLIB)
	 $(LK) $(LKFLAGS1) $(BAS) $(MPI)  $(COMM)  $(GWLIB) $(LKFLAGS2) -o $@

 hvccfp0:	$(MPI) $(VCC)   $(DERFC) $(MPI)  $(COMM) $(GWLIB)
	 $(LK) $(LKFLAGS1) $(VCC) $(DERFC) $(MPI)  $(COMM) $(GWLIB) $(LKFLAGS2) -o $@

 hx0fp0:	$(MPI) $(X0) $(GWLIB) $(MPI)  $(COMM) 
	 $(LK) $(LKFLAGS1) $(X0)     $(GWLIB) $(MPI)  $(COMM) $(LKFLAGS2) -o $@

 hmagnon:	$(MPI) $(GWLIB) $(X0oku) $(MPI)  $(COMM) 
	 $(LK) $(LKFLAGS1) $(GWLIB) $(X0oku) $(MPI)  $(COMM) $(LKFLAGS2) -o $@

 hhomogas:	$(MPI) $(GWLIB) $(X0homo) $(MPI)  $(COMM) 
	 $(LK) $(LKFLAGS1) $(GWLIB) $(X0homo) $(MPI)  $(COMM) $(LKFLAGS2) -o $@

 hx0fp0_sc:	$(MPI) $(X0_SC) $(GWLIB)   $(MPI)  $(COMM) 
	 $(LK) $(LKFLAGS1) $(X0_SC)     $(GWLIB) $(MPI)  $(COMM)  $(LKFLAGS2) -o $@

 hsfp0:		$(MPI) $(SXC) $(GWLIB)   $(MPI)  $(COMM) 
	 $(LK) $(LKFLAGS1) $(SXC)    $(GWLIB) $(MPI)  $(COMM)  $(LKFLAGS2) -o $@

 hsfp0_sc:	$(MPI) $(SXC_SC) $(GWLIB)   $(MPI)  $(COMM) 
	 $(LK) $(LKFLAGS1) $(SXC_SC)    $(GWLIB)  $(MPI)  $(COMM) $(LKFLAGS2) -o $@

 heftet:	$(heftet) $(GWLIB)  $(MPI) $(MPI)  $(COMM) 
	 $(LK) $(LKFLAGS1) $(heftet) $(GWLIB)  $(MPI)  $(COMM) $(LKFLAGS2) -o $@

 hef:		$(hef) $(GWLIB)   $(MPI) $(COMM) 
	 $(LK) $(LKFLAGS1) $(hef)    $(GWLIB)  $(MPI)  $(COMM) $(LKFLAGS2) -o $@

 hqpe:		$(QPE) $(MPI) $(COMM) $(GWLIB)
	 $(LK) $(LKFLAGS1) $(QPE) $(MPI) $(COMM) $(GWLIB) $(LKFLAGS2) -o $@

 hqpe_sc:		$(QPE_SC)  $(MPI) $(COMM) $(GWLIB)
	 $(LK) $(LKFLAGS1) $(QPE_SC) $(GWLIB) $(MPI)  $(COMM) $(LKFLAGS2) -o $@

 hqpe_qsgw:		$(QPE_QSGW)  $(GWLIB) $(MPI) $(COMM) 
	 $(LK) $(LKFLAGS1) $(QPE_QSGW) $(GWLIB) $(MPI)  $(COMM) $(LKFLAGS2) -o $@

 hmergewv:	$(MERGE) $(MPI) $(GWLIB)  $(COMM) 
	 $(LK) $(LKFLAGS1) $(MERGE) $(GWLIB) $(MPI)  $(COMM)  $(LKFLAGS2) -o $@

 convgwin:	$(convg) 
	 $(LK) $(LKFLAGS1) $(convg) $(GWLIB) $(MPI) $(LKFLAGS2) -o $@

 kino_input_test:	$(kino_input_test) $(GWLIB) $(MPI)  $(COMM) 
	 $(LK) $(LKFLAGS1) $(kino_input_test) $(GWLIB)  $(MPI)  $(COMM) $(LKFLAGS2) -o $@

 MLOC = \
 $(maxloc)hmaxloc.o \
 $(maxloc)maxloc0.o \
 $(maxloc)maxloc1.o \
 $(maxloc)maxloc2.o \
 $(maxloc)maxloc3.o \
 $(maxloc)maxloc4.o \
 $(maxloc)m_read_Worb.o

 PHIG = \
 $(maxloc)hphig.o 

 PSIG = \
 $(maxloc)hpsig.o \
 $(maxloc)m_read_Worb.o

 PSIG_MPI = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_qkgroup_mod.o \
 $(maxloc)RSMPI_rotkindex_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)hpsig_MPI.o \
 $(maxloc)m_read_Worb.o

 # PSIPSI = \
 # $(gwt)hpsipsi.o \
 # $(gwsrc)gintxx.o \
 # $(gwsrc)mopen.o \
 # $(gwsrc)u_lat_0.o \
 # $(gwsrc)wronkj.o \
 # $(gwsrc)mklegw.o \
 # $(gwsrc)rxx.o \
 # $(gwsrc)bessl.o \
 # $(gwsrc)cross.o 

 UU2 = \
 $(maxloc)huumat.o 

 UU2_MPI = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)huumat_MPI.o \
 $(maxloc)ent_MPI.o

 WMATK = \
 $(maxloc)hwmatK.o \
 $(maxloc)wmatK.o \
 $(maxloc)maxloc0.o

 WMATK_MPI = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_qkgroup_mod.o \
 $(maxloc)RSMPI_rotkindex_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)hwmatK_MPI.o \
 $(maxloc)wmatK_mpi.o \
 $(maxloc)maxloc0.o \
 $(maxloc)efsimplef_RSMPI.o


 hmaxloc:	$(MLOC)  $(GWLIB) $(MPI)
	$(LK) $(LKFLAGS1) $(MLOC) $(GWLIB) $(MPI) $(LKFLAGS2) -o $@

 hpsig: $(PSIG) $(GWLIB) $(MPI)
	$(LK) $(LKFLAGS1) $(PSIG)     $(GWLIB) $(MPI)   $(LKFLAGS2) -o $@

 hpsig_MPI:	$(PSIG_MPI) $(GWLIB) $(MPI)
	$(LK) $(LKFLAGS1) $(PSIG_MPI)	$(GWLIB) $(MPI)	$(LKFLAGS2) -o $@

 huumat:	$(UU2) $(GWLIB) $(MPI)
	$(LK) -o $@ $(LKFLAGS1) $(UU2)     $(GWLIB) $(MPI)  $(LKFLAGS2) 

 huumat_MPI:	$(UU2_MPI) $(GWLIB) $(MPI)
	$(LK) $(LKFLAGS1) $(UU2_MPI)     $(GWLIB) $(MPI)  $(LKFLAGS2) -o $@

 hwmatK:		$(WMATK) $(GWLIB) $(MPI)
	$(LK) $(LKFLAGS1) $(WMATK) $(COMMX)   $(GWLIB) $(MPI) $(LKFLAGS2) -o $@

 hwmatK_MPI:		$(WMATK_MPI) $(GWLIB) $(MPI)
	$(LK) $(LKFLAGS1) $(WMATK_MPI) $(COMMX)   $(GWLIB) $(MPI) $(LKFLAGS2) -o $@


#huumat.o: $(gwsrc)keyvalue.o $(gwsrc)readqg.o  $(gwsrc)readeigen.o $(gwsrc)rwbzdata.o  \
#          $(gwsrc)genallcf_mod.o 

wanplot : ${OBJS2} $(maxloc)wanplot.o ${GWLIB} $(MPI)
	$(LK) $(LKFLAGS1) $(maxloc)wanplot.o $(OBJS2) $(GWLIB) $(MPI) $(COMM) $(LKFLAGS2) -o $@

htbplot :  $(maxloc)htbplot.o
	$(LK) $(LKFLAGS1) $(maxloc)htbplot.o $(LKFLAGS2) -o $@

################################ test
#
# test_genallcf:	$(test_genallcf) 
#	 $(LK) $(LKFLAGS1) $(test_genallcf) $(LKFLAGS2) -o $@

 $(tags)TAGS: $(progs)
	cd $(tags);etags ./*/*/*.F ./*/*.F

# --- Make rules ---
.SUFFIXES:
.SUFFIXES: .F .o
#.SUFFIXES: .f .o .c1_o .c2_0 .c3_o .c4_o .F

.F.o:
	$(FC) $(FFLAGS) $*.F -c -o $*.o
#	etags $*.f -o $(tags)`echo $*.f| sed 's/..\///' | sed 's/\//-/g'`.tags

#.F.o:
#	$(FC) $(FFLAGS) $*.F -c -o $*.o
#	etags $*.f -o $(tags)`echo $*.f| sed 's/..\///' | sed 's/\//-/g'`.tags

#.f.o:
#	$(FC) $(FFLAGS) $*.f -c -o $*.o
#	etags $*.f -o $(tags)`echo $*.f| sed 's/..\///' | sed 's/\//-/g'`.tags

.f.c1_o:
	$(FC) $(FFLAGS_c1) $*.f -c -o $*.c1_o
	etags $*.f -o $(tags)`echo $*.f| sed 's/..\///' | sed 's/\//-/g'`.tags

.f.c2_o:
	$(FC) $(FFLAGS_c2) $*.f -c -o $*.c2_o
	etags $*.f -o $(tags)`echo $*.f| sed 's/..\///' | sed 's/\//-/g'`.tags

.f.c3_o:
	$(FC) $(FFLAGS_c3) $*.f -c -o $*.c3_o
	etags $*.f -o $(tags)`echo $*.f| sed 's/..\///' | sed 's/\//-/g'`.tags

.f.c4_o:
	$(FC) $(FFLAGS_c4) $*.f -c -o $*.c4_o
	etags $*.f -o $(tags)`echo $*.f| sed 's/..\///' | sed 's/\//-/g'`.tags

$(gwsrc)wintzsg.o :  $(gwsrc)wintzsg.F
	$(FC)  $(FFLAGS) $(gwsrc)wintzsg.F  -c -o $*.o

#check:
#	(cd ../TESTinstallGW;./testgw.py --enforce --all) 

# test for f90 dependency
#../main/hvccfp0.m.o	:	../main/hx0fp0.m.o
#
#../main/hvccfp0.m.o	:	../main/hbasfp0.m.o

include moduledepends.inc





# ################################################################################
# ##### You can comment out these blocks to commnet out memory and time check (verbose output)
# addtime=script/addtime.awk 
# septhen=script/then_separate.awk
# alloclist=script/add_alloclist.awk
# $(main)hsfp0.sc.m.o:  $(main)hsfp0.sc.m.F
# 	gawk -f $(addtime) -vSTART=1 $(main)hsfp0.sc.m.F | gawk -f $(septhen) | gawk -f $(alloclist) > $(main)time_hsfp0.sc.m.F
# 	$(FC) $(FFLAGS) $(main)time_hsfp0.sc.m.F -c -o $*.o

# $(main)hx0fp0.sc.m.o: $(main)hx0fp0.sc.m.F
# 	gawk -f $(addtime) -vSTART=1 $(main)hx0fp0.sc.m.F | gawk -f $(septhen) | gawk -f $(alloclist) > $(main)time_hx0fp0.sc.m.F
# 	$(FC) $(FFLAGS) $(main)time_hx0fp0.sc.m.F -c -o $*.o

# $(gwsrc)sxcf_fal2.sc$(sxcf_g): $(gwsrc)sxcf_fal2.sc.F
# 	gawk -f $(addtime) -vSTART=100 $(gwsrc)sxcf_fal2.sc.F | gawk -f $(septhen) | gawk -f $(alloclist) > $(gwsrc)time_sxcf_fal2.sc.F
# 	$(FC) $(FFLAGS) $(gwsrc)time_sxcf_fal2.sc.F -c -o $*.o

# #$(gwsrc)rppovl.o: $(gwsrc)rppovl.F
# #	gawk -f $(addtime) -vSTART=200 $(gwsrc)rppovl.F | gawk -f $(septhen) | gawk -f $(alloclist) > $(gwsrc)time_rppovl.F
# #	$(FC) $(FFLAGS) $(gwsrc)time_rppovl.F -c -o $*.o

# $(gwsrc)x0kf_v4h$(x0kf_g): $(gwsrc)x0kf_v4h.F 
# 	gawk -f $(addtime) -vSTART=100 $(gwsrc)x0kf_v4h.F | gawk -f $(septhen) | gawk -f $(alloclist) > $(gwsrc)time_x0kf_v4h.F
# 	$(FC) $(FFLAGS) $(gwsrc)time_x0kf_v4h.F -c -o $*.o

# $(gwsrc)ppbafp.fal$(para_g): $(gwsrc)ppbafp.fal.F
# 	gawk -f $(addtime) -vSTART=300 $(gwsrc)ppbafp.fal.F | gawk -f $(septhen) | gawk -f $(alloclist) > $(gwsrc)time_ppbafp.fal.F
# 	$(FC) $(FFLAGS) $(gwsrc)time_ppbafp.fal.F -c -o $*.o
# #$(gwsrc)ppbafp.fal$(para_g): $(gwsrc)ppbafp.fal.F
# #	gawk -f $(addtime) -vSTART=300 $(gwsrc)ppbafp.fal.F | gawk -f $(septhen) | gawk -f $(alloclist) > $(gwsrc)time_ppbafp.fal.F
# #	$(FC) $(FFLAGS) $(gwsrc)time_ppbafp.fal.F -c -o $*.o
# ###################################################################################
