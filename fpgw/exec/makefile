### I think that you don't needs to modify this file. ###
### This file is not machine-dependent. #####
### Machine dependence in make.inc


# ---- Machine-specific compiler flags ---
#include make.inc.ifort_asahi_kino 
#include make.inc.thinkpad_gfortran_tkotani
include make.inc.sora1



#-------------------------------------------------------
# src directories
gwsrc   = ../gwsrc/
main    = ../main/
nfpsrc  = ../nfpsrc/
slatsmlib  = ../slatsmlib/
tote = ../tote/
maxloc = ../Miyake/maxloc/
# tag directory
tags   = ../

#progs  = hbasfp0 hvccfp0 hx0fp0 hsfp0 hef hqpe hchknw qg4gw gwinit heftet hmergewv hparainfo hbndout rdata4gw_v2 convgwin hx0fp0_sc hsfp0_sc hqpe_sc kino_input_test hecor eout eout2 h_uumatrix hsigmconv
# lmf_exec
#progs  = hbasfp0 hvccfp0 hx0fp0 hsfp0 hef hqpe hchknw qg4gw gwinit heftet hmergewv hparainfo hbndout rdata4gw_v2  hx0fp0_fal hx0fp1 

progs  = hbasfp0 hvccfp0 hx0fp0 hsfp0 hef hqpe hchknw qg4gw gwinit heftet hmergewv hparainfo hbndout rdata4gw_v2 convgwin hx0fp0_sc hsfp0_sc hqpe_sc kino_input_test hecor eout eout2 h_uumatrix hsigmconv hwmat hmaxloc huumat qpwf hpsig hnocc_mlw hx0fp0_mlw hphig
# hmaxloc1D
progs2 = $(progs) $(tags)TAGS 
#checkmod 

#script = cleargw* dqpu dtote eps* ex* gw* hqpemetal* inf* lmgw* plotg save* tote_lmfh2 xqp mkG*
script = cleargw* dqpu eps* gw* mkG* 

#### You can choose these options. all is default.

all :$(progs2) 

clean:  
	 rm -f  $(progs)

install:  
	 cp  $(progs)  $(BINDIR)

install2:  
	 cp  $(script) $(BINDIR)

cleanall:  
	 rm -f  $(progs2) $(main)*.o $(gwsrc)*.o  $(nfpsrc)*.o  *.mod  $(slatsmlib)*.o $(tote)*.o  $(maxloc)*.o 

# This is necesaly to compile *.f in right order.
# When you recompile and link, just repeat 'make' (not necessary to repeat 'make init').
# When checkmodule recompile source, you have to repeat 'make'.
init:
	exec ./checkmodule ../*/*.F ../Miyake/*/*.F

checkmod:
	exec ./checkmodule ../*/*.F ../Miyake/*/*.F 
#../../lm7K/subs/m_hamindex.F


####from tote #################################################################################                                     
#LIBLOC  = $(ECAL)/fftw/libfftw.a $(LIBMATH)
##-L/usr/local/ATLAS/lib/Linux_P4SSE2 -llapack -lcblas -lf77blas -latlas                                                             
#LIBSLA  = $(ECAL)/slatsm/slatsm.a
#LIBFP   = $(ECAL)/lm-6.14y/fp/subs.a
#LIBSUBS = $(ECAL)/lm-6.14y/subs/subs.a
#LIBES  = $(LIBSLA) $(LIBLOC)
#lmsrc   = ../../lm-6.14y/
######################################################################################            


ECOR = \
$(tote)hecor.o \
$(tote)rpaq.o 
#eispack.o

NFPLtot  = $(nfpsrc)diagcv2.o 

GW0tot = \
$(gwsrc)rwbzdata.o \
$(gwsrc)keyvalue.o \
$(gwsrc)genallcf_mod.o \
$(gwsrc)rgwinf_mod.o \
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
$(gwsrc)linpackdummy.o \
$(gwsrc)rppovl.o \
$(gwsrc)llnew.o 

#LMFtot= \
#$(lmsrc)lmf.o \
#$(tote)mkpot_exec.o \
#$(tote)locpot_exec.o \
#$(tote)vxcnsp_exec.o \
#$(tote)evxcv_exec.o \
#$(tote)mkehkf_exec.o \
#$(tote)smvxcm_exec.o

EO= \
$(tote)eout.o \
$(gwsrc)rydberg.o

EO2= \
$(tote)eout2.o \
$(gwsrc)rydberg.o

hecor: $(ECOR) $(NFPLtot) $(GW0)   $(MPI) 
	$(LK) $(LKFLAGS1) $(ECOR) $(GW0) $(NFPLtot)  $(MPI) $(LKFLAGS2) -o $@

#	$(LK) $(LKFLAGS1) $(ECORtot) $(GW0tot) $(NFPLtot) $(LKFLAGS2) -o $@

eout: $(EO)   
	$(LK) $(LKFLAGS1) $(EO) $(LKFLAGS2) -o $@

eout2: $(EO2)  
	$(LK) $(LKFLAGS1) $(EO2) $(LKFLAGS2) -o $@

#lmf_exec:	$(LMFtot)
#	$(LK) $(LKC) $(LMFtot) $(LIBFP) $(LIBSUBS) $(LIBES) -o $@
##############################################################################








################
#
# BNDCONN= \
# $(gwsrc)bndconn.o   ### This is not linked but bndconn.o is used in lm/lmfgw. 
# It is now included in lm/gw/
 DERFC=   $(nfpsrc)derfc.o \
          $(nfpsrc)d1mach.o \
          $(nfpsrc)i1mach.o 

# test_genallcf =  \
# $(main)test_genallcf.o \
# $(gwsrc)genallcf_dump.o \
# $(GW0)


 kino_input_test =  \
 $(main)kino_input_test.o 

 convg =  \
 $(main)convgwin.o 

 GWINIT =  \
 $(main)gwinit.m.o \
 $(gwsrc)cross.o \
 $(gwsrc)genqbz.o \
 $(gwsrc)checksymlon.o \
 $(gwsrc)bzmesh.o \
 $(gwsrc)rangedq.o \
 $(gwsrc)iopenxx.o \
 $(gwsrc)iprint.o \
 $(gwsrc)keyvalue.o \
 $(gwsrc)switches.o \
 $(gwsrc)iopen.o 


 QG =  \
 $(gwsrc)conv2gwinput.o \
 $(main)qg4gw.m.o \
 $(gwsrc)getbzdata1.o \
 $(gwsrc)mkqg.o \
 $(gwsrc)q0irre.o \
 $(gwsrc)getgv2.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)tetwt5$(tet5_g) \
 $(gwsrc)zsvd.o \
 $(GW0)

 RDAT_v2 = \
 $(gwsrc)keyvalue.o \
 $(gwsrc)switches.o \
 $(main)rdata4gw_v2.m.o \
 $(gwsrc)rwbzdata.o \
 $(gwsrc)mopen.o \
 $(gwsrc)gintxx.o \
 $(gwsrc)cinvrx.o \
 $(gwsrc)idxk.o \
 $(gwsrc)nword.o \
 $(gwsrc)gwinput_v2.o \
 $(gwsrc)matm.o \
 $(gwsrc)getgv2.o \
 $(gwsrc)iopen.o \
 $(gwsrc)pplmat.o \
 $(gwsrc)bzmesh.o \
 $(gwsrc)ext.o \
 $(gwsrc)ext2.o \
 $(gwsrc)cross.o \
 $(gwsrc)rs.o \
 $(gwsrc)extension.o \
 $(gwsrc)rangedq.o \
 $(gwsrc)llnew.o \
 $(gwsrc)iqindx.o \
 $(gwsrc)polinta.o


 BAS = \
 $(main)hbasfp0.m.o \
 $(gwsrc)reindx.o \
 $(gwsrc)maxnn.o \
 $(gwsrc)icopy.o \
 $(gwsrc)basnfp.o \
 $(gwsrc)rgwinf_mod.o \
 $(gwsrc)keyvalue.o \
 $(gwsrc)switches.o \
 $(gwsrc)gintxx.o \
 $(gwsrc)rs.o \
 $(gwsrc)ext.o \
 $(gwsrc)iopen.o \
 $(gwsrc)excore.o \
 $(gwsrc)rydberg.o \
 $(gwsrc)extension.o \
 $(gwsrc)rangedq.o \
 $(gwsrc)polinta.o \
 $(gwsrc)llnew.o


 VCC= \
 $(main)hvccfp0.m.o \
 $(gwsrc)mkjp.o \
 $(gwsrc)gintxx.o \
 $(gwsrc)extension.o \
 $(gwsrc)rangedq.o \
 $(gwsrc)keyvalue.o \
 $(gwsrc)switches.o \
 $(gwsrc)strxq.o \
 $(gwsrc)iopen.o \
 $(gwsrc)pplmat.o \
 $(gwsrc)matm.o \
 $(gwsrc)getgv2.o \
 $(gwsrc)mopen.o \
 $(gwsrc)cross.o \
 $(gwsrc)llnew.o \
 $(gwsrc)readqg.o \
 $(gwsrc)iqindx.o \
 $(gwsrc)cputid.o

 SXC_SC = \
 $(main)hsfp0.sc.m.o \
 $(gwsrc)wse.o \
 $(gwsrc)sxcf_fal2$(sxcf_g) \
 $(gwsrc)sxcf_fal2.sc$(sxcf_g) \
 $(gwsrc)bzints2.o \
 $(gwsrc)wintzsg.o \
 $(nfpsrc)diagcv2.o

 SXC = \
 $(main)hsfp0.m.o \
 $(gwsrc)wse.o \
 $(gwsrc)wintzsg.o \
 $(gwsrc)sxcf_fal2$(sxcf_g) \
 $(gwsrc)bzints2.o \
 $(gwsrc)genallcf_dump.o \
 $(nfpsrc)diagcv2.o

 WMAT = \
 $(maxloc)hwmat.o \
 $(maxloc)maxloc0.o \
 $(gwsrc)wse.o \
 $(maxloc)wmat.o \
 $(gwsrc)genallcf_dump.o

 MLOC = \
 $(maxloc)hmaxloc.o \
 $(maxloc)maxloc0.o \
 $(maxloc)maxloc1.o \
 $(maxloc)maxloc2.o \
 $(maxloc)maxloc3.o \
 $(gwsrc)wse.o \
 $(gwsrc)genallcf_dump.o

 MLOC1D = \
 $(maxloc)hmaxloc1D.o \
 $(maxloc)maxloc0.o \
 $(maxloc)maxloc1.o \
 $(maxloc)maxloc2.o \
 $(maxloc)maxloc3.o \
 $(gwsrc)wse.o \
 $(gwsrc)genallcf_dump.o

 heftet = \
 $(main)heftet.m.o \
 $(gwsrc)bzints2.o

 hnocc_mlw = \
 $(maxloc)hnocc_mlw.o \
 $(gwsrc)bzints2.o

 hef = \
 $(main)hef.m.o \
 $(gwsrc)wse.o

 CHK = \
 $(main)hchknw.m.o \
 $(gwsrc)genallcf_dump.o

 X0_SC = \
 $(main)hx0fp0.sc.m.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5$(tet5_g) \
 $(gwsrc)x0kf_v2h$(x0kf_g) \
 $(gwsrc)x0kf_v3h$(x0kf_g) \
 $(gwsrc)x0kf_v4h$(x0kf_g) \
 $(nfpsrc)diagcv2.o

 X0 = \
 $(main)hx0fp0.m.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5$(tet5_g) \
 $(gwsrc)x0kf_v2h$(x0kf_g) \
 $(gwsrc)x0kf_v3h$(x0kf_g) \
 $(gwsrc)x0kf_v4h$(x0kf_g) \
 $(nfpsrc)diagcv2.o \
 $(tote)rpaq.o \
 $(gwsrc)cinvrx.o \
 $(gwsrc)zsvd.o

 X0mlw = \
 $(maxloc)hx0fp0.m.o \
 $(maxloc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5$(tet5_g) \
 $(gwsrc)x0kf_v2h$(x0kf_g)\
 $(maxloc)x0kf_v3h$(x0kf_g) \
 $(nfpsrc)diagcv2.o \
 $(tote)rpaq.o \
 $(gwsrc)cinvrx.o

 UU = \
 $(main)h_uumatrix.m.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5$(tet5_g) \
 $(gwsrc)gintxx.o \
 $(gwsrc)pplmat.o \
 $(gwsrc)getgv2.o \
 $(gwsrc)x0kf_v2h$(x0kf_g) \
 $(gwsrc)x0kf_v3h$(x0kf_g) \
 $(gwsrc)x0kf_v4h$(x0kf_g) \
 $(gwsrc)rs.o \
 $(gwsrc)mopen.o \
 $(nfpsrc)u_lat_0.o \
 $(nfpsrc)wronkj.o \
 $(nfpsrc)mklegw.o \
 $(nfpsrc)rxx.o \
 $(nfpsrc)bessl.o \
 $(nfpsrc)cross.o \
 $(nfpsrc)diagcv2.o

 UU2 = \
 $(maxloc)huumat.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5$(tet5_g) \
 $(gwsrc)gintxx.o \
 $(gwsrc)pplmat.o \
 $(gwsrc)getgv2.o \
 $(gwsrc)x0kf_v2h$(x0kf_g) \
 $(gwsrc)rs.o \
 $(gwsrc)mopen.o \
 $(nfpsrc)u_lat_0.o \
 $(nfpsrc)wronkj.o \
 $(nfpsrc)mklegw.o \
 $(nfpsrc)rxx.o \
 $(nfpsrc)bessl.o \
 $(nfpsrc)cross.o

 PSIG = \
 $(maxloc)hpsig.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5$(tet5_g) \
 $(gwsrc)gintxx.o \
 $(gwsrc)pplmat.o \
 $(gwsrc)getgv2.o \
 $(gwsrc)x0kf_v2h$(x0kf_g) \
 $(gwsrc)rs.o \
 $(gwsrc)mopen.o \
 $(nfpsrc)u_lat_0.o \
 $(nfpsrc)wronkj.o \
 $(nfpsrc)mklegw.o \
 $(nfpsrc)rxx.o \
 $(nfpsrc)bessl.o \
 $(nfpsrc)cross.o

 PHIG = \
 $(maxloc)hphig.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5$(tet5_g) \
 $(gwsrc)gintxx.o \
 $(gwsrc)pplmat.o \
 $(gwsrc)getgv2.o \
 $(gwsrc)x0kf_v2h$(x0kf_g) \
 $(gwsrc)rs.o \
 $(gwsrc)mopen.o \
 $(nfpsrc)u_lat_0.o \
 $(nfpsrc)wronkj.o \
 $(nfpsrc)mklegw.o \
 $(nfpsrc)rxx.o \
 $(nfpsrc)bessl.o \
 $(nfpsrc)cross.o

 MPI = $(gwsrc)MPI_fpgw.o

 GW0 = \
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
 $(gwsrc)linpackdummy.o \
 $(gwsrc)switches.o \
 $(gwsrc)rwbzdata.o \
 $(gwsrc)llnew.o  \
 $(gwsrc)readeigen.o \
 $(gwsrc)readqg.o \
 $(gwsrc)iqindx.o  \
 $(gwsrc)alloclist.o 

 QPE_SC = \
 $(gwsrc)switches.o \
 $(gwsrc)keyvalue.o \
 $(main)hqpe.sc.m$(hqpe_g) \
 $(gwsrc)qpe1.sc.o \
 $(gwsrc)icompvv2.o \
 $(gwsrc)iopenxx.o \
 $(gwsrc)iopen.o \
 $(gwsrc)rw.o \
 $(gwsrc)rydberg.o \
 $(gwsrc)iprint.o \
 $(gwsrc)rwbzdata.o \
 $(slatsmlib)dmcpy.o \
 $(slatsmlib)dsifa.o \
 $(slatsmlib)dsisl.o \
 $(slatsmlib)dsidi.o \
 $(slatsmlib)amix.o 

# ../../slatsm/slatsm.a 



 QPE = \
 $(gwsrc)switches.o \
 $(gwsrc)keyvalue.o \
 $(main)hqpe.m$(hqpe_g) \
 $(gwsrc)qpe1.o \
 $(gwsrc)icompvv2.o \
 $(gwsrc)iopenxx.o \
 $(gwsrc)iopen.o \
 $(gwsrc)rw.o \
 $(gwsrc)rydberg.o

 MERGE = \
 $(main)hmergewv.m.o \
 $(gwsrc)switches.o \
 $(gwsrc)keyvalue.o \
 $(gwsrc)iopen.o

 PARAINFO = \
 $(main)hparainfo.m.o \
 $(gwsrc)charext.o


 BNDOUT = \
 $(main)hbndout.m.o \
 $(gwsrc)iqagree.o \
 $(gwsrc)iopenxx.o \
 $(gwsrc)iopen.o \
 $(gwsrc)polinta.o \
 $(gwsrc)rydberg.o \
 $(gwsrc)extension.o \
 $(gwsrc)rangedq.o \
 $(gwsrc)switches.o \
 $(gwsrc)keyvalue.o 


 NFPL  = $(nfpsrc)wronkj.o \
         $(nfpsrc)sylm.o \
         $(nfpsrc)sylmnc.o \
         $(nfpsrc)u_lat_0.o \
         $(nfpsrc)mklegw.o \
         $(nfpsrc)cross.o \
         $(nfpsrc)setpr.o \
         $(nfpsrc)bessl.o \
         $(nfpsrc)rxx.o \
         $(nfpsrc)hsmq.o \
         $(nfpsrc)lgen.o \
         $(nfpsrc)hansr5.o \
         $(nfpsrc)hansr4.o \
         $(nfpsrc)lattc.o \
         $(nfpsrc)ll.o \
         $(nfpsrc)dpcopy.o \
         $(nfpsrc)dpadd.o \
         $(nfpsrc)syscalls.o \
         $(nfpsrc)qdist.o \
         $(nfpsrc)dlmtor.o \
         $(nfpsrc)dpzero.o \
         $(nfpsrc)ropyln.o \
         $(nfpsrc)ropcsm.o \
         $(nfpsrc)dsisl.o \
         $(nfpsrc)dsifa.o \
         $(nfpsrc)diagcv2.o \
         $(gwsrc)scg.o 

 SIGMCONV = \
 $(gwsrc)switches.o \
 $(gwsrc)keyvalue.o \
 $(gwsrc)iopen.o \
 $(main)hsigmconv.m.o 

###########################################

# bndconn.o:	$(BNDCONN)  
#
############### dependency for use ##################



 hsigmconv:	$(SIGMCONV)  $(MPI)  
	 $(LK) $(LKFLAGS1) $(SIGMCONV) $(MPI)  $(LKFLAGS2) -o $@


 gwinit:	$(GWINIT)  $(MPI)  
	 $(LK) $(LKFLAGS1) $(GWINIT) $(MPI)  $(LKFLAGS2) -o $@


 qpwf:		$(maxloc)qpwf.o $(GW0 $(MPI) )
	$(LK) $(LKFLAGS1) $(maxloc)qpwf.o $(GW0) $(MPI)   $(LKFLAGS2) -o $@

 qg4gw:		$(QG)  $(MPI)  
	 $(LK) $(LKFLAGS1) $(QG) $(MPI)   $(LKFLAGS2) -o $@

 rdata4gw_v2:	$(RDAT_v2)  $(NFPL $(MPI) )
	 $(LK) $(LKFLAGS1) $(RDAT_v2) $(NFPL) $(MPI)  $(LKFLAGS2) -o $@

 hbasfp0:	$(BAS)  $(MPI) 
	 $(LK) $(LKFLAGS1) $(BAS) $(MPI)  $(LKFLAGS2) -o $@

 hvccfp0:	$(MPI) $(VCC)  $(NFPL) $(DERFC) $(MPI) 
	 $(LK) $(LKFLAGS1) $(VCC) $(NFPL) $(DERFC) $(MPI) $(LKFLAGS2) -o $@

 hx0fp0:	$(MPI) $(X0) $(GW0) $(MPI) 
	 $(LK) $(LKFLAGS1) $(X0)     $(GW0) $(MPI) $(LKFLAGS2) -o $@

# for maxloc
 hx0fp0_mlw:	$(X0mlw) $(GW0) $(MPI) 
	$(LK) $(LKFLAGS1) $(X0mlw)     $(GW0)  $(MPI)  $(LKFLAGS2) -o $@

 h_uumatrix:	$(UU) $(GW0)   $(MPI) 
	$(LK) $(LKFLAGS1) $(UU)     $(GW0)  $(MPI)  $(LKFLAGS2) -o $@

 huumat:	$(UU2) $(GW0) $(MPI) 
	$(LK) $(LKFLAGS1) $(UU2)     $(GW0)  $(MPI)  $(LKFLAGS2) -o $@

 hphig:	$(PHIG) $(GW0) $(MPI) 
	$(LK) $(LKFLAGS1) $(PHIG)     $(GW0)  $(MPI)  $(LKFLAGS2) $(LIBSLA) -o $@

 hpsig: $(PSIG) $(GW0) $(MPI) 
	$(LK) $(LKFLAGS1) $(PSIG)     $(GW0)  $(MPI)  $(LKFLAGS2) -o $@

 hx0fp0_sc:	$(MPI) $(X0_SC) $(GW0)   $(MPI) 
	 $(LK) $(LKFLAGS1) $(X0_SC)     $(GW0) $(MPI)  $(LKFLAGS2) -o $@

 hwmat:		$(WMAT) $(GW0) $(MPI) 
	$(LK) $(LKFLAGS1) $(WMAT)    $(GW0)  $(MPI)  $(LKFLAGS2) -o $@

 hmaxloc:	$(MLOC)  $(NFPLtot) $(GW0) $(MPI) 
	$(LK) $(LKFLAGS1) $(MLOC) $(NFPLtot) $(GW0)  $(MPI)  $(LKFLAGS2) -o $@

 hmaxloc1D:	$(MLOC1D)  $(NFPLtot) $(GW0) $(MPI) 
	$(LK) $(LKFLAGS1) $(MLOC1D) $(NFPLtot) $(GW0)   $(MPI) $(LKFLAGS2) -o $@

 hsfp0:		$(MPI) $(SXC) $(GW0)   $(MPI) 
	 $(LK) $(LKFLAGS1) $(SXC)    $(GW0) $(MPI)  $(LKFLAGS2) -o $@

 hsfp0_sc:	$(MPI) $(SXC_SC) $(GW0)   $(MPI) 
	 $(LK) $(LKFLAGS1) $(SXC_SC)    $(GW0)  $(MPI) $(LKFLAGS2) -o $@

 hnocc_mlw:	$(hnocc_mlw) $(GW0) $(MPI) 
	$(LK) $(LKFLAGS1) $(hnocc_mlw) $(GW0)  $(MPI)  $(LKFLAGS2) -o $@

 heftet:	$(heftet) $(GW0)  $(MPI) $(MPI) 
	 $(LK) $(LKFLAGS1) $(heftet) $(GW0)  $(MPI) $(LKFLAGS2) -o $@

 hef:		$(hef) $(GW0)   $(MPI)
	 $(LK) $(LKFLAGS1) $(hef)    $(GW0)  $(MPI) $(LKFLAGS2) -o $@

 hchknw:	$(CHK) $(GW0)  $(MPI)
	 $(LK) $(LKFLAGS1) $(CHK)    $(GW0)   $(MPI) $(LKFLAGS2) -o $@

 hqpe:		$(QPE) $(MPI)
	 $(LK) $(LKFLAGS1) $(QPE) $(MPI) $(LKFLAGS2) -o $@

 hqpe_sc:		$(QPE_SC)  $(MPI)
	 $(LK) $(LKFLAGS1) $(QPE_SC) $(MPI) $(LKFLAGS2) -o $@

 hmergewv:	$(MERGE) $(MPI)  
	 $(LK) $(LKFLAGS1) $(MERGE) $(MPI)  $(LKFLAGS2) -o $@

 hparainfo:	$(PARAINFO) $(GW0) $(MPI)  
	 $(LK) $(LKFLAGS1) $(PARAINFO) $(GW0)  $(MPI) $(LKFLAGS2) -o $@

 hbndout:	$(BNDOUT) $(MPI)  
	 $(LK) $(LKFLAGS1) $(BNDOUT) $(MPI)  $(LKFLAGS2) -o $@

 convgwin:	$(convg) $(MPI)  
	 $(LK) $(LKFLAGS1) $(convg) $(MPI)  $(LKFLAGS2) -o $@

 kino_input_test:	$(kino_input_test) $(GW0) $(MPI) 
	 $(LK) $(LKFLAGS1) $(kino_input_test) $(GW0)  $(MPI) $(LKFLAGS2) -o $@

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


check:
	(cd ../TESTinstallGW;./testgw.py --enforce --all) 

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

# #$(gwsrc)x0kf_v4h$(x0kf_g): $(gwsrc)x0kf_v4h.F 
# #	gawk -f $(addtime) -vSTART=100 $(gwsrc)x0kf_v4h.F | gawk -f $(septhen) | gawk -f $(alloclist) > $(gwsrc)time_x0kf_v4h.F
# #	$(FC) $(FFLAGS) $(gwsrc)time_x0kf_v4h.F -c -o $*.o

# $(gwsrc)ppbafp.fal$(para_g): $(gwsrc)ppbafp.fal.F
# 	gawk -f $(addtime) -vSTART=200 $(gwsrc)ppbafp.fal.F > $(gwsrc)time_ppbafp.fal.F
# 	$(FC) $(FFLAGS) $(gwsrc)time_ppbafp.fal.F -c -o $*.o
# ###################################################################################




# DO NOT DELETE
