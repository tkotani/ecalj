yourbin=$(HOME)/bin

#FC=mpiifort
#LK=mpiifort
#LKFLAGS2=-mkl
#FFLAGS = -132 -ffixed-line-length-132 -c -frange-check -I${top}/exec

FC=mpif90
LK=mpif90
LKFLAGS2= /usr/lib/libfftw3.so.3 /usr/lib/liblapack.so.3gf /usr/lib/libblas.so.3gf 
FFLAGS = -ffixed-line-length-132 -c -frange-check -I${top}/exec


### I think you don't need to modify following lines. Not machine-dependent ###########
x0kf_g=.o
para_g=.o

# ---- Machine-specific compiler flags ---
#include make.inc

#-------------------------------------------------------
# src directories
top=../fpgw
gwsrc   = ${top}/gwsrc/
main    = ${top}/main/
nfpsrc  = ${top}/nfpsrc/
slatsmlib  =${top}/../lm7K/slatsm/
tote = ../tote/
# tag directory
tags   = ./

#progs  = hbasfp0 hvccfp0 hx0fp0 hsfp0 hef hqpe hchknw qg4gw gwinit heftet hmergewv hparainfo hbndout rdata4gw_v2 convgwin hx0fp0_sc hsfp0_sc hqpe_sc kino_input_test hecor eout eout2 h_uumatrix hsigmconv
progs= huumat hpsig hmaxloc huumat_MPI hpsig_MPI wanplot #hwmat hwmat_MPI
progs_RSMPI = hvccfp0_RSMPI hx0fp0_RSMPI hsfp0_RSMPI
progs_RSBAND = qg4gw_RSBAND
# lmf_exec
#progs  = hbasfp0 hvccfp0 hx0fp0 hsfp0 hef hqpe hchknw qg4gw gwinit heftet hmergewv hparainfo hbndout rdata4gw_v2  hx0fp0_fal hx0fp1 
#progs2 = $(progs) $(tags)TAGS checkmod 
progs2 = $(progs) 
progs_rs = $(progs_RSMPI) $(progs_RSBAND)

#script = cleargw* dqpu dtote eps* ex* gw* hqpemetal* inf* lmgw* plotg save* tote_lmfh2 xqp

#### You can choose these options. all is default.

all :$(progs2) 
all_rs	:$(progs_rs)

#clean:  
#	 rm -f  $(progs)
clean:  
	 rm -f  $(progs2) *.o 

install:  
	 cp  $(progs) genMLWF $(yourbin)

#install2:  
#	 cp  $(script) $(HOME)/bin

cleanall:  
	 rm -f  $(progs2) $(maxloc)*.o

#	 rm -f  $(progs2) $(main)*.o $(gwsrc)*.o  $(maxloc)*.o  $(gwt)*.o 

# This is necesaly to compile *.f in right order.
# When you recompile and link, just repeat 'make' (not necessary to repeat 'make init').
# When checkmodule recompile source, you have to repeat 'make'.
init:
	exec ./checkmodule

checkmod:
	exec ./checkmodule

wanplot.o: cubeformat.o xsfformat2.o  

####from tote #################################################################################                                     
#LIBLOC  = $(ECAL)/fftw/libfftw.a $(LIBMATH)
##-L/usr/local/ATLAS/lib/Linux_P4SSE2 -llapack -lcblas -lf77blas -latlas                                                             
#LIBSLA  = $(ECAL)/slatsm/slatsm.a
#LIBFP   = $(ECAL)/lm-6.14y/fp/subs.a
#LIBSUBS = $(ECAL)/lm-6.14y/subs/subs.a
#LIBES  = $(LIBSLA) $(LIBLOC)
#lmsrc   = ../../lm-6.14y/
######################################################################################          

OBJS2=m_DATA4GW.o m_QG.o m_FFT3D.o mymath.o myYlm.o wfn2dx_2.o m_MLWF.o cubeformat.o xsfformat2.o  expand_mesh.o  wfn2dx_abc.o wfnrho_abc.o write_mesh2d.o 
 

ECOR = \
$(tote)hecor.o \
$(tote)rpaq.o 
#eispack.o

NFPLtot  = $(nfpsrc)diagcv2.o 

GW0tot = \
$(gwsrc)m_hamindex.o\
$(gwsrc)MPI_fpgw.o \
$(gwsrc)alloclist.o \
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

hecor: $(ECOR) $(NFPLtot) $(GW0)  
	$(LK) $(LKFLAGS1) $(ECOR) $(GW0) $(NFPLtot) $(LKFLAGS2) -o $@

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
          $(nfpsrc)i1mach.o \
          $(nfpsrc)d1mach.o 

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
 $(gwsrc)iopenxx.o \
 $(gwsrc)iprint.o \
 $(gwsrc)keyvalue.o \
 $(gwsrc)switches.o 

 QG =  \
 $(gwsrc)conv2gwinput.o \
 $(main)qg4gw.m.o \
 $(gwsrc)getbzdata1.o \
 $(gwsrc)mkqg.o \
 $(gwsrc)q0irre.o \
 $(gwsrc)getgv2.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)tetwt5.o \
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
 $(gwsrc)polinta.o \
 $(gwsrc)llnew.o


 VCC= \
 $(main)hvccfp0.m.o \
 $(gwsrc)mkjp.o \
 $(gwsrc)gintxx.o \
 $(gwsrc)extension.o \
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
 $(gwsrc)iqindx.o 

# RS: VCC_RSMPI
 VCC_RSMPI= \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(main)hvccfp0_RSMPI.m.o \
 $(gwsrc)mkjp.o \
 $(gwsrc)gintxx.o \
 $(gwsrc)extension.o \
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
 $(gwsrc)iqindx.o 

 SXC_SC = \
 $(main)hsfp0.sc.m.o \
 $(gwsrc)wse.o \
 $(gwsrc)sxcf_fal2$(sxcf_g) \
 $(gwsrc)sxcf_fal2.sc$(sxcf_g) \
 $(gwsrc)bzints2.o \
 $(gwsrc)wintzsg.o

 SXC = \
 $(main)hsfp0.m.o \
 $(gwsrc)wse.o \
 $(gwsrc)wintzsg.o \
 $(gwsrc)sxcf_fal2$(sxcf_g) \
 $(gwsrc)sxcf_fal2.sc$(sxcf_g) \
 $(gwsrc)bzints2.o \
 $(gwsrc)genallcf_dump.o

#RS:
 SXC_RSMPI = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_qkgroup_mod.o \
 $(maxloc)RSMPI_rotkindex_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(main)hsfp0_RSMPI.m.o \
 $(gwsrc)efsimplef_RSMPI.o \
 $(gwsrc)sxcf_fal2_RSMPI.o \
 $(gwsrc)wse.o \
 $(gwsrc)wintzsg.o \
 $(gwsrc)sxcf_fal2$(sxcf_g) \
 $(gwsrc)sxcf_fal2.sc$(sxcf_g) \
 $(gwsrc)bzints2.o \
 $(gwsrc)genallcf_dump.o

 heftet = \
 $(main)heftet.m.o \
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
 $(gwsrc)tetwt5.o \
 $(gwsrc)x0kf_v2h$(x0kf_g) \
 $(gwsrc)x0kf_v3h$(x0kf_g)

 X0 = \
 $(main)hx0fp0.m.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5.o \
 $(gwsrc)x0kf_v2h$(x0kf_g) \
 $(gwsrc)x0kf_v3h$(x0kf_g) \
 $(nfpsrc)diagcv2.o \
 $(tote)rpaq.o \
 $(gwsrc)cinvrx.o

# RS: 
 X0_RSMPI = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_qkgroup_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(gwsrc)x0kf_v3h_RSMPI.o \
 $(main)hx0fp0_RSMPI.m.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5.o \
 $(gwsrc)x0kf_v2h$(x0kf_g) \
 $(gwsrc)x0kf_v3h$(x0kf_g) \
 $(nfpsrc)diagcv2.o \
 $(tote)rpaq.o \
 $(gwsrc)cinvrx.o

 UU = \
 $(main)h_uumatrix.m.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5.o \
 $(gwsrc)gintxx.o \
 $(gwsrc)pplmat.o \
 $(gwsrc)getgv2.o \
 $(gwsrc)x0kf_v2h$(x0kf_g) \
 $(gwsrc)x0kf_v3h$(x0kf_g) \
 $(gwsrc)rs.o \
 $(gwsrc)mopen.o \
 $(nfpsrc)u_lat_0.o \
 $(nfpsrc)wronkj.o \
 $(nfpsrc)mklegw.o \
 $(nfpsrc)rxx.o \
 $(nfpsrc)bessl.o \
 $(nfpsrc)cross.o 

 GW0 = \
 $(gwsrc)m_hamindex.o\
 $(gwsrc)MPI_fpgw.o \
 $(gwsrc)alloclist.o \
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
 $(gwsrc)iqindx.o \
 $(NFPLtot)

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

 RSEC = \
 $(gwsrc)switches.o \
 $(gwsrc)keyvalue.o \
 $(main)read_sec$(hqpe_g) \
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
 $(main)hsigmconv.m.o 

###########################################

# bndconn.o:	$(BNDCONN)  
#
############### dependency for use ##################



 hsigmconv:	$(SIGMCONV)  
	 $(LK) $(LKFLAGS1) $(SIGMCONV) $(LKFLAGS2) -o $@


 gwinit:	$(GWINIT)  
	 $(LK) $(LKFLAGS1) $(GWINIT) $(LKFLAGS2) -o $@


 qg4gw:		$(QG)  
	 $(LK) $(LKFLAGS1) $(QG)  $(LKFLAGS2) -o $@

 rdata4gw_v2:	$(RDAT_v2)  $(NFPL)
	 $(LK) $(LKFLAGS1) $(RDAT_v2) $(NFPL) $(LKFLAGS2) -o $@

 hbasfp0:	$(BAS) 
	 $(LK) $(LKFLAGS1) $(BAS) $(LKFLAGS2) -o $@

 hvccfp0:	$(VCC)  $(NFPL) $(DERFC)
	 $(LK) $(LKFLAGS1) $(VCC) $(NFPL) $(DERFC) $(LKFLAGS2) -o $@

 hvccfp0_RSMPI:	$(VCC_RSMPI)  $(NFPL) $(DERFC)
	 $(LK) $(LKFLAGS1) $(VCC_RSMPI) $(NFPL) $(DERFC) $(LKFLAGS2) -o $@

 hx0fp0:	$(X0) $(GW0)  
	 $(LK) $(LKFLAGS1) $(X0)     $(GW0)  $(LKFLAGS2) -o $@

 hx0fp0_RSMPI:	$(X0_RSMPI) $(GW0)  
	 $(LK) $(LKFLAGS1) $(X0_RSMPI)     $(GW0)  $(LKFLAGS2) -o $@

 h_uumatrix:	$(UU) $(GW0)  
	$(LK) $(LKFLAGS1) $(UU)     $(GW0)  $(LKFLAGS2) -o $@

 hx0fp0_sc:	$(X0_SC) $(GW0)  
	 $(LK) $(LKFLAGS1) $(X0_SC)     $(GW0)  $(LKFLAGS2) -o $@

 hsfp0:		$(SXC) $(GW0)  
	 $(LK) $(LKFLAGS1) $(SXC)    $(GW0)  $(LKFLAGS2) -o $@

 hsfp0_RSMPI:		$(SXC_RSMPI) $(GW0)  
	 $(LK) $(LKFLAGS1) $(SXC_RSMPI)   $(GW0)  $(LKFLAGS2) -o $@

 hsfp0_sc:		$(SXC_SC) $(GW0)  
	 $(LK) $(LKFLAGS1) $(SXC_SC)    $(GW0)  $(LKFLAGS2) -o $@

 heftet:	$(heftet) $(GW0)  
	 $(LK) $(LKFLAGS1) $(heftet) $(GW0)  $(LKFLAGS2) -o $@

 hef:		$(hef) $(GW0)  
	 $(LK) $(LKFLAGS1) $(hef)    $(GW0)  $(LKFLAGS2) -o $@

 hchknw:	$(CHK) $(GW0)  
	 $(LK) $(LKFLAGS1) $(CHK)    $(GW0)  $(LKFLAGS2) -o $@

 hqpe:		$(QPE) 
	 $(LK) $(LKFLAGS1) $(QPE) $(LKFLAGS2) -o $@

 hqpe_sc:		$(QPE_SC) 
	 $(LK) $(LKFLAGS1) $(QPE_SC) $(LKFLAGS2) -o $@

 read_sec:		$(RSEC) 
	 $(LK) $(LKFLAGS1) $(RSEC) $(LKFLAGS2) -o $@

 hmergewv:	$(MERGE) 
	 $(LK) $(LKFLAGS1) $(MERGE) $(LKFLAGS2) -o $@


 hparainfo:	$(PARAINFO) $(GW0) 
	 $(LK) $(LKFLAGS1) $(PARAINFO) $(GW0) $(LKFLAGS2) -o $@

 hbndout:	$(BNDOUT) 
	 $(LK) $(LKFLAGS1) $(BNDOUT) $(LKFLAGS2) -o $@

 convgwin:	$(convg) 
	 $(LK) $(LKFLAGS1) $(convg) $(LKFLAGS2) -o $@

 kino_input_test:	$(kino_input_test) $(GW0)
	 $(LK) $(LKFLAGS1) $(kino_input_test) $(GW0) $(LKFLAGS2) -o $@

################################ test
#
# test_genallcf:	$(test_genallcf) 
#	 $(LK) $(LKFLAGS1) $(test_genallcf) $(LKFLAGS2) -o $@


 $(tags)TAGS: $(progs)
	cd $(tags);etags ./*/*.f ./*/*.F


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


##### Maxloc Wannier fn. and Finite-temperature GW
### Takashi MIYAKE
maxloc = ./
gwt = ../Miyake/gwt/
misc = ../Miyake/misc/

#progs_tm = hwmat hwmatc hmaxloc huumat qpwf hpsig hnocc_mlw hx0fp0_mlw hx0fp0_mlw_ent hphig hx0tau hs0tau hpsipsi hibz2fbz hhmnk hwtmat hnocc_mlwt hwinfo hpade hibz2fbz0 read_gw0 # hdyson # hmaxloc1D
progs_TM= hmaxloc huumat qpwf hpsig hnocc_mlw hx0fp0_mlw hwmat hwmatc hx0tau hs0tau hpsipsi hibz2fbz hhmnk hwtmat hnocc_mlwt hwinfo hwfch hibz2fbz0 read_gw0 qg4gw_TMBAND hqpe_mlw hibz2fbz0wan # hdyson hphig hpade haw
progs_TMMPI= hx0tau_MPI hs0tau_MPI hwtmat_MPI hwmat_MPI hpsig_MPI huumat_MPI hx0fp0_mlw_MPI hwmatc_MPI hsfp0_mlw_MPI hband_ent_MPI hqpband_ent_MPI hupdt_ent_MPI

progs_tm = $(progs_TM) $(progs_TMMPI)
all_tm:	$(progs_tm)

 hnocc_mlw = \
 $(maxloc)hnocc_mlw.o \
 $(gwsrc)bzints2.o

 hnocc_mlwt = \
 $(gwt)hnocc_mlwt.o \
 $(gwsrc)bzints2.o \
 $(gwt)ferdi.o \
 $(gwt)efsimpleft.o

 AW = \
 $(gwsrc)switches.o \
 $(gwsrc)keyvalue.o \
 $(misc)haw$(hqpe_g) \
 $(gwsrc)icompvv2.o \
 $(gwsrc)iopenxx.o \
 $(gwsrc)iopen.o \
 $(gwsrc)rw.o \
 $(gwsrc)rydberg.o

 BAND_ENT = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)hband_ent_MPI.m.o \
 $(maxloc)maxloc0.o \
 $(maxloc)maxloc1.o\
 $(maxloc)ent_MPI.o

 BAND_ENT2 = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)hband_ent_MPI2.m.o \
 $(maxloc)maxloc0.o \
 $(maxloc)maxloc1.o\
 $(maxloc)ent_MPI.o

 DOS_ENT = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)hdos_ent_MPI.m.o \
 $(maxloc)maxloc0.o \
 $(maxloc)maxloc1.o\
 $(maxloc)ent_MPI.o

 DYSON = \
 $(gwt)hdyson.o \
 $(gwt)dysonmu.o \
 $(gwt)geng0.o \
 $(gwt)gentau.o \
 $(gwt)filong.o \
 $(gwsrc)iopen.o \
 $(gwsrc)matm.o

 GW0tm = \
 $(gwsrc)m_hamindex.o\
 $(gwsrc)MPI_fpgw.o \
 $(gwsrc)alloclist.o \
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
 $(gwsrc)iqindx.o \
 $(gwsrc)rangedq.o \
 $(NFPLtot)

 GWTCHK = \
 $(gwt)hgwtchk.o \
 $(gwt)filonx.o \
 $(gwt)filong.o

 HDEBUG = \
 $(gwt)hdebug.o \
 $(gwt)filonx.o

 IBZ2FBZ = \
 $(gwt)hibz2fbz.o \
 $(gwt)filonx.o

 IBZ2FBZ0 = \
 $(maxloc)hibz2fbz0.o \
 $(gwt)filonx.o

 IBZ2FBZ0W = \
 $(maxloc)hibz2fbz0wan.o \
 $(gwt)filonx.o

 IBZ2FBZ0_TEST = \
 $(maxloc)hibz2fbz0_test.o \
 $(gwt)filonx.o

 IBZ2FBZ_TEST = \
 $(gwt)hibz2fbz_test.o \
 $(gwt)filonx.o

 MLOC = \
 $(nfpsrc)rxx.o \
 $(gwsrc)mopen.o \
 $(maxloc)hmaxloc.o \
 $(maxloc)maxloc0.o \
 $(maxloc)maxloc1.o \
 $(maxloc)maxloc2.o \
 $(maxloc)maxloc3.o \
 $(maxloc)maxloc4.o \
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

 PADE = \
 $(gwt)hpade.o \
 $(gwt)pade.o \
 $(gwt)filong.o

 PHIG = \
 $(maxloc)hphig.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5.o \
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
 $(gwsrc)tetwt5.o \
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

 PSIG_MPI = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_qkgroup_mod.o \
 $(maxloc)RSMPI_rotkindex_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)hpsig_MPI.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5.o \
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

 PSIG_TEST = \
 $(maxloc)hpsig_test.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5.o \
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

 PSIPSI = \
 $(gwt)hpsipsi.o \
 $(gwsrc)gintxx.o \
 $(gwsrc)pplmat.o \
 $(gwsrc)getgv2.o \
 $(gwsrc)rs.o \
 $(gwsrc)mopen.o \
 $(nfpsrc)u_lat_0.o \
 $(nfpsrc)wronkj.o \
 $(nfpsrc)mklegw.o \
 $(nfpsrc)rxx.o \
 $(nfpsrc)bessl.o \
 $(nfpsrc)cross.o 

 PSIPSI_TEST = \
 $(gwt)hpsipsi_test.o \
 $(gwsrc)gintxx.o \
 $(gwsrc)pplmat.o \
 $(gwsrc)getgv2.o \
 $(gwsrc)rs.o \
 $(gwsrc)mopen.o \
 $(nfpsrc)u_lat_0.o \
 $(nfpsrc)wronkj.o \
 $(nfpsrc)mklegw.o \
 $(nfpsrc)rxx.o \
 $(nfpsrc)bessl.o \
 $(nfpsrc)cross.o 

 QG_TM =  \
 $(gwsrc)conv2gwinput.o \
 $(maxloc)qg4gw_TMBAND.m.o \
 $(gwsrc)getbzdata1.o \
 $(maxloc)mkqg_TMBAND.o \
 $(gwsrc)q0irre.o \
 $(gwsrc)getgv2.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)tetwt5.o \
 $(GW0)

 QPBAND_ENT = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)hqpband_ent_MPI.m.o \
 $(maxloc)maxloc0.o \
 $(maxloc)maxloc1.o\
 $(maxloc)ent_MPI.o

 QPBAND_ENT_TEST = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)hqpband_ent_MPI_test.m.o \
 $(maxloc)maxloc0.o \
 $(maxloc)maxloc1.o\
 $(maxloc)ent_MPI.o

 QPE_MLW = \
 $(maxloc)hqpe_mlw.m$(hqpe_g) \
 $(maxloc)qpe1_mlw.o \
 $(gwsrc)icompvv2.o \
 $(gwsrc)iopenxx.o \
 $(gwsrc)genallcf_dump.o

 SXC_MPI = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_qkgroup_mod.o \
 $(maxloc)RSMPI_rotkindex_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)hsfp0_MPI.m.o \
 $(gwsrc)efsimplef_RSMPI.o \
 $(maxloc)sxcf_fal2_MPI.o \
 $(gwsrc)wse.o \
 $(gwsrc)wintzsg.o \
 $(gwsrc)sxcf_fal2$(sxcf_g) \
 $(gwsrc)sxcf_fal2.sc$(sxcf_g) \
 $(gwsrc)bzints2.o \
 $(gwsrc)genallcf_dump.o

 SXCTAU = \
 $(gwt)hs0tau.m.o \
 $(gwsrc)wse.o \
 $(gwsrc)wintzsg.o \
 $(gwsrc)sxcf_fal2$(sxcf_g) \
 $(gwsrc)sxcf_fal2.sc$(sxcf_g) \
 $(gwt)sxcf_tau$(sxcf_g) \
 $(gwt)ferdi.o \
 $(gwt)geng0.o \
 $(gwsrc)bzints2.o \
 $(gwsrc)genallcf_dump.o \
 $(gwt)efsimpleft.o

 SXCTAU_MPI = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_qkgroup_mod.o \
 $(maxloc)RSMPI_rotkindex_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(gwt)hs0tau_MPI.m.o \
 $(gwsrc)wse.o \
 $(gwsrc)wintzsg.o \
 $(gwsrc)sxcf_fal2$(sxcf_g) \
 $(gwsrc)sxcf_fal2.sc$(sxcf_g) \
 $(gwt)sxcf_tau$(sxcf_g) \
 $(gwt)sxcf_tau_MPI$(sxcf_g) \
 $(gwt)ferdi.o \
 $(gwt)geng0.o \
 $(gwsrc)bzints2.o \
 $(gwsrc)genallcf_dump.o \
 $(gwsrc)efsimplef_RSMPI.o \
 $(gwt)efsimpleft_MPI.o

 UPDT_ENT = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)hupdt_ent_MPI.m.o \
 $(maxloc)maxloc0.o \
 $(maxloc)maxloc1.o\
 $(maxloc)ent_MPI.o

 UU2 = \
 $(maxloc)huumat.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5.o \
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

 UU2_MPI = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)huumat_MPI.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5.o \
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
 $(nfpsrc)cross.o \
 $(maxloc)ent_MPI.o

 UU2_TEST = \
 $(maxloc)huumat_test.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5.o \
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

 WFCH = \
 $(maxloc)hwfch.o \
 $(gwsrc)gintxx.o \
 $(gwsrc)rs.o \
 $(gwsrc)mopen.o \
 $(nfpsrc)u_lat_0.o \
 $(nfpsrc)wronkj.o \
 $(nfpsrc)mklegw.o \
 $(nfpsrc)rxx.o \
 $(nfpsrc)bessl.o \
 $(nfpsrc)cross.o

 WMAT = \
 $(maxloc)hwmat.o \
 $(gwsrc)wse.o \
 $(maxloc)wmat.o \
 $(maxloc)maxloc0.o \
 $(gwsrc)genallcf_dump.o

 WMAT_MPI = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_qkgroup_mod.o \
 $(maxloc)RSMPI_rotkindex_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)hwmat_MPI.o \
 $(gwsrc)wse.o \
 $(maxloc)wmat_MPI.o \
 $(maxloc)maxloc0.o \
 $(gwsrc)genallcf_dump.o \
 $(gwsrc)efsimplef_RSMPI.o

 WMAT2 = \
 $(maxloc)hwmat2.o \
 $(gwsrc)wse.o \
 $(maxloc)wmat2.o \
 $(maxloc)maxloc0.o \
 $(gwsrc)genallcf_dump.o

 WMATC = \
 $(maxloc)hwmatc.o \
 $(gwsrc)wse.o \
 $(maxloc)wmat.o \
 $(maxloc)maxloc0.o \
 $(gwsrc)genallcf_dump.o

 WMATC_MPI = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_qkgroup_mod.o \
 $(maxloc)RSMPI_rotkindex_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)hwmat_MPI.o \
 $(gwsrc)wse.o \
 $(maxloc)wmat_MPI.o \
 $(maxloc)maxloc0.o \
 $(gwsrc)genallcf_dump.o \
 $(gwsrc)efsimplef_RSMPI.o

 WTMAT = \
 $(gwt)hwtmat.o \
 $(gwt)filong.o \
 $(gwsrc)wse.o \
 $(gwt)wtmat.o \
 $(gwsrc)genallcf_dump.o

 WTMAT_MPI = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_qkgroup_mod.o \
 $(maxloc)RSMPI_rotkindex_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(gwt)hwtmat_MPI.o \
 $(gwt)filong.o \
 $(gwsrc)wse.o \
 $(gwt)wtmat_MPI.o \
 $(gwsrc)genallcf_dump.o \
 $(gwsrc)efsimplef_RSMPI.o

 X0mlw = \
 $(maxloc)hx0fp0.m.o \
 $(maxloc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5.o \
 $(gwsrc)x0kf_v2h$(x0kf_g)\
 $(maxloc)x0kf_v3h$(x0kf_g) \
 $(nfpsrc)diagcv2.o \
 $(tote)rpaq.o \
 $(gwsrc)cinvrx.o

#
 X0mlw_MPI = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_qkgroup_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)hx0fp0_MPI.m.o \
 $(maxloc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5.o \
 $(gwsrc)x0kf_v2h$(x0kf_g) \
 $(gwsrc)x0kf_v3h$(x0kf_g) \
 $(maxloc)x0kf_v3h_MPI$(x0kf_g) \
 $(nfpsrc)diagcv2.o \
 $(tote)rpaq.o \
 $(gwsrc)cinvrx.o \
 $(maxloc)ent_MPI.o

 X0mlw_cfao = \
 $(maxloc)hx0fp0_cfao.m.o \
 $(maxloc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5.o \
 $(gwsrc)x0kf_v2h$(x0kf_g)\
 $(maxloc)x0kf_v3h_cfao$(x0kf_g) \
 $(nfpsrc)diagcv2.o \
 $(tote)rpaq.o \
 $(gwsrc)cinvrx.o

 X0tau = \
 $(gwt)hx0tau.m.o \
 $(maxloc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwt)x0tau$(x0kf_g) \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5.o \
 $(gwsrc)x0kf_v2h$(x0kf_g)\
 $(maxloc)x0kf_v3h$(x0kf_g) \
 $(nfpsrc)diagcv2.o \
 $(tote)rpaq.o \
 $(gwsrc)cinvrx.o \
 $(gwt)ferdi.o \
 $(gwt)filong.o \
 $(gwt)gentau.o \
 $(gwt)wtau.o \
 $(gwt)efsimpleft.o

 X0tau_MPI = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_qkgroup_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(gwt)hx0tau_MPI.m.o \
 $(maxloc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwt)x0tau$(x0kf_g) \
 $(gwt)x0tau_MPI$(x0kf_g) \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5.o \
 $(gwsrc)x0kf_v2h$(x0kf_g)\
 $(maxloc)x0kf_v3h$(x0kf_g) \
 $(maxloc)x0kf_v3h_MPI$(x0kf_g) \
 $(nfpsrc)diagcv2.o \
 $(tote)rpaq.o \
 $(gwsrc)cinvrx.o \
 $(gwt)ferdi.o \
 $(gwt)filong.o \
 $(gwt)gentau.o \
 $(gwt)wtau.o \
 $(gwt)efsimpleft_MPI.o
# $(gwt)efsimpleft.o

 X0tau_test = \
 $(gwt)hx0tau_test.m.o \
 $(maxloc)wcf.o \
 $(gwsrc)tetwt4.o \
 $(gwt)x0tau$(x0kf_g) \
 $(gwsrc)x0kf$(x0kf_g) \
 $(gwsrc)tetwt5.o \
 $(gwsrc)x0kf_v2h$(x0kf_g)\
 $(maxloc)x0kf_v3h$(x0kf_g) \
 $(nfpsrc)diagcv2.o \
 $(tote)rpaq.o \
 $(gwsrc)cinvrx.o \
 $(gwt)ferdi.o \
 $(gwt)filong.o \
 $(gwt)gentau.o \
 $(gwt)wtau.o \
 $(gwt)efsimpleft.o

COMM=$(nfpsrc)rxx.o \
 $(gwsrc)mopen.o \
 $(gwsrc)rangedq.o

COMMX=$(nfpsrc)rxx.o \
 $(gwsrc)mopen.o

 e2h:	$(maxloc)e2h.o  $(NFPLtot)
	$(LK) $(LKFLAGS1) $(maxloc)e2h.o  $(NFPLtot)  $(LKFLAGS2) -o $@

 haw:		$(AW) 
	 $(LK) $(LKFLAGS1) $(AW) $(LKFLAGS2) -o $@

 hband_ent_MPI:		$(BAND_ENT) $(GW0tm)
	 $(LK) $(LKFLAGS1) $(BAND_ENT) $(GW0tm) $(LKFLAGS2) -o $@

 hband_ent_MPI2:		$(BAND_ENT2) $(GW0tm)
	 $(LK) $(LKFLAGS1) $(BAND_ENT2) $(GW0tm) $(LKFLAGS2) -o $@

 hdos_ent_MPI:		$(DOS_ENT) $(GW0tm)
	 $(LK) $(LKFLAGS1) $(DOS_ENT) $(GW0tm) $(LKFLAGS2) -o $@

 hdebug:	$(HDEBUG) $(GW0tm)
	$(LK) $(LKFLAGS1) $(HDEBUG)    $(GW0tm)  $(LKFLAGS2) -o $@

 hdyson:	$(DYSON)
	$(LK) $(LKFLAGS1) $(DYSON)    $(LKFLAGS2) -o $@

 hgwtchk:	$(GWTCHK) $(GW0tm)
	$(LK) $(LKFLAGS1) $(GWTCHK)    $(GW0tm)  $(LKFLAGS2) -o $@

 hhmnk:		$(gwt)hhmnk.o $(GW0tm)
	$(LK) $(LKFLAGS1) $(gwt)hhmnk.o     $(GW0tm)  $(LKFLAGS2) -o $@

 hibz2fbz:	$(IBZ2FBZ) $(GW0tm)
	$(LK) $(LKFLAGS1) $(IBZ2FBZ)    $(GW0tm)  $(LKFLAGS2) -o $@

 hibz2fbz0:	$(IBZ2FBZ0) $(GW0tm)
	$(LK) $(LKFLAGS1) $(IBZ2FBZ0)    $(GW0tm)  $(LKFLAGS2) -o $@

 hibz2fbz0wan:	$(IBZ2FBZ0W) $(GW0tm)
	$(LK) $(LKFLAGS1) $(IBZ2FBZ0W)    $(GW0tm)  $(LKFLAGS2) -o $@

 hibz2fbz0_test:	$(IBZ2FBZ0_TEST) $(GW0tm)
	$(LK) $(LKFLAGS1) $(IBZ2FBZ0_TEST)    $(GW0tm)  $(LKFLAGS2) -o $@

 hibz2fbz_test:	$(IBZ2FBZ_TEST) $(GW0tm)
	$(LK) $(LKFLAGS1) $(IBZ2FBZ_TEST)    $(GW0tm)  $(LKFLAGS2) -o $@

 hmaxloc:	$(MLOC)  $(GW0tm)
	$(LK) $(LKFLAGS1) $(MLOC) $(GW0tm)  $(LKFLAGS2) -o $@

 hmaxloc1D:	$(MLOC1D)  $(GW0tm)
	$(LK) $(LKFLAGS1) $(MLOC1D) $(GW0tm)  $(LKFLAGS2) -o $@

 hnocc_mlw:	$(hnocc_mlw) $(GW0tm)
	$(LK) $(LKFLAGS1) $(hnocc_mlw) $(GW0tm)  $(LKFLAGS2) -o $@

 hnocc_mlwt:	$(hnocc_mlwt) $(GW0tm)
	$(LK) $(LKFLAGS1) $(hnocc_mlwt) $(GW0tm)  $(LKFLAGS2) -o $@

 hpade:	$(PADE) $(GW0tm)
	$(LK) $(LKFLAGS1) $(PADE)    $(GW0tm)  $(LKFLAGS2) -o $@

 hphig:	$(PHIG) $(GW0tm)
	$(LK) $(LKFLAGS1) $(PHIG)     $(GW0tm)  $(LKFLAGS2) $(LIBSLA) -o $@

 hpsig: $(PSIG) $(GW0tm)
	$(LK) $(LKFLAGS1) $(PSIG)     $(GW0tm)  $(LKFLAGS2) -o $@

 hpsig_MPI:	$(PSIG_MPI) $(GW0tm)
	$(LK) $(LKFLAGS1) $(PSIG_MPI)	$(GW0tm)	$(LKFLAGS2) -o $@

 hpsig_test: $(PSIG_TEST) $(GW0tm)
	$(LK) $(LKFLAGS1) $(PSIG_TEST)     $(GW0tm)  $(LKFLAGS2) -o $@

 hpsipsi:	$(PSIPSI) $(GW0tm)
	$(LK) $(LKFLAGS1) $(PSIPSI)     $(GW0tm)  $(LKFLAGS2) -o $@

 hpsipsi_test:	$(PSIPSI_TEST) $(GW0tm)
	$(LK) $(LKFLAGS1) $(PSIPSI_TEST)     $(GW0tm)  $(LKFLAGS2) -o $@

 hqpband_ent_MPI:		$(QPBAND_ENT) $(GW0tm)
	 $(LK) $(LKFLAGS1) $(QPBAND_ENT) $(GW0tm) $(LKFLAGS2) -o $@

 hqpband_ent_MPI_test:		$(QPBAND_ENT_TEST) $(GW0tm)
	 $(LK) $(LKFLAGS1) $(QPBAND_ENT_TEST) $(GW0tm) $(LKFLAGS2) -o $@

 hqpe_mlw:		$(QPE_MLW) $(GW0tm)
	 $(LK) $(LKFLAGS1) $(QPE_MLW) $(GW0tm) $(LKFLAGS2) -o $@

 hsfp0w0:		$(SXCW0) $(GW0tm)  
	 $(LK) $(LKFLAGS1) $(SXCW0)    $(GW0tm)  $(LKFLAGS2) -o $@

 hsfp0_mlw_MPI:		$(SXC_MPI) $(GW0tm)  
	 $(LK) $(LKFLAGS1) $(SXC_MPI)   $(GW0tm)  $(LKFLAGS2) -o $@

 hs0tau:		$(SXCTAU) $(GW0tm)  
	 $(LK) $(LKFLAGS1) $(SXCTAU)    $(GW0tm)  $(LKFLAGS2) -o $@

 hs0tau_MPI:	$(SXCTAU_MPI) $(GW0tm)
	$(LK) $(LKFLAGS1) $(SXCTAU_MPI)	$(GW0tm)	$(LKFLAGS2) -o $@

 hupdt_ent_MPI:		$(UPDT_ENT) $(GW0tm)
	 $(LK) $(LKFLAGS1) $(UPDT_ENT) $(GW0tm) $(LKFLAGS2) -o $@

 hx0fp0_mlw:	$(X0mlw) $(GW0tm)
	$(LK) $(LKFLAGS1) $(X0mlw)     $(GW0tm)  $(LKFLAGS2) -o $@

 hx0fp0_mlw_MPI:	$(X0mlw_MPI) $(GW0tm)
	$(LK) $(LKFLAGS1) $(X0mlw_MPI)     $(GW0tm)  $(LKFLAGS2) -o $@

 hx0fp0_mlw_cfao:	$(X0mlw_cfao) $(GW0tm)
	$(LK) $(LKFLAGS1) $(X0mlw_cfao)     $(GW0tm)  $(LKFLAGS2) -o $@

 hx0tau:	$(X0tau) $(GW0tm)
	$(LK) $(LKFLAGS1) $(X0tau)     $(GW0tm)  $(LKFLAGS2) -o $@

 hx0tau_MPI:	$(X0tau_MPI) $(GW0tm)
	$(LK) $(LKFLAGS1) $(X0tau_MPI)     $(GW0tm)  $(LKFLAGS2) -o $@

 hx0tau_test:	$(X0tau_test) $(GW0tm)
	$(LK) $(LKFLAGS1) $(X0tau_test)     $(GW0tm)  $(LKFLAGS2) -o $@

 huumat:	$(UU2) $(GW0tm)
	$(LK) -o $@ $(LKFLAGS1) $(UU2)     $(GW0tm)  $(LKFLAGS2) 

 huumat_MPI:	$(UU2_MPI) $(GW0tm)
	$(LK) $(LKFLAGS1) $(UU2_MPI)     $(GW0tm)  $(LKFLAGS2) -o $@

 huumat_test:	$(UU2_TEST) $(GW0tm)
	$(LK) $(LKFLAGS1) $(UU2_TEST)     $(GW0tm)  $(LKFLAGS2) -o $@

 hwfch: $(WFCH) $(GW0tm)
	$(LK) $(LKFLAGS1) $(WFCH)     $(GW0tm)  $(LKFLAGS2) -o $@

 hwinfo: $(maxloc)/hwinfo.o $(GW0tm)
	$(LK) $(LKFLAGS1) $(maxloc)/hwinfo.o     $(GW0tm)  $(LKFLAGS2) -o $@

 hwmat:		$(WMAT) $(GW0tm)
	$(LK) $(LKFLAGS1) $(WMAT) $(COMMX)   $(GW0tm)  $(LKFLAGS2) -o $@

 hwmat_MPI:		$(WMAT_MPI) $(GW0tm)
	$(LK) $(LKFLAGS1) $(WMAT_MPI) $(COMM)   $(GW0tm)  $(LKFLAGS2) -o $@

 hwmat2:		$(WMAT2) $(GW0tm)
	$(LK) $(LKFLAGS1) $(WMAT2)    $(GW0tm)  $(LKFLAGS2) -o $@

 hwmatc:		$(WMATC) $(GW0tm)
	$(LK) $(LKFLAGS1) $(WMATC)    $(GW0tm)  $(LKFLAGS2) -o $@

 hwmatc_MPI:		$(WMATC_MPI) $(GW0tm)
	$(LK) $(LKFLAGS1) $(WMATC_MPI)    $(GW0tm)  $(LKFLAGS2) -o $@

 hwtmat:		$(WTMAT) $(GW0tm)
	$(LK) $(LKFLAGS1) $(WTMAT)    $(GW0tm)  $(LKFLAGS2) -o $@

 hwtmat_MPI:		$(WTMAT_MPI) $(GW0tm)
	$(LK) $(LKFLAGS1) $(WTMAT_MPI)    $(GW0tm)  $(LKFLAGS2) -o $@

 read_gw0:	$(maxloc)read_gw0.o  $(NFPLtot)
	$(LK) $(LKFLAGS1) $(maxloc)read_gw0.o  $(NFPLtot)  $(LKFLAGS2) -o $@

 read_gw02:	$(maxloc)read_gw02.o  $(NFPLtot)
	$(LK) $(LKFLAGS1) $(maxloc)read_gw02.o  $(NFPLtot)  $(LKFLAGS2) -o $@

 read_hmnk:	$(gwt)read_hmnk.o  $(GW0tm)
	$(LK) $(LKFLAGS1) $(gwt)read_hmnk.o     $(GW0tm)  $(LKFLAGS2) -o $@

 qg4gw_TMBAND:		$(QG_TM)  
	 $(LK) $(LKFLAGS1) $(QG_TM)  $(LKFLAGS2) -o $@

 qpwf:		$(maxloc)qpwf.o $(GW0tm)
	$(LK) $(LKFLAGS1) $(maxloc)qpwf.o $(GW0tm)  $(LKFLAGS2) -o $@


huumat.o: $(gwsrc)keyvalue.o $(gwsrc)readqg.o  $(gwsrc)readeigen.o $(gwsrc)rwbzdata.o  \
          $(gwsrc)genallcf_mod.o 

wanplot : ${OBJS2} wanplot.o ${GW0}
	$(LK) $(LKFLAGS1) $(maxloc)wanplot.o $(OBJS2) $(GW0) $(COMM) $(LKFLAGS2) -o $@


