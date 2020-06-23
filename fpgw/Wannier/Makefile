############ ifort ##########
#FC=mpiifort
#LK=mpiifort
#LKFLAGS2=-mkl

top=../
PLATFORM=ifort
include ../exec/make.inc.$(PLATFORM)

FFI=-I${top}/exec
BINDIR=${HOME}/ecaljdev/bin
########### gfortran #############
#FC=mpif90
#LK=mpif90
#LKFLAGS2= /usr/lib/libfftw3.so.3 /usr/lib/liblapack.so.3gf /usr/lib/libblas.so.3gf 
#FFLAGS = -ffixed-line-length-132 -c -frange-check -I${top}/exec


### I think you don't need to modify following lines. Not machine-dependent ###########
x0kf_g=.o
para_g=.o

# ---- Machine-specific compiler flags ---
#include make.inc

#-------------------------------------------------------
# src directories
main    = ${top}/main/
gwsrc  = ${top}/gwsrc/
tote = ./
# tag directory
tags   = ./

#progs  = hbasfp0 hvccfp0 hx0fp0 hsfp0 hef hqpe hchknw qg4gw gwinit heftet hmergewv hparainfo hbndout rdata4gw_v2 convgwin hx0fp0_sc hsfp0_sc hqpe_sc kino_input_test hecor eout eout2 h_uumatrix hsigmconv
progs= huumat hpsig hmaxloc huumat_MPI hpsig_MPI wanplot hwmatK hwmatK_MPI htbplot # hx0fp0_mlw 
progs2 = $(progs) 
# lmf_exec
#progs  = hbasfp0 hvccfp0 hx0fp0 hsfp0 hef hqpe hchknw qg4gw gwinit heftet hmergewv hparainfo hbndout rdata4gw_v2  hx0fp0_fal hx0fp1 
#progs2 = $(progs) $(tags)TAGS checkmod 
#progs_rs = $(progs_RSMPI) $(progs_RSBAND)
#progs_RSMPI = hvccfp0_RSMPI hx0fp0_RSMPI hsfp0_RSMPI
#progs_RSBAND = qg4gw_RSBAND

#script = cleargw* dqpu dtote eps* ex* gw* hqpemetal* inf* lmgw* plotg save* tote_lmfh2 xqp

#### You can choose these options. all is default.

all :$(progs2) 
#all_rs	:$(progs_rs)

#clean:  
#	 rm -f  $(progs)
clean:  
	 rm -f  $(progs2) *.o 

install:  
	 cp  $(progs) genMLWF genMLWF_vw genMLWF2 genMLWFmod run_argloc RPAWannier* FLEX_interaction.py $(BINDIR) 

#install2:  
#	 cp  $(script) $(HOME)/bin

cleanall:  
	 rm -f  $(progs2) $(maxloc)/*.o $(maxloc)/*.mod

#	 rm -f  $(progs2) $(main)*.o $(gwsrc)*.o  $(maxloc)*.o  $(gwt)*.o 

# This is necesaly to compile *.f in right order.
# When you recompile and link, just repeat 'make' (not necessary to repeat 'make init').
# When checkmodule recompile source, you have to repeat 'make'.
init:
	exec ../../TOOLS/checkmodule ../gwsrc/*.F *.F 

checkmod:
	exec ../../TOOLS/checkmodule ../gwsrc/*.F *.F 

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
 

NFPLtot  = $(gwsrc)diagcv2.o 

# GW0tot = \
# $(gwsrc)m_hamindex.o\
# $(gwsrc)MPI_fpgw2.o \
# $(gwsrc)alloclist.o \
# $(gwsrc)rwbzdata.o \
# $(gwsrc)keyvalue.o \
# $(gwsrc)genallcf_mod.o \
# $(gwsrc)rgwinf_mod.o \
# $(gwsrc)nocctotg.o \
# $(gwsrc)ppbafp.fal$(para_g) \
# $(gwsrc)psi2b_v2$(para_g) \
# $(gwsrc)psi2b_v3$(para_g) \
# $(gwsrc)wfacx.o \
# $(gwsrc)sortea.o \
# $(gwsrc)rydberg.o \
# $(gwsrc)polinta.o \
# $(gwsrc)efsimplef.o \
# $(gwsrc)extension.o \
# $(gwsrc)nword.o \
# $(gwsrc)scg.o \
# $(gwsrc)matm.o \
# $(gwsrc)rdpp.o \
# $(gwsrc)mptauof.o \
# $(gwsrc)rotdlmm.o \
# $(gwsrc)iopen.o \
# $(gwsrc)cputid.o \
# $(gwsrc)rw.o \
# $(gwsrc)ext.o \
# $(gwsrc)ext2.o \
# $(gwsrc)cross.o \
# $(gwsrc)mate.o \
# $(gwsrc)mate1.o \
# $(gwsrc)icopy.o \
# $(gwsrc)bib1.o \
# $(gwsrc)index.o \
# $(gwsrc)idxk.o \
# $(gwsrc)maxnn.o \
# $(gwsrc)reindx.o \
# $(gwsrc)iprint.o \
# $(gwsrc)bz.o \
# $(gwsrc)bzmesh.o \
# $(gwsrc)genqbz.o \
# $(gwsrc)switches.o \
# $(gwsrc)linpackdummy.o \
# $(gwsrc)rppovl.o \
# $(gwsrc)llnew.o

#LMFtot= \
#$(lmsrc)lmf.o \
#$(tote)mkpot_exec.o \
#$(tote)locpot_exec.o \
#$(tote)vxcnsp_exec.o \
#$(tote)evxcv_exec.o \
#$(tote)mkehkf_exec.o \
#$(tote)smvxcm_exec.o


################
#
# BNDCONN= \
# $(gwsrc)bndconn.o   ### This is not linked but bndconn.o is used in lm/lmfgw. 
# It is now included in lm/gw/
 # DERFC=   $(gwsrc)derfc.o \
 #          $(gwsrc)i1mach.o \
 #          $(gwsrc)d1mach.o 



 # GW0 = \
 # $(gwsrc)m_hamindex.o\
 # $(gwsrc)MPI_fpgw2.o \
 # $(gwsrc)alloclist.o \
 # $(gwsrc)readpomat.o \
 # $(gwsrc)keyvalue.o \
 # $(gwsrc)rppovl.o \
 # $(gwsrc)nocctotg.o \
 # $(gwsrc)ppbafp.fal$(para_g) \
 # $(gwsrc)psi2b_v2$(para_g) \
 # $(gwsrc)psi2b_v3$(para_g) \
 # $(gwsrc)wfacx.o \
 # $(gwsrc)sortea.o \
 # $(gwsrc)rydberg.o \
 # $(gwsrc)polinta.o \
 # $(gwsrc)efsimplef.o \
 # $(gwsrc)extension.o \
 # $(gwsrc)nword.o \
 # $(gwsrc)scg.o \
 # $(gwsrc)matm.o \
 # $(gwsrc)rdpp.o \
 # $(gwsrc)mptauof.o \
 # $(gwsrc)genallcf_mod.o \
 # $(gwsrc)rgwinf_mod.o \
 # $(gwsrc)rotdlmm.o \
 # $(gwsrc)iopen.o \
 # $(gwsrc)cputid.o \
 # $(gwsrc)rw.o \
 # $(gwsrc)ext.o \
 # $(gwsrc)ext2.o \
 # $(gwsrc)cross.o \
 # $(gwsrc)mate.o \
 # $(gwsrc)mate1.o \
 # $(gwsrc)icopy.o \
 # $(gwsrc)bib1.o \
 # $(gwsrc)index.o \
 # $(gwsrc)idxk.o \
 # $(gwsrc)maxnn.o \
 # $(gwsrc)reindx.o \
 # $(gwsrc)iprint.o \
 # $(gwsrc)bz.o \
 # $(gwsrc)bzmesh.o \
 # $(gwsrc)genqbz.o \
 # $(gwsrc)linpackdummy.o \
 # $(gwsrc)switches.o \
 # $(gwsrc)rwbzdata.o \
 # $(gwsrc)llnew.o  \
 # $(gwsrc)readeigen.o \
 # $(gwsrc)readqg.o \
 # $(gwsrc)iqindx.o \
 # $(NFPLtot)


 # QPE = \
 # $(gwsrc)switches.o \
 # $(gwsrc)keyvalue.o \
 # $(main)hqpe.m$(hqpe_g) \
 # $(gwsrc)qpe1.o \
 # $(gwsrc)icompvv2.o \
 # $(gwsrc)iopenxx.o \
 # $(gwsrc)iopen.o \
 # $(gwsrc)rw.o \
 # $(gwsrc)rydberg.o


 # NFPL  = $(gwsrc)wronkj.o \
 #         $(gwsrc)sylm.o \
 #         $(gwsrc)sylmnc.o \
 #         $(gwsrc)u_lat_0.o \
 #         $(gwsrc)mklegw.o \
 #         $(gwsrc)setpr.o \
 #         $(gwsrc)bessl.o \
 #         $(gwsrc)rxx.o \
 #         $(gwsrc)hsmq.o \
 #         $(gwsrc)lgen.o \
 #         $(gwsrc)hansr5.o \
 #         $(gwsrc)hansr4.o \
 #         $(gwsrc)lattc.o \
 #         $(gwsrc)ll.o \
 #         $(gwsrc)dpcopy.o \
 #         $(gwsrc)dpadd.o \
 #         $(gwsrc)syscalls.o \
 #         $(gwsrc)qdist.o \
 #         $(gwsrc)dlmtor.o \
 #         $(gwsrc)dpzero.o \
 #         $(gwsrc)ropyln.o \
 #         $(gwsrc)ropcsm.o \
 #         $(gwsrc)dsisl.o \
 #         $(gwsrc)dsifa.o \
 #         $(gwsrc)diagcv2.o \
 #         $(gwsrc)scg.o 


 $(tags)TAGS: $(progs)
	cd $(tags);etags ./*/*.f ./*/*.F


# --- Make rules ---
.SUFFIXES:
.SUFFIXES: .F .o
#.SUFFIXES: .f .o .c1_o .c2_0 .c3_o .c4_o .F

.F.o:
	$(FC) $(FFI) $(FFLAGS) $*.F -c -o $*.o
#	etags $*.f -o $(tags)`echo $*.f| sed 's/..\///' | sed 's/\//-/g'`.tags

#.F.o:
#	$(FC) $(FFLAGS) $*.F -c -o $*.o
#	etags $*.f -o $(tags)`echo $*.f| sed 's/..\///' | sed 's/\//-/g'`.tags

#.f.o:
#	$(FC) $(FFLAGS) $*.f -c -o $*.o
#	etags $*.f -o $(tags)`echo $*.f| sed 's/..\///' | sed 's/\//-/g'`.tags

.f.c1_o:
	$(FC) $(FFI) $(FFLAGS_c1) $*.f -c -o $*.c1_o
	etags $*.f -o $(tags)`echo $*.f| sed 's/..\///' | sed 's/\//-/g'`.tags

.f.c2_o:
	$(FC) $(FFI) $(FFLAGS_c2) $*.f -c -o $*.c2_o
	etags $*.f -o $(tags)`echo $*.f| sed 's/..\///' | sed 's/\//-/g'`.tags

.f.c3_o:
	$(FC) $(FFI) $(FFLAGS_c3) $*.f -c -o $*.c3_o
	etags $*.f -o $(tags)`echo $*.f| sed 's/..\///' | sed 's/\//-/g'`.tags

.f.c4_o:
	$(FC) $(FFI) $(FFLAGS_c4) $*.f -c -o $*.c4_o
	etags $*.f -o $(tags)`echo $*.f| sed 's/..\///' | sed 's/\//-/g'`.tags

include moduledepends.inc

##### Maxloc Wannier fn. and Finite-temperature GW
### Takashi MIYAKE
maxloc = ./
#gwt = ../Miyake/gwt/
#misc = ../Miyake/misc/

#progs_tm = hwmat hwmatc hmaxloc huumat qpwf hpsig hnocc_mlw hx0fp0_mlw hx0fp0_mlw_ent hphig hx0tau hs0tau hpsipsi hibz2fbz hhmnk hwtmat hnocc_mlwt hwinfo hpade hibz2fbz0 read_gw0 # hdyson # hmaxloc1D
#progs_TM= hmaxloc huumat qpwf hpsig hnocc_mlw hx0fp0_mlw hwmat hwmatc hx0tau hs0tau hpsipsi hibz2fbz hhmnk hwtmat hnocc_mlwt hwinfo hwfch hibz2fbz0 read_gw0 qg4gw_TMBAND hqpe_mlw hibz2fbz0wan # hdyson hphig hpade haw
#progs_TMMPI= hx0tau_MPI hs0tau_MPI hwtmat_MPI hwmatK_MPI hpsig_MPI huumat_MPI hx0fp0_mlw_MPI hwmatc_MPI hsfp0_mlw_MPI hband_ent_MPI hqpband_ent_MPI hupdt_ent_MPI

#progs_tm = $(progs_TM) $(progs_TMMPI)
#all_tm:	$(progs_tm)

 GWLIB =  \
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
 $(gwsrc)m_tetrakbt.o \
 $(gwsrc)tetwt5$(tet5_g) \
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
$(gwsrc)rotgvec.o \
$(gwsrc)m_readVcoud.o\
$(gwsrc)m_itq.o \
$(gwsrc)m_qbze.o \
$(gwsrc)nan.o\
$(gwsrc)m_readgwinput.o\
$(gwsrc)getwemax.o\
$(MPI)

 MPI =  $(gwsrc)MPI_fpgw2.o

 # GW0tm = \
 # $(gwsrc)m_hamindex.o\
 # $(gwsrc)MPI_fpgw2.o \
 # $(gwsrc)alloclist.o \
 # $(gwsrc)readpomat.o \
 # $(gwsrc)keyvalue.o \
 # $(gwsrc)rppovl.o \
 # $(gwsrc)nocctotg.o \
 # $(gwsrc)ppbafp.fal$(para_g) \
 # $(gwsrc)psi2b_v2$(para_g) \
 # $(gwsrc)psi2b_v3$(para_g) \
 # $(gwsrc)wfacx.o \
 # $(gwsrc)sortea.o \
 # $(gwsrc)rydberg.o \
 # $(gwsrc)polinta.o \
 # $(gwsrc)efsimplef.o \
 # $(gwsrc)extension.o \
 # $(gwsrc)nword.o \
 # $(gwsrc)scg.o \
 # $(gwsrc)matm.o \
 # $(gwsrc)rdpp.o \
 # $(gwsrc)mptauof.o \
 # $(gwsrc)genallcf_mod.o \
 # $(gwsrc)rgwinf_mod.o \
 # $(gwsrc)rotdlmm.o \
 # $(gwsrc)iopen.o \
 # $(gwsrc)cputid.o \
 # $(gwsrc)rw.o \
 # $(gwsrc)ext.o \
 # $(gwsrc)ext2.o \
 # $(gwsrc)cross.o \
 # $(gwsrc)mate.o \
 # $(gwsrc)mate1.o \
 # $(gwsrc)icopy.o \
 # $(gwsrc)bib1.o \
 # $(gwsrc)index.o \
 # $(gwsrc)idxk.o \
 # $(gwsrc)maxnn.o \
 # $(gwsrc)reindx.o \
 # $(gwsrc)iprint.o \
 # $(gwsrc)bz.o \
 # $(gwsrc)bzmesh.o \
 # $(gwsrc)genqbz.o \
 # $(gwsrc)linpackdummy.o \
 # $(gwsrc)switches.o \
 # $(gwsrc)rwbzdata.o \
 # $(gwsrc)llnew.o  \
 # $(gwsrc)readeigen.o \
 # $(gwsrc)readqg.o \
 # $(gwsrc)iqindx.o \
 # $(gwsrc)rangedq.o \
 # $(NFPLtot)

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

COMM=

COMMX=

 hmaxloc:	$(MLOC)  $(GWLIB)
	$(LK) $(LKFLAGS1) $(MLOC) $(GWLIB)  $(LKFLAGS2) -o $@

 hpsig: $(PSIG) $(GWLIB)
	$(LK) $(LKFLAGS1) $(PSIG)     $(GWLIB)  $(LKFLAGS2) -o $@

 hpsig_MPI:	$(PSIG_MPI) $(GWLIB)
	$(LK) $(LKFLAGS1) $(PSIG_MPI)	$(GWLIB)	$(LKFLAGS2) -o $@

 huumat:	$(UU2) $(GWLIB)
	$(LK) -o $@ $(LKFLAGS1) $(UU2)     $(GWLIB)  $(LKFLAGS2) 

 huumat_MPI:	$(UU2_MPI) $(GWLIB)
	$(LK) $(LKFLAGS1) $(UU2_MPI)     $(GWLIB)  $(LKFLAGS2) -o $@

 hwmatK:		$(WMATK) $(GWLIB)
	$(LK) $(LKFLAGS1) $(WMATK) $(COMMX)   $(GWLIB)  $(LKFLAGS2) -o $@

 hwmatK_MPI:		$(WMATK_MPI) $(GWLIB)
	$(LK) $(LKFLAGS1) $(WMATK_MPI) $(COMMX)   $(GWLIB)  $(LKFLAGS2) -o $@


huumat.o: $(gwsrc)keyvalue.o $(gwsrc)readqg.o  $(gwsrc)readeigen.o $(gwsrc)rwbzdata.o  \
          $(gwsrc)genallcf_mod.o 

wanplot : ${OBJS2} wanplot.o ${GWLIB}
	$(LK) $(LKFLAGS1) $(maxloc)wanplot.o $(OBJS2) $(GWLIB) $(COMM) $(LKFLAGS2) -o $@

htbplot :  htbplot.o
	$(LK) $(LKFLAGS1) $(maxloc)htbplot.o $(LKFLAGS2) -o $@


