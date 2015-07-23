############ ifort ##########
#FC=mpiifort
#LK=mpiifort
#LKFLAGS2=-mkl

top=../
PLATFORM=gfortran
include ../exec/make.inc.$(PLATFORM)

FFI=-I${top}/exec
BINDIR=${HOME}/bin
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
gwsrc   = ${top}/gwsrc/
main    = ${top}/main/
nfpsrc  = ${top}/nfpsrc/
slatsmlib  =${top}/../lm7K/slatsm/
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
	 cp  $(progs) genMLWF genMLWF2 RPAWannier* $(BINDIR)

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
 

NFPLtot  = $(nfpsrc)diagcv2.o 

GW0tot = \
$(gwsrc)m_hamindex.o\
$(gwsrc)MPI_fpgw2.o \
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


################
#
# BNDCONN= \
# $(gwsrc)bndconn.o   ### This is not linked but bndconn.o is used in lm/lmfgw. 
# It is now included in lm/gw/
 DERFC=   $(nfpsrc)derfc.o \
          $(nfpsrc)i1mach.o \
          $(nfpsrc)d1mach.o 



 GW0 = \
 $(gwsrc)m_hamindex.o\
 $(gwsrc)MPI_fpgw2.o \
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


 GW0tm = \
 $(gwsrc)m_hamindex.o\
 $(gwsrc)MPI_fpgw2.o \
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

 PHIG = \
 $(maxloc)hphig.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt5.o \
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

 PSIG = \
 $(maxloc)hpsig.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt5.o \
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

 PSIG_MPI = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_qkgroup_mod.o \
 $(maxloc)RSMPI_rotkindex_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)hpsig_MPI.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt5.o \
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


 # PSIPSI = \
 # $(gwt)hpsipsi.o \
 # $(gwsrc)gintxx.o \
 # $(gwsrc)pplmat.o \
 # $(gwsrc)getgv2.o \
 # $(gwsrc)rs.o \
 # $(gwsrc)mopen.o \
 # $(nfpsrc)u_lat_0.o \
 # $(nfpsrc)wronkj.o \
 # $(nfpsrc)mklegw.o \
 # $(nfpsrc)rxx.o \
 # $(nfpsrc)bessl.o \
 # $(nfpsrc)cross.o 

 UU2 = \
 $(maxloc)huumat.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt5.o \
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

 UU2_MPI = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)huumat_MPI.o \
 $(gwsrc)wcf.o \
 $(gwsrc)tetwt5.o \
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
 $(nfpsrc)cross.o \
 $(maxloc)ent_MPI.o

 WMATK = \
 $(maxloc)hwmatK.o \
 $(gwsrc)wse.o \
 $(maxloc)wmatK.o \
 $(maxloc)maxloc0.o \
 $(gwsrc)genallcf_dump.o

 WMATK_MPI = \
 $(maxloc)RSMPI_mod.o \
 $(maxloc)RSMPI_qkgroup_mod.o \
 $(maxloc)RSMPI_rotkindex_mod.o \
 $(maxloc)RSMPI_utils.o \
 $(maxloc)hwmatK_MPI.o \
 $(gwsrc)wse.o \
 $(maxloc)wmatK_mpi.o \
 $(maxloc)maxloc0.o \
 $(gwsrc)genallcf_dump.o \
 $(maxloc)efsimplef_RSMPI.o

COMM=$(nfpsrc)rxx.o \
 $(gwsrc)mopen.o \
 $(gwsrc)rangedq.o

COMMX=$(nfpsrc)rxx.o \
 $(gwsrc)mopen.o

 hmaxloc:	$(MLOC)  $(GW0tm)
	$(LK) $(LKFLAGS1) $(MLOC) $(GW0tm)  $(LKFLAGS2) -o $@

 hpsig: $(PSIG) $(GW0tm)
	$(LK) $(LKFLAGS1) $(PSIG)     $(GW0tm)  $(LKFLAGS2) -o $@

 hpsig_MPI:	$(PSIG_MPI) $(GW0tm)
	$(LK) $(LKFLAGS1) $(PSIG_MPI)	$(GW0tm)	$(LKFLAGS2) -o $@

 huumat:	$(UU2) $(GW0tm)
	$(LK) -o $@ $(LKFLAGS1) $(UU2)     $(GW0tm)  $(LKFLAGS2) 

 huumat_MPI:	$(UU2_MPI) $(GW0tm)
	$(LK) $(LKFLAGS1) $(UU2_MPI)     $(GW0tm)  $(LKFLAGS2) -o $@

 hwmatK:		$(WMATK) $(GW0tm)
	$(LK) $(LKFLAGS1) $(WMATK) $(COMMX)   $(GW0tm)  $(LKFLAGS2) -o $@

 hwmatK_MPI:		$(WMATK_MPI) $(GW0tm)
	$(LK) $(LKFLAGS1) $(WMATK_MPI) $(COMMX)   $(GW0tm)  $(LKFLAGS2) -o $@


huumat.o: $(gwsrc)keyvalue.o $(gwsrc)readqg.o  $(gwsrc)readeigen.o $(gwsrc)rwbzdata.o  \
          $(gwsrc)genallcf_mod.o 

wanplot : ${OBJS2} wanplot.o ${GW0}
	$(LK) $(LKFLAGS1) $(maxloc)wanplot.o $(OBJS2) $(GW0) $(COMM) $(LKFLAGS2) -o $@

htbplot :  htbplot.o
	$(LK) $(LKFLAGS1) $(maxloc)htbplot.o $(LKFLAGS2) -o $@


