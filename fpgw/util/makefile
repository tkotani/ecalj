### I think you don't needs to modify this file usually ###
#### This file is not machine-dependent. #####


# ---- Machine-specific compiler flags ---
include make.inc


nfpsrc  = $(HOME)/ecal/fpgw/nfpsrc/
gwsrc  = $(HOME)/ecal/fpgw/gwsrc/

#progs  = calj_interp caljdet2 caljdet2 caljdet calj_LFC calj_noLFC calj_LFC2 calj_search0
progs  = calj_interp_mat calj_nlfc_mat  calj_nlfc_mat_ferro
#calj_search0 ccc calj_s calj_fmat

#readtote readtote2

# lmf_exec
#progs  = hbasfp0 hvccfp0 hx0fp0 hsfp0 hef hqpe hchknw qg4gw gwinit heftet hmergewv hparainfo hbndout rdata4gw_v2  hx0fp0_fal hx0fp1 
progs2 =  $(progs) $(tags)TAGS checkmod 


#### You can choose these options. all is default.

all :$(progs) 

clean:  
	 rm -f  $(progs)

cleanall:  
	 rm -f  $(progs2) $(main)*.o $(gwsrc)*.o 

CALJ_s = \
$(gwsrc)matm.o \
calj_sumruletest.o

CALJ_interp_mat = \
$(gwsrc)matm.o \
calj_interp_mat.o

CALJD = \
$(gwsrc)matm.o \
calj_nlfc_mat.o

CALJDF = \
$(gwsrc)matm.o \
calj_nlfc_mat_ferro.o

calj_fmat = \
$(gwsrc)matm.o \
calj_fmat.o

CCC = \
$(gwsrc)matm.o \
calj_det.org.o 

CALJD2 = \
$(gwsrc)matm.o \
caljdet2.o

CALJS = \
$(gwsrc)matm.o \
calj_search0.o

CALJL = \
calj_LFC.o

CALJL2 = \
calj_LFC.o

CALJNL = \
calj_noLFC.o

READTOTE = \
readtote.o

READTOTE2 = \
readtote2.o

NFPLtot  = $(nfpsrc)diagcv2.o 

calj_nlfc_mat_ferro:	$(CALJDF) $(NFPLtot)
	 $(LK) $(LKFLAGS1) $(CALJDF) $(NFPLtot) $(LKFLAGS2) -o $@

calj_s:	$(CALJ_s) $(NFPLtot)
	 $(LK) $(LKFLAGS1) $(CALJ_s) $(NFPLtot) $(LKFLAGS2) -o $@

calj_interp_mat:	$(CALJ_interp_mat) $(NFPLtot)
	 $(LK) $(LKFLAGS1) $(CALJ_interp_mat) $(NFPLtot) $(LKFLAGS2) -o $@

calj_nlfc_mat:	$(CALJD) $(NFPLtot)
	 $(LK) $(LKFLAGS1) $(CALJD) $(NFPLtot) $(LKFLAGS2) -o $@

calj_fmat:	$(calj_fmat) $(NFPLtot)
	 $(LK) $(LKFLAGS1) $(calj_fmat) $(NFPLtot) $(LKFLAGS2) -o $@

ccc:	$(CCC) $(NFPLtot)
	 $(LK) $(LKFLAGS1) $(CCC) $(NFPLtot) $(LKFLAGS2) -o $@

caljdet2:	$(CALJD2) $(NFPLtot)
	 $(LK) $(LKFLAGS1) $(CALJD2) $(NFPLtot) $(LKFLAGS2) -o $@

calj_search0:	$(CALJS) $(NFPLtot)
	 $(LK) $(LKFLAGS1) $(CALJS) $(NFPLtot) $(LKFLAGS2) -o $@

calj_LFC:	$(CALJL) $(NFPLtot)
	 $(LK) $(LKFLAGS1) $(CALJL) $(NFPLtot) $(LKFLAGS2) -o $@
calj_LFC2:	$(CALJL2) $(NFPLtot)
	 $(LK) $(LKFLAGS1) $(CALJL2) $(NFPLtot) $(LKFLAGS2) -o $@
calj_noLFC:	$(CALJNL) $(NFPLtot)
	 $(LK) $(LKFLAGS1) $(CALJNL) $(NFPLtot) $(LKFLAGS2) -o $@

 readtote:	$(READTOTE) 
	 $(LK) $(LKFLAGS1) $(READTOTE) $(LKFLAGS2) -o $@


 readtote2:	$(READTOTE2) 
	 $(LK) $(LKFLAGS1) $(READTOTE2) $(LKFLAGS2) -o $@

# --- Make rules ---
.SUFFIXES:
.SUFFIXES: .F .o
#.SUFFIXES: .f .o .c1_o .c2_0 .c3_o .c4_o .F

#.F.o:
#	$(FC) $(FFLAGS) $*.F -c -o $*.o
#	etags $*.f -o $(tags)`echo $*.f| sed 's/..\///' | sed 's/\//-/g'`.tags

.F.o:
	$(FC) $(FFLAGS) $*.F -c -o $*.o
#	etags $*.f -o $(tags)`echo $*.f| sed 's/..\///' | sed 's/\//-/g'`.tags








