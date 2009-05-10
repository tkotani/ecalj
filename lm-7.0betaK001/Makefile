###################################
#choose macine dependent PLATFORM
###################################

#PLATFORM=gfortran
PLATFORM=ifort
#PLATFORM=ifort_mpi
#PLATFORM=ifort_mpik

# It will reads MAKEINC/Make.inc.$(PLATFORM).

########################################
#don't edit lines below 
########################################

obj_path=obj.$(PLATFORM)

include slatsm/make.inc
include subs/make.inc
include fp/make.inc
include gwd/make.inc

#################################
.PHONY: slatsm-lib gwd-lib fp-lib subs-lib  slatsm-dir gwd-dir fp-dir subs-dir libs progs clean cleanobj veryclean distclean check checkgw install 
#################################


fplotprogs := fplot plbnds pldos

ALL: libs
	$(MAKE) -f Makefile progs
	-[ -f slatsm/slatsm.a ] && (cd FPLOTdir; $(MAKE) PLATFORM=$(PLATFORM) ; mv -f $(fplotprogs) ../)
#	(cd slatsm; make -f Makefile SLATSM=$(SLATSM))
#	(cd subs; make -f Makefile SUBS=$(SUBS))
#	(cd fp; make -f Makefile FP=$(FP))
#	(cd gwd; make -f Makefile GWD=$(GWD))
#	(cd FPLOTdir; make;cp $(fplotprogs) ../)
#	make -f Makefile $(programs)

#	(cd nc; make -f Makefile )
#	(cd tb; make -f Makefile )


# Generate main routines
libs:
	$(MAKE) -f Makefile slatsm-lib
	$(MAKE) -f Makefile subs-lib
	$(MAKE) -f Makefile fp-lib
	$(MAKE) -f Makefile gwd-lib

slatsm-lib: slatsm-dir
subs-lib: subs-dir
fp-lib: fp-dir subs-lib
gwd-lib: gwd-dir subs-lib



########################################
# general rules
########################################

%.o:%.F
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

fp_obj=$(addprefix fp/$(obj_path)/,$(fp_list))
fp_obj_path=fp/$(obj_path)
fp-lib: $(fp_obj)  
fp-dir:
	[ -d $(fp_obj_path) ] || mkdir  $(fp_obj_path) 
$(fp_obj_path)/%.o: fp/%.F
	$(FC) $(FFLAGS) -c $< -o $@ 

gwd_obj=$(addprefix gwd/$(obj_path)/,$(gwd_list))
gwd_obj_path=gwd/$(obj_path)
gwd-lib: $(gwd_obj)  
gwd-dir:
	[ -d $(gwd_obj_path) ] || mkdir  $(gwd_obj_path) 
$(gwd_obj_path)/%.o: gwd/%.F
	$(FC) $(FFLAGS) -c $< -o $@ 

########################################
# machie dependent rules
# I must load Make.inc here
########################################

include MAKEINC/Make.inc.$(PLATFORM)


########################################
# Generate main routines
########################################
lmchk.o: lmv7.F
	$(FC) $(FFLAGS) -DLMCHK -DFP -c lmv7.F -o $@
#lm.o: lmv7.F
#	$(FC) $(FFLAGS) -DLM -DNC -c lmv7.F -o $@
# non colinear version
#lm.o: lmv7.F
#	$(FC) $(FFLAGS) -DLM -DNC -c lmv7.F -o $@
lmfa.o: lmfav7.F
	$(FC) $(FFLAGS) -DLMFA -c lmfav7.F -o $@
lmf.o: lmv7.F
	$(FC) $(FFLAGS) -DLMF  -c lmv7.F -o $@
lmdos.o: lmv7.F
	$(FC) $(FFLAGS) -DLMDOS -c lmv7.F -o $@
lmfgw.o: lmv7.F
	$(FC) $(FFLAGS) -DLMF -DLMFGWD -c lmv7.F -o $@
lm67.o: lm67.F
	$(FC) $(FFLAGS) -DLM67 -c lm67.F -o $@
lmscell.o: lmv7.F
	$(FC) $(FFLAGS) -DLMSCELL -c lmv7.F -o $@
#tbe.o: lmv7.F
#	$(FC) $(FFLAGS) -DTBE -c lmv7.F -o $@
lmf2gw.o: lmf2gw.F
	$(FC) $(FFLAGS) -c lmf2gw.F -o $@

########################################
# programs
########################################
lmchk:	lmchk.o $(LIBFP) $(LIBSUBS) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o $(LIBFP) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@
 
#lmscell: lmscell.o $(LIBSUBS) $(LIBSLA)
#	$(LK) $(LKFLAGS1) $@.o $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

lmdos:	lmdos.o $(LIBSUBS) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

 
lmfa:	lmfa.o $(LIBSUBS) $(LIBFP) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o $(LIBFP) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@
 
lmf:	lmf.o  $(LIBSUBS) $(LIBFP) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o $(LIBFP) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@
#lm:	lm.o  $(LIBSUBS) 
#	$(LK) $(LKFLAGS1) $@.o $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@
# non-colinear verion
#lm:	lm.o  $(LIBSUBS) 
#	$(LK) $(LKFLAGS1) $@.o $(LIBNC) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@
 
lmfgw:	lmfgw.o $(LIBSUBS) $(LIBFP) $(LIBGWD) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o $(LIBGWD) $(LIBFP) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

lm67:	lm67.o $(LIBSUBS) $(LIBFP) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o  $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

lmf2gw:	lmf2gw.o
	$(FC) $(LKFLAGS1) $@.o $(LKFLAGS2) -o $@

#tbe:	tbe.o $(LIBSUBS) $(LIBTB)
#	$(LK) $(LKFLAGS1) $@.o $(LIBTB) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

rdcmd:	rdcmd.o $(LIBSUBS) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@


####### for MPI
lmf-MPI.o: lmv7.F
	$(FC) $(FFLAGS) -DLMF -c lmv7.F -o $@

lmf-MPI:        lmf-MPI.o $(LIBSUBS) $(LIBFP) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o $(LIBFP) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

lmfgw-MPI.o: lmv7.F
	$(FC) $(FFLAGS) -DLMFGWD -c lmv7.F -o $@

lmfgw-MPI:      lmfgw-MPI.o $(LIBSUBS) $(LIBFP) $(LIBGWD) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o $(LIBGWD) $(LIBFP) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

####### for MPIK
lmf-MPIK.o: lmv7.F
	$(FC) $(FFLAGS) -DLMF -c lmv7.F -o $@

lmf-MPIK:       lmf-MPIK.o $(LIBSUBS) $(LIBFP) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o $(LIBFP) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

lmfgw-MPIK.o: lmv7.F
	$(FC) $(FFLAGS) -DLMFGWD -c lmv7.F -o $@

lmfgw-MPIK:     lmfgw-MPIK.o $(LIBSUBS) $(LIBFP) $(LIBGWD) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o $(LIBGWD) $(LIBFP) $(LIBSUBS) $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@


##########################################
# .PHONY targets
##########################################
progs: $(programs)

clean:
	rm -f $(programs) $(fplotprogs)
	(cd FPLOTdir; $(MAKE) PLATFORM=$(PLATFORM) clean)

cleanobj:
	rm -f $(programs) $(fplotprogs) $(LIBSUBS) $(LIBFP) $(LIBGWD) $(LIBSLA) */obj.*/*.o  *.o  */*.o  */*.a 
	(cd FPLOTdir; $(MAKE) PLATFORM=$(PLATFORM) cleanobj)
	rm -rf MODDIR
	mkdir MODDIR
	mkdir MODDIR/MPI
	mkdir MODDIR/MPIK
	touch MODDIR/dummy
	touch MODDIR/MPI/dummy
	touch MODDIR/MPIK/dummy
 
veryclean:
	for p in $(programs) ; do \
	  rm -f $$p $$p.o ; \
	done
	testing/test.lm --quiet --clean
	fp/test/test.fp --all --quiet --clean
	gwd/test/test.gwd --all --quiet --clean
 
distclean:
	for p in $(programs) ; do \
	  rm -f $$p $$p.o ; \
	done
	testing/test.lm --quiet --clean
	fp/test/test.fp --all --quiet --clean
	gwd/test/test.gwd --all --quiet --clean
	rm -f v7input/subs.a subs/subs.a fp/subs.a gwd/subs.a
	rm -f v7input/{m_gtv,m_rdctrl,m_rdctrlchk,m_toksw}.mod

check:
	(cd fp/test;./test.fp --all --quiet) 
	
checkgw:
	echo 'For GW driver, run samples in TESTinstallGW.'

#	testing/test.lm --quiet $(ADD0)


install: 
	cp $(programs) $(fplotprogs) ~/bin



