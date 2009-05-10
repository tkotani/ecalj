# ---- Machine-specific compiler flags ---
lmdir = ./
include Make.inc

ALL:
	(cd slatsm; make -f Makefile SLATSM=$(SLATSM))
	(cd subs; make -f Makefile SUBS=$(SUBS))
	(cd fp; make -f Makefile FP=$(FP))
	(cd gwd; make -f Makefile GWD=$(GWD))
	(cd FPLOTdir; make;cp fplot plbnds pldos ../)
	make -f Makefile $(programs)

#	(cd nc; make -f Makefile )
#	(cd tb; make -f Makefile )


# Generate main routines axxxxx
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

# Main routines
lmchk:	lmchk.o $(LIBFP) $(LIBSUBS) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o $(LIBFP) $(LIBSUBS) $(LIBES) $(LKFLAGS2) -o $@
 
#lmscell: lmscell.o $(LIBSUBS) $(LIBSLA)
#	$(LK) $(LKFLAGS1) $@.o $(LIBSUBS) $(LIBES) $(LKFLAGS2) -o $@

lmdos:	lmdos.o $(LIBSUBS) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o $(LIBSUBS) $(LIBES) $(LKFLAGS2) -o $@
 
lmfa:	lmfa.o $(LIBSUBS) $(LIBFP) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o $(LIBFP) $(LIBSUBS) $(LIBES) $(LKFLAGS2) -o $@
 
lmf:	lmf.o  $(LIBSUBS) $(LIBFP) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o $(LIBFP) $(LIBSUBS) $(LIBES) $(LKFLAGS2) -o $@
#lm:	lm.o  $(LIBSUBS) 
#	$(LK) $(LKFLAGS1) $@.o $(LIBSUBS) $(LIBES) $(LKFLAGS2) -o $@
# non-colinear verion
#lm:	lm.o  $(LIBSUBS) 
#	$(LK) $(LKFLAGS1) $@.o $(LIBNC) $(LIBSUBS) $(LIBES) $(LKFLAGS2) -o $@
 
lmfgw:	lmfgw.o $(LIBSUBS) $(LIBFP) $(LIBGWD) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o $(LIBGWD) $(LIBFP) $(LIBSUBS) $(LIBES) $(LKFLAGS2) -o $@

lm67:	lm67.o $(LIBSUBS) $(LIBFP) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o  $(LIBSUBS) $(LIBES) $(LKFLAGS2) -o $@

lmf2gw:	lmf2gw.o
	$(FC) $(LKFLAGS1) $@.o $(LKFLAGS2) -o $@

#tbe:	tbe.o $(LIBSUBS) $(LIBTB)
#	$(LK) $(LKFLAGS1) $@.o $(LIBTB) $(LIBSUBS) $(LIBES) $(LKFLAGS2) -o $@

rdcmd:	rdcmd.o $(LIBSUBS) $(LIBSLA)
	$(LK) $(LKFLAGS1) $@.o $(LIBES) $(LKFLAGS2) -o $@

progs:
	make -f Makefile $(programs)
 
clean:
	rm -f $(programs) pldos plbnds fplot 

cleanobj:
	rm -f $(programs) $(LIBSUBS) $(LIBFP) $(LIBGWD) $(LIBSLA) *.o 
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
	cp $(programs) pldos plbnds fplot ~/bin