#lmdir= ../


# Allow configure to specify the install program.
INSTALL = /usr/bin/install -c
INSTALL_PROGRAM = ${INSTALL}
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_INFO = install-info

# Not used now
infodir = ${prefix}/info
mandir = ${prefix}/man

# We don't need high speed here ... avoid compiler bugs by using lower default optimization
FFLAGS=$(FFLAGS_LESS)

PROGS = fplot pldos plbnds
#SUBS = ./subs.a

ALL:	$(PROGS)

fplot:	fplot.o fpsub.o plsub.o contour.o alloc.o
	$(LK) $(LKFLAGS1) $@.o fpsub.o plsub.o contour.o alloc.o $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

pldos:	pldos.o plsub.o alloc.o
	$(LK) $(LKFLAGS1) $@.o plsub.o  alloc.o $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

plbnds:	plbnds.o plsub.o alloc.o
	$(LK) $(LKFLAGS1) $@.o plsub.o  alloc.o $(LIBSLA) $(LIBLOC) $(LKFLAGS2) -o $@

progs:
	make -f Makefile fplot pldos plbnds
 
tar:	
	rm -r -f plot; mkdir plot
	cp -rp README configure Makefile.in etc examples fplot.f plbnds.f pldos.f fpsub.f plsub.f contour.f plot
	tar --gzip -cvf FPLOT.tar.gz plot
	rm -r -f plot
 
clean:	
	rm -f $(PROGS) 
cleanobj:
	rm -f $(PROGS)  *.o 

veryclean:
	for p in  fplot pldos plbnds; do \
	  rm -f $$p $$p.o ; \
	done
	rm -f contour.o fpsub.o plsub.o
	rm -f plot.dos dos.cr ps.dos.cr
	rm -f plot.plbnds bnds.cr ps.cr bnd1.dat bnd2.dat bnd3.dat bnd4.dat bnd5.dat bnd6.dat
	rm -f plot.dos plot.plbnds out.fplot out.plbnds out.pldos

check:
	 examples/test.fplot  --quiet $(ADD0)
	 examples/test.pldos  --quiet $(ADD0)
	 examples/test.plbnds --quiet $(ADD0)

 
include ../MAKEINC/Make.inc.$(PLATFORM)
# use serial version, overwrite LIBSLA 
LIBSLA=../slatsm/slatsm.a

