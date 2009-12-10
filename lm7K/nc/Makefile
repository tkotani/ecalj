lmdir = ../
include ../Make.inc

ALL:	$(NC)

.PRECIOUS:	$(NC)

$(NC):	$(NC)(amagnc.o) \
	$(NC)(bdotsr.o) \
	$(NC)(bndasa.o) \
	$(NC)(bsrhos.o) \
	$(NC)(ham2nc.o) \
	$(NC)(hamss2.o) \
	$(NC)(hmfr2c.o) \
	$(NC)(hmfr3c.o) \
	$(NC)(hml2nc.o) \
	$(NC)(hmladb.o) \
	$(NC)(hmltnc.o) \
	$(NC)(hmltso.o) \
	$(NC)(lmasa.o) \
	$(NC)(magtrq.o) \
	$(NC)(mkbdia.o) \
	$(NC)(mksod.o) \
	$(NC)(mksoph.o) \
	$(NC)(mmag.o) \
	$(NC)(mmdyn.o) \
	$(NC)(mmpair.o) \
	$(NC)(rotspn.o) \
	$(NC)(rotspu.o) \
	$(NC)(secm2c.o) \
	$(NC)(secmat.o) \
	$(NC)(yhmlso.o)
	$(RANLIB) $@
 
clean: 
	-rm -f $(NC) Makefile Names
