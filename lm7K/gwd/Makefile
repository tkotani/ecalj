lmdir = ../
include ../Make.inc

ALL:	$(GWD)

.PRECIOUS:	$(GWD)


$(GWD):	$(GWD)(bndconn_v2.o) \
	$(GWD)(bndfp.o) \
	$(GWD)(chimedit.o) \
	$(GWD)(genqbz.o) \
	$(GWD)(gwcphi.o) \
	$(GWD)(lmfp.o) \
	$(GWD)(polinta.o) \
	$(GWD)(pwmat.o) \
	$(GWD)(stonerpb.o) \
	$(GWD)(stonerrsa.o) \
	$(GWD)(sugw.o) \
	$(GWD)(sugwin.o) \
	$(GWD)(suq0x.o)
	$(RANLIB) $@
 
clean: 
	-rm -f $(GWD) 

