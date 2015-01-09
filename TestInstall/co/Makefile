include ../Makefile.define

TARGET=co

out1=out.lmf.$(TARGET)
out2=out.lmf-dos.$(TARGET)
out2x=pdos.$(TARGET)

all: help

help:
	@echo "Do not call this Makefile alone."

###########################################################
#----- subroutine test1.exec
define test1.exec
	#-------------------------------------------
	# Case co: a hexagonal environment with two equivalent atoms.
	#-------------------------------------------
	#
	# --- Test 1.  Basic check of programs lmfa,lmf ---
	#
	@cat message.test1.$(TARGET)
	$(bindir)/lmfa  $(TARGET) -vmet=4 -vlmf=1 -vnk=8 -vnit=3 --pr31 \
		--no-iactiv > $1 2>&1
	$(LMFP)   $(TARGET) -vmet=4 -vlmf=1 -vnk=8 -vnit=3 -vw2=0 \
		--pr31 --no-iactiv --rs=0 >> $1 2>&1
	rm mixm.$(TARGET)
	$(LMFP)   $(TARGET) -vmet=4 -vlmf=1 -vnk=8 -vnit=3 -vw1=0 \
		--pr31 --no-iactiv >> $1 2>&1
	rm mixm.$(TARGET)
	$(LMFP)   $(TARGET) -vmet=4 -vlmf=1 -vnk=8 -vnit=3 -vw1=1 \
		--pr31 --no-iactiv --time=5 >> $1 2>&1
#	$(LMF)   $(TARGET) -vmet=4 -vlmf=1 -vnk=8 -vnit=3 --pr31 \
#		--no-iactiv -vso=t \
#		--band~col=5:9,14:18~col2=18+5:18+9,18+14:18+18~fn=syml >> $1 2>&1
endef

#----- subroutine test2.exec
define test2.exec
	#
	# --- Test 2.  Core-level spectroscopy (EELS), Mulliken analysis, partial DOS ---
	#
	@cat message.test2.$(TARGET)
	$(bindir)/lmfa  $(TARGET) -vmet=4 -vlmf=1 -vnk=8 -vnit=1 \
		--pr31 --no-iactiv > $1 2>&1
	$(LMF)   $(TARGET) -vmet=4 -vlmf=1 -vnk=8 -vnit=1 \
		--pr31 --no-iactiv -vpdos=t --pdos:mode=2:sites=2 >> $1 2>&1
	mv dos.$(TARGET) tdos.$(TARGET)
	echo 501 -.8 .7 | $(bindir)/lmdos $(TARGET) -vmet=4 -vlmf=1 -vnk=8 \
		-vnit=1  --pr31 --no-iactiv -vpdos=t --pdos:mode=2:sites=2 >> $1 2>&1
	mv dos.$(TARGET) pdos.$(TARGET)
endef

#----- subroutine plot.plbnds
define plot.plbnds
	echo -10,5,5,10 | $(bindir)/plbnds -fplot -ef=0 -scl=13.6 \
		-lt=1,col=0,0,1,colw=1,0,0,colw2=0,1,0 -lbl=M,G,A,L,G,K \
		bnds.$(TARGET) > plot.log 2>&1
	$(bindir)/fplot -disp -f plot.plbnds >> plot.log 2>&1
	mv ps.dat ps.plbnds.dat
	@echo "[note] run make plot-band to plot band"
endef

#----- subroutine plot.pldos
define plot.pldos
	echo 30 15 -.8 .7 | $(bindir)/pldos '-lst=1;3,5,7;9,11,13,15,17' \
		'-lst2=2;4,6,8;10,12,14,16,18' -fplot \
		pdos.$(TARGET) > plot.log 2>&1
	$(bindir)/fplot -disp -pr10 -f plot.dos >> plot.log 2>&1
	mv ps.dat ps.pldos.dat
	@echo "[note] run make plot-dos to plot dos"
endef


###########################################################
# exec default case for $(TARGET).
$(TARGET).exec:
	$(call test1.exec,$(out1))
	$(call test1.check,$(REFDIR)/$(out1),$(out1))
	rm -f atm.* mixm.* rst.* save.* log.* hssn.* wkp.* dos.* tdos.* pdos.* dos-mull.* qpp.* out.lmf-dos*
	$(call test2.exec,$(out2))
	$(call test2.check,$(REFDIR)/$(out2x),$(out2x))

# check
$(TARGET).check:
	$(call test1.check,$(REFDIR)/$(out1),$(out1))
	$(call test2.check,$(REFDIR)/$(out2x),$(out2x))
	$(call plot.plbnds)
#	$(call plot.pldos)

# exec case1 for $(TARGET).
$(TARGET).case1.exec:
	@echo "nothing to exec for $@"

# check case1 for $(TARGET).
$(TARGET).case1.check:
	@echo "nothing to check for $@"

# plot graphs on the display
plot-band:
	ghostscript ps.plbnds.dat

plot-dos:
	ghostscript ps.pldos.dat


