include ../Makefile.define

TARGET=gdn

out1=out.lmf.$(TARGET)
out2=out.lmf-dos.$(TARGET)
out2x=dos-mull.$(TARGET)

all: help

help:
	@echo "Do not call this Makefile alone."

###########################################################
#----- subroutine test1.exec
define test1.exec
	#-------------------------------------------
	# Case GdN: Test of LDA+U, and also LDA with spin polarized 4f core
	#-------------------------------------------
	#
	# --- Test 1.  Basic check of programs lmfa,lmf ---
	#
	rm -f atm.* eula.* moms.* mixm.* rst.* rsta.* save.* log.* hssn.* wkp.* cv.* bsmv.* bnds.* dmats.* dmats-save.* sigm.* qpp.* 
	@cat message.test1.$(TARGET)
	$(bindir)/lmfa $(TARGET) > $1 2>&1
	$(LMFP)  $(TARGET) >> $1 2>&1
	mv rst.$(TARGET) rst-save.$(TARGET)
	mv dmats.$(TARGET) dmats-save.$(TARGET)
	rm -f rst.$(TARGET) wkp.$(TARGET) mixm.$(TARGET)
	cp rst-save.$(TARGET) rst.$(TARGET)
	cp dmats-save.$(TARGET) dmats.$(TARGET)
	$(LMFP)  $(TARGET) --rs=1,0 -vpwmode=1 >> $1 2>&1
	rm -f rst.$(TARGET) wkp.$(TARGET) mixm.$(TARGET)
	$(bindir)/lmfa $(TARGET) -vcoref=t >> $1 2>&1
	$(LMFP)  $(TARGET) -vcoref=t >> $1 2>&1
endef

#----- subroutine test2.exec
define test2.exec
	#
	# --- Test 2.  Core-level spectroscopy (EELS), Mulliken analysis, partial DOS ---
	#
	@cat message.test2.$(TARGET)
	cp rst-save.$(TARGET) rst.$(TARGET)
	cp dmats-save.$(TARGET) dmats.$(TARGET)
	$(LMF) --nosym $(TARGET) --mull:mode=1 -vnk=2 -vnit=1 -vso=t > $1 2>&1
	mv dos.$(TARGET) tdos.$(TARGET)
	$(bindir)/lmdos --nosym -vso=t --mull:mode=1 \
		--dos:fpdos:npts=1001:window=-1,1 \
		-vnk=2 $(TARGET) >> $1 2>&1
	mv dos.$(TARGET) dos-mull.$(TARGET)
endef

#----- subroutine plot.pldos
define plot.pldos
	echo 150 10 -9 10 | $(bindir)/pldos -ef=0 -escl=13.6 -fplot '-lst=8' \
		-ref:fn=tdos.$(TARGET):chan=1:scale dos-mull.$(TARGET) \
		> plot.log 2>&1
	$(bindir)/fplot -disp -pr10 -f plot.dos >> plot.log 2>&1
	mv ps.dat ps.pldos.dat
	@echo "[note] run make plot-dos to plot dos"
endef

###########################################################
# exec default case for $(TARGET).
$(TARGET).exec:
	$(call test1.exec,$(out1))
	$(call test1.check,$(REFDIR)/$(out1),$(out1))
	$(call test2.exec,$(out2))
	$(call test2.check,$(REFDIR)/$(out2x),$(out2x))

# check only
$(TARGET).check:
	$(call test1.check,$(REFDIR)/$(out1),$(out1))
	$(call test2.check,$(REFDIR)/$(out2x),$(out2x))
	$(call plot.pldos)
	ghostscript ps.pldos.dat

# exec case1 for $(TARGET).
$(TARGET).case1.exec:
	@echo "nothing to exec for $@"

# check case1 for $(TARGET).
$(TARGET).case1.check:
	@echo "nothing to check for $@"

# plot graph on the display
plot-dos:
	ghostscript ps.pldos.dat

