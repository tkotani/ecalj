include ../../Makefile
TARGET=gasls
out4c=out.lmf.ls-bands.$(TARGET)
# position of evals of interest
lineeval=3
evalso=5
define test4c.exec
	#-------------------------------------------
	# Case GaAs: GaAs with spin-orbit coupling
	#-------------------------------------------
	#
	# --- Test case 4:  Spin-orbit coupling ---
	#
	@cat message.test4.$(TARGET)
	$(bindir)/lmfa -vso=1 $(TARGET) -vpwmode=0 > $1 2>&1
	$(LMFP) -vso=1 -vnit=1 $(TARGET) --band:fn=syml -vpwmode=0 >> $1 2>&1
endef
###########################################################
$(TARGET).exec:
	$(call test4c.exec,$(out4c))
	$(call test4c.check,$(REFDIR)/$(out4c),$(out4c))
# check default case for $(TARGET).
$(TARGET).check:
	$(call test4c.check,$(REFDIR)/$(out4c),$(out4c))
