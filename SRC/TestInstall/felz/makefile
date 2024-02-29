include ../../Makefile
TARGET=felz
out3=out.lmf.fsmom.$(TARGET)
out4=out.lmf.lzsz.$(TARGET)
out4b=out.lmf.ls.$(TARGET)
###########################################################
define test3.exec
	#-------------------------------------------
	# Case felz: spin-polarized Fe spin-orbit coupling
	#-------------------------------------------
	#
	# --- Test 3.  Check of miscellaneous special features, programs lmfa,lmf ---
	#
	@cat message.test3.$(TARGET)
	mpirun -np 1 $(bindir)/lmfa -vrel=1 -vso=0 $(TARGET) > $1 2>&1
	$(LMFP) -vrel=1 -vnit=3 -vso=2 $(TARGET) \
		-vfsmom=-2 >> $1 2>&1
endef
define test4a.exec
	#
	# --- Test case 4:  Spin-orbit coupling ---
	#
	@cat message.test4.$(TARGET)
	mpirun -np 1 $(bindir)/lmfa -vrel=1 -vnit=1 -vso=0 $(TARGET) > $1 2>&1
	$(LMFP) -vrel=1 -vnit=1 -vso=2 $(TARGET) -vpwmode=11 >> $1 2>1
	rm -f mixm.$(TARGET)
	$(LMFP) -vrel=1 -vnit=1 -vso=2 $(TARGET) >> $1 2>1
endef
define test4b.exec
	$(LMFP) -vrel=1 -vnit=1 -vso=1 $(TARGET) > $1 2>&1
endef
###########################################################
$(TARGET).exec:
	$(call test3.exec,$(out3))
	$(call test3.check,$(REFDIR)/$(out3),$(out3))	
	rm -f atm.* fs.* moms.* mixm.* rst.* save.* log.* hssn.* wkp.* bsmv.* syml.* bnds.*
	$(call test4a.exec,$(out4))
	$(call test4a.check,$(REFDIR)/$(out4),$(out4))
	$(call test4b.exec,$(out4b))
	$(call test4b.check,$(REFDIR)/$(out4b),$(out4b))

# check default case for $(TARGET).
$(TARGET).check:
	$(call test3.check,$(REFDIR)/$(out3),$(out3))	
	$(call test4a.check,$(REFDIR)/$(out4),$(out4))
	$(call test4b.check,$(REFDIR)/$(out4b),$(out4b))
