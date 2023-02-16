include ../../Makefile

TARGET=cu

out1=out.lmf.$(TARGET)
out1x=bnds.$(TARGET)

###########################################################
#----- subroutine test1.exec
define test1.exec
	#-------------------------------------------
	# Case cu: illustration of high-lying local orbitals
	#          and bands of Cu up to ~50 eV.
	#-------------------------------------------
	#
	# --- Test 1.  Basic check of programs lmfa,lmf ---
	#
	@cat message.test1.$(TARGET)
	$(bindir)/lmfa --no-iactiv $(TARGET) -vnk=8 > $1 2>&1
	$(LMFP)  --no-iactiv $(TARGET) -vnk=8 -vbigbas=f >> $1 2>&1
	rm mixm.$(TARGET)
	$(LMFP)  --no-iactiv $(TARGET) -vnk=8 -vbigbas=t -vmetal=3 -vrsm2=1.3 -vrsmd1x=1 -vlmx=4 -vpwmode=0 \
		-voveps=0d-7 >> $1 2>&1
	$(LMFP)  --no-iactiv $(TARGET) -vnk=8 -vbigbas=t -vmetal=3 -vrsm2=1.3 -vrsmd1x=1 -vlmx=4 -vpwmode=0 \
		-voveps=0d-7 --band:fn=syml >> $1 2>&1
endef

###########################################################
# exec default case for $(TARGET).
$(TARGET).exec:
	$(call test1.exec,$(out1))
	$(call test1.check,$(REFDIR)/$(out1),$(out1))
	$(call test2.check,$(REFDIR)/$(out1x),$(out1x))

# check default case for $(TARGET).
$(TARGET).check:
	$(call test1.check,$(REFDIR)/$(out1),$(out1))
	$(call test2.check,$(REFDIR)/$(out1x),$(out1x))

# exec case1 for $(TARGET).
$(TARGET).case1.exec:
	@echo "nothing to exec for $@"

# check case1 for $(TARGET).
$(TARGET).case1.check:
	@echo "nothing to check for $@"
