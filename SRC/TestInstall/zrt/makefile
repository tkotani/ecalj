include ../../Makefile

TARGET=zrt

out1=out.lmf.$(TARGET)

###########################################################
#----- subroutine test1.exec
define test1.exec
	#-------------------------------------------
	# Case zrt: ZrO_2 fluorite in tetragonal setting, with tetragonal distortion
	#-------------------------------------------
	#
	# --- Test 1.  Basic check of programs lmfa,lmf ---
	#
	@cat message.test1.$(TARGET)
	$(bindir)/lmfa -vfp=1 $(TARGET) --no-iactiv > $1 2>&1
	$(LMFP)  -vnitq=1 -vforce=1 -vfp=1 $(TARGET) >> $1 2>&1
	$(LMFP)  -vnitq=1 -vforce=1 -vfp=1 $(TARGET) >> $1 2>&1
endef

###########################################################
# exec default case for $(TARGET).
$(TARGET).exec:
	$(call test1.exec,$(out1))
	$(call test1.check,$(REFDIR)/$(out1),$(out1))

# check default case for $(TARGET).
$(TARGET).check:
	$(call test1.check,$(REFDIR)/$(out1),$(out1))

# exec case1 for $(TARGET).
$(TARGET).case1.exec:
	@echo "nothing to exec for $@"

# check case1 for $(TARGET).
$(TARGET).case1.check:
	@echo "nothing to check for $@"
