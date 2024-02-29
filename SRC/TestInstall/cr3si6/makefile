include ../../Makefile

TARGET=cr3si6

out1=out.lmf.$(TARGET)

###########################################################
#----- subroutine test1.exec
define test1.exec
	#-------------------------------------------
	# Case cr3si6: a hexagonal environment with several atoms, two classes
        #-------------------------------------------
	#
	# --- Test 1.  Basic check of programs lmfa,lmf ---
	#
	@cat message.test1.$(TARGET)
	mpirun -np 1 $(bindir)/lmfa $(TARGET) --pr51 -vnit=2 --time=6 --no-iactiv > $1 2>&1
	$(LMFP)  $(TARGET) --pr51 -vnit=2 --time=6 --no-iactiv >> $1 2>&1
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
