include ../../Makefile

TARGET=copt

out1=out.lmf.$(TARGET)

###########################################################
#----- subroutine test1.exec
define test1.exec
	#-------------------------------------------
	# Case copt: a distorted L12 environment with four atoms.
	#-------------------------------------------
	#
	# --- Test 1.  Basic check of programs lmfa,lmf ---
	#
	@cat message.test1.$(TARGET)
	mpirun -np 1 $(bindir)/lmfa $(TARGET) -vnsp=2 -cstrx=l12 --pr41,41 -vmet=3 \
		-vtetra=0 -vnk1=2 -vlfrce=12 -vdist=1 -vnit=3 \
		--time=5 --no-iactiv > $1 2>&1
	$(LMFP)  $(TARGET) -vnsp=2 -cstrx=l12 --pr41,41 -vmet=3 \
		-vtetra=0 -vnk1=2 -vlfrce=12 -vdist=1 -vnit=3 \
		--time=5 --no-iactiv --cvK:158,790 -vpwmode=0 >> $1 2>&1
endef

###########################################################
# exec default case for $(TARGET).
$(TARGET).exec:
	$(call test1.exec,$(out1))
	$(call test1.check,$(REFDIR)/$(out1),$(out1))

# check case1 for $(TARGET).
$(TARGET).check:
	$(call test1.check,$(REFDIR)/$(out1),$(out1))

# exec case1 for $(TARGET).
$(TARGET).case1.exec:
	@echo "nothing to exec for $@"

# check case1 for $(TARGET).
$(TARGET).case1.check:
	@echo "nothing to check for $@"
