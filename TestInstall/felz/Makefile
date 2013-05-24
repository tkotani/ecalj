include ../Makefile.define

TARGET=felz

out3=out.lmf.fsmom.$(TARGET)
out4=out.lmf.lzsz.$(TARGET)
out4b=out.lmf.ls.$(TARGET)

all: help

help:
	@echo "Do not call this Makefile alone."

###########################################################
#----- subroutine test3.exec
define test3.exec
	#-------------------------------------------
	# Case felz: spin-polarized Fe spin-orbit coupling
	#-------------------------------------------
	#
	# --- Test 3.  Check of miscellaneous special features, programs lmfa,lmf ---
	#
	@cat message.test3.$(TARGET)
	$(bindir)/lmfa -vrel=1 -vnit=1 -vso=0 $(TARGET) > $1 2>&1
	$(bindir)/lmf -vrel=1 -vnit=3 -vso=2 --rs=2 $(TARGET) \
		-vfsmom=-2 >> $1 2>&1
endef

###########################################################
#----- subroutine test4a.exec
define test4a.exec
	#
	# --- Test case 4:  Spin-orbit coupling ---
	#
	@cat message.test4.$(TARGET)
	$(bindir)/lmfa -vrel=1 -vnit=1 -vso=0 $(TARGET) > $1 2>&1
	$(bindir)/lmf -vrel=1 -vnit=1 -vso=2 --rs=2 $(TARGET) -vpwmode=11 >> $1 2>1
	rm -f mixm.$(TARGET)
	$(bindir)/lmf -vrel=1 -vnit=1 -vso=2 --rs=2 $(TARGET) >> $1 2>1
endef

###########################################################
#----- subroutine test4b.exec
define test4b.exec
	$(bindir)/lmf -vrel=1 -vnit=1 -vso=1 --rs=2 $(TARGET) > $1 2>&1
endef

###########################################################
# exec default case for $(TARGET).
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

# exec case1 for $(TARGET).
$(TARGET).case1.exec:
	@echo "nothing to exec for $@"

# check case1 for $(TARGET).
$(TARGET).case1.check:
	@echo "nothing to check for $@"
