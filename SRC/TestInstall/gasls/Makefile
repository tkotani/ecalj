include ../Makefile.define

TARGET=gasls

out4c=out.lmf.ls-bands.$(TARGET)

# position of evals of interest
lineeval=3
evalso=5

all: help

help:
	@echo "Do not call this Makefile alone."

###########################################################
#----- subroutine test4.exec
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
# exec default case for $(TARGET).
$(TARGET).exec:
	$(call test4c.exec,$(out4c))
	$(call test4c.check,$(REFDIR)/$(out4c),$(out4c))

# check default case for $(TARGET).
$(TARGET).check:
	$(call test4c.check,$(REFDIR)/$(out4c),$(out4c))

# exec case1 for $(TARGET).
$(TARGET).case1.exec:
	@echo "nothing to exec for $@"

# check case1 for $(TARGET).
$(TARGET).case1.check:
	@echo "nothing to check for $@"
