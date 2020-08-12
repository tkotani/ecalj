include ../Makefile.define

TARGET=c

out3a=out.lmf.neutral.$(TARGET)
out3b=out.lmf.ionized.$(TARGET)

all: help

help:
	@echo "Do not call this Makefile alone."

###########################################################
#----- subroutine test3.exec
define test3.exec
	#-------------------------------------------
	# Case C: test of homogeneous background
	#-------------------------------------------
	#
	# --- Test 3.  Check of miscellaneous special features, programs lmfa,lmf ---
	#
	@cat message.test3.$(TARGET)
	$(bindir)/lmfa --no-iactiv $(TARGET) -vzbak=$2 >  $1 2>&1
	$(LMFP)  --no-iactiv $(TARGET) -vzbak=$2 >> $1 2>&1
endef

###########################################################
# exec default case for $(TARGET).
$(TARGET).exec:
	$(call test3.exec,$(out3a),0)
	$(call test3.check,$(REFDIR)/$(out3a),$(out3a))
	rm -f mixm.* rst.* save.* log.* hssn.* wkp.* bsmv.* bnds.*
	$(call test3.exec,$(out3b),1)
	$(call test3.check,$(REFDIR)/$(out3b),$(out3b))

# check default case for $(TARGET).
$(TARGET).check:
	$(call test3.check,$(REFDIR)/$(out3a),$(out3a))
	$(call test3.check,$(REFDIR)/$(out3b),$(out3b))

# exec case1 for $(TARGET).
$(TARGET).case1.exec:
	@echo "nothing to exec for $@"

# check case1 for $(TARGET).
$(TARGET).case1.check:
	@echo "nothing to check for $@"
