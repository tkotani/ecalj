include ../Makefile.define

TARGET=tio2

out1=out.lmf.$(TARGET)

all: help

help:
	@echo "Do not call this Makefile alone."

###########################################################
#----- subroutine test1.exec
define test1.exec
	#-------------------------------------------
	# Case tio2: example of relaxation
	#-------------------------------------------
	#
	# --- Test 1.  Basic check of programs lmfa,lmf ---
	#
	@cat message.test1.$(TARGET)
	$(bindir)/lmfa $(TARGET) -vnit=4 $(TARGET) -vlmf=t -vnit=100 \
		-vnk=4 -vnc=0 -vsforce=1 -vsdmod=1 -vbzj=1 -vfile=t \
		-vrelax=t -cimp=Co -vfsmom=0 --rs=0 > $1 2>&1
	$(bindir)/lmf -vnit=4 $(TARGET) -vlmf=t -vnit=100 -vnk=4 \
		-vnc=0 -vsforce=1 -vsdmod=1 -vbzj=1 -vfile=t -vrelax=t \
		-cimp=Co -vfsmom=0 --rs=0 >> $1 2>&1
endef

###########################################################
# commands to execute case1 for $(TARGET).
$(TARGET).exec:
	$(call test1.exec,$(out1))
	$(call test1.check,$(REFDIR)/$(out1),$(out1))

# commands to check result of case1 for $(TARGET).
$(TARGET).check:
	$(call test1.check,$(REFDIR)/$(out1),$(out1))

# commands to execute case1 for $(TARGET).
$(TARGET).case1.exec:
	@echo "nothing to exec for case: $@"

# commands to check result of case1 for $(TARGET).
$(TARGET).case1.check:
	@echo "nothing to check for case: $@"

# commands to execute case2 for $(TARGET).
$(TARGET).case2.exec:
	@echo "nothing to exec for case: $@"

# commands to check result of case2 for $(TARGET).
$(TARGET).case2.check:
	@echo "nothing to check for case: $@"

###########################################################
# Do not edit lines below.

ifneq ($(CHECKONLY),yes)

$(TARGET): $(TARGET).exec
$(TARGET).case1: $(TARGET).case1.exec
$(TARGET).case2: $(TARGET).case2.exec

else

$(TARGET): $(TARGET).check
$(TARGET).case1: $(TARGET).case1.check
$(TARGET).case2: $(TARGET).case2.check

endif

