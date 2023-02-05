include ../Makefile.define

TARGET=te

out1=out.lmf.$(TARGET)

all: help

help:
	@echo "Do not call this Makefile alone."

###########################################################
#----- subroutine test1.exec
define test1.exec
	#-------------------------------------------
	# Case te: molecular statics in an open structure
	#-------------------------------------------
	#
	# --- Test 1.  Basic check of programs lmfa,lmf ---
	#
	@cat message.test1.$(TARGET)
	$(bindir)/lmfa $(TARGET) -vminx=t -vnk=3 -vnit=3 -vlf1=4 \
		-vlmxl=4 -vnk=3 -vngd=20 -vkmx=3 -vconv=1d-4 > $1 2>&1
        $(LMFP) $(TARGET) -vminx=t -vnk=3 -vnit=3 -vlf1=4 \
		-vlmxl=4 -vnk=3 -vngd=20 -vkmx=3 -vconv=1d-4 -vnbas=12 -vnspec=2 >> $1 2>&1
        rm -f mixm.$(TARGET)
	cp rst.$(TARGET) rst.$(TARGET).bk
        $(LMFP) $(TARGET) -vminx=0 -vnk=3 -vnit=3 -vlf1=4 \
		-vlmxl=4 -vnk=3 -vngd=20 -vkmx=5 -vconv=1d-4 -vpwmode=11 >> $1 2>&1
        rm -f mixm.$(TARGET)
	cp rst.$(TARGET).bk rst.$(TARGET)
        $(LMFP)  $(TARGET) -vminx=0 -vnk=3 -vnit=3 -vlf1=4 \
		-vlmxl=4 -vnk=3 -vngd=20 -vkmx=3 -vconv=1d-4 >> $1 2>&1
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
