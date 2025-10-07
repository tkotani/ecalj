include ../Makefile.define

TARGET=er

out1=out.lmf.$(TARGET)

all: help

help:
	@echo "Do not call this Makefile alone."

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
	$(bindir)/lmfa $(TARGET) -vso=0 -vldau=1 -vnsp=2 -vtpp1=3 -vtpf1=4 \
		-vbigbas=3 -vehmax=-.4 -vnk=4 -vnk2=3 -vconvc=1e-5 \
		-vbeta=.3 -vmet=2 -vnl=7c > $1 2>&1
	$(bindir)/lmf  $(TARGET) -vso=0 -vldau=1 -vnsp=2 -vtpp1=3 -vtpf1=4 \
		-vbigbas=3 -vehmax=-.4 -vnk=4 -vnk2=3 -vconvc=1e-5 \
		-vbeta=.3 -vmet=2 -vnl=7 -vnit=20 >> $1 2>&1
	$(bindir)/lmf  $(TARGET) -vso=0 -vldau=1 -vnsp=2 -vtpp1=3 -vtpf1=4 \
		-vbigbas=3 -vehmax=-.4 -vnk=4 -vnk2=3 -vconvc=1e-5 \
		-vbeta=.3 -vmet=2 -vnl=7 -vnit=20 \
		--band~col=10:16,26:32,36:42~fn=syml >> $1 2>&1
	rm -f mixm.$(TARGET)
	$(bindir)/lmf  $(TARGET) -vso=1 -vldau=1 -vnsp=2 -vtpp1=3 -vtpf1=4 \
		-vbigbas=3 -vehmax=-.4 -vnk=4 -vnk2=3 -vconvc=1e-5 \
		-vbeta=.3 -vmet=2 -vnl=7 -vnit=1 >> $1 2>&1
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
