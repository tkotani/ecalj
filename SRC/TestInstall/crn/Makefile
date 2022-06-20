include ../Makefile.define

TARGET=crn

out2=out.lmf-dos.$(TARGET)
out2x=dos-vcdmel.$(TARGET)

all: help

help:
	@echo "Do not call this Makefile alone."

###########################################################
#----- subroutine test2.exec
define test2.exec
	#-------------------------------------------
	# Case CrN: test of CLS with core hole
	#-------------------------------------------
	#
	# --- Test 2.  Core-level spectroscopy (EELS), Mulliken analysis, partial DOS ---
	#
	@cat message.test2.$(TARGET)
	$(bindir)/lmfa $(TARGET) > $1 2>&1
	$(LMFP) $(TARGET) >> $1 2>&1
	#echo 'CLSinput is used for EELS mode (--cls). Works well. But not examined(2022-6-20)'
	$(LMFP) -vnit=1 -vmetal=2 $(TARGET) >> $1 2>&1
	$(LMFP) --cls $(TARGET) >> $1 2>&1
endef


###########################################################
# exec default case for $(TARGET).
$(TARGET).exec:
	$(call test2.exec,$(out2))
	$(call test2.check,$(REFDIR)/$(out2x),$(out2x))
#	$(call plot.pldos)

# check default case for $(TARGET).
$(TARGET).check:
	$(call test2.check,$(REFDIR)/$(out2x),$(out2x))

# # exec case1 for $(TARGET).
# $(TARGET).case1.exec:
# 	@echo "nothing to exec for $@"

# # check case1 for $(TARGET).
# $(TARGET).case1.check:
# 	@echo "nothing to check for $@"

# # plot graph on the display
# plot-dos:
# 	ghostscript ps.pldos.dat
