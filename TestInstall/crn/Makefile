include ../Makefile.define

TARGET=crn

out2=out.lmf-dos.$(TARGET)
out2x=dos-cls.$(TARGET)

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
	$(LMF) --rs=1,0 --cls:5,0,1 -vnit=1 -vmetal=2 $(TARGET) >> $1 2>&1
	$(bindir)/lmdos --dos:cls:window=0,1:npts=101 --cls $(TARGET) >> $1 2>&1
	mv dos.$(TARGET) dos-cls.$(TARGET)
endef

###########################################################
#----- subroutine plot.pldos
define plot.pldos
	echo .25 10 0 1 | $(bindir)/pldos -fplot -lst="1;3;5" -lst2="2;4;6" \
		dos-cls.$(TARGET) > plot.log 2>&1
	$(bindir)/fplot -disp -pr10 -f plot.dos >> plot.log 2>&1
	mv ps.dat ps.pldos.dat
	@echo "[note] run make plot-dos to plot dos"
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

# exec case1 for $(TARGET).
$(TARGET).case1.exec:
	@echo "nothing to exec for $@"

# check case1 for $(TARGET).
$(TARGET).case1.check:
	@echo "nothing to check for $@"

# plot graph on the display
plot-dos:
	ghostscript ps.pldos.dat
