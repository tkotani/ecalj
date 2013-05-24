include ../Makefile.define

TARGET=fe

out2a=out.lmf-dos.$(TARGET)
out2b=out.lmf-cls.$(TARGET)
out2c=out.lmf-mull.$(TARGET)

out2xa=tdos.$(TARGET)
out2xb=dos-cls.$(TARGET)
out2xc=dos-mull.$(TARGET)

all: help

help:
	@echo "Do not call this Makefile alone."

###########################################################
#----- subroutine test2.exec
define test2.exec
	#-------------------------------------------
	# Case Fe: spin-polarized Fe in bcc structure
	#-------------------------------------------
	#
	# --- Test 2.  Core-level spectroscopy (EELS), Mulliken analysis, partial DOS ---
	#
	@cat message.test2.$(TARGET)
	$(bindir)/lmfa $(TARGET) > /dev/null
	$(bindir)/lmf --rs=0 --cls:1,1,2 --mull:mode=2 -vnk=6 -vnit=1 $(TARGET) > $1 2>&1
	mv dos.$(TARGET) tdos.$(TARGET)
	$(bindir)/lmdos --cls --dos:wtfn=cls:npts=1001:window=-.7,.8 -vnk=6 $(TARGET) >> $2 2>&1
	mv dos.$(TARGET) dos-cls.$(TARGET)
	$(bindir)/lmdos --nosym --mull:mode=2 --dos:npts=1001:window=-.7,.8 -vnk=6 $(TARGET) >> $3 2>&1
	mv dos.$(TARGET) dos-mull.$(TARGET)
endef

###########################################################
#----- subroutine plot.pldos
define plot.pldos
	echo 20 10 -0.7 .8 | $(bindir)/pldos -fplot \
		-lst="1:7:2,19:31:2;9,11,15;13,17" -lst2 \
		dos-mull.$(TARGET) > plot.log 2>&1
	$(bindir)/fplot -disp -pr10 -f plot.dos >> plot.log 2>&1
	mv ps.dat ps.pldos.dat
	@echo "[note] run make plot-dos to plot dos"
endef

###########################################################
# exec default case for $(TARGET).
$(TARGET).exec:
	$(call test2.exec,$(out2a),$(out2b),$(out2c))
	$(call test2.check,$(REFDIR)/$(out2xa),$(out2xa))
	$(call test2.check,$(REFDIR)/$(out2xb),$(out2xb))
	$(call test2.check,$(REFDIR)/$(out2xc),$(out2xc))
#	$(call test2.check,$(REFDIR)/$(out2),$(out2))
#	$(call plot.pldos)

# check default case for $(TARGET).
$(TARGET).check:
	$(call test2.check,$(REFDIR)/$(out2xa),$(out2xa))
	$(call test2.check,$(REFDIR)/$(out2xb),$(out2xb))
	$(call test2.check,$(REFDIR)/$(out2xc),$(out2xc))

# exec case1 for $(TARGET).
$(TARGET).case1.exec:
	@echo "nothing to exec for $@"

# check case1 for $(TARGET).
$(TARGET).case1.check:
	@echo "nothing to check for $@"

# plot graph on the display
plot-dos:
	ghostscript ps.pldos.dat
