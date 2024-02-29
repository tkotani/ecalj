include ../../Makefile

TARGET=fe

out2a=out.lmf-dos.$(TARGET)
out2b=out.lmf-cls.$(TARGET)
out2c=out.lmf-mull.$(TARGET)

out2xa=tdos.$(TARGET)
out2xb=dos-cls.$(TARGET)
out2xc=dos-mull.$(TARGET)

outtdos=dos.tot.$(TARGET)

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
	mpirun -np 1 $(bindir)/lmfa $(TARGET) > $1 2>&1
	$(LMF) --cls:1,1,2 --mull:mode=2 -vnk=6 -vnit=1 $(TARGET) >> $1 2>&1
	mv dos.$(TARGET) tdos.$(TARGET)
	$(bindir)/lmdos --cls --dos:wtfn=cls:npts=1001:window=-.7,.8 -vnk=6 $(TARGET) >> $2 2>&1
	mv dos.$(TARGET) dos-cls.$(TARGET)
	$(bindir)/lmdos --nosym --mull:mode=2 --dos:npts=1001:window=-.7,.8 -vnk=6 $(TARGET) >> $3 2>&1
	mv dos.$(TARGET) dos-mull.$(TARGET)

endef

#----- subroutine job_tdos
define testjob1.exec
	#-------------------------------------------
	# Case Fe: spin-polarized Fe in bcc structure
	#-------------------------------------------
	#
	# --- Test 3.  job_tdos ---------
	#
	mpirun -np 1 $(bindir)/lmfa $(TARGET) > /dev/null
	$(bindir)/job_tdos $(TARGET) -np $(mpi_size) -vnk=6 NoGnuplot > ljob_tdos
endef

###########################################################
# exec default case for $(TARGET).
$(TARGET).exec:
	$(call test2.exec,$(out2a),$(out2b),$(out2c))
	$(call test2.check,$(REFDIR)/$(out2xa),$(out2xa))
	$(call test2.check,$(REFDIR)/$(out2xb),$(out2xb))
	$(call test2.check,$(REFDIR)/$(out2xc),$(out2xc))
	$(call testjob1.exec,$(outtdos))
	$(call test2.check,$(REFDIR)/$(outtdos),$(outtdos))

# check default case for $(TARGET).
$(TARGET).check:
	$(call test2.check,$(REFDIR)/$(out2xa),$(out2xa))
	$(call test2.check,$(REFDIR)/$(out2xb),$(out2xb))
	$(call test2.check,$(REFDIR)/$(out2xc),$(out2xc))
	$(call test2.check,$(REFDIR)/$(outtdos),$(outtdos))

# exec case1 for $(TARGET).
$(TARGET).case1.exec:
	@echo "nothing to exec for $@"

# check case1 for $(TARGET).
$(TARGET).case1.check:
	@echo "nothing to check for $@"

# plot graph on the display
plot-dos:
	ghostscript ps.pldos.dat
