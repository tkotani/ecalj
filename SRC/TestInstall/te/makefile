include ../../Makefile
TARGET=te
out1=out.lmf.$(TARGET)
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
	mpirun -np 1 $(bindir)/lmfa $(TARGET) -vdyn='DYN' -vnk=3 -vnit=3 -vlf1=4 \
		-vlmxl=4 -vnk=3 -vngd=20 -vkmx=3 -vconv=1e-4 > $1 2>&1
        $(LMFP) $(TARGET) -vdyn='DYN' -vnk=3 -vnit=3 -vlf1=4 \
		-vlmxl=4 -vnk=3 -vngd=20 -vkmx=3 -vconv=1e-4 -vnbas=12 -vnspec=2 >> $1 2>&1
        rm -f mixm.$(TARGET)
	cp rst.$(TARGET) rst.$(TARGET).bk
        $(LMFP) $(TARGET) -vnk=3 -vnit=3 -vlf1=4 \
		-vlmxl=4 -vnk=3 -vngd=20 -vkmx=5 -vconv=1e-4 -vpwmode=11 >> $1 2>&1
        rm -f mixm.$(TARGET)
	cp rst.$(TARGET).bk rst.$(TARGET)
        $(LMFP)  $(TARGET) -vnk=3 -vnit=3 -vlf1=4 \
		-vlmxl=4 -vnk=3 -vngd=20 -vkmx=3 -vconv=1e-4 >> $1 2>&1
endef
###########################################################
# exec default case for $(TARGET).
$(TARGET).exec:
	$(call test1.exec,$(out1))
	$(call test1.check,$(REFDIR)/$(out1),$(out1))

# check default case for $(TARGET).
$(TARGET).check:
	$(call test1.check,$(REFDIR)/$(out1),$(out1))

