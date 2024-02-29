include ../../Makefile
TARGET=na
out1=out.lmf.$(TARGET)
define test1.exec
	#-------------------------------------------
	# Case na: illustration of low- and high-lying local orbitals
	#-------------------------------------------
	#
	# --- Test 1.  Basic check of programs lmfa,lmf ---
	#
	@cat message.test1.$(TARGET)
	mpirun -np 1 $(bindir)/lmfa --no-iactiv $(TARGET) -vnk=6 -vnapval=0 -vpz1=0 -vp1=3.38  > $1 2>&1
	$(LMFP)  --no-iactiv $(TARGET)       -vnk=6 -vnapval=0 -vpz1=0 -vp1=3.38 >> $1 2>&1
	rm rst.$(TARGET) mixm.$(TARGET)
	mpirun -np 1 $(bindir)/lmfa --no-iactiv $(TARGET) -vnk=6 -vnapval=1 >> $1 2>&1
	$(LMFP)  --no-iactiv $(TARGET) -vnk=6       -vnapval=1 >> $1 2>&1
	rm mixm.$(TARGET) rst.$(TARGET)
	$(LMFP)  --no-iactiv $(TARGET) -vnk=6 -vnapval=2 -vpz1=12.94 >> $1 2>&1
endef
###########################################################
$(TARGET).exec:
	$(call test1.exec,$(out1))
	$(call test1.check,$(REFDIR)/$(out1),$(out1))
# check default case for $(TARGET).
$(TARGET).check:
	$(call test1.check,$(REFDIR)/$(out1),$(out1))
