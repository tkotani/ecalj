include ../../Makefile

# directroy where input/output files are saved.
TARGET=fe_epsPP_lmfh_chipm
# material in the tests.
MATERIAL=fe

# exec default case for $(TARGET).
$(TARGET).exec:
	# exec lmfa command.
	$(bindir)/lmfa $(MATERIAL) > llmfa
	# exec lmf command.
	mpirun -np 4 $(bindir)/lmf $(MATERIAL) > llmf
	# exec epsPP_lmfh_chipm command.
	rm ChiPM*
	$(bindir)/epsPP_lmfh_chipm -np $(mpi_size) $(MATERIAL)
	@echo  >> ../summary.txt
	@-$(bindir)/diffnum ChiPM0001.nlfc.mat $(REFDIR)/ChiPM0001.nlfc.mat && echo "PASSED: $(TARGET) ChiPM0001.nlfc.mat" >> ../summary.txt || echo "FAILED: $(TARGET) ChiPM0001.nlfc.mat" >> ../summary.txt
	@-$(bindir)/diffnum ChiPM0002.nlfc.mat $(REFDIR)/ChiPM0002.nlfc.mat && echo "PASSED: $(TARGET) ChiPM0002.nlfc.mat" >> ../summary.txt || echo "FAILED: $(TARGET) ChiPM0002.nlfc.mat" >> ../summary.txt
	@-$(bindir)/diffnum ChiPM0003.nlfc.mat $(REFDIR)/ChiPM0003.nlfc.mat && echo "PASSED: $(TARGET) ChiPM0003.nlfc.mat" >> ../summary.txt || echo "FAILED: $(TARGET) ChiPM0003.nlfc.mat" >> ../summary.txt
	@-$(bindir)/diffnum ChiPM0004.nlfc.mat $(REFDIR)/ChiPM0004.nlfc.mat && echo "PASSED: $(TARGET) ChiPM0004.nlfc.mat" >> ../summary.txt || echo "FAILED: $(TARGET) ChiPM0004.nlfc.mat" >> ../summary.txt
	@-$(bindir)/diffnum ChiPM0005.nlfc.mat $(REFDIR)/ChiPM0005.nlfc.mat && echo "PASSED: $(TARGET) ChiPM0005.nlfc.mat" >> ../summary.txt || echo "FAILED: $(TARGET) ChiPM0005.nlfc.mat" >> ../summary.txt

# check default case for $(TARGET).
$(TARGET).check:
	# add summary.txt
	@echo  >> ../summary.txt
	@-$(bindir)/diffnum ChiPM0001.nlfc.mat $(REFDIR)/ChiPM0001.nlfc.mat && echo "PASSED: $(TARGET) ChiPM0001.nlfc.mat" >> ../summary.txt || echo "FAILED: $(TARGET) ChiPM0001.nlfc.mat" >> ../summary.txt
	@-$(bindir)/diffnum ChiPM0002.nlfc.mat $(REFDIR)/ChiPM0002.nlfc.mat && echo "PASSED: $(TARGET) ChiPM0002.nlfc.mat" >> ../summary.txt || echo "FAILED: $(TARGET) ChiPM0002.nlfc.mat" >> ../summary.txt
	@-$(bindir)/diffnum ChiPM0003.nlfc.mat $(REFDIR)/ChiPM0003.nlfc.mat && echo "PASSED: $(TARGET) ChiPM0003.nlfc.mat" >> ../summary.txt || echo "FAILED: $(TARGET) ChiPM0003.nlfc.mat" >> ../summary.txt
	@-$(bindir)/diffnum ChiPM0004.nlfc.mat $(REFDIR)/ChiPM0004.nlfc.mat && echo "PASSED: $(TARGET) ChiPM0004.nlfc.mat" >> ../summary.txt || echo "FAILED: $(TARGET) ChiPM0004.nlfc.mat" >> ../summary.txt
	@-$(bindir)/diffnum ChiPM0005.nlfc.mat $(REFDIR)/ChiPM0005.nlfc.mat && echo "PASSED: $(TARGET) ChiPM0005.nlfc.mat" >> ../summary.txt || echo "FAILED: $(TARGET) ChiPM0005.nlfc.mat" >> ../summary.txt

# exec case1 for $(TARGET).
$(TARGET).case1.exec:
	@echo "nothing to exec for $@"

# check case1 for $(TARGET).
$(TARGET).case1.check:
	@echo "nothing to check for $@"
