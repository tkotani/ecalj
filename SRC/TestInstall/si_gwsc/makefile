#include ../Makefile.define

# directroy where input/output files are saved.
TARGET=si_gwsc
# material in the tests.
MATERIAL=si


# exec default case for $(TARGET).
$(TARGET).exec:
	# exec lmfa command.
	mpirun -np 1 $(bindir)/lmfa $(MATERIAL) > llmfa
	# exec gwsc command.
	rm log.$(MATERIAL) QPU
	$(bindir)/gwsc 0 -np $(mpi_size) $(MATERIAL)
	@echo  >> ../summary.txt
	@-$(bindir)/dqpu $(REFDIR)/QPU QPU && echo "PASSED: $(TARGET) QPU" >> ../summary.txt || echo "FAILED: $(TARGET) QPU compared by ./bin/dqpu" >> ../summary.txt
	@-$(bindir)/diffnum $(REFDIR)/log.$(MATERIAL) log.$(MATERIAL) 'fp evl' && echo "PASSED: $(TARGET) diffnum log.$(MATERIAL) fp evl" >> ../summary.txt || echo "FAILED: $(TARGET) diffnum log.$(MATERIAL) fp evl" >> ../summary.txt

# check default case for $(TARGET).
$(TARGET).check:
	# add summary.txt
	@echo  >> ../summary.txt
	@-$(bindir)/dqpu $(REFDIR)/QPU QPU && echo "PASSED: $(TARGET) QPU" >> ../summary.txt || echo "FAILED: $(TARGET) QPU compared by ./bin/dqpu" >> ../summary.txt
	@-$(bindir)/diffnum $(REFDIR)/log.$(MATERIAL) log.$(MATERIAL) 'fp evl' && echo "PASSED: $(TARGET) diffnum log.$(MATERIAL) fp evl" >> ../summary.txt || echo "FAILED: $(TARGET) diffnum log.$(MATERIAL) fp evl" >> ../summary.txt

