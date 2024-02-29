include ../../Makefile
TARGET=gas_epsPP_lmfh
MATERIAL=gas
$(TARGET).exec:
	mpirun -np 1 $(bindir)/lmfa $(MATERIAL) > llmfa
	mpirun -np $(mpi_size) $(bindir)/lmf $(MATERIAL) > llmf
	rm EPS*
	$(bindir)/epsPP_lmfh -np $(mpi_size) $(MATERIAL)
	@echo  >> ../summary.txt
	@-$(bindir)/diffnum EPS0001.nlfc.dat $(REFDIR)/EPS0001.nlfc.dat && echo "PASSED: $(TARGET) EPS0001.nlfc.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0001.nlfc.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0002.nlfc.dat $(REFDIR)/EPS0002.nlfc.dat && echo "PASSED: $(TARGET) EPS0002.nlfc.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0002.nlfc.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0003.nlfc.dat $(REFDIR)/EPS0003.nlfc.dat && echo "PASSED: $(TARGET) EPS0003.nlfc.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0003.nlfc.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0004.nlfc.dat $(REFDIR)/EPS0004.nlfc.dat && echo "PASSED: $(TARGET) EPS0004.nlfc.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0004.nlfc.dat" >> ../summary.txt

# check default case for $(TARGET).
$(TARGET).check:
	@echo  >> ../summary.txt
	@-$(bindir)/diffnum EPS0001.nlfc.dat $(REFDIR)/EPS0001.nlfc.dat && echo "PASSED: $(TARGET) EPS0001.nlfc.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0001.nlfc.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0002.nlfc.dat $(REFDIR)/EPS0002.nlfc.dat && echo "PASSED: $(TARGET) EPS0002.nlfc.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0002.nlfc.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0003.nlfc.dat $(REFDIR)/EPS0003.nlfc.dat && echo "PASSED: $(TARGET) EPS0003.nlfc.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0003.nlfc.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0004.nlfc.dat $(REFDIR)/EPS0004.nlfc.dat && echo "PASSED: $(TARGET) EPS0004.nlfc.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0004.nlfc.dat" >> ../summary.txt

