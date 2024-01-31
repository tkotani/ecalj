TARGET=gas_eps_lmfh
MATERIAL=gas
$(TARGET).exec:
	$(bindir)/lmfa $(MATERIAL) > llmfa
	$(bindir)/lmf $(MATERIAL) > llmf
	rm EPS*
	$(bindir)/eps_lmfh -np $(mpi_size) $(MATERIAL)
	@echo  >> ../summary.txt
	@-$(bindir)/diffnum EPS0001.nlfc.dat $(REFDIR)/EPS0001.nlfc.dat && echo "PASSED: $(TARGET) EPS0001.nlfc.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0001.nlfc.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0002.nlfc.dat $(REFDIR)/EPS0002.nlfc.dat && echo "PASSED: $(TARGET) EPS0002.nlfc.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0002.nlfc.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0003.nlfc.dat $(REFDIR)/EPS0003.nlfc.dat && echo "PASSED: $(TARGET) EPS0003.nlfc.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0003.nlfc.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0004.nlfc.dat $(REFDIR)/EPS0004.nlfc.dat && echo "PASSED: $(TARGET) EPS0004.nlfc.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0004.nlfc.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0001.dat $(REFDIR)/EPS0001.dat && echo "PASSED: $(TARGET) EPS0001.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0001.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0002.dat $(REFDIR)/EPS0002.dat && echo "PASSED: $(TARGET) EPS0002.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0002.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0003.dat $(REFDIR)/EPS0003.dat && echo "PASSED: $(TARGET) EPS0003.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0003.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0004.dat $(REFDIR)/EPS0004.dat && echo "PASSED: $(TARGET) EPS0004.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0004.dat" >> ../summary.txt

$(TARGET).check:
	# add summary.txt
	# check a result by diffnum command.
	@echo  >> ../summary.txt
	@-$(bindir)/diffnum EPS0001.nlfc.dat $(REFDIR)/EPS0001.nlfc.dat && echo "PASSED: $(TARGET) EPS0001.nlfc.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0001.nlfc.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0002.nlfc.dat $(REFDIR)/EPS0002.nlfc.dat && echo "PASSED: $(TARGET) EPS0002.nlfc.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0002.nlfc.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0003.nlfc.dat $(REFDIR)/EPS0003.nlfc.dat && echo "PASSED: $(TARGET) EPS0003.nlfc.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0003.nlfc.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0004.nlfc.dat $(REFDIR)/EPS0004.nlfc.dat && echo "PASSED: $(TARGET) EPS0004.nlfc.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0004.nlfc.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0001.dat $(REFDIR)/EPS0001.dat && echo "PASSED: $(TARGET) EPS0001.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0001.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0002.dat $(REFDIR)/EPS0002.dat && echo "PASSED: $(TARGET) EPS0002.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0002.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0003.dat $(REFDIR)/EPS0003.dat && echo "PASSED: $(TARGET) EPS0003.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0003.dat" >> ../summary.txt
	@-$(bindir)/diffnum EPS0004.dat $(REFDIR)/EPS0004.dat && echo "PASSED: $(TARGET) EPS0004.dat" >> ../summary.txt || echo "FAILED: $(TARGET) EPS0004.dat" >> ../summary.txt

