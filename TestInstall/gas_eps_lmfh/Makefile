include ../Makefile.define

# directroy where input/output files are saved.
TARGET=gas_eps_lmfh
# material in the tests.
MATERIAL=gas

all: help

help:
	@echo "Do not call this Makefile alone."

# exec default case for $(TARGET).
$(TARGET).exec:
	# exec lmfa command.
	$(bindir)/lmfa $(MATERIAL) > llmfa
	# exec lmf command.
	$(bindir)/lmf $(MATERIAL) > llmf
	# exec eps_lmfh command.
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

# check default case for $(TARGET).
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

# exec case1 for $(TARGET).
$(TARGET).case1.exec:
	@echo "nothing to exec for $@"

# check case1 for $(TARGET).
$(TARGET).case1.check:
	@echo "nothing to check for $@"
