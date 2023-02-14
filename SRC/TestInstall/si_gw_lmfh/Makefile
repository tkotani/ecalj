include ../../Makefile

# directroy where input/output files are saved.
TARGET=si_gw_lmfh
# material in the tests.
MATERIAL=si

# exec default case for $(TARGET).
$(TARGET).exec:
	# exec lmfa command.
	$(bindir)/lmfa $(MATERIAL) > llmfa
	# exec lmf command.
	$(bindir)/lmf-MPIK $(MATERIAL) > llmf_lda
	# exec gw_lmfh command.
	rm  QPU
	$(bindir)/gw_lmfh -np $(mpi_size) $(MATERIAL)
	# add summary.txt
	# check a result by dqpu command.
	@echo  >> ../summary.txt
	@-$(bindir)/dqpu $(REFDIR)/QPU QPU && echo "PASSED: $(TARGET) QPU" >> ../summary.txt || echo "FAILED: $(TARGET) QPU compared by ./bin/dqpu" >> ../summary.txt


# check default case for $(TARGET).
$(TARGET).check:
	# add summary.txt
	# check a result by dqpu command.
	@echo  >> ../summary.txt
	@-$(bindir)/dqpu $(REFDIR)/QPU QPU && echo "PASSED: $(TARGET) QPU" >> ../summary.txt || echo "FAILED: $(TARGET) QPU compared by ./bin/dqpu" >> ../summary.txt

# exec case1 for $(TARGET).
$(TARGET).case1.exec:
	@echo "nothing to exec for $@"

# check case1 for $(TARGET).
$(TARGET).case1.check:
	@echo "nothing to check for $@"
