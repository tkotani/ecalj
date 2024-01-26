ecalroot=../..
export mpi_size = 4
export bindir = $(CURDIR)/bin
export work=work
LM_TEST  = copt te zrt co cr3si6 felz gasls eras c crn cu na
GW_TEST  = gas_eps_lmfh gas_epsPP_lmfh fe_epsPP_lmfh_chipm si_gw_lmfh
GW_TEST += gas_pw_gw_lmfh si_gwsc gas_gwsc nio_gwsc fe_gwsc ni_crpa srvo3_crpa
#GW_TEST += nio_gwsc444 yh3fcc_gwsc666 gas_gwsc666 pdo_gwsc443 cugase2_gwsc222
TL_DIR   = $(ecalroot)/TOOLS
TL_PROG  = diffnum diffnum2 diffnum0.py comp compall comp.eval
LM_DIR   = $(ecalroot)/SRC/exec
LM_PROG  = lmfa lmf-MPIK run_arg job_pdos job_tdos ctrl2ctrlp.py a2vec.py
GW_DIR   = $(ecalroot)/SRC/exec
GW_PROG  = gwutil.py gwsc qg4gw hvccfp0 hsfp0_sc hqpe_sc hmaxloc hpsig_MPI huumat_MPI hwmatK_MPI hx0init hrcxq 
GW_PROG += rdata4gw_v2 heftet hbasfp0 dqpu
GW_PROG += gw_lmfh hx0fp0 hsfp0 hqpe eps_lmfh epsPP_lmfh epsPP_lmfh_chipm 
MLWF_DIR = $(ecalroot)/SRC/exec
MLWF_PROG  = genMLWFx 
###########################################################
help:
	@cat README.txt
show-lmtarget:
	@echo "the following tests are provided for lm7k:"
	@echo "  $(LM_TEST)."
show-gwtarget:
	@echo "the following tests are provided for fpgw:"
	@echo "  $(GW_TEST)."
show-target: show-lmtarget show-gwtarget
# for each tool program, copy it from tool dir to bin dir if it was updated
TL_BINPROG = $(addprefix $(bindir)/,$(TL_PROG))
$(TL_BINPROG): $(bindir)/%: $(TL_DIR)/%
	cp -f $< $@
# for each lm7k program, copy it from tool dir to bin dir if it was updated
LM_BINPROG = $(addprefix $(bindir)/,$(LM_PROG))
$(LM_BINPROG): $(bindir)/%: $(LM_DIR)/%
	cp -f $< $@
# for each fpgw program, copy it from tool dir to bin dir if it was updated
GW_BINPROG = $(addprefix $(bindir)/,$(GW_PROG))
$(GW_BINPROG): $(bindir)/%: $(GW_DIR)/%
	cp -f $< $@
# for each fpgw program, copy it from tool dir to bin dir if it was updated
MLWF_BINPROG = $(addprefix $(bindir)/,$(MLWF_PROG))
$(MLWF_BINPROG): $(bindir)/%: $(MLWF_DIR)/%
	cp -f $< $@
# install all programs to bindir
install-program: $(TL_BINPROG) $(LM_BINPROG) $(GW_BINPROG) $(MLWF_BINPROG)
clear-summary:
	rm -f $(work)/summary.txt
show-summary:
	@cat   $(work)/summary.txt
.PHONY: $(LM_TEST)
.PHONY: $(GW_TEST)
lmall: $(LM_TEST)
gwall: $(GW_TEST)
all: lmall gwall
$(LM_TEST) $(GW_TEST): TARGET=$(basename $@)
$(LM_TEST) $(GW_TEST): install-program clear-summary

########### main part of test. call makefile at */LM_TEST #########################
# test: $@ target: $(TARGET) work: $(work)/$(TARGET)
ifneq ($(checkonly),yes)
	mkdir -p $(work)
	rm -rf $(work)/$(TARGET)
	cp -r $(TARGET) $(work)
	@cd $(work)/$(TARGET) # only effective in shell script
	# exec Makefile in the work target directory.
	$(MAKE) -C $(work)/$(TARGET) $@.exec REFDIR=$(CURDIR)/$(TARGET)
	@cd $(CURDIR) # only effective in shell script
else
	@cd $(work)/$(TARGET) # only effective in shell script
	# exec Makefile in the work target directory.
	$(MAKE) -C $(work)/$(TARGET) $@.check REFDIR=$(CURDIR)/$(TARGET)
	@cd $(CURDIR) # only effective in shell script
endif	
	@echo 
	cat  $(work)/summary.txt
	@echo 
	@-grep 'FAILED' $(work)/summary.txt > /dev/null && echo "   FAILED in some tests! See " $(work)/summary.txt  || echo "   OK! All passed! "
####################################################### subroutines
bindir=../../bin/
LMF=$(bindir)/lmf-MPIK
LMFP=mpirun -np $(mpi_size) $(bindir)/lmf-MPIK
comp=test.comp
# test1
defatol1=1e-5
dehf1tol1=1e-5
dfmax1tol1=0.1
dmom1tol1=1e-4
dehf1toln=1e-5
dmomntol1=1e-4
drmsqtol1=1e-4
bndstol=1e-4
dfmaxntol1=0.1
# test2
dosclstol=0.003 #1e-3
dosmulltol=1e-3
pdostol=0.01
# test3
dmom1tol3=1e-4
# test4
dorbmtol=1e-5
# test4c
gmtol=1e-4
define test1.check
	@echo "compare $1 and $2"
	@$(bindir)/comp $1 $2 'FA etot (last species)  ' $(defatol1) 'etot='  > $(comp)
	@$(bindir)/comp $1 $2 '1st  iter ehf.eV       ' $(dehf1tol1) 'ehf\(eV\)=' -v '^h' >> $(comp)
	@$(bindir)/comp $1 $2 '1st  iter ehk.eV       ' $(dehf1tol1) 'ehk\(eV\)=' -v '^h' >> $(comp)
	@$(bindir)/comp $1 $2 '1st  iter mmom          ' $(dmom1tol1) 'mmom=' -v '^h' >> $(comp)
	@$(bindir)/comp $1 $2 '2nd  iter ehf.Ry       ' $(dehf1toln) 'ehf=' -v 'it  2'  >> $(comp)
	@$(bindir)/comp $1 $2 '2nd  iter ehk.Ry       ' $(dehf1toln) 'ehk=' -v 'it  2'  >> $(comp)
	@$(bindir)/comp $1 $2 '9th  iter ehf.Ry       ' $(dehf1toln) 'ehf=' -v 'it  9'  >> $(comp)
	@$(bindir)/comp $1 $2 '9th  iter ehk.Ry       ' $(dehf1toln) 'ehk=' -v 'it  9'  >> $(comp)
	@$(bindir)/comp $1 $2 'last iter x ehf.eV       ' $(dehf1toln) 'ehf\(eV\)=' -v '^x'  >> $(comp)
	@$(bindir)/comp $1 $2 'last iter x ehk.eV       ' $(dehf1toln) 'ehk\(eV\)=' -v '^x'  >> $(comp)
	@$(bindir)/comp $1 $2 'last iter c ehf.eV       ' $(dehf1toln) 'ehf\(eV\)=' -v '^c'  >> $(comp)
	@$(bindir)/comp $1 $2 'last iter c ehk.eV       ' $(dehf1toln) 'ehk\(eV\)=' -v '^c'  >> $(comp)
	@$(bindir)/comp $1 $2 'last iter E(LDA+U).eV   ' $(dehf1toln) 'Etot\(LDA+U\)='  >> $(comp)
	@$(bindir)/comp $1 $2 'last iter max force     ' $(dfmax1tol1) 'Maximum Harris force ='  >> $(comp)
	@$(bindir)/comp $1 $2 'last iter mmom          ' $(dmom1tol1) 'mmom='  >> $(comp)
	@$(bindir)/comp $1 $2 'chk1ch last iter RMS dq ' $(drmsqtol1) 'RMS DQ=' >> $(comp)
	@cat $(comp)
	@-grep 'ERR' $(comp) > /dev/null && echo '$(TARGET) TEST 1 FAILED comparison $1 $2' || echo '$(TARGET) TEST 1 PASSED $2'
	@-grep 'ERR' $(comp) > /dev/null && echo " FAILED: $(TARGET) TEST1 $2" >> ../summary.txt || echo "PASSED: $(TARGET) TEST1 $2" >> ../summary.txt
	@echo "[note] run make show-summary to confirm the result."
endef
define test2.check
	@rm -f $(comp)
	@-$(bindir)/compall $1 $2 $(dosclstol) >> $(comp) 
	@tail -1 $(comp)
	@-grep 'ERR' $(comp) > /dev/null && echo '$(TARGET) TEST2 FAILED comparison $1 $2' || echo '$(TARGET) TEST2 PASSED $2'
	@-grep 'ERR' $(comp) > /dev/null && echo " FAILED: $(TARGET) TEST2 $2" >> ../summary.txt || echo "PASSED: $(TARGET) TEST2 $2" >> ../summary.txt
	@echo "[note] run make show-summary to confirm the result."
endef
define test3.check
	@echo "compare $1 and $2"
	@$(bindir)/comp $1 $2 'FA etot (last species)  ' $(defatol1) 'etot=' -t > $(comp)
	@$(bindir)/comp $1 $2 '1st  iter ehf.eV       ' $(dehf1tol1) 'ehf\(eV\)=' -v '^h'  >> $(comp)
	@$(bindir)/comp $1 $2 '1st  iter ehk.eV       ' $(dehf1tol1) 'ehk\(eV\)=' -v '^h'  >> $(comp)
	@$(bindir)/comp $1 $2 '1st  iter unconst. mmom ' $(dmom1tol3) 'Mag. moment:' >> $(comp)
	@$(bindir)/comp $1 $2 '1st  iter mmom          ' $(dmom1tol1) 'mmom=' >> $(comp)
	@$(bindir)/comp $1 $2 'last iter ehf.eV       ' $(dehf1toln) 'ehf\(eV\)=' -v '^c' -t >> $(comp)
	@$(bindir)/comp $1 $2 'last iter ehk.eV       ' $(dehf1toln) 'ehk\(eV\)=' -v '^c' -t >> $(comp)
	@$(bindir)/comp $1 $2 'last iter max force     ' $(dfmax1tol1) 'Maximum Harris force =' -t >> $(comp)
	@$(bindir)/comp $1 $2 'last iter mmom          ' $(dmom1tol1) 'mmom=' -t >> $(comp)
	@$(bindir)/comp $1 $2 'chk1ch last iter RMS dq ' $(drmsqtol1) 'RMS DQ=' >> $(comp)
	@cat $(comp)
	@-grep 'ERR' $(comp) > /dev/null && \
	echo '$(TARGET) TEST 3 FAILED comparison $1 $2' || \
	echo '$(TARGET) TEST 3 PASSED $2'
	@-grep 'ERR' $(comp) > /dev/null && \
	echo " FAILED: $(TARGET) TEST3 $2" >> ../summary.txt || \
	echo "PASSED: $(TARGET) TEST3 $2" >> ../summary.txt
	@echo "[note] run make show-summary to confirm the result."
endef
define test4a.check
	@echo "compare $1 and $2"
	@$(bindir)/comp $1 $2 'Orbital moment          ' $(dorbmtol) 'total orbital moment   1:' -t > $(comp)
	@$(bindir)/comp $1 $2 'last iter ehf.eV       ' $(dehf1toln) 'ehf\(eV\)=' -v '^   it' -t >> $(comp)
	@$(bindir)/comp $1 $2 'last iter ehk.eV       ' $(dehf1toln) 'ehk\(eV\)=' -v '^   it' -t >> $(comp)
	@$(bindir)/comp $1 $2 'last iter ehf.eV E(MTO+PW)' $(dehf1toln) 'pwmode=[^0].*ehf\(eV\)=' -t >> $(comp)
	@$(bindir)/comp $1 $2 'last iter ehk E(MTO+PW)  ' $(dehf1toln) 'pwmode=[^0].*ehk\(eV\)=' -t >> $(comp)
	@cat $(comp)
	@-grep 'ERR' $(comp) > /dev/null && echo '$(TARGET) TEST 4a FAILED comparison $1 $2' || echo '$(TARGET) TEST 4a PASSED $2'
	@-grep 'ERR' $(comp) > /dev/null && echo " FAILED: $(TARGET) TEST4a $2 " >> ../summary.txt || echo "PASSED: $(TARGET) TEST4a $2 " >> ../summary.txt
	@echo "[note] run make show-summary to confirm the result."
endef
define test4b.check
	@echo "compare $1 and $2"
	@$(bindir)/comp $1 $2 'Orbital moment          ' $(dorbmtol) 'total orbital moment   1:' -t > $(comp)
	@cat $(comp)
	@-grep 'ERR' $(comp) > /dev/null && echo '$(TARGET) TEST 4b FAILED comparison $1 $2' || echo '$(TARGET) TEST 4b PASSED $2'
	@-grep 'ERR' $(comp) > /dev/null && echo " FAILED: $(TARGET) TEST4b $2" >> ../summary.txt || echo "PASSED: $(TARGET) TEST4b $2" >> ../summary.txt
	@echo "[note] run make show-summary to confirm the result."
endef
define test4c.check
	$(bindir)/comp.eval $1 $2 ' 0.00000  0.00000  0.00000' $(lineeval) $(evalso) $(gmtol) > $(comp)
	@cat $(comp)
	@-grep 'ERR' $(comp) > /dev/null && echo '$(TARGET) TEST 4c FAILED comparison $1 $2' || echo '$(TARGET) TEST 4c PASSED $2'
	@-grep 'ERR' $(comp) > /dev/null && echo " FAILED: $(TARGET) TEST4c $2" >> ../summary.txt || echo "PASSED: $(TARGET) TEST4c $2" >> ../summary.txt
	@echo "[note] run make show-summary to confirm the result."
endef
