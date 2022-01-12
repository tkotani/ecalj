ecalroot=../..

# default mpi size, it is overwritten by envirionmen variable
export mpi_size = 4
# path to the executable programs or scripts
export bindir = $(CURDIR)/bin

# default work directory, it is overwritten by envirionmen variable
export work=work


# all lm tests provided by this makefile.
LM_TEST  = copt te zrt co cr3si6 felz gasls eras c crn cu na
# fe gdn

# all gw tests provided by this makefile.
GW_TEST  = gas_eps_lmfh gas_epsPP_lmfh fe_epsPP_lmfh_chipm si_gw_lmfh
GW_TEST += gas_pw_gw_lmfh si_gwsc gas_gwsc nio_gwsc fe_gwsc ni_crpa srvo3_crpa
#GW_TEST += nio_gwsc444 yh3fcc_gwsc666 gas_gwsc666 pdo_gwsc443 cugase2_gwsc222

# programs under tools dir, required in the tests
TL_DIR   = $(ecalroot)/TOOLS
TL_PROG  = diffnum diffnum2 diffnum0.py comp compall comp.eval

# programs under lm7k dir, required in the tests
LM_DIR   = $(ecalroot)/SRC/exec
LM_PROG  = lmfa lmf-MPIK lmfgw-MPIK run_arg job_pdos job_tdos

# programs under fpgw dir, required in the tests
GW_DIR   = $(ecalroot)/SRC/exec
GW_PROG  = gwsc qg4gw hvccfp0 hsfp0_sc hx0fp0_sc hqpe_sc hmaxloc hpsig_MPI huumat_MPI hwmatK_MPI hx0init hhilbert hrcxq hx0zmel
GW_PROG += rdata4gw_v2 heftet hbasfp0 dqpu
GW_PROG += gw_lmfh hx0fp0 hsfp0 hqpe eps_lmfh epsPP_lmfh epsPP_lmfh_chipm 

MLWF_DIR = $(ecalroot)/SRC/exec
MLWF_PROG  = genMLWFx 
###########################################################
# Do not edit lines below.

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
	############################################
	# test: $@ target: $(TARGET) work: $(work)/$(TARGET)
	############################################
ifneq ($(checkonly),yes)
	# make a work directory.
	mkdir -p $(work)
	# clear the work target directory if it exists.
	rm -rf $(work)/$(TARGET)
	# copy Makefile.define to the work.
	cp -f Makefile.define $(work)
	# copy target directory to the work.
	cp -r $(TARGET) $(work)
	@cd $(work)/$(TARGET) # only effective in shell script
	# exec Makefile in the work target directory.
	$(MAKE) -C $(work)/$(TARGET) $@.exec \
		REFDIR=$(CURDIR)/$(TARGET)
	@cd $(CURDIR) # only effective in shell script
else
	@cd $(work)/$(TARGET) # only effective in shell script
	# exec Makefile in the work target directory.
	$(MAKE) -C $(work)/$(TARGET) $@.check \
		REFDIR=$(CURDIR)/$(TARGET)
	@cd $(CURDIR) # only effective in shell script
endif	
	@echo 
	cat  $(work)/summary.txt
	@echo 
	@-grep 'FAILED' $(work)/summary.txt > /dev/null && \
	echo "   FAILED in some tests! See " $(work)/summary.txt  || \
	echo "   OK! All passed! "
