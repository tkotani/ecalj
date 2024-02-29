include ../../Makefile

TARGET=co

out1=out.lmf.$(TARGET)
out2=out.lmf-dos.$(TARGET)
outbnds=bnds.$(TARGET)
pdos1001=dos.isp1.site001.$(TARGET)
pdos1002=dos.isp1.site002.$(TARGET)
pdos2001=dos.isp2.site001.$(TARGET)
pdos2002=dos.isp2.site002.$(TARGET)

###########################################################
#----- subroutine test1.exec
define test1.exec
	#-------------------------------------------
	# Case co: a hexagonal environment with two equivalent atoms.
	#-------------------------------------------
	#
	# --- Test 1.  Basic check of programs lmfa,lmf ---
	#
	@cat message.test1.$(TARGET)
	mpirun -np 1 $(bindir)/lmfa  $(TARGET) -vmet=3 -vlmf=1 -vnk=8 -vnit=3 --pr31 > $1 2>&1
	$(LMFP)   $(TARGET) -vmet=3 -vlmf=1 -vnk=8 -vnit=3 -vw2=0 --pr31 >> $1 2>&1
	rm mixm.$(TARGET)
	$(LMFP)   $(TARGET) -vmet=3 -vlmf=1 -vnk=8 -vnit=3 -vw1=0 --pr31 >> $1 2>&1
	rm mixm.$(TARGET)
	$(LMFP)  $(TARGET) -vmet=3 -vlmf=1 -vnk=8 -vnit=3 --pr31 --time=5 >> $1 2>&1
	$(LMFP)  $(TARGET) -vmet=3 -vnk=8 -vnit=3 --pr31  -vso=t --band:fn=syml >> $1 2>&1 
endef

#----- subroutine test2.exec
define test2.exec
	#
	# --- Test 2.  Core-level spectroscopy (EELS), Mulliken analysis, partial DOS ---
	#
	@cat message.test2.$(TARGET)
	mpirun -np 1 $(bindir)/lmfa  $(TARGET) -vmet=3 -vlmf=1 -vnk=8 -vnit=1 \
		--pr31 --no-iactiv > $1 2>&1
        $(bindir)/job_pdos $(TARGET)  -np $(mpi_size) -vmet=3 -vlmf=1 -vnk=8 -vnit=1 \
		--pr31 --no-iactiv ---NoGnuplot > $1 2>&1
endef

#----- subroutine plot.plbnds
# define plot.plbnds
# 	echo -10,5,5,10 | $(bindir)/plbnds -fplot -ef=0 -scl=13.6 \
# 		-lt=1,col=0,0,1,colw=1,0,0,colw2=0,1,0 -lbl=M,G,A,L,G,K \
# 		bnds.$(TARGET) > plot.log 2>&1
# 	$(bindir)/fplot -disp -f plot.plbnds >> plot.log 2>&1
# 	mv ps.dat ps.plbnds.dat
# 	@echo "[note] run make plot-band to plot band"
# endef

# #----- subroutine plot.pldos
# define plot.pldos
# 	echo 30 15 -.8 .7 | $(bindir)/pldos '-lst=1;3,5,7;9,11,13,15,17' \
# 		'-lst2=2;4,6,8;10,12,14,16,18' -fplot \
# 		pdos.$(TARGET) > plot.log 2>&1
# 	$(bindir)/fplot -disp -pr10 -f plot.dos >> plot.log 2>&1
# 	mv ps.dat ps.pldos.dat
# 	@echo "[note] run make plot-dos to plot dos"
# endef


###########################################################
# exec default case for $(TARGET).
$(TARGET).exec:
	$(call test1.exec,$(out1))
	$(call test1.check,$(REFDIR)/$(out1),$(out1))
	rm -f atm.* mixm.* rst.* save.* log.* hssn.* wkp.* dos.* tdos.* pdos.* dos-mull.* qpp.* out.lmf-dos*
	$(call test2.exec,$(out2))
	$(call test2.check,$(REFDIR)/$(outbnds),$(outbnds))
	$(call test2.check,$(REFDIR)/$(pdos2002),$(pdos2002))

#	$(call test2.check,$(REFDIR)/$(pdos1001),$(pdos1001))
#	$(call test2.check,$(REFDIR)/$(pdos1002),$(pdos1002))
#	$(call test2.check,$(REFDIR)/$(pdos2001),$(pdos2001))

# check
$(TARGET).check:
	$(call test1.check,$(REFDIR)/$(out1),$(out1))
	$(call test2.check,$(REFDIR)/$(outbnds),$(outbnds))
	$(call test2.check,$(REFDIR)/$(pdos2002),$(pdos2002))

#	$(call test2.check,$(REFDIR)/$(pdos1001),$(pdos1001))
#	$(call test2.check,$(REFDIR)/$(pdos1002),$(pdos1002))
#	$(call test2.check,$(REFDIR)/$(pdos2001),$(pdos2001))
#	$(call plot.plbnds)
#	$(call plot.pldos)


# exec case1 for $(TARGET).
#$(TARGET).case1.exec:
#	@echo "nothing to exec for $@"

# check case1 for $(TARGET).
#$(TARGET).case1.check:
#	@echo "nothing to check for $@"

# # plot graphs on the display
# plot-band:
# 	ghostscript ps.plbnds.dat

# plot-dos:
# 	ghostscript ps.pldos.dat


