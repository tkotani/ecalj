Usage:
>make target mpi_size=8
1.Default mpi_size=4

>make target mpi_size=8
2.targets are name of directories as 
 LM_TEST  = copt te zrt co cr3si6 felz gasls eras c crn cu na
 GW_TEST  = gas_eps_lmfh gas_epsPP_lmfh fe_epsPP_lmfh_chipm si_gw_lmfh
 GW_TEST += gas_pw_gw_lmfh si_gwsc gas_gwsc nio_gwsc fe_gwsc ni_crpa srvo3_crpa
or choose one of
all
lmall
gwall

When you run
>make cr3si6
This Makefile automatically copy all binaries and scripts to
TestInstall/bin/ at first (only new ones).
Then run calculations with binaries at TestInstall/bin/ 

>make foobar 
foobar are
 show-target
 show-gwtarget
 show-lmtarget

options are
 checkonly=yes
 
 
=======================================
To add new sample, MgO, for example,
perform a calculation at first. Then
 1. Make directory MgO/
 2. Set ctrl file and so on, which are required for the calculation.
 3. To set test run schedule and comparion test, write MgO/makefile.
    Look into other samples as template. This makefile is invoked at work/MgO
    by TestInstall/Makefile.
 4. Add name of the directory to LM_TEST or GW_TEST in Makefile here.

That is, a test is controlled by two files
* TestInstall/Makefile
* MgO/makefile
 
