==============================
How to add new test ?
1. make a test calculation at a directory.
2. Make a test directory ecalj/TestInstall/foobar/.
   Then set keep ctrl.* and so on in it (See other test).
3. Write Makefile at ecalj/TestInstall/foobar
   (see other test).
4. Add test in the main Makefile at ecalj/TestInstall.
===========================

NAME
	Makefile - makefile script to test and check ecalj programs.

SYNOPSIS
	make [targets] [mpi_size=size] [checkonly=yes] [work=dir] [bindir=dir]

PREPARATION
	Make sure that fpgw TOOLS have been compiled,
	and its directory has been listed by TL_DIR in this Makefile.

	Make sure that lm7K programs have been compiled,
	and its directory has been listed by LM_DIR in this Makefile.

	Make sure that fpgw programs have been compiled,
	and its directory has been listed by GW_DIR in this Makefile.

DESCRIPTION
	To show all available test targets,
		make show-target

	To perform a test, for example, copt under lm7k,
		make copt

	To perform two tests successively, for example,
        copt under lm7k and si_gwsc under fpgw,
		make copt si_gwsc

	To perform all tests under lm7k,
		make lmall

	To perform all tests under fpgw,
		make gwall

	To perform all tests under lm7k and fpgw,
		make all


RESULT
	To show summary of results,
		make show-summary

OPTIONS
	To perform programs parallelly using MPI,
		make copt mpi_size=2

	To perform a check without execution,
		make copt checkonly=yes

	To install programs under a specified directory instead of the default,
		make copt bindir=$(HOME)/mybin

	To use a specified work directory instead of the default,
		make copt work=/var/tmp/mywork

NOTE:
	If perl cause 'warning: Falling back to...', and so on,
        set "export PERL_BADLANG=0" to suprese them.

        FA ---> free atom 
        ehf ---> Harris energy 
	ehk --->  K-Sham energy
        max force ---> maximum force 
	mmon --->  magnetic moment
        RMS dq --->  RMS drho
 
