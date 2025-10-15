from comp import test1_check,test2_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir):
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    message1='''
    # Case te: molecular statics in an open structure
    # --- Test 1.  Basic check of programs lmfa,lmf ---
    The te test also checks and illustrates the following:

    1.  a simple molecular relaxation with one symmetry-allowed degree of freedom
    2.  an spd-sp basis augmented by floating orbitals, and also by q-dependent APWs
    3.  Comparison of forces and total energy with and without floating orbitals, and with and without APWs.

    lmf will first relax the atoms with the basis including floating orbitals.

    After relaxation, the a new calculation is performed that remove floating orbitals
    but adding plane waves to the basis, so the total energy and forces may be compared. 
    The basis size is variable, but averages ~80 orbitals, a little more than the floating
    orbitals case (~70 orbitals). About 3 mRy is gained relative to the floating orbitals case.

    Note that KMXA=5 is used with the PW calculation.  It isn't necessary in this case; still the
    user is cautioned to monitor this parameter when using energy cutoffs higher than 3 Ry or so.

    As a last step the calculation is repeated at the relaxed position with only atom-centered
    MTO's (neither floating orbitals nor plane waves).  Elimination of floating orbitals reduces 
    the basis from 75 to 39 orbitals, and reduces the LDA total energy by about 4 mRy.

    The forces are not strongly affected by the local orbitals (or APWs), as can be seen by
    looking at the maximum force after the last step.
'''
    print(message1)
    outfile='out.lmf.te'
    rmfiles(workdir,[outfile])
    runprogs([\
        lmfa+"te -vdyn='DYN' -vnk=3 -vnit=3 -vlf1=4 -vlmxl=4 -vnk=3 -vngd=20 -vkmx=3 -vconv=1e-4 > "+outfile, 
        lmf+ "te -vdyn='DYN' -vnk=3 -vnit=3 -vlf1=4 -vlmxl=4 -vnk=3 -vngd=20 -vkmx=3 -vconv=1e-4 -vnbas=12 -vnspec=2 >> "+outfile,
        f"rm -f *mixm.te",
        f"cp rst.te rst.te.bk",
        lmf+"te -vnk=3 -vnit=3 -vlf1=4 -vlmxl=4 -vnk=3 -vngd=20 -vkmx=5 -vconv=1e-4 -vpwmode=11 >> "+outfile,
        f"rm -f *mixm.te",
        f"cp rst.te.bk rst.te",
        lmf+"te -vnk=3 -vnit=3 -vlf1=4 -vlmxl=4 -vnk=3 -vngd=20 -vkmx=3 -vconv=1e-4 >> "+outfile 
    ])
    tall=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    return tall
