from comp import test1_check,test2_check,runprogs
def test(args,bindir,testdir,workdir):
    np4=' -np '+args.np+' '
    job_pdos=bindir+'/job_pdos '
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    outfile='out.lmf.co'
    message1='''
    # Case co: a hexagonal environment with two equivalent atoms.
    # --- Test 1.  Basic check of programs lmfa,lmf ---
        The co test also checks and illustrates the following:
        1.  a spin-polarized case
        2.  tetrahedron method (METAL=3)
        3.  Constrained mixing (first spin is kept frozen, then charge, then both are allowed to change)
        4.  Broyden mixing
        5.  bands mode (see command-line argument --band in last lmf invocation)
            SO coupling and color weights for projection onto Co d channels are included.
            Two separate weights are made (one each for d majority and minority bands.)
            This script will prompt you to see whether it should create 
            a figure from the bands file (The fplot plotting package needs
            to be installed for this function).
        '''       
    print(message1)
    tall=''
    runprogs([
                lmfa+" co -vmet=3 -vlmf=1 -vnk=8 -vnit=3 --pr31 > "+ outfile,
                lmf+ " co -vmet=3 -vlmf=1 -vnk=8 -vnit=3 -vw2=0 --pr31 >> "+outfile,
                "rm *mixm.co",
                lmf+ " co -vmet=3 -vlmf=1 -vnk=8 -vnit=3 -vw1=0 --pr31 >> "+outfile,
                "rm *mixm.co",
                lmf+ " co -vmet=3 -vlmf=1 -vnk=8 -vnit=3 --pr31 --time=5 >> "+outfile,
                lmf+ " co -vmet=3 -vnk=8 -vnit=3 --pr31  -vso=t --band:fn=syml >> "+outfile,
                "rm -f atm.* *mixm.* rst.* save.* log.* *hssn.* wkp.* dos.* tdos.* pdos.* dos-mull.* qpp.* out.lmf-dos*"
    ])
    tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    outbnds='bnds.co'
    tall+=test2_check(testdir+'/'+outbnds, workdir+'/'+outbnds)
    message1='''
    # --- Test 2.  Core-level spectroscopy (EELS), Mulliken analysis, partial DOS ---
        The Co test case illustrates partial dos resolved by both l and m.
        Note: because the routine generating partial DOS within MT spheres
        does not properly symmetrize it, symmetry operations must be
        suppressed when resolving DOS by m.
    '''
    print(message1)
    runprogs([
                lmfa+" co -vmet=3 -vlmf=1 -vnk=8 -vnit=1 --pr31 > out.lmf-dos.co",\
                job_pdos+" co "+ np4 +" -vmet=3 -vlmf=1 -vnk=8 -vnit=1 --pr31 ---NoGnuplot > out.lmf-dos.co"\
    ])
    pdos2002='dos.isp2.site002.co'
    tall+=test2_check(testdir+'/'+pdos2002, workdir+'/'+pdos2002)
    return tall
