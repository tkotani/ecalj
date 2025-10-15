from comp import test1_check,test2_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir):
    np4=' -np '+args.np+' '
    job_pdos=bindir+'/job_pdos '
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    out1='out.lmf.co'
    out2='bnd001.spin1 bnd002.spin1  bnd003.spin1  bnd004.spin1  bnd005.spin1'
    out2=out2.split()
    out3='dos.isp2.site002.co'
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
    rmfiles(workdir,[out1,out3]+out2)
    runprogs([
                lmfa+" co -vmet=3 -vlmf=1 -vnk=8 -vnit=3 --pr31 > "+ out1,
                lmf+ " co -vmet=3 -vlmf=1 -vnk=8 -vnit=3 -vw2=0 --pr31 >> "+out1,
                "rm *mixm.co",
                lmf+ " co -vmet=3 -vlmf=1 -vnk=8 -vnit=3 -vw1=0 --pr31 >> "+out1,
                "rm *mixm.co",
                lmf+ " co -vmet=3 -vlmf=1 -vnk=8 -vnit=3 --pr31 --time=5 >> "+out1,
                lmf+ " co -vmet=3 -vnk=8 -vnit=3 --pr31  -vso=t --band:fn=syml >> "+out1,
                "rm -f atm.* *mixm.* rst.* save.* log.* *hssn.* wkp.* dos.* tdos.* pdos.* dos-mull.* qpp.* out.lmf-dos*"
    ])
    tall+=test1_check(testdir+'/'+out1, workdir+'/'+out1)
    for out in out2:
        print(out,end=" ")
        tall+=test2_check(testdir+'/'+out, workdir+'/'+out)
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
    tall+=test2_check(testdir+'/'+out3, workdir+'/'+out3)
    return tall
