import os
from comp import test1_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir):
    outfile='out.lmf.copt'
    message1='''
    # Case copt: a distorted L12 environment with four atoms. 
    # --- Basic check of programs lmfa,lmf ---
    '''
    print(message1)
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    rmfiles(workdir,[outfile])
    runprogs([
        lmfa+'copt --ctrlg:ham.nspin=2 --ctrlg:verbose=41 --ctrlg:bz.metal=3 --ctrlg:bz.tetra=0 --ctrlg:bz.nkabc=[2,2,2] --ctrlg:ham.forces=12  --ctrlg:iter.nit=3 --ctrlg:time=[5,999] > '+outfile,
        lmf +'copt --ctrlg:ham.nspin=2 --ctrlg:verbose=41 --ctrlg:bz.metal=3 --ctrlg:bz.tetra=0 --ctrlg:bz.nkabc=[2,2,2] --ctrlg:ham.forces=12  --ctrlg:iter.nit=3 --ctrlg:time=[5,999] > '+outfile 
    ])
    result= test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    return result
