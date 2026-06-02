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
        lmfa+'copt --toml.ham.nspin=2 --toml.verbose=41 --toml.bz.metal=3 --toml.bz.tetra=0 --toml.bz.nkabc=[2,2,2] --toml.ham.forces=12  --toml.iter.nit=3 --toml.time=[5,999] > '+outfile,
        lmf +'copt --toml.ham.nspin=2 --toml.verbose=41 --toml.bz.metal=3 --toml.bz.tetra=0 --toml.bz.nkabc=[2,2,2] --toml.ham.forces=12  --toml.iter.nit=3 --toml.time=[5,999] > '+outfile 
    ])
    result= test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    return result
