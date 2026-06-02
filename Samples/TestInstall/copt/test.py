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
        lmfa+'copt -v[ham.nspin]=2 --pr=41-v[bz.metal]=3 -v[bz.tetra]=0 -v[bz.nkabc]=[2,2,2] -v[ham.forces]=12  -v[iter.nit]=3 --time=5 > '+outfile,
        lmf +'copt -v[ham.nspin]=2 --pr=41-v[bz.metal]=3 -v[bz.tetra]=0 -v[bz.nkabc]=[2,2,2] -v[ham.forces]=12  -v[iter.nit]=3 --time=5 > '+outfile 
    ])
    result= test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    return result
