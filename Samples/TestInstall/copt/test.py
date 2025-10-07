import os
from comp import test1_check,test2_check,runprogs
def test(args,bindir,testdir,workdir):
    outfile='out.lmf.copt'
    message1='''
    # Case copt: a distorted L12 environment with four atoms. 
    # --- Basic check of programs lmfa,lmf ---
    '''
    print(message1)
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    runprogs([
        lmfa+'copt -vnsp=2 -cstrx=l12 --pr41 -vmet=3 -vtetra=0 -vnk1=2 -vlfrce=12 -vdist=1 -vnit=3 --time=5 > '+outfile,
        lmf +'copt -vnsp=2 -cstrx=l12 --pr41 -vmet=3 -vtetra=0 -vnk1=2 -vlfrce=12 -vdist=1 -vnit=3 --time=5 > '+outfile 
    ])
    result= test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    return result
