import os
from comp import test2_check,runprogs,diffnum,rmfiles
def test(args,bindir,testdir,workdir):
    lmfa   = f'mpirun -np 1 {bindir}/lmfa '
    lmf    = f'mpirun -np {args.np} {bindir}/lmf '
    job_eps= f'{bindir}/job_eps -np {args.np} --nognuplot '
    epsfile="EPS0001.nlfc.dat EPS0002.nlfc.dat EPS0003.nlfc.dat EPS0004.nlfc.dat"
    tall=''
    rmfiles(workdir,epsfile.split())
    runprogs([
                lmfa+" gas  > llmfa" ,
                lmf+ " gas  > llmf",
                job_eps+" gas"
    ])
    for outfile in epsfile.split():
        tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=3e-3,comparekeys=[])
    return tall
