import os
from comp import test1_check,test2_check,runprogs,diffnum
def test(args,bindir,testdir,workdir):
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    eps_lmfh= f'{bindir}/eps_lmfh -np  {args.np} '
    epsPP_lmfh= f'{bindir}/epsPP_lmfh -np  {args.np} '
    tall=''
    runprogs([
                lmfa+" gas  > llmfa" ,
                lmf+ " gas  > llmf",
                "rm -rf EPS*",
                eps_lmfh+" gas" 
    ])
    epsfile="EPS0001.nlfc.dat EPS0002.nlfc.dat EPS0003.nlfc.dat EPS0004.nlfc.dat EPS0001.dat EPS0002.dat EPS0003.dat EPS0004.dat"
    for outfile in epsfile.split(): tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=1e-3,comparekeys=[])
    return tall
