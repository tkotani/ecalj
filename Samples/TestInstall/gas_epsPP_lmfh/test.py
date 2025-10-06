import os
from comp import test1_check,test2_check,runprogs,diffnum
def test(args,bindir,testdir,workdir):
    lmfa= bindir +'/lmfa '
    lmf = f'mpirun -np {args.np} '+ bindir +'/lmf '
    eps_lmfh= bindir + f'/eps_lmfh -np {args.np} '
    epsPP_lmfh= bindir + f'/epsPP_lmfh -np {args.np} '
    tall=''
    runprogs([
                lmfa+" gas  > llmfa" ,
                lmf+ " gas  > llmf",
                "rm -rf EPS*",
                epsPP_lmfh+" gas" 
    ])
    epsfile="EPS0001.nlfc.dat EPS0002.nlfc.dat EPS0003.nlfc.dat EPS0004.nlfc.dat"
    for outfile in epsfile.split(): tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=3e-3,comparekeys=[])
    return tall
