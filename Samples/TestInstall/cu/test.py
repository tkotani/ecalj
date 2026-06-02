import os
from comp import test1_check,test2_check,runprogs,compeval,rmfiles
def test(args,bindir,testdir,workdir):
        lmfa= f'mpirun -np 1 {bindir}/lmfa '
        lmf = f'mpirun -np {args.np} {bindir}/lmf '
        tall=''
        out1='out.lmf.cu'
        out2='bnd001.spin1 bnd002.spin1 bnd003.spin1 bnd004.spin1'
        out2=out2.split()
        message1='''
        # Case cu: illustration of high-lying local orbitals and bands of Cu up to ~50 eV.
        # --- Test 1.  Basic check of programs lmfa,lmf ---
         The cu test also illustrates the following:
         1.  high-lying local orbitals (Cu 5s,5p,4d are included as local orbitals)
         2.  METAL=3 for BZ integration
         3.  bands mode (see command-line argument --band in last lmf invocation)
        '''
        print(message1)
        rmfiles(workdir,[out1]+out2)
        runprogs([
                 lmfa+" cu  > "+out1 ,
                 lmf+ " cu --ctrlg:bz.nkabc=[8,8,8] >>"+out1,
                 "rm *mixm.cu",
                 lmf+ " cu --ctrlg:bz.nkabc=[8,8,8] --ctrlg:spec.1.rsmh2=[1.3,0.0,1.0,1.3,0.0] --ctrlg:spec.1.lmxa=4 >>"+out1,
                 lmf+ " cu --ctrlg:bz.nkabc=[8,8,8] --ctrlg:spec.1.rsmh2=[1.3,0.0,1.0,1.3,0.0] --ctrlg:spec.1.lmxa=4 --band >>"+out1
        ])
        tall+=test1_check(testdir+'/'+out1, workdir+'/'+out1)
        for out in out2:
                tall+=test2_check(testdir+'/'+out, workdir+'/'+out)
        return tall
