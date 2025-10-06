import os
from comp import test1_check,test2_check,runprogs,compeval
def test(args,bindir,testdir,workdir):
        lmfa= f'mpirun -np 1 {bindir}/lmfa '
        lmf = f'mpirun -np {args.np} {bindir}/lmf '
        tall=''
        outfile='out.lmf.cu'
        message1='''
        # Case cu: illustration of high-lying local orbitals and bands of Cu up to ~50 eV.
        # --- Test 1.  Basic check of programs lmfa,lmf ---
         The cu test also illustrates the following:
         1.  high-lying local orbitals (Cu 5s,5p,4d are included as local orbitals)
         2.  METAL=3 for BZ integration
         3.  bands mode (see command-line argument --band in last lmf invocation)
        '''
        print(message1)
        runprogs([
                 lmfa+" cu  > "+outfile ,
                 lmf+ " cu -vnk=8 -vbigbas=f >>"+outfile,
                 "rm *mixm.cu",
                 lmf+ " cu -vnk=8 -vbigbas=t -vmetal=3 -vrsm2=1.3 -vrsmd1x=1 -vlmx=4 -vpwmode=0 -voveps=0d-7 >>"+outfile,
                 lmf+ " cu -vnk=8 -vbigbas=t -vmetal=3 -vrsm2=1.3 -vrsmd1x=1 -vlmx=4 -vpwmode=0 -voveps=0d-7 --band:fn=syml >>"+outfile
        ])
        tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
        outfile='bnds.cu'
        tall+=test2_check(testdir+'/'+outfile, workdir+'/'+outfile)
        return tall
