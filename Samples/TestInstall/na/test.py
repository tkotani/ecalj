import os
from comp import test1_check,runprogs,compeval,rmfiles
def test(args,bindir,testdir,workdir):
        lmfa= f'mpirun -np 1 {bindir}/lmfa '
        lmf = f'mpirun -np {args.np} {bindir}/lmf '
        tall=''
        outfile='out.lmf.na'
        message1='''
        # Case na: illustration of low- and high-lying local orbitals
        # --- Test 1.  Basic check of programs lmfa,lmf ---
         The na test also illustrates the following:
         1.  compare the total energy for conventional and extended Na 2p orbitals
             After the test finishes, compare the three energies with out.lmf.na
        '''
        print(message1)
        rmfiles(workdir,[outfile])
        runprogs([
                 lmfa+" na -vnk=6 -vnapval=0 -vpz1=0 -vp1=3.38  > "+outfile ,
                 lmf+ " na -vnk=6 -vnapval=0 -vpz1=0 -vp1=3.38 >>"+outfile,
                 "rm *mixm.na rst.na",
                 lmfa+" na -vnk=6  -vnapval=1 >>"+outfile,
                 lmf+ " na -vnk=6  -vnapval=1>>"+outfile,
                 "rm *mixm.na rst.na",
                 lmf+ " na -vnk=6 -vnapval=2 -vpz1=12.94 >>"+outfile
        ])
        tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
        return tall
    
