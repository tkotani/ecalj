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
                 lmfa+" na  > "+outfile ,
                 lmf+ " na >>"+outfile,
                 "rm *mixm.na rst.na",
                 lmfa+" na --ctrlg:spec.1.p=[3.7,3.5,3.2,4.12,5.1] --ctrlg:spec.1.pz=[0.0,2.9] >>"+outfile,
                 lmf+ " na --ctrlg:spec.1.p=[3.7,3.5,3.2,4.12,5.1] --ctrlg:spec.1.pz=[0.0,2.9]>>"+outfile,
                 "rm *mixm.na rst.na",
                 lmf+ " na --ctrlg:spec.1.p=[3.7,3.5,3.2,4.12,5.1] --ctrlg:spec.1.pz=[0.0,12.94] >>"+outfile
        ])
        tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
        return tall
    
