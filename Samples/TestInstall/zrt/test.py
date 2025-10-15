import os
from comp import test1_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir):
        lmfa= f'mpirun -np 1 {bindir}/lmfa '
        lmf = f'mpirun -np {args.np} {bindir}/lmf '
        outfile='out.lmf.zrt'
        message1='''
        # Case zrt: ZrO_2 fluorite in tetragonal setting, with tetragonal distortion
        # --- Test 1.  Basic check of programs lmfa,lmf ---
        The zrt test also checks and illustrates the following:
                1.  the use of restart files in both binary and ascii forms
                2.  two-kappa basis
        '''
        
        print(message1)
        rmfiles(workdir,[outfile])
        runprogs([
                        lmfa+" -vfp=1 zrt --no-iactiv > out.lmf.zrt > "+outfile,
                        lmf+"  -vnitq=1 -vforce=1 -vfp=1 zrt >> "+outfile,
                        lmf+"  -vnitq=1 -vforce=1 -vfp=1 zrt >> "+outfile 
        ])
        tall=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
        return tall
