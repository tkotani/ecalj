from comp import test2_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir):
        lmfa= f'mpirun -np 1 {bindir}/lmfa '
        lmf = f'mpirun -np {args.np} {bindir}/lmf '
        tall=''
        out1='out.lmf-dos.crn'
        out2='dos-vcdmel.crn'
        message1='''
        # Case CrN: test of CLS with core hole
        # --- Test 2.  Core-level spectroscopy (EELS), Mulliken analysis, partial DOS ---
         The CrN test case generates core-level spectroscopy for the
         1s state in N.  The self-consistent calculation proceeds with
         an electron missing from the N 1s core, which corresponds to
         the 'sudden approximation' (system relaxes instantanously
         from electron exited out of hole).
        '''
        print(message1)
        rmfiles(workdir,[out1,out2])
        runprogs([
                 lmfa+" crn > "+out1 ,
                 lmf+ " crn >>"+out1,
                 lmf+ " -vnit=1 -vmetal=2 crn >>"+out1,
                 lmf+ " --cls crn >>"+out1
        ])
        tall+=test2_check(testdir+'/'+out2, workdir+'/'+out2)
        return tall
