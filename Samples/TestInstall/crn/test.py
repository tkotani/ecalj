from comp import test2_check,runprogs
def test(args,bindir,testdir,workdir):
        lmfa= f'mpirun -np 1 {bindir}/lmfa '
        lmf = f'mpirun -np {args.np} {bindir}/lmf '
        tall=''
        outfile='out.lmf-dos.crn'
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
        runprogs([
                 lmfa+" crn > "+outfile ,
                 lmf+ " crn >>"+outfile,
                 lmf+ " -vnit=1 -vmetal=2 crn >>"+outfile,
                 lmf+ " --cls crn >>"+outfile
        ])
        outfile='dos-vcdmel.crn'
        tall+=test2_check(testdir+'/'+outfile, workdir+'/'+outfile)
        return tall
