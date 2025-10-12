from comp import test1_check,test2_check,runprogs,compeval,rmfiles
def test(args,bindir,testdir,workdir):
        lmfa= f'mpirun -np 1 {bindir}/lmfa '
        lmf = f'mpirun -np {args.np} {bindir}/lmf '
        outfile='out.lmf.eras'
        message1='''
        # Case ErAs: Test of LDA+U
        # --- Test 1.  Basic check of programs lmfa,lmf ---
         The ErAs test also illustrates the LDA+U implementation for:
         1. a case when the LDA puts f orbitals at the Fermi energy.
         2. demonstration of U on both d and f orbitals
         (3. convergence to a metastable solution with a reasonable spin moment but wrong orbital moment.)
        '''
        print(message1)
        rmfiles(workdir,[outfile])
        runprogs([
                 lmfa+" eras  > "+outfile,
                 lmf+" -vnit=1 --pr51 eras >> "+outfile,
                 lmf+" -vnit=3 eras        >> "+outfile 
        ])
        tall=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
        return tall
