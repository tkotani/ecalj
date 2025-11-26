from comp import test1_check,test2_check,runprogs,compeval,rmfiles
def test(args,bindir,testdir,workdir):
        lmfa= f'mpirun -np 1 {bindir}/lmfa '
        lmf = f'mpirun -np {args.np} {bindir}/lmf '
        outfile='out.lmf.nio'
        message1='''
        # Case NiO: Test of AFsymmetry
        '''
        print(message1)
        rmfiles(workdir,[outfile])
        runprogs([
                 lmfa+" nio > "+outfile,
                 lmf+" -vnit=3 nio >> "+outfile,
        ])
        tall=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
        return tall
