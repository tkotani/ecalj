from comp import runprogs,dqpu,rmfiles
def test(args,bindir,testdir,workdir):
        lmfa= f'mpirun -np 1 {bindir}/lmfa '
        lmf = f'mpirun -np {args.np} {bindir}/lmf '
        gw_lmfh= f'{bindir}/gw_lmfh -np {args.np} '
        tall=''
        out= ["QPU"]
        rmfiles(workdir,out)
        runprogs([
                lmfa+" si  > llmfa" ,
                lmf+ " si  > llmf",
                gw_lmfh+" si" 
        ])
        for outfile in out:   
                tall+=dqpu(testdir+'/'+outfile, workdir+'/'+outfile)
        return tall
