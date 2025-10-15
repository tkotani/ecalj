from comp import runprogs,dqpu,rmfiles
def test(args,bindir,testdir,workdir):
        lmfa= f'mpirun -np 1 {bindir}/lmfa '
        lmf = f'mpirun -np {args.np} {bindir}/lmf '
        gw_lmfh= f'{bindir}/gw_lmfh -np {args.np} '
        out1=["QPU"]
        tall=''
        rmfiles(workdir,out1)
        runprogs([
                        lmfa+" gas  > llmfa" ,
                        lmf+ " gas  > llmf",
                        gw_lmfh+" gas" 
        ])
        for outfile in out1: 
                tall+=dqpu(testdir+'/'+outfile, workdir+'/'+outfile)
        return tall
