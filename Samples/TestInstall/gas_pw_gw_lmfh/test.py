from comp import runprogs,dqpu
def test(args,bindir,testdir,workdir):
        lmfa= f'mpirun -np 1 {bindir}/lmfa '
        lmf = f'mpirun -np {args.np} {bindir}/lmf '
        gw_lmfh= f'{bindir}/gw_lmfh -np {args.np} '
        tall=''
        runprogs([
                        lmfa+" gas  > llmfa" ,
                        lmf+ " gas  > llmf",
                        "rm QPU",
                        gw_lmfh+" gas" 
        ])
        epsfile="QPU"
        for outfile in epsfile.split(): 
                tall+=dqpu(testdir+'/'+outfile, workdir+'/'+outfile)
        return tall
