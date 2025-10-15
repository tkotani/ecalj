from comp import runprogs,diffnum,rmfiles
def test(args,bindir,testdir,workdir):
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    epsPP_lmfh_chipm= f'{bindir}/epsPP_lmfh_chipm -np {args.np}'
    tall=''
    epsfile="ChiPM0001.nlfc.mat ChiPM0002.nlfc.mat ChiPM0003.nlfc.mat ChiPM0004.nlfc.mat ChiPM0005.nlfc.mat"
    rmfiles(workdir,epsfile.split())
    runprogs([
                 lmfa+" fe  > llmfa" ,
                 lmf+ " fe  > llmf",
                 epsPP_lmfh_chipm+" fe" 
    ])
    for outfile in epsfile.split(): 
        tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=3e-3,comparekeys=[])
    return tall
