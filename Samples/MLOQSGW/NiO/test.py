from comp import runprogs,diffnum,dqpu,rmfiles
def test(args,bindir,testdir,workdir):
        # MLO-QSGW: two iterations from LDA; the 2nd iteration runs on Sigma
        # interpolated through the MLO (ECALJ_MLO_MIX=1, mixbeta=0.5). See ../README.md.
        gwsc = f'ECALJ_MLO_MIX=1 {bindir}/gwsc 2 -np {args.np} --mlo '
        out1= ["QPU","QPD"]
        out2='log.nio'
        tall=''
        rmfiles(workdir,out1+[out2])
        runprogs([
                 gwsc + " nio" + f' {args.run_args}',
        ])
        if not args.mp:
                for outfile in out1:
                        tall+=dqpu(testdir+'/'+outfile, workdir+'/'+outfile)
        tol_log = 5e-3 if args.mp else 3e-3
        tall+=diffnum(testdir+'/'+out2, workdir+'/'+out2,tol=tol_log,comparekeys=['fp evl'])
        return tall
