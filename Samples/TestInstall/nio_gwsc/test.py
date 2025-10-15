from comp import runprogs,diffnum,dqpu,rmfiles
def test(args,bindir,testdir,workdir):
        gwsc0= bindir + f'/gwsc 0 -np {args.np} '
        tall=''
        out1="QPU"
        out2='log.nio'
        rmfiles(workdir,[out1,out2])
        runprogs([
                 gwsc0+ " nio",
        ])
        for outfile in out1.split():   
                tall+=dqpu(testdir+'/'+outfile, workdir+'/'+outfile)
        tall+=diffnum(testdir+'/'+out2, workdir+'/'+out2,tol=3e-3,comparekeys=['fp evl'])
        return tall
