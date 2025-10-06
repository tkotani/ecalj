from comp import runprogs,diffnum,dqpu
def test(args,bindir,testdir,workdir):
        gwsc0= bindir + f'/gwsc 0 -np {args.np} '
        tall=''
        runprogs([
                 "rm -f log.nio QPU QPD",
                 gwsc0+ " nio",
        ])
        dfile="QPU"
        for outfile in dfile.split():   
                tall+=dqpu(testdir+'/'+outfile, workdir+'/'+outfile)
        outfile='log.nio'
        tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=3e-3,comparekeys=['fp evl'])
        return tall