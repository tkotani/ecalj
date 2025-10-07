from comp import runprogs,diffnum,dqpu
def test(args,bindir,testdir,workdir):
        gwsc0= bindir + f'/gwsc 0 -np {args.np} '
        tall=''
        runprogs([
                 "rm -f log.fe QPU QPD",
                 gwsc0+ " fe",
        ])
        dfile="QPU QPD"
        for outfile in dfile.split():   
                tall+=dqpu(testdir+'/'+outfile, workdir+'/'+outfile)
        outfile='log.fe'
        tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=3e-3,comparekeys=['fp evl'])
        return tall
