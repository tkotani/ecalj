from comp import runprogs,diffnum,dqpu
def test(args,bindir,testdir,workdir):
        gwsc0= bindir + f'/gwsc 0 -np {args.np} '
        tall=''
        runprogs([
                 gwsc0+ " gas",
        ])
        epsfile="QPU"
        for outfile in epsfile.split(): 
                tall+=dqpu(testdir+'/'+outfile, workdir+'/'+outfile)
        outfile='log.gas'
        tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=3e-3,comparekeys=['fp evl'])
        return tall
