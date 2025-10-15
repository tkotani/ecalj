from comp import runprogs,diffnum,dqpu,rmfiles
def test(args,bindir,testdir,workdir):
        gwsc0= bindir + f'/gwsc 0 -np {args.np} '
        tall=''
        out1= ["QPU"]
        out2='log.si'
        rmfiles(workdir,out1+[out2])
        runprogs([
                 gwsc0+ " si"
        ])
        for out in out1:   
                tall+=dqpu(testdir+'/'+out, workdir+'/'+out)
        tall+=diffnum(testdir+'/'+out2, workdir+'/'+out2,tol=3e-3,comparekeys=['fp evl'])
        return tall
