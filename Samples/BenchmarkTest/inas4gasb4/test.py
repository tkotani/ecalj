from comp import runprogs,diffnum,dqpu,rmfiles
def test(args,bindir,testdir,workdir):
        gwsc0= bindir + f'/gwsc 0 -np {args.np} '
        tall=''
        out1=["QPU"]
        out2='log.inas4gasb4'
        rmfiles(workdir,out1+[out2])
        runprogs([
                 gwsc0+ " inas4gasb4" + f' {args.run_args}',
        ])
        if not args.mp:
                for outfile in out1:
                        tall+=dqpu(testdir+'/QPU.1run', workdir+'/'+outfile)
        return tall
