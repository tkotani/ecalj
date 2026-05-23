from comp import runprogs,diffnum,dqpu,rmfiles
def test(args,bindir,testdir,workdir):
        gwsc0= bindir + f'/gwsc 0 -np {args.np} '
        tall=''
        out1=["QPU.1run"]
        rmfiles(workdir,out1)
        runprogs([
                 gwsc0+ " inas2gasb2" + f' {args.run_args}',
        ])
        if not args.mp:
                for outfile in out1:
                        tall+=dqpu(testdir+'/QPU.1run', workdir+'/'+outfile)
        return tall
