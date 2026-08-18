from comp import runprogs,diffnum,dqpu,rmfiles
def test(args,bindir,testdir,workdir):
        gwsc0= bindir + f'/gwsc 0 -np {args.np} '
        tall=''
        out1=["QPU"]
        rmfiles(workdir,out1)
        runprogs([
                 gwsc0+ " inas4gasb4" + f' {args.run_args}',
        ])
        if not args.mp:
                for outfile in out1:
                        tall+=dqpu(testdir+'/QPU.1run', workdir+'/QPU')
        return tall
