from comp import runprogs,diffnum,dqpu,rmfiles
def test(args,bindir,testdir,workdir):
        # NiO, QSGW with the antiferromagnetic symmetry (symgrpaf): one iteration from LDA on the 2x2x2 mesh.
        # The same input without symgrpaf and af is Samples/TestInstall/nio_gwsc; the two agree (../README.md, table 1).
        gwsc1= bindir + f'/gwsc 1 -np {args.np} '
        tall=''
        out1=["QPU"]
        out2='log.nio'
        rmfiles(workdir,out1+[out2])
        runprogs([
                 gwsc1+ " nio" + f' {args.run_args}',
        ])
        if not args.mp:
                for outfile in out1:
                        tall+=dqpu(testdir+'/'+outfile, workdir+'/'+outfile)
        tol_log = 5e-3 if args.mp else 3e-3
        tall+=diffnum(testdir+'/'+out2, workdir+'/'+out2,tol=tol_log,comparekeys=['fp evl'])
        return tall
