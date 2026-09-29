from comp import runprogs,diffnum,dqpu,rmfiles
def test(args,bindir,testdir,workdir):
        # Si at the electron temperature 3000 K (t_tetrakbt = t_sigmaw = 3000): one QSGW iteration from LDA on the 4x4x4 mesh.
        # A semiconductor at finite temperature: the Fermi level EFERMI_kbt lies in the gap, and carriers excited across the
        # gap screen (README.md: the gap is 0.16 eV smaller than at 300 K after five iterations).
        gwsc1= bindir + f'/gwsc 1 -np {args.np} '
        tall=''
        out1= ["QPU"]
        out2='log.si'
        rmfiles(workdir,out1+[out2])
        runprogs([
                 "sed -i -e 's/^t_tetrakbt .*/t_tetrakbt    = 3000/' -e 's/^t_sigmaw .*/t_sigmaw      = 3000/' ctrlg.si.toml",
                 gwsc1+ " si" + f' {args.run_args}'
        ])
        if not args.mp:
                for out in out1:
                        tall+=dqpu(testdir+'/'+out, workdir+'/'+out)
        tall+=diffnum(testdir+'/'+out2, workdir+'/'+out2,tol=3e-3,comparekeys=['fp evl'])
        return tall
