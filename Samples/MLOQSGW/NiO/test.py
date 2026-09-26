from comp import runprogs,diffnum,dqpu,rmfiles
def test(args,bindir,testdir,workdir):
        # MLO-QSGW: two iterations from LDA; the 2nd iteration runs on Sigma
        # interpolated through the MLO (ECALJ_MLO_MIX=1, mixbeta=0.5). See ../README.md.
        gwsc = f'ECALJ_MLO_MIX=1 {bindir}/gwsc 2 -np {args.np} --mlo '
        out1= ["QPU","QPD"]
        out2='log.nio'
        tall=''
        rmfiles(workdir,out1+[out2])
        runprogs([
                 gwsc + " nio" + f' {args.run_args}',
        ])
        if not args.mp:
                for outfile in out1:
                        tall+=dqpu(testdir+'/'+outfile, workdir+'/'+outfile)
        # --mp without --fp32 runs the GPU products in TF32 (10-bit mantissa): SEx moves by up to 0.03 eV
        # already in iteration 1 and the MLO route carries that into both lmf runs; kt1 2026-09-26 gave
        # a MaxDiff of 1.35e-2 in fp evl.  Loosen only for --mp; the CPU and GPU checks stay tight.
        tol_log = 2e-2 if args.mp else 3e-3
        tall+=diffnum(testdir+'/'+out2, workdir+'/'+out2,tol=tol_log,comparekeys=['fp evl'])
        return tall
