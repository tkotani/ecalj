from comp import test1_check,runprogs
def test(args,bindir,testdir,workdir):
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    tall=''
    message1='''
    # Case C: test of homogeneous background
    # zbak is overridden via --ctrlg:bz.zbak=N (N=0 neutral, N=1 ionized).
'''
    print(message1)
    outfile='out.lmf.neutral.c'
    runprogs([
        lmfa+" c --ctrlg:bz.zbak=0 > "+outfile,
        lmf+ " c --ctrlg:bz.zbak=0 >>"+outfile,
        "rm -f *mixm.* rst.* save.* log.* *hssn.* wkp.* bsmv.* bnds.*"
    ])
    tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    outfile='out.lmf.ionized.c'
    runprogs([
        lmfa+" c --ctrlg:bz.zbak=1 > "+outfile,
        lmf+ " c --ctrlg:bz.zbak=1 >>"+outfile,
        "rm -f *mixm.* rst.* save.* log.* *hssn.* wkp.* bsmv.* bnds.*"
    ])
    tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    return tall
