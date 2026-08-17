from comp import test1_check, test2_check, runprogs, rmfiles
def test(args, bindir, testdir, workdir):
    M = 'gdn'
    lmfa = f'mpirun -np 1 {bindir}/lmfa '
    lmf  = f'mpirun -np {args.np} {bindir}/lmf '
    tall = ''
    out1 = f'out.lmf.{M}'
    dos  = f'dos.tot.{M}'
    message = '''
    # Case gdn: GdN (rocksalt), LDA+U on Gd 4f + total DOS.
    #  1. LDA+U (idu=[0,0,0,2], U=0.515 Ry on 4f), mmom=7 (Gd f^7)
    #  2. short fixed run (nit=2, by design) starting from lmfa
    #  3. total DOS via job_tdos (dos.tot.gdn)
    '''
    print(message)
    rmfiles(workdir, [out1, dos])
    if args.checkonly:
        runprogs(["rm -rf summary.txt"], quiet=True)
    else:
        runprogs([
            lmfa + f"{M} > llmfa",
            lmf  + f"{M} > " + out1,
            f"{bindir}/job_tdos {M} -np {args.np} > ljob_tdos",
        ])
    tall += test1_check(testdir + '/' + out1, workdir + '/' + out1)
    print(dos, end=': ')
    tall += test2_check(testdir + '/' + dos, workdir + '/' + dos)
    return tall
