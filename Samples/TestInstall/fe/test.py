from comp import test1_check, test2_check, runprogs, rmfiles
def test(args, bindir, testdir, workdir):
    M = 'fe'
    lmfa = f'mpirun -np 1 {bindir}/lmfa '
    lmf  = f'mpirun -np {args.np} {bindir}/lmf '
    tall = ''
    out1 = f'out.lmf.{M}'
    dos  = f'dos.tot.{M}'
    message = '''
    # Case fe: bcc Fe, spin-polarized LDA + total DOS.
    #  1. ferromagnetic metal, METAL=3 tetrahedron BZ integration (nkabc=10)
    #  2. self-consistency from scratch (lmfa start), conv=1e-4 Ry
    #  3. total DOS via job_tdos (dos.tot.fe)
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
