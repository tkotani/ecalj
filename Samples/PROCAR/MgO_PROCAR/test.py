from comp import test2_check, runprogs, rmfiles
def test(args, bindir, testdir, workdir):
    MATERIAL = "mgo"
    ncore = args.np
    lmfa = f'mpirun -np 1 {bindir}/lmfa '
    lmf  = f'mpirun -np {args.np} {bindir}/lmf '
    outfile = f'out.lmf.{MATERIAL}'
    dat = 'bw.dat'
    tall = ''
    rmfiles(workdir, [outfile, dat])
    if args.checkonly:
        runprogs(["rm -rf summary.txt"], quiet=True)
    else:
        runprogs([
            lmfa + f"{MATERIAL} > llmfa",
            lmf  + f"{MATERIAL} > llmf",
            f"{bindir}/job_band {MATERIAL} -np {ncore} --NoGnuplot > ljob_band",
            "rm -rf PROCAR*",
            lmf + f"--mkprocar --band:fn=syml {MATERIAL} > lbandW",
            "cat PROCAR.UP.* >> PROCAR.UP",
            "rm PROCAR.UP.*",
            f"{workdir}/BandWeight.py > bw.dat",
            "gnuplot bnds.gnu.mgoW",
        ])
    print(dat, end=': ')
    tall += test2_check(testdir + '/' + dat, workdir + '/' + dat)
    print(f'''
    ==========================================================================
    Fat band PDF: {workdir}/mgoWeight.pdf (O-2p weight)
    To view: evince {workdir}/mgoWeight.pdf
    ==========================================================================
    ''')
    return tall
