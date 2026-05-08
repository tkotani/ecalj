from comp import test2_check, runprogs, rmfiles
def test(args, bindir, testdir, workdir):
    MATERIAL = "ni2mnga"
    ncore = args.np
    lmfa = f'mpirun -np 1 {bindir}/lmfa '
    lmf  = f'mpirun -np {args.np} {bindir}/lmf '
    dats = ['bandweight_atom1.spin1', 'bandweight_atom2.spin1',
            'bandweight_atom1.spin2', 'bandweight_atom2.spin2']
    tall = ''
    rmfiles(workdir, dats)
    if args.checkonly:
        runprogs(["rm -rf summary.txt"], quiet=True)
    else:
        # rst.ni2mnga is committed in source dir and copied into workdir by testecalj.
        # lmfa + lmf converge in ~1 iter from the saved rst, then job_band --fatband.
        runprogs([
            lmfa + f"{MATERIAL} > llmfa",
            lmf  + f"{MATERIAL} > llmf",
            f"{bindir}/job_band {MATERIAL} -np {ncore} --fatband --emin=-5 --emax=5 --NoGnuplot > ljob_band",
            "gnuplot fatband.glt",
        ])
    for dat in dats:
        print(dat, end=': ')
        tall += test2_check(testdir + '/' + dat, workdir + '/' + dat, abs_tol=0.0001)
    print(f'''
    ==========================================================================
    Fat band PDF: {workdir}/fatband.pdf
    To view: evince {workdir}/fatband.pdf
    ==========================================================================
    ''')
    return tall
