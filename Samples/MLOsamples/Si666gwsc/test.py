from comp import test2_check, runprogs, rmfiles
def test(args, bindir, testdir, workdir):
    job_band = f'{bindir}/job_band '
    npflag = f'-np {args.np} '
    bndfiles = [f'bnd00{i}.spin1' for i in range(1, 7)]
    tall = ''
    rmfiles(workdir, ['llmf_ef', 'llmf_band'] + bndfiles)
    runprogs([
        f'{job_band} si {npflag} --NoGnuplot',
    ])
    for out in bndfiles:
        print(out, end=' ')
        tall += test2_check(testdir + '/' + out, workdir + '/' + out)
    return tall
