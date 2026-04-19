from comp import test2_check, runprogs, rmfiles
def test(args, bindir, testdir, workdir):
    job_mlo = f'{bindir}/job_mlo '
    npflag = f'-np {args.np} '
    out = 'band_MLO_spin1.dat'
    tall = ''
    rmfiles(workdir, ['lwriteham', 'lmlo', 'HamRsMLO', out])
    runprogs([
        'rm -f HamiltonianPMT.* PROCAR.UP.* PROCAR.DN.* PROCAR.UP PROCAR.DN __HamiltonianPMT*',
        f'{job_mlo} gaas {npflag} --NoGnuplot',
    ])
    tall += test2_check(testdir + '/' + out, workdir + '/' + out, abs_tol=7.4e-5)  # 0.001 eV in Ry
    return tall
