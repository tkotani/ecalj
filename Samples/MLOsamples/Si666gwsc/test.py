from comp import test2_check, runprogs, rmfiles
import glob, os

def test(args, bindir, testdir, workdir):
    job_band = f'{bindir}/job_band '
    job_mlo  = f'{bindir}/job_mlo '
    npflag = f'-np {args.np} '
    bndfiles = [f'bnd00{i}.spin1' for i in range(1, 7)]
    tall = ''

    # --- Test 1: DFT bands via job_band ---
    rmfiles(workdir, ['llmf_ef', 'llmf_band'] + bndfiles)
    runprogs([f'{job_band} si {npflag} --NoGnuplot'])
    for out in bndfiles:
        print(out, end=' ')
        tall += test2_check(testdir + '/' + out, workdir + '/' + out)

    # --- Test 2: MLO bands via job_mlo (non-magnetic Si, single spin) ---
    def cleanup_mlo():
        for pat in ('HamiltonianPMT.*', '__HamiltonianPMT*', 'PROCAR.UP.*', 'PROCAR.DN.*'):
            for p in glob.glob(os.path.join(workdir, pat)):
                try: os.remove(p)
                except OSError: pass
        rmfiles(workdir, ['HamRsMLO', 'band_MLO_spin1.dat',
                          'PROCAR.UP', 'PROCAR.DN', 'lwriteham', 'lmlo',
                          'bandplot_MLO.isp1.glt'])
    cleanup_mlo()
    runprogs([f'{job_mlo} si {npflag} --NoGnuplot'])
    tall += test2_check(os.path.join(testdir, 'band_MLO_spin1.dat'),
                        os.path.join(workdir, 'band_MLO_spin1.dat'),
                        abs_tol=7.4e-5)  # 0.001 eV in Ry

    return tall
