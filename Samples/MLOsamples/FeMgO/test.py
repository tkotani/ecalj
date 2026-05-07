from comp import test2_check, runprogs, rmfiles
import glob, os

def test(args, bindir, testdir, workdir):
    job_mlo = f'{bindir}/job_mlo '
    npflag = f'-np {args.np} '
    tall = ''
    def cleanup_mlo():
        for pat in ('HamiltonianPMT.*', '__HamiltonianPMT*', 'PROCAR.UP.*', 'PROCAR.DN.*'):
            for p in glob.glob(os.path.join(workdir, pat)):
                try: os.remove(p)
                except OSError: pass
        rmfiles(workdir, ['HamRsMLO', 'band_MLO_spin1.dat', 'band_MLO_spin2.dat',
                          'PROCAR.UP', 'PROCAR.DN', 'lwriteham', 'lmlo',
                          'bandplot_MLO.isp1.glt', 'bandplot_MLO.isp2.glt'])
    cleanup_mlo()
    runprogs([f'{job_mlo} femgo {npflag} --NoGnuplot'])
    tall += test2_check(os.path.join(testdir, 'band_MLO_spin1.dat'),
                        os.path.join(workdir, 'band_MLO_spin1.dat'),
                        abs_tol=7.4e-5)
    tall += test2_check(os.path.join(testdir, 'band_MLO_spin2.dat'),
                        os.path.join(workdir, 'band_MLO_spin2.dat'),
                        abs_tol=7.4e-5)
    return tall
