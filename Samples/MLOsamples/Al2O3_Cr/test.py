from comp import test2_check, runprogs, rmfiles
import glob, os

def test(args, bindir, testdir, workdir):
    job_mlo = f'{bindir}/job_mlo '
    npflag = f'-np {args.np} '
    tall = ''

    # Al2O3_Cr is QSGW80 (ssig=0.8) with rst + sigm pre-converged.
    # mlo_emax=7 eV: Cr 3d localised orbitals reach ~7 eV above E_F.
    # nspin=2, so=2 (Lz.Sz coupling baked into self-energy/density).
    def cleanup_mlo():
        for pat in ('HamiltonianPMT.*', '__HamiltonianPMT*', 'PROCAR.UP.*', 'PROCAR.DN.*'):
            for p in glob.glob(os.path.join(workdir, pat)):
                try: os.remove(p)
                except OSError: pass
        rmfiles(workdir, ['HamRsMLO', 'band_MLO_spin1.dat', 'band_MLO_spin2.dat',
                          'PROCAR.UP', 'PROCAR.DN', 'lwriteham', 'lmlo',
                          'bandplot_MLO.isp1.glt', 'bandplot_MLO.isp2.glt'])
    cleanup_mlo()
    runprogs([f'{job_mlo} al2o3_cr {npflag} --NoGnuplot'])
    tall += test2_check(os.path.join(testdir, 'band_MLO_spin1.dat'),
                        os.path.join(workdir, 'band_MLO_spin1.dat'),
                        abs_tol=7.4e-5)
    tall += test2_check(os.path.join(testdir, 'band_MLO_spin2.dat'),
                        os.path.join(workdir, 'band_MLO_spin2.dat'),
                        abs_tol=7.4e-5)
    return tall
