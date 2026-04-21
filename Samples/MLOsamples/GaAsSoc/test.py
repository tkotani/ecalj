from comp import test2_check, runprogs, rmfiles
import glob, os

def test(args, bindir, testdir, workdir):
    job_mlo     = f'{bindir}/job_mlo '
    job_mlo_soc = f'{bindir}/job_mlo_soc '
    npflag = f'-np {args.np} '
    tall = ''

    def cleanup():
        for pat in ('HamiltonianPMT.*', '__HamiltonianPMT*', 'PROCAR.UP.*', 'PROCAR.DN.*'):
            for p in glob.glob(os.path.join(workdir, pat)):
                try: os.remove(p)
                except OSError: pass
        rmfiles(workdir, ['HamRsMLO', 'band_MLO_spin1.dat', 'band_MLO_spin2.dat',
                          'PROCAR.UP', 'PROCAR.DN', 'lwriteham', 'lmlo', 'llmf_ef_soc',
                          'bandplot_MLO.isp1.glt', 'bandplot_MLO.isp2.glt'])

    # --- Test 1: non-SOC MLO bands ---
    cleanup()
    runprogs([f'{job_mlo} gaas {npflag} --NoGnuplot'])
    tall += test2_check(os.path.join(testdir, 'band_MLO_spin1.dat'),
                        os.path.join(workdir, 'band_MLO_spin1.dat'),
                        abs_tol=7.4e-5)  # 0.001 eV in Ry

    # --- Test 2: SOC MLO bands via job_mlo_soc ---
    cleanup()
    runprogs([f'{job_mlo_soc} gaas {npflag} --NoGnuplot'])
    tall += test2_check(os.path.join(testdir, 'band_MLO_spin1.soc.dat'),
                        os.path.join(workdir, 'band_MLO_spin1.dat'),
                        abs_tol=7.4e-5)

    print()
    print('=' * 70)
    print(f'To view the SOC band plot (last run is job_mlo_soc, 2N=36 spinor):')
    print(f'  cd {workdir}')
    print(f'  gnuplot -p bandplot_MLO.isp1.glt')
    print('  (red points = MLO-SOC bands, black lines = non-SOC DFT bands)')
    print(f'SOC Fermi energy saved in: {workdir}/efermi_soc')
    print('=' * 70)
    return tall
