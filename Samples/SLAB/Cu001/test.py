from comp import test1_check, test2_check, runprogs, rmfiles
def test(args, bindir, testdir, workdir):
    M = 'cu001'
    lmfa = f'mpirun -np 1 {bindir}/lmfa '
    lmf  = f'mpirun -np {args.np} {bindir}/lmf '
    tall = ''
    out1 = f'save.{M}'               # total energies of the iterations
    wf   = f'workfunction.{M}.txt'   # vacuum level, Fermi energy, work function (eV)
    bnd  = 'bnd001.spin1'            # bands along Gamma-X of the surface Brillouin zone (eV from Ef)
    dos  = f'dos.tot.{M}'            # total DOS
    message = '''
    # Case SLAB/Cu001: Cu(001) slab of 4 layers (4 atoms in the 1x1 surface cell) with vacuum, LDA
    #  1. lmfa, lmf: self-consistency (nkabc = 6 6 1)
    #  2. workfunction.py: vacuum level of the electrostatic potential - Fermi energy
    #  3. job_band along Gamma-X-M-Gamma of the surface Brillouin zone (syml.cu001), job_tdos
    '''
    print(message)
    if args.checkonly:
        runprogs(["rm -rf summary.txt"], quiet=True)
    else:
        rmfiles(workdir, [out1, wf, bnd, dos])
        runprogs([
            lmfa + f"{M} > llmfa",
            lmf  + f"{M} > llmf",
            f"python3 {workdir}/workfunction.py {M} > lworkfunction",
            f"{bindir}/job_band {M} -np {args.np} --NoGnuplot > ljob_band",
            f"{bindir}/job_tdos {M} -np {args.np} --NoGnuplot > ljob_tdos",
        ])
    tall += test1_check(testdir + '/' + out1, workdir + '/' + out1)
    print(wf, end=': ')
    tall += test2_check(testdir + '/' + wf, workdir + '/' + wf, abs_tol=2e-3)
    print(bnd, end=': ')
    tall += test2_check(testdir + '/' + bnd, workdir + '/' + bnd, abs_tol=2e-3)
    print(dos, end=': ')
    tall += test2_check(testdir + '/' + dos, workdir + '/' + dos, abs_tol=2e-2)
    print(f'''
    ==========================================================================
    Figure: cd {workdir}; python3 plot_slab.py   (slab_cu001.png)
    ==========================================================================
    ''')
    return tall
