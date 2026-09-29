from comp import test2_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir):
    M = 'gaas'
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    message1='''
    # Case EffectiveMass/GaAs: effective masses of GaAs from the QSGW bands with spin-orbit coupling
    #  1. self-consistency with the stored QSGW self-energy sigm.gaas (nspin=1, so=0)
    #  2. bands with SOC (so=1) along three short lines from Gamma: syml.gaas in the mass mode
    #  3. massfit.py: fit of k^2 = a|E| + b E^2 for each band (numpy) -> mass.txt
    '''
    print(message1)
    if args.checkonly:   # only the fit and the comparison, in the work directory of an earlier run
        runprogs([f'python3 {testdir}/massfit.py > lmassfit', 'rm -f summary.txt'], quiet=True)
    else:
        rmfiles(workdir,['mass.txt','mass_detail.txt','massfit.npz','massfit.png'])
        runprogs([
            lmfa + f'{M} > llmfa',
            lmf  + f'{M} > llmf',
            f'{bindir}/job_band {M} -np {args.np} --ctrlg:ham.nspin=2 --ctrlg:ham.so=1 --NoGnuplot > ljob_band',
            f'python3 {testdir}/massfit.py > lmassfit',
        ])
    # masses (in units of the electron mass), gap and spin-orbit splitting at Gamma (eV)
    print('mass.txt', end=': ')
    tall = test2_check(testdir+'/mass.txt', workdir+'/mass.txt', abs_tol=1e-3)
    return tall
