from comp import test1_check,test2_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir):
    M = 'fept'
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    message1='''
    # Case FePt_MAE: magnetic anisotropy energy (MAE) of L1_0 FePt by the force theorem
    #  1. SCF without SOC (so=0), spin-averaged radial functions (phispinsym = true in ctrlg.fept.toml)
    #  2. one-shot band energies sev without symmetry (--nosym --quit=band): so=0, so=1 axis 001, so=1 axis 110
    #  3. mae.py: MAE = sev(110) - sev(001)
    # The k mesh 6x6x4 is for a test; the MAE is not converged.
    '''
    print(message1)
    save = f'save.{M}'
    oneshot = f' {M} --nosym --quit=band '
    if args.checkonly:   # only the comparison, in the work directory of an earlier run
        runprogs([f'python3 {testdir}/mae.py > lmae', 'rm -f summary.txt'], quiet=True)
    else:
        rmfiles(workdir,[save,'energies.txt','mae.txt'])
        runprogs([
            lmfa + f'{M} > llmfa',
            lmf  + f'{M} > llmf_scf',
            lmf  + oneshot + '--ctrlg:ham.so=0 > llmf_so0',
            lmf  + oneshot + '--ctrlg:ham.so=1 --ctrlg:ham.socaxis=[0,0,1] > llmf_so001',
            lmf  + oneshot + '--ctrlg:ham.so=1 --ctrlg:ham.socaxis=[1,1,0] > llmf_so110',
            f'python3 {testdir}/mae.py > lmae',
        ])
    # save.fept: first iteration (h line) and the converged SCF (c line)
    tall = test1_check(testdir+'/'+save, workdir+'/'+save)
    # energies (eV): sev moves by about 1e-3 eV in the last SCF iterations, so 2e-3
    print('energies.txt', end=': ')
    tall+= test2_check(testdir+'/energies.txt', workdir+'/energies.txt', abs_tol=2e-3)
    # MAE (meV): a difference of two runs on the same potential
    print('mae.txt', end=': ')
    tall+= test2_check(testdir+'/mae.txt', workdir+'/mae.txt', abs_tol=0.02)
    return tall
