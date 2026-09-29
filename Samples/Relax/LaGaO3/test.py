from comp import test2_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir):
    M = 'lagao3'
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    message1='''
    # Case Relax/LaGaO3: relaxation of the atomic positions of LaGaO3 (20 atoms), [dyn] mode=5 (Fletcher-Powell)
    #  A short run for the test: two geometries (--ctrlg:dyn.nit=2) with five iterations of the
    #  self-consistency each (--ctrlg:iter.nit=5; not converged, the lines of save.lagao3 end with x).
    #  Checked: energies and forces at the end of each geometry, and the positions of AtomPos.lagao3
    '''
    print(message1)
    if args.checkonly:   # only the comparison, in the work directory of an earlier run
        runprogs([f'python3 {testdir}/relax_result.py {M} > lrelax', 'rm -f summary.txt'], quiet=True)
    else:
        rmfiles(workdir,[f'save.{M}', f'AtomPos.{M}', 'relax_energy.txt', 'relax_force.txt'])
        runprogs([
            lmfa + f'{M} > llmfa',
            lmf  + f'{M} --ctrlg:iter.nit=5 --ctrlg:dyn.nit=2 > llmf',
            f'python3 {testdir}/relax_result.py {M} > lrelax',
        ])
    print('relax_energy.txt', end=': ')
    tall = test2_check(testdir+'/relax_energy.txt', workdir+'/relax_energy.txt', abs_tol=2e-3)   # eV
    print('relax_force.txt', end=': ')
    tall+= test2_check(testdir+'/relax_force.txt', workdir+'/relax_force.txt', abs_tol=0.05)     # mRy/bohr
    print(f'AtomPos.{M}', end=': ')
    # the reference is AtomPos.lagao3.ref: a file AtomPos.lagao3 in the run directory is read by lmf as the positions
    tall+= test2_check(testdir+f'/AtomPos.{M}.ref', workdir+f'/AtomPos.{M}', abs_tol=2e-5)       # units of alat
    return tall
