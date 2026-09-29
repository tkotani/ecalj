from comp import test2_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir):
    M = 'fe'
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    message1='''
    # Case AHC/Fe: anomalous Hall conductivity (AHC) of bcc Fe, magnetization along z
    #  1. SCF with the spin-orbit coupling LzSz (ham.so=2, ham.phispinsym=true in ctrlg.fe.toml)
    #  2. job_band (job_AHC asks for bnd001.spin1), job_AHC: <u_k|u_k+b> on the mesh gw.n1n2n3 (huumat_MPI --ahc)
    #     and the AHC at the Fermi energy (hahc)
    #  3. hx0ahc.py: AHC for shifts of the Fermi energy from -1 to 1 eV, 11 points, one after another
    #  4. ahc_table.py: sigma_xy of the two spin channels and their sum -> ahc.txt
    # The mesh 4x4x4 of the AHC is for a test; the AHC is far from converged.
    '''
    print(message1)
    table = f'python3 {testdir}/ahc_table.py {M} > lahc_table'
    if args.checkonly:   # only the table and the comparison, in the work directory of an earlier run
        runprogs([table, 'rm -f summary.txt'], quiet=True)
    else:
        rmfiles(workdir,['ahc.txt','scf.txt'])
        runprogs([
            lmfa + f'{M} > llmfa',
            lmf  + f'{M} > llmf',
            f'{bindir}/job_band {M} -np {args.np} --NoGnuplot > ljob_band',
            f'{bindir}/job_AHC {M} -np {args.np} > ljob_AHC',
            # without mpirun: hx0ahc.py divides the shifts over the ranks only when mpi4py and mpirun are of the same MPI
            f'python3 {bindir}/hx0ahc.py -1. 1. 11 > lhx0ahc',
            table,
        ])
    # magnetic moment (1e-3 mu_B) and total energies (eV) of the SCF
    print('scf.txt', end=': ')
    tall = test2_check(testdir+'/scf.txt', workdir+'/scf.txt', abs_tol=1e-3)
    # sigma_xy (Ohm^-1 cm^-1), values of 20 to 500. As a function of the shift it is a staircase on this mesh:
    # 0.5 allows for the change of the matrix elements, not for a level that crosses the Fermi energy.
    print('ahc.txt', end=': ')
    tall+= test2_check(testdir+'/ahc.txt', workdir+'/ahc.txt', abs_tol=0.5)
    return tall
