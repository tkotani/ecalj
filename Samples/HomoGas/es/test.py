from comp import test2_check, runprogs, rmfiles
def test(args, bindir, testdir, workdir):
    M = 'es'
    lmfa  = f'mpirun -np 1 {bindir}/lmfa '
    lmf1  = f'mpirun -np 1 {bindir}/lmf '
    lmf   = f'mpirun -np {args.np} {bindir}/lmf '
    qg4gw = f'mpirun -np 1 {bindir}/qg4gw '
    hhomogas = f'mpirun -np 1 {bindir}/hhomogas '
    tall = ''
    dat = 'x0homo.dat'   # chi0(q,omega)*volume of the electron gas for 10 q points
    message = '''
    # Case HomoGas/es: Lindhard function of the homogeneous electron gas by the tetrahedron method
    #  1. lmfa, lmf --jobgw=0, qg4gw --job=1, lmf --jobgw=1 --skipCPHI: lattice, symmetry and the k mesh
    #     [gw] n1n2n3 with its tetrahedra (no self-consistent run is needed)
    #  2. hhomogas: tetrahedron weights for the free-electron band, Im chi0, and Re chi0 by the
    #     Hilbert transformation -> x0homo.dat
    '''
    print(message)
    if args.checkonly:
        runprogs(["rm -rf summary.txt"], quiet=True)
    else:
        rmfiles(workdir, [dat])
        runprogs([
            lmfa  + f"{M} > llmfa",
            lmf1  + f"{M} --jobgw=0 > llmfgw00",
            qg4gw + f"{M} --job=1 > lqg4gw",
            lmf   + f"{M} --jobgw=1 --skipCPHI > llmfgw01",
            hhomogas + f"{M} > lhhomogas",
        ])
    print(dat, end=': ')
    tall += test2_check(testdir + '/' + dat, workdir + '/' + dat, abs_tol=1e-3, rel_tol=1e-3)
    print(f'''
    ==========================================================================
    Figure: cd {workdir}; python3 plot_lindhard.py   (lindhard.png)
    ==========================================================================
    ''')
    return tall
