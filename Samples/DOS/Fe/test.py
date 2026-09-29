from comp import test2_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir):
    M = 'fe'
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    message1='''
    # Case DOS/Fe: total DOS and partial DOS of bcc Fe (ferromagnetic, GGA)
    #  1. lmfa, lmf: self-consistency
    #  2. job_tdos: total DOS by the tetrahedron method (lmf --tdos) -> dos.tot.fe, dosi.tot.fe
    #  3. job_pdos: DOS in the spheres for each lm (lmf --mkprocar --fullmesh, lmf --writepdos) -> dos.isp*.site*.fe
    #  4. dos_table.py: tables in eV, sums over m -> tdos.txt, pdos_l.txt, dos.png
    '''
    print(message1)
    table = f'python3 {testdir}/dos_table.py {M} --emin -15 --emax 15 --title "bcc Fe" --sites Fe > ldos_table'
    if args.checkonly:   # only the tables and the comparison, in the work directory of an earlier run
        runprogs([table, 'rm -f summary.txt'], quiet=True)
    else:
        rmfiles(workdir,['tdos.txt','pdos_l.txt','dos.png'])
        runprogs([
            lmfa + f'{M} > llmfa',
            lmf  + f'{M} > llmf',
            f'{bindir}/job_tdos {M} -np {args.np} --NoGnuplot > ljob_tdos',
            f'{bindir}/job_pdos {M} -np {args.np} --emin=-15 --emax=15 --ndos=1000 --NoGnuplot > ljob_pdos',
            table,
        ])
    # DOS in states/eV (peaks of 5 to 30) and numbers of states; energies in eV
    print('tdos.txt', end=': ')
    tall = test2_check(testdir+'/tdos.txt', workdir+'/tdos.txt', abs_tol=5e-3)
    print('pdos_l.txt', end=': ')
    tall+= test2_check(testdir+'/pdos_l.txt', workdir+'/pdos_l.txt', abs_tol=5e-3)
    return tall
