from comp import test1_check, test2_check, runprogs, rmfiles
def test(args, bindir, testdir, workdir):
    M = 'cu'
    lmfa = f'mpirun -np 1 {bindir}/lmfa '
    lmf  = f'mpirun -np {args.np} {bindir}/lmf '
    tall = ''
    out1 = f'save.{M}'       # total energies of the iterations (h, i, c lines)
    bxsf = 'fermiup.bxsf'    # eigenvalues (Ry) on the 11x11x11 grid for xcrysden
    message = '''
    # Case FermiSurface/Cu: Fermi surface of fcc Cu for xcrysden.
    #  1. lmfa, lmf: LDA self-consistency (nkabc=12, tetrahedron)
    #  2. job_fermisurface: Ef on the 10x10x10 mesh (lmf --quit=band), then the eigenvalues on
    #     the 11x11x11 grid that includes both ends (lmf --fermisurface) -> fermiup.bxsf
    '''
    print(message)
    if args.checkonly:
        runprogs(["rm -rf summary.txt"], quiet=True)
    else:
        rmfiles(workdir, [out1, bxsf])
        runprogs([
            lmfa + f"{M} > llmfa",
            lmf  + f"{M} > llmf",
            f"{bindir}/job_fermisurface {M} -np {args.np} -vnk1=10 -vnk2=10 -vnk3=10 > ljob_fermisurface",
        ])
    tall += test1_check(testdir + '/' + out1, workdir + '/' + out1)
    print(bxsf, end=': ')
    tall += test2_check(testdir + '/' + bxsf, workdir + '/' + bxsf, abs_tol=2e-4)
    print(f'''
    ==========================================================================
    Fermi surface: xcrysden --bxsf {workdir}/fermiup.bxsf
    ==========================================================================
    ''')
    return tall
