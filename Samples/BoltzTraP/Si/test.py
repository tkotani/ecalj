from comp import test1_check, test2_check, diffnum, runprogs, rmfiles
def test(args, bindir, testdir, workdir):
    M = 'si'
    lmfa = f'mpirun -np 1 {bindir}/lmfa '
    lmf  = f'mpirun -np {args.np} {bindir}/lmf '
    tall = ''
    out1 = f'save.{M}'
    intrans = f'{M}.intrans_template.boltztrap'   # Ef (Ry) and the number of electrons
    energy  = f'{M}.energy.isp1.boltztrap'        # eigenvalues (Ry) at the irreducible k points
    struct  = f'{M}.struct.boltztrap'             # lattice vectors and the 48 symmetry operations
    message = '''
    # Case BoltzTraP/Si: input files of BoltzTraP (generic format) written by lmf --boltztrap
    #  1. lmfa, lmf: LDA self-consistency (nkabc=8)
    #  2. lmf --boltztrap on the 12x12x12 mesh
    '''
    print(message)
    if args.checkonly:
        runprogs(["rm -rf summary.txt"], quiet=True)
    else:
        rmfiles(workdir, [out1, intrans, energy, struct])
        runprogs([
            lmfa + f"{M} > llmfa",
            lmf  + f"{M} > llmf",
            lmf  + f"{M} --boltztrap '--ctrlg:bz.nkabc=[12,12,12]' > llmf_boltztrap",
        ])
    tall += test1_check(testdir + '/' + out1, workdir + '/' + out1)
    print(intrans, end=': ')
    tall += test2_check(testdir + '/' + intrans, workdir + '/' + intrans, abs_tol=2e-4)
    # diffnum, not test2_check: the eigenvalues are written as 0.22D+00, which test2_check (compall) skips
    tall += diffnum(testdir + '/' + energy, workdir + '/' + energy, tol=2e-4, comparekeys=[])
    print(struct, end=': ')
    tall += test2_check(testdir + '/' + struct, workdir + '/' + struct, abs_tol=1e-6)
    return tall
