from comp import test1_check, diffnum, runprogs, rmfiles
def test(args, bindir, testdir, workdir):
    M = 'si'
    lmfa = f'mpirun -np 1 {bindir}/lmfa '
    lmf  = f'mpirun -np {args.np} {bindir}/lmf '
    tall = ''
    zbakonly = "'--ctrlg:spec.1.z=14' '--ctrlg:spec.1.q=[2.0,2.0]' '--ctrlg:bz.zbak=0.2'"
    cases = [  # name, options given to both lmfa and lmf
        ('caseA', ""),                               # Z=14.2 (ctrlg.si.toml as it is): 8.4 valence electrons
        ('caseB', "'--ctrlg:bz.zbak=0.3'"),          # Z=14.2 and background 0.3       : 8.1
        ('caseC', zbakonly),                         # Z=14 and background 0.2         : 7.8
        ('caseD', zbakonly + " '--ctrlg:ham.nspin=2' '--ctrlg:bz.fsmom=0.2'"),  # case C with the fixed spin moment 0.2
    ]
    message = '''
    # Case Doping/Si: fractional number of electrons in Si (2 atoms per cell)
    #  caseA: fractional nuclear charge Z=14.2 with Q=2,2.2 (written in ctrlg.si.toml)
    #  caseB: Z=14.2 and the homogeneous background charge bz.zbak=0.3
    #  caseC: Z=14 and bz.zbak=0.2
    #  caseD: caseC, spin polarized, with the fixed spin moment bz.fsmom=0.2
    # Checked: total energies (save.si.<case>), and Z, the valence charge, the background
    #          charge and the Fermi energy of every iteration (log.si.<case>)
    '''
    print(message)
    if args.checkonly:
        runprogs(["rm -rf summary.txt"], quiet=True)
    else:
        rmfiles(workdir, [f'save.{M}.{c}' for c, _ in cases] + [f'log.{M}.{c}' for c, _ in cases])
        for c, opt in cases:
            runprogs([
                f"rm -f rst.{M} __mixm.{M} mixm.{M} save.{M} log.{M}",   # every case starts from the free atoms
                lmfa + f"{M} {opt} > llmfa.{c}",
                lmf  + f"{M} {opt} > llmf.{c}",
                f"mv save.{M} save.{M}.{c}",
                f"mv log.{M} log.{M}.{c}",
            ])
    for c, _ in cases:
        tall += test1_check(testdir + f'/save.{M}.{c}', workdir + f'/save.{M}.{c}')
        # '===' (the START lines, no number) is there because diffnum does not compare the first line it picks:
        # without it the line 'fa Z' of lmfa, which shows the fractional Z, would not be compared
        tall += diffnum(testdir + f'/log.{M}.{c}', workdir + f'/log.{M}.{c}', tol=1e-4,
                        comparekeys=['===', 'fa Z', 'fp qvl', 'bzmet'])
    return tall
