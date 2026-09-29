from comp import test2_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir):
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    message1='''
    # Case ReN: rare-earth nitrides GdN (cgdn) and PrN (cprn), rocksalt
    #  LDA+U on the 4f shell (idu=2: FLL) together with spin-orbit coupling so=2 (LzSz);
    #  the 4f density matrix starts from occnum.<sname> (Hund's rule).
    #  Checked: result.<sname>.txt = total spin moment and energy (first iteration and converged),
    #           spin and orbital moments of the sites, occupations of the 4f orbitals
    '''
    print(message1)
    tall=''
    for M in ['cgdn','cprn']:
        result = f'result.{M}.txt'
        if args.checkonly:   # only the comparison, in the work directory of an earlier run
            runprogs([f'python3 {testdir}/ldau_result.py {M} > lresult.{M}'] + (['rm -f summary.txt'] if M=='cgdn' else []), quiet=True)
        else:
            rmfiles(workdir,[result, f'save.{M}', f'dmats.{M}', f'orbitalmom.{M}.chk'])
            runprogs([
                lmfa + f'{M} > llmfa.{M}',
                lmf  + f'{M} > llmf.{M}',
                f'python3 {testdir}/ldau_result.py {M} > lresult.{M}',
            ])
        # energies (eV) and moments: the self-consistency of the density matrix is slow in PrN
        # (ehk still moves by 1e-4 eV per iteration at the end), so 1e-3
        print(result, end=': ')
        tall+= test2_check(testdir+'/'+result, workdir+'/'+result, abs_tol=1e-3)
    return tall
