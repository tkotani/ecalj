from comp import test2_check,dqpu,runprogs,rmfiles
def test(args,bindir,testdir,workdir):
    M = 'c'
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    message1='''
    # Case IIR/C: impact ionization rate of diamond = width 2 Z Im Sigma of the states in one-shot GW
    #  1. self-consistency with the stored QSGW80 self-energy sigm.c (ham.scaledsigma = 0.8)
    #  2. gw_lmfh: diagonal Sigma with its imaginary part at all the irreducible points of the mesh gw.n1n2n3 = 6x6x6,
    #     for the bands up to gw.EMAXforGW = 30 eV above the Fermi energy -> QPU, QPU_life
    #  3. iir_rate.py: rate = |FWHM|/hbar -> iir.txt, iir.png
    # The mesh 6x6x6 is for a test; the rates are not converged.
    '''
    print(message1)
    rate = f'python3 {testdir}/iir_rate.py --title "diamond, 6x6x6" > liir_rate'
    if args.checkonly:   # only the comparison, in the work directory of an earlier run
        runprogs([rate, 'rm -f summary.txt'], quiet=True)
    else:
        rmfiles(workdir,['QPU','QPU_life','iir.txt','iir.png'])
        runprogs([
            lmfa + f'{M} > llmfa',
            lmf  + f'{M} > llmf',
            f'{bindir}/gw_lmfh {M} -np {args.np} > lgw_lmfh',
            rate,
        ])
    # QPU: self-energies and quasiparticle energies (eV), tolerance 1.1e-2 of dqpu
    tall = dqpu(testdir+'/QPU', workdir+'/QPU')
    # QPU_life: energy (eV, three decimals: a change of the last digit is 1e-3) and FWHM = 2 Z Im Sigma (eV, up to 1.4)
    print('QPU_life', end=': ')
    tall+= test2_check(testdir+'/QPU_life', workdir+'/QPU_life', abs_tol=2e-3)
    return tall
