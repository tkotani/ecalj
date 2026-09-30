from comp import test2_check,peak_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir):
    M = 'fe'
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    message1='''
    # Case Fe_mlo_magnon: magnon spectrum of bcc Fe from the MLO model (job_mlo_magnon), the MLO counterpart of Fe_magnon
    #  1. lmfa, lmf (LDA, spin polarized), job_band (the bands that the MLO fit uses)
    #  2. job_mlo_magnon: the MLO Hamiltonian on mlo_nkabc, W in the MLO basis, the overlaps at k and k+q, chi^+-
    #  Checked: MagSuscep.syml001 = K(q,omega) and R(q,omega) along syml.fe (values, and the positions of the peaks of Im R)
    '''
    print(message1)
    dat = 'MagSuscep.syml001'
    if args.checkonly:
        runprogs(['rm -f summary.txt'], quiet=True)
    else:
        rmfiles(workdir,[dat, 'MagSpec.syml001', 'HamRsMLO', 'HamiltonianPMTInfo'])
        runprogs([
            lmfa + f'{M} > llmfa',
            lmf  + f'{M} > llmf',
            f'{bindir}/job_band {M} -np {args.np} --NoGnuplot > ljob_band',
            f'{bindir}/job_mlo_magnon {M} -np {args.np} > ljob_mlo_magnon',
        ])
    # R(q,omega) has sharp poles: a shift of a pole by one frequency step fails a comparison point by point although the
    # spectra agree, so the values are compared loosely and the positions of the peaks of Im R (column 9) within 5 %
    print(dat, end=': ')
    skip_w0 = lambda line: len(line.split()) >= 5 and abs(float(line.split()[4])) < 1e-6
    tall = test2_check(testdir+'/'+dat, workdir+'/'+dat, abs_tol=0.5, rel_tol=0.1, skipcond=skip_w0)
    tall+= peak_check(testdir+'/'+dat, workdir+'/'+dat, qcol=3, wcol=4, vcol=8, wmin=1e-6, wmax=1.6, rel_tol=0.05, abs_tol=0.003)
    return tall
