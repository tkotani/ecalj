from comp import test2_check, runprogs, rmfiles
# SiC (Materials Project SiC, DFT level): lmfa -> lmf -> job_band -> job_mlo, then compare the MLO bands.
# Every s,p,d channel of every atom is in the model (mlo_lm); mlo_delta = mlo_w = 2 eV (defaults).
def test(args, bindir, testdir, workdir):
    material = "sic"
    lmfa = f'mpirun -np 1 {bindir}/lmfa '
    lmf  = f'mpirun -np {args.np} {bindir}/lmf '
    dats = ['band_MLO_spin1.dat']
    rmfiles(workdir, ['out.lmf.' + material] + dats)
    if not args.checkonly:
        runprogs([
            lmfa + f'{material} > out.lmf.{material}',
            lmf  + f'{material} >> out.lmf.{material}',
            f'{bindir}/job_band {material} -np {args.np} --nognuplot',
            f'{bindir}/job_mlo -np {args.np} {material} --nognuplot',
        ])
    tall = ''
    for dat in dats:
        print(dat, end=': ')
        tall += test2_check(testdir + '/' + dat, workdir + '/' + dat)
    return tall
