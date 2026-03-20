from comp import test2_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir): #Fixed. called as >testecalj Fe_magnon
    MATERIAL="fe"
    ncore=args.np
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    outfile='out.lmf.fe'
    dats=['band_MLO_spin1.dat', 'band_MLO_spin2.dat']
    tall=''

    rmfiles(workdir,[outfile]+dats)
    if(args.checkonly): runprogs([
            "rm -rf summary.txt"
            ],quiet=True)
    else: runprogs([
        lmfa + f" {MATERIAL} > "+ outfile,
        lmf  + f" {MATERIAL} > "+ outfile,
        f"{bindir}/job_band {MATERIAL} -np {ncore} --nognuplot",
        f"{bindir}/job_mlo -np {ncore} {MATERIAL} --nognuplot",
    ])
    for dat in dats:
        print(dat,end=': ')
        tall+=test2_check(testdir+'/'+dat, workdir+'/'+dat)
    return tall
