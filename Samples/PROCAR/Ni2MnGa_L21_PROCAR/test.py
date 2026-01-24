from comp import test2_check,runprogs, rmfiles
def test(args,bindir,testdir,workdir): #Fixed. called as >testecalj Fe_magnon
    MATERIAL="ni2mnga"
    ncore=args.np
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    dats = ['bandweight_atom1.spin1', 'bandweight_atom2.spin1', 'bandweight_atom1.spin2', 'bandweight_atom2.spin2']
    tall=''

    rmfiles(workdir, dats)
    if(args.checkonly): runprogs([
            "rm -rf summary.txt"
            ],quiet=True)
    runprogs([
            lmfa +f"{MATERIAL} >llmfa",
            lmf  +f"{MATERIAL} > llmf",
            f"{bindir}/job_band {MATERIAL} -np {ncore} --fatband --emin=-5 --emax=5 > ljob_band",
            "gnuplot fatband.glt",
            f"evince {workdir}/fatband.pdf &"
    ])
    for dat in dats:
        print(dat, end=': ')
        tall+=test2_check(testdir+'/'+dat, workdir+'/'+dat, abs_tol=0.0001) #numerical agreement check
    return tall
