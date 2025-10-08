from comp import test2_check,runprogs
def test(args,bindir,testdir,workdir): #Fixed. called as >testecalj Fe_magnon
    tall=''
    MATERIAL="fe"
    ncore=args.np
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    outfile='out.lmf.fe'
    if(args.checkonly): runprogs([
            "rm -rf summary.txt"
            ],quiet=True)
    else: runprogs([
        lmfa + f" {MATERIAL} > "+ outfile,
        lmf  + f" {MATERIAL} > "+ outfile,
        f"{bindir}/job_band   {MATERIAL} -np {ncore}",              # band plot
        f"{bindir}/genMLWF_vw {MATERIAL} -np {ncore}",              # Wannier
        "echo --- Go into epsPP_magnon. It may take several minutes ---",
        "date",
        f"{bindir}/epsPP_magnon_chipm_mpi -np {ncore} {MATERIAL}",  # magnon calculation
        "date",
        "gnuplot fbplot.glt" ,
        "gnuplot wanplot.glt",
        "gnuplot mag3d.glt"
    ])
    dat='wan_ChiPMz.mat.syml1'
    print(dat,end=': ')
    tall+=test2_check(testdir+'/'+dat, workdir+'/'+dat) #numerical agreement check
    dat='wan_ChiPMr.mat.syml1'
    print(dat,end=': ')
    tall+=test2_check(testdir+'/'+dat, workdir+'/'+dat)
    message1='''
     ======================================================
     Magnon calculation finished                           
     'wan_ChiPMr.dat' <--- R(q,omega)   
     'wan_ChiPMz.dat' <--- K(q,omega)
     '*.eps' are genereted!
     Compare the results to the prepared eps file in ./eps/
     ======================================================
    '''
    print(message1)
    return tall
