from comp import test2_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir): #Fixed. called as >testecalj Fe_magnon
    MATERIAL="fe"
    ncore=args.np
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    outfile='out.lmf.fe'
    dat1='wan_ChiPMz.mat.syml1'
    dat2='wan_ChiPMr.mat.syml1'
    tall=''

    rmfiles(workdir,[outfile,dat1,dat2])
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
    print(dat1,end=': ')
    tall+=test2_check(testdir+'/'+dat1, workdir+'/'+dat1) #numerical agreement check
    print(dat2,end=': ')
    tall+=test2_check(testdir+'/'+dat2, workdir+'/'+dat2)
    message1=f'''
     ======================================================
     Magnon calculation finished                           
     'wan_ChiPMr.dat' <--- R(q,omega)   
     'wan_ChiPMz.dat' <--- K(q,omega)
     '*.eps' are genereted!
    
     Compare the results pdf ./pdf/
       >evince {workdir}/magnon3d_100.pdf
       >evince {testdir}/pdf/magnon3d_100.pdf
     ======================================================
    '''
    print(message1)
    return tall
