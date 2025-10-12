from comp import test1_check,test2_check,runprogs
def test(args,bindir,testdir,workdir):
    MATERIAL="ni"
    NSLOTS=args.np
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    outfile=f'out.lmf.{MATERIAL}'
    dat1='wan_ChiPMz.mat.syml1'
    dat2='wan_ChiPMz.mat.syml2'
    dat3='wan_ChiPMr.mat.syml1'
    dat4='wan_ChiPMr.mat.syml2'
    tall=''
    
    rmfiles(workdir,[outfile,dat1,dat2])
    if(args.checkonly): runprogs([
            "rm -rf summary.txt"
            ],quiet=True)
    else: runprogs([
        lmfa + f" {MATERIAL} > "+ outfile,
        lmf  + f" {MATERIAL} > "+ outfile,
        f"{bindir}/job_band   {MATERIAL} -np {NSLOTS}",
        f"{bindir}/genMLWF_vw {MATERIAL} -np {NSLOTS}",              # Wannier
        "echo --- Go into epsPP_magnon. It may take several minutes ---",
        "date",
        f"{bindir}/epsPP_magnon_chipm_mpi -np {NSLOTS} {MATERIAL}",  # magnon calculation
        "date",
        "gnuplot fbplot.glt" ,
        "gnuplot wanplot.glt",
            "gnuplot mag3d.glt",
        f"evince {workdir}/magnon3d.pdf &"
    ])
    tol=0.001
    tall+=test2_check(testdir+'/'+dat1, workdir+'/'+dat1,tol)
    tall+=test2_check(testdir+'/'+dat2, workdir+'/'+dat2,tol)
    tall+=test2_check(testdir+'/'+dat3, workdir+'/'+dat3,tol)
    tall+=test2_check(testdir+'/'+dat4, workdir+'/'+dat4,tol)
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
