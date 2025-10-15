from comp import test1_check,test2_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir):
    MATERIAL="ni"
    ncore=args.np
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    outfile=f'out.lmf.{MATERIAL}'
    dat1='wan_ChiPMz.mat.syml1'
    dat2='wan_ChiPMz.mat.syml2'
    dat3='wan_ChiPMr.mat.syml1'
    dat4='wan_ChiPMr.mat.syml2'
    tall=''
    
    rmfiles(workdir,[outfile,dat1,dat2,dat3,dat4])
    if(args.checkonly): runprogs([
            "rm -rf summary.txt"
            ],quiet=True)
    else: runprogs([
        lmfa + f" {MATERIAL} > "+ outfile,
        lmf  + f" {MATERIAL} > "+ outfile,
        f"{bindir}/job_band   {MATERIAL} -np {ncore}",
        f"{bindir}/genMLWF_vw {MATERIAL} -np {ncore}",              # Wannier
        "echo --- Go into epsPP_magnon. It may take several minutes ---",
        "date",
        f"{bindir}/epsPP_magnon_chipm_mpi -np {ncore} {MATERIAL}",  # magnon calculation
        "date",
        "gnuplot fbplot.glt" ,
        "gnuplot wanplot.glt",
        "gnuplot mag3d.glt",
        f"evince {workdir}/magnon3d.pdf &"
    ])
    tol=0.001
    skipcond = lambda line: len(line.split()) >= 5 and all(float(x) == 0.0 for x in line.split()[:5])
    tall+=test2_check(testdir+'/'+dat1, workdir+'/'+dat1, tol, rel_tol=1e-3, skipcond=skipcond)
    tall+=test2_check(testdir+'/'+dat2, workdir+'/'+dat2, tol, rel_tol=1e-3, skipcond=skipcond)
    tall+=test2_check(testdir+'/'+dat3, workdir+'/'+dat3, tol, rel_tol=1e-3, skipcond=skipcond)
    tall+=test2_check(testdir+'/'+dat4, workdir+'/'+dat4, tol, rel_tol=1e-3, skipcond=skipcond)
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
