from comp import test1_check,test2_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir):
    MATERIAL="ni"
    ncore=args.np
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    outfile=f'out.lmf.{MATERIAL}'
    dat1='TrKpm.syml001'
    dat2='TrRpm.syml001'
    dat3='TrRpm.syml002'
    dat4='TrRpm.syml002'
    tall=''
    
    rmfiles(workdir,[outfile,dat1,dat2,dat3,dat4])
    if(args.checkonly): runprogs([
            "rm -rf summary.txt"
            ],quiet=True)
    else: runprogs([
        lmfa + f" {MATERIAL} > "+ outfile,
        lmf  + f" {MATERIAL} > "+ outfile,
        f"{bindir}/job_band   {MATERIAL} -np {ncore}",
        "date",
        f"{bindir}/job_magnon -np {ncore} {MATERIAL}",  # magnon calculation
        "date",
        "gnuplot mag3d.glt",
        "gnuplot r_k.glt",
        "gnuplot wan_bandplot.glt"
    ])
    tol=0.001
    skipcond = lambda line: len(line.split()) >= 5 and all(float(x) == 0.0 for x in line.split()[:4])
    tall+=test2_check(testdir+'/'+dat1, workdir+'/'+dat1, tol, rel_tol=1e-3, skipcond=skipcond)
    tall+=test2_check(testdir+'/'+dat2, workdir+'/'+dat2, tol, rel_tol=1e-3, skipcond=skipcond)
    tall+=test2_check(testdir+'/'+dat3, workdir+'/'+dat3, tol, rel_tol=1e-3, skipcond=skipcond)
    tall+=test2_check(testdir+'/'+dat4, workdir+'/'+dat4, tol, rel_tol=1e-3, skipcond=skipcond)
    message1='''
     ======================================================
     Magnon calculation finished                           
     'TrRpm.dat' <--- R(q,omega)   
     'TrKpm.dat' <--- K(q,omega)
     '*.pdf' are genereted!
     Compare the results to the prepared eps file in ./ref/
     ======================================================
    '''
    print(message1)
    return tall
