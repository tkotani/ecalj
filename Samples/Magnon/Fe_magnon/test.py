from comp import test2_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir): #Fixed. called as >testecalj Fe_magnon
    MATERIAL="fe"
    ncore=args.np
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    outfile='out.lmf.fe'
    dat1='TrKpm.syml001'
    dat2='TrRpm.syml001'
    tall=''

    rmfiles(workdir,[outfile,dat1,dat2])
    if(args.checkonly): runprogs([
            "rm -rf summary.txt"
            ],quiet=True)
    else: runprogs([
        lmfa + f" {MATERIAL} > "+ outfile,
        lmf  + f" {MATERIAL} > "+ outfile,
        f"{bindir}/job_band   {MATERIAL} -np {ncore}",              # band plot
        "date",
        f"{bindir}/job_magnon -np {ncore} {MATERIAL}",  # magnon calculation
        "date",
        "gnuplot mag3d.glt",
        "gnuplot r_k.glt",
        "gnuplot wan_bandplot.glt"
    ])
    print(dat1,end=': ')
    tall+=test2_check(testdir+'/'+dat1, workdir+'/'+dat1) #numerical agreement check
    print(dat2,end=': ')
    # skip condition: q=0 and ω=0 point
    skipcond = lambda line: len(line.split()) >= 5 and all(float(x) == 0.0 for x in line.split()[:5])
    tall+=test2_check(testdir+'/'+dat2, workdir+'/'+dat2, abs_tol =1e-3, rel_tol=1e-3, skipcond=skipcond)
    message1=f'''
     ======================================================
     Magnon calculation finished                           
     'TrRpm.dat' <--- R(q,omega)   
     'TrKpm.dat' <--- K(q,omega)
     '*.pdf' are genereted!
    
     Compare the prevous results ./ref/
       >evince {workdir}/magnon3d_100.pdf
       >evince {testdir}/eps/magnon3d_100.eps
     ======================================================
    '''
    print(message1)
    return tall
