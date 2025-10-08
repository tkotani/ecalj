from comp import test2_check,runprogs
def test(args,bindir,testdir,workdir): #Fixed. called as >testecalj Fe_magnon
    tall=''
    MATERIAL="cu"
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
        f"{bindir}/getsyml   {MATERIAL}",
        f"{bindir}/job_band  {MATERIAL} -np {ncore} ",        
        f"{bindir}/epsPP0 {MATERIAL} -np {ncore} ",            
        "gnuplot -p epsinter.glt" ,
        "gnuplot -p epsintra.glt",
        "gnuplot -p epsall.glt"
    ])
    import re
    skipcond = lambda line: any( abs(float(tok)) >= 1e4
                                 for tok in re.split(r'\s+', line.strip())[4:]
                                 if tok != '' )
    for file in ['EPS0001','EPS0002','EPS0003']:
        dat= file+'.nlfc.dat.interbandonly'
        print(dat,end=': ')
        tall+=test2_check(testdir+'/'+dat, workdir+'/'+dat) 
        dat= file+'.nlfc.dat.intrabandonly'
        print(dat,end=': ')
        tall+=test2_check(testdir+'/'+dat, workdir+'/'+dat,skipcond=skipcond) 
    message1='''
     ======================================================
     See plots
        "gnuplot -p epsinter.glt"  interband
        "gnuplot -p epsintra.glt"  intraband
        "gnuplot -p epsall.glt"    all
     ======================================================
    '''
    print(message1)
    return tall
