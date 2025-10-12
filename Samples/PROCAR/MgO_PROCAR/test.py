from comp import test2_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir): #Fixed. called as >testecalj Fe_magnon
    tall=''
    MATERIAL="mgo"
    ncore=args.np
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    outfile='out.lmf.fe'
    dat= 'bw.dat'
    rmfiles(workdir,[outfile,dat])
    if(args.checkonly): runprogs([
            "rm -rf summary.txt"
            ],quiet=True)
    else: runprogs([
            lmfa +f"{MATERIAL} >llmfa",
            lmf  +f"{MATERIAL} > llmf",
            f"{bindir}/job_band {MATERIAL} -np {ncore} > ljob_band",
            "rm -rf PROCAR*",
            lmf + f"--mkprocar --band:fn=syml mgo >lbandW", # This is for pdos mode
            "cat PROCAR.UP.* >>PROCAR.UP",
            "rm PROCAR.UP.*",
            f"{workdir}/BandWeight.py > bw.dat",
            "gnuplot bnds.gnu.mgoW",
            f"evince {workdir}/mgoWeight.pdf"
    ])
    print(dat,end=': ')
    tall+=test2_check(testdir+'/'+dat, workdir+'/'+dat) 
    message1=f'''
    ======================================================
    We show fat band for O(2p) in {workdir}/mgoWeight.pdf
    ======================================================
    '''
    print(message1)
    return tall
