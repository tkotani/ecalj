from comp import test2_check,runprogs
def test(args,bindir,testdir,workdir): #Fixed. called as >testecalj Fe_magnon
    tall=''
    MATERIAL="mgo"
    ncore=args.np
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    outfile='out.lmf.fe'
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
            "rm PROCAR.UP.* bw.dat",
            f"{workdir}/BandWeight.py > bw.dat",
            "gnuplot bnds.gnu.mgoW"
    ])
    dat= 'bw.dat'
    print(dat,end=': ')
    tall+=test2_check(testdir+'/'+dat, workdir+'/'+dat) 
    message1=f'''
     ======================================================
     See {workdir}/mgoWeight.pdf file for fat band
     ======================================================
    '''
    print(message1)
    return tall
