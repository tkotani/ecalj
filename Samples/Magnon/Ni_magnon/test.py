from comp import test1_check,test2_check,runprogs
def test(args,bindir,testdir,workdir):
    tall=''
    MATERIAL="ni"
    NSLOTS=args.np
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    outfile=f'out.lmf.{MATERIAL}'
    runprogs([
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
        "gnuplot mag3d.glt"
    ])
    dat='wan_ChiPMz.mat.syml1'
    tall+=test2_check(testdir+'/'+dat, workdir+'/'+dat)
    dat='wan_ChiPMz.mat.syml2'
    tall+=test2_check(testdir+'/'+dat, workdir+'/'+dat)
    dat='wan_ChiPMr.mat.syml1'
    tall+=test2_check(testdir+'/'+dat, workdir+'/'+dat)
    dat='wan_ChiPMr.mat.syml2'
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
