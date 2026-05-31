import re
from comp import test2_check,runprogs
def test(args,bindir,testdir,workdir):
    MATERIAL="cu"
    ncore=args.np
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    outfile=f'out.lmf.{MATERIAL}'
    tall=''
    if(args.checkonly): runprogs(["rm -rf summary.txt"],quiet=True)
    else: runprogs([
        lmfa + f" {MATERIAL} > "+ outfile,
        lmf  + f" {MATERIAL} > "+ outfile,
        f"{bindir}/getsyml   {MATERIAL}",
        f"{bindir}/job_band  {MATERIAL} -np {ncore}",
        f"{bindir}/job_eps   {MATERIAL} -np {ncore} --decompose",
    ])
    skipcond = lambda line: any( abs(float(tok)) >= 1e4
                                 for tok in re.split(r'\s+', line.strip())[4:]
                                 if tok != '' )
    for file in ['EPS0001','EPS0002','EPS0003']:
        dat= file+'.nlfc.dat.interbandonly'
        print(dat,end=': ')
        tall+=test2_check(testdir+'/'+dat, workdir+'/'+dat)
        dat= file+'.nlfc.dat.intrabandonly'
        print(dat,end=': ')
        tall+=test2_check(testdir+'/'+dat, workdir+'/'+dat, skipcond=skipcond, rel_tol=0.01)
    print(f'''
     ======================================================
     See plots
        "gnuplot -p eps_interbandonly_{MATERIAL}.glt"  interband
        "gnuplot -p eps_intrabandonly_{MATERIAL}.glt"  intraband
        "gnuplot -p eps_total_{MATERIAL}.glt"          total
     ======================================================
    ''')
    return tall
