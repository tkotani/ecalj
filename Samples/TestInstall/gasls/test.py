import os
from comp import test1_check,test2_check,runprogs,compeval
def test(args,bindir,testdir,workdir):
        lmfa= f'mpirun -np 1 {bindir}/lmfa '
        lmf = f'mpirun -np {args.np} {bindir}/lmf '
        message1='''
        # Case GaAs: GaAs with spin-orbit coupling
        # --- Test case 4:  Spin-orbit coupling ---
         The GaAs test computes the energy bands at (0,0,0) (Gamma point),
         (1/4,1/4,1/4) and (1/2,1/2,1/2) (L point).
	 The spin-orbit splitting of the valence states is tested.
         This test checks SO coupling in conjunction with conventional local orbitals.
        '''
        outfile='out.lmf.ls-bands.gasls'
        runprogs([
                 lmfa+" -vso=1 gasls -vpwmode=0 >"+outfile,
                 lmf+ " -vso=1 -vnit=1 gasls --band:fn=syml -vpwmode=0 >>"+outfile
        ])
        print(message1)
        tall=compeval(testdir+'/'+outfile, workdir+'/'+outfile,' 0.00000  0.00000  0.00000',lineeval=3,evalso=5,tol=1e-4) 
        return tall
