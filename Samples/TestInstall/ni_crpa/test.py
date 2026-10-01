from comp import runprogs,diffnum,rmfiles
# cRPA of the MLO model (2026-10-02; until then the Wannier path, genMLWFx, git tag last-wannier):
# DFT -> bands (job_band) -> MLO model (job_mlo) -> v, W-v and the cRPA W-v in the MLO basis (job_mloW --crpa).
# The MLO subspace is [mlo] mlo_lm (ni: Ni d, 5 orbitals); U, J at R=0, omega=0 in the first rows.
def test(args,bindir,testdir,workdir):
        np= f'-np {args.np} '
        tall=''
        out1=["Coulomb_v.h","Screening_W-v.h","Screening_W-v_crpa.h"]
        rmfiles(workdir,out1+['HamRsMLO'])
        runprogs([
                 f'mpirun -np 1 {bindir}/lmfa ni > llmfa',
                 f'mpirun {np} {bindir}/lmf ni > llmf',
                 f'{bindir}/job_band ni {np} --NoGnuplot > ljob_band',
                 f'{bindir}/job_mlo ni {np} --NoGnuplot > ljob_mlo',
                 f'{bindir}/job_mloW ni {np} --crpa > ljob_mloW',
                 "head -300 Coulomb_v.UP > Coulomb_v.h",
                 "head -1000 Screening_W-v.UP > Screening_W-v.h",
                 "head -1000 Screening_W-v_crpa.UP > Screening_W-v_crpa.h"
        ])
        for outfile in out1:
                tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=1e-3,comparekeys=[])
        return tall
