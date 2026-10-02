from comp import runprogs,diffnum,rmfiles
# cRPA of the MLO model (2026-10-02; until then the Wannier path, genMLWFx, git tag last-wannier):
# DFT -> bands (job_band) -> MLO model (job_mlo) -> v, W-v and the cRPA W-v in the MLO basis (job_mloW --crpa).
# Since 2026-10-02 11:44 job_mloW takes the Loewdin orthonormalized MLOs (projected Wannier functions of the MLO subspace;
# --mlo_raw for the raw MLOs); the reference files are in that basis.
# The MLO subspace is [mlo] mlo_lm (srvo3: V t2g, 3 orbitals); U, J at R=0, omega=0 in the first rows.
def test(args,bindir,testdir,workdir):
        np= f'-np {args.np} '
        tall=''
        out1=["Coulomb_v.h","Screening_W-v.h","Screening_W-v_crpa.h"]
        rmfiles(workdir,out1+['HamRsMLO'])
        runprogs([
                 f'mpirun -np 1 {bindir}/lmfa srvo3 > llmfa',
                 f'mpirun {np} {bindir}/lmf srvo3 > llmf',
                 f'{bindir}/job_band srvo3 {np} --NoGnuplot > ljob_band',
                 f'{bindir}/job_mlo srvo3 {np} --NoGnuplot > ljob_mlo',
                 f'{bindir}/job_mloW srvo3 {np} --crpa > ljob_mloW',
                 "head -300 Coulomb_v.UP > Coulomb_v.h",
                 "head -1000 Screening_W-v.UP > Screening_W-v.h",
                 "head -1000 Screening_W-v_crpa.UP > Screening_W-v_crpa.h"
        ])
        for outfile in out1:
                tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=1e-3,comparekeys=[])
        return tall
