from comp import runprogs,diffnum
def test(args,bindir,testdir,workdir):
        genm= bindir + f'/genMLWFx '
        np4= f'-np {args.np} '
        tall=''
        dfile= ["Screening_W-v.h", "Screening_W-v_crpa.h"]
        runprogs([
                 *[f"rm -rf {workdir}/{fname}" for fname in dfile],
                 genm + " srvo3 "+ np4,
                 "head -1000 Screening_W-v.UP > Screening_W-v.h",
                 "head -1000 Screening_W-v_crpa.UP > Screening_W-v_crpa.h"
        ])
        for outfile in dfile:   
                tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=1e-3,comparekeys=[])
        return tall
