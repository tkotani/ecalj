from comp import runprogs,diffnum
def test(args,bindir,testdir,workdir):
        genm= bindir + f'/genMLWFx '
        np4= f'-np {args.np} '
        tall=''
        runprogs([
                 genm + " ni "+ np4,
                 "head -1000 Screening_W-v.UP > Screening_W-v.h",
                 "head -1000 Screening_W-v_crpa.UP > Screening_W-v_crpa.h"
        ])
        dfile="Screening_W-v.h Screening_W-v_crpa.h"
        for outfile in dfile.split():   
                tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=1e-3,comparekeys=[])
        return tall
