import os
from comp import test1_check,test2_check,runprogs
def test(args,bindir,testdir,workdir):
        lmfa= f'mpirun -np 1 {bindir}/lmfa '
        lmf = f'mpirun -np {args.np} {bindir}/lmf '
        tall=''
        message1='''
        # Case felz: spin-polarized Fe spin-orbit coupling
        # --- Test 3.  Check of miscellaneous special features, programs lmfa,lmf ---
         The Fe test tests the code's implementation of fixed-spin moment
         with and without spin-orbit coupling
        '''
        print(message1)
        outfile='out.lmf.fsmom.felz'
        runprogs([
                 lmfa+ " -vrel=1 -vso=0 felz > "+outfile,
                 lmf + " -vrel=1 -vnit=3 -vso=2 felz -vfsmom=-2 >> "+outfile ,
	         "rm -f atm.* fs.* moms.* *mixm.* rst.* save.* log.* *hssn.* wkp.* bsmv.* syml.* bnds.*"
        ])
        tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
        message1='''
        # --- Test case 4:  Spin-orbit coupling ---
         The felz test computes the orbital moment in Fe.
         lmf calculates the orbital magnetic moment.
           * In the first part of this test only LzSz is used.
           * The APW basis with LzSz is also checked.
	   * In the second part the FULL SPIN ORBIT is used.
           * Symmetry operations must be suppressed at present.
           * Only 4x4x4 k points are used in this test.
        '''
        outfile='out.lmf.lzsz.felz'
        runprogs([
                 lmfa + " -vrel=1 -vnit=1 -vso=0 felz > "+ outfile,
                 lmf  + " -vrel=1 -vnit=1 -vso=2 felz -vpwmode=11 >> "+outfile,
                 "rm -f *mixm.felz",
                 lmf  + " -vrel=1 -vnit=1 -vso=2 felz >> "+outfile 
        ])
        tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
        outfile='out.lmf.ls.felz'
        runprogs([
                 lmf+" -vrel=1 -vnit=1 -vso=1 felz > "+outfile 
        ])
        tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
        return tall