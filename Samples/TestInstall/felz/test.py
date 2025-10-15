import os
from comp import test1_check,test2_check,runprogs,rmfiles
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
        out1='out.lmf.fsmom.felz'
        out2='out.lmf.lzsz.felz'
        out3='out.lmf.ls.felz'
        rmfiles(workdir,[out1,out2,out3])
        runprogs([
                 lmfa+ " -vrel=1 -vso=0 felz > "+out1,
                 lmf + " -vrel=1 -vnit=3 -vso=2 felz -vfsmom=-2 >> "+out1 ,
	         "rm -f atm.* fs.* moms.* *mixm.* rst.* save.* log.* *hssn.* wkp.* bsmv.* syml.* bnds.*"
        ])
        tall+=test1_check(testdir+'/'+out1, workdir+'/'+out1)
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
        runprogs([
                 lmfa + " -vrel=1 -vnit=1 -vso=0 felz > "+ out2,
                 lmf  + " -vrel=1 -vnit=1 -vso=2 felz -vpwmode=11 >> "+out2,
                 "rm -f *mixm.felz",
                 lmf  + " -vrel=1 -vnit=1 -vso=2 felz >> "+out2 
        ])
        tall+=test1_check(testdir+'/'+out2, workdir+'/'+out2)
        runprogs([
                 lmf+" -vrel=1 -vnit=1 -vso=1 felz > "+out3 
        ])
        tall+=test1_check(testdir+'/'+out3, workdir+'/'+out3)
        return tall
