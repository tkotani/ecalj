from comp import test1_check,runprogs
def test(args,bindir,testdir,workdir):
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    tall=''
    message1='''
    # CAUTION! 2024-3-12: this documents is old. Need to be examined. Only for memo 
    # Case C: test of homogeneous background
    # --- Test 3.  Check of miscellaneous special features, programs lmfa,lmf ---
    The C tests the code's implementation of homogeneous background mode
    It checks that program lmf generates the correct total energy and 
    ionization potential for a single C atom.
    *The total energy of the neutral atom computed by lmf (-74.996 Ry) is very close
    to the free-atom energy computed by lmfa (-74.995 Ry), and that the free-atom 
    potential is approximately self-consistent; 
    *Also the SUM of total energy of the ionized system (-74.443 (newcode2024 gives -74.545, right?)) 
    and the interaction of a point charge with a homogeneous background E^I,
    estimated from E^I=9/5R, with 4*pi*R^3/3=vol => R=6.20 and E^I = 0.290 Ry
    evaluates to -74.443+0.290 = -74.153 Ry, is close to the total energy of
    the positively charged free ion as computed by lmfa (-74.171 Ry).
    We should have; ---------------------
    Energy of charged system         = .5524587
    Estat energy q*q/9/5/<r>         = .2902
      Corrected charged system energy  =  0.842659  <---
      Energy of neutral system         =  -.0012168  <---
    ----------------------------------------------         
    difference                        = 0.843876
    '''
    print(message1)
    outfile='out.lmf.neutral.c'
    runprogs([
        lmfa+" c -vzbak=0 > "+outfile,
        lmf+ " c -vzbak=0 >>"+outfile,
        "rm -f *mixm.* rst.* save.* log.* *hssn.* wkp.* bsmv.* bnds.*"
    ])
    tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    message1='''# Case C: test of homogeneous background \n continue'''
    outfile='out.lmf.ionized.c'
    runprogs([
        lmfa+" c -vzbak=1 > "+outfile,
        lmf+ " c -vzbak=1 >>"+outfile,
        "rm -f *mixm.* rst.* save.* log.* *hssn.* wkp.* bsmv.* bnds.*"
    ])
    tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    return tall
