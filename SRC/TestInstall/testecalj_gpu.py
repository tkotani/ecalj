#!/usr/bin/env python3
''' test system for ecalj '''
import os, datetime, shutil, glob,sys
from comp import comp,compall,compeval,diffnum,dqpu,test1_check,test2_check

"Set directory"
testroot=os.getcwd()
ecaljroot= testroot+'/../../'
bindir=testroot+'/bin'
workroot=testroot+'/work'
shutil.rmtree(workroot,ignore_errors=True)
os.mkdir(workroot)

"Set arguments"
ttall=''
npsize='4' #default
np=False
usegpu=False
showt=False
usemp=False
option=''
for arg in sys.argv[1:]:
    if(arg=='--help'):
        ttall=''
        showt=True
        break
    if(np):
        npsize= arg
        np=False
        continue
    if(arg=='-np'):
        np=True
        continue
    if(arg=='--gpu'):
        usegpu=True
        continue
    if(arg=='--mp'):
        usemp=True
        continue
    if(arg=='gwall'):
        ttall="gas_eps_lmfh gas_epsPP_lmfh fe_epsPP_lmfh_chipm si_gw_lmfh gas_pw_gw_lmfh si_gwsc gas_gwsc nio_gwsc fe_gwsc ni_crpa srvo3_crpa"
        continue
    #if(arg[0:2]=='--'):
    #    option = option+" "+arg
    ttall=ttall+" "+arg
#print('option=',option)
#sys.exit()
np4 = "-np "+npsize+" " # mpisize
if(len(ttall)==0): #when no tests are specified
    testall="copt te zrt co cr3si6 felz gasls eras c crn cu na "\
        +"gas_eps_lmfh gas_epsPP_lmfh fe_epsPP_lmfh_chipm si_gw_lmfh gas_pw_gw_lmfh si_gwsc gas_gwsc nio_gwsc fe_gwsc ni_crpa srvo3_crpa"
else:
    testall=ttall
print('mpi size= ',np4)
print('run tests for ',testall)
if(showt):
    print('usage: testecalj.py -np 4 [test names]')
    sys.exit()

"Install to bin"
os.makedirs(bindir,exist_ok=True)
exec='lmfa lmf run_arg job_pdos job_tdos ctrl2ctrlp.py a2vec.py \
 gwsc qg4gw hvccfp0 hsfp0_sc hqpe_sc hmaxloc hpsig_MPI huumat_MPI hwmatK_MPI hrcxq \
 hsfp0_sc_gpu hrcxq_gpu hx0fp0_gpu \
 hsfp0_sc_mp hrcxq_mp hx0fp0_mp \
 hsfp0_sc_mp_gpu hrcxq_mp_gpu hx0fp0_mp_gpu \
 heftet hbasfp0 gw_lmfh hx0fp0 hsfp0 hqpe eps_lmfh epsPP_lmfh epsPP_lmfh_chipm genMLWFx'
for ex in exec.split():
    if os.path.exists(ecaljroot+'/SRC/exec/'+ex):
        shutil.copy(ecaljroot+'/SRC/exec/'+ex,bindir)
        print ('cp ' + ecaljroot+'/SRC/exec/' +ex+ ' to ',bindir)
#sys.exit()

"Alias binaries combined with np "
lmfa = 'mpirun -np 1 '+testroot+'/bin/lmfa ' 
lmf  = 'mpirun '+np4+ testroot+'/bin/lmf '
eps_lmfh = testroot+'/bin/eps_lmfh '+np4
epsPP_lmfh = testroot+'/bin/epsPP_lmfh '+np4
epsPP_lmfh_chipm = testroot+'/bin/epsPP_lmfh_chipm '+np4
gw_lmfh = testroot+'/bin/gw_lmfh '+np4
gwsc0  = testroot+'/bin/gwsc 0 '+np4
genm   = testroot+'/bin/genMLWFx '
job_pdos= testroot+'/bin/job_pdos '
if usegpu:
    gwsc0 = testroot+'/bin/gwsc --gpu 0 '+np4
if usemp:
    gwsc0 = testroot+'/bin/gwsc --mp 0 '+np4
if usemp and usegpu:
    gwsc0 = testroot+'/bin/gwsc --gpu --mp 0 '+np4

def runprogs(runlist):
    for irun in runlist:
        print(irun)
        err = os.system(irun)
        if(err): print('Error exit!')
        if(err): sys.exit(-1)
                
"Start test as tname from testall"
tall=''
#try:
#    from tqdm import tqdm
#except:
#    def tqdm(aaa):
#        aaa
#for tname in tqdm(testall.split()):
for tname in testall.split():
    testdir=testroot+'/'+tname
    workdir=workroot+'/'+tname
    print()
    print('=== Test ',tname,' at ',workdir)
    try:
        shutil.copytree(testdir, workdir)
    except:
        print('Error! No such test: testdir=',testdir)
        sys.exit(-1)
    os.chdir(workdir)
    outfile='out.lmf.'+tname
    if(tname=='copt'):
        message1='''
        # Case copt: a distorted L12 environment with four atoms. 
        # --- Basic check of programs lmfa,lmf ---
        '''
        print(message1)
        runprogs([
                 lmfa+'copt -vnsp=2 -cstrx=l12 --pr41 -vmet=3 -vtetra=0 -vnk1=2 -vlfrce=12 -vdist=1 -vnit=3 --time=5 > '+outfile,
                 lmf +'copt -vnsp=2 -cstrx=l12 --pr41 -vmet=3 -vtetra=0 -vnk1=2 -vlfrce=12 -vdist=1 -vnit=3 --time=5 > '+outfile 
        ])
        tall+= test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    elif(tname=='te'):
        message1='''
        # Case te: molecular statics in an open structure
        # --- Test 1.  Basic check of programs lmfa,lmf ---
         The te test also checks and illustrates the following:

          1.  a simple molecular relaxation with one symmetry-allowed degree of freedom
          2.  an spd-sp basis augmented by floating orbitals, and also by q-dependent APWs
          3.  Comparison of forces and total energy with and without floating orbitals, and with and without APWs.

         lmf will first relax the atoms with the basis including floating orbitals.
	 
         After relaxation, the a new calculation is performed that remove floating orbitals
         but adding plane waves to the basis, so the total energy and forces may be compared. 
         The basis size is variable, but averages ~80 orbitals, a little more than the floating
         orbitals case (~70 orbitals). About 3 mRy is gained relative to the floating orbitals case.

         Note that KMXA=5 is used with the PW calculation.  It isn't necessary in this case; still the
         user is cautioned to monitor this parameter when using energy cutoffs higher than 3 Ry or so.

         As a last step the calculation is repeated at the relaxed position with only atom-centered
         MTO's (neither floating orbitals nor plane waves).  Elimination of floating orbitals reduces 
         the basis from 75 to 39 orbitals, and reduces the LDA total energy by about 4 mRy.
	 
         The forces are not strongly affected by the local orbitals (or APWs), as can be seen by
         looking at the maximum force after the last step.
	'''
        print(message1)
        runprogs([\
            lmfa+"te -vdyn='DYN' -vnk=3 -vnit=3 -vlf1=4 -vlmxl=4 -vnk=3 -vngd=20 -vkmx=3 -vconv=1e-4 > "+outfile, 
            lmf+ "te -vdyn='DYN' -vnk=3 -vnit=3 -vlf1=4 -vlmxl=4 -vnk=3 -vngd=20 -vkmx=3 -vconv=1e-4 -vnbas=12 -vnspec=2 >> "+outfile,
            "rm -f mixm.te",
            "cp rst.te rst.te.bk",
            lmf+"te -vnk=3 -vnit=3 -vlf1=4 -vlmxl=4 -vnk=3 -vngd=20 -vkmx=5 -vconv=1e-4 -vpwmode=11 >> "+outfile,
            "rm -f mixm.te",
            "cp rst.te.bk rst.te",
            lmf+"te -vnk=3 -vnit=3 -vlf1=4 -vlmxl=4 -vnk=3 -vngd=20 -vkmx=3 -vconv=1e-4 >> "+outfile 
        ])
        tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    elif(tname=='zrt'):
        message1='''
        # Case zrt: ZrO_2 fluorite in tetragonal setting, with tetragonal distortion
        # --- Test 1.  Basic check of programs lmfa,lmf ---
        The zrt test also checks and illustrates the following:
         1.  the use of restart files in both binary and ascii forms
         2.  two-kappa basis
        '''
        print(message1)
        runprogs([
                 lmfa+" -vfp=1 zrt --no-iactiv > out.lmf.zrt > "+outfile,
                 lmf+"  -vnitq=1 -vforce=1 -vfp=1 zrt >> "+outfile,
                 lmf+"  -vnitq=1 -vforce=1 -vfp=1 zrt >> "+outfile 
        ])
        tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    elif(tname=='co'):
        message1='''
        # Case co: a hexagonal environment with two equivalent atoms.
        # --- Test 1.  Basic check of programs lmfa,lmf ---
         The co test also checks and illustrates the following:
         1.  a spin-polarized case
         2.  tetrahedron method (METAL=3)
         3.  Constrained mixing (first spin is kept frozen, then charge, then both are allowed to change)
         4.  Broyden mixing
         5.  bands mode (see command-line argument --band in last lmf invocation)
             SO coupling and color weights for projection onto Co d channels are included.
             Two separate weights are made (one each for d majority and minority bands.)
             This script will prompt you to see whether it should create 
             a figure from the bands file (The fplot plotting package needs
             to be installed for this function).
         '''       
        print(message1)
        runprogs([
                 lmfa+" co -vmet=3 -vlmf=1 -vnk=8 -vnit=3 --pr31 > "+ outfile,
                 lmf+ " co -vmet=3 -vlmf=1 -vnk=8 -vnit=3 -vw2=0 --pr31 >> "+outfile,
                 "rm mixm.co",
                 lmf+ " co -vmet=3 -vlmf=1 -vnk=8 -vnit=3 -vw1=0 --pr31 >> "+outfile,
                 "rm mixm.co",
                 lmf+ " co -vmet=3 -vlmf=1 -vnk=8 -vnit=3 --pr31 --time=5 >> "+outfile,
                 lmf+ " co -vmet=3 -vnk=8 -vnit=3 --pr31  -vso=t --band:fn=syml >> "+outfile,
                 "rm -f atm.* mixm.* rst.* save.* log.* hssn.* wkp.* dos.* tdos.* pdos.* dos-mull.* qpp.* out.lmf-dos*"
        ])
        tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
        outbnds='bnds.co'
        tall+=test2_check(testdir+'/'+outbnds, workdir+'/'+outbnds)
        message1='''
        # --- Test 2.  Core-level spectroscopy (EELS), Mulliken analysis, partial DOS ---
         The Co test case illustrates partial dos resolved by both l and m.
         Note: because the routine generating partial DOS within MT spheres
         does not properly symmetrize it, symmetry operations must be
         suppressed when resolving DOS by m.
        '''
        print(message1)
        runprogs([
                 lmfa+" co -vmet=3 -vlmf=1 -vnk=8 -vnit=1 --pr31 > out.lmf-dos.co",\
                 job_pdos+" co "+ np4 +" -vmet=3 -vlmf=1 -vnk=8 -vnit=1 --pr31 ---NoGnuplot > out.lmf-dos.co"\
        ])
        pdos2002='dos.isp2.site002.co'
        tall+=test2_check(testdir+'/'+pdos2002, workdir+'/'+pdos2002)
    elif(tname=='cr3si6'):
        message1='''
        # Case cr3si6: a hexagonal environment with several atoms, two classes
        # --- Test 1.  Basic check of programs lmfa,lmf ---
         Case cr3si6: a hexagonal environment with several atoms, two classes
           Other checks:       verbose output, insulator, forces(mode 1)
         --- Test 1.  Basic check of programs lmfa,lmf ---
         Checks that program lmfa produces a sensible atm file
         and that program lmf iterates to the proper energy.
         The cr3si6 test also checks and illustrates the following:
         1.  insulator
         2.  forces with correction mode 1 (FORCES=1)
         3.  verbose output
         4.  A case with low kmxa (kmxa=2)
         '''       
        print(message1)
        runprogs([
                 lmfa+" cr3si6 --pr51 -vnit=2 --time=6 >  "+outfile,
                 lmf +" cr3si6 --pr51 -vnit=2 --time=6 >> "+outfile
        ])
        tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    elif(tname=='felz'):
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
	         "rm -f atm.* fs.* moms.* mixm.* rst.* save.* log.* hssn.* wkp.* bsmv.* syml.* bnds.*"
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
                 "rm -f mixm.felz",
                 lmf  + " -vrel=1 -vnit=1 -vso=2 felz >> "+outfile 
        ])
        tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
        outfile='out.lmf.ls.felz'
        runprogs([
                 lmf+" -vrel=1 -vnit=1 -vso=1 felz > "+outfile 
        ])
        tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    elif(tname=='gasls'):
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
        compeval(testdir+'/'+outfile, workdir+'/'+outfile,' 0.00000  0.00000  0.00000',lineeval=3,evalso=5,tol=1e-4) 
    elif(tname=='eras'):
        message1='''
        # Case ErAs: Test of LDA+U
        # --- Test 1.  Basic check of programs lmfa,lmf ---
         The ErAs test also illustrates the LDA+U implementation for:
         1. a case when the LDA puts f orbitals at the Fermi energy.
         2. demonstration of U on both d and f orbitals
         (3. convergence to a metastable solution with a reasonable spin moment but wrong orbital moment.)
        '''
        print(message1)
        runprogs([
                 lmfa+" eras  > "+outfile,
                 lmf+" -vnit=1 --pr51 eras >> "+outfile,
                 lmf+" -vnit=3 eras        >> "+outfile 
        ])
        tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    elif(tname=='c'):
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
                 "rm -f mixm.* rst.* save.* log.* hssn.* wkp.* bsmv.* bnds.*"
        ])
        tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
        message1='''
        # Case C: test of homogeneous background
        continue
        '''
        outfile='out.lmf.ionized.c'
        runprogs([
                 lmfa+" c -vzbak=1 > "+outfile,
                 lmf+ " c -vzbak=1 >>"+outfile,
                 "rm -f mixm.* rst.* save.* log.* hssn.* wkp.* bsmv.* bnds.*"
        ])
        tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    elif(tname=='crn'):
        message1='''
        # Case CrN: test of CLS with core hole
        # --- Test 2.  Core-level spectroscopy (EELS), Mulliken analysis, partial DOS ---
         The CrN test case generates core-level spectroscopy for the
         1s state in N.  The self-consistent calculation proceeds with
         an electron missing from the N 1s core, which corresponds to
         the 'sudden approximation' (system relaxes instantanously
         from electron exited out of hole).
        '''
        print(message1)         #outfile='out.lmf-dos.crb'
        runprogs([
                 lmfa+" crn > "+outfile ,
                 lmf+ " crn >>"+outfile,
                 lmf+ " -vnit=1 -vmetal=2 crn >>"+outfile,
                 lmf+ " --cls crn >>"+outfile
        ])
        outfile='dos-vcdmel.crn'
        tall+=test2_check(testdir+'/'+outfile, workdir+'/'+outfile)
    elif(tname=='cu'):
        message1='''
        # Case cu: illustration of high-lying local orbitals and bands of Cu up to ~50 eV.
        # --- Test 1.  Basic check of programs lmfa,lmf ---
         The cu test also illustrates the following:
         1.  high-lying local orbitals (Cu 5s,5p,4d are included as local orbitals)
         2.  METAL=3 for BZ integration
         3.  bands mode (see command-line argument --band in last lmf invocation)
        '''
        print(message1)
        runprogs([
                 lmfa+" cu  > "+outfile ,
                 lmf+ " cu -vnk=8 -vbigbas=f >>"+outfile,
                 "rm mixm.cu",
                 lmf+ " cu -vnk=8 -vbigbas=t -vmetal=3 -vrsm2=1.3 -vrsmd1x=1 -vlmx=4 -vpwmode=0 -voveps=0d-7 >>"+outfile,
                 lmf+ " cu -vnk=8 -vbigbas=t -vmetal=3 -vrsm2=1.3 -vrsmd1x=1 -vlmx=4 -vpwmode=0 -voveps=0d-7 --band:fn=syml >>"+outfile
        ])
        tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
        outfile='bnds.cu'
        tall+=test2_check(testdir+'/'+outfile, workdir+'/'+outfile)
    elif(tname=='na'):
        message1='''
        # Case na: illustration of low- and high-lying local orbitals
        # --- Test 1.  Basic check of programs lmfa,lmf ---
         The na test also illustrates the following:
         1.  compare the total energy for conventional and extended Na 2p orbitals
             After the test finishes, compare the three energies with out.lmf.na
        '''
        print(message1)
        runprogs([
                 lmfa+" na -vnk=6 -vnapval=0 -vpz1=0 -vp1=3.38  > "+outfile ,
                 lmf+ " na -vnk=6 -vnapval=0 -vpz1=0 -vp1=3.38 >>"+outfile,
                 "rm mixm.na rst.na",
                 lmfa+" na -vnk=6  -vnapval=1 >>"+outfile,
                 lmf+ " na -vnk=6  -vnapval=1>>"+outfile,
                 "rm mixm.na rst.na",
                 lmf+ " na -vnk=6 -vnapval=2 -vpz1=12.94 >>"+outfile
        ])
        tall+=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    elif(tname=='gas_eps_lmfh'):
        runprogs([
                 lmfa+" gas  > llmfa" ,
                 lmf+ " gas  > llmf",
                 "rm EPS*",
                 eps_lmfh+" gas" 
        ])
        epsfile="EPS0001.nlfc.dat EPS0002.nlfc.dat EPS0003.nlfc.dat EPS0004.nlfc.dat EPS0001.dat EPS0002.dat EPS0003.dat EPS0004.dat"
        for outfile in epsfile.split(): tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=1e-3,comparekeys=[])
    elif(tname=='gas_epsPP_lmfh'):
        runprogs([
                 lmfa+" gas  > llmfa" ,
                 lmf+ " gas  > llmf",
                 "rm EPS*",
                 epsPP_lmfh+" gas" 
        ])
        epsfile="EPS0001.nlfc.dat EPS0002.nlfc.dat EPS0003.nlfc.dat EPS0004.nlfc.dat"
        for outfile in epsfile.split(): tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=3e-3,comparekeys=[])
    elif(tname=='fe_epsPP_lmfh_chipm'):
        runprogs([
                 lmfa+" fe  > llmfa" ,
                 lmf+ " fe  > llmf",
                 "rm ChiPM*",
                 epsPP_lmfh_chipm+" fe" 
        ])
        epsfile="ChiPM0001.nlfc.mat ChiPM0002.nlfc.mat ChiPM0003.nlfc.mat ChiPM0004.nlfc.mat ChiPM0005.nlfc.mat"
        for outfile in epsfile.split(): tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=3e-3,comparekeys=[])
    elif(tname=='si_gw_lmfh'):
        runprogs([
                 lmfa+" si  > llmfa" ,
                 lmf+ " si  > llmf",
                 "rm QPU",
                 gw_lmfh+" si" 
        ])
        dfile="QPU"
        for outfile in dfile.split():   tall+=dqpu(testdir+'/'+outfile, workdir+'/'+outfile)
    elif(tname=='gas_pw_gw_lmfh'):
        runprogs([
                 lmfa+" gas  > llmfa" ,
                 lmf+ " gas  > llmf",
                 "rm QPU",
                 gw_lmfh+" gas" 
        ])
        epsfile="QPU"
        for outfile in epsfile.split(): tall+=dqpu(testdir+'/'+outfile, workdir+'/'+outfile)
    elif(tname=='si_gwsc'):
        runprogs([
                 lmfa+" si  > llmfa",
                 "rm log.si QPU",
                 gwsc0+ " si"
        ])
        dfile="QPU"
        for outfile in dfile.split():   tall+=dqpu(testdir+'/'+outfile, workdir+'/'+outfile)
        outfile='log.si'
        tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=3e-3,comparekeys=['fp evl'])
    elif(tname=='gas_gwsc'):
        runprogs([
                 lmfa+" gas  > llmfa" ,
                 "rm -f log.gas QPU",
                 gwsc0+ " gas",
        ])
        epsfile="QPU"
        for outfile in epsfile.split(): tall+=dqpu(testdir+'/'+outfile, workdir+'/'+outfile)
        outfile='log.gas'
        tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=3e-3,comparekeys=['fp evl'])
    elif(tname=='nio_gwsc'):
        runprogs([
                 lmfa+" nio  > llmfa" ,
                 "rm -f log.nio QPU QPD",
                 gwsc0+ " nio",
        ])
        dfile="QPU"
        for outfile in dfile.split():   tall+=dqpu(testdir+'/'+outfile, workdir+'/'+outfile)
        outfile='log.nio'
        tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=3e-3,comparekeys=['fp evl'])
    elif(tname=='fe_gwsc'):
        runprogs([
                 lmfa+" fe  > llmfa" ,
                 "rm -f log.fe QPU QPD",
                 gwsc0+ " fe",
        ])
        dfile="QPU QPD"
        for outfile in dfile.split():   tall+=dqpu(testdir+'/'+outfile, workdir+'/'+outfile)
        outfile='log.fe'
        tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=3e-3,comparekeys=['fp evl'])
    elif(tname=='ni_crpa'):
        runprogs([
                 genm + " ni "+ np4,
                 "head -1000 Screening_W-v.UP > Screening_W-v.h",
                 "head -1000 Screening_W-v_crpa.UP > Screening_W-v_crpa.h"
        ])
        dfile="Screening_W-v.h Screening_W-v_crpa.h"
        for outfile in dfile.split():   tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=1e-3,comparekeys=[])
    elif(tname=='srvo3_crpa'):
        runprogs([
                 genm + " srvo3 "+ np4,
                 "head -1000 Screening_W-v.UP > Screening_W-v.h",
                 "head -1000 Screening_W-v_crpa.UP > Screening_W-v_crpa.h"
        ])
        dfile="Screening_W-v.h Screening_W-v_crpa.h"
        for outfile in dfile.split():   tall+=diffnum(testdir+'/'+outfile, workdir+'/'+outfile,tol=1e-3,comparekeys=[])
    else:
        print('such test id have not given yet')
        sys.exit()
        
    print('=== EndOf',tname,' at ',workdir)
    print()
    os.system('cat ../summary.txt')
    if('err' in tall):
        print('FAILED at some tests ===')
    else:    
        print('OK! ALL PASSED ===')
print("   See work/summary.txt")
