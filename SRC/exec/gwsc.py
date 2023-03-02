#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
QSGW calculation with MPI.
'''
import os, datetime, shutil, glob
from gwutil import *
'''
    main program of gwsc.py
'''
#get argments
target,nloop,ncore,ncore2=Obtain_args()
epath=os.path.dirname(os.path.abspath(__file__)) #if gwsc.py is in ecalj bin dir
tmp=[int(s.strip('RUN.ITER')) for s in os.listdir() if 'RUN.ITER' in s]
Iter0=max(tmp) if len(tmp)!=0 else 0
initxt=("### START gwsc: ITER= %d, MPI size=  %d, TARGET= %s"%(nloop,ncore,target)
        if ncore==ncore2 else
        "### START gwsc: ITER= %d, MPI size=  %d, %d, TARGET= %s"%(nloop,ncore,ncore2,target))
print(initxt,flush=True)
#initial directry config
gen_dir('SEBK')
gen_dir('STDOUT')
#set target fine name
sigm_name='sigm.%s'%target
rst_name='rst.%s'%target
#serch sigm and sigm.target files
if os.path.isfile('sigm') and os.path.isfile(sigm_name):
    shutil.move(sigm_name, sigm_name+'.bakup')
    os.symlink('sigm',sigm_name)
    print('--- sigm is used. sigm.$TARGET is softlink to it  ---',flush=True)
elif os.path.isfile(sigm_name):
    shutil.move(sigm_name, sigm)
    os.symlink('sigm',sigm_name)
    print('--- sigm.$TARGET is moved to sigm. sigm.$TARGET is softlink now.  ---',flush=True)
else:
    print('--- No sigm nor sigm.$TARGET files for starting ---',flush=True)
#main iteration
for i in range(nloop+1):
    niter=i+Iter0
    if not (Iter0!=0 and i==0):
        print(" ---- goto sc calculation for given sigma-vxc --- ix=%d"%niter,flush=True)
        #lmf calculation
    if Iter0!=0 and i==0:
        '''
        if continue previous calc, we do not need initial lda calculation
        skip meaning is keeping consistency of iteration number
        '''
        continue
    elif i==0:
        if os.path.isfile(sigm_name):
            print(" we have sigm already, skip iter=0",flush=True)
        else:
            print("No sigm ---> LDA caculation for eigenfunctions ",flush=True)
            remove('llmf')
            run_program(100,epath,'lmf-MPIK','llmf_lda',ncore,target)
            shutil.copy(rst_name,rst_name+'.lda')
    else:
        run_program(100,epath,'lmf-MPIK','llmf',ncore,target)
    #gw initialize
    run_program(0,epath,'lmf-MPIK','llmfgw00',0,target+' --jobgw=0')
    run_program(1,epath,'qg4gw','lqg4gw --job=1')
    run_program(1,epath,'lmf-MPIK','llmfgw01',ncore,target+' --jobgw=1')
    ##### main stage of gw ####
    run_program(1,epath,'rdata4gw_v2','lrdata4gw_v2') #prepare files
    mv_files(glob.glob('norm*'),'STDOUT')
    rm_files(['gwa']+glob.glob('gwb*'))
    run_program(1,epath,'heftet','leftet --job=1')       # A file EFERMI for hx0fp0
    ### Core part of the self-energy (exchange only) ###
    run_program(3,epath,'hbasfp0','lbasC --job=3')       # Product basis generation
    run_program(3,epath,'hvccfp0','lvccC --job=3',ncore) # Coulomb matrix for lbasC
    run_program(3,epath,'hsfp0_sc','lsxC --job=3',ncore2) # Sigma from core1
    mv_files(glob.glob('stdout*'),'STDOUT')
    ### Valence part of the self-energy Sigma ###
    run_program(0,epath,'hbasfp0','lbas --job=0')        # Product basis generation
    run_program(0,epath,'hvccfp0','lvcc --job=0',ncore) # Coulomb matrix for lbas
    run_program(1,epath,'hsfp0_sc','lsx --job=1',ncore2) # Exchange Sigma
    mv_files(glob.glob('stdout*'),'STDOUT')
    run_program(11,epath,'hx0init','lx0init --job=11',ncore2)    #x0 part
    run_program(100,epath,'hx0zmel','lx0zmel',ncore2)   #x0 part
    run_program(100,epath,'hrcxq','lrcxq',ncore2)       #x0 part
    run_program(100,epath,'hhilbert','lhilbert',ncore2) #x0 part
    mv_files(glob.glob('stdout*'),'STDOUT')
    run_program(2,epath,'hsfp0_sc','lsc --job=2',ncore2) #correlation Sigma
    mv_files(glob.glob('stdout*'),'STDOUT')
    run_program(0,epath,'hqpe_sc','lqpe')        #all Sigma are combined.
    ### final part of iteration loop. Manupulate files ###
    remove(sigm_name)
    os.symlink('sigm',sigm_name)
    mv_files(glob.glob('SEX*')+glob.glob('SEC*')+glob.glob('XC*'),'SEBK')
    for f in ['sigm','QPU','QPD','TOTE.UP','TOTE.DN','lqpe','lsc','lsx','lx0','llmfgw01','evecfix.chk','llmf','ESEAVR']:
        if os.path.isfile(f):
            shutil.copy(f,f+'.%drun'%niter)
    if i==0 and Iter0!=0 and nloop!=0:
        run0='RUN0'
        gen_dir(run0)
        run_program(100,epath,'lmf-MPIK','llmf_oneshot',ncore,target)
        cp_files(['ctrl.%s'%target,rst_name,sigm_name,'llmf_oneshot','save.%s'%target],run0)
    rundir='RUN.ITER%d'%niter
    gen_dir(rundir)
    cp_files(['ctrl.%s'%target,rst_name,sigm_name,'GWinput','save.%s'%target],rundir)
    print('OK! --> == %d iteration over =='%niter,flush=True)
    ##### end of loop #####

### finally we have llmf_gwscend ###
run_program(100,epath,'lmf-MPIK','llmf_gwscemd.%d'%niter,ncore,target)
rm_files(glob.glob('ewindow.%s*')+glob.glob('qbyl.%s*'%target)+glob.glob('eigze*.%s*'%target)+['_IN_'])
if nloop==0:
    run0='RUN0'
    gen_dir(run0)
    cp_files(['ctrl.%s'%target,rst_name,sigm_name,'llmf_gwscend.%d'%niter,'save.%s'%target],run0)
    argv=''
    print('OK! ==== All calclation finished for  gwsc %s ===='%argv)

    
