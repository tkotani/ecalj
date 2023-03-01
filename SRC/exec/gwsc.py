#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
QSGW calculation with MPI.
'''
import os, datetime, shutil, glob
def Obtain_args():
    '''
    arguments settings
    Returns
      target: target material name
       nloop: iteration number of QSGW
       ncore: number of MPI thereads in lmf
      ncore2: number of MPI thereads in lxsC etc.
    '''
    import argparse
    parser=argparse.ArgumentParser(prog='gwsc.py',description='QSGW calculation script')
    parser.add_argument("nloop",help='iteration number of QSGW loop')
    parser.add_argument("mat_name",help='material name')
    parser.add_argument("-np",help='number of mpi core in lmf',action='store')
    parser.add_argument("-np2",help='number of mpi core in lxc etc.',action='store') 
    args=parser.parse_args()
    target=args.mat_name
    nloop=int(args.nloop)
    if args.np!=None:
        ncore=int(args.np)
        if args.np2!=None:
            ncore2=int(args.np2)
        else:
            ncore2=ncore
    else:
        ncore=1
        ncore2=ncore
    return(target,nloop,ncore,ncore2)

def gen_dir(dirname):
    '''
    serch directry and else generate that one
    Arguments
       dirname: the name of serch or generate directry
    '''
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

def cp_files(files,cp_dir):
    for fname in files:
        if os.path.isfile(fname):
            shutil.copy(fname,cp_dir+'/'+fname)

def mv_files(files,mv_dir):
    for fname in files:
        shutil.move(fname,mv_dir+'/'+fname)

def remove(fname):
    if os.path.isfile(fname):
        os.remove(fname)

def rm_files(files):
    for fname in files:
        remove(fname)

def run_codes(argin,epath,command,output,mpi_size=0,target=''):
    '''
    run ecalj code and write std output
    Arguments
          argin: argument of ecalj code 
          epath: path of ecalj bin
        command: run command
         output: output file
    Optional arguments
       mpi_size: number of MPI parallel
         target: target name
    '''
    import subprocess
    echo_arg='%d'%argin if argin !=100 else '---'
    mpirun='mpirun -np %d '%mpi_size if mpi_size!=0 else ''
    run_command = mpirun+epath+'/'+command+' '+target+' > '+output
    print('OK! --> Start echo '+run_command,flush=True)
    info=subprocess.run(run_command,shell=True)
    if info.returncode!=0: #if return error
        print('Error in '+command+' input_arg='+echo_arg+'. See OutputFile='+output,flush=True)
        exit()

def main():
    '''
    main program of gwsc.py
    '''
    #get argments
    target,nloop,ncore,ncore2=Obtain_args()
    epath=os.path.dirname(os.path.abspath(__file__)) #if gwsc.py is in ecalj bin dir
    #get previous iteration number
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
                run_codes(100,epath,'lmf-MPIK','llmf_lda',ncore,target)
                shutil.copy(rst_name,rst_name+'.lda')
        else:
            run_codes(100,epath,'lmf-MPIK','llmf',ncore,target)
        #gw initialize
        run_codes(0,epath,'lmf-MPIK','llmfgw00',0,target+' --jobgw=0')
        run_codes(1,epath,'qg4gw','lqg4gw --job=1')
        run_codes(1,epath,'lmf-MPIK','llmfgw01',ncore,target+' --jobgw=1')

        ##### main stage of gw ####
        run_codes(1,epath,'rdata4gw_v2','lrdata4gw_v2') #prepare files
        mv_files(glob.glob('norm*'),'STDOUT')
        rm_files(['gwa']+glob.glob('gwb*'))

        run_codes(1,epath,'heftet','leftet --job=1')       # A file EFERMI for hx0fp0
        ### Core part of the self-energy (exchange only) ###
        run_codes(3,epath,'hbasfp0','lbasC --job=3')       # Product basis generation
        run_codes(3,epath,'hvccfp0','lvccC --job=3',ncore) # Coulomb matrix for lbasC
        run_codes(3,epath,'hsfp0_sc','lsxC --job=3',ncore2) # Sigma from core1
        mv_files(glob.glob('stdout*'),'STDOUT')
        ### Valence part of the self-energy Sigma ###
        run_codes(0,epath,'hbasfp0','lbas --job=0')        # Product basis generation
        run_codes(0,epath,'hvccfp0','lvcc --job=0',ncore) # Coulomb matrix for lbas
        run_codes(1,epath,'hsfp0_sc','lsx --job=1',ncore2) # Exchange Sigma
        mv_files(glob.glob('stdout*'),'STDOUT')
        if os.path.isfile('WV.d'):
            rm_files(glob.glob('WV*'))
        # following two runs are most expensive #
        try:
            lx0_para_option
        except NameError:
            lx0_para_option='' #set lx0_para_option='-nq 4 -ns 1'
        run_codes(11,epath,'hx0init','lx0init --job=11',ncore2,lx0_para_option)    #x0 part
        run_codes(100,epath,'hx0zmel','lx0zmel',ncore2,lx0_para_option)   #x0 part
        run_codes(100,epath,'hrcxq','lrcxq',ncore2,lx0_para_option)       #x0 part
        run_codes(100,epath,'hhilbert','lhilbert',ncore2,lx0_para_option) #x0 part

        mv_files(glob.glob('stdout*'),'STDOUT')
        run_codes(2,epath,'hsfp0_sc','lsc --job=2',ncore2) #correlation Sigma
        mv_files(glob.glob('stdout*'),'STDOUT')
        run_codes(0,epath,'hqpe_sc','lqpe')        #all Sigma are combined.
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
            run_codes(100,epath,'lmf-MPIK','llmf_oneshot',ncore,target)
            cp_files(['ctrl.%s'%target,rst_name,sigm_name,'llmf_oneshot','save.%s'%target],run0)

        rundir='RUN.ITER%d'%niter
        gen_dir(rundir)
        cp_files(['ctrl.%s'%target,rst_name,sigm_name,'GWinput','save.%s'%target],rundir)
        print('OK! --> == %d iteration over =='%niter,flush=True)
    ##### end of loop #####

    ### finally we have llmf_gwscend ###
    run_codes(100,epath,'lmf-MPIK','llmf_gwscemd.%d'%niter,ncore,target)
    rm_files(glob.glob('ewindow.%s*')+glob.glob('qbyl.%s*'%target)+glob.glob('eigze*.%s*'%target)+['_IN_'])
    if nloop==0:
        run0='RUN0'
        gen_dir(run0)
        cp_files(['ctrl.%s'%target,rst_name,sigm_name,'llmf_gwscend.%d'%niter,'save.%s'%target],run0)
    argv=''
    print('OK! ==== All calclation finished for  gwsc %s ===='%argv)

if __name__=="__main__":
    main()
