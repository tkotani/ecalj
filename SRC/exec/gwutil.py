#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
QSGW calculation with MPI.
'''
import os, datetime, shutil, glob
def gwsc_args():
    '''
    arguments settings
    Returns
      target: target material name ctrl.target
      nloop: numbdr of QSGW iterations starting from current result (nloop=0 is replaced by nloop=1 internally)
      ncore: number of MPI thereads in lmf
      ncore2: number of MPI thereads in lxsC etc.
      option: options
    '''
    import argparse
    parser=argparse.ArgumentParser(prog='gwsc',description='QSGW calculation')
    parser.add_argument("-np",     help='number of mpi cores in lmf',action='store')
    parser.add_argument("-np2",    help='number of mpi cores in lxc etc.',action='store') 
    parser.add_argument("nloop",   help='iteration number of QSGW loop')
    parser.add_argument("material_name",help='material name')
    parser.add_argument('--phispinsym',action='store_true',help='spin-symmetrized augmentation')
    parser.add_argument('--afsym',action='store_true',help='AF symmetry mode')
    args=parser.parse_args()
    print(args)
    target=args.material_name
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
    option=''
    if args.phispinsym==True: option=' --phispinsym'
    if args.afsym==True: option=option+' --afsym'
    return(target,nloop,ncore,ncore2,option)

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

def run_program(commandline,ncore=0):
    '''
    Run codes:
        mpi_size:
        command line: 
    '''
    import subprocess,datetime
    xdate=datetime.datetime.today().isoformat()
    mpirun='mpirun -np %d '%ncore if ncore!=0 else ''
    run_command = mpirun + commandline
    print(xdate+'  '+run_command,flush=True)
    info=subprocess.run(run_command,shell=True)
    if info.returncode!=0: #if return error
        print('Error in '+run_command,flush=True)
        exit(-1)
