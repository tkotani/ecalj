#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
QSGW calculation with MPI.
'''
import os, datetime, shutil, glob
def gw_args(pname,note):
    '''
    arguments settings
    Returns
      target: target material name ctrl.target
      ncore: number of MPI thereads in lmf
      option: options
    '''
    import argparse
    parser=argparse.ArgumentParser(prog=pname,description=note)
    parser.add_argument("-np",     help='number of mpi cores in lmf',action='store')
    parser.add_argument("material_name",help=':name of extension')
    parser.add_argument('--phispinsym',action='store_true',help='spin-symmetrized augmentation')
#    parser.add_argument('--afsym',action='store_true',help='AF symmetry mode')
    args=parser.parse_args()
    print(args)
    target=args.material_name
    ncore=int(args.np)
    option=''
    if args.phispinsym==True: option=' --phispinsym'
    return(target,ncore,option)

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
