#!/usr/bin/env python3
import os, subprocess, sys, shutil,time,argparse,re

parser = argparse.ArgumentParser(prog='InstallAll',description=
'''
Install ecalj and tests. Instead of InstallAll, we will use InstallAll.py
''')
parser.add_argument("-np",    help='number of mpi cores',default=8,type=int)
parser.add_argument('--clean',help='Clean CMakeCache CMakeFiles before make',action='store_true')
parser.add_argument('--gpu'  ,help='nvfortran for GPU',action='store_true')
parser.add_argument('--bindir' ,help='ecalj binaries and scripts',type=str,default='bin')
parser.add_argument('--fc'   ,help='fortran compilar  gfortran/ifort/nvfortran',type=str,required=True)
args=parser.parse_args()

def main():
    BUILD_TYPE = "Release"    # = "Debug"
    CWD =os.getcwd()
    HOME=os.getenv('HOME')
    BINDIR = os.path.join(HOME,args.bindir) #os.path.join(HOME, 'bin')  # Make directory for ecalj binaries and scripts.
    ncore=args.np
    #FC = os.getenv('FC')
    #if not FC:
    FC=args.fc
    #    if(FC==''):
    #        print('Usage: >FC=gfortran ./InstallAll [options]. Run ./InstallAll -h for help.')
    #        sys.exit()
    if not os.path.exists(BINDIR): os.makedirs(BINDIR)
    print(f"Going to install required binaries and scripts to {BINDIR}")
    start0_time = time.time()
    # Make links
    for scr in ['StructureTool/viewvesta', 'StructureTool/ctrl2vasp', 'StructureTool/vasp2ctrl','GetSyml/getsyml']:
        src   = os.path.join(CWD, scr+'.py')
        slink = os.path.join(BINDIR, scr.split('/')[-1])
        if os.path.exists(os.path.join(BINDIR, slink)):
            os.remove(os.path.join(BINDIR, slink))
        print(f"ln -s {src} {slink}")
        if os.path.islink(slink): os.remove(slink)
        os.symlink(src, slink)
    # Make executables
    os.chdir(f'{CWD}/SRC/exec')
    if(args.clean): os.system('rm -rf CMakeFiles CMakeCache.txt')
    if(args.gpu): #Obata for nvfortran
        if os.system(f'FC={FC} cmake . -DBUILD_MP=ON -DBUILD_GPU=ON -DBUILD_MP_GPU=ON -DCMAKE_BUILD_TYPE={BUILD_TYPE}') != 0:sys.exit(1)
        if os.system('make -j 32') != 0: 
            if os.system('make -j 32') != 0: sys.exit(1)
    elif(FC in ["gfortran", "ifort", "nvfortran"]):
        if os.system(f'FC={FC} cmake .') != 0: sys.exit(1)
        if os.system('make -j')          != 0: sys.exit(1)
    else:
        print('Check InstallAll')
        sys.exit(-1)
    # Copy executables to BINDIR
    executables = [f for f in os.listdir('.') if os.path.isfile(f) and os.access(f, os.X_OK)]
    print(executables)
    for exe in executables:
        shutil.copy(exe, BINDIR)
    # Install test
    print('=== goto test ===')
    os.chdir(f'{CWD}/SRC/TestInstall')
    start_time = time.time()
    os.system(f'./testecalj.py -np {ncore} ')
    end_time = time.time()
    elapsed_time = end_time - start_time
    elapsed0_time = start_time-start0_time
    print(f"Elapsed time for make        : {elapsed0_time:.0f} seconds")
    print(f"Elapsed time for testecalj.py: {elapsed_time:.0f} seconds")

if __name__ == "__main__":
    main()
