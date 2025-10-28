#!/usr/bin/env python3
import os, subprocess, sys, shutil,time,argparse,re

parser = argparse.ArgumentParser(prog='InstallAll',description=
'''
Install ecalj and run tests. Instead of InstallAll, we will use InstallAll.py
''')
parser.add_argument("-np",    help='number of mpi cores for install test',default=8,type=int)
parser.add_argument('--clean',help='Clean CMakeCache CMakeFiles before make',action='store_true')
parser.add_argument('--gpu'  ,help='nvfortran for GPU',action='store_true')
parser.add_argument('--bindir' ,help='ecalj binaries and scripts',type=str,default='bin')
parser.add_argument('--fc'   ,help='fortran compilar  gfortran/ifort/ifx/nvfortran',type=str,required=True)
parser.add_argument('--notest' ,help='no test. only compile',action='store_true')
parser.add_argument('--verbose' ,help='verbose on for debug',action='store_true')
parser.add_argument('--debug' ,help='debug',action='store_true')
args=parser.parse_args()

def main():
    if(args.debug):
        BUILD_TYPE = "Debug"    # = "Debug"
    else:
        BUILD_TYPE = "Release"    # = "Debug"
    CWD =os.getcwd()
    HOME=os.getenv('HOME')
    BINDIR = os.path.join(HOME,args.bindir) #os.path.join(HOME, 'bin')  # Make directory for ecalj binaries and scripts.
    ncore=args.np
    #FC = os.getenv('FC')
    #if not FC:
    FC=args.fc
    verbose=''
    if(args.verbose): verbose='VERBOSE=1 '
    #    if(FC==''):
    #        print('Usage: >FC=gfortran ./InstallAll [options]. Run ./InstallAll -h for help.')
    #        sys.exit()
    if not os.path.exists(BINDIR): os.makedirs(BINDIR)
    print(f"Going to install required binaries and scripts to {BINDIR}")
    start0_time = time.time()
    # Make links
    EXECDIR = os.path.join(CWD,'SRC/exec')
    BUILDIR = os.path.join(CWD,'SRC/exec/build')
    print(EXECDIR)
    for scr in ['StructureTool/viewvesta', 'StructureTool/ctrl2vasp', 'StructureTool/vasp2ctrl','GetSyml/getsyml']:
        src   = os.path.join(CWD, scr+'.py')
        slink  = os.path.join(BINDIR,  scr.split('/')[-1])
        slink2 = os.path.join(EXECDIR, scr.split('/')[-1])
        if os.path.exists(os.path.join(BINDIR, slink)):
            os.remove(os.path.join(BINDIR,  slink))
        if os.path.exists(os.path.join(BINDIR, slink2)):
            os.remove(os.path.join(EXECDIR, slink2))
        print(f"ln -s {src} {slink}",f"ln -s {src} {slink2}")
        os.symlink(src, slink)
        os.symlink(src, slink2)
    if(args.clean):
        # Clean CMakeCache & CMakeFiles in exec
        os.path.exists(f'{EXECDIR}/CMakeCache.txt') and os.remove(f'{EXECDIR}/CMakeCache.txt')
        shutil.rmtree(f'{EXECDIR}/CMakeFiles', ignore_errors=True)
        shutil.rmtree(f'{BUILDIR}', ignore_errors=True)
    os.makedirs(f'{BUILDIR}', exist_ok=True)
    os.chdir(f'{BUILDIR}')
    if(args.gpu): #Obata for nvfortran
        if os.system(f'FC={FC} cmake .. -DBUILD_MP=ON -DBUILD_GPU=ON -DBUILD_MP_GPU=ON -DCMAKE_BUILD_TYPE={BUILD_TYPE}') != 0:sys.exit(1)
    elif(FC in ["gfortran", "ifort", "ifx", "nvfortran"]):
        if os.system(f'FC={FC} cmake .. -DCMAKE_BUILD_TYPE={BUILD_TYPE}') != 0: sys.exit(1)
    else:
        print('Check InstallAll')
        sys.exit(-1)
    jobs = min(os.cpu_count(), 32)
    if os.system(f'{verbose}make -j {jobs}') != 0: sys.exit(1)
    # Copy executables to BINDIR (but not soft link)
    executables = [
        os.path.join(d, f)
        for d in [BUILDIR, EXECDIR]
        for f in os.listdir(d)
        if os.path.isfile(os.path.join(d, f))
        and not os.path.islink(os.path.join(d, f))
        and os.access(os.path.join(d, f), os.X_OK)
    ]
    print('COPY to BINDIR',executables)
    for exe in executables:
        shutil.copy(exe, BINDIR)
    if(args.notest):
        print('No test. Only compile.')
        return
    
    # Install test
    print()
    print('=== goto test ===')
    os.chdir(f'{CWD}/Samples/TestInstall')
    start_time = time.time()
    os.system(f'{BINDIR}/testecalj -np {ncore} --all')
    end_time = time.time()
    elapsed_time = end_time - start_time
    elapsed0_time = start_time-start0_time
    print(f"Elapsed time for make        : {elapsed0_time:.0f} seconds")
    print(f"Elapsed time for testecalj.py: {elapsed_time:.0f} seconds")

if __name__ == "__main__":
    main()
