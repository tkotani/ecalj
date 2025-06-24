import os,shutil
'''
util for QSGW scripts
'''
def gen_dir(dirname):
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

nxdict = {}
def run_program(cmd, ncore=0, x0=0):
    import subprocess, datetime
    global nxdict
    nx = ncore
    while nx > 0:
        t = datetime.datetime.now()
        if cmd in nxdict and nx==ncore:
            run_cmd = nxdict[cmd]
        else:
            run_cmd = (f'mpirun -np {nx} ' if ncore != 0 else '') + cmd
        print(t if x0 == 0 else t - x0, ' ', run_cmd, flush=True)
        result = subprocess.run(run_cmd, shell=True)
        if result.returncode == 0:
            nxdict[cmd] = run_cmd
            return t
        print(f'Error for ncore={nx} in {run_cmd}', flush=True)
        nx //= 2
    exit(-1)
