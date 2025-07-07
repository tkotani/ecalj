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

const_b={}
def run_program_breduction(cmd, ncore=0, x0=0, ext=''):
    """Run a program with b-mix reduction.
       Decreasing b until lmf converges."""
    import subprocess, datetime,shutil,glob,re,os
    global const_b
    with open(f'ctrl.{ext}', 'r') as f:
        text = f.read()
        bval = float(re.search(r'b=0?\.\d+', text).group().split('=')[1])
        
    while ( bval > 0.05): #bmix reducing
        t = datetime.datetime.now()
        if os.path.isfile(f'rst.{ext}'):
            shutil.copy(f'rst.{ext}', f'rst.{ext}.bk')
        if cmd in const_b: bval= const_b[cmd]
        run_cmd = (f'mpirun -np {ncore} ' if ncore != 0 else '') + cmd
        print(t if x0 == 0 else t - x0, ' ', run_cmd, flush=True)

        #result = subprocess.run(run_cmd, shell=True)
        result = subprocess.run(run_cmd, shell=True, stderr=subprocess.PIPE)
        if result.returncode != 0:
            err_lines = result.stderr.decode().splitlines()
            if err_lines: print(err_lines[0], flush=True)
        with open(f'save.{ext}', 'r') as f:
            lines = f.readlines()
            last_line=lines[-1].strip()
            first_word = last_line.split()[0]
            if first_word == 'c' or first_word =='x': 
                const_b[cmd] = bval
                return t
        print(f'Error for const_b={bval} in {run_cmd}', flush=True)
        for file in glob.glob("mix*"): os.remove(file)
        os.remove(f'rst.{ext}')
        shutil.copy(f'rst.{ext}.bk',f'rst.{ext}')

        bval = round(bval - 0.05, 2)
        print(f'Set b={bval} in {run_cmd}', flush=True)
        with open(f'ctrl.{ext}', 'r') as f:
            text = f.read()
        text_new = re.sub(r'b=0?\.\d+', f'b={bval}', text)
        with open(f'ctrl.{ext}', 'w') as f:
            f.write(text_new)
        
    exit(-1)
