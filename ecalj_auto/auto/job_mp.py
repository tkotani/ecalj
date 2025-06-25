#originally shota takano 2025
import os,sys,shutil,re
from pathlib import Path
import pandas as pd
import datetime
from contextlib import contextmanager
import argparse, configparser

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)
sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', buffering=1)
sys.stdin = os.fdopen(sys.stdin.fileno(), 'r', buffering=1)

def get_float(value):
    if '/' in value:
        numerator, denominator = value.split('/')
        return float(numerator) / float(denominator)
    else:
        return float(value)

#parser = argparse.ArgumentParser(description="Run QSGW")
#parser.add_argument('--config', type=Path, default=Path(__file__).resolve().parent/'config.ini', help='Path to the config.ini')
#print('ppp path=',Path(__file__).parent.resolve()/'config.ini')
#parser.add_argument('--config', type=Path, default=opath/'config.ini')
#args, remaining_argv = parser.parse_known_args()

'''python script for qsub. jobx.foobar script calls this script.
   This is called once from each qsub process. As shown in the schedule in joblist.*,
   we run calculations for each material.
'''

class Args:
    pass
args = Args()

# Readin from augments
parser = argparse.ArgumentParser(description="Run QSGW")
parser.add_argument('--dir', type=Path, default=Path(__file__).resolve().parent, help='Output directory')
parser.add_argument('--file', type=str, default=None, help='Listfile of mpid. only 1 file is allowed')
aarg=parser.parse_known_args()
args.file=aarg[0].file
args.dir=aarg[0].dir
print('args.file=',args.file)
print('args.dir=',args.dir)

# Readin from config.ini
config = configparser.ConfigParser()
config_path = args.dir / 'config.ini'
config.read(config_path)
config = config['RUNjob'] 
args.epath = Path(config.get('epath'))
args.poscar = Path(config.get('ppath'))
args.autopath = Path(config.get('autopath'))
args.apikey = config.get('apikey')
args.niter = config.getint('niter')
args.ncore = config.getint('ncore')
args.bnd4all = config.getboolean('bnd4all')
args.gw80 = config.getboolean('gw80')
args.koption = config.getint('koption')
args.kratio = round(float(config.get('kratio')), 8)
args.kkmesh = config.get('kkmesh')
args.ldaonly = config.getboolean('ldaonly') 
import ast
if(config.get('kkmesh') is not None):args.kkmesh=ast.literal_eval(config.get('kkmesh'))
print('kkmesh=',args.kkmesh,type(args.kkmesh))
print(config.get('epath'))

# parser.add_argument('--epath', type=Path, default=Path(config.get('epath')), help='Path of ecalj package')
# parser.add_argument('--dir', type=Path, default=Path(__file__).resolve().parent, help='Output directory')
# parser.add_argument('--poscar', type=Path, default=Path(config.get('ppath')))
# parser.add_argument('--file', type=str, default=None, help='Listfile of mpid. only 1 file is allowed')
# parser.add_argument('--autopath', type=Path, default=config.get('autopath'))
# parser.add_argument('--apikey', type=str, default=config.get('apikey'))
# parser.add_argument('--niter', type=int, default=config.getint('niter'))
# parser.add_argument('--ncore', type=int, default=config.getint('ncore'))
# parser.add_argument('--bnd4all', type=bool, default=config.getboolean('bnd4all'))
# parser.add_argument('--gw80', type=bool, default=config.getboolean('gw80'))
# parser.add_argument('--koption', nargs='+', type=int, default=eval(config.get('koption')))
# parser.add_argument('--kratio', type=float,        default=get_float(config.get('kratio')))
# parser.add_argument('--kkmesh', type=int, nargs=6, default=kkmesh)
# #parser.add_argument('--mpid', type=str, nargs='+', default=config.getint('mpid'))
# #parser.add_argument('--lmxa6', type=bool, default=False)
# #args = parser.parse_args(sys.argv[1:]) #read auto directory.
# print('args=',sys.argv[1:])
# print(args)

sys.path.append(str(args.autopath))
import creplot

args.kratio = round(float(args.kratio), 8)
args.dir = args.dir.resolve()

#if args.mpid:
#    args.mpid = [f for f in args.mpid if re.fullmatch(r'\d+', f)]
if args.file:
    path_log = args.dir / f'log.{args.file}'
    path_lst = args.dir / args.file
else:
    path_log = args.dir / 'log.joblist'


### Set joblist
job_dict = {}
if args.file:
    with path_lst.open('r') as f:
        for line in f:
            if len(line.strip()) == 0:  continue
            if line.split()[0]=='#': continue #skip comment lines
            mag = line.split()[-1].rstrip()
            if mag not in ['NM','FM','AFM','FiM']:
                mag = ''
            job_dict[line.split()[0]] = mag
                
    if path_log.exists():
        try:
            log_df = pd.read_csv(path_log, header=None, delim_whitespace=True, usecols=[0,1,8])
            log_fin = log_df[log_df[1] == 'GW'][0].astype(str).tolist()  # GW calculations have already been done, regardless of convergence
            log_error = log_df[log_df[8] != 'c'][0].astype(str).tolist() # errors during the calculation of LDA/GW

            num_remove = [num for num in job_dict if (num in log_fin) or (num in log_error)]
            for num in num_remove:
                del job_dict[num]
        except pd.errors.EmptyDataError:
            print('Error: log file is empty! skip reading')
#if args.mpid:
#    for mpid in args.mpid:
#        job_dict[mpid] = 'NM'
joblist = list(job_dict.keys())
print('joblist=',joblist)

#print('lmxa6= ',args.lmxa6)
print('job_mp core=',args.ncore)
fout=path_log.open('a')

### errcodes for various outputs
# recalculate if code = 1
dict_errcode = {
    "conv=x": 1,    # not convereged, sometimes because of the lack of k-points
    "kmesh": 1,     # for error in lmf; kmesh-incompatible
    "job_band": 1,  # for error in job_band (maybe lack of k-points)
    "conv=c": 0,    # convereged
    "conv=i": 0,    # normal exit is "c" or "x". "i" means something strange
    "ctrlgen": 0,   # for error in ctrlgenM1.py
    "others": 0     # for other errors
}


def save_result(i, outc, calc, start):
    endtime=datetime.datetime.now()
    difft=int((endtime-start).total_seconds())
    print(i,calc,'start@',starttime.strftime('%x %X'),':',difft,' sec',outc)
    print(i,calc,'start@',starttime.strftime('%x %X'),':',difft,' sec',outc, file=fout,end='\n',flush=True)

@contextmanager
def change_directory(path):
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(args.dir)
    
for i in joblist:
    num = i #    num = 'mp-'+i
    print()
    print('MATERIAL=', num)
    num_dir = args.dir / num
    num_dir.mkdir(exist_ok=True)
    starttime=datetime.datetime.now()
    ordering = job_dict.get(i, '')
    koption = args.koption
    print(' koption=',koption)    #kitmx=3    #for kadd in range(kitmx): # k point choices. Need fixing.
    calc = creplot.Calc(num,args.epath,args.ncore)
    nitmax=40
    with change_directory(num_dir): #go into material directory
        ### LDA
        if Path('LDA').exists():
            print('LDA calc. is already done! skip LDA calc. and go to QSGW calc.')
            pass
        else:
            out_LDA, errcode = calc.run_LDA(args,args.apikey,koption,nitmax,ordering,args.poscar,dict_errcode)
            # lmxa6 removed.
            if errcode == 1: # something wrong in k-mesh. run calc. again
                #if kadd== kitmx: # if the last k-mesh, save the result and go next
                #    save_result(i, out_LDA, 'LDA', starttime)
                continue
            elif errcode == 0: # LDA converged or something wrong in calc. code
                save_result(i, out_LDA, 'LDA', starttime)
                if not out_LDA.startswith('c'):  # something wrong in the code. go next
                    break
        if(args.ldaonly): continue
        
        ### QSGW
        starttime = datetime.datetime.now()
        out_GW = calc.run_QSGW(args,args.niter, args.bnd4all, args.gw80, args.kratio)
        save_result(i, out_GW, 'GW', starttime)
        gap_GW = calc.gap_GW
        if gap_GW is None: pass
        elif gap_GW < 1e-5:
            if calc.gap_LDA > 1e-4:
                print('Run gwsc 3 more loop')
                starttime = datetime.datetime.now()
                out_GW = calc.run_QSGW(3, args.bnd4all, args.gw80, args.kratio)
                save_result(i, out_GW, 'RepeatGW', starttime)
        if not out_GW.startswith('x'): continue

        ### error handling with x ### Not checked well
        # 'x': maybe lack of k-mesh, then run LDA/GW with more k-volume
        dir_save = '_'.join(map(str, calc.k_points))
        if Path(dir_save).exists():
            shutil.rmtree(dir_save)
            shutil.copytree('.', dir_save)
        if Path('osgw.out').exists():
            shutil.copy('osgw.out', dir_save)
            print(f'information are stored in {dir_save}')
        # clear files which are already stored in {dir_save}
        for f in Path('.').glob('*'):
            if f.is_file():
                f.unlink()
        for f in Path('LDA').glob('*'):
            if not f.match('rst.*'):
                shutil.copy(f, Path('.'))
                shutil.rmtree('LDA')
        for run_iter_dir in Path('.').glob('QSGW.*run'):
            if run_iter_dir.is_dir():
                shutil.rmtree(run_iter_dir)
        continue
            
fout.close()
