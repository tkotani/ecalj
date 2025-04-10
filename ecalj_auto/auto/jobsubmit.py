import sys
import argparse, configparser
from pathlib import Path
#from Job import Job
import os, re, shutil, sys
from pathlib import Path

def chk_exist(path):
    if not Path(path).exists():
        print(f"Error: The epath '{path}' does not exist.")
        sys.exit()
def get_float(value):
    if '/' in value:
        numerator, denominator = value.split('/')
        return float(numerator) / float(denominator)
    else:
        return float(value)

class JobQUE:
    '''  Divide joblist file and submit jobs using qsub.
    '''
    def __init__(self, config):
        import datetime
        self.xdate = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        self.fname = 'joblist'
        self.autopath = Path(__file__).resolve().parent
        self.lpath = Path(config.inpath) / self.fname
        self.ppath = Path(config.inpath) / 'POSCARALL'
        parts = config.inpath.split('INPUT', 1)
        self.opath = Path(parts[0] + 'OUTPUT' + parts[1], f'start@{self.xdate}') #OUTPUT path
        self.pythonbin = sys.executable

        config.ppath = self.ppath
        self.config = config
        print('date       =',self.xdate)    # datetime stamp
        print('input  path=',config.inpath)
        print('output path=',self.opath)
        print('epath=',Path(config.epath))

    def run(self):
        self.divide_file()
        self.qsuball()

    def divide_file(self): #joblist file is divided into nqsub files
        nqsub = self.config.nqsub
        fpath = self.opath / self.fname
        (self.opath).mkdir(parents=True, exist_ok=True)

        lines = (self.lpath).read_text().splitlines()
        print(f"Reading {self.lpath}")
        
        lines = [line for line in lines if line.strip()]
        lines = [l for l in lines if not (l.startswith('#') or l[0].isspace())]
        if len(lines) < nqsub:
            nqsub = len(lines)
            self.config.nqsub = nqsub

        print(f"{self.fname} is divided into {nqsub} files")
        if nqsub == 1: return

        if len(lines) % nqsub == 0:
            nsize = len(lines) // nqsub
        else:
            nsize = len(lines) // nqsub + 1

        i, fnum = 0, 1
        ffw = Path(f"{fpath}.{fnum}").open("w")
        for line in lines:
            i += 1
            print(line,file=ffw)
            if i % nsize==0 and i != len(lines):
                fnum += 1
                ffw.close()
                ffw = Path(f"{fpath}.{fnum}").open("w")
        ffw.close()

    def qsuball(self): # qsub jobx.1...   We use sed jobtemplate with actual arguments
        shutil.copy(self.lpath, self.opath)
        config = self.config
        job_mp = self.autopath / 'job_mp.py'
        
        max_attr_length = max(len(attr) for attr in config.__dict__.keys())
        config.kratio = round(float(config.kratio),8)
        with (self.opath / 'config.ini').open('w') as f:  # Here we write opath/config.ini from auto/config.ini
            f.write('[RUNjob]\n')
            f.write(f'autopath : {self.autopath} \n')      
            for attr, value in config.__dict__.items():
                if(value == None): continue
                f.write(f'{attr.ljust(max_attr_length)} : {value} \n') 
                
        if config.nqsub == 1:  pattern = self.fname
        if config.nqsub != 1:  pattern = self.fname + '.[0-9]+'

        for i in (self.opath).glob('*'):
            f = i.name
            if not re.fullmatch(pattern, f): continue
            if re.fullmatch(r'\S+\.\d+', f):
                job_num = f.split('.')[-1]
            else:
                job_num = '0'
            job_name = 'job' + job_num
            
            jobx = self.opath / f"jobx.{job_num}"
            fout = self.opath / f"{job_name}.out"
            ftouch = self.opath / f"{job_name}.finished"
            
            ccc = (f'sed -e "s/job_name/{job_name}/" '
                   f'-e "s/ncore/{config.ncore}/" '
                   f'-e "s#pythonbin#{self.pythonbin}#" '
                   f'-e "s#job_mp#{job_mp}#" '
                   f'-e "s/infile/{f}/" '
                   f'-e "s#dir_out#{self.opath}#" '
                   f'-e "s#outfile#{fout}#" '
                   f'-e "s#file_to_touch#{ftouch}#" jobtemplate > {jobx}') 
                   # We obtain qsub scricpt jobx.* from jobtemplate
            print()
            print(f, 'ncore=', config.ncore)
            print(f' >>>  {jobx}')
            os.system(ccc)  #jobx.* is created
            os.system(f'qsub {jobx}') #qsub here
            
    # def replace_all(self,ffile, dic):
    #     text = Path(ffile).read_text().splitlines()
    #     text0 = []
    #     for line in text:
    #         for i, j in dic.items():
    #             if i in line:
    #                 default = line.split('default=')[1]
    #                 if default[-1] == ')': default = default[:-1]
    #                 else: default = default.split(',')[0]
    #                 line = line.replace(default, j)
    #         text0.append(line)
    #     with Path(ffile).open('w') as f:
    #         f.write('\n'.join(text0))


if __name__=="__main__":
    config = configparser.ConfigParser()
    config_path = Path(__file__).resolve().parent.parent / 'config.ini'
    config.read(config_path)
    config = config['DEFAULT']
    parser = argparse.ArgumentParser(description="Run QSGW")
    parser.add_argument('--inpath', type=Path, default=Path(config.get('inpath')), help='Path to the input directory. We automatically have correspoining OUTPUT directory')
    parser.add_argument('--epath', type=Path, default=Path(config.get('epath')), help='Path of ecalj package')
    parser.add_argument('--apikey', type=str,   default=config.get('apikey'), help='API key for Materials Project')
    parser.add_argument('--nqsub', type=int,    default=config.getint('nqsub'), help='Number of qsub process to divide the input file')
    parser.add_argument('--niter', type=int,    default=config.getint('niter'), help='Number of iterations for QSGW calculation')
    parser.add_argument('--ncore', type=int,    default=config.getint('ncore'), help='Number of cores per qsub process ')
    parser.add_argument('--bnd4all', type=bool, default=config.getboolean('bnd4all'), help='Generate band plots or not')
    parser.add_argument('--gw80' , type=bool,    default=False, help='Run QSGW80 or not')
    parser.add_argument('--koption', nargs='+', type=int, default=eval(config.get('koption')), help='number of k-points option for QSGW calculation')
    parser.add_argument('--kratio', type=str,   default=get_float(config.get('kratio')), help='ratio of k point for self energy/band calculation')
    parser.add_argument('--kkmesh', type=int, nargs=6, default=None, help='Given k mesh LDA and QSGW 3+3 int numbers ')
    args = parser.parse_args(sys.argv[1:])
    print(args)
    args.inpath = str(args.inpath.resolve())
    chk_exist(args.inpath)
    chk_exist(args.epath)     #print('job_mp: inpath= ',args.inpath,' epath= ',args.epath)
    jfile= Path(args.inpath) / './joblist' # joblist file exists, or joblist is created, which contains names of POSCAR in POSCARALL/
    pcar = Path(args.inpath) / 'POSCARALL'
    print('joblist exists?  =',jfile.exists()) 
    print('POSCARALL exists?=',pcar.exists()) 
    if not jfile.exists(): # Write joblist fie joblist is not exist.
        print('joblist not exists, create joblist')
        if not pcar.exists():
            print('POSCARALL not exists, please check your input path')
            sys.exit()
        print('POSCARALL exists, create joblist')
        jobl=[item.name for item in pcar.iterdir() if item.is_file()]
        with jfile.open('w') as f:
            for item in jobl:
                f.write(f"{item.strip('POSCAR.')}\n")
    ### Run calc.
    job = JobQUE(args)
    job.run()
