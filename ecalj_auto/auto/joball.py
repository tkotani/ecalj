import os, re, shutil, sys
from pathlib import Path


def replace_all(ffile, dic):
    text = Path(ffile).read_text().splitlines()
    text0 = []
    for line in text:
        for i, j in dic.items():
            if i in line:
                default = line.split('default=')[1]
                if default[-1] == ')': default = default[:-1]
                else: default = default.split(',')[0]
                line = line.replace(default, j)
        text0.append(line)
    with Path(ffile).open('w') as f:
        f.write('\n'.join(text0))

class Job:

    def __init__(self, config):
        
        import datetime
        self.xdate = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        self.fname = 'joblist'
        self.lpath = Path(config.inpath) / self.fname
        self.ppath = Path(config.inpath) / 'POSCARALL'
        parts = config.inpath.split('INPUT', 1)
        self.opath = Path(parts[0] + 'OUTPUT' + parts[1], f'start@{self.xdate}')
        self.pythonbin = sys.executable

        config.ppath = self.ppath
        self.config = config
        print(self.xdate)
        print(config.inpath)
        print(self.opath)

    def run(self):
        self.divide_file()
        self.qsub()

    def divide_file(self):

        nqsub = self.config.nqsub
        fpath = self.opath / self.fname
        (self.opath).mkdir(parents=True, exist_ok=True)

        lines = (self.lpath).read_text().splitlines()
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

    def qsub(self):

        shutil.copy(self.lpath, self.opath)
        config = self.config
        
        job_mp = self.opath / 'job_mp.py'
        auto = Path(__file__).resolve().parent
        os.system(f'cp auto/job_mp.py {job_mp}')
        defaults = {
            '--config': "Path(__file__).parent.resolve()/'config.ini'",
            '--poscar': "Path(config.get('ppath'))",
            '--auto': f"Path('{auto}')"
        }
        replace_all(job_mp, defaults)

        max_attr_length = max(len(attr) for attr in config.__dict__.keys())
        config.kratio = round(float(config.kratio), 3)
        with (self.opath / 'config.ini').open('w') as f:
            f.write('[DEFAULT]\n')
            for attr, value in config.__dict__.items():
                f.write(f'{attr.ljust(max_attr_length)} : {value} \n')

        if config.nqsub == 1:
            pattern = self.fname
        if config.nqsub != 1:
            pattern = self.fname + '.[0-9]+'

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
            print()
            print(f, 'ncore=', config.ncore)
            print(f' >>>  {jobx}')
            os.system(ccc)
            os.system(f'qsub {jobx}')
