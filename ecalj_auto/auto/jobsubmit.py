import sys
import argparse, configparser
from pathlib import Path
from Job import Job

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
#parser.add_argument('--lmxa6', type=bool,    default=False, help='Run QSGW80 or not')
parser.add_argument('--koption', nargs='+', type=int, default=eval(config.get('koption')), help='number of k-points option for QSGW calculation')
parser.add_argument('--kratio', type=str,   default=get_float(config.get('kratio')), help='ratio of k point for self energy/band calculation')
parser.add_argument('--kkmesh', type=int, nargs=6, default=None, help='Given k mesh LDA and QSGW 3+3 int numbers ')
args = parser.parse_args(sys.argv[1:])
print(args)
if __name__=="__main__":
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
    job = Job(args)
    job.run()
