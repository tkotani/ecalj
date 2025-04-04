import sys
import argparse, configparser
from pathlib import Path
from joball import Job

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
parser.add_argument('--inpath', type=Path, default=Path(config.get('inpath')), help='Path to the input directory')
parser.add_argument('--epath', type=Path, default=Path(config.get('epath')), help='Path of ecalj package')
parser.add_argument('--apikey', type=str, default=config.get('apikey'))
parser.add_argument('--nqsub', type=int, default=config.getint('nqsub'))
parser.add_argument('--niter', type=int, default=config.getint('niter'))
parser.add_argument('--ncore', type=int, default=config.getint('ncore'))
parser.add_argument('--bnd4all', type=bool, default=config.getboolean('bnd4all'))
parser.add_argument('--gw80', type=bool, default=config.getboolean('gw80'))
parser.add_argument('--koption', nargs='+', type=int, default=eval(config.get('koption')))
parser.add_argument('--kratio', type=str, default=get_float(config.get('kratio')))
args = parser.parse_args(sys.argv[1:])
print(args)


if __name__=="__main__":

    args.inpath = str(args.inpath.resolve())
    chk_exist(args.inpath)
    chk_exist(args.epath)
    
    ### Run calc.
    job = Job(args)
    job.run()
