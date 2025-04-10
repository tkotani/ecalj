#!/usr/bin/env python3
from pathlib import Path
from itertools import groupby
import os, tarfile, sys,argparse

img_form = 'png'
#path_calc = Path('/home/stakano/gw1000/start@20241107-121739/')
path_calc = Path(sys.argv[1]).resolve()
fpath = Path(__file__).absolute().parent

parser = argparse.ArgumentParser(prog='bandpng',description='''summarize band plot to md file''')
parser.add_argument("--qsgwnum",    help='QSGW.?run',default=2,type=int)
args,unknown=parser.parse_known_args()
qsgwnum=args.qsgwnum
print(args)

def exe():
    LDA = GNU('LDA',img_form, path_calc)
    gw = GNU('GW', img_form, path_calc)
    gen_md(Path(f'{path_calc}_band{img_form}'), img_form)


class GNU:
    def __init__(self, calc, img_form, path_calc):
        ### calc is "LDA" or "GW"
        if calc == 'LDA':
            self.dir_calc = 'LDA'
            self.ti = 'LDA'
        if calc == 'GW':
            if(qsgwnum==1): self.dir_calc = 'QSGW.1run'
            if(qsgwnum==2): self.dir_calc = 'QSGW.2run'
            self.ti = 'QSGW'
        self.img_form = img_form
        self.calc = calc
        self.path_calc = path_calc
            
        self.glt = Path(f'{path_calc}_band{img_form}/gnu_{calc}')
        self.img_dir = Path(f'{path_calc}_band{img_form}/{calc}')
        self.glt.mkdir(parents=True ,exist_ok=True)
        self.img_dir.mkdir(parents=True, exist_ok=True)
        
        if img_form == 'eps': self.term = 'postscript'
        elif img_form in ['pdf','png']: self.term = f'{img_form}cairo'
        else:
            print(f'This image format {img_form} is not supported by this code')
            return
        self.main()
        

    def main(self):
        mpids = []
        for f in self.path_calc.glob('mp-*'):
            mpids.append(f)

        for mp in mpids:
            if not (mp/'LDA').exists():
                print(f'{mp.name}: lmf not complete')
                continue
            self.bandplot = mp/f'{self.dir_calc}'/'bandplot.isp1.glt'
            if not self.bandplot.exists():
                print(str(self.bandplot), ' not exists !!!!')
                continue
            print(self.bandplot)
            lines = (mp/'SiteInfo.lmchk').read_text().split('\n')
            atoms = []
            for line in lines:
                if not line: break
                atoms.append(line.split()[3])
            comp = compound(atoms)
            self.mk_glt(mp, comp)
            os.system(f'gnuplot {self.glt}/{mp.name}.glt')
            #print(mp.name)
            #return
        #self.compress_png_files()

    def mk_glt(self, mp, comp):

        mpid = mp.name
        fout = self.glt/f'{mpid}.glt'
        bnd_lines = [f'"{bnd}" u ($2):($3) w l lt 1 lw 1 lc rgb "red" ti "",\\' for bnd in (mp/f'{self.dir_calc}').glob('bnd*.spin*')]

        #TERM OUTPUT TITLE KPATH1 KPATH2 BNDxxx.SPINx TDOS
        dict_replace = {
            'TERM': self.term,
            'OUTPUT': "{}/{}.{}".format(self.img_dir, mpid, self.img_form),
            'TITLE': "[{}] {} {}".format(self.ti, mpid, comp),
            'KPATH': '\n'.join(extract_kpath(self.bandplot)),
            'BND': '\n'.join(bnd_lines),
            'TDOS': (mp/f'{self.dir_calc}')/f'dos.tot.{mpid}'
        }
        gen_gnu(*dict_replace.values(), fout)

    def compress_png_files(self):
        with tarfile.open(f"{self.path_calc}_band{self.img_form}/{self.calc}_{self.img_form}.tar.gz", "w:gz") as tar:
            for _file in Path(self.dir_calc).glob(f"*.{self.img_form}"):
                tar.add(_file, arcname=_file.name)


def compound(lst):
    result = ""
    for key, group in groupby(lst):
        count = len(list(group))
        result += "{}{}".format(key, f"_{count}" if count > 1 else '')
    return result

def extract_kpath(bandplot):
    bandplot = bandplot.read_text().split('\n')
    flag = False
    kpath = []
    for line in bandplot:
        if 'xtics' in line:
            line = line.replace('set xtics', '')
            flag = True
        if 'tpia' in line: break
        if flag:
            if 'GAMMA' in line: line = line.replace('GAMMA', 'Γ')
            if 'SIGMA' in line: line = line.replace('SIGMA', 'Σ')
            if 'DELTA' in line: line = line.replace('DELTA', 'Δ')
            kpath.append(line)
    return kpath
    #kpathを読み取る

def gen_md(img_dir, img_form):
    LDA_dir = img_dir / "LDA"
    content = "## Comparison of LDA & QSGW with respect to band structure \n"
    content += "### Table of contents \n - [Band structure & total DOS](#band-structure--total-dos) \n - [Jump links](#jump-links)\n\n"
    content += "### Band structure & total DOS\n"

    def get_mpid(x):
        return x.stem.split("-")[1].split(".")[0]
    sorted_png_files = sorted(LDA_dir.glob(f"*.{img_form}"), key=lambda x: int(get_mpid(x)))

    for png_file in sorted_png_files:
        fname = png_file.name
        img_number = get_mpid(png_file)
        content += f'<a id="img_{img_number}"></a>'
        content += f'<img src="LDA/{fname}" width="50%">'
        if (img_dir / "GW" / (fname)).exists():
            content += f'<img src="GW/{fname}" width="50%">\n'
        else:
            content += '\n'

    content += "\n### jump links\n"
    sorted_png_files_by_number = sorted(LDA_dir.glob(f"*.{img_form}"), key=lambda x: int(get_mpid(x)))
    for png_file in sorted_png_files_by_number:
        img_number = get_mpid(png_file)
        content += f'[{img_number}](#img_{img_number}) / '

    output_path = img_dir / "bandpng.md"
    output_path.write_text(content, encoding="utf-8")

def gen_gnu(term, output, title, kpath, bnd, tdos, fout):
    content = f"""set term {term} enhanced size 9.5in, 5in
set output "{output}"
unset key
set multiplot
set xzeroaxis
set yzeroaxis
set grid

# # # Variables # # #
# # # setting font of letters # # #
set title font "Helvetica,25" offset 0,-0.7
set xlabel font "Helvetica,20"
set ylabel font "Helvetica,20" offset 1,0
set xtics font "Helvetica,15"
set ytics font "Helvetica,15"

# # # band plot majority # # # 
set border lw 1
set ylabel "Energy (eV)"
set ytics 4
set title "{title}"
set size 0.7 , 0.5
set yrange [   -10.00000:    14.00000]
set tmargin screen 0.85
set bmargin screen 0.15
set lmargin screen 0.10
set rmargin screen 0.80
set xtics {kpath}
tpia=2*3.1415926/     1.88973

plot \
{bnd}


unset xtics
unset ytics

# # # DOS plot # # #
set border lw 1
set size 0.15 , 0.5
set title "DOS"
unset y2label
unset grid
set xtics 20
set yrange [   -10.00000:    14.00000]
set xrange [0:15]
set tmargin screen 0.85
set bmargin screen 0.15
set lmargin screen 0.8
set rmargin screen 0.95
unset ylabel

# # # setting font of letters # # #
set title font "Helvetica,25" offset 0,-0.7
set xlabel font "Helvetica,25"
set ylabel font "Helvetica,25"
set xtics font "Helvetica,25"
set noytics
set noxtics

set style fill transparent solid 0.5 noborder

rydberg=13.605
plot \
"{tdos}" u ($2/rydberg):($1*rydberg) w l lw 1 lc rgb "black"
"""
    Path(fout).write_text(content, encoding="utf-8")


if __name__=='__main__':
    exe()
