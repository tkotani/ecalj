import sys
ppp=sys.executable

import os, subprocess, re, shutil, time
from pathlib import Path
import change_k
import pandas as pd
from contextlib import contextmanager


@contextmanager
def change_directory(path):
    cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(cwd)

def replace(f, str1, str2, count=1):
    file_path = Path(f)
    content = file_path.read_text()
    new_content = re.sub(str1, str2, content, count=count)
    file_path.write_text(new_content)

def check_save(save, n=1):
    outc, returncode = run_command(['tail',f'-n{n}',f'{save}'], sw_return=True, flush=False)
    if returncode != 0:
        return f'Error! {save} not exist'
    return outc
    
def run_command(command, out='', mode='w', sw_print=False, sw_return=False, flush=True):
    command = [str(c) for c in command]
    if flush: print(f"{' '.join(command)} {f'> {out}' if out else ''}")
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = result.stdout.decode('utf-8')
    if out:
        with Path(out).open(mode) as f:
            f.write(output)
    elif sw_print and output:
        print(output)
        
    if result.returncode != 0:
        print(f"Command '{' '.join(command)}' failed with exit code {result.returncode}")
        if result.stderr.decode('utf-8'):
            print(result.stderr.decode('utf-8'))
    if sw_return:
        return output, result.returncode
        
def run_popen(main_command, sub_command, out, oute, mode):
    main_command = [str(c) for c in main_command]
    sub_command = [str(c) for c in sub_command]
    print(' '.join(main_command), f'> {out}')
    with Path(out).open(mode) as f, Path(oute).open(mode) as ff
        main_process = subprocess.Popen(main_command, stdout=f, stderr=ff)
        sub_process = subprocess.Popen(sub_command)
        try:
            return_code = main_process.wait()
            if return_code != 0:
                print(f"Main command '{' '.join(main_command)}' failed with exit code {return_code}")
                return
        finally:
            time.sleep(1)
            sub_process.terminate()

def run_with_save(command, out, mode):
    command = [str(c) for c in command]
    print(' '.join(command), f'> {out}')
    with subprocess.Popen(command, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as process:
        with Path(out).open(mode, buffering=1) as fout:
            for line in process.stdout:
                print(line, end='')
                fout.write(line)
        return_code = process.wait()
        if return_code != 0:
            print(f"Command '{' '.join(command)}' failed with exit code {return_code}")
            return

class Calc:

    def __init__(self, num, epath, ncore, so=False):
        self.num = num
        self.epath = Path(epath)
        self.ncore = str(ncore)
        self.k_points = None
        self.const_b = 0.2
        self.gap_LDA = None
        self.option_lmf = []
        self.option_bnd = []
        if so: self.option_bnd.extend(['-vnspin=2', '-vso=1'])

    def k_points_from_ctrl(self):
        line, code = run_command(['grep', 'nk1=', f'ctrl.{self.num}'], sw_return=True, flush=False)
        if code != 0:
            print(f'something is wrong with ctrl.{self.num}')
            return
        match = re.search(r'nk1=(\d+) nk2=(\d+) nk3=(\d+)', line)
        nk1, nk2, nk3 = map(int, match.groups())
        self.k_points = [nk1, nk2, nk3]

    def set_bmix(self, path_rst, save='llmf_bmix'):

        def run_lmf_bmix():
            for f in Path('.').glob('mix*'):
                f.unlink()
            if path_rst != '.':
                shutil.copy(f'{path_rst}/rst.{self.num}', f'rst.{self.num}')
            self.const_b = round(self.const_b - 0.05, 2)
            self.option_lmf = [f'-vnb={self.const_b}']
            return self.run_lmf(save,save+'.err')
        
        outc = 'i'
        while (not outc.startswith('c')) and self.const_b > 0.05:
            outc = run_lmf_bmix()
            if outc.split()[0] == 'c':
                replace(f'ctrl.{self.num}', r'b=0?\.\d+', f'b={self.const_b}')
                return outc
        print('Not converged with small bmix')
        return outc

    def run_lmf(self, fout,foute):
        command1 = ['mpirun', '-np', self.ncore, self.epath/'lmf', self.num] + self.option_lmf
        command2 = ['tail', '-f', f'save.{self.num}']
        run_popen(command1, command2, fout, foute, 'a')
        return check_save(f'save.{self.num}')

    def run_gwsc(self, niter, mode):
        gwsc_command = [self.epath/'gwsc', niter, '-np', self.ncore, self.num]
        run_with_save(gwsc_command, 'osgw.out', mode)
        return check_save('osgw.out')
    
    def run_gwsc_re(self, niter):
        iter_files = sorted(Path('.').glob('QSGW.*run'), key=lambda f: int(re.match(r'QSGW.(\d+)run', f.name).group(1)))
        x = len(iter_files)
        if iter_files:
            max_run_iter = iter_files[-1]
            path_rst = max_run_iter
        else:
            path_rst = 'LDA'
        outc = self.set_bmix(path_rst)
        if outc.split()[0] == 'c':
            return self.run_gwsc(niter - x, 'a')
        else:
            return 'ERROR: not converged with smaller bmix'

    def run_plot(self, path_save):

        self.mkplot_and_save(path_save)
        read_gap = ReadBND(path_save)
        if read_gap.Eg is None:
            ## loosen symmetry and calc. again
            if not Path('outlmchk').exists():
                run_command([self.epath/'lmchk', self.num], out='outlmchk')
            lines = Path('outlmchk').read_text().splitlines()
            for i,line in enumerate(lines):
                if 'ig group ops' in line:
                    lsymg = lines[i+1].split()[-1].rstrip('\n') # the loosest symmetry-group
                    break
            replace(f'ctrl.{self.num}', 'SYMGRP find', f'SYMGRP {lsymg}')
            self.mkplot_and_save(path_save)
            replace(f'ctrl.{self.num}', f'SYMGRP {lsymg}', 'SYMGRP find')
            read_gap = ReadBND(path_save)
        return read_gap.Eg

    def mkplot_and_save(self, path_save):

        print(path_save)
        if path_save == 'LDA':
            Path(path_save).mkdir(parents=True, exist_ok=True)
            for f in Path('.').glob(f'*.{self.num}'):
                try:
                    if f.exists():
                        shutil.copy(f, path_save)
                except:  # maybe sigm.{num}
                    continue
            shutil.copy('POSCAR', path_save)

        with change_directory(path_save):
            #run_command([self.epath/'lmfa', self.num], out='llmfa', mode='a')
            #run_command(['mpirun', '-np', self.ncore, self.epath/'lmf', self.num] + self.option_bnd, out='llmf', mode='a')
            run_command([ppp, self.epath/'getsyml', self.num, '--nobzview'], out='lgetsyml')
            for job in ['job_band', 'job_tdos', 'job_pdos']:
                file_out = 'l' + job.replace('_', '')
                jjob = [self.epath/job, self.num, '-np', self.ncore, 'NoGnuplot'] + self.option_bnd
                run_command(jjob, out=file_out)

    def run_LDA(self, key, lmxa6, option, ordering, path_poscar, errcode):
        num = self.num
        epath = self.epath
        print(Path.cwd())
        
        ### lmfa
        if Path(f'atm.{num}').exists(): #if already lmfa is finished
            pass
        else:
            ### get POSCAR, ctrls.{num}
            try:
                path_poscar = path_poscar / f'POSCAR.{num}'
                shutil.copy(path_poscar, 'POSCAR')
            except: 
                from pymatgen.ext.matproj import MPRester
                with MPRester(key) as mpr:
                    struc = mpr.get_structure_by_material_id(num)
                    struc.to(fmt='poscar', filename='POSCAR')
            run_command([epath/'vasp2ctrl', 'POSCAR'], out='llmf')
            shutil.copy('ctrls.POSCAR.vasp2ctrl', f'ctrls.{num}')

            ### Set magnetic conditions if ferro, anti-ferro or ferri
            print(ordering)
            if ordering == ('NM' or '') : pass
            else:
                import mag
                mag.SetMag(key, num, 'POSCAR', 3, ordering)

            ### Get ctrl.{num}
            nsp = '1' if ordering=='NM' else '2'
            run_command([epath/'ctrlgenM1.py', num, f'--nspin={nsp}'], out='llmf', mode='a')
            if not Path(f'ctrlgenM1.ctrl.{num}').exists():
                print(f'ctrl.{num} doesn\'t exist !!!')
                return 'ERROR: ctrlgen', errcode['ctrlgen']
            shutil.copy(f'ctrlgenM1.ctrl.{num}', f'ctrl.{num}')
            if lmxa6:
                replace(f'ctrl.{num}', 'LMXA=4', 'LMXA=6')
            
            ### Run lmchk
            run_command([epath/'lmchk', num], out='llmchk')
            run_command(['grep', 'conf', 'llmchk'], sw_print=True, flush=False)

            ### Run lmfa
            run_command([epath/'lmfa', num], out='llmfa')

        ### Get k_mesh from input option
        self.k_points = change_k.get_kpoints([option[1]] * 3)
        if self.k_points is None:
            return 'ERROR: PlatQplat.chk in lmfa', errcode['ctrl']
        self.option_lmf = '-vnit={} -vnk1={} -vnk2={} -vnk3={}'.format(option[0], *self.k_points).split() + option[2:]
        kkk = 'nk1={} nk2={} nk3={}'.format(*self.k_points)

        ### Run lmf
        outc = self.run_lmf('llmf','llmf.err')
        print('lmf finished')

        ### Check convergence and if k-mesh is appropreate or not
        conv = outc.split()[0]
        if conv == 'x': return f'x {kkk}', errcode['conv=x']
        elif conv == 'i': return f'i {kkk}', errcode['conv=i']
        elif conv != 'c':
            #errmsgs = check_save('llmf', n=3)
            #!for msg in errmsgs.split('\n'):
            #   if ('incompatible with this mesh' in msg) or ('qp mapped to is not on k-mesh' in msg):
            if(os.path.exists('bzmesh.err')):
                os.remove('bzmesh.err')
                print('Run lmf')
                self.k_points = [option[1]] * 3
                kkk = 'nk1={} nk2={} nk3={}'.format(*self.k_points)
                self.option_lmf = '-vnit={} -vnk1={} -vnk2={} -vnk3={}'.format(option[0], *self.k_points).split() + option[2:]
                outc = self.run_lmf('llmf','llmf.err')
                conv = outc.split()[0]
                if conv == 'c':
                    self.kmesh_cube = True
                    break
                return f'{conv} {kkk}', errcode['kmesh']
            else:
                return 'ERROR: lmf', errcode['others']

        ### Make plots
        replace(f'ctrl.{num}', r'nk1=(\d+) nk2=(\d+) nk3=(\d+)', kkk)
        self.gap_LDA = self.run_plot('LDA')
        if self.gap_LDA is None:
            return 'ERROR: job_band', errcode['job_band']
    
        return f'c {kkk} gap={self.gap_LDA}', errcode['conv=c']


    def run_QSGW(self, niter, bnd4all, gw80, k_gw):
        num = self.num
        epath = self.epath
        self.gap_GW = None
        if gw80:
            self.option_bnd.append('-vssig=0.8')
        if not self.k_points:
            self.k_points_from_ctrl()
            if not self.k_points:
                return 'ERROR: something wrong in ctrl'
        if not self.gap_LDA:
            self.gap_LDA = self.run_plot('LDA')
        q_points = change_k.get_q(self.k_points, k_gw)   #k_points is int list, like [6, 6, 6]
        print('q-mesh', q_points)
        
        ### Prepare input for QSGW calc.
        replace(f'ctrl.{num}', r'nit=\d+', 'nit=80')
        run_command([epath/'mkGWinput', num], out='lgwin')
        shutil.copy('GWinput.tmp', 'GWinput')
        replace('GWinput', r'n1n2n3\s+\d+\s+\d+\s+\d+', 'n1n2n3   {}   {}   {}'.format(*q_points))

        ### Run n-shot QSGW
        out_gwsc = self.run_gwsc(niter, 'w')
        outc = check_save(f'save.{num}')

        def check_gw_out(osgw):
            log = osgw.split('>')[-1].strip()
            return f'ERROR: {log} in gwsc'
        
        ### Check if QSGW calc. has been successfully completed or not
        if 'Error' in out_gwsc:
            err_msg = check_gw_out(out_gwsc)
            if 'lmf' not in out_gwsc:
                return err_msg
            if outc.split()[0] == 'i':   # lmf error in gwsc loop
                while (outc.split()[0] != 'c') and (self.const_b > 0.05):
                    out_gwsc = self.run_gwsc_re(niter)
                    outc = check_save(f'save.{num}')
                    if 'ERROR' in out_gwsc: return out_gwsc
                    if (outc.split()[0] == 'c' ) and ('Error' not in out_gwsc):
                        break
                else:
                    return check_gw_out(out_gwsc)
                
        ### Check if calc. is converged or not
        if outc.split()[0] == 'x':
            outc = self.set_bmix('.', save=f'llmf_gwscemd.{niter-1}')
        if outc and (outc.split()[0] not in ['c','x']):
            if 'Start' in outc:
                xxx = outc.split()
                return f"ERROR: {xxx[xxx.index('Start') + 1]} in gwsc"
            return outc

        ### Save QPU / Plot Band for all iterations of QSGW calc., swtich is "bnd4all"
        dirs = [d for d in Path('.').glob('QSGW.*run') if d.is_dir()]
        dirs = sorted(dirs, key=lambda d: int(re.match(r'QSGW.(\d+)run', d.name).group(1)))
        for n, d in enumerate(dirs):
            n = n+1
            if (d/'ljobband').exists(): continue
            shutil.copy(f'QPU.{n}run', d)
            print(f'creplot: QPU.{n}run to {d}')
            if bnd4all or (n == niter - 1):
                self.gap_GW = self.run_plot(d)
                if self.gap_GW is None:
                    return 'ERROR: job_band'

        run_command([epath/'cleargw','.'], out='clear')
        return f'{outc.split()[0]} gap={self.gap_GW}'



#### Now this is only for nonmagnetic system...
class ReadBND:
    
    def __init__(self, path, nr=4, tolerance=1e-4, nsp=1, gap0=True):
        "Read Effective Mass from bnd*** files"
        self.path = Path(path)
        self.nr = nr
        self.tol = tolerance
        self.gap0 = gap0
        self.nsp = nsp
        self.bnds = list(self.path.glob('bnd*.spin*'))
        self.onsite = None
        self.Eg = None
        self.Ecbm = None
        self.Evbm = None
        self.main()

    def main(self):
        ### check bnd***.spin*
        if not self.bnds: return

        ### get the number of onsite electrons                                                                                                                                                                    
        self.parse_efermi_lmf(self.path / "efermi.lmf")
        if self.onsite is None: return

        ### read and concatenate bnd files
        self.df_conduction = self.concatenate_bnd(self.onsite+1)
        self.df_valence = self.concatenate_bnd(self.onsite)

        ### get k-points of CBM from number of electrons
        df_CBM, self.Ecbm = self.get_df_extremum("min")
        df_VBM, self.Evbm = self.get_df_extremum("max")
        self.k_points_CBM = self.get_k_dict(df_CBM)
        self.k_points_VBM = self.get_k_dict(df_VBM)

        ### get band gap and transition type
        self.Eg = round(self.Ecbm - self.Evbm, self.nr)
        if self.Eg <= 0: self.metal = True
        if self.gap0: self.Eg = max(self.Eg, 0)
        self.transition = "D" if any(k_v in self.k_points_CBM.values() for k_v in self.k_points_VBM.values()) else "I"

    def parse_efermi_lmf(self, efermi_path):
        def get_energy(line):
            return float(re.sub('D', 'E', line.split()[0]))

        lines = efermi_path.read_text().splitlines()
        if len(lines) < 4: return

        efermi = get_energy(lines[0])
        valtop = get_energy(lines[1])
        conbtm = get_energy(lines[2])
        self.metal = (valtop > efermi) or ((conbtm - valtop) <= 0)

        ele = re.findall(r'\d+\.\w+\+\d+', lines[3])[0]
        onsite = round(float(re.sub('D', 'E', ele)) / 2, 3)
        if not onsite.is_integer():
            print('Number of electrons is odd.')
            print('This should be calculated as a magnetic compound')
            self.onsite = None
        else:
            self.onsite = int(onsite)
            if self.nsp == 2: self.onsite *= 2

    def concatenate_bnd(self, num):
        df_list = [self.load_data(bnd, num) for bnd in self.bnds]
        df_list = [df for df in df_list if df is not None]

        for i, df in enumerate(df_list):
            df['index'] = i
        df = pd.concat(df_list, ignore_index=True)
        return df

    def load_data(self, bnd_file, num):
        try:
            df = pd.read_csv(bnd_file, header=None, delimiter="\s+", dtype=str, skiprows=1).dropna()
            df[0] = df[0].astype(int)
            df[1] = df[1].astype(float).round(decimals=5)
            df[2] = df[2].astype(float).round(decimals=5)
            df = df[df[0] == num]
            return df if not df.empty else None
        except:
            return

    def get_df_extremum(self, min_max):
        if min_max == "min":
            df = self.df_conduction
            index = df.iloc[:, 2].idxmin()
        elif min_max == "max":
            df = self.df_valence
            index = df.iloc[:, 2].idxmax()
        E_ext = round(df.iloc[index, 2], self.nr+1)
        df0 = df[abs(df.iloc[:, 2] - E_ext) <= self.tol]
        return df0, E_ext

    @staticmethod
    def get_k_dict(df):
        k_points_dict = {}
        for index, row in df.iterrows():
            x = row[1]
            k1, k2, k3 = row[4], row[5], row[6]
            k_points_dict[x] = (k1, k2, k3)
        return k_points_dict
