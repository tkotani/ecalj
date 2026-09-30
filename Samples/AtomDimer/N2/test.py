from comp import test2_check,runprogs,rmfiles
import os, shutil
DS = [1.05, 1.10, 1.15]      # bond lengths (A) of the dimer; the atom in the same box
def test(args,bindir,testdir,workdir):
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    message1='''
    # Case AtomDimer/N2: the N2 molecule and the N atom in a cubic box of 15 A (PBE, spin polarized, fixed spin moment,
    #  Gamma point only, Fermi-Dirac occupation of 157 K).  The dimer at three bond lengths (subdirectories d_<d>/, the
    #  positions given by --ctrlg:site.<n>.pos), the atom in atom/.  dimer_curve.py: the parabola through the three
    #  energies -> r_e, E_min, D_e = 2 E(N) - E_min.  Checked: dimer.txt (energies 2e-3 eV; r_e and D_e).
    '''
    print(message1)
    if not args.checkonly:
        rmfiles(workdir,['dimer.txt','dimer.npz','dimer.png'])
        cmds = []
        os.makedirs(f'{workdir}/atom', exist_ok=True); shutil.copy(f'{testdir}/ctrlg.n.toml', f'{workdir}/atom/')
        cmds += [f'cd atom && {lmfa} n > llmfa && {lmf} n > llmf']
        prev = None
        for d in DS:
            dd = f'd_{d:.3f}'; os.makedirs(f'{workdir}/{dd}', exist_ok=True); shutil.copy(f'{testdir}/ctrlg.n2.toml', f'{workdir}/{dd}/')
            z = d/2                       # pos is in units of alat = 1 A here
            cp_rst = f'cp ../{prev}/rst.n2 . && ' if prev else ''   # start from the density of the previous bond length
            cmds += [f'cd {dd} && {lmfa} n2 > llmfa && {cp_rst}{lmf} n2 --ctrlg:site.1.pos=[0,0,{-z:.4f}] --ctrlg:site.2.pos=[0,0,{z:.4f}] > llmf']
            prev = dd
        cmds += ['python3 '+testdir+'/dimer_curve.py n2 n '+' '.join(str(d) for d in DS)+' > ldimer_curve']
        runprogs(cmds)
    print('dimer.txt', end=': ')
    return test2_check(testdir+'/dimer.txt', workdir+'/dimer.txt', abs_tol=2e-3)
