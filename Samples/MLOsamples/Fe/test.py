"""Fe MLO test: job_mloW + diagonal V, W-V check (no reference files)."""
import glob, os, subprocess, sys
from comp import rmfiles

# Expected on-site diagonal <i i | V/W-V | i i> at R=(0,0,0), omega=0  [eV]
# Baseline: 2026-05-08 job_mloW fe -np 8 (Fe bcc, 1 Fe atom, s+3p+5d = 9 orbitals)
EXPECTED = {
    'UP': {1:(10.6977,-9.7736), 2:(10.7015,-9.3768), 3:(10.7015,-9.3774),
           4:(10.7015,-9.3768), 5:(23.9480,-22.3433), 6:(23.9480,-22.3431),
           7:(23.9239,-22.1601), 8:(23.9480,-22.3429), 9:(23.9239,-22.1602)},
    'DN': {1:(10.6868,-9.7572), 2:(10.7161,-9.3825), 3:(10.7161,-9.3831),
           4:(10.7161,-9.3825), 5:(23.1946,-21.6553), 6:(23.1946,-21.6552),
           7:(23.0995,-21.4155), 8:(23.1946,-21.6550), 9:(23.0995,-21.4156)},
}
TOL = 0.05  # eV; tighter than physical changes, looser than numerical noise


def _parse_diag(workdir, spin):
    cv = os.path.join(workdir, f'Coulomb_v.{spin}')
    sw = os.path.join(workdir, f'Screening_W-v.{spin}')
    v, wmv = {}, {}
    with open(cv) as fh:
        for line in fh:
            f = line.split()
            if (len(f) >= 12 and f[0] == 'Wannier'
                    and float(f[3]) == 0 and float(f[4]) == 0 and float(f[5]) == 0
                    and f[7] == f[8] == f[9] == f[10]):
                v[int(f[7])] = float(f[11])
    with open(sw) as fh:
        for line in fh:
            f = line.split()
            if (len(f) >= 14 and f[0] == 'Wannier'
                    and float(f[3]) == 0 and float(f[4]) == 0 and float(f[5]) == 0
                    and f[7] == f[8] == f[9] == f[10]
                    and float(f[11]) == 0):
                wmv[int(f[7])] = float(f[13])
    return {i: (v[i], wmv[i]) for i in v if i in wmv}


def test(args, bindir, testdir, workdir):
    if args.checkonly:
        from comp import runprogs
        runprogs(['rm -rf summary.txt'], quiet=True)
        return ''

    # cleanup prior outputs (keep bnd00*.spin? as job_mloW prerequisite)
    for pat in ('HamiltonianPMT.*', '__HamiltonianPMT*',
                'PROCAR.UP.*', 'PROCAR.DN.*',
                'Coulomb_v.*', 'Screening_W-v.*'):
        for p in glob.glob(os.path.join(workdir, pat)):
            try: os.remove(p)
            except OSError: pass
    rmfiles(workdir, ['HamRsMLO', 'band_MLO_spin1.dat', 'band_MLO_spin2.dat',
                      'PROCAR.UP', 'PROCAR.DN', 'lwriteham', 'lmlo',
                      'bandplot_MLO.isp1.glt', 'bandplot_MLO.isp2.glt'])

    # Run job_mloW with live-streamed stdout/stderr to terminal
    cmd = [f'{bindir}/job_mloW', 'fe', '-np', str(args.np)]
    print('### exec:', ' '.join(cmd), flush=True)
    proc = subprocess.run(cmd, cwd=workdir)
    if proc.returncode != 0:
        msg = f'FAILED! job_mloW returncode={proc.returncode}'
        print(msg)
        with open(os.path.join(workdir, 'summary.txt'), 'a') as fh: print(msg, file=fh)
        return 'err! '

    # Compare diagonal V, W-V vs EXPECTED (no file kept as reference)
    tall = ''
    for sp in ('UP', 'DN'):
        try:
            actual = _parse_diag(workdir, sp)
        except FileNotFoundError as e:
            print(f'FAILED! parse spin {sp}: {e}')
            tall += 'err! '
            continue
        print(f'\n=== Fe MLO diagonal V, W-V spin {sp} (tol={TOL} eV) ===')
        print(f'  {"i":>3}  {"V_exp":>9} {"V_act":>9} {"dV":>7}  {"WmV_exp":>9} {"WmV_act":>9} {"dWmV":>7}')
        spin_ok = True
        for i, (vexp, wexp) in sorted(EXPECTED[sp].items()):
            va_wa = actual.get(i)
            if va_wa is None:
                print(f'  {i:3d}  MISSING in run')
                spin_ok = False
                continue
            vact, wact = va_wa
            dv, dw = abs(vact - vexp), abs(wact - wexp)
            ok = dv < TOL and dw < TOL
            mark = '   ' if ok else ' X'
            print(f'  {i:3d}  {vexp:9.4f} {vact:9.4f} {dv:7.4f}  {wexp:9.4f} {wact:9.4f} {dw:7.4f}{mark}')
            spin_ok = spin_ok and ok
        out = 'ok! ' if spin_ok else 'err! '
        tag = 'PASSED' if spin_ok else 'FAILED'
        msg = f'{tag}! Fe MLO diagonal V/W-V spin {sp}'
        print(msg)
        with open(os.path.join(workdir, 'summary.txt'), 'a') as fh: print(msg, file=fh)
        tall += out
    return tall
