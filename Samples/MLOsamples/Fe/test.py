"""Fe MLO test: job_mloW + diagonal V, W-V check (no reference files)."""
import glob, os, subprocess, sys
from comp import rmfiles

# Expected on-site diagonal <i i | V/W-V | i i> at R=(0,0,0), omega=0  [eV]
# Baseline: 2026-09-16 job_mloW fe -np 4 (Fe bcc, 1 Fe atom, s+3p+5d = 9 orbitals)
#           with mlo_method = 4 (Delta = w = 2.0 eV).
#
# The previous baseline (2026-05-08) used mlo_method = 0. Switching to 4 leaves
# s and p almost unchanged but moves the five d orbitals by 0.9 eV (UP) and
# 1.7-2.1 eV (DN): V 23.95 -> 22.98 and 23.19 -> 21.45. That is a real change,
# not noise -- with the default w = 2.0 eV the floor at ecbot+Delta sits above
# E_F for a metal, so more PMT weight enters and the d MLOs come out more
# extended, lowering the on-site U. Fe is one of the systems whose band fit
# improves with a much wider w (measured optimum ~11 eV); see
# https://ecalj.github.io/ecaljdoc/manual/mlo section 2.
#
# Measured with mlo_w = 11 eV instead of the 2.0 default:
#     band error   dE_win 19.8 -> 9.7 meV,  dv/v 0.057 -> 0.024
#                  (the gain is ABOVE E_F, where the floor acts; below E_F the
#                   two are mixed -- see the Fe_w_error figure in the doc)
#     d-orbital W  UP 1.5815 -> 1.6308 eV,  DN 1.4248 -> 1.5546 eV
#
# Which W is right is NOT settled here. All that is established is that W
# depends on mlo_w at the 10% level. This sample keeps the default w = 2.0 so
# that all 18 samples share one setting; if you use W itself, check its w
# dependence for your own system rather than taking either value on faith.
EXPECTED = {
    'UP': {1:(10.6227,-9.7058), 2:(10.7465,-9.3756), 3:(10.7465,-9.3762),
           4:(10.7465,-9.3756), 5:(22.9839,-21.4666), 6:(22.9839,-21.4664),
           7:(23.0628,-21.3852), 8:(22.9839,-21.4663), 9:(23.0628,-21.3853)},
    'DN': {1:(10.5700,-9.6560), 2:(10.7341,-9.3615), 3:(10.7341,-9.3621),
           4:(10.7341,-9.3615), 5:(21.4495,-20.0617), 6:(21.4495,-20.0616),
           7:(20.9761,-19.4960), 8:(21.4495,-20.0614), 9:(20.9762,-19.4962)},
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

    cmd = [f'{bindir}/job_mlo', 'fe', '-np', str(args.np), '--mlo_diagnorm']
    print('### exec:', ' '.join(cmd), flush=True)
    proc = subprocess.run(cmd, cwd=workdir)
    # Run job_mloW with live-streamed stdout/stderr to terminal
    cmd = [f'{bindir}/job_mloW', 'fe', '-np', str(args.np), '--mlo_diagnorm']
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
