#!/usr/bin/env python3
"""Grade of an MLO model of gw1500_mlo.sh (2026-10-05, user: a deviation at the top of the window is not bad; steep bands
there deviate easily).

    gw1500_mlo_grade.py <dir> [bindir]

Uses bandcheck.json (window [VBM - 8, CBM + 2] eV, metals E_F) and runs mlo_bandcheck.py --delta 1.5 for the inner window
[VBM - 8, CBM + 1.5] eV (bandcheck_in.json). Writes grade.json and prints the grade:
  PASS  the check of the window passes (largest deviation <= 0.1 eV, no jump > 0.1 eV, no broken band)
  OK    the inner window passes, and in the whole window the largest deviation is <= 0.2 eV and no band is broken (spike < 2 eV)
  FAIL  otherwise
"""
import json, os, subprocess, sys
d = sys.argv[1]; B = sys.argv[2] if len(sys.argv) > 2 else os.environ.get('ECALJ_BIN', os.path.expanduser('~/bin'))
c = next(iter(json.load(open(os.path.join(d, 'bandcheck.json'))).values()))
fi = os.path.join(d, 'bandcheck_in.json')
if not os.path.exists(fi) or os.path.getmtime(fi) < os.path.getmtime(os.path.join(d, 'bandcheck.json')):
    subprocess.run([sys.executable, f'{B}/mlo_bandcheck.py', d, '--delta', '1.5', '--json', fi], capture_output=True, text=True)
ci = next(iter(json.load(open(fi)).values()))
if c['check'] == 'PASS':
    g = 'PASS'
elif ci['check'] == 'PASS' and c['max'] <= 0.2 and c['spike'] < 2.0:
    g = 'OK'
else:
    g = 'FAIL'
json.dump(dict(grade=g, max=c['max'], max_in=ci['max'], jump=c['jump'], spike=c['spike']), open(os.path.join(d, 'grade.json'), 'w'))
print(g)
