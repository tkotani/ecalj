#!/usr/bin/env python3
"""Cross-check m_GWinput Fortran reader against Python tomllib.

For each Samples/*/GWinput:
  1. Convert to GWinput.toml via gwinput2toml.py
  2. Load via Fortran gwinput_dump  -> capture key=value lines
  3. Load via tomllib (Python)      -> compute expected key=value lines
  4. Diff the two; flag any mismatch.

This catches type mismatches, missing keys, default-value drift.
"""
import sys, os, subprocess, tempfile, shutil, tomllib
from pathlib import Path

ECALJ = Path(os.environ.get('ECALJ_DIR', os.path.expanduser('~/ecaljdeveloper')))
CONV = ECALJ / 'SRC/exec/gwinput2toml.py'
DUMP = ECALJ / 'SRC/exec/build/gwinput_dump'
SAMPLES = ECALJ / 'Samples'

# Keys we cross-check (subset of [gw] section that has an exact 1:1 mapping)
CHECK_INT  = ['BZmesh', 'iSigMode', 'niw', 'nband_chi0', 'nband_sigm']
CHECK_REAL = ['QpGcut_psi', 'QpGcut_cou', 'alpha_OffG', 'emax_sigm', 'esmr', 'delta',
              'deltaw', 'HistBin_ratio', 'HistBin_dw', 'ecut_p', 'ecuts_p', 'gauss_img',
              'EMINforGW', 'EMAXforGW']
CHECK_BOOL = ['unit_2pioa', 'GaussSmear', 'KeepPositiveCou']
CHECK_IVEC3 = ['n1n2n3', 'n1n2n3eps', 'multitet']

def parse_dump(text):
    """Parse 'key = value' lines from gwinput_dump output."""
    vals = {}
    for line in text.splitlines():
        if '=' not in line: continue
        k, _, v = line.partition('=')
        k = k.strip()
        v = v.strip()
        vals[k] = v
    return vals

def fortran_int(v): return int(v)
def fortran_real(v):
    if '*' in v: return float('inf')   # f10.5 overflow shows as ******
    return float(v)
def fortran_bool(v): return v.upper().startswith('T')
def fortran_ivec3(v):
    return [int(x) for x in v.split()][:3]

def fmt_real(x, tol=1e-5):
    return x

def close_real(a, b, rel=1e-4, abs_=1e-9):
    if a == b: return True
    if a == float('inf') and b >= 1e9: return True   # overflow case
    if abs(a) < abs_ and abs(b) < abs_: return True
    if a == 0 or b == 0: return abs(a-b) < abs_
    return abs(a-b) / max(abs(a), abs(b)) < rel

def check_one(gwinput_path, work):
    name = str(gwinput_path).replace(str(SAMPLES) + '/', '').replace('/GWinput', '').replace('/', '_')
    shutil.copy(gwinput_path, work / 'GWinput')
    # 1. Convert
    r = subprocess.run(['python3', str(CONV), 'GWinput'], cwd=work, capture_output=True, text=True)
    if r.returncode != 0:
        return name, ['convert failed: ' + r.stderr.strip()]
    # 2. Fortran dump
    r = subprocess.run([str(DUMP), 'GWinput.toml'], cwd=work, capture_output=True, text=True)
    if r.returncode != 0:
        return name, ['fortran dump failed: ' + r.stderr.strip()[:200]]
    fdump = parse_dump(r.stdout)
    # 3. Python tomllib
    with open(work / 'GWinput.toml', 'rb') as f:
        toml = tomllib.load(f)
    gw = toml.get('gw', {})

    diffs = []
    for k in CHECK_INT:
        f_v = fortran_int(fdump.get(k, '0'))
        if k not in gw: continue
        raw = gw[k]
        if isinstance(raw, list):
            # legacy reads first int from list (e.g. nband_sigm)
            p_v = int(raw[0]) if raw else None
        else:
            p_v = int(raw)
        if p_v is not None and f_v != p_v:
            diffs.append(f'{k}: fortran={f_v} python={p_v}')
    for k in CHECK_REAL:
        f_v = fortran_real(fdump.get(k, '0'))
        if k not in gw: continue
        raw = gw[k]
        if isinstance(raw, list):
            p_v = float(raw[0]) if raw else None
        else:
            p_v = float(raw)
        if p_v is not None and not close_real(f_v, p_v):
            diffs.append(f'{k}: fortran={f_v} python={p_v}')
    for k in CHECK_BOOL:
        f_v = fortran_bool(fdump.get(k, 'F'))
        p_v = bool(gw[k]) if k in gw else None
        if p_v is not None and f_v != p_v:
            diffs.append(f'{k}: fortran={f_v} python={p_v}')
    for k in CHECK_IVEC3:
        f_v = fortran_ivec3(fdump.get(k, '0 0 0'))
        if k not in gw: continue
        raw = list(gw[k])
        # Fortran reads at most 3; truncate Python side for fair comparison
        p_v = raw[:3]
        if f_v != p_v:
            diffs.append(f'{k}: fortran={f_v} python={p_v} (raw={raw})')
    return name, diffs

def main():
    inputs = sorted(p for p in SAMPLES.rglob('GWinput') if '_work' not in str(p))
    work_root = Path(tempfile.mkdtemp(prefix='gwinput_xcheck_'))
    try:
        n_pass = n_fail = 0
        fails = []
        for i, gw in enumerate(inputs):
            work = work_root / f'case_{i}'
            work.mkdir()
            name, diffs = check_one(gw, work)
            if diffs:
                n_fail += 1
                fails.append((name, diffs))
                print(f'FAIL: {name}')
                for d in diffs: print(f'  {d}')
            else:
                n_pass += 1
        print()
        print(f'===== Cross-check (Fortran toml-f vs Python tomllib) =====')
        print(f'PASS: {n_pass}')
        print(f'FAIL: {n_fail}')
        if fails:
            sys.exit(1)
        print('ALL OK')
    finally:
        shutil.rmtree(work_root, ignore_errors=True)

if __name__ == '__main__':
    main()
