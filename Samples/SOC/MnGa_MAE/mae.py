#!/usr/bin/env python3
'''Magnetic anisotropy energy by the force theorem, from the logs of lmf.
usage: python3 mae.py [log_scf log_so0 log_so001 log_so110]   (default: llmf_scf llmf_so0 llmf_so001 llmf_so110)
Reads the last line with sev(eV)= (band energy = sum of the occupied eigenvalues) of each log and writes
  energies.txt : mmom, ehf, ehk of the SCF and sev of the three one-shot runs (eV)
  mae.txt      : MAE(meV) = 1000*[ sev(so=1, axis 110) - sev(so=1, axis 001) ]   (positive: 001 is the easy axis)
'''
import re, sys

def lastvalues(logfile):
    'mmom, ehf, ehk, sev of the last line of logfile that has sev(eV)='
    found = None
    for line in open(logfile):
        m = re.search(r'mmom=\s*(\S+)\s+ehf\(eV\)=\s*(\S+)\s+ehk\(eV\)=\s*(\S+)\s+sev\(eV\)=\s*(\S+)', line)
        if m: found = [float(x) for x in m.groups()]
    if found is None: sys.exit(f'mae.py: no line with sev(eV)= in {logfile}')
    return found

logs = sys.argv[1:5] if len(sys.argv) >= 5 else ['llmf_scf', 'llmf_so0', 'llmf_so001', 'llmf_so110']
scf, so0, so001, so110 = [lastvalues(f) for f in logs]
mae = 1000*(so110[3] - so001[3])
with open('energies.txt', 'w') as f:
    f.write('# energies (eV) and the magnetic moment (Bohr magneton) from ' + ' '.join(logs) + '\n')
    f.write('# run               mmom         ehf(eV)          ehk(eV)         sev(eV)\n')
    for name, v in zip(['SCF_so0', 'oneshot_so0', 'oneshot_so1_001', 'oneshot_so1_110'], [scf, so0, so001, so110]):
        f.write(f'{name:16s} {v[0]:8.4f} {v[1]:16.6f} {v[2]:16.6f} {v[3]:14.6f}\n')
with open('mae.txt', 'w') as f:
    f.write('# MAE(meV) = 1000*[sev(so=1,110) - sev(so=1,001)]; positive: 001 is the easy axis\n')
    f.write('# sev(so=1,001)-sev(so=0) and sev(so=1,110)-sev(so=0) are the band-energy gains by SOC (meV)\n')
    f.write(f'MAE_meV          {mae:10.3f}\n')
    f.write(f'gain001_meV      {1000*(so001[3]-so0[3]):10.3f}\n')
    f.write(f'gain110_meV      {1000*(so110[3]-so0[3]):10.3f}\n')
print(open('energies.txt').read(), end='')
print(open('mae.txt').read(), end='')
