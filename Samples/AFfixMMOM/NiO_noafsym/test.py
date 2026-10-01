import re
from comp import runprogs, rmfiles
# AFfixMMOM (2026-10-01): AF NiO with the moment of the Ni pair held at mmtarget.aftest (1.6).
# Compared with the reference out.lmf.nio: the last 'mmaftest:' line (field uhx, moments of sites 1 and 2) and the
# energies of the converged 'c' line. The field is a feedback loop over some 20 iterations, so the tolerances are looser
# than those of test1_check (1e-5 eV), which a compiler can move. m1 + m2 = 0 checks the AF equivalence of the two Ni.
TOL_M, TOL_U, TOL_E = 2e-3, 2e-3, 1e-3   # moments, uhx (Ry), energies (eV)

def last_aftest(f):
    v = None
    for l in open(f, errors='replace'):
        if l.startswith('mmaftest:'):
            v = [float(x) for x in l.split()[2:5]]   # uhx, m1, m2
    return v

def last_c(f):
    v = None
    for l in open(f, errors='replace'):
        m = re.match(r'^c .*ehf\(eV\)=\s*(\S+)\s+ehk\(eV\)=\s*(\S+)', l)
        if m: v = [float(m.group(1)), float(m.group(2))]
    return v

def test(args, bindir, testdir, workdir):
    lmfa = f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    out = 'out.lmf.nio'
    rmfiles(workdir, [out, 'mmagfield.aftest', 'mixmag.aftest'])
    runprogs([lmfa + ' nio > ' + out,
              lmf + ' nio --ctrlg:iter.nit=80 >> ' + out])
    ref, new = testdir + '/' + out, workdir + '/' + out
    a, b = last_aftest(ref), last_aftest(new)
    ca, cb = last_c(ref), last_c(new)
    lines = []
    def chk(label, x, y, tol):
        ok = x is not None and y is not None and abs(x - y) < tol
        lines.append(f"{'PASSED!' if ok else 'FAILED!'} AFfixMMOM {label:22s} ref={x} new={y} tol={tol}")
        return ok
    ok = a is not None and b is not None and ca is not None and cb is not None
    if ok:
        ok &= chk('uhx (Ry)', a[0], b[0], TOL_U)
        ok &= chk('moment site 1', a[1], b[1], TOL_M)
        ok &= chk('moment site 2', a[2], b[2], TOL_M)
        ok &= chk('m1 + m2 (AF)', 0.0, b[1] + b[2], 1e-4)
        ok &= chk('ehf (eV)', ca[0], cb[0], TOL_E)
        ok &= chk('ehk (eV)', ca[1], cb[1], TOL_E)
    else:
        lines.append(f'FAILED! AFfixMMOM no converged mmaftest/c line in {new}')
    with open('summary.txt', 'a') as s:
        for l in lines: print(l); print(l, file=s)
    return 'ok! ' if ok else 'err! '
