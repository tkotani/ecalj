#!/usr/bin/env python3
"""損失をフェルミ面近傍の分散で測る。
金属: E_F を横切るバンドの、E_F 近傍 (|E-E_F|<Ew) での
      (a) エネルギー誤差 (b) 傾き dE/dk の誤差 (=フェルミ速度) (c) E_F 交差位置 k_F の誤差
絶縁体: VBM/CBM 近傍で同じ量(band edge の分散)
"""
import glob, os, re, sys, bisect
import numpy as np

RY = 13.605
BASE = os.path.expanduser('~/ecalj/Samples/MLOsamples')

def read_ef(g):
    for l in open(g):
        m = re.search(r'ef\s*=\s*([-+0-9.eEdD]+)', l)
        if m: return float(m.group(1).replace('D', 'E'))

def load_mlo(f, ef):
    P = []
    for l in open(f):
        t = l.split()
        if len(t) < 2 or l.lstrip().startswith('#'): continue
        try: P.append((float(t[0]), (float(t[1]) - ef) * RY))
        except ValueError: pass
    return P

def load_bnd(fs):
    P = []
    for f in fs:
        for l in open(f):
            t = l.split()
            if len(t) < 3 or l.lstrip().startswith('#'): continue
            try: P.append((float(t[1]), float(t[2])))
            except ValueError: pass
    return P

def by_x(P):
    d = {}
    for x, e in P: d.setdefault(round(x, 5), []).append(e)
    for k in d: d[k].sort()
    return d

def near_ef(P, Ew):
    """|E| < Ew の点だけ"""
    return [(x, e) for x, e in P if abs(e) < Ew]

def slope_err(M, B, Ew=1.0):
    """E_F 近傍のバンドの傾き dE/dx を数値微分で比べる。
    各 x で |E|<Ew の準位を集め、隣接 x との差分で傾きを取り、
    MLO と DFT で最も近い準位同士を対応させて傾き差の rms を返す。"""
    mb, bb = by_x(M), by_x(B)
    xs = sorted(set(mb) & set(bb))
    if len(xs) < 3: return None, None, 0
    de, ds, n, sb_all = [], [], 0, []
    for i in range(1, len(xs) - 1):
        x0, xm, xp = xs[i], xs[i-1], xs[i+1]
        dxm, dxp = x0 - xm, xp - x0
        if dxm <= 0 or dxp <= 0 or dxm > 0.08 or dxp > 0.08: continue
        for eb in bb[x0]:
            if abs(eb) >= Ew: continue
            # DFT 側の傾き
            def nearest(lst, v): return min(lst, key=lambda z: abs(z - v)) if lst else None
            ebm, ebp = nearest(bb[xm], eb), nearest(bb[xp], eb)
            if ebm is None or ebp is None: continue
            sb = (ebp - ebm) / (dxm + dxp)
            # MLO 側で対応する準位
            em = nearest(mb[x0], eb)
            if em is None or abs(em - eb) > 0.5: continue
            emm, emp = nearest(mb[xm], em), nearest(mb[xp], em)
            if emm is None or emp is None: continue
            sm = (emp - emm) / (dxm + dxp)
            de.append(em - eb); ds.append(sm - sb); n += 1; sb_all.append(abs(sb))
    if not de: return None, None, None, 0
    de, ds = np.array(de), np.array(ds)
    sref = float(np.sqrt((np.array(sb_all)**2).mean())) if sb_all else float('nan')
    rms_s = float(np.sqrt((ds**2).mean()))
    return float(np.sqrt((de**2).mean())), rms_s, rms_s/sref if sref else float('nan'), n

def kf_err(M, B, Ew=1.0):
    """E_F 交差位置 (k_F) の誤差: |E|<Ew の準位の符号変化点を x で線形補間"""
    def crossings(P):
        d = by_x(P); xs = sorted(d); out = []
        for i in range(len(xs) - 1):
            if xs[i+1] - xs[i] > 0.08: continue
            for e0 in d[xs[i]]:
                if abs(e0) >= Ew: continue
                e1 = min(d[xs[i+1]], key=lambda z: abs(z - e0))
                if e0 * e1 < 0:
                    out.append(xs[i] + (xs[i+1]-xs[i]) * abs(e0)/(abs(e0)+abs(e1)))
        return sorted(out)
    cm, cb = crossings(M), crossings(B)
    if not cb: return None, 0
    errs = []
    for x in cb:
        if not cm: break
        errs.append(min(abs(x - y) for y in cm))
    if not errs: return None, 0
    return float(np.sqrt((np.array(errs)**2).mean())), len(cb)

def report(cases, Ew=1.0):
    print(f"{'case':16s} {'ΔE_rms(eV)':>11s} {'Δv/v(相対)':>11s} {'Δk_F':>8s} {'n_pts':>6s}")
    for w in cases:
        p = os.path.join(BASE, w)
        try:
            ef = read_ef(os.path.join(p, 'bandplot_MLO.isp1.glt'))
            M = load_mlo(os.path.join(p, 'band_MLO_spin1.dat'), ef)
            B = load_bnd(sorted(glob.glob(os.path.join(p, 'bnd0*.spin1'))))
        except Exception as e:
            print(f'{w:16s} (読めず: {e})'); continue
        de, ds, dsr, n = slope_err(M, B, Ew)
        kf, nc = kf_err(M, B, Ew)
        name = w.replace('_work', '')
        print(f"{name:16s} {de if de is not None else float('nan'):11.4f} "
              f"{dsr if dsr is not None else float('nan'):11.3f} "
              f"{kf if kf is not None else float('nan'):8.4f} {n:6d}")

if __name__ == '__main__':
    report(sys.argv[1:] or ['Fe__m3_work','FeCo__m3_work','Cu__spd_work','RuO2__m3_work','FeMgO__m3_work'])
