#!/usr/bin/env python3
"""CBM 近傍の質だけを見る簡易評価: ギャップ誤差 / m*比 / m*のフィット窓依存性(滑らかさ)"""
import glob, os, re, sys
import numpy as np

RY = 13.605
BASE = os.path.expanduser('~/ecalj/Samples/MLOsamples')
DXS = (0.10, 0.15, 0.20)

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

def cbm_info(P, thr=0.25):
    """CBM の (位置, エネルギー) と各 dx での曲率"""
    cand = [(x, e) for x, e in P if e >= thr]
    if not cand: return None
    x0, e0 = min(cand, key=lambda p: p[1])
    cs = []
    for dx in DXS:
        sel = {}
        for x, e in P:
            if abs(x - x0) <= dx and e >= thr:
                sel[round(x, 5)] = min(sel.get(round(x, 5), 1e9), e)
        xs = np.array(sorted(sel)); es = np.array([sel[x] for x in xs])
        # 頂点固定フィット E-e0 = c2 (x-x0)^2 : 片側しか点が無い Gamma 中心 CBM でも成立
        if len(xs) >= 3:
            u = (xs - x0)**2
            cs.append(float(np.dot(u, es - e0) / np.dot(u, u)) if np.dot(u, u) > 0 else np.nan)
        else:
            cs.append(np.nan)
    return x0, e0, np.array(cs)

def vbm(P, thr=0.25):
    occ = [e for _, e in P if e < thr]
    return max(occ) if occ else None

def evaluate(workdir):
    ef = read_ef(os.path.join(workdir, 'bandplot_MLO.isp1.glt'))
    M = load_mlo(os.path.join(workdir, 'band_MLO_spin1.dat'), ef)
    B = load_bnd(sorted(glob.glob(os.path.join(workdir, 'bnd0*.spin1'))))
    mi, bi = cbm_info(M), cbm_info(B)
    if mi is None or bi is None: return None
    gap_m = mi[1] - vbm(M); gap_b = bi[1] - vbm(B)
    ratio = bi[2] / mi[2]                      # m*_MLO/m*_DFT = c2_DFT/c2_MLO
    r_mid = ratio[DXS.index(0.15)]
    rough = float(np.nanstd(ratio) / np.nanmean(ratio))   # フィット窓依存性 = 非放物線性
    return dict(gap_err=gap_m - gap_b, mstar=r_mid, rough=rough,
                dpos=mi[0] - bi[0], ratios=ratio)

def loss(r):
    """3項のみ: ギャップ / m* / 滑らかさ。単位は「この誤差で 1」"""
    return ((r['gap_err']/0.010)**2 + (np.log(r['mstar'])/0.05)**2
            + (r['rough']/0.05)**2 + (r['dpos']/0.02)**2)

if __name__ == '__main__':
    print(f"{'case':16s} {'gap誤差':>8s} {'m*比':>7s} {'粗さ':>7s} {'Δpos':>7s} {'loss':>9s}   m*(dx=.10/.15/.20)")
    for w in sys.argv[1:]:
        p = os.path.join(BASE, w)
        try: r = evaluate(p)
        except Exception as e:
            print(f'{w:16s} (失敗: {e})'); continue
        if r is None:
            print(f'{w:16s} (CBM 検出できず)'); continue
        print(f"{w.replace('_work',''):16s} {r['gap_err']:+8.3f} {r['mstar']:7.2f} {r['rough']:7.3f} "
              f"{r['dpos']:+7.3f} {loss(r):9.1f}   " + " / ".join(f"{x:.2f}" for x in r['ratios']))
