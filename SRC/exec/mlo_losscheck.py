#!/usr/bin/env python3
"""MLO バンドの損失関数 (縮退に強い版)。

窓        : VBM-2 eV .. CBM+2 eV (金属は E_F±2 eV)
対応づけ  : 順序を保つ最適割当 (DP)。縮退・交差・バンド数の差に耐える。
            MLO 模型は DFT より本数が少ないので、順位で揃えると偽の誤差が出る。
損失      : L = (ΔE_rms/σ_E)^2 + (Δv/v / σ_v)^2      σ_E = 20 meV, σ_v = 5%
            第1項がマッチング、第2項が微分 (スムーズネス)。

使い方: mlo_losscheck.py <workdir> [<workdir> ...]
        work ディレクトリの band_MLO_spin1.dat (MLO) と bnd0*.spin1 (DFT) を比べる。
"""
import glob, os, re, sys
import numpy as np

RY = 13.605
BASE = os.path.expanduser('~/ecalj/Samples/MLOsamples')

S_E, V_REL, PAD = 0.020, 0.05, 2.0


def read_ef(g):
    for l in open(g):
        m = re.search(r'ef\s*=\s*([-+0-9.eEdD]+)', l)
        if m:
            return float(m.group(1).replace('D', 'E'))


def load_mlo(fn, ef):
    d = {}
    for l in open(fn):
        if l.lstrip().startswith('#'):
            continue
        t = l.split()
        if len(t) < 2:
            continue
        try:
            d.setdefault(round(float(t[0]), 5), []).append((float(t[1]) - ef) * RY)
        except ValueError:
            pass
    for k in d:
        d[k].sort()
    return d


def load_bnds(fs):
    d = {}
    for f in fs:
        for l in open(f):
            if l.lstrip().startswith('#'):
                continue
            t = l.split()
            if len(t) < 3:
                continue
            try:
                d.setdefault(round(float(t[1]), 5), []).append(float(t[2]))
            except ValueError:
                pass
    for k in d:
        d[k].sort()
    return d


def window(B):
    occ = [e for v in B.values() for e in v if e < 0.0]
    uno = [e for v in B.values() for e in v if e >= 0.0]
    if not occ or not uno:
        return -PAD, PAD
    vbm, cbm = max(occ), min(uno)
    return (-PAD, PAD) if (cbm - vbm) < 0.15 else (vbm - PAD, cbm + PAD)


def align(m, b):
    """順序を保つ最適対応 (DP)。m[i] を b[j] に、i,j の順序を保って割り当て
    総 |Δε| を最小化する。返り値は (i, j) のリスト。"""
    n, M = len(m), len(b)
    if n == 0 or M == 0 or n > M:
        if n > M:
            m, b = b, m
            n, M = M, n
            swap = True
        else:
            return []
    else:
        swap = False
    INF = 1e18
    dp = np.full((n + 1, M + 1), INF)
    dp[0, :] = 0.0
    bk = np.zeros((n + 1, M + 1), dtype=np.int8)
    for i in range(1, n + 1):
        for j in range(i, M + 1):
            skip = dp[i, j - 1]
            take = dp[i - 1, j - 1] + abs(m[i - 1] - b[j - 1])
            if take <= skip:
                dp[i, j] = take; bk[i, j] = 1
            else:
                dp[i, j] = skip; bk[i, j] = 0
    if dp[n, M] >= INF:
        return []
    pairs, i, j = [], n, M
    while i > 0 and j > 0:
        if bk[i, j] == 1:
            pairs.append((i - 1, j - 1)); i -= 1; j -= 1
        else:
            j -= 1
    pairs = pairs[::-1]
    return [(b_, a_) for a_, b_ in pairs] if swap else pairs


def evaluate(M, B):
    lo, hi = window(B)
    xs = sorted(set(M) & set(B))
    if len(xs) < 3:
        return None
    win = lambda d, x: [e for e in d[x] if lo <= e <= hi]
    dE, dV, vref = [], [], []
    for i, x in enumerate(xs):
        b, m = win(B, x), win(M, x)
        pr = align(m, b)
        if not pr:
            continue
        dE.extend(m[a] - b[c] for a, c in pr)
        if 0 < i < len(xs) - 1:
            xm, xp = xs[i - 1], xs[i + 1]
            if (x - xm) <= 0.08 and (xp - x) <= 0.08:
                bm, bp = win(B, xm), win(B, xp)
                mm, mp = win(M, xm), win(M, xp)
                dm = dict(align(mm, bm)); dpp = dict(align(mp, bp))
                dx = xp - xm
                for a, c in pr:
                    if a in dm and a in dpp:
                        sb = (bp[dpp[a]] - bm[dm[a]]) / dx
                        sm = (mp[a] - mm[a]) / dx
                        dV.append(sm - sb); vref.append(abs(sb))
    if not dE:
        return None
    dE = np.array(dE)
    dV = np.array(dV) if dV else np.array([0.0])
    vr = float(np.sqrt(np.mean(np.array(vref) ** 2))) if vref else 1.0
    rmsE = float(np.sqrt(np.mean(dE ** 2)))
    rmsV = float(np.sqrt(np.mean(dV ** 2))) / vr if vr > 1e-9 else np.nan
    return dict(rmsE=rmsE, rmsV=rmsV, loss=(rmsE / S_E) ** 2 + (rmsV / V_REL) ** 2,
                npt=len(dE), win=(lo, hi))


def eval_workdir(workdir, mlodat=None, glt=None):
    """work ディレクトリ一つを評価する。
    mlodat/glt を与えれば MLO 側だけ別ファイルから読む(パラメタ走査用)。"""
    glt = glt or os.path.join(workdir, 'bandplot_MLO.isp1.glt')
    mlodat = mlodat or os.path.join(workdir, 'band_MLO_spin1.dat')
    M = load_mlo(mlodat, read_ef(glt))
    B = load_bnds(sorted(glob.glob(os.path.join(workdir, 'bnd0*.spin1'))))
    return evaluate(M, B)


if __name__ == '__main__':
    args = sys.argv[1:]
    if not args:
        sys.exit(f'usage: {os.path.basename(sys.argv[0])} <workdir> [<workdir> ...]')
    print(f"{'系':20s} {'窓 (eV)':>16s} {'ΔE rms':>8s} {'Δv/v':>8s} {'loss':>9s} {'n':>6s}")
    tot = []
    for w in args:
        p = w if os.path.isdir(w) else os.path.join(BASE, w)
        name = os.path.basename(p.rstrip('/')).replace('_work', '')
        try:
            v = eval_workdir(p)
        except Exception as e:
            print(f'{name:20s} (評価不能: {e})'); continue
        if v is None:
            print(f'{name:20s} (データ不足)'); continue
        print(f"{name:20s} [{v['win'][0]:6.2f},{v['win'][1]:6.2f}] {v['rmsE']:8.4f} "
              f"{v['rmsV']:8.3f} {v['loss']:9.1f} {v['npt']:6d}")
        tot.append(v['loss'])
    if len(tot) > 1:
        print(f"{'平均':20s} {'':16s} {'':8s} {'':8s} {np.mean(tot):9.1f}")
