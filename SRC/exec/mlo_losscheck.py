#!/usr/bin/env python3
"""MLO バンドの損失関数 (縮退に強い版)。

窓        : VBM-2 eV .. CBM+2 eV (金属は E_F±2 eV)
対応づけ  : 順序を保つ最適割当 (DP)。縮退・交差・バンド数の差に耐える。
            MLO 模型は DFT より本数が少ないので、順位で揃えると偽の誤差が出る。
誤差      : 一方向では決まらないので、ベクトルで返す。用途によってどれを
            重んじるかが変わるため、スカラーに潰さない。規格化も 2 乗もせず、
            **生の誤差のまま**出す (どれくらいずれているかが直接読めるように)。

  dE_win   窓 VBM-2..CBM+2 eV (金属 EF±2) の ΔE rms   [meV]
  dE_occ   EF-12..EF-2 eV の ΔE rms [meV]。窓だけ見ると見落とす
  dv/v     分散 (経路方向の勾配) の相対誤差 [無次元]
  dE_gap   バンドギャップ誤差 [meV]   (絶縁体のみ、符号付き)
  m*比     m*_MLO / m*_DFT [無次元]   (絶縁体のみ、1.00 が一致)

  m*比 は「比」であって誤差ではない。1.18 なら MLO の質量が 18% 重い。

注意: m* は補間の質に敏感で、粗い k メッシュでは雑音が乗る。Si では 6^3 の m* が
      1.64 なのに 8^3 で 1.05、10^3 で 0.93 になる (w=0.13)。しかも 6^3 は系統的に
      大きい w を良く見せる。L_edge を使う最適化は 8^3 以上、かつ 2 つのメッシュで
      一致する範囲でのみ意味を持つ。

使い方: mlo_losscheck.py <workdir> [<workdir> ...]
        work ディレクトリの band_MLO_spin1.dat (MLO) と bnd0*.spin1 (DFT) を比べる。
"""
import glob, os, re, sys
import numpy as np

RY = 13.605
BASE = os.path.expanduser('~/ecalj/Samples/MLOsamples')

PAD = 2.0                            # eV: 窓の VBM 下 / CBM 上への張り出し
DXS = (0.10, 0.15, 0.20)             # m* フィット窓 (経路座標)。中央値を採用
OCC_LO, OCC_HI = -12.0, -2.0         # eV: 深い占有側の範囲


def read_ef(g):
    for l in open(g):
        m = re.search(r'ef\s*=\s*([-+0-9.eEdD]+)', l)
        if m:
            return float(m.group(1).replace('D', 'E'))


def load_mlo(fn, ef):
    """band_MLO_spin*.dat を読む。列は x, E(Ry), スピン, バンド番号。
    経路セグメントの継ぎ目の k は 2 回書かれる (バンド番号 1..n が 2 周する) ので、
    (k, バンド番号) で重複を落とす。落とさないと Si は 92 k 点中 5 点だけ
    18 本でなく 36 本になり、準位の対応づけが壊れる。"""
    d, seen = {}, set()
    for l in open(fn):
        if l.lstrip().startswith('#'):
            continue
        t = l.split()
        if len(t) < 2:
            continue
        try:
            x, e = round(float(t[0]), 5), (float(t[1]) - ef) * RY
        except ValueError:
            continue
        if len(t) >= 4:
            key = (x, t[3])
            if key in seen:
                continue
            seen.add(key)
        d.setdefault(x, []).append(e)
    for k in d:
        d[k].sort()
    return d


def load_bnds(fs):
    """bnd00*.spin1 を読む。各ファイルは経路の 1 セグメントで、隣り合う
    セグメントは端点の k を共有する。素朴に追記すると共有点だけ準位が二重に
    並ぶので (Si: 総本数 75 のところ 5 点だけ 150)、既出の k は読み飛ばす。"""
    d = {}
    for f in fs:
        seg = {}
        for l in open(f):
            if l.lstrip().startswith('#'):
                continue
            t = l.split()
            if len(t) < 3:
                continue
            try:
                seg.setdefault(round(float(t[1]), 5), []).append(float(t[2]))
            except ValueError:
                pass
        for k, v in seg.items():
            if k not in d:
                d[k] = v
    for k in d:
        d[k].sort()
    return d


def is_insulator(B):
    """真の絶縁体なら占有本数が k によらず一定。大域の VBM/CBM の差だけで判定すると
    E_F 交差が経路の標本から抜けた金属 (Cu) を絶縁体と誤る。"""
    occ = [e for v in B.values() for e in v if e < 0.0]
    uno = [e for v in B.values() for e in v if e >= 0.0]
    if not occ or not uno:
        return False
    if len({sum(1 for e in v if e < 0.0) for v in B.values()}) != 1:
        return False
    return min(uno) - max(occ) >= 0.15


def window(B):
    if not is_insulator(B):
        return -PAD, PAD
    occ = [e for v in B.values() for e in v if e < 0.0]
    uno = [e for v in B.values() for e in v if e >= 0.0]
    return max(occ) - PAD, min(uno) + PAD


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
    r = dict(rmsE=rmsE, rmsO=rms_range(M, B, OCC_LO, OCC_HI), rmsV=rmsV,
             npt=len(dE), win=(lo, hi), gap=None, mstar=None, rough=None)
    e = edge_terms(M, B)
    if e is not None:
        r['gap'], r['mstar'], r['rough'] = e
    return r


def rms_range(M, B, lo, hi):
    """エネルギー範囲 [lo, hi] に限った ΔE rms。窓の外 (深い占有側) を見るのに使う。"""
    d = []
    for x in sorted(set(M) & set(B)):
        b = [e for e in B[x] if lo <= e <= hi]
        m = [e for e in M[x] if lo <= e <= hi]
        for a, c in align(m, b):
            d.append(m[a] - b[c])
    return float(np.sqrt(np.mean(np.array(d) ** 2))) if d else np.nan


def cbm_curvature(P, thr):
    """CBM の (位置, エネルギー, 曲率) 。頂点固定フィット E-e0 = c2 (x-x0)^2 を
    3 つのフィット窓 DXS で行う。頂点固定なので Gamma 中心 CBM でも 3 点で立つ。
    m* ∝ 1/c2 なので、比 c2_DFT/c2_MLO が m*_MLO/m*_DFT になる。"""
    cand = [(x, e) for x, e in P if e >= thr]
    if not cand:
        return None
    x0, e0 = min(cand, key=lambda p: p[1])
    cs = []
    for dx in DXS:
        sel = {}
        for x, e in P:
            if abs(x - x0) <= dx and e >= thr:
                sel[round(x, 5)] = min(sel.get(round(x, 5), 1e9), e)
        xs = np.array(sorted(sel)); es = np.array([sel[x] for x in xs])
        u = (xs - x0) ** 2
        cs.append(float(np.dot(u, es - e0) / np.dot(u, u))
                  if len(xs) >= 3 and np.dot(u, u) > 0 else np.nan)
    return x0, e0, np.array(cs)


def edge_terms(M, B):
    """絶縁体なら (ギャップ誤差 eV, m*比, m* のフィット窓依存性) を返す。金属は None。"""
    flat = lambda d: [(x, e) for x, v in d.items() for e in v]
    fb, fm = flat(B), flat(M)
    occ = [e for _, e in fb if e < 0.0]; uno = [e for _, e in fb if e >= 0.0]
    if not occ or not uno:
        return None
    vbm, cbm = max(occ), min(uno)
    # 金属判定は「大域の VBM/CBM が離れているか」では不十分。経路の刻みが粗いと
    # E_F 交差が標本から抜け、金属が絶縁体に見える (Cu は大域ギャップ 0.20 eV なのに
    # 92 k 点中 1 点しか |E-E_F|<0.1 eV を持たない)。真の絶縁体なら占有本数が k に
    # よらず一定なので、それで判定する (Cu 5/6, Fe 5-9, RuO2 18-21 に対し
    # Si/GaAs/Al2O3/NiO は一定)。
    if not is_insulator(B):
        return None
    thr = 0.5 * (vbm + cbm)
    mi, bi = cbm_curvature(fm, thr), cbm_curvature(fb, thr)
    if mi is None or bi is None:
        return None
    vm = max((e for _, e in fm if e < thr), default=None)
    if vm is None:
        return None
    ratio = bi[2] / mi[2]           # m*_MLO / m*_DFT
    return ((mi[1] - vm) - (cbm - vbm), float(ratio[1]),
            float(np.nanstd(ratio) / np.nanmean(ratio)))


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
    print(f"{'系':20s} {'窓 (eV)':>15s} | {'dE_win':>8s} {'dE_occ':>8s} "
          f"{'dv/v':>7s} {'dE_gap':>8s} {'m*比':>7s}")
    print('-' * 80)
    for w in args:
        path = w if os.path.isdir(w) else os.path.join(BASE, w)
        name = os.path.basename(path.rstrip('/')).replace('_work', '')
        try:
            v = eval_workdir(path)
        except Exception as e:
            print(f'{name:20s} (評価不能: {e})'); continue
        if v is None:
            print(f'{name:20s} (データ不足)'); continue
        occ = f"{v['rmsO']*1000:8.1f}" if v['rmsO'] == v['rmsO'] else f"{'—':>8s}"
        gap = f"{v['gap']*1000:+8.0f}" if v['gap'] is not None else f"{'—':>8s}"
        ms = f"{v['mstar']:7.2f}" if v['mstar'] is not None else f"{'—':>7s}"
        print(f"{name:20s} [{v['win'][0]:6.2f},{v['win'][1]:6.2f}] | "
              f"{v['rmsE']*1000:8.1f} {occ} {v['rmsV']:7.3f} {gap} {ms}")
    print()
    print('dE_win 窓内 / dE_occ 占有側 EF-12..-2 eV   [meV]')
    print('dv/v   分散の相対誤差 [無次元]')
    print('dE_gap ギャップ誤差 [meV, 符号付き] / m*比 = m*_MLO/m*_DFT [1.00 が一致]')
    print('金属には CBM が無いので後ろ 2 つは —。合計は取らない。')
