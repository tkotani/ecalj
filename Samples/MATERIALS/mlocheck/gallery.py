#!/usr/bin/env python3
# 2026-10-01 mlocheck page, version 3 (16:1x, user's layout): the baseline is the model with the semicore local orbitals
# added automatically (m_HamPMT: band top above EF - 17 eV, weight over the species). Other models: 2. EH2 s,p on the
# cations, 3. empty spheres (SiO2 only). The old default (until 2026-10-01 14:4x) is shown struck through.
#   Fig. 1, 2 and 3: the best model per material (smallest error); Table 1: baseline and EH2 per material, SiO2 with and
#   without ES in two rows; Fig. 4: the baseline against models 2 and 3.
#   python3 gallery.py <OUT>        (writes OUT/index.html and OUT/img/; the plotted numbers go to page_data/ here)
# Inputs: JSON of SRC/exec/mlo_bandcheck.py --json for each set (made here), the band files of each work directory.
# Version 2 (15:xx) is in the git history of Samples/MATERIALS/mlocheck/gallery.py; version 1 is gallery_v1.py.
import json, os, sys, math, datetime, subprocess, glob, importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': ['DejaVu Sans', 'Noto Sans CJK JP'], 'hatch.linewidth': 0.6})

HERE = os.path.dirname(os.path.abspath(__file__))
W = os.path.dirname(HERE)                                      # ~/work
SETS = {'old': f'{W}/mlocheck_20261001', 'base': f'{W}/mlocheck_auto_20261001', 'eh2': f'{W}/mlocheck_eh2cat'}
ES = {'SiO2c': f'{W}/mlocheck_es/SiO2c_ES'}
OUT = sys.argv[1]
DATA = f'{HERE}/page_data'                                     # the numbers of each figure (npz), next to the scripts
os.makedirs(os.path.join(OUT, 'img'), exist_ok=True); os.makedirs(DATA, exist_ok=True)
mats = json.load(open(f'{HERE}/materials.json'))
EXE = os.path.expanduser('~/ecalj/SRC/exec')
def mod(name):
    s = importlib.util.spec_from_file_location(name, f'{EXE}/{name}.py'); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
BP, BC = mod('mlo_bandplot'), mod('mlo_bandcheck')
rev = subprocess.run(['git', '-C', os.path.expanduser('~/ecalj'), 'log', '-1', '--format=%h'], capture_output=True, text=True).stdout.strip()
NOW = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')

def has_mlo(d): return os.path.exists(f'{d}/band_MLO_spin1.dat') and os.path.getsize(f'{d}/band_MLO_spin1.dat') > 0
def evaluate(dirs, out):
    dirs = [d for d in dirs if has_mlo(d)]
    if dirs: subprocess.run(['python3', f'{EXE}/mlo_bandcheck.py', *dirs, '--json', out], capture_output=True, text=True)
    return json.load(open(out)) if dirs and os.path.exists(out) else {}

R = {k: evaluate([f'{p}/{m}' for m in mats], f'{OUT}/{k}.json') for k, p in SETS.items()}
r_es = evaluate(list(ES.values()), f'{OUT}/es.json')
R['es'] = {m: r_es[os.path.basename(d)] for m, d in ES.items() if os.path.basename(d) in r_es}
json.dump(R, open(f'{OUT}/all.json', 'w'), indent=1)
DIR = {'base': lambda m: f"{SETS['base']}/{m}", 'eh2': lambda m: f"{SETS['eh2']}/{m}", 'es': lambda m: ES[m]}
LAB = {'base': '1. 基準', 'eh2': '2. 陽イオンに EH2 (s,p)', 'es': '3. 空格子球 2 つ (s,p)'}

fail_eh2 = {}
st = f"{SETS['eh2']}/status.log"
if os.path.exists(st):
    for l in open(st):
        t = l.split()
        if len(t) > 3 and t[3] == 'FAIL': fail_eh2[t[2]] = ' '.join(t[4:])

def worst(v):   # the error of a model: max(rms MLO->DFT, rms DFT->MLO, |gap error|)
    if not v: return None
    w = [v['rms_m2d'], v['rms_d2m']]
    if v.get('insulator'): w.append(abs(v['gapM'] - v['gapD']))
    return max(w)
def sev(x):     # good <= 0.02, fair <= 0.05, marginal <= 0.1, poor > 0.1 eV (user 2026-10-01 14:3x)
    if x is None: return 'none'
    return 'good' if x <= 0.02 else ('fair' if x <= 0.05 else ('marginal' if x <= 0.1 else 'poor'))
CL = ['good', 'fair', 'marginal', 'poor']

GROUPS = [
 ('sp', 's,p の半導体', 'C Si Ge Sn 3cSiC 2hSiC 4hSiC AlN AlNzb AlP AlAs AlSb GaN GaNzb GaP GaAs GaSb InN InNzb InP InAs InSb MgS MgSe MgTe'),
 ('d10', 'Zn・Cd・Hg・Pb の化合物', 'ZnO ZnS ZnSe ZnTe wZnS CdO CdS CdSe CdTe wCdS HgO HgS HgSe HgTe PbS PbTe'),
 ('ox', '酸化物・ペロブスカイト', 'MgO SiO2c ZrO2 HfO2 SrTiO3 SrVO3 BaTiO3 LaGaO3 La2CuO4'),
 ('mag', '金属と磁性体', 'Li Cu Fe Ni YMn2 MnO NiO'),
 ('f', '4f（s,p,d,f の模型）', 'Ce EuO EuS EuSe EuTe'),
 ('soc', 'スピン軌道（job_mlo_soc）', 'GaAs_so Bi2Te3'),
 ('sl', '超格子', 'InAsGaSb_n4'),
]
allm = [m for _, _, ms in GROUPS for m in ms.split() if m in R['base']]
missing = [m for _, _, ms in GROUPS for m in ms.split() if m not in R['base']]
W_ = {k: {m: worst(R[k].get(m)) for m in allm} for k in ('old', 'base', 'eh2', 'es')}
BEST = {m: min([(W_['base'][m], 'base')] + [(W_[k][m], k) for k in ('eh2', 'es') if W_[k][m] is not None]) for m in allm}
def counts(k):
    c = {x: 0 for x in CL}
    for m in allm:
        w = BEST[m][0] if k == 'best' else W_[k][m]
        if w is not None: c[sev(w)] += 1
    return c
CNT = {k: counts(k) for k in ('old', 'base', 'eh2', 'best')}

def fmt(x, n=3, sign=False):
    if x is None or (isinstance(x, float) and x != x): return '—'
    if abs(x) < 0.5 * 10 ** -n: x = 0.0                       # no "-0.000"
    return f'{x:+.{n}f}'.replace('-', '−') if sign else f'{x:.{n}f}'
def gap_small(v): return min(v['gap_mesh'], v['gapD']) if v.get('insulator') else None
def dgap(v): return (v['gapM'] - v['gapD']) if v and v.get('insulator') else None

# ---------- band figures ----------
# DFT: grey lines; MLO: red x; hatched: the window of the error [VBM-8, CBM+mlo_delta] (both directions, 17:54);
# orange dots: DFT points with no MLO band within 0.1 eV (a band missing from the
# model); dark rings: MLO points with no DFT band within 0.1 eV (a wrong band). The matching is that of mlo_bandcheck.py.
TOL = 0.1
def spins_of(d): return [s for s in (1, 2) if os.path.exists(f'{d}/band_MLO_spin{s}.dat') and os.path.getsize(f'{d}/band_MLO_spin{s}.dat') > 0]
def far(P, Q, a, b):
    """points (x, E) of P in [a, b] whose nearest band of Q at the same x is more than TOL away"""
    xs = np.array(sorted(Q)); out = []
    for x, e in P:
        if not (a <= e <= b): continue
        i = np.searchsorted(xs, x); c = [xs[k] for k in (i - 1, i) if 0 <= k < len(xs)]
        if not c: continue
        xq = min(c, key=lambda t: abs(t - x))
        if abs(xq - x) <= 0.01 and np.min(np.abs(np.array(Q[xq]) - e)) > TOL: out.append((x, e))
    return np.array(out).reshape(-1, 2)
def panel(ax, d, isp, res, title, store):
    ef = BP.read_ef(d)
    segs = BP.load_dft(d, isp); m = BP.load_mlo(d, isp, ef)
    lo, hi = res['emin'], res['emax']
    ax.axhspan(lo, hi, facecolor='none', edgecolor='#7fa3c7', hatch='////', lw=0, zorder=0)
    for s in segs: ax.plot(s[:, 0], s[:, 1], '-', color='0.40', lw=0.9, zorder=2)
    D, M = {}, {}                                              # both spins together, as mlo_bandcheck.py
    for s in spins_of(d):
        for x, es in BC.load(sorted(glob.glob(f'{d}/bnd*.spin{s}')), lambda e: e, bnd=True).items(): D.setdefault(x, []).extend(es)
        for x, es in BC.load([f'{d}/band_MLO_spin{s}.dat'], lambda e: (e - ef) * BC.RY).items(): M.setdefault(x, []).extend(es)
    dpts = np.concatenate(segs) if segs else np.zeros((0, 2))
    miss = far([(round(x, 5), e) for x, e in dpts], M, lo, hi)
    wrong = far([(round(x, 5), e) for x, e in m], D, lo, hi) if m is not None else np.zeros((0, 2))
    if len(miss): ax.plot(miss[:, 0], miss[:, 1], 'o', color='#e08a1e', ms=4.2, mew=0, alpha=0.9, zorder=3)
    if m is not None: ax.plot(m[:, 0], m[:, 1], 'x', color='#c0392b', ms=3.6, mew=0.9, zorder=4)
    if len(wrong): ax.plot(wrong[:, 0], wrong[:, 1], 'o', mfc='none', mec='#1d2733', ms=6.5, mew=0.9, zorder=5)
    ax.axhline(0.0, color='#2c3e50', lw=0.7, ls='--', zorder=1)
    tics = BP.read_xtics(d)
    if tics:
        ax.set_xticks([x for x, _ in tics]); ax.set_xticklabels([l for _, l in tics])
        for x, _ in tics[1:-1]: ax.axvline(x, color='0.82', lw=0.6, zorder=1)
        ax.set_xlim(tics[0][0], tics[-1][0])
    ax.set_title(title, fontsize=11)
    store.update({f'{title}|dft_x': dpts[:, 0], f'{title}|dft_E': dpts[:, 1], f'{title}|mlo_x': m[:, 0] if m is not None else [],
                  f'{title}|mlo_E': m[:, 1] if m is not None else [], f'{title}|window': np.array([lo, hi]),
                  f'{title}|missing': miss, f'{title}|wrong': wrong, f'{title}|dir': d})
    return lo - 1.2, hi + 2.5

def figure(m, keys, fname):
    """one material; keys = models (columns); spins as columns for one model, as rows for several"""
    pans = [(k, DIR[k](m), R[k][m]) for k in keys if m in R[k]]
    ns = max(len(spins_of(d)) for _, d, _ in pans)
    nrow, ncol = (1, ns) if len(pans) == 1 else (ns, len(pans))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.2 * ncol, 4.5 * nrow), squeeze=False, sharey=True)
    store = {}; ylo, yhi = [], []
    for j, (k, d, res) in enumerate(pans):
        for i, s in enumerate(spins_of(d)):
            ax = axes[0][i] if len(pans) == 1 else axes[i][j]
            a, b = panel(ax, d, s, res, f'{m}  ' + LAB[k] + (f'  spin {s}' if ns > 1 else ''), store); ylo.append(a); yhi.append(b)
    for ax in axes.flat: ax.set_ylim(min(ylo), max(yhi))
    for ax in axes[:, 0]: ax.set_ylabel(r'$E-E_F$  (eV)')
    fig.tight_layout(); fig.savefig(f'{OUT}/img/{fname}', dpi=110); plt.close(fig)
    np.savez_compressed(f'{DATA}/{fname[:-4]}.npz', **{k: np.asarray(v) for k, v in store.items()},
                        _made=f'{NOW} by {os.path.abspath(__file__)} (ecalj {rev})')
    return f'img/{fname}', ns

FIG3 = {m: figure(m, [BEST[m][1]], f'{m}_best.png') for m in allm}
cmp4 = [m for m in allm if m in ES or sev(W_['base'][m]) != 'good' or (W_['eh2'][m] is not None and sev(W_['eh2'][m]) != sev(W_['base'][m]))]
cmp4 += [m for m in ('Cu', 'Ni') if m in allm and m not in cmp4]
cmp4.sort(key=lambda m: (m not in ES, allm.index(m)))
FIG4 = {m: figure(m, ['base', 'eh2'] + (['es'] if m in ES else []), f'{m}_cmp.png') for m in cmp4}

# ---------- chart 1: gap error of the best model; arrows from the baseline where another model is better ----------
X0, X1, Y0, Y1 = 0.0, 6.0, -0.10, 0.45
PW, PH, ML, MR, MT, MB = 720, 360, 56, 18, 16, 44
def sx(x): return ML + (min(max(x, X0), X1) - X0) / (X1 - X0) * (PW - ML - MR)
def sy(y): return MT + (Y1 - min(max(y, Y0), Y1)) / (Y1 - Y0) * (PH - MT - MB)
svg = [f'<svg viewBox="0 0 {PW} {PH}" role="img" aria-labelledby="c1t" class="chart"><title id="c1t">一番良い模型のバンドギャップの誤差</title>']
for yt in (-0.1, 0, 0.1, 0.2, 0.3, 0.4):
    svg.append(f'<line x1="{ML}" x2="{PW-MR}" y1="{sy(yt):.1f}" y2="{sy(yt):.1f}" class="{"zero" if yt == 0 else "grid"}"/>'
               f'<text x="{ML-8}" y="{sy(yt)+4:.1f}" class="tick" text-anchor="end">{yt:+.1f}</text>')
for xt in range(0, 7):
    svg.append(f'<line x1="{sx(xt):.1f}" x2="{sx(xt):.1f}" y1="{MT}" y2="{PH-MB}" class="grid"/>'
               f'<text x="{sx(xt):.1f}" y="{PH-MB+18}" class="tick" text-anchor="middle">{xt}</text>')
svg.append(f'<text x="{(ML+PW-MR)/2}" y="{PH-6}" class="axis" text-anchor="middle">DFT のギャップ（対称線の上、eV）</text>')
svg.append(f'<text x="14" y="{(MT+PH-MB)/2}" class="axis" text-anchor="middle" transform="rotate(-90 14 {(MT+PH-MB)/2})">MLO − DFT（eV）</text>')
for m in allm:
    w0, k = BEST[m]; v = R[k][m]
    if not v.get('insulator'): continue
    x, y = sx(v['gapD']), sy(dgap(v))
    vb = R['base'][m]; d0 = dgap(vb) if vb.get('insulator') else 0.0
    if k != 'base' and vb.get('insulator'):
        y0 = sy(d0)
        if abs(y0 - y) > 6:
            svg.append(f'<line x1="{x:.1f}" y1="{y0:.1f}" x2="{x:.1f}" y2="{y-5 if y0 < y else y+5:.1f}" class="arrow" marker-end="url(#ah)"/>')
        svg.append(f'<circle cx="{x:.1f}" cy="{y0:.1f}" r="4" class="pt was"><title>{m} 1. 基準: {d0:+.3f} eV</title></circle>')
    svg.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" class="pt {sev(w0)}"><title>{m} {LAB[k]}: DFT {v["gapD"]:.3f}, MLO {v["gapM"]:.3f}, 誤差 {dgap(v):+.3f} eV</title></circle>')
    if abs(dgap(v)) > 0.04 or (k != 'base' and abs(d0) > 0.06):
        yl = sy(d0) if k != 'base' else y
        right = x > 0.75 * PW
        svg.append(f'<text x="{x-8 if right else x+8:.1f}" y="{yl+4:.1f}" class="lab" text-anchor="{"end" if right else "start"}">{m}{" (+%.2f eV、枠外)" % d0 if d0 > Y1 else ""}</text>')
svg.append('<defs><marker id="ah" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" class="ahead"/></marker></defs></svg>')
chart1 = '\n'.join(svg)

# ---------- chart 2: error of the best model, sorted, log scale; the bar colour is the model (1, 2, 3; user 2026-10-01 18:0x),
# the verdict is read from the dashed lines; the values of model 2 (EH2 on cations) as a polyline over all materials
# (user 18:2x: a line, not rings) ----------
srt = sorted(allm, key=lambda m: BEST[m][0])
BW, BH, bl, bb, bt = max(720, 11 * len(srt) + 70), 300, 56, 92, 14
L0, L1 = math.log10(0.0005), math.log10(5.0)
def by(val): return bt + (L1 - math.log10(min(max(val, 0.0005), 5.0))) / (L1 - L0) * (BH - bt - bb)
sv = [f'<svg viewBox="0 0 {BW} {BH}" role="img" aria-labelledby="c2t" class="chart"><title id="c2t">物質ごとの一番良い模型の誤差の最大値</title>']
for val, lab in ((0.001, '0.001'), (0.01, '0.01'), (0.02, '0.02'), (0.05, '0.05'), (0.1, '0.1'), (1, '1')):
    cls = 'thr' if val in (0.02, 0.05, 0.1) else 'grid'
    sv.append(f'<line x1="{bl}" x2="{BW-8}" y1="{by(val):.1f}" y2="{by(val):.1f}" class="{cls}"/><text x="{bl-6}" y="{by(val)+4:.1f}" class="tick" text-anchor="end">{lab}</text>')
step = (BW - bl - 8) / len(srt)
for i, m in enumerate(srt):
    w, k = BEST[m]; x = bl + i * step
    sv.append(f'<rect x="{x+1:.1f}" y="{by(w):.1f}" width="{max(step-2,2):.1f}" height="{BH-bb-by(w):.1f}" rx="2" class="bar m-{k}"><title>{m} {LAB[k]}: {w:.3f} eV ({sev(w)})</title></rect>')
    sv.append(f'<text x="{x+step/2:.1f}" y="{BH-bb+8}" class="blab" transform="rotate(60 {x+step/2:.1f} {BH-bb+8})">{m}</text>')
# polylines of model 1 (grey, user 18:3x) and model 2 (blue); a line breaks where the model has no value (EuO stopped
# in model 2). SiO2c is left out of both lines: its best is model 3 and its 1 and 2 lack the conduction band (user 18:3x).
for k, cls in (('base', 'line1'), ('eh2', 'eh2line')):
    segs, cur = [], []
    for i, m in enumerate(srt):
        if m in ES: continue
        v = W_[k][m]
        if v is None:
            if cur: segs.append(cur); cur = []
            continue
        cur.append((bl + (i + 0.5) * step, by(v), m, v))
    if cur: segs.append(cur)
    for s in segs:
        sv.append(f'<polyline class="{cls}" points="' + ' '.join(f'{x:.1f},{y:.1f}' for x, y, _, _ in s) + '"/>')
        for x, y, m, v in s:
            sv.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" class="hit"><title>{m} {LAB[k]}: {v:.3f} eV</title></circle>')
sv.append(f'<text x="14" y="{(bt+BH-bb)/2}" class="axis" text-anchor="middle" transform="rotate(-90 14 {(bt+BH-bb)/2})">eV</text></svg>')
chart2 = '\n'.join(sv)

# ---------- the models and their counts (one table: what is added, where it is written, the verdicts) ----------
CNT['es'] = counts('es')
def nrun(k): return sum(1 for m in allm if (BEST[m][0] if k == 'best' else W_[k][m]) is not None)
def mrow(no, name, what, where, k, cls=''):
    c = CNT[k]; dl = (lambda s: f'<del>{s}</del>') if cls == 'old' else (lambda s: s)
    return (f'<tr class="{cls}"><td class="no">{dl(no)}</td><td>{dl(name)}</td><td class="what">{dl(what)}</td><td class="what">{dl(where)}</td>'
            f'<td class="num">{dl(nrun(k))}</td>' + ''.join(f'<td class="num">{dl(c[x])}</td>' for x in CL) + '</tr>')
counts_tbl = ('<table class="cnt"><thead><tr><th></th><th>模型</th><th>MLO のシード</th><th>書く所</th><th class="num">物質</th>'
              '<th class="num">good<br>≤ 0.02</th><th class="num">fair<br>≤ 0.05</th><th class="num">marginal<br>≤ 0.1</th><th class="num">poor<br>&gt; 0.1 eV</th></tr></thead><tbody>'
  + mrow('', '旧既定（〜2026-10-01 14:4x）', '原子ごと・lm ごとに EH 1 本をシードにする。浅い局所軌道は EH と入れ替え（E_F − 10 eV より上）', '<code>mlo_lm</code>', 'old', 'old')
  + mrow('1', '<b>基準</b>', 'EH に加えて、半内殻の局所軌道（帯の上端が E_F − 17 eV より上）を MLO のシードとして加える', '自動（<code>m_HamPMT</code>）', 'base')
  + mrow('2', '1 + 陽イオンに EH2', '陽イオンの s,p に第 2 の smooth Hankel 関数（EH2）を MLO のシードとして加える', '<code>mlo_lm2</code>', 'eh2')
  + mrow('3', '1 + 空格子球', '空隙に置いた z = 0 の球の s,p を MLO のシードとして加える', '<code>[[site]]</code>・<code>[[spec]]</code>・<code>mlo_lm</code>（DFT から回し直す）', 'es')
  + mrow('', '<b>物質ごとに一番良いもの</b>', '1・2・3 のうち誤差の最大値が一番小さい模型（図 1〜3）', '', 'best', 'best')
  + '</tbody></table>')

# ---------- table 1 ----------
def cells(v, w, best):
    if v is None: return '<td class="num">—</td><td class="num">—</td>'
    return (f'<td class="num">{fmt(dgap(v), sign=True) if v.get("insulator") else "—"}</td>'
            f'<td class="num sev-{sev(w)}{" best" if best else ""}">{fmt(w)}</td>')
def window(v):
    return f'{fmt(v["emin"], 1, True)} 〜 {fmt(v["emax"], 1, True)}'
trs = []
for m in allm:
    v = R['base'][m]; k = BEST[m][1]; e2 = R['eh2'].get(m); g = gap_small(v)
    e2c = cells(e2, W_['eh2'][m], k == 'eh2') if e2 else ('<td class="num" colspan="2">止まった</td>' if m in fail_eh2 else cells(None, None, False))
    trs.append(f'<tr><td class="name"><a href="#fig-{m}">{m}</a></td><td class="num">{fmt(g) if g is not None else "金属"}</td>'
               f'<td class="num win">{window(v)}</td>{cells(v, W_["base"][m], k == "base")}{e2c}</tr>')
    if m in R['es']:
        ve = R['es'][m]
        trs.append(f'<tr class="sub"><td class="name"><a href="#cmp-{m}">{m} + 空格子球</a></td><td class="num">{fmt(gap_small(ve))}</td>'
                   f'<td class="num win">{window(ve)}</td>{cells(ve, W_["es"][m], k == "es")}<td class="num">—</td><td class="num">—</td></tr>')
table = '\n'.join(trs)

# ---------- fig. 3: the best model per material ----------
def capnums(v, w):
    g = gap_small(v)
    return (f'ギャップ {g:.3f} eV · 誤差 {fmt(dgap(v), sign=True)} · ' if v.get('insulator') else '金属 · ') + f'最大 {w:.3f}'
gal = []
for g, gname, ms in GROUPS:
    items = []
    for m in ms.split():
        if m not in R['base']: continue
        w, k = BEST[m]; src, ns = FIG3[m]
        items.append(f'<figure class="card{" wide" if ns > 1 else ""}" id="fig-{m}" data-sev="{sev(w)}">'
                     f'<button class="zoom" aria-label="{m} を拡大" data-src="{src}" data-cap="{m} {LAB[k]}"><img loading="lazy" src="{src}" alt="{m}: DFT（灰の線）と MLO（赤の ×）"></button>'
                     f'<figcaption><span class="nm">{m}</span><span class="mdl{" alt" if k != "base" else ""}">{LAB[k]}</span><span class="chip {sev(w)}">{sev(w)}</span>'
                     f'<span class="meta">{capnums(R[k][m], w)}</span></figcaption></figure>')
    if items: gal.append(f'<section class="grp"><h3>{gname}</h3><div class="grid">{"".join(items)}</div></section>')

# ---------- fig. 4: models 2 and 3 against the baseline ----------
NOTE4 = {
 'SiO2c': 'クリストバライトの Si は隙間の多いダイヤモンド網で、伝導帯の底は網の空隙に広がった状態。原子の上の関数（1 の EH、2 の EH2）だけでは表しきれず、'
          '伝導帯に橙の点（模型に無い DFT の帯）が残る。空隙 2 か所（立方体の単位で ½(111) と ¾(111)、r = 2.6 a.u.）に z = 0 の球を置き、その s,p を模型に入れると（3）、伝導帯まで DFT に重なる。',
}
def w3(m, k): return fmt(W_[k].get(m))
_e = evaluate([f'{W}/mlocheck_eh2main/sp/EuO'], f'{OUT}/euo_eh2o.json').get('EuO')
EUO_O = worst(_e) if _e else None
NOTE4['C'] = f'2 で {w3("C", "base")} → {w3("C", "eh2")} eV に良くなるが、good（0.02 以下）には届かない。'
NOTE4['Bi2Te3'] = (f'1 も 2 も fair の下の方（{w3("Bi2Te3", "base")}、{w3("Bi2Te3", "eh2")} eV）。抜けた帯・余計な帯は無く（橙の点・黒丸はほとんど無い）、'
                   '窓の中の帯全体に小さなずれが広がるだけで、問題は無い（user 2026-10-01 18:1x）。スピン軌道は <code>job_mlo_soc</code>（摂動）。')
# the problem cases (user 2026-10-01 18:1x: "Cu 以下は問題のあるケース", "ボトムに"), fig. 5 at the bottom of the page
PROBLEM = {
 'Cu': f'<b>2 で壊れる</b>（{w3("Cu", "base")} → {w3("Cu", "eh2")} eV）。単体の金属なので、陽イオンの規則では全原子に EH2 が入る。'
       'Γ–X の 2〜3 点の k だけで MLO の帯が E_F + 0.3 eV に集まり（黒丸）、その k の DFT の帯が模型から抜ける（橙の点）。ほかの k は重なる。'
       '決めた既定（遷移金属・4f・5f 以外に EH2 の s,p）では Cu に EH2 は入らない。原因は未確認（同じ原子の EH と EH2 がほぼ一次従属になっているのではと疑っている）。',
 'Ni': f'<b>2 で壊れる</b>（{w3("Ni", "base")} → {w3("Ni", "eh2")} eV）。Cu と同じ壊れ方で、Γ–X の 2〜3 点の k だけ（両方のスピン）。決めた既定では Ni に EH2 は入らない。',
}
def cmpfig(m, prob=False):
    src, ns = FIG4[m]
    nums = ' / '.join(f'{LAB[k].split(".")[0]}: {fmt(W_[k][m])}' for k in ('base', 'eh2', 'es') if W_[k][m] is not None)
    stop = '（2 は止まった）' if m in fail_eh2 and W_['eh2'][m] is None else ''
    note = PROBLEM.get(m) if prob else NOTE4.get(m)
    cls = ' lead' if m in ES else (' prob' if prob else '')
    return (f'<figure class="cmp{cls}" id="cmp-{m}"><figcaption><b>{m}</b> · 誤差の最大値 {nums} eV{stop}'
            + (f'<p>{note}</p>' if note else '') + '</figcaption>'
            f'<button class="zoom" aria-label="{m} を拡大" data-src="{src}" data-cap="{m}"><img loading="lazy" src="{src}" alt="{m}: 模型 1・2・3 の比較"></button></figure>')
cmp = [cmpfig(m) for m in cmp4 if m not in PROBLEM]
cmp_prob = [cmpfig(m, True) for m in PROBLEM if m in FIG4]
# a problem case with no figure: model 2 could not be made (user 2026-10-01 18:2x: "EuO も書いておいて、図は無くてもいい")
PROBLEM_NOFIG = {
 'EuO': f'<b>2 で模型が作れない</b>（1: {w3("EuO", "base")} eV、2: 止まった）。陽イオンの規則では Eu（4f）の s,p に EH2 が入り、<code>job_mlo</code> の '
        '<code>Hreduction</code> が止まる: 「PMT completeness loss too large: band 26 dev= −0.012970」（2 番目の q 点、spin 1）。'
        '同じ原子の EH と EH2 をどちらもシードにすると、PMT の固有ベクトルが一次従属に近い向きを落とす（<code>zhev_tk4</code> の oveps）ので、'
        'シードのノルムが減る。減りが 1 % を超えると止める（<code>m_hreduction</code> の NormalizationCheck、1 % 未満なら規格化し直して進む）。'
        f'決めた既定（遷移金属・4f・5f 以外に EH2 の s,p）では Eu に EH2 は入らない。O だけに入れた模型は {fmt(EUO_O)} eV（good）で、止まらない。',
}
cmp_prob += [f'<div class="cmp prob" id="cmp-{m}"><p><b>{m}</b> · 図なし</p><p>{s}</p></div>' for m, s in PROBLEM_NOFIG.items()]

# materials whose best model is 2 or 3 by less than 0.001 eV (the extra seeds lower the error a little anyway)
tiny = [m for m in allm if BEST[m][1] != 'base' and W_['base'][m] - BEST[m][0] < 0.001]
TINY = '・'.join(tiny) + f'（{len(tiny)} 物質）'
page = open(f'{HERE}/page_template.html').read()
REP = {'@@CHART1@@': chart1, '@@CHART2@@': chart2, '@@TABLE@@': table, '@@CMP@@': '\n'.join(cmp), '@@CMPPROB@@': '\n'.join(cmp_prob), '@@NPROB@@': str(len(cmp_prob)), '@@GALLERY@@': '\n'.join(gal),
       '@@COUNTS@@': counts_tbl, '@@TINY@@': TINY, '@@SIO2_1@@': fmt(W_['base'].get('SiO2c'), 2), '@@SIO2_2@@': fmt(W_['eh2'].get('SiO2c'), 2), '@@N@@': str(len(allm)), '@@NCMP@@': str(len(cmp)), '@@REV@@': rev, '@@DATE@@': NOW}
for k, val in REP.items(): page = page.replace(k, val)
open(f'{OUT}/index.html', 'w').write(page)
print('wrote', OUT, len(allm), 'materials; base', CNT['base'], 'eh2', CNT['eh2'], 'best', CNT['best'], 'missing', missing,
      'fig4', cmp4, 'eh2 fail', fail_eh2)
