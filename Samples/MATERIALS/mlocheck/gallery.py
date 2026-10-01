#!/usr/bin/env python3
# 2026-10-01 mlocheck page, version 2 (15:xx): the baseline is the model with the semicore local orbitals added
# automatically (m_HamPMT: band top above EF - 17 eV, weight over the species). Columns per material: the old default
# (until 2026-10-01 14:4x), the baseline, the baseline + EH2 s,p on the cations, and empty spheres (SiO2 only).
#   python3 gallery.py <OUT>        (writes OUT/index.html and OUT/img/)
# Inputs: JSON of SRC/exec/mlo_bandcheck.py --json for each set (made here), the figures of each set.
# The first version (old default and hand-picked variants) is gallery_v1.py / page_template_v1.html.
import json, os, shutil, sys, html, math, datetime, subprocess, re, glob
HERE = os.path.dirname(os.path.abspath(__file__))
W = os.path.dirname(HERE)                                      # ~/work
SETS = {'old': f'{W}/mlocheck_20261001', 'base': f'{W}/mlocheck_auto_20261001', 'eh2': f'{W}/mlocheck_eh2cat'}
ES = {'SiO2c': f'{W}/mlocheck_es/SiO2c_ES'}
OUT = sys.argv[1]
os.makedirs(os.path.join(OUT, 'img'), exist_ok=True)
mats = json.load(open(f'{HERE}/materials.json'))
CHECK = os.path.expanduser('~/ecalj/SRC/exec/mlo_bandcheck.py')
rev = subprocess.run(['git', '-C', os.path.expanduser('~/ecalj'), 'log', '-1', '--format=%h'], capture_output=True, text=True).stdout.strip()

def evaluate(dirs, out):
    dirs = [d for d in dirs if os.path.exists(os.path.join(d, 'band_MLO_spin1.dat')) and os.path.getsize(os.path.join(d, 'band_MLO_spin1.dat')) > 0]
    if dirs: subprocess.run(['python3', CHECK, *dirs, '--json', out], capture_output=True, text=True)
    return json.load(open(out)) if dirs and os.path.exists(out) else {}

R = {k: evaluate([f'{p}/{m}' for m in mats], f'{OUT}/{k}.json') for k, p in SETS.items()}
r_es = evaluate(list(ES.values()), f'{OUT}/es.json')
R['es'] = {m: r_es[os.path.basename(d)] for m, d in ES.items() if os.path.basename(d) in r_es}
json.dump(R, open(f'{OUT}/all.json', 'w'), indent=1)

def nmlo(d):
    p = os.path.join(d, 'lmlo')
    if not os.path.exists(p): return None
    m = re.findall(r'HamRsMTO=\s*(\d+)', open(p, errors='replace').read())
    return int(m[-1]) if m else None
NM = {'base': {m: nmlo(f"{SETS['base']}/{m}") for m in mats}, 'eh2': {m: nmlo(f"{SETS['eh2']}/{m}") for m in mats},
      'es': {m: nmlo(d) for m, d in ES.items()}}
fail_eh2 = {}
st = f"{SETS['eh2']}/status.log"
if os.path.exists(st):
    for l in open(st):
        t = l.split()
        if len(t) > 3 and t[3] == 'FAIL': fail_eh2[t[2]] = ' '.join(t[4:])

def worst(v):
    if not v: return None
    w = [v['rms_m2d'], v['rms_d2m']]
    if v.get('insulator'): w.append(abs(v['gapM'] - v['gapD']))
    return max(w)
def sev(x):   # good <= 0.02, fair <= 0.05, marginal <= 0.1, poor > 0.1 eV (user 2026-10-01 14:3x)
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
def best_of(m):
    c = [(W_['base'][m], 'base')] + [(W_[k][m], k) for k in ('eh2', 'es') if W_[k][m] is not None]
    return min(c)
BEST = {m: best_of(m) for m in allm}
def counts(k):
    c = {x: 0 for x in CL}
    for m in allm:
        w = BEST[m][0] if k == 'best' else W_[k][m]
        if w is not None: c[sev(w)] += 1
    return c
CNT = {k: counts(k) for k in ('old', 'base', 'eh2', 'best')}

def fmt(x, n=3, sign=False):
    if x is None or (isinstance(x, float) and x != x): return '—'
    return f'{x:+.{n}f}' if sign else f'{x:.{n}f}'
def chip(x): return f'<span class="chip {sev(x)}">{sev(x)}</span>' if x is not None else '—'

# ---------- figures ----------
def cp(src, dst):
    if os.path.exists(src): shutil.copy(src, f'{OUT}/img/{dst}'); return f'img/{dst}'
    return ''
FIG = {m: cp(f"{SETS['base']}/{m}/mlo_{m}.png", f'{m}.png') for m in allm}
FIG_EH2 = {m: cp(f"{SETS['eh2']}/{m}/mlo_{m}_eh2cat.png", f'{m}_eh2cat.png') for m in allm}
FIG_ES = {m: cp(f'{d}/mlo_{os.path.basename(d)}.png', f'{m}_es.png') for m, d in ES.items()}

# ---------- chart 1: gap error of the baseline, arrows to the best option ----------
X0, X1, Y0, Y1 = 0.0, 6.0, -0.10, 0.45
PW, PH, ML, MR, MT, MB = 720, 360, 56, 18, 16, 44
def sx(x): return ML + (min(max(x, X0), X1) - X0) / (X1 - X0) * (PW - ML - MR)
def sy(y): return MT + (Y1 - min(max(y, Y0), Y1)) / (Y1 - Y0) * (PH - MT - MB)
svg = [f'<svg viewBox="0 0 {PW} {PH}" role="img" aria-labelledby="c1t" class="chart"><title id="c1t">基準の模型のギャップの差</title>']
for yt in (-0.1, 0, 0.1, 0.2, 0.3, 0.4):
    svg.append(f'<line x1="{ML}" x2="{PW-MR}" y1="{sy(yt):.1f}" y2="{sy(yt):.1f}" class="{"zero" if yt == 0 else "grid"}"/>'
               f'<text x="{ML-8}" y="{sy(yt)+4:.1f}" class="tick" text-anchor="end">{yt:+.1f}</text>')
for xt in range(0, 7):
    svg.append(f'<line x1="{sx(xt):.1f}" x2="{sx(xt):.1f}" y1="{MT}" y2="{PH-MB}" class="grid"/>'
               f'<text x="{sx(xt):.1f}" y="{PH-MB+18}" class="tick" text-anchor="middle">{xt}</text>')
svg.append(f'<text x="{(ML+PW-MR)/2}" y="{PH-6}" class="axis" text-anchor="middle">DFT のギャップ（対称線の上、eV）</text>')
svg.append(f'<text x="14" y="{(MT+PH-MB)/2}" class="axis" text-anchor="middle" transform="rotate(-90 14 {(MT+PH-MB)/2})">MLO − DFT（eV）</text>')
for m in allm:
    v = R['base'][m]
    if not v.get('insulator'): continue
    d = v['gapM'] - v['gapD']; x, y = sx(v['gapD']), sy(d)
    w0, k = BEST[m]
    if k != 'base' and R[k][m].get('insulator'):
        vv = R[k][m]; y2 = sy(vv['gapM'] - vv['gapD'])
        if abs(y2 - y) > 6:
            svg.append(f'<line x1="{x:.1f}" y1="{y:.1f}" x2="{x:.1f}" y2="{y2+5 if y2 > y else y2-5:.1f}" class="arrow" marker-end="url(#ah)"/>')
        svg.append(f'<circle cx="{x:.1f}" cy="{y2:.1f}" r="4.5" class="pt fixed"><title>{m} {"+ 陽イオン EH2" if k == "eh2" else "+ 空格子球"}: {vv["gapM"]-vv["gapD"]:+.3f} eV</title></circle>')
    svg.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" class="pt {sev(W_["base"][m])}"><title>{m}: DFT {v["gapD"]:.3f}, MLO {v["gapM"]:.3f}, 差 {d:+.3f} eV</title></circle>')
    if abs(d) > 0.06:
        right = x > 0.75 * PW
        svg.append(f'<text x="{x-8 if right else x+8:.1f}" y="{y+4:.1f}" class="lab" text-anchor="{"end" if right else "start"}">{m}{" (+%.2f eV、枠外)" % d if d > Y1 else ""}</text>')
svg.append('<defs><marker id="ah" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" class="ahead"/></marker></defs></svg>')
chart1 = '\n'.join(svg)

# ---------- chart 2: worst value of the baseline, sorted, log scale ----------
srt = sorted(allm, key=lambda m: W_['base'][m])
BW, BH, bl, bb, bt = max(720, 11 * len(srt) + 70), 300, 56, 92, 14
L0, L1 = math.log10(0.0005), math.log10(5.0)
def by(val): return bt + (L1 - math.log10(min(max(val, 0.0005), 5.0))) / (L1 - L0) * (BH - bt - bb)
sv = [f'<svg viewBox="0 0 {BW} {BH}" role="img" aria-labelledby="c2t" class="chart"><title id="c2t">基準の模型の物質ごとの最悪値</title>']
for val, lab in ((0.001, '0.001'), (0.01, '0.01'), (0.02, '0.02'), (0.05, '0.05'), (0.1, '0.1'), (1, '1')):
    cls = 'thr' if val in (0.02, 0.05, 0.1) else 'grid'
    sv.append(f'<line x1="{bl}" x2="{BW-8}" y1="{by(val):.1f}" y2="{by(val):.1f}" class="{cls}"/><text x="{bl-6}" y="{by(val)+4:.1f}" class="tick" text-anchor="end">{lab}</text>')
step = (BW - bl - 8) / len(srt)
for i, m in enumerate(srt):
    w = W_['base'][m]; x = bl + i * step
    sv.append(f'<rect x="{x+1:.1f}" y="{by(w):.1f}" width="{max(step-2,2):.1f}" height="{BH-bb-by(w):.1f}" rx="2" class="bar {sev(w)}"><title>{m}: {w:.3f} eV</title></rect>')
    sv.append(f'<text x="{x+step/2:.1f}" y="{BH-bb+8}" class="blab" transform="rotate(60 {x+step/2:.1f} {BH-bb+8})">{m}</text>')
sv.append(f'<text x="14" y="{(bt+BH-bb)/2}" class="axis" text-anchor="middle" transform="rotate(-90 14 {(bt+BH-bb)/2})">eV</text></svg>')
chart2 = '\n'.join(sv)

# ---------- the counts table ----------
def crow(label, c): return f'<tr><td>{label}</td>' + ''.join(f'<td class="num">{c[x]}</td>' for x in CL) + '</tr>'
counts_tbl = ('<table class="cnt"><thead><tr><th>模型</th><th class="num">good<br>≤ 0.02</th><th class="num">fair<br>≤ 0.05</th>'
              '<th class="num">marginal<br>≤ 0.1</th><th class="num">poor<br>&gt; 0.1 eV</th></tr></thead><tbody>'
              + crow('旧既定（〜2026-10-01 14:4x）', CNT['old']) + crow('<b>基準</b>: 自動の半内殻入り', CNT['base'])
              + crow('基準 + 陽イオンに EH2 の s,p', CNT['eh2']) + crow('物質ごとに一番良いもの（基準・陽イオン EH2・空格子球）', CNT['best'])
              + '</tbody></table>')

# ---------- pairs: materials whose verdict an option improves ----------
pairs = []
for m in allm:
    w0, k = BEST[m]
    if k == 'base' or sev(w0) == sev(W_['base'][m]): continue
    f2 = FIG_EH2[m] if k == 'eh2' else FIG_ES.get(m, '')
    what = '陽イオンに EH2 の s,p' if k == 'eh2' else '空隙に空格子球 2 つ（s,p）'
    pairs.append(f'<figure class="pair"><figcaption><b>{m}</b> + {what}: 最悪値 {W_["base"][m]:.3f} → {w0:.3f} eV</figcaption>'
                 f'<div class="two"><img loading="lazy" src="{FIG[m]}" alt="{m} 基準"><img loading="lazy" src="{f2}" alt="{m} {what}"></div></figure>')

# ---------- table 1 ----------
trs = []
for m in allm:
    v = R['base'][m]; md = mats[m]
    kind = '絶縁体' if v.get('insulator') else '金属'
    extra = []
    if md['nspin'] == 2 and not md['so']: extra.append('スピン')
    if md['so']: extra.append('SOC')
    if md['f']: extra.append('4f')
    gd = v.get('gapD'); gm = v.get('gapM')
    e2 = R['eh2'].get(m); es = R['es'].get(m)
    e2s = f'{fmt(W_["eh2"][m])} {chip(W_["eh2"][m])}' if e2 else ('止まった' if m in fail_eh2 else '—')
    ess = f'{fmt(W_["es"][m])} {chip(W_["es"][m])}' if es else ''
    trs.append(f'<tr data-sev="{sev(W_["base"][m])}"><td class="name"><a href="#fig-{m}">{m}</a></td><td class="num">{md["nsite"]}</td>'
               f'<td>{kind}{"・"+"・".join(extra) if extra else ""}</td>'
               f'<td class="num">{fmt(v.get("gap_mesh")) if v.get("insulator") else "—"}</td><td class="num">{fmt(gd)}</td>'
               f'<td class="num">{fmt(W_["old"].get(m))}</td>'
               f'<td class="num">{fmt((gm-gd) if gd is not None else None, sign=True)}</td>'
               f'<td class="num">{fmt(v["rms_m2d"])}</td><td class="num">{fmt(v["max_m2d"],2)}</td><td class="num">{fmt(v["rms_d2m"])}</td>'
               f'<td class="num">{fmt(W_["base"][m])}</td><td>{chip(W_["base"][m])}</td><td class="num">{NM["base"][m] or "—"}</td>'
               f'<td class="num">{e2s}</td><td class="num">{NM["eh2"][m] or "—"}</td><td class="num">{ess}</td></tr>')
table = '\n'.join(trs)

# ---------- gallery ----------
gal = []
for g, gname, ms in GROUPS:
    items = []
    for m in ms.split():
        if m not in R['base']: continue
        v = R['base'][m]; s = sev(W_['base'][m])
        gap = f'ギャップ {v["gapD"]:.2f} → {v["gapM"]:.2f} eV' if v.get('insulator') else '金属'
        items.append(f'<figure class="card" id="fig-{m}" data-sev="{s}"><button class="zoom" aria-label="{m} を拡大" data-src="{FIG[m]}" data-cap="{m}">'
                     f'<img loading="lazy" src="{FIG[m]}" alt="{m}: DFT（灰）と MLO（赤）"></button>'
                     f'<figcaption><span class="nm">{m}</span><span class="chip {s}">{s}</span><span class="meta">{gap} · rms {v["rms_m2d"]:.3f} / {v["rms_d2m"]:.3f}</span></figcaption></figure>')
    if items: gal.append(f'<section class="grp"><h3>{gname}</h3><div class="grid">{"".join(items)}</div></section>')

page = open(f'{HERE}/page_template.html').read()
REP = {'@@CHART1@@': chart1, '@@CHART2@@': chart2, '@@TABLE@@': table, '@@PAIRS@@': '\n'.join(pairs), '@@GALLERY@@': '\n'.join(gal),
       '@@COUNTS@@': counts_tbl, '@@N@@': str(len(allm)), '@@REV@@': rev, '@@DATE@@': datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}
for x in CL:
    REP[f'@@B_{x}@@'] = str(CNT['base'][x]); REP[f'@@S_{x}@@'] = str(CNT['best'][x])
for k, val in REP.items(): page = page.replace(k, val)
open(f'{OUT}/index.html', 'w').write(page)
print('wrote', OUT, len(allm), 'materials; base', CNT['base'], 'eh2', CNT['eh2'], 'best', CNT['best'], 'missing', missing,
      'eh2 done', len(R['eh2']), 'eh2 fail', fail_eh2)
