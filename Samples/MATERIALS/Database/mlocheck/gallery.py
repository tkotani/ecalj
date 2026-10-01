#!/usr/bin/env python3
# 2026-10-01 mlocheck: build the page of the MLO test (index.html + img/) in OUT from
#   mlocheck.json (mlo_bandcheck.py --json over the material dirs), variants.json (same over ../mlocheck_variants),
#   materials.json (prep.py). Figures: <Material>/mlo_<Material>.png and ../mlocheck_variants/<v>/mlo_<v>.png.
import json, os, shutil, sys, html, math, datetime, subprocess
W = os.path.dirname(os.path.abspath(__file__))
V = os.path.join(os.path.dirname(W), 'mlocheck_variants')
OUT = sys.argv[1]
os.makedirs(os.path.join(OUT, 'img'), exist_ok=True)
res = json.load(open(os.path.join(W, 'mlocheck.json')))
var = json.load(open(os.path.join(W, 'variants.json')))
mats = json.load(open(os.path.join(W, 'materials.json')))
rev = subprocess.run(['git', '-C', os.path.expanduser('~/ecalj'), 'log', '-1', '--format=%h'], capture_output=True, text=True).stdout.strip()

GROUPS = [
 ('sp', 'sp semiconductors', 'C Si Ge Sn 3cSiC 2hSiC 4hSiC AlN AlNzb AlP AlAs AlSb GaN GaNzb GaP GaAs GaSb InN InNzb InP InAs InSb MgS MgSe MgTe'),
 ('d10', 'Zn, Cd, Hg and Pb compounds', 'ZnO ZnS ZnSe ZnTe wZnS CdO CdS CdSe CdTe wCdS HgO HgS HgSe HgTe PbS PbTe'),
 ('ox', 'oxides and perovskites', 'MgO SiO2c ZrO2 HfO2 SrTiO3 SrVO3 BaTiO3 LaGaO3 La2CuO4'),
 ('mag', 'metals and magnets', 'Li Cu Fe Ni YMn2 MnO NiO'),
 ('f', '4f (spdf model)', 'Ce EuO EuS EuSe EuTe'),
 ('soc', 'spin-orbit (job_mlo_soc)', 'GaAs_so Bi2Te3'),
 ('sl', 'superlattices', 'InAsGaSb_n4 InAsGaSb_n10'),
]
group_of = {m: g for g, _, ms in GROUPS for m in ms.split()}

def worst(v):
    w = [v.get('rms_m2d', 0), v.get('rms_d2m', 0)]
    if v.get('insulator'): w.append(abs(v['gapM'] - v['gapD']))
    return max(x for x in w if x == x)

def sev(x):   # good <= 0.02, fair <= 0.05, marginal <= 0.1, poor > 0.1 eV (user 2026-10-01 14:3x: four classes)
    return 'good' if x <= 0.02 else ('fair' if x <= 0.05 else ('marginal' if x <= 0.1 else 'poor'))

def fmt(x, n=3, sign=False):
    if x is None or (isinstance(x, float) and x != x): return '—'
    return (f'{x:+.{n}f}' if sign else f'{x:.{n}f}')

# copy figures
figs = {}
for m in res:
    src = os.path.join(W, m, f'mlo_{m}.png')
    if os.path.exists(src):
        shutil.copy(src, os.path.join(OUT, 'img', f'{m}.png')); figs[m] = f'img/{m}.png'
vfigs = {}
for v in var:
    src = os.path.join(V, v, f'mlo_{v}.png')
    if os.path.exists(src):
        shutil.copy(src, os.path.join(OUT, 'img', f'v_{v}.png')); vfigs[v] = f'img/v_{v}.png'

allm = [m for g, _, ms in GROUPS for m in ms.split()]
missing = [m for m in allm if m not in res]
rows = [(m, res[m]) for m in allm if m in res]
nsev = {'good': 0, 'fair': 0, 'marginal': 0, 'poor': 0}
for m, v in rows: nsev[sev(worst(v))] += 1

# ---------- chart 1: gap error against the DFT gap (insulators), arrows to the fixed variants ----------
FIX = {  # base -> variant that fixes it (and what was added)
 'GaN': ('GaN_lm3d', 'Ga 3d LO'), 'GaNzb': ('GaNzb_lm3d', 'Ga 3d LO'), 'MgS': ('MgS_lm2sp', 'EH2 s,p'), 'MgSe': ('MgSe_lm2sp', 'EH2 s,p'),
 'MgTe': ('MgTe_lm2sp', 'EH2 s,p'), 'CdTe': ('CdTe_lm2sp', 'EH2 s,p'), 'ZnTe': ('ZnTe_lm2sp', 'EH2 s,p'), 'AlN': ('AlN_lm2sp', 'EH2 s,p'),
 'SiO2c': ('SiO2c_lm2sp', 'EH2 s,p'), 'LaGaO3': ('LaGaO3_lm3semi', 'Ga 3d, La 5p LO'), 'SrTiO3': ('SrTiO3_both', 'Sr 4p LO + EH2 s,p'),
 'InN': ('InN_lm3d', 'In 4d LO'), 'InNzb': ('InNzb_lm3d', 'In 4d LO'), 'La2CuO4': ('La2CuO4_lm3semi', 'La 5p LO'),
 'EuO': ('EuO_lm3semi', 'Eu 5p LO'), 'SrVO3': ('SrVO3_lm3semi', 'Sr 4s,4p and V 3p LO'),
}
X0, X1, Y0, Y1 = 0.0, 6.0, -0.10, 0.45
PW, PH, ML, MR, MT, MB = 720, 360, 56, 18, 16, 44
def sx(x): return ML + (min(max(x, X0), X1) - X0) / (X1 - X0) * (PW - ML - MR)
def sy(y): return MT + (Y1 - min(max(y, Y0), Y1)) / (Y1 - Y0) * (PH - MT - MB)
svg = [f'<svg viewBox="0 0 {PW} {PH}" role="img" aria-labelledby="c1t" class="chart"><title id="c1t">Gap of the MLO model minus the DFT gap, against the DFT gap</title>']
for yt in (-0.1, 0, 0.1, 0.2, 0.3, 0.4):
    svg.append(f'<line x1="{ML}" x2="{PW-MR}" y1="{sy(yt):.1f}" y2="{sy(yt):.1f}" class="{"zero" if yt == 0 else "grid"}"/>'
               f'<text x="{ML-8}" y="{sy(yt)+4:.1f}" class="tick" text-anchor="end">{yt:+.1f}</text>')
for xt in range(0, 7):
    svg.append(f'<line x1="{sx(xt):.1f}" x2="{sx(xt):.1f}" y1="{MT}" y2="{PH-MB}" class="grid"/>'
               f'<text x="{sx(xt):.1f}" y="{PH-MB+18}" class="tick" text-anchor="middle">{xt}</text>')
svg.append(f'<text x="{(ML+PW-MR)/2}" y="{PH-6}" class="axis" text-anchor="middle">DFT gap along the path (eV)</text>')
svg.append(f'<text x="14" y="{(MT+PH-MB)/2}" class="axis" text-anchor="middle" transform="rotate(-90 14 {(MT+PH-MB)/2})">MLO gap − DFT gap (eV)</text>')
for m, v in rows:
    if not v.get('insulator'): continue
    d = v['gapM'] - v['gapD']; x, y = sx(v['gapD']), sy(d)
    if m in FIX and FIX[m][0] in var and var[FIX[m][0]].get('insulator'):
        vv = var[FIX[m][0]]; y2 = sy(vv['gapM'] - vv['gapD'])
        svg.append(f'<line x1="{x:.1f}" y1="{y:.1f}" x2="{x:.1f}" y2="{y2+5 if y2 > y else y2-5:.1f}" class="arrow" marker-end="url(#ah)"/>')
        svg.append(f'<circle cx="{x:.1f}" cy="{y2:.1f}" r="4.5" class="pt fixed"><title>{m} + {FIX[m][1]}: {vv["gapM"]-vv["gapD"]:+.3f} eV</title></circle>')
    clip = ' clipped' if d > Y1 else ''
    svg.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" class="pt {sev(worst(v))}{clip}"><title>{m}: DFT {v["gapD"]:.3f} eV, MLO {v["gapM"]:.3f} eV, diff {d:+.3f} eV</title></circle>')
    if abs(d) > 0.06:
        right = x > 0.75 * PW   # near the right edge the label goes to the left of the point
        svg.append(f'<text x="{x-8 if right else x+8:.1f}" y="{y+4:.1f}" class="lab" text-anchor="{"end" if right else "start"}">{m}{" (+%.2f eV, off scale)" % d if d > Y1 else ""}</text>')
svg.append('<defs><marker id="ah" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" class="ahead"/></marker></defs></svg>')
chart1 = '\n'.join(svg)

# ---------- chart 2: rms MLO->DFT, sorted, log scale ----------
srt = sorted(rows, key=lambda r: worst(r[1]))
BW, BH, bl, bb, bt = max(720, 11 * len(srt) + 70), 300, 56, 92, 14
L0, L1 = math.log10(0.0005), math.log10(2.0)
def by(val): return bt + (L1 - math.log10(min(max(val, 0.0005), 2.0))) / (L1 - L0) * (BH - bt - bb)
sv = [f'<svg viewBox="0 0 {BW} {BH}" role="img" aria-labelledby="c2t" class="chart"><title id="c2t">Worst of the band rms and the gap error, per material (log scale)</title>']
for val, lab in ((0.001, '0.001'), (0.01, '0.01'), (0.02, '0.02'), (0.05, '0.05'), (0.1, '0.1'), (1, '1')):
    cls = 'thr' if val in (0.02, 0.05, 0.1) else 'grid'
    sv.append(f'<line x1="{bl}" x2="{BW-8}" y1="{by(val):.1f}" y2="{by(val):.1f}" class="{cls}"/><text x="{bl-6}" y="{by(val)+4:.1f}" class="tick" text-anchor="end">{lab}</text>')
step = (BW - bl - 8) / len(srt)
for i, (m, v) in enumerate(srt):
    w = worst(v); x = bl + i * step
    sv.append(f'<rect x="{x+1:.1f}" y="{by(w):.1f}" width="{max(step-2,2):.1f}" height="{BH-bb-by(w):.1f}" rx="2" class="bar {sev(w)}"><title>{m}: {w:.3f} eV</title></rect>')
    sv.append(f'<text x="{x+step/2:.1f}" y="{BH-bb+8}" class="blab" transform="rotate(60 {x+step/2:.1f} {BH-bb+8})">{m}</text>')
sv.append(f'<text x="14" y="{(bt+BH-bb)/2}" class="axis" text-anchor="middle" transform="rotate(-90 14 {(bt+BH-bb)/2})">eV</text></svg>')
chart2 = '\n'.join(sv)

# ---------- radial functions added: list (placed above the figures) ----------
import re as _re
def nmlo(d):
    p = os.path.join(d, 'lmlo')
    if not os.path.exists(p): return None
    m = _re.search(r'HamRsMTO=\s*(\d+)', open(p, errors='replace').read())
    return int(m.group(1)) if m else None
KIND = {'lm2sp': ('<code>mlo_lm2</code>', '第 2 の smooth Hankel 関数（EH2）の s,p を全原子に', '(b)'),
        'lm2s':  ('<code>mlo_lm2</code>', 'EH2 の s を全原子に', '(b)'),
        'lm3d':  ('<code>mlo_lm3</code>', '陽イオンの半内殻 d の局所軌道', '(a)'),
        'lm3semi': ('<code>mlo_lm3</code>', '半内殻の局所軌道（s,p,d）すべて', '(a)'),
        'both':  ('<code>mlo_lm2</code> と <code>mlo_lm3</code>', 'EH2 の s,p と半内殻の局所軌道', '(a)+(b)')}
rad = []
for base, (vn, what) in FIX.items():
    if base not in res or vn not in var: continue
    k = vn.split('_', 1)[1] if '_' in vn else ''
    key, desc, typ = KIND.get(k, ('', '', ''))
    n0, n1 = nmlo(os.path.join(W, base)), nmlo(os.path.join(V, vn))
    rad.append(f'<tr><td class="name">{base}</td><td>{typ}</td><td>{key}</td><td>{html.escape(what.replace(" LO", " の局所軌道").replace(" and ", " と ").replace(" + ", " と ").replace("EH2 s,p", "EH2 の s,p（全原子）"))}</td>'
               f'<td class="num">{n0} → {n1}</td><td class="num">{worst(res[base]):.3f} → {worst(var[vn]):.3f}</td></tr>')
bad_rows = []
for vn, note in (('Cu_both', '止まらずに模型が壊れる'), ('EuO_both', '<code>Hreduction: PMT completeness loss too large</code> で止まる'),
                 ('Si_both', ''), ('GaAs_both', ''), ('Fe_both', ''), ('NiO_both', ''), ('ZnO_both', '')):
    base = vn.split('_')[0]
    if base not in res: continue
    after = f'{worst(var[vn]):.3f}' if vn in var else '止まった'
    n0, n1 = nmlo(os.path.join(W, base)), nmlo(os.path.join(V, vn))
    bad_rows.append(f'<tr><td class="name">{base}</td><td class="num">{n0} → {n1 if n1 else "—"}</td>'
                    f'<td class="num">{worst(res[base]):.3f} → {after}</td><td>{note}</td></tr>')
radial = ('<table><thead><tr><th>物質</th><th>型</th><th>キー</th><th>足したもの</th><th>MLO の数</th><th>最悪値 (eV)</th></tr></thead><tbody>'
          + ''.join(rad) + '</tbody></table>')
radial_bad = ('<table><thead><tr><th>物質</th><th>MLO の数</th><th>最悪値 (eV)</th><th></th></tr></thead><tbody>'
              + ''.join(bad_rows) + '</tbody></table>')

# ---------- after the radial functions: the verdicts with the fixed models in place of the default ones ----------
best = {m: (worst(var[FIX[m][0]]) if (m in FIX and FIX[m][0] in var) else worst(v)) for m, v in res.items()}
nbest = {'good': 0, 'fair': 0, 'marginal': 0, 'poor': 0}
for w in best.values(): nbest[sev(w)] += 1
left = sorted(((w, m) for m, w in best.items() if w > 0.02), reverse=True)
after = ('<table><thead><tr><th></th><th>good（≤ 0.02 eV）</th><th>fair（≤ 0.05 eV）</th><th>marginal（≤ 0.1 eV）</th><th>poor（&gt; 0.1 eV）</th></tr></thead><tbody>'
         f'<tr><td>既定の模型</td><td class="num">{nsev["good"]}</td><td class="num">{nsev["fair"]}</td><td class="num">{nsev["marginal"]}</td><td class="num">{nsev["poor"]}</td></tr>'
         f'<tr><td>足した後</td><td class="num">{nbest["good"]}</td><td class="num">{nbest["fair"]}</td><td class="num">{nbest["marginal"]}</td><td class="num">{nbest["poor"]}</td></tr>'
         '</tbody></table>'
         '<p class="note">足した後も good にならないもの（最悪値 eV）: ' + '、'.join(f'{m} {w:.3f}' for w, m in left) +
         '。SiO₂ のほかは何も足して試していない（C・AlSb・InSb・Sn・GaSb・SiC は s,p の半導体で、(b) の型が軽く出ている可能性がある）。'
         '何を足すかは型を見て物質ごとに選んだもので、自動ではない。</p>')
afterline = f'動径関数を足した後（下の「動径関数を足すとは」）: good {nbest["good"]}、fair {nbest["fair"]}、poor {nbest["poor"]}。'

# ---------- table ----------
trs = []
for m, v in rows:
    s = sev(worst(v)); md = mats[m]
    kind = 'insulator' if v.get('insulator') else 'metal'
    extra = []
    if md['nspin'] == 2 and not md['so']: extra.append('spin')
    if md['so']: extra.append('SOC')
    if md['f']: extra.append('4f')
    gd = v.get('gapD'); gm = v.get('gapM')
    trs.append(f'<tr data-sev="{s}"><td class="name"><a href="#fig-{m}">{m}</a></td><td>{group_of.get(m,"")}</td><td>{md["nsite"]}</td>'
               f'<td>{kind}{" · "+", ".join(extra) if extra else ""}</td>'
               f'<td class="num">{fmt(v.get("gap_mesh")) if v.get("insulator") else "—"}</td><td class="num">{fmt(gd)}</td><td class="num">{fmt(gm)}</td>'
               f'<td class="num">{fmt((gm-gd) if gd is not None else None, sign=True)}</td>'
               f'<td class="num">{fmt(v.get("dVBM"), sign=True)}</td><td class="num">{fmt(v.get("dCBM"), sign=True)}</td>'
               f'<td class="num">{fmt(v["rms_m2d"])}</td><td class="num">{fmt(v["max_m2d"],2)}</td><td class="num">{fmt(v["rms_d2m"])}</td>'
               f'<td><span class="chip {s}">{s}</span></td></tr>')
table = '\n'.join(trs)

# ---------- fixes ----------
fx = []
for base, (vn, what) in FIX.items():
    if base not in res or vn not in var: continue
    b, a = res[base], var[vn]
    fx.append(f'<figure class="pair"><figcaption><b>{base}</b> + {html.escape(what)}: worst {worst(b):.3f} → {worst(a):.3f} eV</figcaption>'
              f'<div class="two"><img loading="lazy" src="{figs.get(base,"")}" alt="{base}, default model"><img loading="lazy" src="{vfigs.get(vn,"")}" alt="{base} with {html.escape(what)}"></div></figure>')
if 'HfO2_spd' in var and 'HfO2' in res:
    fx.append(f'<figure class="pair"><figcaption><b>HfO2</b>: s,p,d only (left, the first run) → Hf 4f in the model (right, in the table): '
              f'worst {worst(var["HfO2_spd"]):.3f} → {worst(res["HfO2"]):.3f} eV</figcaption>'
              f'<div class="two"><img loading="lazy" src="{vfigs.get("HfO2_spd","")}" alt="HfO2 spd"><img loading="lazy" src="{figs.get("HfO2","")}" alt="HfO2 spdf"></div></figure>')
bad = []
for vn, what in (('Cu_both', 'EH2 s,p on Cu: the model breaks without a stop'),):
    if vn in var and vn in vfigs:
        bad.append(f'<figure class="pair"><figcaption><b>Cu</b> + EH2 s,p: worst {worst(res["Cu"]):.3f} → {worst(var[vn]):.3f} eV. {html.escape(what)}</figcaption>'
                   f'<div class="two"><img loading="lazy" src="{figs["Cu"]}" alt="Cu default"><img loading="lazy" src="{vfigs[vn]}" alt="Cu with EH2 s,p"></div></figure>')

# ---------- gallery ----------
gal = []
for g, gname, ms in GROUPS:
    items = []
    for m in ms.split():
        if m not in res: continue
        v = res[m]; s = sev(worst(v))
        gap = f'gap {v["gapD"]:.2f} → {v["gapM"]:.2f} eV' if v.get('insulator') else 'metal'
        items.append(f'<figure class="card" id="fig-{m}" data-sev="{s}"><button class="zoom" aria-label="Enlarge {m}" data-src="{figs.get(m,"")}" data-cap="{m}">'
                     f'<img loading="lazy" src="{figs.get(m,"")}" alt="{m}: DFT bands (grey) and MLO model (red)"></button>'
                     f'<figcaption><span class="nm">{m}</span><span class="chip {s}">{s}</span><span class="meta">{gap} · rms {v["rms_m2d"]:.3f} / {v["rms_d2m"]:.3f}</span></figcaption></figure>')
    if items: gal.append(f'<section class="grp"><h3>{gname}</h3><div class="grid">{"".join(items)}</div></section>')

page = open(os.path.join(W, 'page_template.html')).read()
for k, val in {'@@AGOOD@@': str(nbest['good']), '@@AFAIR@@': str(nbest['fair']), '@@AMARG@@': str(nbest['marginal']), '@@APOOR@@': str(nbest['poor']), '@@NMARG@@': str(nsev['marginal']), '@@AFTER@@': after, '@@AFTERLINE@@': afterline, '@@RADIAL@@': radial, '@@RADIALBAD@@': radial_bad, '@@CHART1@@': chart1, '@@CHART2@@': chart2, '@@TABLE@@': table, '@@FIXES@@': '\n'.join(fx), '@@BAD@@': '\n'.join(bad),
               '@@GALLERY@@': '\n'.join(gal), '@@N@@': str(len(rows)), '@@NGOOD@@': str(nsev['good']), '@@NFAIR@@': str(nsev['fair']),
               '@@NPOOR@@': str(nsev['poor']), '@@REV@@': rev, '@@MISSING@@': ', '.join(missing) if missing else 'none',
               '@@DATE@@': datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}.items():
    page = page.replace(k, val)
open(os.path.join(OUT, 'index.html'), 'w').write(page)
print('wrote', OUT, len(rows), 'materials', nsev, 'missing', missing, 'figs', len(figs), 'variant figs', len(vfigs))
