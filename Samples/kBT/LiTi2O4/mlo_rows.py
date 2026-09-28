#!/usr/bin/env python3
"""LiTi2O4 6^3: conventional (MTO) QSGW beside MLO-QSGW (chi~ frozen).
One ROW per QSGW iteration, growing downward as the chains proceed.
usage: mlo_rows.py ROOT OUT.png
ROOT holds  ref/bnd_lda.dat, ref/bnd_iter<N>.dat  and  frozen/bnd_lda.dat, frozen/bnd_iter<N>.dat"""
import sys, os, datetime, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

ROOT, OUT = sys.argv[1], sys.argv[2]
# argv[3] and argv[4] are comma-separated: one MLO chain per column, left to right
MLODIR = (sys.argv[3] if len(sys.argv) > 3 else 'frozen').split(',')
MLOLAB = (sys.argv[4] if len(sys.argv) > 4 else 'MLO-QSGW').split(',')
MLOCOL = (sys.argv[5].split(',') if len(sys.argv) > 5
          else ['tab:green', 'tab:red', 'tab:purple', 'tab:orange'])
NOREF = len(sys.argv) > 6 and sys.argv[6] == 'noref'   # no conventional chain for this mesh
COLS = ([] if NOREF else [('conventional QSGW (MTO), pwmode=1', f'{ROOT}/ref', 'tab:blue')]) + \
       [(MLOLAB[i] if i < len(MLOLAB) else d, f'{ROOT}/{d}', MLOCOL[i % len(MLOCOL)])
        for i, d in enumerate(MLODIR)]

def rd(f):
    xs, es, cx, ce = [], [], [], []
    for ln in open(f):
        if ln.startswith('#'): continue
        t = ln.split()
        if len(t) < 3:
            if cx: xs.append(cx); es.append(ce); cx = []; ce = []
            continue
        cx.append(float(t[1])); ce.append(float(t[2]))
    if cx: xs.append(cx); es.append(ce)
    return np.array(xs[0]), np.array(es)

## HILITE=1: colour only the two bands that carry the oscillation.  The plotted lines are energy-sorted
## at each k, so a band that rises through flat ones is split over several sorted lines; to follow it
## as one band, track by continuity: predict E(k_i) linearly from the two previous k and assign the
## sorted energies at k_i to the tracks with the Hungarian method.
HILITE = os.environ.get('HILITE', '0') == '1'
def track(S):
    from scipy.optimize import linear_sum_assignment
    nb, nk = S.shape; T = np.empty_like(S); T[:, 0] = S[:, 0]; T[:, 1] = S[:, 1]
    for i in range(2, nk):
        pred = 2*T[:, i-1] - T[:, i-2]
        r, c = linear_sum_assignment(np.abs(pred[:, None] - S[None, :, i]))
        T[r, i] = S[c, i]
    return T
def bulger(T, x):
    """index of the tracked band that bulges most in the middle of Gamma->X (None if < 50 meV)"""
    mid = (x > 0.40) & (x < 0.70)
    e3 = np.interp(0.30, x, T[0]) if False else None
    best, bi = 0.05, None
    for b in range(T.shape[0]):
        ends = max(np.interp(0.30, x, T[b]), np.interp(0.80, x, T[b]))
        bulge = T[b, mid].max() - ends
        if bulge > best: best, bi = bulge, b
    return bi

## Sigma q-mesh points on Gamma->X: (0,0,x) = (n/N)(b1+b2) with (0,0,1) = (b1+b2)/2, so x = 2n/N.
## N=6: 0, 1/3, 2/3, 1.  N=9: 0, 2/9, 4/9, 6/9, 8/9 -- X itself is not on an odd mesh.
_N = int(os.environ.get('MESH', '6').split(',')[0])
def qmesh(N): return [2*n/N for n in range(N) if 2*n/N <= 1 + 1e-9]
QMESH = qmesh(_N); LO, HI = 32, 44
## MESHCOLS="6,6,9,9": a mesh per column, when chains with different meshes sit side by side
_MC = [int(v) for v in os.environ['MESHCOLS'].split(',')] if 'MESHCOLS' in os.environ else None
def path(d, it): return f'{d}/bnd_lda.dat' if it == 0 else f'{d}/bnd_iter{it}.dat'

# rows are driven by the MLO chain (last column); the MTO column is shown alongside
iters = [0] + [i for i in range(1, 101) if any(os.path.exists(path(c[1], i)) for c in (COLS if NOREF else COLS[1:]))]
## ROWS="0,1,3,5,10": keep only these iterations (0 = LDA), for a figure that has to stay readable when small
if 'ROWS' in os.environ: iters = [i for i in iters if i in {int(v) for v in os.environ['ROWS'].split(',')}]
n = len(iters)
## A column whose directory name ends in _eg (e.g. a link to the t2g column's directory) shows the eg bands b45-52 (grey),
## b53 (purple; the dispersive band above them, 4-7 eV, which keeps rising at the 9^3 mesh point 4/9) and b54-58 (light grey)
## on 2.6-EGMAX eV, so that the t2g and the eg of a run can stand side by side (2026-09-28, user: 12 columns)
def is_eg(d): return os.path.basename(d.rstrip('/')).endswith('_eg')
YL_EG = (2.6, float(os.environ.get('EGMAX', '8.0')))   # 2026-09-28: up to 8 eV (user: a bit higher), bands b45-58
ROWH = float(os.environ.get('ROWH', '2.5'))   # height of a row (inch); ROWH=5 for a figure of one or two rows (2026-09-28)
MINW = float(os.environ.get('MINW', '9'))     # the title lines need about 9 inch; one column was clipped (2026-09-28)
COLW = float(os.environ.get('COLW', '4.9'))   # width of a column (inch); COLW=3.6 for the 12-column figure
## HEADS='big line|small line;...': one header per calculation (GROUP columns each) above the panels, in large type
## (';' and '|' separate the headers and their lines: the texts must not contain them)
## (2026-09-28, user: the settings of the 7 calculations written large on top)
HEADS = [h.split('|') for h in os.environ['HEADS'].split(';')] if os.environ.get('HEADS') else []
HEADH = 0.9 if HEADS else 0.0
FIGH = ROWH * n + 1.5 + HEADH   # 1.5 inch on top for the title of up to 5 lines (+ the headers), whatever the number of rows
fig, AX = plt.subplots(n, len(COLS), figsize=(max(COLW * len(COLS) + 0.4, MINW), FIGH),
                       sharex=True, sharey=('col' if any(is_eg(c[1]) for c in COLS) else True), squeeze=False)
store = {}
## OCC=1: only the two lowest t2g bands (the occupied pair near Gamma) on 0 < x < 0.5, with the dip of the lowest band below its
## Gamma value (x < 0.3) and their roughness (mean |E(k-1) - 2E(k) + E(k+1)| on x < 0.5, sorted per k) in each panel (2026-09-28)
OCC = os.environ.get('OCC', '0') == '1'
XMAX, YL = (0.5, (-0.58, -0.18)) if OCC else (1.0, (-0.9, 1.5))
DATA = {}   # the numbers of every panel, written to DATAOUT (2026-09-28, user: the numbers behind a plot must be kept)
## GROUP=2: two columns belong to one calculation (the 14-column figure); every other calculation gets a light background
## and a black line separates the calculations (2026-09-28, user: easier to see which panels go together)
GROUP = int(os.environ.get('GROUP', '1'))
TINT = '#e9ecf4'
LINE = dict(lw=0.6, color='k')      # the bands that are not highlighted: thin black (2026-09-28, user: not grey)
for r, it in enumerate(iters):
    for c, (lab, d, col) in enumerate(COLS):
        ax = AX[r][c]; f = path(d, it)
        if GROUP > 1 and (c // GROUP) % 2 == 1:
            ax.set_facecolor(TINT)
        if not os.path.exists(f):
            ## EMPTY=frame: an empty frame with the same axes, so that the columns stay aligned when runs
            ## have different iterations (2026-09-28, user: the figure looked shifted); default: no frame
            if os.environ.get('EMPTY', 'off') == 'frame':
                ax.set_ylim(*(YL_EG if is_eg(d) else YL)); ax.set_xlim(0, 1.0 if is_eg(d) else XMAX)
                ax.set_title(f'{lab}   iter {it}', fontsize=10, color='0.6')
                ax.text(0.5, 0.5, 'no band (not run yet / not kept)', transform=ax.transAxes, ha='center', va='center',
                        color='0.6', fontsize=9)
                ax.yaxis.set_major_locator(MultipleLocator(0.1 if OCC else 0.5))   # the y axis is shared: same ticks
            else:
                ax.axis('off')
            continue
        EG = is_eg(d)
        x, E = rd(f); S = np.sort(E, axis=0)[44:58] if EG else np.sort(E[LO:HI], axis=0); store[(c, it)] = S
        DATA[f'{os.path.basename(d)}_{"lda" if it == 0 else f"iter{it}"}'] = (x, S, os.path.realpath(f))
        if EG:
            for b in S[9:]: ax.plot(x, b, '-', lw=0.5, color='k')            # b54-58
            for b in S[:8]: ax.plot(x, b, '-', **LINE)                       # eg b45-52
            ax.plot(x, S[8], '-', lw=1.6, color='tab:purple')                 # b53
        elif OCC:
            ax.plot(x, S[1], '-', lw=0.9, color='k')
            ax.plot(x, S[0], '-', lw=1.6, color='tab:green')
            h = x < 0.5
            dip = (S[0][0] - S[0][x < 0.3].min()) * 1e3
            rough = np.mean([np.abs(b[h][:-2] - 2 * b[h][1:-1] + b[h][2:]).mean() for b in S[:2]]) * 1e3
            ax.text(0.02, 0.96, f'dip {dip:.1f} meV,  roughness {rough:.3f} meV', transform=ax.transAxes, va='top',
                    fontsize=9, backgroundcolor='white')
        elif HILITE and it > 0:
            T = track(S); bb = bulger(T, x)
            for k, b in enumerate(T):
                ax.plot(x, b, '-', **LINE)
            ax.plot(x, S[0], '-', lw=2.0, color='tab:green')                 # lowest band (energy-sorted)
            if bb is not None: ax.plot(x, T[bb], '-', lw=2.0, color='tab:blue')  # the band bulging mid Gamma-X
        else:
            for b in S: ax.plot(x, b, '-', **(LINE if it == 0 else dict(lw=1.2, color=col)))
        ## <band file>.label (one line) marks a panel whose data is not from the column's own run, e.g. an older code
        ## standing in until the run gets there (2026-09-28); the title gets it in orange
        tag = open(f + '.label').read().strip() if os.path.exists(f + '.label') else ''
        ## <band file>.stamp: when the band was made (the date of the source file), added in parentheses (2026-09-28)
        stamp = open(f + '.stamp').read().strip() if os.path.exists(f + '.stamp') else ''
        ax.set_title(f'{lab}   {"LDA" if it == 0 else f"iter {it}"}' + (f'  [{tag}]' if tag else '')
                     + (f'  ({stamp})' if stamp else ''), fontsize=float(os.environ.get('TITLEFS', '9.5')),
                     color='tab:orange' if tag else 'k')
        ax.set_ylim(*(YL_EG if EG else YL)); ax.set_xlim(0, 1.0 if EG else XMAX)
        if not EG: ax.axhline(0, color='k', lw=0.9)
        ## Sigma q-mesh points: a small red cross on every band (where the interpolation is exact).
        ## On 9^3 the mesh x = 2n/9 is not on the 211-point grid, so the band is interpolated there.
        _qm = [q for q in (qmesh(_MC[c]) if _MC else QMESH) if q <= (1.0 if EG else XMAX) + 1e-9]
        for b in (S[:2] if OCC and not EG else S):
            ax.plot(_qm, np.interp(_qm, x, b), 'x', color='red', ms=4.5, mew=1.1, zorder=5)
        if OCC and not EG: ax.yaxis.set_major_locator(MultipleLocator(0.1)); ax.yaxis.set_minor_locator(MultipleLocator(0.02))
        else:   ax.yaxis.set_major_locator(MultipleLocator(0.5)); ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.grid(axis='y', which='major', color='0.78', lw=0.6)
        ax.grid(axis='y', which='minor', color='0.92', lw=0.4)
    AX[r][0].set_ylabel('$E-E_F$ [eV]')
for a in AX[-1]: a.set_xlabel('$\\Gamma \\to X$')
stamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
MESH = os.environ.get('MESH', '6')
_ML = ' / '.join(f'{m}$^3$' for m in dict.fromkeys(MESH.split(',')))
fig.suptitle(f'LiTi$_2$O$_4$  {_ML} (nkabc = n1n2n3 = mesh)   '
             + ('the occupied t$_{{2g}}$ pair (the two lowest of b33-44) on the first half of $\\Gamma\\to X$\n' if OCC else
                f't$_{{2g}}$ (b33-44) along $\\Gamma\\to X$, 211 points\n')
             + 'red x = $\\Sigma$ q-mesh points (interpolation is exact there)\n'
             + (os.environ.get('MLOINFO', '') + '\n' if os.environ.get('MLOINFO') else '')
             + ('green = the lowest band, black = the second (energy-sorted at each k)\n' if OCC else
                ('green = lowest band, blue = the band that bulges mid Gamma-X (tracked through crossings), black = the rest\n' if HILITE else ''))
             + ('' if NOREF else
                'CAUTION: different APW cutoffs - MTO chain pwmode=1 (|G|), MLO chains pwmode=11 (|q+G|).\n'
                'Each column states its own damping; the MTO chain always Anderson-mixes sigm at [gw] mixbeta\n') +
             
             f'generated {stamp}', fontsize=10.5, y=1 - 0.1 / FIGH, va='top')
plt.tight_layout(rect=[0, 0, 1, 1 - (1.4 + HEADH) / FIGH])   # the title sits in the top 1.4 inch, the headers below it
if HEADS:
    G = max(GROUP, 1)
    top = max(a.get_position().y1 for a in AX[0])
    for g, h in enumerate(HEADS[:len(COLS) // G]):
        xc = 0.5 * (AX[0][g * G].get_position().x0 + AX[0][g * G + G - 1].get_position().x1)
        fig.text(xc, top + 0.62 / FIGH, h[0], ha='center', va='center', fontsize=22, fontweight='bold')
        if len(h) > 1:
            fig.text(xc, top + 0.3 / FIGH, h[1], ha='center', va='center', fontsize=11, color='0.25')
if GROUP > 1:                                        # a black line between the calculations, over the full height
    from matplotlib.lines import Line2D
    top = max(a.get_position().y1 for a in AX[0]); bot = min(a.get_position().y0 for a in AX[-1])
    for g in range(1, len(COLS) // GROUP):
        xs = 0.5 * (AX[0][g * GROUP - 1].get_position().x1 + AX[0][g * GROUP].get_position().x0)
        fig.add_artist(Line2D([xs, xs], [bot - 0.3 / FIGH, top + 0.25 / FIGH], transform=fig.transFigure, color='k', lw=1.2))
plt.savefig(OUT, dpi=float(os.environ.get('DPI', '115')))
## DATAOUT=<file>.npz: per panel <column dir>_<lda|iterN>: x (Gamma-X = 0..1) and the t2g energies (b33-44 sorted per k, eV from
## E_F) as plotted, with <key>__src (the band file: its .src if given, e.g. where it lies on the machine that ran it),
## __label and __stamp (from the .label/.stamp files)
if os.environ.get('DATAOUT'):
    arrs = {}
    for k, (x, S, src) in DATA.items():
        arrs[k + '__x'] = x; arrs[k + '__E'] = S; arrs[k + '__src'] = np.array(src)
        base = [c[1] for c in COLS if os.path.basename(c[1]) == k.rsplit('_', 1)[0]][0]   # the .label/.stamp/.src sit
        it = k.rsplit('_', 1)[1]; link = f'{base}/bnd_lda.dat' if it == 'lda' else f'{base}/bnd_{it}.dat'   # by the panel link
        for ext in ('label', 'stamp', 'src'):
            if os.path.exists(link + '.' + ext): arrs[k + '__' + ext] = np.array(open(link + '.' + ext).read().strip())
    np.savez_compressed(os.environ['DATAOUT'], **arrs)
    print('wrote', os.environ['DATAOUT'], f'({len(DATA)} panels)')
print('wrote', OUT, f'({n} rows x {len(COLS)} cols)')
