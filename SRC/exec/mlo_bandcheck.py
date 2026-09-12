#!/usr/bin/env python3
"""v6.2 統合評価: [occupied .. CBM+3eV] 窓 rms / ギャップ / 有効質量(曲率比)"""
import glob, os, re, math, bisect
import numpy as np

RY = 13.605
BASE = os.path.expanduser('~/ecalj/Samples/MLOsamples')

def read_ef(g):
    for l in open(g):
        m = re.search(r'ef\s*=\s*([-+0-9.eEdD]+)', l)
        if m: return float(m.group(1).replace('D','E'))

def load_mlo(f, ef):
    P=[]
    for l in open(f):
        t=l.split()
        if len(t)<2 or l.lstrip().startswith('#'): continue
        try: P.append((float(t[0]),(float(t[1])-ef)*RY))
        except: pass
    return P

def load_bnd(fs):
    P=[]
    for f in fs:
        for l in open(f):
            t=l.split()
            if len(t)<3 or l.lstrip().startswith('#'): continue
            try: P.append((float(t[1]),float(t[2])))
            except: pass
    return P

def edges(P, thr=0.25):
    occ=[e for _,e in P if e<thr]; uno=[e for _,e in P if e>=thr]
    return (max(occ) if occ else None, min(uno) if uno else None)

def window_rms(M, B, elo, ehi):
    bx={}
    for x,e in B: bx.setdefault(round(x,5),[]).append(e)
    xs=sorted(bx); dev=[]
    for x,e in M:
        if not (elo<=e<=ehi): continue
        i=bisect.bisect_left(xs,x); best=None
        for k in (i-1,i,i+1):
            if 0<=k<len(xs) and (best is None or abs(xs[k]-x)<abs(best-x)): best=xs[k]
        if best is None or abs(best-x)>0.02: continue
        dev.append(min(abs(e-b) for b in bx[best]))
    dev.sort()
    if not dev: return None,None,0
    r=math.sqrt(sum(d*d for d in dev)/len(dev))
    p95=dev[int(0.95*len(dev))-1] if len(dev)>1 else dev[-1]
    return r,p95,len(dev)

def edge_curv(P, edge, thr=0.25, dx=0.15):
    if edge=='vbm':
        cand=[(x,e) for x,e in P if e<thr]; x0,e0=max(cand,key=lambda p:p[1])
        sel={}
        for x,e in P:
            if abs(x-x0)<=dx and e<thr: sel[round(x,5)]=max(sel.get(round(x,5),-1e9),e)
    else:
        cand=[(x,e) for x,e in P if e>=thr]; x0,e0=min(cand,key=lambda p:p[1])
        sel={}
        for x,e in P:
            if abs(x-x0)<=dx and e>=thr: sel[round(x,5)]=min(sel.get(round(x,5),1e9),e)
    xs=np.array(sorted(sel)); es=np.array([sel[x] for x in xs])
    if len(xs)<5: return None
    return np.polyfit(xs-np.float64(x0),es-e0,2)[0]

INSUL={'C','GaAs','SrTiO3','NiO666lda','Al2O3_Cr','Si666gwsc','C.sp'}
ALL=['Fe','Cu','C','C.sp','GaAs','SrTiO3','NiO666lda','Al2O3_Cr','RuO2','SmP','FeCo','FeMgO','Si666gwsc','GdCo5','GdION']
print(f"{'system':12s} | {'rms[occ,CBM+3]':>14s} {'p95':>7s} {'n':>5s} | {'DFTgap':>7s} {'v6.2gap':>8s} {'dVBM':>7s} {'dCBM':>7s} | {'m*VBM':>6s} {'m*CBM':>6s}")
for s in ALL:
    w=os.path.join(BASE,s+'__m3_work')
    try:
        ef=read_ef(os.path.join(w,'bandplot_MLO.isp1.glt'))
        M=load_mlo(os.path.join(w,'band_MLO_spin1.dat'),ef)
        B=load_bnd(sorted(glob.glob(os.path.join(w,'bnd0*.spin1'))))
        assert M and B
    except Exception as e:
        print(f"{s:12s} | (missing)"); continue
    vd,cd_=edges(B); vm,cm_=edges(M)
    cbm_ref = cd_ if (s in INSUL and cd_ is not None) else 0.0
    r,p95,n = window_rms(M,B,-25.0, cbm_ref+3.0)
    gapstr=f"{'—':>7s} {'—':>8s} {'—':>7s} {'—':>7s}"
    if s in INSUL and None not in (vd,cd_,vm,cm_):
        gapstr=f"{cd_-vd:7.3f} {cm_-vm:8.3f} {vm-vd:+7.3f} {cm_-cd_:+7.3f}"
    ms=[]
    for edge in ('vbm','cbm'):
        try:
            a,b=edge_curv(B,edge),edge_curv(M,edge)
            ms.append(f"{a/b:6.2f}" if (a and b and abs(b)>1e-9) else f"{'—':>6s}")
        except Exception:
            ms.append(f"{'—':>6s}")
    print(f"{s:12s} | {r if r else float('nan'):14.3f} {p95 if p95 else float('nan'):7.3f} {n:5d} | {gapstr} | {ms[0]} {ms[1]}")
