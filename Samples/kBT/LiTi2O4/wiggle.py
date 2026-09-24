import numpy as np, os
def rd(f):
    xs,es,cx,ce=[],[],[],[]
    for ln in open(f):
        if ln.startswith('#'): continue
        t=ln.split()
        if len(t)<3:
            if cx: xs.append(cx); es.append(ce); cx=[];ce=[]
            continue
        cx.append(float(t[1])); ce.append(float(t[2]))
    if cx: xs.append(cx); es.append(ce)
    return np.array(xs[0]), np.array(es)
IDX=[0,70,140,210]; LO,HI=32,44
def chord(y,x):
    o=[]
    for a,b in zip(IDX[:-1],IDX[1:]):
        c=np.interp(x[a:b+1],[x[a],x[b]],[y[a],y[b]]); o.append(np.abs(y[a:b+1]-c).max())
    return np.array(o)*1e3
def stats(f):
    x,E=rd(f); ch=[];tp=[];res=[]
    for ib in range(LO,HI):
        y=E[ib]; ch.append(chord(y,x).mean())
        d=np.gradient(y,x)
        # turning points: sign changes of slope, ignoring numerical noise below 0.05 eV/x
        s=np.sign(d); s[np.abs(d)<0.05]=0
        s=s[s!=0]; tp.append(int(np.sum(s[1:]!=s[:-1])))
        # residual from a smooth cubic fit over the whole line
        c=np.polyfit(x,y,3); res.append(np.abs(y-np.polyval(c,x)).max()*1e3)
    return np.mean(ch), np.mean(tp), np.max(tp), np.mean(res), np.max(res)
cases=[('LDA','liti_mlo/bnd_lda.dat'),('MTO iter1','liti_ref/bnd_iter1.dat'),
       ('MLO76 iter1','liti_mlo/bnd_iter1.dat'),('MLO154 iter1','liti_mlo154/bnd_iter1.dat'),
       ('MLO76 iter4','liti_mlo/bnd_iter4.dat'),('MLO76 iter10','liti_mlo/bnd_iter10.dat')]
print(f"{'case':14s} {'chord':>7s} {'turn/band':>10s} {'turnmax':>8s} {'cubicres mean':>14s} {'max':>8s}")
for l,f in cases:
    if os.path.exists(f):
        a,b,c,d,e=stats(f); print(f"{l:14s} {a:7.1f} {b:10.2f} {c:8d} {d:14.1f} {e:8.1f}")
