#!/usr/bin/env python3
import pandas as pd
import pymatgen.core as mg
from pymatgen.ext.matproj import MPRester
import sys,os
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)
sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', buffering=1)
sys.stdin  = os.fdopen(sys.stdin.fileno(), 'r', buffering=1)

apikey = "FRaVT8iIl5JXLimgwLg"
maxsite=4
ffile='list.job_mp'
f=open(ffile,'w')

RG='He Ne Ar Xe Kr Xe Rn '
LN='Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu ' #La can be contained.
AC='Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr '  #Ac can be contained.
noelements=(RG+LN+AC).split()
print(noelements)
with MPRester(apikey) as mpr:
    c={"nsites":{"$lte":maxsite},'elements':{'$nin':noelements}}
    p=['task_id','formula',"band_gap","elasticity","nsites"]
    material = mpr.query(criteria=c,properties=p)
    print('number of files query returned=',len(material))
    count=0
    for m in material:
        if(m['task_id'][0:3]!='mp-'): continue #I found mvc-15384. is it mp file?
        
        try:
            bm = m["elasticity"]['K_Voigt_Reuss_Hill']
        except:
            bm = None
        bg =m['band_gap']
        nn =m['nsites']
        count = count+1
        print(m['task_id'].strip('mp-'),nn,m['formula'],bg,bm,file=f)
f.close()
print('### number of total mp file (no mvc file=)',count)
sys.exit()
