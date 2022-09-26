#!/usr/bin/env python3
import pandas as pd
import pymatgen.core as mg
from pymatgen.ext.matproj import MPRester
import sys,os
def getlist(avoided,required,minsite,maxsite,ffile):
    RG='He Ne Ar Xe Kr Xe Rn '
    LN='La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu ' 
    AC='Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg' 
    nn=''
    if 'RG' in avoided: nn=nn+RG
    if 'LN' in avoided: nn=nn+LN
    if 'AC' in avoided: nn=nn+AC
    noelements=nn.split()
    dicelements={'$nin':noelements}
    if('ANY' not in required):
        nn=''
        if 'RG' in required: nn=nn+RG
        if 'LN' in required: nn=nn+LN
        if 'AC' in required: nn=nn+AC
        dicelements['$in']=nn.split()
    apikey = "FRaVT8iIl5JXLimgwLg"
    with MPRester(apikey) as mpr:
        c={"nsites":{"$lte":maxsite,"$gte":minsite},'elements':dicelements}
        p=['task_id','formula',"band_gap","elasticity","nsites"]
        material = mpr.query(criteria=c,properties=p)
        #print('number of files query returned=',len(material))
        
        nxx= 500 #divided to ~nxx jobs per file 
        GetPOSCAR=False
        ndiv=(len(material)//nxx+1) 
        nsize=len(material)//ndiv+1 
        print(ndiv,nsize)
        count=0
        ff=open(ffile,'w')
        fnum=0
        for m in material:
            if count % nsize==0:
                fnum=fnum+1
                ff.close()
                ff=open(ffile+'.'+str(fnum),'w')
                print(count)
            if(m['task_id'][0:3]!='mp-'): continue #I found mvc-15384. is it mp file?
            try:
                bm = m["elasticity"]['K_Voigt_Reuss_Hill']
            except:
                bm = None
            bg =m['band_gap']
            nn =m['nsites']
            count = count+1
            print(m['task_id'].strip('mp-'),nn,m['formula'],bg,bm,file=ff)
            #print(m['task_id'])
            if(GetPOSCAR):
                struc = mpr.get_structure_by_material_id(m['task_id'])
                struc.to(fmt='poscar', filename='POSCARALL/POSCAR.'+m['task_id'])
        ff.close()
        print('### file    =',ffile)
        print('  maxsite =',maxsite)
        print('  required=',required)
        print('  avoided =',avoided)
        print('  numbder of mp (no mvc file=)',count)
    return count
#-------------------------------------
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)
sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', buffering=1)
sys.stdin  = os.fdopen(sys.stdin.fileno(), 'r', buffering=1)

maxsite=8
minsite=1
avoided='RG AC LN'
required='ANY'
ffile='list.job_mp.no4f'
count=getlist(avoided,required,minsite,maxsite,ffile)

avoided='RG AC'
required='LN'
ffile='list.job_mp.4f'
count=getlist(avoided,required,minsite,maxsite,ffile)

sys.exit()
