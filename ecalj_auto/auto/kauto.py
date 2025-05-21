import numpy as np
import re, os,sys,math
def readPlatQlat():
    qvec, length = [], []
    Plat=[None]*3
    Qlat=[None]*3
    with open("PlatQlat.chk","r") as f:
        txt = f.readlines()
        for i,t in enumerate(txt[0:3]):
            Plat[i]= [float(x) for x in t.split()[0:3]]
            Qlat[i]= [float(x) for x in t.split()[3+0:3+3]]
        Alat=float(txt[4].split()[0])
        Blat=float(txt[5].split()[0])
        Plat=np.array(Plat)
        Qlat=np.array(Qlat) 
        print('Plat=',Plat)
        print('Qlat=',Qlat)
        print('Alat Blat=',Alat,Blat) # bohr
    return Plat,Qlat,Alat,Blat
def kauto(n): #https://www.vasp.at/wiki/index.php/KPOINTS
    # n=4 returns 4 4 4 for Si
    Plat,Qlat,Alat,Blat=readPlatQlat()
    absQlat= [np.sum(x**2)**.5*Blat for x in Qlat] # bohr-1
    R_k=  6.7* .529177*n/4  
    nnn= [int(max(np.floor(R_k*x+.5),1)) for x in absQlat]
    #nnn= [max(        (R_k*x+.5),1) for x in absQlat]
    return(nnn)
if __name__=='__main__':
    n=8
    nnn=kauto(n)
    print(nnn)
    sys.exit()
