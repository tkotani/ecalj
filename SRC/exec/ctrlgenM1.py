#!/usr/bin/env python3
# CAUTION; R here is in Angstrom, converted to a.u. in ctrlgenM1.ctrl.* 
#########################################################################
# Generate a temprate of ctrl file from ctrls.
# Work with python 2.5 or so.
#  Takao Kotani and Hiori Kino
# ---  a ctrls is ---
# STRUC   ALAT=7.37 
#         PLAT=1 0 0  0 1 0  0 0 1
# SITE
#   ATOM=Sr POS=1/2 1/2 1/2
#   ATOM=Ti POS= 0   0   0+0
#   ATOM=O  POS=1/2  0   0
#   ATOM=O  POS= 0  1/2  0
#   ATOM=O  POS= 0   0  1/2
# ---- end here ------
# (HEADER and so are at the begining of lines)
# take Kino's change.
#########################################################################
import os, sys, string, re, locale

atomlist="""

### !!! CAUTION!!! R= in this table is in Angstrome. !!!!!

### eh.3 eh2.2 works OK?
H=   " atomz=1@  pz=''@       p=''@          eh=-1*3@ eh2=-2*2@    R=0.38@" 
He=  " atomz=2@  pz='PZ=2.5'@ p=''@          eh=-1*3@ eh2=-2*2@    R=1.17@"  
Li=  " atomz=3@  pz='PZ=1.9'@ p=''@          eh=-1*3@ eh2=-2*2@    R=1.37@"        
Be=  " atomz=4@  pz='PZ=1.9'@ p=''@          eh=-1*3@ eh2=-2*2@    R=1.22@"        
B=   " atomz=5@  pz='PZ=1.9'@ p=''@          eh=-1*3@ eh2=-2*2@    R=0.77@"      
C=   " atomz=6@  pz=''@       p=''@          eh=-1*3@ eh2=-2*2@    R=0.66@"       
N=   " atomz=7@  pz=''@       p=''@          eh=-1*3@ eh2=-2*2@    R=0.55@"      
O=   " atomz=8@  pz=''@       p=''@          eh=-1*3@ eh2=-2*2@    R=0.61@"      
F=   " atomz=9@  pz=''@       p=''@          eh=-1*3@ eh2=-2*2@    R=0.71@"      
Ne=  " atomz=10@ pz='PZ=3.5,3.5'@  p=''@     eh=-1*3@ eh2=-2*2@    R=1.28@"       

Na=  " atomz=11@ pz='PZ=0,2.8'@  p=''@       eh=-1*4@ eh2=-2*3@    R=1.55@"   
Mg=  " atomz=12@ pz=''@          p=''@       eh=-1*4@ eh2=-2*3@    R=1.75@"                 
Al=  " atomz=13@ pz=''@          p=''@       eh=-1*4@ eh2=-2*3@    R=1.25@"
#Na=  " atomz=11@ pz='PZ=2.8,2.8'@ p=''@     eh=-1*4@ eh2=-2*3@    R=1.55@"                 
#Mg=  " atomz=12@ pz='PZ=2.8,2.8'@ p=''@     eh=-1*4@ eh2=-2*3@    R=1.75@"                 
#Al=  " atomz=13@ pz='PZ=2.9,2.9'@ p=''@     eh=-1*4@ eh2=-2*3@    R=1.25@"
Si=  " atomz=14@ pz=''@            p=''@     eh=-1*4@ eh2=-2*3@    R=1.15@"
P=   " atomz=15@ pz=''@            p=''@     eh=-1*4@ eh2=-2*3@    R=0.96@"       
S=   " atomz=16@ pz=''@            p=''@     eh=-1*4@ eh2=-2*3@    R=0.97@"       
Cl=  " atomz=17@ pz=''@            p=''@     eh=-1*4@ eh2=-2*3@    R=1.03@"       
Ar=  " atomz=18@ pz='PZ=4.3,4.3'@  p=''@     eh=-1*4@ eh2=-2*3@    R=1.72@"       

#################################################################################
##### For Fe2 dimer, using 3.9 (like PZ=0,0,3.9) is better for stable convergence
##### than PZ=13.9, which means "automatic determination of EH for 3rd MTO"
##### can introduce some unsystematic behevior.
#################################################################################

#K=   " atomz=19@ pz='PZ=3.9,3.9'@ p=''@ eh=-1*4@ eh2=-2*3@             R=2.00@"       
#Ca=  " atomz=20@ pz='PZ=3.9,3.9'@ p=''@ eh=-1*4@ eh2=-2*3@             R=2.08@"       
#Sc=  " atomz=21@ pz='PZ=3.9,3.9'@ p=''@ eh=-1*4@ eh2=-2*3@             R=1.31@"       
### 3p is deep, but it may affect to band gaps and so on. ###
#Fe= "  atomz=26@ pz='PZ=0,3.9,4.5'@ p=''@ eh=-1*4@ eh2=-2*3@ R=1.00@"
#Co= "  atomz=27@ pz='PZ=0,3.9,4.5'@ p=''@ eh=-1*4@ eh2=-2*3@ R=1.00@"
#Ni= "  atomz=28@ pz='PZ=0,3.9,4.5'@ p=''@ eh=-1*4@ eh2=-2*3@ R=1.06@"

K=   " atomz=19@ pz='PZ=0,3.9'@   p=''@ eh=-1*4@ eh2=-2*3@          R=2.00@"       
Ca=  " atomz=20@ pz='PZ=0,3.9'@   p=''@ eh=-1*4@ eh2=-2*3@          R=2.08@"       
Sc=  " atomz=21@ pz='PZ=0,3.9'@   p=''@ eh=-1*4@ eh2=-2*3@          R=1.27@"       
Ti=  " atomz=22@ pz='PZ=0,3.9'@   p=''@ eh=-1*4@ eh2=-2*3@          R=0.95@"       
V=   " atomz=23@ pz='PZ=0,3.9'@   p=''@ eh=-1*4@ eh2=-2*3@          R=0.87@"       
Cr=  " atomz=24@ pz='PZ=0,3.9'@   p=''@ eh=-1*4@ eh2=-2*3@          R=0.80@"       
Mn=  " atomz=25@ pz='PZ=0,3.9'@   p=''@ eh=-1*4@ eh2=-2*3@          R=0.82@"       
Fe= "  atomz=26@ pz='PZ=0,3.9'@   p=''@ eh=-1*4@ eh2=-2*3@   R=1.00@"
Co= "  atomz=27@ pz='#PZ=0,3.9'@   p=''@ eh=-1*4@ eh2=-2*3@   R=1.00@"
Ni= "  atomz=28@ pz='#PZ=0,3.9'@   p=''@ eh=-1*4@ eh2=-2*3@   R=1.06@"
Cu= "  atomz=29@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@   R=1.13@"
Zn= "  atomz=30@ pz='PZ=0,0,3.9'@ p='P=0,0,4.5'@ eh=-1*4@ eh2=-2*3@   R=1.60@"
Ga= "  atomz=31@ pz='PZ=0,0,3.9'@ p='P=0,0,4.5'@ eh=-1*4@ eh2=-2*3@   R=1.37@"
Ge= "  atomz=32@ pz='PZ=0,0,3.9'@ p='P=0,0,4.5'@ eh=-1*4@ eh2=-2*3@   R=1.21@"
As= "  atomz=33@ pz='#PZ=0,0,3.9'@ p='#P=0,0,4.5'@ eh=-1*4@ eh2=-2*3@   R=1.06@"
Se= "  atomz=34@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@                      R=1.10@"
Br= "  atomz=35@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@                      R=1.16@"
Kr= "  atomz=36@ pz='PZ=5.3,5.3'@  eh=-1*4@ eh2=-2*3@                 R=1.88@"

## When we use large R, it may cause basis-independency error due to s-orbital
### Followings are not tested well ###
#Sr  =" atomz=38@ pz='PZ=4.9,4.9'@ p=''@ eh=-1*4@ eh2=-2*3@   R=2.08@"       
#Y   =" atomz=39@ pz='PZ=4.9,4.9'@ p=''@ eh=-1*4@ eh2=-2*3@   R=1.31@" (same R with Sc) OK???
Rb  =" atomz=37@ pz='PZ=0,4.9'@ p=''@ eh=-1*4@ eh2=-2*3@   R=2.00@"
### PZ=4.9,4.9 is for GW
Sr  =" atomz=38@ pz='PZ=4.9,4.9'@ p=''@ eh=-1*4@ eh2=-2*3@   R=2.0@"       
### PZ=4.9,4.9 is for GW
Y   =" atomz=39@ pz=''@ p='PZ=4.9,4.9'@ eh=-1*4@ eh2=-2*3@   R=1.3@" (same R with Sc) OK???
Zr  =" atomz=40@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@   R=?@"
Nb  =" atomz=41@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@   R=?@" 
Mo  =" atomz=42@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@   R=?@" 
Tc  =" atomz=43@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@   R=?@" 
Ru  =" atomz=44@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Rh  =" atomz=45@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Pd  =" atomz=46@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Ag  =" atomz=47@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Cd  =" atomz=48@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
In= "  atomz=49@ pz='PZ=0,0,4.9'@ p='P=0,0,5.2'@ eh=-1*4@ eh2=-2*3@ R=1.37@"
Sn  =" atomz=50@ pz='PZ=0,0,4.9'@ p='P=0,0,5.2'@ eh=-1*4@ eh2=-2*3@ R=?@" 
Sb  =" atomz=51@ pz='PZ=0,0,4.9'@ p='P=0,0,5.2'@ eh=-1*4@ eh2=-2*3@ R=?@" 
#Sb  =" atomz=51@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Te  =" atomz=52@ pz='#PZ=0,0,4.9'@ p='#P=0,0,5.5'@ eh=-1*4@ eh2=-2*3@ R=?@" 
I   =" atomz=53@ pz='#PZ=0,0,4.9'@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Xe  =" atomz=54@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 

# We need to give P= for p channel for Cs excplicitly, because default is P=0,5 for p channel (due to historical reason...)
Cs  =" atomz=55@ pz='PZ=5.9,5.9'@ p='P=6.2,6.2'@ eh=-1*4@ eh2=-2*3@ R=?@" 
Ba  =" atomz=56@ pz='PZ=0,5.9'@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
#La  =" atomz=57@ pz='PZ=0,5.9'@ p=''@ eh=-1*4@ eh2=-2*3@  R=1.6@" ( R is just given for test) 
La  =" atomz=57@ pz='PZ=0,5.9'@ p=''@ eh=-1*4@ eh2=-2*3@  R=1.6@" ( R is just given for test) 
Ce  =" atomz=58@ pz='PZ=0,5.9'@ p=''@ eh=-1*4@ eh2=-2*4@ R=2.80@" 
Pr  =" atomz=59@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Nd  =" atomz=60@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Pm  =" atomz=61@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Sm  =" atomz=62@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Eu  =" atomz=63@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Gd  =" atomz=64@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Tb  =" atomz=65@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Dy  =" atomz=66@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Ho  =" atomz=67@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Er  =" atomz=68@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Tm  =" atomz=69@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Yb  =" atomz=70@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Lu  =" atomz=71@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
#Hf  =" atomz=72@  pz='PZ=0,0,0,4.9'@ p=''@ eh=-1*5@ eh2=-2*4@  R=?@" 
Hf  =" atomz=72@  pz='PZ=0,5.5,0,0'@ p='P=0,6.5,0,4.5 Q=2,6,2,14'@ eh=-1*5@ eh2=-2*4@
#PZ=0,0,0,4.9'@ p=''@ eh=-1*5@ eh2=-2*4@  R=?@" 
Ta  =" atomz=73@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
W   =" atomz=74@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Re  =" atomz=75@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Os  =" atomz=76@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Ir  =" atomz=77@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Pt  =" atomz=78@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Au  =" atomz=79@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Hg  =" atomz=80@ pz='PZ=0,0,5.9'@ p='#P=0,0,6.9'@ eh=-1*4@ eh2=-2*3@ R=?@" 
Tl  =" atomz=81@ pz='PZ=0,0,5.9'@ p='#P=0,0,6.9'@ eh=-1*4@ eh2=-2*3@ R=?@" 
Pb  =" atomz=82@ pz='PZ=0,0,5.9'@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Bi  =" atomz=83@ pz='PZ=0,0,5.9'@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
#Bi  =" atomz=83@ pz='PZ=0,0,5.9'@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 

Po  =" atomz=84@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
At  =" atomz=85@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Rn  =" atomz=86@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Fr  =" atomz=87@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Ra  =" atomz=88@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Ac  =" atomz=89@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Th  =" atomz=90@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Pa  =" atomz=91@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
U   =" atomz=92@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Np  =" atomz=93@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Pu  =" atomz=94@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Am  =" atomz=95@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Cm  =" atomz=96@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Bk  =" atomz=97@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@" 
Cf  =" atomz=98@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@"
Es  =" atomz=99@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@"
Fm  =" atomz=100@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@"
Md  =" atomz=101@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@"
No  =" atomz=102@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@"
Lr  =" atomz=103@ pz=''@ p=''@ eh=-1*4@ eh2=-2*3@ R=?@"
### CAUTION!!! R= in this table is in Angstrome.
### This gives smaller basis than those used in our dimer paper. ###
"""


#---------------------------------------------------
def manip_argset(argset):
#    global nspin_val, xcfun_val, mmom_val,  systype_val,nk_val,showatomlist,showhelp #rlmchk,r_mul_val,
    error_title="ARGUMENT ERROR"
    ierror=0
# defalut setting
    nspin_val="1"
    so_val="0"
    xcfun_str="vwn"
    xcfun_val=103
    mmom_val="#MMOM=0,0,0,0"
    systype_val="bulk"
    nk_val1="4"
    nk_val2="-999999"
    nk_val3="-999999"
    mmom_val='#MMOM=0 0 1 0'
    showatomlist=0
    showhelp=0
    metali=3
    fsmom_val=0.0
    ssig_val=1.0
    touchingratio=.97 #default value was -1.0 
    eh1set=1
#    readrmt=0
#    rlmchk=0
#    r_mul_val="1.0"
#    systype_val="molecule"
    for arg in argset:
        if re.match("--nspin",arg)!=None :
            nspinlist=arg.split("=")
            if len(nspinlist)==2:
                nspin_val=nspinlist[1]
        elif re.match("--so",arg)!=None:
            nsolist=arg.split("=")
            if len(nsolist)==2:
                so_val=nsolist[1]
            nspin_val='2'
        elif re.match("--xcfun",arg)!=None:
            xclist=arg.split("=")
            if len(xclist)==2:
                xcfun_str=xclist[1]
        elif re.match("--mmom",arg)!=None:
            mmomlist=arg.split("=")
            if len(mmomlist)==2:
                mmom_val=mmomlist[1]
        elif re.match("--insulator",arg)!=None:
            metali=0
#    elif re.match("--r_mul",arg)!=None:
#                rlist=arg.split("=")
#                if len(rlist)==2:
#                        r_mul_val=rlist[1]
        elif re.match("--nk1",arg)!=None:
            nklist=arg.split("=")
            if len(nklist)==2:
                nk_val1=nklist[1]
        elif re.match("--nk2",arg)!=None:
            nklist=arg.split("=")
            if len(nklist)==2:
                nk_val2=nklist[1]
        elif re.match("--nk3",arg)!=None:
            nklist=arg.split("=")
            if len(nklist)==2:
                nk_val3=nklist[1]
        elif re.match("--systype",arg)!=None:
            syslist=arg.split("=")
            if len(syslist)==2:
                systype_val=syslist[1]
        elif re.match("--fsmom",arg)!=None:
            nklist=arg.split("=")
            if len(nklist)==2:
                fsmom_val=nklist[1]
        elif re.match("--ssig",arg)!=None:
            nklist=arg.split("=")
            if len(nklist)==2:
                ssig_val=nklist[1]
        elif re.match("--tratio",arg)!=None:
            ttt=arg.split("=")
            if len(ttt)==2:
                touchingratio=ttt[1]
                touchingratio=float(touchingratio)
        elif re.match("--ehmol",arg)!=None:
                eh1set=0
#    elif arg=="--rlmchk":
#        rlmchk=1
        elif arg=="--showatomlist":
            showatomlist=1
        elif arg=="--help":
            showhelp=1
#    elif arg=="--readrmt":
#        readrmt=1
        elif re.match("-",arg):
            sys.stderr.write( error_title + ", unknown arg:  "+arg+"\n")
            ierror+=1
        
    if xcfun_str.upper()=="PBE":
        xcfun_val="103"
    elif xcfun_str.upper()=="VWN":
        xcfun_val="1"
    elif xcfun_str.upper()=="BH":
        xcfun_val="2"
    else:
        sys.stderr.write(error_title+", --xc="+xcfun_str+" : unknown\n")
        ierror+=1

    if systype_val.upper()=="MOLECULE":
        do_nothing=0
        nk_val1="1"
        nk_val2="1"
        nk_val3="1"
    elif systype_val.upper()=="BULK":
        do_nothing=0
    else:
        sys.stderr.write(error_title+", --systype="+systype_val+" : unknown\n")
        ierror+=1

    if ierror!=0:
        print ("ABORT. Check names of args. Some are not defined.")
        sys.exit(-1)
    return nspin_val,so_val,xcfun_val,xcfun_str,mmom_val,systype_val,nk_val1,nk_val2,nk_val3, showatomlist, showhelp,metali,fsmom_val,ssig_val,touchingratio,eh1set


#-------------------------------------------------------
def  line2Token(linein):
    """ convert the result of readline to token """
#    """ input=readlines() output=token""
    listout = []
    for s in linein:
        if s=="":
            continue
        s = s.replace('\n','')
        s = s.replace(',',' ')
        s = s.replace('=',"= ")
        s = s.replace(':',": ")

        lista=s.split()
        for x in lista:
            if x!="":
                listout.append(x)

    return listout

#---------------------------------------------------
def lineReadfile(filename):
    """ read file and make list of the content of the file, and \n->'' """
#    "input:filename output=readlines() "
    f = open(filename)
    list1 =[]
    while 1:
        s = f.readline()
        if s=="":
            break
        s=s.replace("\n","")
        if s=="":
            continue
        list1.append(s)
    f.close()
    return list1

#----------------------------------------------------------------------------
def RemoveCat(listctrl,key):
    """ remove a category (key) from listctrl input:devided line, output: combined line"""
    res=''
    ix=0
    sss ='^('+key.upper()+'|'+key.lower()+')'+'(\s|\Z)'
    for x in listctrl:
        if(ix==1):
            if(re.match('^\w',x)): ix=0
        if(re.match(sss,x)): ix=1
        if(ix==0): res = res+'\n'+ x

    res=res+'\n'
#    print res
#    sys.exit()

    return res

def RemoveCat2(listctrl,key):
    """ remove a category (key) from listctrl input:devided line, output: devided line"""
    res=[]
    ix=0
    sss ='^('+key.upper()+'|'+key.lower()+')'+'(\s|\Z)'
    for x in listctrl:
        if(ix==1):
            if(re.match('^\w',x)): ix=0
        if(re.match(sss,x)): ix=1
        if(ix==0): res.append(x)

#    res=res+'\n'
#    print res
#    sys.exit()
    return res

def GetCat(listctrl,key):
    """ get a category (key) from listctrl This returns lines. Not good correspondence to RemoveCat """
    res=[]
    ix=0
    sss ='^('+key.upper()+'|'+key.lower()+')'+'(\s|\Z)'
#    print sss
    for x in listctrl:
#        print x
        if(ix==1):
            if(re.match('^\w',x)): break
        if(re.match(sss, x )): ix=1
        if(ix==1): res.append(x)

#    res=res+'\n'
#    print res
#    sys.exit()
    return res

def countnum(mmm,key):
    xx2=re.split(key+"\s*",mmm)
    #print 'xx2=',key, xx2
    try:
        xx=re.split(' *',xx2[1])
    except:
        return 0
    num=0
    for i in xx:
        try:
            yy = float(i)
            #print yy
            num=num+1
        except:
            break
    return num

def getsitename(listsite):
    ddd=[]
    for x in listsite:
        xx=re.split('\WATOM=\W*',x)
        ddd=ddd+xx[1:]
    rrr=[]
    for i in ddd:
        rrr.append(re.split(' ',i)[0])
    return rrr

def glist(list):
    aaa=''
    for i in list:
        aaa=aaa+i+'\n'
    return aaa

def uniq(list):
    result = []
    for l in list:
        if not l in result:
            result.append(l)
    return result
def getdataa(aaa,key):
    return (aaa.split(key)[1].split('@')[0]).strip()+' '

def getdataa2(aaa,key):
    return (aaa.split(key+"'")[1].split("'@")[0]).strip()+' '



#===============================================================================
#          MAIN 
#===============================================================================

### readin args and set them in global variables ###
nargv = len(sys.argv) -1
argset= set(sys.argv[1:])
nspin_val, so_val,xcfun_val,xcfun_str, mmom_val, systype_val, nk_val1,nk_val2,nk_val3, showatomlist, showhelp,metali,fsmom_val,ssig_val,touchingratio,eh1set \
= manip_argset(argset) #set global argument defined at the top of mainp_argset.

if(so_val !='0' and nspin_val !='2'): sys.exit("Error: so is not zero but with nspin=1")
manip_argset(argset) #set global argument defined at the top of mainp_argset.
if nk_val2=="-999999":nk_val2=nk_val1
if nk_val3=="-999999":nk_val3=nk_val1

showhelpw='Not exist'
if(showhelp): showhelpw='given'
showatomlistw='Not exist'
if(showatomlist): showatomlistw='given'
insulatorw='Not exist'
if(metali==0): insulatorw='given'


### help sections ###
if(showhelp==1):
    print( \
""" 
ctrlgenM1.py. tkotani and h.kino aug_2013 version :
---------------
 Purpose: 
     Generate a template of ctrl file names as ctrlgenM1.ctrl.{ext}."
     Before you run lmfa, you have to copy ctrlgenM1.ctrl.{ext} to ctrl.{ext} and edit it.
 Usage  : ctrlgenM1 {extension of ctrl file} [option]
          [options] = INPUT arguments in the followings.
          Your given options (also defaults when not specified) are shown at the begining of console output.
 Example: 
       After you write ctrls.si, run
       >ctrlgenM1.py si --nk1=8 --nk2=8 --nk3=8 --tratio=1.0  --xcfun=vwn
""")
print() 
print( " === INPUT options (shown values are default) === ")
print( "  --help  %s"      % showhelpw)
print( "  --showatomlist  %s"      % showatomlistw)
print( "  --nspin=%s"  % nspin_val)
print( "  --so=%s"     % so_val)
print( "  --nk1=%s Division for BZ integral along a-axis"   % nk_val1)
print( "  --nk2=%s   (if not give, nk2=nk1) along b-axis"   % nk_val2)
print( "  --nk3=%s   (if not give, nk3=nk1) along c-axis"   % nk_val3)
print( "  --xcfun=%s   !(bh,vwn,pbe)"      % xcfun_str )#,xcfun_val
#print( "  --mmom='%s\' ! mmom is the initial magnetic moment for each spec,l-channel. Effective for --nspin=2" % mmom_val
print( " +++ Followings are for experts to change +++")
print( "  --tratio=%s (for MT radius: we use touching MT radius \\times this ratio. lmf --getwsr is called." % touchingratio)
print( "               if negative, we use use defalut MT radius in ctrlgenM1.py)")
print( "  --systype=%s !(bulk,molecule)" % systype_val)
print( "  --insulator  %s !not set this if you are not expert. (do not set for --systype=molecule)"    %  insulatorw)
print( "  --fsmom=%s ! (only for FSMOM mode. --systype=molecule automatically set this)"    %  fsmom_val)
print( "  --ssig=%s ! ScaledSigma(experimental =1.0 is the standard QSGW"    %  ssig_val)
#print( "  --ehmol ! if this exists, set EH used for a molecule paper (Not for PMT-QSGW. --ehmol may give better total energy in LDA)")
print( )
if(showhelp==1): sys.exit('--- end of help ---')


### readin atomlist ###
alist=atomlist.split("\n")
aused=open('atomlist.used','w')
dicatom={}
for line in alist:
    if(line[0:1]=='#' or len(line)==0): continue
    keya=(line.split("=")[0]).strip()
    dicatom[keya]=''.join(line)   #line.split("@")[0:-1]
    if(showatomlist): print( ''.join(line) )#line.split("@")[0:-1]
    aaa= '%s' % ''.join(line) #line.split("@")[0:-1]
    aused.write(aaa+'\n')
#for ikey in dicatom.keys():
#    print ikey,dicatom[ikey]
#sys.exit()
#
#line.split("atomz=")[1].split('@')[0]
#    mat=re.search('\".*\"',line)
#    aaa=mat.group().split('"')[1]
#    print keya,line.split("@")[2:-1]
if(showatomlist==1):
    print( "--- This is atomlist in ctrlgenM1.py (not yet set after Kr) ---.")
    sys.exit()
specstd=dicatom.keys()
#print specstd
#print dicatom
#keya='Fe'
#print dicatom[keya]
#rmt = float(dicatom[keya][-1].split('R=')[1])/.529177 * Rfac
#'rmt=',rmt,dicatom[keya]
#rmt = eval(dicatom[keya+'2dis'].split('@')[0].split('discenter=')[1])/2.*Rfac

z2dicatom={}
for ils in dicatom.keys():
    #print ils
    if(re.search('\.',ils)): # skip SPECKEY.# (such as Co.2) in z2dicatomlist.
        continue
    #print ils
    z2dicatom[dicatom[ils].split('atomz=')[1].split('@')[0]]=ils
#print dicatom.keys()
#sys.exit()


#### Read in ctrls #####
try:
    ext=sys.argv[1]
except:
    sys.exit()
print 
print( "... Generate ctrlgenM1.ctrl."+ext+" from ctrls."+ext + " ...")
ctrls = "ctrls." + ext
f=open(ctrls,'rt')
ctrlsdat = f.read() 
f.close()

listctrls  = lineReadfile(ctrls) 
listspec   = GetCat(listctrls,"SPEC")  # SPEC section only
print ('readin SPEC and #=',listspec,len(listspec))

listsite   = GetCat(listctrls,"SITE")  # SITE section only
liststruc  = GetCat(listctrls,"STRUC")  # SPEC section only 
listno     = RemoveCat2(RemoveCat2(RemoveCat2(listctrls,"SITE"),"SPEC"),"STRUC")  # sections except SITE and SPEC and STRUC

sitename = getsitename(listsite)
speclist = getsitename(listspec)
print( '### SITE  ', sitename)
print( '### SPEC  ', speclist)
print( '### other ', listno)

########### obtain mapping spec to Z dictionary spec2z ####
spec2z={}
specextra={}
zspec=False
if(len(listspec)==0):
    #listspec = re.split('\n',specstd)  # SPEC standard if no SPEC is in ctrls.*
    #print specstd
    #sys.exit()
    print( " NO SPEC is found in "+ctrls+". USE standard SPEC; try to see; ctrlgenM1.py --showatomlist")
#    for ils in sitename:
#        print 'site=',ils
#         spec2z[ils]=dicatom[ils].split('atomz=')[1].split('@')[0]
#         specname=ils.split('ATOM=')[1].split(' ')[0]
#         specz=ils.split('Z=')[1].split(' ')[0]
#         spec2z[specname]=specz
    specextradata=False
else:
    zspec=True
    for ils in listspec:
        atomsss=ils.split('ATOM=')
        print( 'readin spec=',len(atomsss),atomsss)
        if(len(atomsss)>1):
            specname=atomsss[1].split(' ')[0]
            specz=ils.split('Z=')[1].split(' ')[0]
            spec2z[specname]=specz
            sss =re.split(r'Z\s*=\s*[0-9]+\s',' '.join(atomsss))[-1].strip()
            specextra[specname]= sss
            print( ils,specname,specz)
    specextradata=True

#print spec2z
#sys.exit()

ansite = '%i' % len(sitename)
anspec = '%i' % len(uniq(sitename))

#### specsec0 for --getwsr (calculate touching MT radius) ###
specsec0=''
for ispec in uniq(sitename):
    aaa=''
#    aaa=aaa+ ' R='+ '%6.3f' %
    if(zspec):
        z=spec2z[ispec]
        speckey=z2dicatom[z]
    else:
        speckey=ispec
        z=dicatom[speckey].split('atomz=')[1].split('@')[0]

    aaa= '    ATOM='+ispec +' Z='+ z 
    specsec0= specsec0 + aaa +'\n'
    #print aaa

#print 'specsec0=\n',specsec0

################################################################
os.system("rm -rf llmchk_getwsr llmfa.tmp2")
head="### This is generated by ctrlgenM1.py from ctrls \n"
head= head+ """
###  Try "lmf --input", to see all effective category and token.
###
###  CAUTION: For GW calculations, you have to generate GWinput by mkGWIN_lmf2.
###          Repeat it when you change MTO settings ---> <PRODUCTBASIS> section is affected
###
VERS    LM=7 FP=7        # version check. Fixed.

IO      SHOW=T VERBOS=35 TIM=0,0 #Use TIM=3,3 or more for debug (to show which routines go through).
             # SHOW=T shows readin data (and default setting at the begining of console output)
         # It is useful to check ctrl is read in correctly or not (equivalent with --show option).
         # larger VERBOSE gives more detailed console output.

SYMGRP find   # 'find' evaluate space-group symmetry (just from lattice) automatically.
              #
              # Usually 'find is OK', but lmf may use lower symmetry
          # if you did not supply accurate structure.
              # 'lmchk foobar --pr60' shows, what symmetry is recognized.
              # See file://Document/Manual/CaterogyAndToken.org
              # To read its results.
              # If symmetry of electronic system has lower symmetry than the lattice symmetry,
              # we need to choose lower symmetry than 'find' gives.

%const kmxa=5  # =radial degree of freedom to expand eigenfuncitons in tail sites.
               #  kmxa=5 is good for pwemax \sim 4 or less.
               # larger kmxa is better but time-consuming. A rule of thumb: kmxa>pwemax in Ry.\n
               # Enlarge this, when your enlarge pwemax, and check little dependence on kmxa. 
"""
alltmp = head + glist(listno) \
     + glist(liststruc)  + "     NL=4  NBAS= "+ ansite + "  NSPEC="+ anspec +'\n' \
     + glist(listsite) \
     + 'SPEC\n' #specsec
f = open("ctrl.tmp",'wt')
f.write(alltmp+specsec0)
f.close()



#### Get touching MT radius and make rdic ########################
rlmchk=0
#print type(touchingratio)
if touchingratio>0: rlmchk=1
# ### Get R= by lmchk ###
if(rlmchk==1):
    os.system("lmchk --getwsr tmp > llmchk_getwsr; echo $? >exitcode")
    f=open("exitcode",'rt')
    iexit=int(f.read())
    f.close()
    try:
        listr = lineReadfile("rmt.tmp")
    except:
        print( ' Error: Can not readin rmt.tmp! ')
        sys.exit()
    print( listr)
    rdic={}
    for i in listr:
        xx=re.split(' +',i)
        print(xx)
        #rdic[xx[0]]=xx[1]
        rdic[xx[0]]= float(xx[1]) #*string.atof(touchingratio) #r_mul_val)
        rdic[xx[0]]= str(rdic[xx[0]])

#print " Rmax is taken from lmchk --getwsr. See llmchk_getwsr "
    print()
    print( ' rmt.tmp: --getwsr gives  R= -->', rdic)
    print( '  note: we use R=3.0 if R is larger than 3.0')
    print( 'rdic',rdic)
################################################################### 

#print dicatom
print ('zspec=',zspec)
#specdat = re.split('\WATOM=\W*',glist(listspec))[1:]
#print  'specdat',specdat
#specdic={}
#for i in specdat:
#    xx=re.split(' *',i)
#    ii=re.sub("\n","",i)
#    specdic[xx[0]]= '  ATOM='+ii +'\n'
###############
#print specdic
#print sitename,'yyy',uniq(sitename)
#print 'xxxx',sitename    
specsec=''
for ispec in uniq(sitename):
    aaa=''
#    aaa=aaa+ ' R='+ '%6.3f' %
    if(zspec):
        z=spec2z[ispec]
        speckey=z2dicatom[z]
    else:
        speckey=ispec
        z=dicatom[speckey].split('atomz=')[1].split('@')[0]

    try:
#### in the case of touching MT
#        print ispec,speckey,'tratio=',touchingratio
        print( 'from atomlist= ', dicatom[speckey])
        if touchingratio >0:
#            rrr = string.atof(rdic[speckey]) * touchingratio
            rrr = float(rdic[ispec]) * touchingratio
        else:
            rrr = float(getdataa(dicatom[speckey],'R='))/0.529177 #*r_mul_val
        print( rrr)
        if(rrr>3.0): rrr=3.0 #upper limit of R
        rrrh = rrr/2.0
        if(rrrh<0.5): rrrh=0.5 #lower limit of RSMH
    except:
        print ('ERROR: R are not set. need bug fix')
        sys.exit(-1)

    try:
        print('skey=',speckey)
        eh1data= getdataa( dicatom[speckey],'eh=')
        eh1value= eh1data.split('*')[0]+' '
        eh1count = int(eh1data.split('*')[1])
        print('   ',eh1data,eh1value,eh1count)

        eh2data= getdataa( dicatom[speckey],'eh2=')
        eh2value= eh2data.split('*')[0]+' '
        eh2count= int(eh2data.split('*')[1])
        print('   ',eh2data,eh2value,eh2count)

        if eh1set==1: eh1value='-1 '
        rsize= '%6.2f' % rrr
        rsizeh= '%6.2f' % rrrh

        print('  ',rsize,rsizeh)

        aaa= '    ATOM='+ispec +' Z='+ z + ' R='+rsize.strip()
        aaa=aaa+  ' '+getdataa2( dicatom[speckey],'pz=')
        aaa=aaa+  ' '+getdataa2( dicatom[speckey],'p=')+'\n'
        if specextradata: aaa=aaa+ ' '*6+specextra[ispec] +'\n'  #mar2013
        aaa=aaa+ '      EH='+  eh1value*eh1count
        aaa=aaa+ ' RSMH='+(rsizeh.strip()+' ')*eh1count+'\n'
        aaa=aaa+ '      EH2='+ eh2value*eh2count
        aaa=aaa+ ' RSMH2='+(rsizeh.strip()+' ')*eh2count+'\n' \
                +'      KMXA={kmxa}  LMX=3 LMXA=4 NMCORE=1\n' \
                            +'      '+mmom_val+' #s,p,d,f initial condition\n' \
                            +'      #NOTE: lmfa(rhocor) generates spin-averaged rho for any MMOM,jun2012\n'\
                            +'      #Q=0 0.5 1 0 #s,p,d,f initial condition \n' \
                            +'      #MMOM and Q are to set electron population. grep conf: in lmfa output\n'
    except:
        print( 'ERROR this is probably because we have not yet set default values in ctrlgenM1.py for a spec')
        sys.exit(-1)
    specsec= specsec + aaa +'\n'


### Read in "ctrls."+ext ###################################################
# try:
#     listctrl  = lineReadfile("ctrl.tmp") 
# except:
#     print '---> no ctrl or some problem'
#     sys.exit()
#print listctrl
### ctrl.tmp contains R=. Do lmfa to get mtopara.*
#f.write('\n'.join(listctrl) +'\nHAM  XCFUN=')
f = open("ctrl.tmp2",'wt')
f.write(alltmp+specsec +'\nHAM  XCFUN=')
f.write(xcfun_val+'\n')
f.close()


### check lmfa works OK or not #############################
os.system("lmfa tmp2 > llmfa.tmp2; echo $? >exitcode")
f=open("exitcode",'rt')
iexit=int(f.read())
f.close()
if (iexit != 0):
    print ('! Exit -1: not go through lmfa. You may need to modify SPEC for atoms not in atomlist.')
else:
    print( ' ------tail of llmf.tmp2 ----------------------------------')
    os.system("tail llmfa.tmp2")
    print( ' ----  lmfa has done! --------------------------------------')
    print ()

# mmmx = mtodic[ikey]
# mmm = re.sub(","," ",mmmx)
# pz=re.split("PZ",mmmx)
# #Over ride by new setting
# rsmh= max(string.atof(rdic[ikey])*1.0/2.0,0.5)
# rsmh= '%6.3f' % rsmh
# mmm= 'RSMH='+4*rsmh+' EH= -0.5 -0.5 -0.5 -0.5 \n'
# mmm= mmm+ '     RSMH2='+4*rsmh+' EH2= -2 -2 -2 -2\n'
# if(len(pz)==2):    mmm= mmm +'     PZ'+pz[1]

# il1 = countnum(mmm,'RSMH=')
# il2 = countnum(mmm,'PZ=')
# lx = max(il1,il2,3)
# lll = "%i" % lx
# #print il1,il2,lx
# aaa = aaa+' R='+rdic[ikey]+'\n'+' '*5+mmm+'\n' \
#     + ' '*5+'KMXA={kmxa} LMXA='+lll+'\n'+' '*5+'MMOM=0 0 0 0'


#########################################
metali_val= '%i' % metali
tail="""
\n"""
tail = tail+ "% const pwemax=3 nk1="+nk_val1+" nk2="+nk_val2+" nk3="+nk_val3+" nit=30  gmax=12  nspin="+nspin_val+ " metal="+ metali_val +" so=" +so_val +" xcfun="+xcfun_val+" ssig="+str(ssig_val)+"\n"
tail = tail + "BZ    NKABC={nk1} {nk2} {nk3} # division of BZ for q points.\n"\
            + "      METAL={metal}"\
"""
                # METAL=3 is safe setting (double path method but not repeat diagonalization), 
                # no problem even for insulator.
                # For insulator, METAL=0 may be a little faster.

              DOSMAX=1.5 NPTS=2001 SAVDOS=T
                # These are used to plot total dos, and pdos.
                # DOSMAX: is the maximum of the total dos plot. It is relative to the Fermi energy.
                #   Corresponding mimimum is automatically chosen to cover all valence states.
                # NPTS: division of plots. To get a high-energy resolution plot, use large NPTS.
                #
        # To get better plot, Use large NKABC (even after converged).
                # Larger NKABC may give smoother plot. Use large enough NKABC.
                #
                # NOTE: current version automatically enfoce TETRA=T and MEAL=3 for --pdos and --tdos
                # ---TOTAL DOS plot---      
                # you can use job_tdos_nspin*.
                # It just do lmf --tdos foobar; then you have dos.tot.foobar, and we can plot it.
                #  
                # --- PDOS plot --
                # you can do folloings by job_pdos_nspin*
                # run lmf --pdos foobar  (this stops quickly, or write SYMOPS (only e) in it by hand.)
                #     echo 1| qg4gw      (this is needed to generate QGpsi for no symmetry).
                #     lmf --pdos foobar  (it is a little time-consuming since we do not use symmetry).
                #     lmdos --pdos foobar

             # KNOWN BUG: For a hydrogen in a large cell, METAL=0 for (NSPIN=2 MMOM=1 0 0) 
             # results in non-magnetic solution. Use METAL=3 for a while in this case.

      # TETRA=0 N=-1 W=0.001 FSMOM below
      #  are for molecules. No tetrahedron integration. (Smearing)
"""

tail =  tail + "%const bzw=1e-4  fsmom="+str(fsmom_val)
if (systype_val.upper()=="BULK") :
    tail =  tail + \
"""
      #TETRA=0 
      #N=-1    #Negative is the Fermi distribution function W= gives temperature.
      #W=0.001 #W=0.001 corresponds to T=157K as shown in console output
               #W=0.01 is T=1573K. It makes stable convergence for molecule. 
               #Now you don't need to use NEVMX in double band-path method,
               #which obtain only eigenvalues in first-path to obtain integration weights
               #, and accumulate eigenfunctions in second path.

      #FSMOM={fsmom} real number (fixed moment method)
      #  FSMOM is for he fixed-spin moment method. 
      #  a spin-dependent potential shift is added to constrain the total magnetic moment to value 
      #  assigned by FSMOM=. Default is NULL (no FSMOM). FSMOM=0 works now (takao Dec2010)
      #  NOTE: current version is for ferro magnetic case (total mag moment) only.
      #FSMOMMETHOD=0 #only effective when FSMOM exists. #Added by t.kotani on Dec8.2010
      #  =0: original mode suitable for solids.(default)
      #  =1: discrete eigenvalue case. Calculate bias magnetic field from LUMO-HOMO gap for each spins.
      #      Not allowed to use together with HAM_SO=1 (L.S). 
      #      It seems good enough to use W=0.001. Smaller W= may cause instability.
"""
else :
    tail =  tail + \
"""
      TETRA=0 
      N=-1     #Negative is Fermi distribution function W= gives temperature.
      W={bzw}  #This corresponds to T=157K as shown in console output
               #W=0.01 is T=1573K. It makes stable nonvergence for molecule. 
               #Now you don't need to use NEVMX in double band-path method,
               #which obtain only eigenvalues in first-path to obtain integration weights
               #, and accumulate eigenfunctions in second path.
      FSMOM={fsmom} (fixed moment method)
      #  Set the global magnetic moment (collinear magnetic case). In the fixed-spin moment method, 
      #  a spin-dependent potential shift is added to constrain the total magnetic moment to value 
      #  assigned by FSMOM=. Default is NULL (no FSMOM). FSMOM=0 works now (takao Dec2010)
      # NOTE: current version is for ferro magnetic case (total mag moment) only.
      #
      FSMOMMETHOD=1 #only effective when FSMOM exists. #Added by t.kotani on Dec8.2010
      #  =0: original mode suitable for solids.(default)
      #  =1: discrete eigenvalue case. Calculate bias magnetic field from LUMO-HOMO gap for each spins.
      #      Not allowed to use together with HAM_SO=1 (L.S). 
      #      It seems good enough to use W=0.001. Smaller W= may cause instability.
"""


tail = tail + """      #For Molecule, you may also need to set FSMOM=n_up-n_dn, and FSMOMMETHOD=1 below.

      # NOINV=F (If F, we enfoce |psi_\sigma^\bfk|^2=|psi_\sigma^{-\bfk}|^2)
      #    usually judged automatically. So this is for debugging).
      # If NOINV=F, we assume H=H^* for each spin (k <-> -k symmetry). 
      #   This is the case in DFT without spin-orbit coupling.
      #   Then we use |psi_\sigma^\bfk|^2=|psi_\sigma^{-\bfk}|^2 to reduce comutational time.
      # If so/=0, lmf automatically sets NOINV=T.
      # You must get correct results evenif NOINV=T in the case of DFT without SO (for debug).
      #
      # WARN: Even if sigm exist, we use NOINV=F; exactly speaking, this is not correct.
      #       sigm without SO can cause orbital moments, that is, breaking time-reversal symmetry.
      # NOTE: because of inversion in space-group symmetry, we may have 
      #       |phi_sigm^\bfk|^2 = |phi_sigm^{-\bfk}|^2. This is not for NOINV.

ITER MIX=A2,b=.3,n=3 CONV=1e-5 CONVC=1e-5 NIT={nit}
#ITER MIX=B CONV=1e-6 CONVC=1e-6 NIT={nit}
                # MIX=A: Anderson mixing.
                # MIX=B: Broyden mixing (default). 
                #        Unstable than Anderson mixing. But faseter. It works fine for sp bonded systems.
                #  See file://Document/Manual/CaterogyAndToken.org

HAM   NSPIN={nspin}   # Set NSPIN=2 for spin-polarize case; then set SPEC_MMOM (initial guess of magnetic polarization).
      FORCES=0  # 0: no force calculation, 1: forces calculaiton 
      GMAX={gmax}   # this is for real space mesh. See GetStarted. (Real spece mesh for charge density).
                # Instead of GMAX, we can use FTMESH.
                # You need to use large enough GMAX to reproduce smooth density well.
                # Look into sugcut: shown at the top of console output. 
                # It shows required gmax for given tolelance HAM_TOL.
      REL=T     # T:Scaler relativistic, F:non rela.
"""
tail = tail + "      XCFUN={xcfun}"+ """
          # =1 for VWN.
                # =2 Birth-Hedin (if this variable is not set).
          #   (subs/evxc.F had a problem when =2 if rho(up)=0 or rho(down)=0).
                # =103 PBE-GGA

      PWMODE=1  # 0:  MTO basis only (LMTO) !2021feb. I think PWMODE=1 is usually better because of the smoothness.
                # 11: APW+MTO        (PMT)
                # 12: APW basis only (LAPW)
                #
      PWEMAX={pwemax} # (in Ry). When you use larger pwemax more than 5, be careful
                      # about overcompleteness. In cases, e.g.Fe, you may need larger KMXA for larger PWEMAX.
"""
if (False): #systype_val.upper()=="BULK") :
    tail = tail + """      ELIND=-1    # this is to accelarate convergence. Not affect to the final results.
"""
else :
    tail =  tail + """      ELIND=0    # this is to accelarate convergence. Not affect to the final results.
"""

tail = tail + """                 # For sp-bonded solids, ELIND=-1 may give faster convergence.
                 # For O2 molecule, Fe, and so on, use ELIND=0(this is default).

      FRZWF=F #If T, fix augmentation function. This is worth to test in future.
      #  See HAM_FRZWF, file://Document/Manual/CaterogyAndToken.org

      #For LDA+U calculation, see ecalj manual.

      #For QSGW. you have to set them. Better to get some samples.
      RDSIG=12
      SIGP[MODE=3]
      # Now we use no cutoff procedure for Sigma-Vxc in lmf. (only in emax_sigm is effective in GWinput). 
      # default:SIGP_EMAX=9999.
      RSRNGE=15.0  #If you see Exit -1 rdsigm: Bloch sum derivates mor..., Set this large enough.
      # This occurs when you set large n1n2n3 in GWinput, and/or large cell.
      # Reducing RSRNGE makes speed up a little ---> no effect to final results, but have chance 
      # to show the above message.

      ScaledSigma={ssig} # ScaledSigma* \Sigma + (1-ScaledSigma)*Vxc^LDA
      # This is a RPA-level hybridyzation method, in contrast to the B3LYP (hartree-fock level hybridyzation)

      SO={so}   #default = 0 
      #Spin-orbit coupling (for REL=1)
      #0 : no SO coupling
      #1 : Add L.S to hamiltonian (but non-colinear density yet).
      #2 : Add Lz.Sz only to hamiltonian
      #
      # For QSGW with SO=1 pertubation, run SO=1 with --rs=1,0 with NIT=1,
      #  starting from converged rst.* sigm.*  ESEAVR, and
      #  QpGpsi (it is easily generated by echo 1|qg4gw if you lost it).

      OVEPS=1d-7
      # nov2015 added. This is for removing poor linear-dependent basis in the diagonalization (H-eO)z=0
      # See cal zhev_tk3 in ecalj/lm7K/fp/bndfp.F. For using large pwemax, you may need to set OVEPS=1d-6 ~ 1d-8.
      # If you use larger OVEPS, you have smaller number of basis (APW+MTO) for expanding eigenfunctions.

OPTIONS PFLOAT=1 #not need to change this. Just for backward compatibility.

      # Q=band (this is quit switch if you like to add)


# Relaxiation sample
#DYN     MSTAT[MODE=5 HESS=T XTOL=.001 GTOL=0 STEP=.015]  NIT=20
# See file://Document/Manual/CaterogyAndToken.org

"""

### Write ctrl.ext
#g = open("ctrl."+ext,'wt')
#g.write(ctrlnospec+aaa+tail)
#g.close()
g = open("ctrlgenM1.ctrl."+ext,'wt')
#g.write('\n'.join(listctrl)+tail)
g.write(alltmp+specsec+tail)
g.close()
print( "OK! A template of ctrl file, ctrlgenM1.ctrl."+ext+", is generated.")
#print "Copy it to ctrl."+ext+"; then Do lmchk "+ext" to check crystal structure (option --pr60 show more info)."
sys.exit()
