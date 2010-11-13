#!/usr/bin/env python
#########################################################################
# Generate a temprate of ctrl file from ctrls.
# Work with python 2.5 or so.
#  Takao Kotani and Hiori Kino
# ---  a ctrls is ---
# HEADER  SrTiO3 cubic 
# STRUC   NBAS=5 NSPEC=3 ALAT=7.37 
#         DALAT=0 PLAT=1 0 0  0 1 0  0 0 1
# SITE
#   ATOM=Sr POS=1/2 1/2 1/2
#   ATOM=Ti POS= 0   0   0+0
#   ATOM=O  POS=1/2  0   0
#   ATOM=O  POS= 0  1/2  0
#   ATOM=O  POS= 0   0  1/2
# SPEC
#   ATOM=Sr Z=38 
#   ATOM=Ti Z=22 
#   ATOM=O  Z=8
# ---- end here ------
# (HEADER and so are at the begining of lines)
#########################################################################

import os
import sys
import string
import re

#-------------------------------------------------------
def  line2Token(linein):
	""" convert the result of readline to token """
#	""" input=readlines() output=token""
	listout = []
	for s in linein:
		if s=="":
			continue
		s = string.replace(s,'\n','')
		s = string.replace(s,',',' ')
                s = string.replace(s,'=',"= ")
		s = string.replace(s,':',": ")

		lista=string.split(s)
		for x in lista:
			if x<>"":
				listout.append(x)

	return listout

#---------------------------------------------------
def lineReadfile(filename):
	""" read file and make list of the content of the file, and \n->'' """
#	"input:filename output=readlines() "
	f = open(filename)
	list1 =[]
	while 1:
		s = f.readline()
		if s=="":
			break
		s=string.replace(s,"\n","")
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
#	print res
#	sys.exit()

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

#	res=res+'\n'
#	print res
#	sys.exit()
	return res

def GetCat(listctrl,key):
	""" get a category (key) from listctrl This returns lines. Not good correspondence to RemoveCat """
	res=[]
	ix=0
	sss ='^('+key.upper()+'|'+key.lower()+')'+'(\s|\Z)'
#	print sss
	for x in listctrl:
#		print x
		if(ix==1):
			if(re.match('^\w',x)): break
		if(re.match(sss, x )): ix=1
		if(ix==1): res.append(x)

#	res=res+'\n'
#	print res
#	sys.exit()
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

specstd=\
"""
SPEC (standard setting)
  ATOM=H   Z=1  
  ATOM=He  Z=2  
  ATOM=Li  Z=3  
  ATOM=Be  Z=4  
  ATOM=B   Z=5  
  ATOM=C   Z=6  
  ATOM=N   Z=7  
  ATOM=O   Z=8  
  ATOM=F   Z=9  
  ATOM=Ne  Z=10 
  ATOM=Na  Z=11 
  ATOM=Mg  Z=12 
  ATOM=Al  Z=13 
  ATOM=Si  Z=14 
  ATOM=P   Z=15 
  ATOM=S   Z=16 
  ATOM=Cl  Z=17 
  ATOM=Ar  Z=18 
  ATOM=K   Z=19 
  ATOM=Ca  Z=20 
  ATOM=Sc  Z=21 
  ATOM=Ti  Z=22 
  ATOM=V   Z=23 
  ATOM=Cr  Z=24 
  ATOM=Mn  Z=25 
  ATOM=Fe  Z=26 
  ATOM=Co  Z=27 
  ATOM=Ni  Z=28 
  ATOM=Cu  Z=29 
  ATOM=Zn  Z=30 
  ATOM=Ga  Z=31 
  ATOM=Ge  Z=32 
  ATOM=As  Z=33 
  ATOM=Se  Z=34 
  ATOM=Br  Z=35 
  ATOM=Kr  Z=36 
  ATOM=Rb  Z=37 
  ATOM=Sr  Z=38 
  ATOM=Y   Z=39 
  ATOM=Zr  Z=40 
  ATOM=Nb  Z=41 
  ATOM=Mo  Z=42 
  ATOM=Tc  Z=43 
  ATOM=Ru  Z=44 
  ATOM=Rh  Z=45 
  ATOM=Pd  Z=46 
  ATOM=Ag  Z=47 
  ATOM=Cd  Z=48 
  ATOM=In  Z=49 
  ATOM=Sn  Z=50 
  ATOM=Sb  Z=51 
  ATOM=Te  Z=52 
  ATOM=I   Z=53 
  ATOM=Xe  Z=54 
  ATOM=Cs  Z=55 
  ATOM=Ba  Z=56 
  ATOM=La  Z=57 
  ATOM=Ce  Z=58 
  ATOM=Pr  Z=59 
  ATOM=Nd  Z=60 
  ATOM=Pm  Z=61 
  ATOM=Sm  Z=62 
  ATOM=Eu  Z=63 
  ATOM=Gd  Z=64 
  ATOM=Tb  Z=65 
  ATOM=Dy  Z=66 
  ATOM=Ho  Z=67 
  ATOM=Er  Z=68 
  ATOM=Tm  Z=69 
  ATOM=Yb  Z=70 
  ATOM=Lu  Z=71 
  ATOM=Hf  Z=72 
  ATOM=Ta  Z=73 
  ATOM=W   Z=74 
  ATOM=Re  Z=75 
  ATOM=Os  Z=76 
  ATOM=Ir  Z=77 
  ATOM=Pt  Z=78 
  ATOM=Au  Z=79 
  ATOM=Hg  Z=80 
  ATOM=Tl  Z=81 
  ATOM=Pb  Z=82 
  ATOM=Bi  Z=83 
  ATOM=Po  Z=84 
  ATOM=At  Z=85 
  ATOM=Rn  Z=86 
  ATOM=Fr  Z=87 
  ATOM=Ra  Z=88 
  ATOM=Ac  Z=89 
  ATOM=Th  Z=90 
  ATOM=Pa  Z=91 
  ATOM=U   Z=92 
  ATOM=Np  Z=93 
  ATOM=Pu  Z=94 
  ATOM=Am  Z=95 
  ATOM=Cm  Z=96 
  ATOM=Bk  Z=97 
  ATOM=Cf  Z=98 
  ATOM=Es  Z=99 
  ATOM=Fm  Z=100
  ATOM=Md  Z=101
  ATOM=No  Z=102
  ATOM=Lr  Z=103
"""

#===============================================================================
#          MAIN 
#===============================================================================
nargv = len(sys.argv) -1
argset= set(sys.argv[1:])
help=0
helpatomname=0
if (nargv==0 or  '--help' in  argset): help=1
if ('--helpatomname' in  argset): helpatomname=1
version=" tkotani Nov12_2010"
if(help==1):
	print 
	print " Purpose: Generate ctrl.{ext} file from ctrls.{ext} ;"+version
	print
	print " Usage  : ctrlgen {extension of ctrl file} [option]"
	print "        : [options] = --helpatomname, --readrmt"
	print 
	print " You should have no SPEC if you use ATOM name shown by --helpatomname"
	print 
	print " If we supply rmt.tmp by hand with --readrmt, the give rmt.tmp is used."
        print " The rmt(Muffin-tin radius) are not calculated."
        print "      rmt.tmp consists of specname R with rmt for each line. For example,---"
        print '       --- rmt.tmp for SrTiO3 --- '
	print '       Sr          3.616323'
	print '       Ti          2.089960'
	print '       O           1.595007'
        print '       --- end of rmt.tmp ------- '
        print ' After you write rmt.tmp, do ctrlgen.py again'
	sys.exit()
if(helpatomname==1):
	print "--- This is a standard name when no SPEC is specified."
	print   specstd
	sys.exit()


#########################
ext=sys.argv[1]
print "Generate ctrl."+ext+" from ctrls."+ext + "..."

ctrls = "ctrls." + ext
if ('--readrmt' in  argset):
	readrmt=1
else:
	readrmt=0

#### Read in ctrls #####
f=open(ctrls,'rt')
ctrlsdat = f.read() 
f.close()

listctrls  = lineReadfile(ctrls) 
listspec   = GetCat(listctrls,"SPEC")  # SPEC section only
print 'listspec=',listspec
if(len(listspec)==0):
	listspec = re.split('\n',specstd)  # SPEC standard if no SPEC is in ctrls.*
	#print specstd
	#sys.exit()
	print " NO SPEC is found in " + ctrls + ". ---> USE standard SPEC; it is shown by ctrlgen.py without argument"
	print 
listsite   = GetCat(listctrls,"SITE")  # SITE section only
liststruc  = GetCat(listctrls,"STRUC")  # SPEC section only 
listno     = RemoveCat2(RemoveCat2(RemoveCat2(listctrls,"SITE"),"SPEC"),"STRUC")  # sections except SITE and SPEC and STRUC

sitename = getsitename(listsite)
#speclist = getsitename(listspec)
#print '### SITE  ', sitelist,len(sitelist)
#print '### SPEC  ', speclist,len(speclist)
#print '### other ', listno
ansite = '%i' % len(sitename)
anspec = '%i' % len(uniq(sitename))
specdat = re.split('\WATOM=\W*',glist(listspec))[1:]
specdic={}
for i in specdat:
	xx=re.split(' *',i)
	ii=re.sub("\n","",i)
	specdic[xx[0]]= '  ATOM='+ii +'\n'
###############33
#print specdic
#print sitename,'yyy',uniq(sitename)
#print 'xxxx',sitename	
specsec=''
for i in uniq(sitename):
	specsec= specsec + specdic[i]
#print specsec
#sys.exit()

##################
os.system("rm -rf llmchk_getwsr llmfa.tmp2")
head="### This is generated by ctrlgen.py from ctrls ;" + version +"\n"
head= head+ """
### For tokens, See http://titus.phy.qub.ac.uk/packages/LMTO/tokens.html. 
### However, lm7K is now a little different from Mark's lmf package in a few points.
### Do lmf --input to see all effective category and token ###
### It will be not so difficult to edit ctrlge.py for your purpose ###
VERS    LM=7 FP=7
             # version check. Fixed.

IO      SHOW=T VERBOS=35
             # SHOW=T shows readin data (and default setting at the begining of console output)
	     # It is useful to check ctrl is read in correctly or not (equivalent with --show option).
	     #
	     # lerger VERBOSE gives more detailed console output.

SYMGRP find  # 'find' evaluate space-group symmetry automatically.
             # Usually 'find is OK', but lmf may use lower symmetry
	     # because of numerical problem.
             # Do lmchk to check how it is evaluated.
             # See http://titus.phy.qub.ac.uk/packages/LMTO/tokens.html#SYMGRPcat
             
%const kmxa=5  # kmxa=5 is good for pwemax=3 or lower.
               # larger kmxa is better but time-consuming. A rule of thumb: kmxa>pwemax in Ry.\n
"""

alltmp = head + glist(listno) \
	 + glist(liststruc)  + "  NBAS= "+ ansite + "  NSPEC="+ anspec +'\n' \
	 + glist(listsite) \
	 + 'SPEC\n'+specsec

#print alltmp

f = open("ctrl.tmp",'wt')
f.write(alltmp)
f.close()

### Get R= by lmchk
if(readrmt==0):
	os.system("lmchk --getwsr tmp > llmchk_getwsr; echo $? >exitcode")
	f=open("exitcode",'rt')
	iexit=int(f.read())
	f.close()
else:
	iexit=2

if(iexit==0):
	print ' -----tail of llmchk_getwsr ----------------------------'
	os.system("tail llmchk_getwsr")
	print ' --- lmchk --getwsr has done successfully! rmt.tmp is generated -----'
elif(iexit==2):
	print ' readin rmt.tmp! because of --readrmt '
else:
	print ' lmchk --getwsr can not find the muffin-tin radius SPEC_ATOM_R.'
	print ' Wrong ctrls? or bug in ctrlgen.py? (see llmchk_getwsr)'
        print '   or you have to write rmt.tmp by yourself.'
	print ' rmt.tmp consists of "specname R" for each line. For example,---'
        print '       --- rmt.tmp for SrTiO3 --- '
	print '       Sr          3.616323'
	print '       Ti          2.089960'
	print '       O           1.595007'
        print '       --- end of rmt.tmp ------- '
        print ' After you write rmt.tmp, repeat ctrlgen.py --readrmt'
        print ' Then, check that "R= -->" below shows the content of rmt.tmp'
        print ' If so, neglect this error.'
	print 

### Read rmax, and make rdic ##############################################
try:
	listr = lineReadfile("rmt.tmp")
except:
	print ' Error: Can not readin rmt.tmp! '
	sys.exit()
	
rdic={}
for i in listr:
	xx=re.split(' *',i)
	#print xx
	rdic[xx[0]]=xx[1]
#print " Rmax is taken from lmchk --getwsr. See llmchk_getwsr "
print 
print ' rmt.tmp: gives  R= -->', rdic
print ' we use R=3.0 if R is larger than 3.0'
#sys.exit()



### Read in "ctrls."+ext ######################################################
try:
	listctrl  = lineReadfile("ctrl.tmp") 
except:
	print '---> no ctrl or some problem'
	sys.exit()
ictrlnospec = RemoveCat2(listctrl,"SPEC")
ctrlnospec=''
for i in ictrlnospec:
	ctrlnospec=ctrlnospec+i+'\n'

### spec section
listspec =  GetCat(listctrl,"SPEC")
tokenspec=line2Token(listspec)
#print '### ', tokenspec
if(tokenspec[0]!='SPEC'):
	print 'tokenspec[0]!=SPEC :ctrls is wrong or bug?'
	sys.exit()


### Generate ctrl.tmp2 including R=
aaa=''
ispec=0
ix=0
aaa='SPEC\n  '
for ii in tokenspec[1:]:
	if(ix==1):
		ikey=ii
		#print ikey
		ix=0
	elif(ii=='ATOM='):
		if(ispec>0): aaa = aaa+' R= '+rdic[ikey]+'\n      '
		ispec=ispec+1
		ix=1
	aaa=aaa+' '+ii
aaa = aaa+' R='+rdic[ikey]+'\n'

### ctrl.tmp2 contains R=. Do lmfa to get mtopara.*
f = open("ctrl.tmp2",'wt')
f.write(ctrlnospec+aaa +'\nHAM  XCFUN=1\n')
f.close()
os.system("lmfa tmp2 > llmfa.tmp2; echo $? >exitcode")
f=open("exitcode",'rt')
iexit=int(f.read())
f.close()
if (iexit != 0):
	print '! Exit -1: not go through lmfa. Wrong ctrls?'
else:
	print ' ------tail of llmf.tmp2 ----------------------------------'
	os.system("tail llmfa.tmp2")
	print ' ----  lmfa has done! --------------------------------------'
	print 

#NOTE: "lmfa tmp2 >& llmfa.tmp2" caused a  bug in python. "Bad fd number"

### Readin mtopara
ext2='tmp2'
try:
	listmto = lineReadfile("mtopara."+ext2)
except:
	print "---> No mtopara."+ext2, "or some problem"
	sys.exit()
mtodic={}
for i in listmto:
	i=re.sub("KMXA=","\n    KMXA={kmxa} LMXA=5",i)
	xx=i.split('@')
	mtodic[xx[0]]= xx[1]

### Get new SPEC section (taken setting from mtopara).
ispec=0
ix=0
aaa='SPEC\n'
tokenspec.append("ENDATOM")
for ii in tokenspec[1:]:
	#print 'xxxxxxxxxxx',ii
	if(ix==1):
		ikey=ii
		ix=0
	elif(ii=='ATOM=' or ii=='ENDATOM'):
		if(ispec>0): 

			mmmx=mtodic[ikey]
			mmm=re.sub(","," ",mmmx)
			#il1 = countnum(mmm,'RSMH=')
			il1=4
			pz=re.split("PZ",mmmx)
                        #Over ride by new setting
			rmt = min(string.atof(rdic[ikey]),3.0)
			rsmh = max(rmt*1.0/2.0,0.5)
			rsmh= '%6.3f' % rsmh
			mmm= 'RSMH=  '          +il1*rsmh+' EH= '+il1*' -1.0'+'\n' #-1 and -0.5 which is better? 
			mmm= mmm+ '     RSMH2= '+il1*rsmh+' EH2='+il1*' -2.0'+'\n'
			if(len(pz)==2):	mmm= mmm +'     PZ'+pz[1]
			
			il2 = countnum(mmm,'PZ=')
			lx = max(il1,il2)
			lll = "%i" % lx
			il1m=il1-1
			lmx = "%i" % il1m
			rmt= '%6.3f' % rmt
			#print il1,il2,lx
			aaa = aaa+' R='+rmt +'\n'+' '*5+mmm \
			    +' '*5+'KMXA={kmxa} '+' LMX='+lmx+' LMXA='+lll+'\n'+'     MMOM=0 0 0 0'+'\n\n'
		ispec=ispec+1
		ix=1
	if(ii=='ENDATOM'): break
	aaa=aaa+' '+ii

# mmmx = mtodic[ikey]
# mmm = re.sub(","," ",mmmx)
# pz=re.split("PZ",mmmx)
# #Over ride by new setting
# rsmh= max(string.atof(rdic[ikey])*1.0/2.0,0.5)
# rsmh= '%6.3f' % rsmh
# mmm= 'RSMH='+4*rsmh+' EH= -0.5 -0.5 -0.5 -0.5 \n'
# mmm= mmm+ '     RSMH2='+4*rsmh+' EH2= -2 -2 -2 -2\n'
# if(len(pz)==2):	mmm= mmm +'     PZ'+pz[1]

# il1 = countnum(mmm,'RSMH=')
# il2 = countnum(mmm,'PZ=')
# lx = max(il1,il2,3)
# lll = "%i" % lx
# #print il1,il2,lx
# aaa = aaa+' R='+rdic[ikey]+'\n'+' '*5+mmm+'\n' \
#     + ' '*5+'KMXA={kmxa} LMXA='+lll+'\n'+' '*5+'MMOM=0 0 0 0'



tail="""
\n
% const pwemax=2 nk=2 nit=30 gmax=12 nspin=1
BZ    NKABC={nk} {nk} {nk}  # division of BZ for q points.
      METAL=3   # METAL=3 is safe setting. For insulator, METAL=0 is good enough.
		# When you plot dos, set SAVDOS=T and METAL=3, and with DOS=-1 1 (range) NPTS=2001 (division) even for insulator.
		#   (SAVDOS, DOS, NPTS gives no side effect for self-consitency calculaiton).
                # 
                #BUG: For a hydrogen in a large cell, I(takao) found that METAL=0 for
                #(NSPIN=2 MMOM=1 0 0) results in non-magnetic solution. Use METAL=3 for a while in this case.
                # 

      BZJOB=0	# BZJOB=0 (including Gamma point) or =1 (not including Gamma point).
		#  In cases , BZJOB=1 makes calculation efficient.


      #Setting for molecules. No tetrahedron integration. (Smearing))
      #TETRA=0 
      #N=-1
      #W=0.001


      #For Total DOS.   DOS:range, NPTS:division. We need to set METAL=3 with default TETRA (no TETRA).
      #SAVDOS=T DOS=-1 1 NPTS=2001

      #EFMAX= (not implemented yet, but maybe not so difficult).            


      #  See http://titus.phy.qub.ac.uk/packages/LMTO/tokens.html#HAMcat for tokens below.

      #NOINV=T (for inversion symmetry)
      #  Suppress the automatic addition of the inversion to the list of point group operations. 
      #  Usually the inversion symmetry can be included in the determination of the irreducible 
      #  part of the BZ because of time reversal symmetry. There may be cases where this symmetry 
      #  is broken: e.g. when spin-orbit coupling is included or when the (beyond LDA) 
      #  self-energy breaks time-reversal symmetry. In most cases, lmf program will automatically 
      #  disable this addition in cases that knows the symmetry is broken
      #

      #FSMOM=real number (fixed moment method)
      #  Set the global magnetic moment (collinear magnetic case). In the fixed-spin moment method, 
      #  a spin-dependent potential shift is added to constrain the total magnetic moment to value 
      #  assigned by FSMOM=. No constraint is imposed if this value is zero (the default).
      #

      #INVIT=F
      #  Enables inverse iteration generate eigenvectors (this is the default). 
      #  It is more efficient than the QL method, but occasionally fails to find all the vectors. 
      #   When this happens, the program stops with the message:
      #     DIAGNO: tinvit cannot find all evecs
      #   If you encounter this message set INVIT=F.
      #  T.Kotani think (this does not yet for lm7K).
    
ITER MIX=A2,b=.5,n=3 CONV=1e-6 CONVC=1e-6 NIT={nit}
#ITER MIX=B CONV=1e-6 CONVC=1e-6 NIT={nit}
                # MIX=A: Anderson mixing.
                # MIX=B: Broyden mixing (default). 
                #        Unstable than Anderson mixing. But faseter. It works fine for sp bonded systems.
                #  See http://titus.phy.qub.ac.uk/packages/LMTO/tokens.html#ITERcat

HAM   NSPIN={nspin}   # Set NSPIN=2 for spin-polarize case; then set SPEC_MMOM (initial guess of magnetic polarization).
      FORCES=0  # 0: no force calculation, 1: forces calculaiton 
      GMAX={gmax}   # this is for real space mesh. See GetStarted. (Real spece mesh for charge density).
                # Instead of GMAX, we can use FTMESH.
                # You need to use large enough GMAX to reproduce smooth density well.
                # Look into sugcut: shown at the top of console output. 
                # It shows required gmax for given tolelance HAM_TOL.
      REL=T     # T:Scaler relativistic, F:non rela.

      XCFUN=1   # =1 for VWN.
                # =2 Birth-Hedin (if this variable is not set).
		#    (subs/evxc.F had a problem when =2 if rho(up)=0 or rho(down)=0).
                # =103 PBE-GGA

      PWMODE=11 # 10: MTO basis only (LMTO) PW basis is not used.
                # 11: APW+MTO        (PMT)
                # 12: APW basis only (LAPW) MTO basis is not used.

      PWEMAX={pwemax} # (in Ry). When you use larger pwemax more than 5, be careful
                      # about overcompleteness. See GetStarted.

      ELIND=0    # this is to accelarate convergence. Not affect to the final results.
                 # For sp-bonded solids, ELIND=-1 may give faster convergence.
                 # For O2 molecule, Fe, and so on, use ELIND=0(this is default).
  
      #STABILIZE=1e-10 #!!! Test option for convergence check. Not tested well.
                       # default is negative, then STABILIZER in diagonalization is not effective 
                       # (See slatsm/zhev.F delta_stabilize).
                       # I am not sure wether this stabilizer works OK or not(in cases this gives little help).
                       # STABILIZE=1e-10 may make convergence stabilized 
                       # (by pushing up poorly-linear-dependent basis to high eigenvalues).
                       # STABILIZE=1e-8 may give more stable convergence. 
                       # If STABILIZE is too large, it may affect to low eigenvalues around E_Fermi

      #FRZWF=T #to fix augmentation function. 
      #  See http://titus.phy.qub.ac.uk/packages/LMTO/tokens.html#HAMcat

      #For LDA+U calculation, see http://titus.phy.qub.ac.uk/packages/LMTO/fp.html#ldaplusu

      #For QSGW. you have to set them. Better to get some samples.
      #RDSIG=
      #RSRNGE=
               
OPTIONS PFLOAT=1 
        # Q=band (this is quit switch if you like to add)

# Relaxiation sample
#DYN     MSTAT[MODE=5 HESS=T XTOL=.001 GTOL=0 STEP=.015]  NIT=20
# See http://titus.phy.qub.ac.uk/packages/LMTO/tokens.html#DYNcat

"""

### Write ctrl.ext
#g = open("ctrl."+ext,'wt')
#g.write(ctrlnospec+aaa+tail)
#g.close()
g = open("ctrlgen.ctrl."+ext,'wt')
g.write(ctrlnospec+aaa+tail)
g.close()
print "OK! A template of ctrl file, ctrlgen.ctrl."+ext+", is generated."
sys.exit()
