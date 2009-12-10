#!/usr/bin/python
#########################################################################
# Generate a temprate of ctrl file from ctrls.
#
#  T.Kotani. March.2009.
#  this worked by python 2.5.2
#
#  By Hiori Kino, Aug, 2004
#  modified: Sep, 24, 2004 by H.Kino
#  
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

# import textwrap
# version 2.3 or higher is necessary to import textwrap


# listctrlkey=["HEADER","VERS","IO","CONST","STRUC","OPTIONS","BZ","SYMGRP","SPEC","SITE",
#                 "EWALT","MIX","HAM","START","GW"]

# #----------------------------------------------------------------------------
# def token2Line(tokena):
# 	""" connect a tokana list and return strings to write to the file """
# 	keylist=tokenMakeKeylist(tokena)
# 	n=len(tokena)
# 	listout=[]
# 	str1 = "%-8s" % (tokena[0])
# 	newone=0
# 	for i in range(1,n):
# 		str1 = str1 + tokena[i]+" "
# 		newone=1
# 		if len(str1)> 50:
# 			iskey=0
# 			for x in keylist:
# 				if x==tokena[i]:
# 					iskey=1
# 			if iskey==0:
# 				listout.append(str1)
# 				str1="        "
# 				newone=0
# 	if newone==1:
# 		listout.append(str1)
# 	return listout
		


# #------------------------------------------------------------------------------
# def tokenExpandMtoq(atom,specatom,listmtofp):
# 	""" replace MTOQ secction with the content of listmtofp """
# 	speckey=tokenMakeKeylist(specatom)
# 	mtoqlist=tokenExtract(specatom,"MTOQ=","",speckey)
#         specatom=tokenRemove(specatom,mtoqlist[0])
# 	mtoqlist.remove("MTOQ=")
# 	speclist=lineExtract(listmtofp,"SPEC",atom,[])
# 	l=0
# 	casedata=[]
# 	for x in mtoqlist:
# 		caselist = lineExtract(speclist,"CASE",x,[])
# 		if caselist==[]:
# 			print
# 			print "Error: Invalid CASE ",x, " is not found at MTOQ= section of ATOM=",atom
# 			print
# 			sys.exit(10)
# 		casedata.append( caseExtract(caselist,l) )
# 		l=l+1
	
# 	return tokenArrangeMtoq(casedata)
	
	
# #------------------------------------------------------------------------------
# def tokenArrangeMtoq(casedata):
# 	""" make EH= -1.0 -1.0 0.0,,, RSMH= 2.0 2.4 1.0 ... """
# 	""" from RH=1.0 RSMH=2.0.. EH=-1.0 RSMH=2.4 ... """
# 	flageh=0
# 	flagrsmh=0
# 	flageh2=0
# 	flagrsmh2=0
# 	flagpz=0
# 	flagp=0
# 	flagidmod=0
# 	for x in casedata:
# 		for st in x:
# 			if st=="EH=":	
# 				flageh=1
# 			if st=="RSMH=":
# 				flagrsmh=1
# 			if st=="EH2=":
# 				flageh2=1
# 			if st=="RSMH2=":
# 				flagrsmh2=1
# 			if st=="P=":
# 				flagp=1
# 			if st=="PZ=":
# 				flagpz=1
# 			if st=="IDMOD=":
# 				flagidmod=1
# 			if x==[]:
# 				break

# 	lmx = len(casedata)-1
	

# 	list1=[]
# 	defaultvalueeh="0.00"
# 	defaultvaluersmh="-1.00"
# 	defaultvaluep="0.00"
# 	defaultvaluepz="0.00"
# 	defaultvalueidmod=["0","0","0","1","1","1","1"]
	
# # EH
# 	if flageh==1:
# 		list1.append("EH=")
# 		for x in casedata:
# 			value=tokenKeyvalue(x,"EH=")	
# 			if value<>[]:
# 				list1.append(value[0])
# 			else:
# 				list1.append(defaultvalueeh)
# 		str= defaultvalueeh
# 		while str==defaultvalueeh:
# 			str=list1.pop()
# 		list1.append(str)
		
# # RSMH
# 	if flagrsmh==1:
# 		list1.append("RSMH=")
# 		for x in casedata:
# 			value=tokenKeyvalue(x,"RSMH=")	
# 			if value<>[]:
# 				list1.append(value[0])
# 			else:
# 				list1.append(defaultvaluersmh)
# 		str=defaultvaluersmh
# 		while str==defaultvaluersmh:
# 			str=list1.pop()
# 		list1.append(str)
		
# #EH2
# 	if flageh2==1:
# 		list1.append("EH2=")
# 		for x in casedata:
# 			value=tokenKeyvalue(x,"EH2=")	
# 			if value<>[]:
# 				list1.append(value[0])
# 			else:
# 				list1.append(defaultvalueeh)
# 		str=defaultvalueeh
# 		while str==defaultvalueeh:
# 			str=list1.pop()
# 		list1.append(str)

# # RSMH
# 	if flagrsmh2==1:
# 		list1.append("RSMH2=")
# 		for x in casedata:
# 			value=tokenKeyvalue(x,"RSMH2=")	
# 			if value<>[]:
# 				list1.append(value[0])
# 			else:
# 				list1.append(defaultvaluersmh)
# 		str=defaultvaluersmh
# 		while str==defaultvaluersmh:
# 			str=list1.pop()
# 		list1.append(str)
# # P
# 	if flagp==1:
# 		list1.append("P=")
# 		for x in casedata:
# 			value=tokenKeyvalue(x,"P=")	
# 			if value<>[]:
# 				list1.append(value[0])
# 			else:
# 				list1.append(defaultvaluep)

# # IDMOD
# 	if flagidmod==1:
# 		list1.append("IDMOD=")
# 		l=0
# 		for x in casedata:
# 			value=tokenKeyvalue(x,"IDMOD=")	
# 			l=l+1
# 			if value<>[]:
# 				list1.append(value[0])
# 			else:
# 				list1.append(defaultvalueidmod[l])

# #PZ
# 	if flagpz==1:
# 		list1.append("PZ=")
# 		for x in casedata:
# 			value=tokenKeyvalue(x,"PZ=")
# 			if value<>[]:
# 				list1.append(value[0])
# 			else:
# 				list1.append(defaultvaluepz)


	
	
# 	return list1

#-------------------------------------------------------------------------------
#def tokenRemoveList(list1,listremove):
#	""" remove listremove from list1, list1 itself is changed, nothing returns """
#
#	n =len(list1)
#	pos=-1
#	for i in range(0,n):
#		if list1[i]==listremove[0]:
#			pos=i	
#			break
#	if pos>=0:
#		n=len(listremove)
#		for i in range(0,n):
#			list1.pop(pos)		
#
#-------------------------------------------------------------------------------
#def tokenRemoveN(list1,name,n):
#	""" remove n tokens staring name , list1 is changed, nothing returns """
#	n0 = len(list1)
#	pos=-1
#	for i in range (0,n0):
#		if list1[i]==name:
#			pos=i
#			break
#	if pos>=0:
#		for i in range(0,n+1):
#			list1.pop(pos)
#
#-------------------------------------------------------------------------------
# def tokenRemove(list1,name):
# 	""" remove a block starting name , list1 is changed, nothing returns"""
# 	listkey=tokenMakeKeylist(list1)
# 	listname=tokenExtract(list1,name,"",listkey)
# 	n0 = len(list1)
# 	pos=-1
# 	for i in range (0,n0):
# 		if string.upper(list1[i])==string.upper(name):
# 			pos=i
# 			break
# 	if pos>=0:
# 		for i in range(0,len(listname)):
# 			list1.pop(pos)

	
		
# #------------------------------------------------------------------------------
# def tokenPrintSpec(list1):
# 	"""  print token list1  which is spec section of a ctrl file """
# 	NLlist=["EH=","RSMH=","EH2=", "RSMH2=", "IDMOD=", "P=","PZ="]	
# 	clearlist=["ATOM="]
# 	nlflag=0
# 	for x in list1:
# 		if x=="":
# 			continue
# 		for y in NLlist:
# 			if  y == x:
# 				nlflag=1
# 				sys.stdout.write("\n           ")
# 		for y in clearlist:
# 			if y==x:
# 				nlflag=0
# 				sys.stdout.write("        ")
				

# 		if nlflag==1:
# 			sys.stdout.write("%6s" % (x))
# 			sys.stdout.write(" ")
# 		else:
# 			sys.stdout.write(x)
#                         sys.stdout.write(" ")


# 	print


# #------------------------------------------------------------------------------
# def tokenExpandQuality(specatom,listmtofp):
# 	""" change specatom, replace quality=??? with the result of listmtofp """
# 	qualitylist=tokenKeyvalue(specatom,"quality=")
# 	qualityvalue=""
# 	if qualitylist<>[]:
# 		qualityvalue=qualitylist.pop()
# 		qualitydetail=findQuality(listmtofp, atom,qualityvalue)
# 		tokenRemove(specatom,"quality=")
# 		specatom=tokenMergelist(specatom,qualitydetail)
# 	return specatom


# #------------------------------------------------------------------------------
# def tokenKeyvalue(token,key):
# 	""" get token followed by key """
# 	compare=0
# 	listout=[]
# 	for x in token:
# 		if compare==1:
# 			listout.append(x)
# 			compare=0
# 		if string.upper(x)==string.upper(key):
# 			compare=1

# 	return listout



# #---------------------------------------------------

# def tokenMakeKeylist(tokena):
# 	""" list up keywords ending at = """
# 	listout=[]
# 	for x in tokena:
# 		n = len(x)
# 		if x[n-1]=="=" or x[n-1]==":":
# 			listout.append(x)
# 	return listout

# #---------------------------------------------------------
# def tokenMergelist(tokena,tokenb):
# 	""" merge a tokena list  and a tokenb list, if the key of the tokanb list does not exist in the tokena list """
# 	keya=tokenMakeKeylist(tokena)
# 	keyb=tokenMakeKeylist(tokenb)
# 	for a in keya:
# 		for b in keyb:
# 			if string.upper(a)==string.upper(b):
# 				listb=tokenExtract(tokenb, b,"", keyb)
# 				tokenRemove(tokenb,listb[0])	

# 	return tokena+tokenb


#----------------------------------------------------
# def tokenExtract(tokenin,name,key,namelist):
# 	""" extract token staring name following key and ending one of namelist or name """
# 	""" key can be null and namelist can be null list """
# #	""" input=list output=list """
# 	listout=[]
# 	start=0
# 	compare=0
# 	s1=""
# 	for s in tokenin:
# 		if compare==1:
# 			if string.upper(s)==string.upper(key):
# 				start=1
# 				listout.append(s1)
# 				listout.append(s)
# 				compare=0
# 				continue
# 		if string.upper(s)==string.upper(name):
# 			if start==1:
# 				return listout
# 			else:
# 				if len(key)>0:
# 					compare=1	
# 					s1=s
# 					continue
# 				else:
# 					start=1
# 					listout.append(s)
# 		else:
# 			if start==1 and namelist<>[]:
# 				for x in namelist:
# 					if string.upper(x)==string.upper(s):
# 						return listout
# 			if start==1:
# 				listout.append(s)
# 	return listout


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

#-----------------------------------------------------

#-----------------------------------------------------

# def lineExtract(linein,name,key,namelist):
# 	""" extract list staring name following key and ending name of namelist """
# 	""" key can be null and namelist can be null """
# #	""" input=readlines() output=lines"""
# 	listout=[]
# 	start=0
# 	for s in linein:
# 		ret=string.split(s)
# 		if ret==[]:
# 			continue
# 		st=string.upper(ret.pop(0))
# 		if st==string.upper(name):
# 			if start==1:
# 				return listout
# 			else:
# 				if len(key)>0:
# 					st = string.upper(ret.pop(0))
# 					if st==string.upper(key):
# 						start=1
# 						listout.append(s)	
# 				else:
# 					start=1
# 					listout.append(s)
# 		else:
# 			if start==1 and namelist<>[]:
# 				for x in namelist:
# 					if st==string.upper(x):
# 						return listout
# 			if start==1:
# 				listout.append(s)

# 	return listout

# def lineRemove(linein,name,key,namelist):
# 	""" extract list staring name following key and ending name of namelist """
# 	""" key can be null and namelist can be null """
# #	""" input=readlines() output=lines"""
# 	listout=[]
# 	start=0
# 	for s in linein:
# 		ret=string.split(s)
# 		if ret==[]:
# 			continue
# 		st=string.upper(ret.pop(0))
# 		if st==string.upper(name):
# 			if start==1:
# 				return listout
# 			else:
# 				if len(key)>0:
# 					st = string.upper(ret.pop(0))
# 					if st==string.upper(key):
# 						start=1
# 					else:	
# 						listout.append(s)	
# 				else:
# 					start=1
# 		else:
# 			if start==1 and namelist<>[]:
# 				for x in namelist:
# 					if st==string.upper(x):
# 						return listout
# 			if start==1:
# 				pass
# 			else:	
# 				listout.append(s)

# 	return listout
#------------------------------------------------------------------------------

# def listCopy(listin):
# 	""" copy list """
# 	listout=[]
# 	for x in listin:
# 		listout.append(x)

# 	return listout


			
#-----------------------------------------------------------
# def findQuality(linemtofp, atom,qualityvalue):
# 	""" get qualist section of atom , quality = qualityvalue """
# 	listspec=lineExtract(linemtofp,"spec", atom,[])
# 	listquality=lineExtract(listspec,"quality",qualityvalue,["quality","case"])
# 	if listquality==[]:
# 		print
# 		print "Error: quality=",qualityvalue, " is not found in ATOM=",atom 
# 		print 
# 		sys.exit(10)
# 	tokenquality=line2Token(listquality)
# 	listout=[]
# 	for x in tokenquality:
# 		if string.upper(x)<>"QUALITY" and string.upper(x)<>string.upper(qualityvalue):
# 			listout.append(x)

# 	return listout


#-------------------------------------------
#def RemoveKeyandValue(list1,key,n):
#	""" make a new list where key and n token after key is removed """
#	flag=0
#	listout=[]
#	for x in list1:
#		if flag==1:
#			n=n-1
#			if n==0:
#				flag=0
#			continue
#		if string.upper(x)==string.upper(key):
#			flag=1
#			continue
#		listout.append(x)
#				
#	return listout

#------------------------------------------------------------------------------
# def caseExtract(caselist,l):
# 	" make a list such as EH= -1.0 RSMH= 3.0 EH2=... , for l"
# 	listout=[]
# 	casetoken=line2Token(caselist)
# 	casekey=tokenMakeKeylist(casetoken)
# 	for token in casekey:
# 		token1=tokenExtract(casetoken,token,"",casekey)
# 		if token1<>[]:
# 			token1.remove(token)
# 			l0=0
# 			value=""
# 			for x in token1:
# 				if l0==l:
# 					value=x
# 					break
# 				l0=l0+1
# 			listout.append(token)
# 			listout.append(value)
# 	return listout


# #-------------------------------------------------------------------------
# def LineRemoveSection(listctrl,key,namelist):
# 	""" remove key section from listctrl """
# 	start=-1
# 	end=len(listctrl)
# 	i=-1
# 	for x in listctrl:
# 		print 'aaa ', x
# 		i=i+1
# 		x1=[x]
# 		tok = line2Token(x1)
# 		if tok==[]:
# 			continue
# 		if start==-1 and string.upper(tok[0])== string.upper(key):
# 			start=i	
# 			continue
# 		if start>=0:
# 			for y in namelist:
# 				if string.upper(y)==string.upper(tok[0]):
# 					end=i
# 					break
# 		if start>=0 and end>=0:
# 			break
# 	if start>=0 and end>=0:
# 		for i in range(start,end):
# 			listctrl.pop(start)


# #---------------------------------------------------------------------------
# def  lineWriteCtrl(listctrl,listspecatom):
# 	""" write listctrl and listspecatom """
# 	for x in listctrl:
# 		print x
# 	print
# 	print "SPEC"
# 	for x in listspecatom:
# 		tokenPrintSpec(x)
# 	print
	

# #--------------------------------------------------------------------------
# def lineAddCtrlFeature(listctrl):
# 	""" add keyword lists to ctrl """
# 	listkeyfeature=[
# 		"IO      SHOW=F HELP=F VERBOS=30 ",
# 		"OPTIONS NSPIN=2 REL=T XCFUN=2",
# 		"BZ      NKABC=2 2 2 TETRA=1 BZJOB=0 TETRA=1 METAL=3 ",
# 		"EWALD   TOL=1D-8 NKDMX=1999 NKRMX=1999",
# 		"MIX     CONV=1d-4 CONVC=1d-4",
# 		"HAM     ELIND=-.7 TOL=1d-6  GMAX=12.5                   "+
# 		"        RDSIG=12 SIGP:3,0,0,0,2.,0,.06,0",
# 		"START   NIT=100",
# 		"GW      NKABC= GCUTB= GCUTX="]
# 	listkey=["IO","OPTIONS","BZ","EWALD","MIX","HAM","START","GW"]

# 	listout=[]
# 	for x in listctrl:
# 		t = line2Token([x])
# 		if t==[]:
# 			listout.append(x)
# 			continue
# 		found=0
# 		for y in listkey:
# 			if string.upper(y)==string.upper(t[0]):
# 				listy=lineExtract(listctrl,y,"",listctrlkey)
# 				tokeny=line2Token(listy)
# 				listadd=listkeyfeature[ listkey.index(y) ]
# 				tokenadd=line2Token([listadd])
# 				tokeny.pop(0)
# 				tokenadd.pop(0)
# 				tokennewy= tokenMergelist(tokeny,tokenadd)
# 				tokennewy.insert(0,y)
# 				found=1
# 				break
# 		if found==0:
# #			list1= textwrap.wrap(x)
# 			list1= textwrap_wrap(x)
# 			for i in range(0,len(list1)):
# 				if i==0:
# 					listout.append(list1[i])
# 				else:
# 					listout.append("        "+list1[i])
# 		else:
# 			listout = listout+ token2Line(tokennewy)
# 			listkeyfeature[ listkey.index(y) ] = ""

# 	for x in listkeyfeature:
# 		if x<>"":
# #			list1=textwrap.wrap(x)
# 			list1=textwrap_wrap(x)
# 			for i in range(0,len(list1)):
# 				if i==0:
# 					listout.append(list1[i])
# 				else:
# 					listout.append(7*" "+list1[i])
		
# 	return listout	
		
# def textwrap_wrap(x): #Word Wrap at 70 colmuns	
# 	nwrap=70
# 	ic=0
# 	n =0
# 	b =[]
# 	for i in range(0,len(x)):
# #		print '3333  ', i, len(x)
# 		if(x[i]==' '): ise=i
# 		n=n+1
# 		if(n>nwrap):
# 			b.append( x[ic:ise] )
# 			ic= ise
# 			n=0
# 		elif(i==len(x)-1):
# 			b.append( x[ic:len(x)] )
# #			print '222 ', x[ic:len(x)]
# #	print b
# 	return b


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
try:
	ext=sys.argv[1]
	print "Generate ctrl."+ext+" from ctrls."+ext + "..."
except:
	print "--- This is a standard setting when no SPEC is specified"
	print   specstd
	print "--- SPEC shown above is standard setting when no SPEC is specified"
	print 
	print " Purpose: Generate ctrl.{ext} file from ctrls.{ext}"
	print
	print " Usage  : ctrlgen {extension of ctrl file}"
	print
	sys.exit()
ctrls = "ctrls." + ext


#### Read in ctrls #####
f=open(ctrls,'rt')
ctrlsdat = f.read() 
f.close()

listctrls  = lineReadfile(ctrls) 
listspec   = GetCat(listctrls,"SPEC")  # SPEC section only
if(len(listspec)==0):
	listspec = re.split('\n',specstd)  # SPEC standard if no SPEC is in ctrls.*
	print 
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
	specdic[xx[0]]= '  ATOM='+ii+'\n'
#print specdic
#print sitename,'yyy',uniq(sitename)
	
specsec=''
for i in uniq(sitename):
	specsec= specsec + specdic[i]
#print specsec

##################
os.system("rm -rf llmchk_getwsr llmfa.tmp2")
head= """
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
%const kmxa=7  # kmxa=7 is good for pwemax=5 or lower.
               # larger kmxa is better but time-consuming (maybe not the critical part for large systems).\n
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
os.system("lmchk --getwsr tmp > llmchk_getwsr; echo $? >exitcode")
f=open("exitcode",'rt')
iexit=int(f.read())
f.close()

if(iexit==0):
	print ' -----tail of llmchk_getwsr ----------------------------'
	os.system("tail llmchk_getwsr")
	print ' --- lmchk --getwsr has done successfully! rmt.tmp is generated -----'
else:
	print 
	print ' lmchk --getwsr can not find the muffin-tin radius SPEC_ATOM_R.'
	print ' Wrong ctrls? or bug in ctrlgen.py? (see llmchk_getwsr)'
        print '   or you have to write rmt.tmp by yourself.'
	print ' rmt.tmp consists of "specname R" for each line. For example,---'
        print '       --- rmt.tmp for SrTiO3 --- '
	print '       Sr          3.616323'
	print '       Ti          2.089960'
	print '       O           1.595007'
        print '       --- end of rmt.tmp ------- '
        print ' After you write rmt.tmp, repeat ctrlgen.py'
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
print 
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
for ii in tokenspec[1:]:
	if(ix==1):
		ikey=ii
		ix=0
	elif(ii=='ATOM='):
		if(ispec>0): 
			mmmx=mtodic[ikey]
			mmm=re.sub(","," ",mmmx)
			il1 = countnum(mmm,'RSMH=')
			il2 = countnum(mmm,'PZ=')
			lx = max(il1,il2)+1
			lll = "%i" % lx
			#print il1,il2,lx
			aaa = aaa+' R='+rdic[ikey]+'\n'+' '*5+mmm+'\n' \
			    +' '*6+'KMXA={kmxa} LMXA='+lll+'\n'+'     MMOM=0 0 0 0'+'\n\n'
		ispec=ispec+1
		ix=1
	aaa=aaa+' '+ii
mmmx = mtodic[ikey]
mmm = re.sub(","," ",mmmx)
il1 = countnum(mmm,'RSMH=')
il2 = countnum(mmm,'PZ=')
lx = max(il1,il2,3)+1
lll = "%i" % lx
#print il1,il2,lx
aaa = aaa+' R='+rdic[ikey]+'\n'+' '*5+mmm+'\n' \
    + ' '*6+'KMXA={kmxa} LMXA='+lll+'\n'+' '*6+'MMOM=0 0 0'

tail="""
\n
% const pwemax=3 nk=2
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

ITER  CONV=1e-6 CONVC=1e-6 NIT=30
                # An other choice is
                # ITER MIX=A2,b=.5,n=3 CONV=1e-6 CONVC=1e-6 NIT=20
                # Practically results are independent from mixing procedure.
		
HAM   NSPIN=1   # Set NSPIN=2 for spin-polarize case; then set SPEC_MMOM (initial guess of magnetic polarization).
      FORCES=0  # 0: no force calculation, 1: forces calculaiton 
      GMAX=9    # this is for real space mesh. See GetStarted.
      REL=T     # T:Scaler relativistic, F:non rela.

      XCFUN=1   # =1 for VWN; GGA is not yet.
                # XCFUN=2 shows a bug for Hydrogen atom. 
		# (subs/evxc.F works only for XCFUN=1 if rho(up)=0 or rho(down)=0).

      PWMODE=11 # 10: MTO basis only (LMTO) PW basis is not used.
                # 11: APW+MTO        (PMT)
                # 12: APW basis only (LAPW) MTO basis is not used.

      PWEMAX={pwemax} # (in Ry). When you use larger pwemax more than 5, be careful
                      # about overcompleteness. See GetStarted.
      ELIND=-1  # this is only for accelaration of convergence. Not need to change.
      
OPTIONS PFLOAT=1 # Q=band (this is quit switch if you like to add)
                 # 
"""

### Write ctrl.ext
g = open("ctrl."+ext,'wt')
g.write(ctrlnospec+aaa+tail)
g.close()
print " Check ctrl."+ext
sys.exit()

