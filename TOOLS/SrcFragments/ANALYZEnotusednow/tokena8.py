#   By Hiori Kino, Aug, 2004
#  modified: Sep, 24, 2004 by H.Kino
#
import sys
import string
#import textwrap
# version 2.3 or higher is necessary to import textwrap

listctrlkey=["HEADER","VERS","IO","CONST","STRUC","OPTIONS","BZ","SYMGRP","SPEC","SITE",
                "EWALT","MIX","HAM","START","GW"]


#----------------------------------------------------------------------------
def token2Line(tokena):
	""" connect a tokana list and return strings to write to the file """
	keylist=tokenMakeKeylist(tokena)
	n=len(tokena)
	listout=[]
	str1 = "%-8s" % (tokena[0])
	newone=0
	for i in range(1,n):
		str1 = str1 + tokena[i]+" "
		newone=1
		if len(str1)> 50:
			iskey=0
			for x in keylist:
				if x==tokena[i]:
					iskey=1
			if iskey==0:
				listout.append(str1)
				str1="        "
				newone=0
	if newone==1:
		listout.append(str1)
	return listout
		


#------------------------------------------------------------------------------
def tokenExpandMtoq(atom,specatom,listmtofp):
	""" replace MTOQ secction with the content of listmtofp """
	speckey=tokenMakeKeylist(specatom)
	mtoqlist=tokenExtract(specatom,"MTOQ=","",speckey)
        specatom=tokenRemove(specatom,mtoqlist[0])
	mtoqlist.remove("MTOQ=")
	speclist=lineExtract(listmtofp,"SPEC",atom,[])
	l=0
	casedata=[]
	for x in mtoqlist:
		caselist = lineExtract(speclist,"CASE",x,[])
		if caselist==[]:
			print
			print "Error: Invalid CASE ",x, " is not found at MTOQ= section of ATOM=",atom
			print
			sys.exit(10)
		casedata.append( caseExtract(caselist,l) )
		l=l+1
	
	return tokenArrangeMtoq(casedata)
	
	
#------------------------------------------------------------------------------
def tokenArrangeMtoq(casedata):
	""" make EH= -1.0 -1.0 0.0,,, RSMH= 2.0 2.4 1.0 ... """
	""" from RH=1.0 RSMH=2.0.. EH=-1.0 RSMH=2.4 ... """
	flageh=0
	flagrsmh=0
	flageh2=0
	flagrsmh2=0
	flagpz=0
	flagp=0
	flagidmod=0
	for x in casedata:
		for st in x:
			if st=="EH=":	
				flageh=1
			if st=="RSMH=":
				flagrsmh=1
			if st=="EH2=":
				flageh2=1
			if st=="RSMH2=":
				flagrsmh2=1
			if st=="P=":
				flagp=1
			if st=="PZ=":
				flagpz=1
			if st=="IDMOD=":
				flagidmod=1
			if x==[]:
				break

	lmx = len(casedata)-1
	

	list1=[]
	defaultvalueeh="0.00"
	defaultvaluersmh="-1.00"
	defaultvaluep="0.00"
	defaultvaluepz="0.00"
	defaultvalueidmod=["0","0","0","1","1","1","1"]
	
# EH
	if flageh==1:
		list1.append("EH=")
		for x in casedata:
			value=tokenKeyvalue(x,"EH=")	
			if value<>[]:
				list1.append(value[0])
			else:
				list1.append(defaultvalueeh)
		str= defaultvalueeh
		while str==defaultvalueeh:
			str=list1.pop()
		list1.append(str)
		
# RSMH
	if flagrsmh==1:
		list1.append("RSMH=")
		for x in casedata:
			value=tokenKeyvalue(x,"RSMH=")	
			if value<>[]:
				list1.append(value[0])
			else:
				list1.append(defaultvaluersmh)
		str=defaultvaluersmh
		while str==defaultvaluersmh:
			str=list1.pop()
		list1.append(str)
		
#EH2
	if flageh2==1:
		list1.append("EH2=")
		for x in casedata:
			value=tokenKeyvalue(x,"EH2=")	
			if value<>[]:
				list1.append(value[0])
			else:
				list1.append(defaultvalueeh)
		str=defaultvalueeh
		while str==defaultvalueeh:
			str=list1.pop()
		list1.append(str)

# RSMH
	if flagrsmh2==1:
		list1.append("RSMH2=")
		for x in casedata:
			value=tokenKeyvalue(x,"RSMH2=")	
			if value<>[]:
				list1.append(value[0])
			else:
				list1.append(defaultvaluersmh)
		str=defaultvaluersmh
		while str==defaultvaluersmh:
			str=list1.pop()
		list1.append(str)
# P
	if flagp==1:
		list1.append("P=")
		for x in casedata:
			value=tokenKeyvalue(x,"P=")	
			if value<>[]:
				list1.append(value[0])
			else:
				list1.append(defaultvaluep)

# IDMOD
	if flagidmod==1:
		list1.append("IDMOD=")
		l=0
		for x in casedata:
			value=tokenKeyvalue(x,"IDMOD=")	
			l=l+1
			if value<>[]:
				list1.append(value[0])
			else:
				list1.append(defaultvalueidmod[l])

#PZ
	if flagpz==1:
		list1.append("PZ=")
		for x in casedata:
			value=tokenKeyvalue(x,"PZ=")
			if value<>[]:
				list1.append(value[0])
			else:
				list1.append(defaultvaluepz)


	
	
	return list1

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
def tokenRemove(list1,name):
	""" remove a block starting name , list1 is changed, nothing returns"""
	listkey=tokenMakeKeylist(list1)
	listname=tokenExtract(list1,name,"",listkey)
	n0 = len(list1)
	pos=-1
	for i in range (0,n0):
		if string.upper(list1[i])==string.upper(name):
			pos=i
			break
	if pos>=0:
		for i in range(0,len(listname)):
			list1.pop(pos)

	
		
#------------------------------------------------------------------------------
def tokenPrintSpec(list1):
	"""  print token list1  which is spec section of a ctrl file """
	NLlist=["EH=","RSMH=","EH2=", "RSMH2=", "IDMOD=", "P=","PZ="]	
	clearlist=["ATOM="]
	nlflag=0
	for x in list1:
		if x=="":
			continue
		for y in NLlist:
			if  y == x:
				nlflag=1
				sys.stdout.write("\n           ")
		for y in clearlist:
			if y==x:
				nlflag=0
				sys.stdout.write("        ")
				

		if nlflag==1:
			sys.stdout.write("%6s" % (x))
			sys.stdout.write(" ")
		else:
			sys.stdout.write(x)
                        sys.stdout.write(" ")


	print


#------------------------------------------------------------------------------
def tokenExpandQuality(specatom,listmtofp):
	""" change specatom, replace quality=??? with the result of listmtofp """
	qualitylist=tokenKeyvalue(specatom,"quality=")
	qualityvalue=""
	if qualitylist<>[]:
		qualityvalue=qualitylist.pop()
		qualitydetail=findQuality(listmtofp, atom,qualityvalue)
		tokenRemove(specatom,"quality=")
		specatom=tokenMergelist(specatom,qualitydetail)
	return specatom


#------------------------------------------------------------------------------
def tokenKeyvalue(token,key):
	""" get token followed by key """
	compare=0
	listout=[]
	for x in token:
		if compare==1:
			listout.append(x)
			compare=0
		if string.upper(x)==string.upper(key):
			compare=1

	return listout



#---------------------------------------------------

def tokenMakeKeylist(tokena):
	""" list up keywords ending at = """
	listout=[]
	for x in tokena:
		n = len(x)
		if x[n-1]=="=" or x[n-1]==":":
			listout.append(x)
	return listout

#---------------------------------------------------------
def tokenMergelist(tokena,tokenb):
	""" merge a tokena list  and a tokenb list, if the key of the tokanb list does not exist in the tokena list """
	keya=tokenMakeKeylist(tokena)
	keyb=tokenMakeKeylist(tokenb)
	for a in keya:
		for b in keyb:
			if string.upper(a)==string.upper(b):
				listb=tokenExtract(tokenb, b,"", keyb)
				tokenRemove(tokenb,listb[0])	

	return tokena+tokenb


#----------------------------------------------------
def tokenExtract(tokenin,name,key,namelist):
	""" extract token staring name following key and ending one of namelist or name """
	""" key can be null and namelist can be null list """
#	""" input=list output=list """
	listout=[]
	start=0
	compare=0
	s1=""
	for s in tokenin:
		if compare==1:
			if string.upper(s)==string.upper(key):
				start=1
				listout.append(s1)
				listout.append(s)
				compare=0
				continue
		if string.upper(s)==string.upper(name):
			if start==1:
				return listout
			else:
				if len(key)>0:
					compare=1	
					s1=s
					continue
				else:
					start=1
					listout.append(s)
		else:
			if start==1 and namelist<>[]:
				for x in namelist:
					if string.upper(x)==string.upper(s):
						return listout
			if start==1:
				listout.append(s)
	return listout


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

def lineExtract(linein,name,key,namelist):
	""" extract list staring name following key and ending name of namelist """
	""" key can be null and namelist can be null """
#	""" input=readlines() output=lines"""
	listout=[]
	start=0
	for s in linein:
		ret=string.split(s)
		if ret==[]:
			continue
		st=string.upper(ret.pop(0))
		if st==string.upper(name):
			if start==1:
				return listout
			else:
				if len(key)>0:
					st = string.upper(ret.pop(0))
					if st==string.upper(key):
						start=1
						listout.append(s)	
				else:
					start=1
					listout.append(s)
		else:
			if start==1 and namelist<>[]:
				for x in namelist:
					if st==string.upper(x):
						return listout
			if start==1:
				listout.append(s)

	return listout

#------------------------------------------------------------------------------

def listCopy(listin):
	""" copy list """
	listout=[]
	for x in listin:
		listout.append(x)

	return listout


			
#-----------------------------------------------------------
def findQuality(linemtofp, atom,qualityvalue):
	""" get qualist section of atom , quality = qualityvalue """
	listspec=lineExtract(linemtofp,"spec", atom,[])
	listquality=lineExtract(listspec,"quality",qualityvalue,["quality","case"])
	if listquality==[]:
		print
		print "Error: quality=",qualityvalue, " is not found in ATOM=",atom 
		print 
		sys.exit(10)
	tokenquality=line2Token(listquality)
	listout=[]
	for x in tokenquality:
		if string.upper(x)<>"QUALITY" and string.upper(x)<>string.upper(qualityvalue):
			listout.append(x)

	return listout


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
def caseExtract(caselist,l):
	" make a list such as EH= -1.0 RSMH= 3.0 EH2=... , for l"
	listout=[]
	casetoken=line2Token(caselist)
	casekey=tokenMakeKeylist(casetoken)
	for token in casekey:
		token1=tokenExtract(casetoken,token,"",casekey)
		if token1<>[]:
			token1.remove(token)
			l0=0
			value=""
			for x in token1:
				if l0==l:
					value=x
					break
				l0=l0+1
			listout.append(token)
			listout.append(value)
	return listout


#-------------------------------------------------------------------------
def LineRemoveSection(listctrl,key,namelist):
	""" remove key section from listctrl """
	start=-1
	end=len(listctrl)
	i=-1
	for x in listctrl:
		i=i+1
		x1=[x]
		tok = line2Token(x1)
		if tok==[]:
			continue
		if start==-1 and string.upper(tok[0])== string.upper(key):
			start=i	
			continue
		if start>=0:
			for y in namelist:
				if string.upper(y)==string.upper(tok[0]):
					end=i
					break
		if start>=0 and end>=0:
			break
	if start>=0 and end>=0:
		for i in range(start,end):
			listctrl.pop(start)


#---------------------------------------------------------------------------
def  lineWriteCtrl(listctrl,listspecatom):
	""" write listctrl and listspecatom """
	for x in listctrl:
		print x
	print
	print "SPEC"
	for x in listspecatom:
		tokenPrintSpec(x)
	print
	

#--------------------------------------------------------------------------
def lineAddCtrlFeature(listctrl):
	""" add keyword lists to ctrl """
	listkeyfeature=[
		"IO      SHOW=F HELP=F VERBOS=30 ",
		"OPTIONS NSPIN=2 REL=T XCFUN=2",
		"BZ      NKABC=2 2 2 TETRA=1 BZJOB=0 TETRA=1 METAL=3 ",
		"EWALD   TOL=1D-8 NKDMX=1999 NKRMX=1999",
		"MIX     CONV=1d-4 CONVC=1d-4",
		"HAM     ELIND=-.7 TOL=1d-6  GMAX=12.5                   "+
		"        RDSIG=12 SIGP:3,0,0,0,2.,0,.06,0",
		"START   NIT=100",
		"GW      NKABC= GCUTB= GCUTX="]
	listkey=["IO","OPTIONS","BZ","EWALD","MIX","HAM","START","GW"]

	listout=[]
	for x in listctrl:
		t = line2Token([x])
		if t==[]:
			listout.append(x)
			continue
		found=0
		for y in listkey:
			if string.upper(y)==string.upper(t[0]):
				listy=lineExtract(listctrl,y,"",listctrlkey)
				tokeny=line2Token(listy)
				listadd=listkeyfeature[ listkey.index(y) ]
				tokenadd=line2Token([listadd])
				tokeny.pop(0)
				tokenadd.pop(0)
				tokennewy= tokenMergelist(tokeny,tokenadd)
				tokennewy.insert(0,y)
				found=1
				break
		if found==0:
#			list1= textwrap.wrap(x)
			list1= textwrap_wrap(x)
			for i in range(0,len(list1)):
				if i==0:
					listout.append(list1[i])
				else:
					listout.append("        "+list1[i])
		else:
			listout = listout+ token2Line(tokennewy)
			listkeyfeature[ listkey.index(y) ] = ""

	for x in listkeyfeature:
		if x<>"":
#			list1=textwrap.wrap(x)
			list1=textwrap_wrap(x)
			for i in range(0,len(list1)):
				if i==0:
					listout.append(list1[i])
				else:
					listout.append(7*" "+list1[i])
		
	return listout	
		
def textwrap_wrap(x): #Word Wrap at 70 colmuns	
	nwrap=70
	ic=0
	n =0
	b =[]
	for i in range(0,len(x)):
#		print '3333  ', i, len(x)
		if(x[i]==' '): ise=i
		n=n+1
		if(n>nwrap):
			b.append( x[ic:ise] )
			ic= ise
			n=0
		elif(i==len(x)-1):
			b.append( x[ic:len(x)] )
#			print '222 ', x[ic:len(x)]
#	print b
	return b

#----------------------------------------------------------------------------



#===============================================================================
#          MAIN 
#===============================================================================




ext=sys.argv[1]
print "% made by  python", sys.argv[0], sys.argv[1]

listmtofp = lineReadfile("mtofp."+ext)

#---

listctrl = lineReadfile("ctrl."+ext) #ctrl file 


### to get nspec from ctrl
liststruc=lineExtract(listctrl,"struc","",listctrlkey)
#  STRUCT section of the ctrl file

tokenstruc=line2Token(liststruc)
structkey=tokenMakeKeylist(tokenstruc)
listnspec=tokenExtract(tokenstruc,"nspec=","",structkey)
nspec = listnspec[1]


### to get SPEC section of ctrl
listspec=lineExtract(listctrl,"spec","",listctrlkey)
tokenspec=line2Token(listspec)
#   tokenspec= tokens of SPEC section of the ctrl file


atomlist=tokenKeyvalue(tokenspec,"ATOM=")
# e.g.  atomlist=["Mn", "La", "O"] 

listspecatom=[]

for atom in atomlist:
	# to extract ATOM= atom in SPEC section of ctrl 
	specatom=tokenExtract(tokenspec,"atom=",atom,[])
#        tokenRemove(specatom,"LMX=")

	specatom=tokenExpandQuality(specatom,listmtofp)	
	st      =tokenExpandMtoq(atom,specatom,listmtofp)
	specatom=specatom+st
	listspecatom.append(specatom)


LineRemoveSection(listctrl,"spec", listctrlkey)
	

listctrl2=lineAddCtrlFeature(listctrl)

lineWriteCtrl(listctrl2,listspecatom)


