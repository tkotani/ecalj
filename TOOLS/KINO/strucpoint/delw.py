#/bin/python
import os
import re
import sys
from types import *
import string

import mFLine
import mStrucList

from mToken import *

spc10='          '
spc30=spc10+spc10+spc10
#spc6='123456'
spc6='      '


thisprogram='Cdelw2'
strinfo='...info...'+spc10

class t_struc_elem:
	struc=''
	elem=''
	target=''
	source=''
	funcname=''
	def __init__(self,ll):
		struc,elem,funcname,target,source=ll
		self.struc=struc
		self.elem=elem
		self.funcname=funcname
		self.target=target
		self.source=source	
	def __str__(self):
		p=' '
		s1=list2str(self.target)
		if len(s1)==0:
			s1='(None)'
		s2=list2str(self.source)
		if len(s2)==0:
			s2='(None)'
		name=self.struc+'%'+self.elem
		s=name.ljust(12)+" target="+s1.ljust(12)+p+"source="+s2.ljust(12)+p+self.funcname
		return s

def struc_elem_cmp(x,y):
    if x.funcname>y.funcname:
	return 1
    elif x.funcname<y.funcname:
	return -1
    else:
	if x.struc>y.struc:
		return 1
	elif x.struc<y.struc:
                return -1
	else:
		if x.elem > y.elem:
			return 1
		elif x.elem==y.elem:
			return 0
		else:
			return -1

class t_reason:
	reason=[]
	name=''
	funcname=''
	def __init__(self,ll):
		reason,name,funcname=ll
		self.name=name
		self.reason=reason
		self.funcname=funcname
	def __str__(self,p=' '):
		return self.name+p+list2string(self.reason)+p+self.funcname
	def makelist(self):
		return [self.reason,self.name,self.funcname]

class t_var7:
        name=''
        nameorg=''
        type=''
        indent=0
        size=''
        prevar=''
        funcname=''
        def __init__(self,ll):
		name,indent,nameorg,size,type,prevar,funcname =ll
                self.name=name
		self.nameorg=nameorg
		self.indent=indent
		self.size=size
                self.type=type
		self.prevar=prevar
                self.funcname=funcname
        def __str__(self):
		p=' '
		name=self.name
		nameorg=self.nameorg
		type=self.type
		prevar=self.prevar
                s=name+p+str(self.indent)+p+nameorg+p+str(self.size)+p+type+p+prevar+p+self.funcname
		return s
	def makelist(self):
		return [ self.name,self.indent,self.nameorg,self.size,self.type,self.prevar,self.funcname]

def t_var7_cmp(x,y):
	if x.funcname>y.funcname:
		return 1
	elif x.funcname<y.funcname:
                return -1
	else:
		if x.name > y.name:
			return 1
		elif x.name < y.name:
                        return -1
		else:
			return 0


class t_arraydef:
	name=''
	arraysize=[]
	def __init__(self,name,arraysize=[]):
		self.name=name
		self.arraysize=arraysize
	def __str__(self,p=','):
		s=name+p+"["+list2string(arraysize,p)+"]"
		return s


def insert_tok(orig,rep,change):
	thisfunc="insert_tok"
	if orig[0]=='if':
		if 'then' in orig:
			print thisfunc,"uknown grammar", orig,rep,change
			print >>sys.stderr, thisfunc,"uknown grammar", orig,rep,change
			sys.exit(10)
		else:
			#search rep
			for i in range(len(orig)):
				if orig[i]==rep:
					break
			newtok=[]
			for j in range(len(orig)):
				if j==i:
					break
				newtok.append(orig[j])
			newtok.append('then')
			for nn in change:
				newtok.append(nn)
			if rep=='return':
				newtok.append([rep])
			newtok.append('endif')
			if i!=len(orig)-1:
				print thisfunc,"uknown grammar", orig,rep,change
				print >>sys.stderr, thisfunc,"uknown grammar", orig,rep,change
	                        sys.exit(10)
			return newtok
	elif orig[0]=='return':
		newtok=[]
		for nn in change:
			newtok.append(nn)
		newtok.append('return')
		return newtok

	print thisfunc,"uknown grammar", orig,rep,change
	print >>sys.stderr, thisfunc,"uknown grammar", orig,rep,change
	sys.exit(10)
				
					


class Keyword:
    list = []
    def __init__(self):
        self.list.append("end")
        self.list.append("if")
        self.list.append("then")
	self.list.append("else")
	self.list.append("elseif")
	self.list.append("endif")
        self.list.append("continue")
        self.list.append("do")
        self.list.append("enddo")
	self.list.append("defrr")
	self.list.append("defr")
	self.list.append("defdr")
	self.list.append("defi")
	self.list.append("defdc")
	self.list.append("defcc")
	self.list.append("redfrr")
	self.list.append("redfi")
	self.list.append("rlse")
	self.list.append("return")
	self.list.append("go")
	self.list.append("goto")
	self.list.append("select")

#---------------------------------------------------------

class DoParse:
    printit=0

    printsource=0

    execfuncname=['','']

    var_dic={};

    newvar=[]
    newvarfile='__newvar.dat'
    newvarfile_desc=0

    keyword=[]
    struclist=[]
    uselist=[]
    use_toadd=[]

    indentstack=[]
    w_varlist=[]
    w_varlistall=[]
    w_varlistundel=[]

    w_varlistundel_fix=0

    printmap=0
    mapfile='__file.map'

    callvarlist=[]

    idisable_change=0

    struc_table=[]
    w_pointerlist=[]

    struc_elem=[]

    def __init__(self,varfile0,mapfile0,printsource0,printmap0):
	thisfunc='__init__'
	self.keyword=Keyword()
        self.struclist=mStrucList.StrucList()

	self.newvarfile=varfile0
	self.mapfile=mapfile0
    	self.printsource=printsource0
    	self.printmap=printmap0

        self.read_struc_table()
        self.w_pointerlist=[]
	self.w_varlistundel_fix=0

	struc_elem=[]

	if printmap0==1:
	    #print thisfunc,"clear", self.mapfile
	#clear self.
	    try:
		fo=open(self.mapfile,'w')
	    except:
		print 'ERROR: failed to open ',thisfunc, self.mapfile
		print>>sys.stderr, 'ERROR: failed to open ',thisfunc, self.mapfile
		sys.exit(10)
	    fo.close()

	self.read_newvar()

	try:
		sys.unlink(self.newvarfile)
	except:
		idummy=0

	# reopen for writing
        try:
                self.newvarfile_desc=open(self.newvarfile,'w')
        except:
		print 'ERROR: failed to open ',thisfunc, self.newvarfile,self.newvarfile_desc
		print >> sys.stderr, 'ERROR: failed to open ',thisfunc, self.newvarfile,self.newvarfile_desc
		sys.exit(10)

	#print "FILE,open",self.newvarfile,thisfunc
	

    def __del__(self):
	    thisfunc="__del__"
	    #print "FILE,close",self.newvarfile,thisfunc
	    try:
		self.newvarfile_desc.close()
	    except:
		print 'failed to close for writing',self.newvarfile
		print>>sys.stderr, 'failed to close for writing',self.newvarfile
		sys.exit(10)

    def read_struc_table(self):
	thisfunc="read_struc_table"
	name="struc_types"
	try:
		fr=open(name,"r")
	except:
		print thisfunc,"failed to open",name
		sys.exit(10)
        #print "FILE,open",name,thisfunc
	self.struc_table=[]
        for line in fr:
                aline= line[:-1].split(' ')
		if len(aline)==0 or aline[0][0]=="#":
			continue
		type = aline[2]
		if len(aline)>4:
			if aline[4]=='unchange':
				type='unchange'
		self.struc_table.append([aline[0],aline[1],type])

	fr.close()
	#print "FILE,close",name,thisfunc

    def struc_element_type(self,struc,element):
	for nn in self.struc_table:
		if nn[0]==struc and nn[1]==element:
			return nn[2]
	print thisfunc,"ERROR failed to find struc,element",struc,element
	sys.exit(10)

    def save_wvar(self,name):
	thisfunc='save_wvar'
	try:
		fw=open('#delw.wlist.'+name,'w')
	except:
		print thisfunc,'failed to open',name
		sys.exit(10)

	for nn0 in self.w_varlist:
		fw.write(nn0+'\n')

	fw.close()
	#print "FILE,close",'#delw.wlist.'+name,thisfunc
	if self.printsource==1:
		print thisprogram,"save w_varlist=",
		self.w_varlist_show()

    def load_wvar(self,name):
        thisfunc='restore_wvar'
        try:
                fr=open('#delw.wlist.'+name,'r')
        except:
                print thisfunc,'failed to open',name
                sys.exit(10)
	#print "FILE,open",'#delw.wlist.'+name,thisfunc

	self.w_varlist=[]
	for line in fr:
		aline= line[:-1].split(' ')
		if len(aline)<7:
			print thisfunc,"ERROR, field!=7"
			print line
			sys.exit(10)	
		aline=aline[0:7]
		for i in range(len(aline)):
			if aline[i]=='***':
				aline[i]=''
		aline[1]=int(aline[1])
		self.w_varlist.append(t_var7(aline[0:7]))
		
	fr.close()
	#print "FILE,close",'#delw.wlist.'+name,thisfunc

	if self.printsource==1:
		print thisprogram,"load w_varlist=",
		self.w_varlist_show()


    def read_newvar(self):
# store in var_toadd
	thisfunc='read_newvar'
        try:
                fr=open(self.newvarfile,'r')
        except:
		return 
	#print "FILE,open",self.newvarfile,thisfunc

	self.w_varlistall=[]
	self.w_varlistundel=[]
	for line in fr:
		aline= line[:-1].split(' ')
		#print thisfunc, aline
		if aline[0]=='wvar':
			var, name,level,nameorg,namesize,type,prevar,execfuncname = aline[0:8]
		elif aline[0]=='change':
			var,ok,execfuncname = aline[:3]
		elif aline[0]=='undel':
			#var,reason,name,execfuncname= aline[:4]
			var,name,execfuncname=aline[:3]
			listreason=aline[3:]
		elif aline[0]=='strucelem_t' or aline[0]=='strucelem_s':
			var,struc,elem,execfuncname=aline[:4]
			listreason=aline[4:]
		elif aline[0]=='pointer':
			var, name,level,nameorg,namesize,type,prevar,execfuncname = aline[0:8]
		else:
			print thisfunc,"ERROR ",aline
			print >>sys.stderr, thisfunc,"ERROR ",aline
			sys.exit(10)

		if var=='wvar':
			self.w_varlistall_append_uniq(t_var7([name,level,nameorg,namesize,type,prevar,execfuncname]),0,[])
		elif var=='undel':
    			self.w_varlistundel_append_uniq(t_reason([listreason,name,execfuncname]),0,[])
		elif var=='strucelem_t':
			self.struc_elem_append_uniq(t_struc_elem([struc,elem,execfuncname,listreason,'']))
		elif var=='strucelem_s':
			self.struc_elem_append_uniq(t_struc_elem([struc,elem,execfuncname,'',listreason]))
		elif var=='pointer':
			self.w_pointerlist_append_uniq(t_var7([name,level,nameorg,namesize,type,prevar,execfuncname]))

	fr.close()
	#except:
	#	print thisfunc,"failed to close for reading, ",self.newvarfile
	#	sys.exit(10)
	#print "FILE,close",self.newvarfile,thisfunc

	# delete undel if it is in pointer
	ltmp=[]
	for nn in self.w_varlistundel:
		wp=  self.w_pointerlist_find(nn.name,nn.funcname)
		if isinstance(wp,t_var7)==1:
			reason=nn.reason
			if wp.type=='undef':
				ltmp.append(nn)
			elif 'call' in reason or 'pass' in reason or 'o_eq' in reason or 'eq' in reason:
				ltmp.append(nn)
		else:
			ltmp.append(nn)
	self.w_varlistundel=ltmp	

        self.struc_elem.sort(struc_elem_cmp)

	if self.printsource==1:
	    vlist=[]
	    for nn in self.w_varlistall:
		vlist.append(nn.name)
	    for nn in self.w_varlistundel:
		vlist.append(nn.name)
	    for nn in self.w_pointerlist:
		vlist.append(nn.name)
	    vlist=list(set(vlist))
	    for name in vlist:
	      for nn in self.w_varlistall:
		if name==nn.name:	
			msg="wvar".ljust(8)
			print thisprogram,msg,nn
	      for nn in self.w_varlistundel:
		if name==nn.name:	
			msg="undel".ljust(8)
			print thisprogram,msg,nn
	      for nn in self.w_pointerlist:
		if name==nn.name:	
			msg="pointer".ljust(8)
			print thisprogram,msg,nn

	    for nn in self.struc_elem:
		print thisprogram,"struc_elem",nn


    def write_newvar(self):
# write newvar
	thisfunc='write_newvar'
	fw=self.newvarfile_desc
	if fw==0:
		print thisfunc,'ERROR newvarfile is not opened'
		print >>sys.stderr, thisfunc,'ERROR newvarfile is not opened'
		sys.exit(10)


	for aline in self.w_varlistall:
		#name,level,nameorg,namesize,type,prevar,execfuncname = aline
	    if aline.name in self.callvarlist:
		name=aline.name
		execfuncname=aline.funcname
		fw.write('undel '+name+' '+execfuncname+' call\n')
	    else:
		linew='wvar '+ aline.__str__()
		fw.write(linew+'\n')

	for aline in self.w_varlistundel:
		linew='undel '
		reason,name,execfuncname = aline.makelist()
		if name in self.callvarlist:
			reason.append('call')
		reason=list(set(reason))
		#fw.write(linew+reason+' '+name+' '+execfuncname+'\n')
		fw.write(linew+name+' '+execfuncname+' '+list2string(reason,'','')+'\n')

        if self.idisable_change==1:
		ok='disable'
	else:
		ok='enable'
	#print 'change '+ok+ ' '+self.execfuncname[0]
	fw.write('change '+ok+ ' '+self.execfuncname[0] + '\n')

	for nn in self.w_pointerlist:
		fw.write('pointer '+nn.__str__()+'\n')
	p=' '
	for nn in self.struc_elem:
		fw.write('strucelem_t '+nn.struc+p+nn.elem+p+nn.funcname+p+list2str(nn.target)+'\n')
	for nn in self.struc_elem:
		fw.write('strucelem_s '+nn.struc+p+nn.elem+p+nn.funcname+p+list2str(nn.source)+'\n')


    def do_stringoption__is_ol(self,cmd0):
# cmd=list
# 'a b:15' -> 'a b' and and range 
#input string: cmd 
	thisfunc="do_stringoption__is_ol"
	tok=[]
	cmd=cmd0.split(' ')
	for name in cmd:
		if len(name)==0: 
			continue
		namel=[]
		namel=name.split(';')
		namel=do_split0__ils_ol(namel,':')
		namel=do_split0__ils_ol(namel,'.')
		namel=do_split0__ils_ol(namel,',')
		tok.append(namel)
	tok2=[]
	for elem in tok:
		name=elem[0]
		alias=''
		mask=''
		range1=''
		range2=''
		n=len(elem)
		i=0
		while i<n:
			a=divmod(i,2)
			if elem[i]==';':
				if a[1]==0:
					print 'ERROR',thisfunc,'delm=; ', elem
				i=i+1
				alias=elem[i]
			if elem[i]==',':
				if a[1]==0:
					print 'ERROR',thisfunc,'delm=, ', elem
				i=i+1
				mask=elem[i]
			if elem[i]==':':
				if a[1]==0:
					print 'ERROR',thisfunc,'delm=: ', elem
				i=i+1
				range1=elem[i]
			if elem[i]=='.':
				if [1]==0:
					print 'ERROR',thisfunc,'delm=. ', elem
				i=i+1
				range2=elem[i]
			i=i+1
		if range1!='' and range2=='':
			range2=range1
		tok2.append([name,alias,mask,range1,range2])
	return tok2


    def add_struc(self,strstruc,struc):
        alist=[self.execfuncname[0],strstruc,struc]
	for nn in self.newvar:
		if alist==nn:
			return
	self.newvar.append(alist)

	
    def merge_variableinline__il_ol(self,tok):
	tok2=[]
	lasttype='sep'
	for name in tok:
		if name=='(' or name==',' or name==')':
			currenttype='sep'		
		else:
			currenttype=''
		if lasttype==currenttype and lasttype=='':
			a=tok2.pop()
			name=a+name
		tok2.append(name)
		lasttype=currenttype
	return tok2

    def get_nextindex(a):
	n=len(a)
	return n


    def do_nothing(self):
	i=1

    def w_varlist_release(self,name,level,linenumber,printit=1):
	thisfunc="w_varlist_release"
	if printit==1:
		print thisprogram,"rlse name=",name,"old_list=",
		self.w_varlist_show()

	if len(self.w_varlist)==0:

	    	if printit==1:
			print ''
			print '#error ERROR, try to release name=',name, ',but list=null at linenumber=',linenumber,
			print 'list=',
			self.w_varlist_show()
                        print >>sys.stderr,'ERROR, try to release name=',name, ',but list=null, at linenumber=',linenumber

		self.idisable_change=1
		return

	idel=-1
	for i in range(len(self.w_varlist)):
		a=self.w_varlist[i].name
		if a == name:
			idel=i
			break

	list_rel=[]
	if idel>=0:
		list_rel=self.w_varlist[idel:]
		tok2= []
		for i in range(idel):
			tok2.append(self.w_varlist[i])
		self.w_varlist=tok2	
	else:
		if printit==1:
			print ''
			print '#error ERROR, try to release name=',name, ',but list does not have ',name,'at linenumber=',linenumber,
			print "list=",
			self.w_varlist_show()
			print ''
			print >>sys.stderr, 'ERROR, try to release name=',name, ',but list does not have ',name,'at linenumber=',linenumber

	if printit==1:
		print thisprogram,"rlse name=",name,"new_list=",
		self.w_varlist_show()

	return list_rel


    def indentstack_release(self,fnlabel,searchw):
	thisfunc="indentstack_release"
	if self.printit>9:
		print thisfunc,"fnlabel=",fnlabel,"search=",searchw,"stack=",self.indentstack
        if len(fnlabel)>0:	
		nfirst=len(self.indentstack)
		idel=-1
		for i in range(len(self.indentstack)):
			word,label=self.indentstack[i]
			if word==searchw and label==fnlabel:
				# find the first one
				idel=i
				break
		if idel>=0:
			tok2=[]
			for i in range(idel):
				tok2.append(self.indentstack[i])
			self.indentstack=tok2
		nlast=len(self.indentstack)
	else:
		if self.printit>9:
			print thisfunc,'search unlabel'
		nfirst=len(self.indentstack)
                idel=-1
                for i in range(len(self.indentstack)):
                        word,label=self.indentstack[i]
			if self.printit>9:
				print thisfunc,"word,label=",word,label
                        if word==searchw:
				# find the last one
                                idel=i
                if idel>=0:
                        tok2=[]
                        for i in range(idel):
                                tok2.append(self.indentstack[i])
                        self.indentstack=tok2
                nlast=len(self.indentstack)

	if self.printit>9:
		print thisfunc,"fnlabel=",fnlabel,"idel=",idel,"new stack=",self.indentstack,"released=",nfirst-nlast
	return nfirst-nlast


    def get_labellist(self,tok):
	thisfunc='get_labellist'
	#print thisfunc,tok
	if tok[0].isdigit():
		return [tok[0]]
	labels=[]
	if tok[0]=='(':
		ic=0
		for nn in tok:
			if nn=='(':
				ic=ic+1
			elif nn==')':
				ic=ic-1
				if ic<=0:
					return labels
			elif nn==',':
				continue
			elif nn.isdigit():
				labels.append(nn)
			else:
				print thisfunc,'ERROR, unknown label',tok
	else:
		print thisfunc,'ERROR, strange label',tok
		
	return []


    def get_gotolabel(self,tok):
	thisfunc='get_gotolabel'
	#print thisfunc,tok
	if 'goto' in tok:
		for i in range(len(tok)):
			if tok[i]=='goto':
				return self.get_labellist(tok[i+1:])
	elif 'go' in tok:
		for i in range(len(tok)-2):
			if tok[i]=='go' and  tok[i+1]=='to':
				return self.get_labellist(tok[i+2:])
	else:
		print thisfunc,"ERROR ",tok
	print thisfunc,"ERROR ",tok
	print >>sys.stderr, thisfunc,"ERROR ",tok
	sys.exit(10)
	return []

    def struc_elem_append_uniq(self,w):
	thisfunc="struc_elem_append_uniq"
	if isinstance(self.struc_elem,list)==0:
		print thisfunc,"self.struc_elem must be a list"
		print >> sys.stderr, thisfunc,"self.struc_elem must be a list"
		sys.exit(10)
		
	if len(self.struc_elem)==0:
		self.struc_elem.append(w)
		return 
	for i in range(len(self.struc_elem)):
		nn=self.struc_elem[i]
		if nn.struc==w.struc and nn.elem==w.elem and nn.funcname==w.funcname:
			self.struc_elem[i].source=add_list_uniq(nn.source,w.source)
			self.struc_elem[i].target=add_list_uniq(nn.target,w.target)
			return 
	self.struc_elem.append(w)


    def w_varlist_append_uniq(self,w_var,linen,tok,printit=0):
	thisfunc="w_varlist_append_uniq"
	#print thisfunc,len(w_var),w_var
	[name,level,nameorg,sizeorig,type,prevar,routinename] = w_var.makelist()
	for w in self.w_varlist:
		if name==w.name and routinename==w.funcname:
		    	#if type!=w[4] and printit==1:
			#	print thisprogram,thisfunc,"ERROR, define again with anther type:",name
			#	print>>sys.stderr, thisprogram,thisfunc,"ERROR, define again with anther type:",name
			#if printit==1:
			#	print thisprogram,thisfunc,"warning, duplicated symbol:",name
			#print thisprogram,"line=",linen,tok
			return 
	self.w_varlist.append(w_var)

    def w_varlistall_append_uniq(self,w_var,linen,tok,printit=0):
        thisfunc="w_varlistall_append_uniq"
        [name,level,nameorg,sizeorig,type,prevar,routinename] = w_var.makelist()
        for w in self.w_varlistall:
                if name==w.name and routinename==w.funcname:
                        if type!=w.type and type!='undef':
			    self.w_varlistundel_append_uniq(t_reason([['redef'],name,routinename]),linen,tok)
			    if printit==1:
                                print thisprogram,thisfunc,"ERROR, define again with another type:",name,w.type,"and",type
                                print>>sys.stderr, thisprogram,thisfunc,"ERROR, define again with another type:",name,w.type,"and",type
			#if printit==1:
                        #	print thisprogram,thisfunc,"warning, duplicated symbol:",name
                        return
        self.w_varlistall.append(w_var)
	if type=='undef':
		self.w_varlistundel_append_uniq(t_reason([type,name,routinename]),linen,tok,printit)

    def w_varlistundel_append_uniq(self,ll,linen,tok,printit=0):
        thisfunc="w_varlistundel_append_uniq"
        if self.w_varlistundel_fix==1:
		return 
	reason,name,routinename=ll.makelist()
	for i in range(len(self.w_varlistundel)):
		if name==self.w_varlistundel[i].name and routinename==self.w_varlistundel[i].funcname:
			listreason=self.w_varlistundel[i].reason
			if isinstance(reason,list)==1:
				for areason in reason:
					listreason=list_adduniq(listreason,areason)
			elif isinstance(reason,basestring)==1:
				listreason=list_adduniq(listreason,reason)
			else:
				print thisfunc,"error, unknown reason type, reason=",reason
				sys.exit(10)
			self.w_varlistundel[i].reason=listreason
                        return
	if isinstance(reason,basestring)==1:
        	self.w_varlistundel.append(t_reason([[reason],name,routinename]))
	elif isinstance(reason,list)==1:
                self.w_varlistundel.append(t_reason([reason,name,routinename]))


    def get_defname(self,tok,istart):
	thisfunc='get_defname'
	name=tok[istart]
	namelist=[name]
	idone=istart
	for i in range(istart+1,len(tok)):
		if tok[i]=='(':
			# search ')'
			ic=0
			for j in range(i,len(tok)):
				namelist.append(tok[j])
				if tok[j]=='(':
					ic=ic+1
				elif tok[j]==')':
					ic=ic-1
					if ic==0:
						idone=j
						break
						
		else:
			break

		
		if idone>=0:
			break


	nameout=''
	for n in namelist:
		nameout = nameout+ n

	# size
	sizelist=[]
	idone=idone+1
	if tok[idone]!=',':
		print thisfunc,'ERROR not [,], i=',idone
		print thisfunc,tok,istart
		print >>sys.stderr,thisfunc,'ERROR not [,], i=',idone
		print >>sys.stderr,thisfunc,tok,istart
		sys.exit(10)
	idone=idone+1
	ic=0
	for i in range(idone,len(tok)):
		sizelist.append(tok[i])
		if tok[i]=='(':
			ic=ic+1
		elif tok[i]==')':
			ic=ic-1
			if ic<0:
				sizelist.pop()
				break
	sizeout=''
	for n in sizelist:
		sizeout = sizeout+n

	if self.printit>9:
		print thisfunc,tok,nameout,sizeout
	#print thisfunc,nameout,sizeout
	return [nameout,sizeout]

    def get_prevar(self,type,p="w"):
	thisfunc='get_prevar'
	prevar=''
	if type=='real(8)':
		#prevar='rv_w_'
		prevar='rv_'+p+'_'
	elif type=='integer':
		#prevar='iv_w_'
		prevar='iv_'+p+'_'
	elif type=='complex(8)':
		#prevar='zv_w_'
		prevar='zv_'+p+'_'
	else:
		print thisfunc,'unknown type=<'+type+'>'
		print >>sys.stderr, thisfunc,'unknown type=<'+type+'>'
		sys.exit(10)
	return prevar


    def get_kind(self,callname,p='w'):
	if callname in ['defrr','defdr','defr','redfrr']:
        	type='real(8)'
       	elif callname in ['defi','redfi']:
               	type='integer'
       	elif callname in ['defcc','defdc']:
               	type='complex(8)'
       	else:
               	print thisfunc, "unknown type=",callname
               	print >> sys.stderr, thisfunc, "unknown type=",callname
		sys.exit(10)

	prevar=self.get_prevar(type,p)

	return [type,prevar]


    def make_statement_alloc(self,toks,alist,routinename):
        thisfunc='make_statement_alloc'
        [name,level,nameorg,namesize,type,prevar]= alist

	tok,istart=toks
	tok_out=[]
	for i in range(len(tok)):
		if i<istart:
			tok_out.append(tok[i])
		else:
			break
	# search end of call def*
	ic=0
	for i in range(istart,len(tok)):
		if tok[i]=='(':
			ic=ic+1
		elif tok[i]==')':
			ic=ic-1
			if ic==0:
				iend=i
				break

	if namesize[0]=='-':
        	tok1='allocate('+prevar+name+'(abs('+namesize+')))'
	else:
        	tok1='allocate('+prevar+name+'('+namesize+'))'

	tok_out.append([tok1])

	tok1='if ('+namesize+'<0) ' 
	if prevar=='iv_w_':
		tok1=tok1+prevar+name+'(:)=0'
	else:
		tok1=tok1+prevar+name+'(:)=0.0d0'

       	tok_out.append([tok1])

	for i in range(iend+1,len(tok)):
		tok_out.append(tok[i])

	#print thisfunc,"original=",tok
	#print thisfunc,"changed =",tok_out
        return tok_out


    def make_statement_rlse_dealloc(self,toks,lists):
        thisfunc='make_statement_rlse_dealloc'
	tok,istart=toks

	# 'if ( ) call rlse' -> [if () then] [dealloc...] [endif]
	ihave_if=-1
	ihave_then=-1
	if 'if' in tok[:istart]:
		ihave_if=tok[:istart].index('if')
		if 'then' in tok[:istart]:
			ihave_then=tok[:istart].index('then')
	if ihave_if==-1 and ihave_then>=0:
		print thisfunc,"found then, not failed to find if",tok,istart
		print >> sys.stderr, thisfunc,"found then, not failed to find if",tok,istart
		sys.exit(10)
	if ihave_if>=0 and ihave_if!=0:
		print thisfunc,"analysis ERROR (if) ",tok,istart
		print >>sys.stderr, thisfunc,"analysis ERROR (if) ",tok,istart
		sys.exit(10)
	if ihave_then>=0 and ihave_then!=istart-1:
		print thisfunc,"analysis ERROR (then) ",tok,istart
		print >>sys.stderr, thisfunc,"analysis ERROR (then) ",tok,istart
		sys.exit(10)

	iadd_endif=0
	if ihave_if>=0 and ihave_then==-1:
		iadd_endif=1
		
	# start copying
        tok_out=[]
        for i in range(len(tok)):
                if i<istart:
                        tok_out.append(tok[i])
                else:
                        break

	if ihave_if>=0 and ihave_then==-1:
		tok_out.append('then')
        # search end of call def*
        ic=0
        for i in range(istart,len(tok)):
                if tok[i]=='(':
                        ic=ic+1
                elif tok[i]==')':
                        ic=ic-1
                        if ic==0:
                                iend=i
                                break

	tok2=self.make_statement_dealloc(lists,1) #=1 means also call rlse()
	tok_out.append(tok2)

	if iadd_endif==1:
		tok_out.append('endif')
	#print thisfunc,tok_out
        return tok_out


    def make_statement_dealloc(self,lists,alsorlse=0):
	thisfunc="make_statement_dealloc"
	if isinstance(lists,list)==0:
		return []
	if len(lists)==0:
		return []
        tok2=[]
	listrlse=[]
        lists.reverse()
        for  alist in lists:
                [name,level,nameorg,namesize,type,prevar,routinename]= alist.makelist()
		id= self.w_varlistundel_contains(name,self.execfuncname[0])
		if id>=0:
			listrlse.append(name)
			if self.printsource==1:
				print thisprogram,"not deallocate",name+" because of "+list2string(self.w_varlistundel[id].reason)
			continue
		nn = self.w_pointerlist_find(name,self.execfuncname[0])
		if isinstance(nn,t_var7)==1:
			#listrlse.append(name)
                        if self.printsource==1:
                                print thisprogram,"not deallocate",name+" because of pointer"
                        continue
                realname=prevar+name
                tok='if (allocated('+realname+')) deallocate('+realname+')'
                tok2.append([tok])
	if len(listrlse)>0 and alsorlse==1:
		n=len(listrlse)-1
		name=listrlse[n]
		namedic=self.var_dic[listrlse[n]]
		if namedic[0].find("(")>=0:
		    if self.printsource==1:
			print "#error, rlse? ",name,"=",namedic
			print >>sys.stderr,"ERROR, rlse? ",name,"=",namedic
		tok='call rlse('+name+')'
		tok2.append(tok)
	return tok2

    def make_statement_dealloc_uselist(self,lists):
        thisfunc="make_statement_dealloc_uselist"
        if isinstance(lists,list)==0:
                return []
        if len(lists)==0:
                return []
        tok2=[]
        lists.reverse()
        for  alist in lists:
                [name,level,nameorg,namesize,type,prevar,routinename]= alist
                id= self.w_varlistundel_contains(name,self.execfuncname[0])
                if id>=0:
                        print thisprogram,"not deallocate",name+"("+self.w_varlistundel[id][0]+")"
                        continue
                realname=prevar+name
                tok='if (allocated('+realname+')) deallocate('+realname+')'
                tok2.append([tok])
        return tok2




    def get_rlsename(self,tok,istart):
	thisfunc='get_rlsename'
	#tok[istart]=='rlse'
	if tok[istart]!='rlse':
		print thisfunc,"ERROR statement mismatch, not rlse",tok[istart]
		print >>sys.stderr,thisfunc,"ERROR statement mismatch, not rlse",tok[istart]
		sys.exit(10)
	istart=istart+1
	if tok[istart]!='(':
                print thisfunc,"ERROR statement mismatch, not (",tok[istart]
                print >>sys.stderr, thisfunc,"ERROR statement mismatch, not (",tok[istart]
                sys.exit(10)
	istart=istart+1
	name1=tok[istart]
	namelist=[]
	ic=0
	for i in range(istart,len(tok)):
		namelist.append(tok[i])
		if tok[i]=='(':
			ic=ic+1
		elif tok[i]==')':
			ic=ic-1
			if ic<0:
				namelist.pop()
				break
	name2=''
	for tok in namelist:
		name2=name2+tok
	return [name1,name2]

 
    def analyze_wvar(self,tok,linenumber,fnlabel):
	thisfunc="analyze_wvar"
	level=-1
	size=0
	type='undef'
	prevar='x'
	nameorg='undef'
	for i in range(len(tok)-2):
		if tok[i]=="w" and tok[i+1]=='(' :
			icmax,var=get_var(tok[i:])
			namenext=var[3]
			# reason if there is 
			if namenext=='(':
				reason='array'
			else:
				reason='wvar'
			#      w ( oa ( i ) ) 
			#      w ( oa + i )
			if len(var)>4:
			    for nn in var[2:]:
				name=nn
				if name=='(' or name==')' or name=='+' or name=='-' or name=='*' or name=='/':
					continue
				if name.isdigit()==1:
					if self.printsource==1:
						print thisprogram,name,"is digit, dropped"
					continue
				if name[0]=='o':
					self.w_varlistundel_append_uniq(t_reason([[reason],name,self.execfuncname[0]]),linenumber,tok,self.printsource)
					if '+' in var or '-' in var:
						print "#error  w(o+i) type access",list2string(var),"in",self.execfuncname[0]
						print >> sys.stderr ,  "ERROR , w(o+i) type access",list2string(var),"in",self.execfuncname[0]
			elif len(var)==4:
			# w ( i ) = 4 
				name=var[2]
				if name.isdigit()==1: 
			  	    if self.printsource==1:
					print thisprogram,name,"is digit, dropped"
				else:
				    self.w_varlistall_append_uniq(t_var7([name,level,nameorg,size,type,prevar,self.execfuncname[0]]),linenumber,tok,self.printsource)
			else:
				print thisfunc,"ERROR strange var=",var 
				sys.exit(10)


    def connect_element(self,tok):
	thisfunc="connect_element"
	print thisfunc,tok
	newvar=''
	id=0
	if tok[id].isalpha==1 and tok[id+1]=='(':
		ic=1
		newvar=tok[id]+tok[id+1]
		g_i
		for i in range(id+2,len(tok)):
			g_i=i
			newvar=newvar+tok[i]	
			if tok[i]=='(':
				ic=ic+1
			if tok[i]==')':
				ic=ic-1
			if ic==0:
				break	
		if tok[g_i+1][0]=='%':
			newvar=newvar+tok[g_i+1]
		print thisfunc,"var=",newvar
	elif tok[id].find('%')>=0:
		newvar=tok[id]
		print thisfunc,"var=",newvar
	else:
		print thisfunc,"ERROR1",tok


    def element_connect(self,tok):
	thisfunc="element_connect"
	id=0
	l_tok=[]
	while id<len(tok):
		if tok[id][0].isalpha() and id+1<len(tok) and tok[id+1]=='(':
			ic = 1
			newvar=tok[id]+tok[id+1]
			id=id+2
			ifound=0
			for i in range(id,len(tok)):
				newvar=newvar+tok[i]
				if tok[i]=='(':
					ic=ic+1
				elif tok[i]==')':
					ic=ic-1
				if ic==0:
					id=i
					ifound=1
					break
			if ifound==0:
				print thisfunc,"ERROR",tok,ic
				sys.exit(10)

			if id+1<len(tok) and tok[id+1][0]=='%':
				id=id+1
				newvar=newvar+tok[id]
			l_tok.append(newvar)
			newvar=''
		else:
			l_tok.append(tok[id])		
		id=id+1
	return l_tok
			

    def w_pointerlist_append_uniq(self,w7):
	thisfunc="w_pointerlist_append_uniq"
	if isinstance(w7,t_var7)==0:
		print thisfunc,"ERROR in type",t_var7
		sys.exit(10)
	for nn in self.w_pointerlist:
		if nn.name==w7.name and nn.funcname==w7.funcname:
			if nn.type != w7.type:
				print thisfunc,"ERROR , different type",nn,w7
				print >> sys.stderr,thisfunc,"ERROR , different type",nn,w7
				sys.exit(10)
			else:
				return
	self.w_pointerlist.append(w7)

    def change_wpointer_eq(self,eqtype,tok,wvar7,var,vardic,varinfo,struc,strucdic,linenum):
	thisfunc="change_wpointer_eq"
	#print thisfunc,
	#print eqtype,tok,"/wvar7=",wvar7,"/var=",var,"/vardic=",vardic,
	#print "/varinfo=",varinfo,"/struc=",struc,"/strucdic=",strucdic
	elementtype=wvar7.type

	if tok[1]!='=':
		print thisfunc,"ERROR tok[1] must be  =. tok=", tok
		print >>sys.stderr, thisfunc,"ERROR tok[1] must be  =. tok=", tok
		exit(10)

        struc_real=strucdic[1]
        struc_real= struc_real.replace("type(","")
        struc_real= struc_real.replace(")","")
        struc_type=self.struc_element_type(struc_real,struc.element)
        if struc_type=="unchange":
		if self.printsource==1:
			print thisprogram,"unchange", struc.name,"%",struc.element
		self.w_varlistundel_append_uniq(t_reason(["o_eq",var.name,self.execfuncname[0]]),linenum,tok)
		return [0,tok]

	if eqtype=='vs':
		ivar=0
		istruc=2
	else:
		ivar=2
		istruc=0

	p_var=tok[ivar]
	p_var_array=var.array
	p_struc=tok[istruc]
	p_struc_array=struc.array
	prevar=wvar7.prevar

	wvar= self.w_varlistall_find(p_var,self.execfuncname[0])
	if isinstance(wvar,t_var7)==1 and wvar.type!='undef' and eqtype=='vs' :
	    if self.printsource==1:
		print thisprogram,"warning, pointer and allocate",p_var,tok, 
		print ", call def*=",wvar
		

		

	id= self.w_varlistundel_contains(p_var,self.execfuncname[0])
	if id>=0:
		if self.printsource==1:
			print "#error ",p_var,"is changed, but check if I can change it because of",list2string(self.w_varlistundel[id].reason)


	if p_var_array[0]==0:
		p_var2=prevar + p_var 
	else:
		print thisfunc,"ERROR array var is not supported",p_var
		sys.exit(10)

	l_tmp=p_struc
	l_tmp=do_split0__ils_ol(l_tmp,'(')
	l_tmp=do_split0__ils_ol(l_tmp,')')
	l_tmp=do_split0__ils_ol(l_tmp,',')
	l_tmp=do_split0__ils_ol(l_tmp,'%')
	if p_struc_array[1]==0:
		l_struc2=[]
		for i in range(len(l_tmp)):
			if l_tmp[i]=='%':
				l_struc2.append(l_tmp[i])
				l_struc2.append(prevar)
			else:
				l_struc2.append(l_tmp[i])
	else:
		print thisfunc,"ERROR array element is not supported"
                sys.exit(10)
				
	p_struc2=list2str(l_struc2,'')
	if eqtype=='vs':
		l_tok=[p_var2,"=>",p_struc2]
		self.struc_elem_append_uniq(t_struc_elem([struc.name,struc.element,self.execfuncname[0],p_var2,'']))
        else:
		l_tok=[p_struc2,"=>",p_var2]
		self.struc_elem_append_uniq(t_struc_elem([struc.name,struc.element,self.execfuncname[0],'',p_var2]))

	return [1,l_tok]

	
    def do_element_in_type(self,tok,linenumber,fnlabel):
	thisfunc="do_element_in_type"
	ifound_perc=0
	if len(tok)>0:
	    if isinstance(tok[0],basestring)==1:
	    	for nn in tok:
		    if isinstance(nn,basestring)==0:
			print thisprogram,thisfunc,"possible ERROR, tok is a list, tok=",nn,tok
			continue

		    percpos=nn.find('%')
		    if len(nn)>0 and percpos>=0 and nn[0]!="'":
			var = nn.split('%')
			# ask kind of structure
			if len(var)>=2 and var[1][0]=='o':
				ifound_perc=1
				break
		pos_eq=-1
		for  i in range(len(tok)):
			if tok[i]=='=':
				pos_eq=i
				break

		if ifound_perc==1 and pos_eq>=0:
			l_tok=self.element_connect(tok)
			if l_tok[0].find('%')>=0:
				# source =l_tok[0]
				# target= l_tok[2]
				#targetname,dummy1,targetname_orig,targetarray= separate_var([l_tok[2]])
                                var = separate_var([l_tok[2]])
                                #strucvar[0],strucvar[1],strucname_orig,strucarray= separate_var([l_tok[0]])
                                struc = separate_var([l_tok[0]])
				eqtype="sv"
			elif l_tok[2].find('%')>=0:
                                var = separate_var([l_tok[0]])
                                struc = separate_var([l_tok[2]])
				eqtype="vs"
			if not (struc.name in ["size","int"]):
                                varinfo=self.w_varlistall_find(var.name,self.execfuncname[0])
                                if isinstance(varinfo,t_var7)>0:
                                        vartype=varinfo.type
                                else:   
                                        vartype='None'
				if vartype!="None" and vartype !="undef":
					varinfo.prevar=self.get_prevar(vartype,"p")
				vardic=[]
				strucdic=[]
				if len(var.name)>0 and var.name[0].isalpha()==1:
				    try:
					vardic=self.var_dic[var.name]
				    except:
					print >> sys.stderr, "var.name=",var.name
				if len(struc.name)>0 and struc.name[0].isalpha()==1:
				    try:
					strucdic=self.var_dic[struc.name]
				    except:
					print >> sys.stderr, "struc.name=",struc.name
                                #print thisfunc," var=",eqtype,var.name,var.tok,vardic,var.array,varinfo,l_tok
				#print thisfunc,"",struc
                                #print thisfunc," struc=",eqtype,struc.name,struc.tok,strucdic,"source=",varinfo,struc.array,l_tok
				# struc=strucdic[1], element=struc.element
				# type of var = 
				struc_real=strucdic[1]
				struc_real= struc_real.replace("type(","")
				struc_real= struc_real.replace(")","")
				type_of_var = self.struc_element_type(struc_real,struc.element)
				#print thisfunc," pointer=",var.name,type_of_var
				if type_of_var!="undef" and type_of_var !='unchange':
					prevar=self.get_prevar(type_of_var,"p")
				else:
					prevar=""
				wvar7=t_var7([var.name,0,var.name,'*',type_of_var,prevar,self.execfuncname[0]])
				self.w_pointerlist_append_uniq(wvar7)
				newtok=self.change_wpointer_eq(eqtype,l_tok,wvar7,var,vardic,varinfo,struc,strucdic,linenumber)
				return newtok
				
		
	    elif isinstance(tok[0],list)==1:
		listtok=[]
		ichanged=0
		for nn in tok:
			rettok=self.do_element_in_type(nn,linenumber,fnlabel)
			ichanged = ichanged and rettok[0]
			listtok.append(rettok[1])
		tok=[ichanged,listtok]
	    else:
		print thisfunc,"ERROR, strange?",tok
		sys.exit(10)

	return [0,tok]


    def do_function(self,tok,linenumber,fnlabel):	
	thisfunc="do_function"
	ichanged=0


	if self.printit>9:
		print thisfunc,linenumber,tok,fnlabel
	if tok[0]=="if":
		if "then" in tok:
			self.indentstack.append(['if',''])
	if 'endif' in tok or (len(tok)>=2 and tok[0]=="end" and tok[1]=="if"):
		self.indentstack_release(fnlabel,'if')

	if tok[0]=="do":
		# do
		# do i=...
		# do 100 i=...
		if len(tok)>1:
			a=re.match("[0-9]+",tok[1])
			if a != None:
				if self.printit>9:
					print "do number=",tok[1];
				self.indentstack.append(["do",tok[1]])
			else:
				self.indentstack.append(["do",''])
		else:
			self.indentstack.append(["do",""])
	if 'enddo' in tok or 'end' in tok:
		endlist=[]
		for i in range(len(tok)):
			if tok[i]=='enddo' or (i+1<len(tok) and  tok[i]=='end' and tok[i+1]=='do'):
				endlist.append(i)
		if self.printit>9:
			print thisfunc,tok
			print thisfunc,"# of enddo=",len(endlist),endlist
		if len(endlist)>0:
		    for i in range(len(endlist)):
			n=self.indentstack_release('','do')
	if len(fnlabel)>0:
		if self.printit>9:
			print thisfunc,'search do',fnlabel
                n=self.indentstack_release(fnlabel,'do')
	if tok[0]=="select":
		self.indentstack.append(["select",fnlabel])
	if tok[0]=="endselect" or (len(tok)>=2 and tok[0]=="end" and tok[1]=="select"):
		self.indentstack.pop()
	if  'call' in tok:
	    iloopagain=1
	    l_tok=tok
	    while iloopagain==1:	
		iloopagain=0
		for i in range(len(l_tok)-3):
			if l_tok[i]=='call' and (l_tok[i+1] in ["defdr","defrr","defr","defcc",'defdc',"defi",'redfi','redfrr']):
				name=l_tok[i+3]
				id= self.w_varlistundel_contains(name,self.execfuncname[0])
				level=len(self.indentstack)
				nameorg,namesize=self.get_defname(l_tok,i+3)

				nn= self.w_pointerlist_find(name,self.execfuncname[0])
				if isinstance(nn,t_var7)==1:
					type,prevar=self.get_kind(l_tok[i+1],'p')
				else:
					type,prevar=self.get_kind(l_tok[i+1])

				self.w_varlist_append_uniq(t_var7([name,level,nameorg,namesize,type,prevar,self.execfuncname[0]]),linenumber,tok,self.printsource)
				self.w_varlistall_append_uniq(t_var7([name,level,nameorg,namesize,type,prevar,self.execfuncname[0]]),linenumber,tok,self.printsource)
				# new statement 
				if id<0:
				    if name==nameorg:
					l_tok=self.make_statement_alloc([l_tok,i],[name,level,nameorg,namesize,type,prevar],self.execfuncname[0])
					ichanged=1
					iloopagain=1
					#print thisfunc,"def",l_tok,ichanged
				    else:
					self.w_varlistundel_append_uniq(t_reason(['array',name,self.execfuncname[0]]),linenumber,tok,0)
					ichange_ERROR=1
			elif l_tok[i]=='call' and l_tok[i+1]=='rlse':
				#name=tok[i+3]
				name,nameorg=self.get_rlsename(l_tok,i+1)
				level=len(self.indentstack)
				# new statement
				if name==nameorg:
					list_rel=self.w_varlist_release(name,level,linenumber,self.printsource)
					l_tok=self.make_statement_rlse_dealloc([l_tok,i],list_rel)
					ichanged=1
					iloopagain=1
					#print thisfunc,"rlse",l_tok,ichanged
				else:
					ichange_ERROR=1
			if iloopagain==1:
				break
				# for i loop again

	    if ichanged==1 and self.printit==1:
			print thisfunc,"changed",l_tok

	    tok=l_tok

	#if tok[0]=="continue":
	#	if self.printit==1:
	#		print thisfunc,"continue(enddo?)",linenumber,tok,fnlabel
	#	if len(fnlabel)>0:
	#	    n=self.indentstack_release(fnlabel,'do')
	if "goto" in tok:
		label=self.get_gotolabel(tok)
		if self.printit>9:
			print "goto"
		for nn in label:
			if self.printsource==1:
				print >>sys.stderr , "ERROR, possible spagettie code at line",linenumber,"label=", nn
	if "go" in tok:
		#search go and to
		for n in range(len(tok)-1):
			if tok[n]=="go" and tok[n+1]=="to":
				label=self.get_gotolabel(tok)
				for nn in label:
					if self.printsource==1:
						print >>sys.stderr , "ERROR, possible spagettie code at",linenumber,"label=", nn
 				break
					
				
	if "return" in tok:
		for name in tok:
			if self.w_varlistall_contains(name,self.execfuncname[0]):
				print "#error tok contains return and",name
				print >>sys.stderr, "#error tok contains return and",name
				sys.exit(10)
		if self.printit>9:
			print "return"
		#if self.printsource and len(self.w_varlist)>0:
		#	print "#error, have return with len(w_varlist)>0 at line",linenumber
		#	print >>sys.stderr ,"ERROR, have return with len(w_varlist)>0",linenumber
		l_tok=self.deallocate_var(linenumber)
		#if len(fnlabel)>0 and len(l_tok)>0:
                #        print thisfunc,"ERROR in analysis return, have_fnlabel", tok,l_tok,fnlabel
                #        print >> sys.stderr, thisfunc,"ERROR in analysis return,", tok,l_tok,fnlabel
                #        sys.exit(10)
		if len(l_tok)>0:
			l_tok=insert_tok(tok,'return',l_tok)
			return [1,l_tok]
		else:
			l_tok=tok
			return [0,[]]

#        self.analyze_wvar(tok,linenumber,fnlabel)

        if ichanged==1:
               return [1,tok]
        else:
               return [0,tok]



    def typesize__is_ol(self,tok,istart0):
	# *4 *8 *(2) (4) 
	i=istart0
	size=0
	success=1
	if tok[i]=='*':
		try:
			iint=int(tok[i+1])
		except ValueError:
			iint=0
		if tok[i+1]=='(':
			if tok[i+3]==')':
				size=tok[i+2]
				i=i+4
			else:
				success=-1
		elif tok[i+1].isdigit()==1:
			size=tok[i+1]
			i=i+2
		else:
			success=-1
	elif tok[i]=='(' :
	    if  tok[i+2]==')':
		size=tok[i+1]
		i=i+3
	    else:
		success=-1
	return [success,i,size]


    def do_typevar__il_on(self,tok):
        thisfunc='do_typevar__il_on'
	typestr0=''
	typesize='0'
	istart=-1
	ERROR=0
	if tok[0]=='type':
		istart=1
		typesize=-1	
		typestr0=tok[0]
		l_out=self.typesize__is_ol(tok,istart)
                if l_out[0]:
                        istart=l_out[1]
                        typesize=l_out[2]
                else:
                        ERROR=1
	elif tok[0]=='integer' or tok[0]=='real':
		istart=1
		typesize='0'
		typestr0=tok[0]
		l_out=self.typesize__is_ol(tok,istart)
		if l_out[0]:
			istart=l_out[1]
			typesize=l_out[2]
		else:
			ERROR=1
	elif tok[0]=='character':
		istart=1
		typesize='0'
		typestr0=tok[0]
		l_out=self.typesize__is_ol(tok,istart)
                if l_out[0]:
                        istart=l_out[1]
                        typesize=l_out[2]
                else:
                        ERROR=1
	elif tok[0]=='double' and tok[1]=='precision':
		istart=2
		typesize='8'
		typestr0='real'
	elif tok[0]=='logical':
		istart=1
		typesize='0'
		typestr0='logical'
		l_out=self.typesize__is_ol(tok,istart)
                if l_out[0]:
                        istart=l_out[1]
                        typesize=l_out[2]
                else:
                        ERROR=1
	elif tok[0]=='complex':
		istart=1
		typesize='4'
		typestr0='complex'
		l_out=self.typesize__is_ol(tok,istart)
                if l_out[0]:
                        istart=l_out[1]
                        typesize=l_out[2]
                else:
                        ERROR=1
	elif tok[0]=='double' and tok[1]=='complex':
		istart=2
		typesize='8'
		typestr0='complex'
	elif tok[0]=='external' or tok[0]=='parameter' or tok[0]=='data':
		tok2=[]
		return tok2	
	else:
		ERROR=1
	if ERROR==1:
		print thisfunc,"ERROR",tok
		print >> sys/stderr , thisfunc,"ERROR",tok
		sys.exit(10)

	tok=merge_w__il_ol(tok[istart:])
        tok=merge_stringslash__il_ol(tok)
	level=1
        tok=merge_variable__il_pi_ol(tok,level)
        tok=merge_ast__il_ol(tok)
	found=0
	for i in range(len(tok)):
		if tok[i]=='::':
			istart=i+1
			found=1
			break
	supp0=[]
	if found==1:
		if istart-1>=0:
			supp0=tok[:istart-1]
		tok=tok[istart:]
	supp=''
	for nn in supp0:
		supp=supp+nn
	

	#make type
	if typesize==0:
		typedef=typestr0
	else:
		typedef=typestr0+'('+typesize+')'

	tok2=[]
	for nn in tok:
		if nn==',': 
			continue
		ll=separate_array__il_ol(nn)
		alist=[ll[0],ll[1],typedef,supp]
		tok2.append(alist)

	return tok2

    def manip_varlist(self,varlist):
	thisfunc='manip_varlist'
	tok2=[]
	# var_toadd
	#print thisfunc,"varlist=",varlist

	if len(varlist)==0:
		return tok2
	[var,var0,typestr,supp]=varlist[0]
		# use only typestr
	deftypeorig=typestr+supp+"::"
	deforigadd=[]
	for nn in varlist:
		[struc,strucorig,struc_typedef,supp]=nn
		# check whether it is in the list, var_toadd
		changed=0
		wpointer=self.w_pointerlist_find(struc,self.execfuncname[0])
		if isinstance(wpointer,t_var7)==1 and wpointer.type!="unchange":
			if self.printsource==1:
				print thisprogram,"POINTER change",wpointer
			id= self.w_varlistundel_contains(struc,self.execfuncname[0])
			if id==-1:
				aline = wpointer.type + ",pointer :: "+wpointer.prevar+ wpointer.name+ '(:)'
				changed=1
				tok2.append(aline)
	
		if changed==0:
		    id= self.w_varlistundel_contains(struc,self.execfuncname[0])
		    if id>=0:
		    	if self.printsource==1:
				print thisprogram,"do not change ",self.w_varlistundel[id]
		    else:
		    	for n2 in self.w_varlistall:
			    name,level,nameorg,namesize,type,prevar,execfuncname0 = n2.makelist()
			    if (execfuncname0==self.execfuncname[0] or execfuncname0=='*') and  name==struc:
				#print thisfunc,"1>",var,var0,typestr,supp,self.execfuncname
				#print thisfunc,"2>",name,level,nameorg,namesize,type,prevar,execfuncname0
				aline = type+' ,allocatable :: '+prevar+name+'(:)'
				changed=1
				tok2.append(aline)
				break

		if changed==0:
			deforigadd.append(strucorig)
			deforigadd.append(',')

	n=len(deforigadd)
	deforig=[]
	if n>1 and deforigadd[n-1]==',':
		deforigadd.pop()
		deforig.append(deftypeorig)
		for nn in deforigadd:
			deforig.append(nn)	
	if len(tok2)>0:
		tok3=deforig
		tok3.append(tok2)
		return [1,tok3]
	else:
		return [0,'']
			
	

    def do_use__il_on(self,tok):
	thisfunc="do_use__il_on"
	nn=tok[1]
	if nn[0].isalpha()==1:
		self.uselist.append([self.execfuncname,nn])

    def add_use__in_ol(self):
	#search the same execfuncname[0] in  uselist
	thisfunc='add_use__in_ol'
	found=0
	for nn in self.use_toadd:
		if nn[0]==self.execfuncname[0]:
			found=1
	if found==1:
		# search the same execfuncname[0] in varlist
		tok2=[]
		for nn in self.var_toadd:
			if self.execfuncname[0]==nn[0]:
				tok2.append(nn[1])
		if len(tok2)==0:
			print thisfunc,"ERROR tok2==[]",self.execfuncname
			print >>sys.stderr, thisfunc,"ERROR tok2==[]",self.execfuncname
			sys.exit(10)
		l_tok=[]
		l_tok.append(['','use m_struc_def  !'+thisprogram])
		#l_tok.append('use m_struc_def , only:')
		#for i in range(len(tok2)):
		#	l_tok.append('s_'+tok2[i])
		#	if i!=len(tok2)-1:
		#		l_tok.append(',')
		return l_tok
	else:
		return []


    def do_w_var(self,w_var,tok,linenumber,fnlabel):
	thisfunc="do_w_var"
	name=w_var.name
	prefix=w_var.prevar
	tok0=tok

	if '=>' in tok:
		return [tok,0]
	
	# w(var) or var itself
	idone=0
	ichanged=1
	while ichanged==1:
	    ichanged=0
	    for i in range(len(tok)):
		if tok[i]==name:
			if i-2>=0 and tok[i-2]=='w' and tok[i-1]=='(':
				ic=0
				# indent ( )
				for j in range(i-2,len(tok)):
					if tok[j]=='(':
						ic=ic+1
					elif tok[j]==')':
						ic=ic-1
						if ic<=0:
							break
				# copy
				l_tok=tok[:i-2]
				l_tok.append(prefix+name)
				for k in range(j+1,len(tok)):
					l_tok.append(tok[k])
				tok=l_tok
				ichanged=1
				idone=1
				break

			elif i-2>=0 and tok[i-2] in ["defdr","defrr","defr","defcc",'defdc',"defi"] and tok[i-1]=='(':
				idonothing=1
				return [tok0,0]

	if idone==0:
		if '=' in tok:
			if self.printsource==1:
				print thisprogram,"eq",w_var.name,tok
			self.w_varlistundel_append_uniq(t_reason(['eq',w_var.name,w_var.funcname]),linenumber,tok,0)
                	return [tok,0]

		self.idisable_change=1
		if self.printsource==1:
			print thisprogram,"pass",w_var.name,tok
		self.w_varlistundel_append_uniq(t_reason(['pass',w_var.name,w_var.funcname]),linenumber,tok,0)
		return [tok,0]

	#print thisfunc,'changed:',tok
	return [tok,1]



    def  w_varlist_addall(self):
	nn=''
        for atom in self.w_varlist:
                nn=nn+ atom.name+' '
	return nn


    def w_varlist_show(self):
	if len(self.w_varlist)>0:
	    for atom in self.w_varlist:
		print atom.name,
	    print ''
	else:
	    print '(None)'


    def w_varlistall_addall(self):
	nn=''
        for atom in self.w_varlistall:
                nn=nn+atom.name+' '
	return nn


    def w_varlistall_show(self):
	for atom in self.w_varlistall:
		print atom.name,
	print ''

    def w_varlistundel_addall(self):
        nn=''
        for name in self.w_varlistundel:
	    if name.funcname==self.execfuncname[0]:
                nn=nn+name.name+'('+list2string(name.reason,'','')+') '

	for name in self.w_varlistall:
		if name.name in self.callvarlist and name.funcname==self.execfuncname[0]:
			nn=nn+name.name+'(call) '
	return nn

    def w_varlistnotundel_addall(self):
	listundel=[]
	for name in self.w_varlistundel:
	    if name.funcname==self.execfuncname[0]:
		listundel.append(name.name)
	errorlist=[]
	for name in self.w_varlist:
		if not name.name in listundel and name.funcname==self.execfuncname[0]:
			errorlist.append(name.name)
	nn=''
	if len(errorlist)==0:
		return nn
	for name in errorlist:
		nn=nn+name+' '
	return nn	

    def w_varlist_contains(self,name,funcname): 
	thisfunc='w_varlist_contains'
	for w_var in self.w_varlist:
		if name==w_var.name and (w_var.funcname==funcname or funcname=='*'):
			return 1
	return 0

    def w_varlistall_contains(self,name,funcname): 
	thisfunc='w_varlistall_contains'
	for w_var in self.w_varlistall:
		if name==w_var.name and (w_var.funcname==funcname or funcname=='*'):
			return 1
	return 0

    def w_varlistall_find(self,name,funcname):
        thisfunc='w_varlistall_contains'
        for w_var in self.w_varlistall:
                if name==w_var.name and (w_var.funcname==funcname or funcname=='*'):
                        return w_var
        return []

    def w_pointerlist_find(self,name,funcname):
	thisfunc="w_pointerlist_find"
	for nn in self.w_pointerlist:
		if nn.name==name and nn.funcname==funcname:
			return nn
	return []


    def indentstring(self,level):
	spc='    '
	s=''
        for i in range(level):
		s=s+spc
	return s	


    def deallocate_var(self,linenumber):
	thisfunc='deallocate_var'
        # self.w_varlistundel
	#if (self.w_varlistall)>0:
	#	return []

        if  self.printsource==1:
                print thisprogram,"w_varlist remains:",
                self.w_varlist_show()
        nn=[]
        for n2 in self.w_varlistundel:
	    if n2.funcname == self.execfuncname[0]:
                nn.append(n2.name)
        if len(nn)>0 and self.printsource==1:
                print thisprogram,"w_varlistundel:",list2string(nn)
	for n2 in self.w_pointerlist:
		if n2.funcname==self.execfuncname[0]:
			nn.append(n2.name)
        if len(nn)>0 and self.printsource==1:
                print thisprogram,"w_varlistundel(+pointer):",list2string(nn)
        # make list in self.w_varlist ,but not nn
        listshoulddel=[]
	listshouldnotdel=[]
        for n2 in self.w_varlist:
                if not n2.name in nn:
                        listshoulddel.append(n2.name)
		else:
			listshouldnotdel.append(n2.name)
        if self.printsource==1:
		if len(listshouldnotdel)>0:
			str=list2string(listshouldnotdel)
		else:
			str='(None)'
                print thisprogram,"w_varlist (undel), remains:",str
        if self.printsource==1:
		if len(listshoulddel)>0:
			str=list2string(listshoulddel)
		else:
			str='(None)'
                print thisprogram,"w_varlist (del), remains:",str

        #print thisprogram,"dealloc them for safty"
        l_tok=[]
        l_tok=self.make_statement_dealloc(self.w_varlist)

        return l_tok


    def do_at_end_of_subroutine(self,linenumber):
	thisfunc='do_at_end_of_subroutine'
	self.write_newvar()

	l_tok= self.deallocate_var(linenumber)

        #if len(self.w_varlistall)>0:
        #       print "w_varlist used (all):",self.w_varlistall

	return l_tok

    def  do_at_beginning_of_subroutine(self,tok):
		thisfunc="do_at_beginning_of_subroutine"
                self.var_dic.clear()
                self.indentstack=[]
                self.w_varlist=[]
                #self.w_varlistall=[]
                #self.w_varlistundel=[]
		self.idisable_change=0
		self.callvarlist=[]

		if len(tok)==0:
			return
		is0=0
		ie=len(tok)-1
		#print thisfunc,tok
		if tok[is0]=='(' and tok[ie]==')':
			for i in range(is0+1,ie):
				if tok[i] ==',':
					continue
				self.callvarlist.append(tok[i])	
		else:
			print thisfunc,"ERROR position of ( and/or )", tok
			print >>sys.stderr, thisfunc,"ERROR position of ( and/or )", tok
			sys.exit(10)
		#print thisfunc,'call var list=',self.callvarlist

    def w_varlistundel_contains(self,name,funcname):
	thisfunc='w_varlistundel_contains'
	for i in range(len(self.w_varlistundel)):
		nn=self.w_varlistundel[i]
		if nn.name==name and nn.funcname==funcname:
			return i
	return -1


#input,string: line
    def do_line__is_ol(self,line,linenumber,fnlabel):
	thisfunc='do_line__is_ol'
	tok2=line2tok__is_ol(line)

	# lower case
	tok=[]
	for nn in tok2:
		if nn[0]=='\'':
			tok.append(nn)
		else:
			tok.append(nn.lower())

	if len(tok)==0:
		return [0,'']

	if (len(tok)==1 and tok[0]=='end'):
		if len(self.execfuncname)==0:
			print "ERROR,",thisfunc,", found 'end', but function name is unknown"
			print tok
			print >>sys.stderr, "ERROR,",thisfunc,", found 'end', but function name is unknown"
			print >>sys.stderr,tok
			sys.exit(10)
		l_tok=self.do_at_end_of_subroutine(linenumber)
		if len(l_tok)>0:
			return [3,l_tok]
		else:
			return [0,[]]


	if len(tok)>=2 and (tok[0]=='end' and (tok[1]=='subroutine' or tok[1]=='function')):
                if len(self.execfuncname)==0:
                        print "ERROR,",thisfunc,", found 'end', but function name is unknown"
                        print tok
                        print >>sys.stderr,"ERROR,",thisfunc,", found 'end', but function name is unknown"
                        print >>sys.stderr,tok
                        sys.exit(10)
		l_tok=self.do_at_end_of_subroutine(linenumber)
		if len(l_tok)>0:
			return [3,l_tok]
		else:
                	return [0,[]]


	if tok[0]=='subroutine':
		#self.do_at_end_of_subroutine()

		self.execfuncname=[tok[1],'subroutine']
		#print
		self.do_at_beginning_of_subroutine(tok[2:])
		#l_tok2=self.add_use__in_ol()
		return [0,'']

	if 'function' in tok:
	    def_func=0
	    for i in range(len(tok)):
		if tok[i]=='function':
			j=i-1
			if tok[j]=='integer' or tok[j]=='precision' or tok[j]=='logical':
				#self.do_at_end_of_subroutine()
				#print "function",tok[i+1]
				self.execfuncname=[tok[i+1],'function']
				def_func=1
				break
	    if def_func==1:
		# function=i foo=i+1 (=i+2
		self.do_at_beginning_of_subroutine(tok[i+2:])
		#l_tok2=self.add_use__in_ol()
		if i+1 >= len(tok):
			print thisfunc,"ERROR try to access tok[i+1], i+1=",i+1, "len=",len(tok)
			print thisfunc, linenumber,tok 
			print >>sys.stderr, thisfunc,"ERROR try to access tok[i+1], i+1=",i+1, "len=",len(tok)
			print  >>sys.stderr,thisfunc, linenumber,tok 
			# don't exit to know linenumber 
			#sys.exit(10) 
		return [0,'']

	if tok[0]=='use':
		self.do_use__il_on(tok)
			
	if tok[0]=='integer' or tok[0]=='double' or tok[0]=='logical' or tok[0]=='real' or  tok[0]=='type' or tok[0]=='character' :
		varlist=self.do_typevar__il_on(tok)
		# make dictionary
		for var1 in varlist:
			self.var_dic[var1[0]]=[var1[1],var1[2]]

		l_tok2=self.manip_varlist(varlist)
		#print thisfunc,"typedef",l_tok2
		return l_tok2
	if tok[0]=='external' or tok[0]=='data' or tok[0]=='parameter':
		return [0,[]]

	if len(tok)>=2 and tok[0]=='end' and (tok[1]=='module' or tok[1]=='subroutine'):
		return [0,[]]


	# in keyword
	have_name=0
	for name in self.keyword.list:
		if name in tok:
			have_name=1

	if len(fnlabel)>0:
		have_name=1

	ichanged=0
	if have_name==1:
		l_tok=self.do_function(tok,linenumber,fnlabel)
		ichanged,tok= l_tok


	ichanged_element,l_tok=self.do_element_in_type(tok,linenumber,fnlabel)
	tok=l_tok
	if ichanged_element>0:
		ichanged=1


	self.analyze_wvar(tok,linenumber,fnlabel)


	iloop=1
        while iloop>0:
            iloop=0
            for w_var in self.w_pointerlist:
                if w_var.name in tok and self.execfuncname[0]==w_var.funcname :
                        if self.w_varlist_contains(w_var.name,self.execfuncname[0])==0 and self.printsource==1:
                            if tok[0]=="common":
                                print
                                print "#error w variable,",w_var.name,"is defined in common block"
                                print
                                print >> sys.stderr, "ERROR w variable,",w_var.name,"is defined in common block"
                        # check w_var[0] is in self.w_varlistundel
                        id= self.w_varlistundel_contains(w_var.name,self.execfuncname[0])
			if id>=0:
				continue
                        nn = self.w_pointerlist_find( w_var.name,self.execfuncname[0] )
			if isinstance(nn,t_var7)==0:
				continue
                        if self.printsource==1:
                        	print thisprogram,"POINTER, change ","[",nn,"]",w_var.name,self.execfuncname[0]

                       	tok,iloop0=self.do_w_var(w_var, tok,linenumber,fnlabel)
                       	ichanged=ichanged+iloop0
                       	iloop=iloop+iloop0

        iloop=1
	while iloop>0:
            iloop=0
	    for w_var in self.w_varlistall:
		ipos_w_var=-1
		for i in range(len(tok)):
			if w_var.name== tok[i]:
				ipos_w_var = i
				break;
		if ipos_w_var>=0 and self.execfuncname[0]==w_var.funcname :
			if self.w_varlist_contains(w_var.name,self.execfuncname[0])==0 and self.printsource==1:
			    if self.printsource==1:
				print thisprogram,'warning(1) , probably ',w_var.name,' is not defined yet at linenumber=',linenumber
			    if tok[0]=="common":
				print
				print "#error w variable,",w_var.name,"is defined in common block"
				print
				print >> sys.stderr, "ERROR w variable,",w_var.name,"is defined in common block"
			# check w_var[0] is in self.w_varlistundel
			id= self.w_varlistundel_contains(w_var.name,self.execfuncname[0])
			if id>=0:
				# check w ( foo ...
				if ipos_w_var-2>=0:
					if tok[ipos_w_var-2] in ["w","defdr","defrr","defr","defcc",'defdc',"defi"] and tok[ipos_w_var-1]=="(":
						continue
				if '=' in tok:
					if self.printsource==1:
						print thisprogram,"eq",w_var.name,tok
					self.w_varlistundel_append_uniq(t_reason(['eq',w_var.name,w_var.funcname]),linenumber,tok)
					continue
				if self.printsource==1:
					print thisprogram,"pass",w_var.name,tok
				self.w_varlistundel_append_uniq(t_reason(['pass',w_var.name,w_var.funcname]),linenumber,tok)
				if self.printsource==1:
					print thisprogram,'(1)do not change',w_var.name,'because of',self.w_varlistundel[id].reason
				continue

			tok,iloop0=self.do_w_var(w_var,tok,linenumber,fnlabel)
			ichanged=ichanged+iloop0
			iloop=iloop+iloop0

	if ichanged>0:
		return [1,tok]
	else:
		return [0,'']




#----------------------------------------------------------

def print_original__isl_on(pad,packline):
	line=packline[0]
	if isinstance(line[0],list)==1:
	    for s in line:
		if len(pad)==0:
			print s[0]
		else:
			print pad,s[0]
	else:
		s=line
		if len(pad)==0:
			print s[0]
		else:
			print pad,s[0]

	
def make_line_recur__ili_ol(lines,level=0):
	thisfunc='make_line_recur__ili_ol'
	nmax=60
	cmd2=''
	cmd=[]
	first=0
	for line in lines:
		if isinstance(line,basestring)==1:
			if level==0:
				cmd2=cmd2+line+' '
				if len(cmd2)>nmax:
					cmd.append([first,cmd2])
					cmd2=''
					first=1
			else:
				cmd.append([first,line])
				first=0
		elif isinstance(line,list)==1:
			if len(cmd2)>0:
				cmd.append([first,cmd2])
				cmd2=''
				first=0
			cmdnew=make_line_recur__ili_ol(line,level+1)
			for n2 in cmdnew:
				cmd.append(n2)
			first=0
			cmd.append('')
	if len(cmd2)>0:
		cmd.append([first,cmd2])
	return cmd


def make_list_recur__ili_ol(nspc,list,fnlabel=''):
	thisfunc='make_list_recur__ili_ol'
	if nspc<6: 
		nspc=6
	tok=[]
	pad=''
	for i in range(nspc):
		pad=pad+' '
	if len(fnlabel)>0:
		pad1 =fnlabel.ljust(nspc-1)
		pad1=' '+pad1
	else:
		pad1=pad

	padcont=''
	for i in range(nspc):
		if i==5:
			padcont=padcont+'.'
		else:
			padcont=padcont+' '
	cmd=make_line_recur__ili_ol(list)
	ifirst=1
        for name in cmd:
	    if len(name)>0:
		if name[0]==0:
                	#print pad+name[1]
		    if ifirst==1:
			tok.append(pad1+name[1])
			ifirst=0
		    else:
			tok.append(pad+name[1])
		else:
			#print padcont+name[1]
			tok.append(padcont+name[1])
	#print
	tok.append('')
	return tok

def get_fnlabel(line,printit=0):
	thisfunc="get_fnlabel"
	if printit==1:
		print thisfunc,"line=",line,isinstance(line,list)
	if isinstance(line,list)==1:
		aline=line[0]
		if isinstance(aline,list)==1:
			aline=aline[0]
	elif isinstance(line,basestring)==1:
		aline=line
	else:
		print thisfunc,"ERROR"
		print >>sys.stderr,thisfunc,"ERROR"
		sys.exit(10)
	if printit==1:
		print thisfunc,"aline=",aline
	if isinstance(aline,basestring)==0:
		print thisfunc,"type ERROR",aline
		print >>sys.stderr,thisfunc,"type ERROR",aline
		sys.exit(10)
		
	if len(aline)>5:
		a=''
		if aline[0]==' ':
			str=aline[1:5]
			a= str.strip()
		elif aline[0].isdigit():
			str=aline[0:5]
			a= str.strip()
		if len(a)>0:
			return a
	return ''

def get_linenumber(line):
        thisfunc="get_linenumber"
        if isinstance(line[0],list)==1:
                aline=line[0]
		aline= aline[1]
        elif isinstance(line[0],basestring)==1:
                aline=line[1]
        else:
                print thisfunc,"ERROR"
                print>>sys.stderr, thisfunc,"ERROR"
                sys.exit(10)
	return aline


def analyze_option0(packline,printsource=0):
    thisfunc='analize_option0'
# CDELW savewlist name
# CDELW loadwlist name
    if len(packline)<=15:
	return []
    if packline[0:6]!='CDELW2':
	return []
    s=packline.split(' ')
    s2=[]
    for n in s:
	if n=='':
		continue
	s2.append(n)
    s=s2
    if len(s)>=3:
#	opt,cmd,name=s[0:3]
#	if opt=='CDELW1':
#		if cmd=='savewlist' or cmd=='loadwlist':
#			#print thisfunc,cmd,name
#			return [cmd,name]
	opt=s[0]
	if opt=='CDELW2':
		if s[1]=='allocate':
			name=s[2]
			if s[3].isdigit()==1:
				dimension=int(s[3])
			else:
				print thisfunc,"option error",packline
				print >>sys.stderr, thisfunc,"option error",packline
				sys.exit(10)
			arraysize=s[4:]
			if len(arraysize)!=dimension:
				print thisfunc,"ERROR, dimension!=arraysize",packline
				print >>sys.stderr,thisfunc,"ERROR, dimension!=arraysize",packline
				sys.exit(10)
			return t_arraydef(s[1],name,arraysize)
		if s[1]=='deallocate':
			name=s[2]
			return t_arraydef(s[1],name,[])
		if s[1]=='activate':
			if printsource==1:
				print spc6,
				for nn in s[2:]:
					print nn,
				print 
			
    return []


def analyze_option(packline,printsource=0):
    thisfunc='analyze_option'
    ret=[]
    if isinstance(packline[0],list)==1:
	if isinstance(packline[0][0],basestring)==1:
		for nn in packline:
			ret2=analyze_option0(nn[0])
			if len(ret2)==2:
				return ret2
	else:
		print thisfunc,'ERROR, strange packline'
		print packline[0]
		print packline[0][0]
		sys.exit(10)
    elif isinstance(packline[0],basestring)==1:
	ret=analyze_option0(packline[0],printsource)
    return ret 

def analyze_cpp(packline,printsource):
	thisfunc='analyze_cpp'
	
	if isinstance(packline[0],list)==1:
        	a=re.match("#if.*F90",packline[0][0])
	elif isinstance(packline[0],basestring)==1:
        	a=re.match("#if.*F90",packline[0])
	else:
		print thisfunc,'ERROR, strange packline'
                print packline[0]
                print packline[0][0]
                sys.exit(10)
	if a!=None:
		if printsource==1:
			print "#error have #if F90 directive"




def do_main(contline,varfile0,mapfile0,fix,printsource0,printmap0):
    thisfunc='do_main'
#parse constructor
    doparse=DoParse(varfile0,mapfile0,printsource0,printmap0)


    doparse.w_varlistundel_fix=fix

#print
    printcont=1

    have_module=0
    #for packline in contline:
#	tok=line2tok__is_ol(packline[1])
#	if len(tok)>0 and tok[0]=='module':
#		have_module=1
#		break
#
#    if len(doparse.var_toadd)>0 and have_module==0 and doparse.self.printit==1:
#	print spc6+"module func_"+doparse.var_toadd[0][0]
#	print spc6+"contains"

    for packline in contline:
	if printcont==1:
		# parse mode
	    fnlabel=get_fnlabel(packline[0])
	    if len(packline[1])==0:
		fnlabel=''
	    linenumber=get_linenumber(packline[0])

	    analyze_cpp(packline[0],printsource0)

	    wopt=analyze_option(packline[0],doparse.printsource)
	    if len(wopt)==2:
		if wopt[0]=='loadwlist':
			doparse.load_wvar(doparse.execfuncname[0]+'.'+wopt[1])
		elif wopt[0]=='savewlist':
			doparse.save_wvar(doparse.execfuncname[0]+'.'+wopt[1])

	    listchanged=doparse.do_line__is_ol(packline[1],linenumber,fnlabel)

	    if doparse.printsource==1:
		if listchanged[0]==1:
			print_original__isl_on(thisprogram,packline)
	                nspc=count_spc__il_oi(packline[0])
       	         	cmd=make_list_recur__ili_ol(nspc,listchanged[1],fnlabel)
                    	for nn in cmd:
                        	print nn
		elif listchanged[0]==2:
			print_original__isl_on('',packline)
                        nspc=count_spc__il_oi(packline[0])
                        cmd=make_list_recur__ili_ol(nspc,listchanged[1])
                        for nn in cmd:
                                print nn
		elif listchanged[0]==3:
                        nspc=count_spc__il_oi(packline[0])
                        cmd=make_list_recur__ili_ol(nspc,listchanged[1],fnlabel)  # may cause error
                        for nn in cmd:
                                print nn
                        print_original__isl_on('',packline)
		elif listchanged[0]==0:
			print_original__isl_on('',packline)


	elif printcont==2 and doparse.printsource==1:
		# print original line mode
	    for list in packline[0]:
		print "3>", list
	elif printcont==3 and doparse.printsource==1:
		print packline

#    if len(doparse.var_toadd)>0 and have_module==0 and doparse.self.printit==1:
#	print spc6+"end module func_"+doparse.var_toadd[0][0]
			
#print doparse.var_dic



#---main---
thisfunc="main"
#constructure
fline=mFLine.FLine()
contline=fline.get__in_ol()

try:
	os.unlink('newvar.dat')
except:
	idummy=0

do_main(contline,'newvar.dat','file.map',0,0,1)
do_main(contline,'newvar.dat','file.map',0,0,0)
do_main(contline,'newvar.dat','file.map',1,1,0)

