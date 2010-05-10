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


thisprogram='Cdelw1'
strinfo='...info...'+spc10


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
    indentmap=[]

    printmap=0
    mapfile='__file.map'

    callvarlist=[]

    idisable_change=0

    def __init__(self,varfile0,mapfile0,printsource0,printmap0):
	thisfunc='__init__'
	self.keyword=Keyword()
        self.struclist=mStrucList.StrucList()

	self.newvarfile=varfile0
	self.mapfile=mapfile0
    	self.printsource=printsource0
    	self.printmap=printmap0

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
	

    def __del__(self):
	    thisfunc="__del__"
	    try:
		self.newvarfile_desc.close()
	    except:
		print 'failed to close ',self.newvarfile
		print>>sys.stderr, 'failed to close ',self.newvarfile
		sys.exit(10)

    def save_wvar(self,name):
	thisfunc='save_wvar'
	try:
		fw=open('#delw.wlist.'+name,'w')
	except:
		print thisfunc,'failed to open',name
		sys.exit(10)

	for nn0 in self.w_varlist:
		nn=nn0
		for i in range(len(nn)):
			if nn[i]=='':
				nn[i]='***'
		s=''
		for n2 in nn:
			s=s+str(n2)+' '
		fw.write(s+'\n')

	fw.close()
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
		self.w_varlist.append(aline[0:7])
		
	fr.close()

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

	self.w_varlistall=[]
	self.w_varlistundel=[]
	for line in fr:
		aline= line[:-1].split(' ')
		#print thisfunc, aline
		if aline[0]=='var':
			var, name,level,nameorig,namesize,type,prevar,execfuncname = aline[0:8]
		elif aline[0]=='change':
			var,ok,execfuncname = aline[:3]
		elif aline[0]=='undel':
			#var,reason,name,execfuncname= aline[:4]
			var,name,execfuncname=aline[:3]
			listreason=aline[3:]
		else:
			print thisfunc,"ERROR ",aline
			print >>sys.stderr, thisfunc,"ERROR ",aline
			sys.exit(10)

		if var=='var':
			#self.w_varlistall.append([name,level,nameorig,namesize,type,prevar,execfuncname])
			self.w_varlistall_append_uniq([name,level,nameorig,namesize,type,prevar,execfuncname],0,[])
		elif var=='undel':
			#self.w_varlistundel.append([reason,name,execfuncname])
    			self.w_varlistundel_append_uniq(listreason,name,execfuncname,0,[])
			#if self.printsource==1:
			#	print thisprogram,"undel list",[name,reason,execfuncname]
	try:
		fr.close()
	except:
		print thisfunc,"failed to close",self.newvarfile


	if self.printsource==1:
	    for nn in self.w_varlistall:
		print thisprogram,"var",nn
	    for nn in self.w_varlistundel:
		print thisprogram,"undel", nn 


    def write_newvar(self):
# write newvar
	thisfunc='write_newvar'
	fw=self.newvarfile_desc
	if fw==0:
		print thisfunc,'ERROR newvarfile is not opened'
		print >>sys.stderr, thisfunc,'ERROR newvarfile is not opened'
		sys.exit(10)


	for aline in self.w_varlistall:
		#name,level,nameorig,namesize,type,prevar,execfuncname = aline
	    if aline[0] in self.callvarlist:
		name=aline[0]
		execfuncname=aline[6]
		fw.write('undel '+name+' '+execfuncname+' call\n')
		#print thisfunc,'undel call',name,execfuncname
	    else:
		linew='var '
		for nn in aline:
			if isinstance(nn,int)==1:
				linew=linew+str(nn)+' '	
			else:
				linew=linew+nn+' '	
		fw.write(linew+'\n')

	for aline in self.w_varlistundel:
		linew='undel '
		reason,name,execfuncname = aline
		#fw.write(linew+reason+' '+name+' '+execfuncname+'\n')
		fw.write(linew+name+' '+execfuncname+' '+list2string(reason,'','')+'\n')

        if self.idisable_change==1:
		ok='disable'
	else:
		ok='enable'
	#print 'change '+ok+ ' '+self.execfuncname[0]
	fw.write('change '+ok+ ' '+self.execfuncname[0] + '\n')


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
		a=self.w_varlist[i][0]
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


    def w_varlist_append_uniq(self,w_var,linen,tok,printit=0):
	thisfunc="w_varlist_append_uniq"
	#print thisfunc,len(w_var),w_var
	[name,level,nameorig,sizeorig,type,prevar,routinename] = w_var
	for w in self.w_varlist:
		if name==w[0] and routinename==w[6]:
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
        [name,level,nameorig,sizeorig,type,prevar,routinename] = w_var
        for w in self.w_varlistall:
                if name==w[0] and routinename==w[6]:
                        if type!=w[4] and type!='undef':
			    self.w_varlistundel_append_uniq(['redef'],name,routinename,linen,tok)
			    if printit==1:
                                print thisprogram,thisfunc,"ERROR, define again with another type:",name,w[4],"and",type
                                print>>sys.stderr, thisprogram,thisfunc,"ERROR, define again with another type:",name,w[4],"and",type
			#if printit==1:
                        #	print thisprogram,thisfunc,"warning, duplicated symbol:",name
                        return
        self.w_varlistall.append(w_var)
	if type=='undef':
		self.w_varlistundel_append_uniq(type,name,routinename,linen,tok,printit)

    def w_varlistundel_append_uniq(self,reason,name,routinename,linen,tok,printit=0):
        thisfunc="w_varlistundel_append_uniq"
        #[name,level,nameorig,sizeorig,type,prevar,routinename] = w_var

	for i in range(len(self.w_varlistundel)):
		if name==self.w_varlistundel[i][1] and routinename==self.w_varlistundel[i][2]:
			listreason=self.w_varlistundel[i][0]
			if isinstance(reason,list)==1:
				for areason in reason:
					listreason=list_adduniq(listreason,areason)
			elif isinstance(reason,basestring)==1:
				listreason=list_adduniq(listreason,reason)
			else:
				print thisfunc,"error, unknown reason type, reason=",reason
				sys.exit(10)
			self.w_varlistundel[i][0]=listreason
                        return
	if isinstance(reason,basestring)==1:
        	self.w_varlistundel.append([[reason],name,routinename])
	elif isinstance(reason,list)==1:
                self.w_varlistundel.append([reason,name,routinename])


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

    def get_prevar(self,type):
	thisfunc='get_prevar'
	prevar=''
	if type=='real(8)':
		prevar='rv_w_'
	elif type=='integer':
		prevar='iv_w_'
	elif type=='complex(8)':
		prevar='zv_w_'
	else:
		print thisfunc,'unknown type=<'+kind+'>'
		print >>sys.stderr, thisfunc,'unknown type=<'+kind+'>'
		sys.exit(10)
	return prevar


    def get_kind(self,callname):
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

	prevar=self.get_prevar(type)

	return [type,prevar]


    def make_statement_alloc(self,toks,alist,routinename):
        thisfunc='make_statement_alloc'
        [name,level,nameorig,namesize,type,prevar]= alist

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
                [name,level,nameorig,namesize,type,prevar,routinename]= alist
		id= self.w_varlistundel_contains(name,self.execfuncname[0])
		if id>=0:
			listrlse.append(name)
			if self.printsource==1:
				print thisprogram,"not deallocate",name+" because of "+list2string(self.w_varlistundel[id][0])
			continue
                realname=prevar+name
                tok='if (allocated('+realname+')) deallocate('+realname+')'
                tok2.append([tok])
	if len(listrlse)>0 and alsorlse==1:
		n=len(listrlse)-1
		tok='call rlse('+listrlse[n]+')'
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
                [name,level,nameorig,namesize,type,prevar,routinename]= alist
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
	nameorig='undef'
	for i in range(len(tok)-2):
		if tok[i]=="w" and tok[i+1]=='(' :
			icmax,var=get_var(tok[i:])
			namenext=var[3]
			# reason if there is 
			if namenext=='(':
				reason='call'
			else:
				reason='wref'
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
					self.w_varlistundel_append_uniq([reason],name,self.execfuncname[0],linenumber,tok,self.printsource)
					print >> sys.stderr ,  "ERROR , w(o+i) type access",list2string(var),"in",self.execfuncname[0]
			elif len(var)==4:
			# w ( i ) = 4 
				name=var[2]
				if name.isdigit()==1: 
			  	    if self.printsource==1:
					print thisprogram,name,"is digit, dropped"
				else:
				    self.w_varlistall_append_uniq([name,level,nameorig,size,type,prevar,self.execfuncname[0]],linenumber,tok,self.printsource)
			else:
				print thisfunc,"ERROR strange var=",var 
				sys.exit(10)



    def do_function(self,tok,linenumber,fnlabel):	
	thisfunc="do_function"
	ichanged=0


	if self.printit>9:
		print thisfunc,linenumber,tok,fnlabel
	if len(fnlabel)>0:
		self.indentmap.append([linenumber,['label',fnlabel],len(self.indentstack),fnlabel,tok])
	if tok[0]=="if":
		if "then" in tok:
			self.indentstack.append(['if',''])
			self.indentmap.append([linenumber,['if'],len(self.indentstack),fnlabel,tok])
	if tok[0]=="elseif" or (len(tok)>=2 and tok[0]=="else" and tok[1]=="if"):
		if "then" in tok:
			self.indentmap.append([linenumber,['elseif'],len(self.indentstack),fnlabel,tok])
	if tok[0]=="else":
		self.indentmap.append([linenumber,['else'],len(self.indentstack),fnlabel,tok])
	if 'endif' in tok or (len(tok)>=2 and tok[0]=="end" and tok[1]=="if"):
		self.indentstack_release(fnlabel,'if')
		self.indentmap.append([linenumber,['endif'],len(self.indentstack),fnlabel,tok])

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
				self.indentmap.append([linenumber,['do',tok[1]],len(self.indentstack),fnlabel,tok])
			else:
				self.indentstack.append(["do",''])
				self.indentmap.append([linenumber,['do',''],len(self.indentstack),fnlabel,tok])
		else:
			self.indentstack.append(["do",""])
			self.indentmap.append([linenumber,['do',''],len(self.indentstack),fnlabel,tok])
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
			self.indentmap.append([linenumber,['enddo',''],len(self.indentstack),fnlabel,tok])
	if len(fnlabel)>0:
		if self.printit>9:
			print thisfunc,'search do',fnlabel
                n=self.indentstack_release(fnlabel,'do')
                for i in range(n):
                        self.indentmap.append([linenumber,['enddo',fnlabel],len(self.indentstack),fnlabel,tok])
	if tok[0]=="select":
		self.indentstack.append(["select",fnlabel])
		self.indentmap.append([linenumber,['select'],len(self.indentstack),fnlabel,tok])
	if tok[0]=="endselect" or (len(tok)>=2 and tok[0]=="end" and tok[1]=="select"):
		self.indentstack.pop()
		self.indentmap.append([linenumber,['endselect'],len(self.indentstack),fnlabel,tok])
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
				nameorig,namesize=self.get_defname(l_tok,i+3)

				type,prevar=self.get_kind(l_tok[i+1])

				self.w_varlist_append_uniq([name,level,nameorig,namesize,type,prevar,self.execfuncname[0]],linenumber,tok,self.printsource)
				self.w_varlistall_append_uniq([name,level,nameorig,namesize,type,prevar,self.execfuncname[0]],linenumber,tok,self.printsource)
				self.indentmap.append([linenumber,[l_tok[i+1],name],level,fnlabel,tok])
				# new statement 
				if id<0:
				    if name==nameorig:
					l_tok=self.make_statement_alloc([l_tok,i],[name,level,nameorig,namesize,type,prevar],self.execfuncname[0])
					ichanged=1
					iloopagain=1
					#print thisfunc,"def",l_tok,ichanged
				    else:
					self.w_varlistundel_append_uniq('array',name,self.execfuncname[0],linenumber,tok,0)
					ichange_ERROR=1
			elif l_tok[i]=='call' and l_tok[i+1]=='rlse':
				#name=tok[i+3]
				name,nameorg=self.get_rlsename(l_tok,i+1)
				level=len(self.indentstack)
				# new statement
				if name==nameorg:
					list_rel=self.w_varlist_release(name,level,linenumber,self.printsource)
					self.indentmap.append([linenumber,[l_tok[i+1],name],level,fnlabel,tok])
					l_tok=self.make_statement_rlse_dealloc([l_tok,i],list_rel)
					ichanged=1
					iloopagain=1
					#print thisfunc,"rlse",l_tok,ichanged
				else:
					#if self.printsource==1:
					#	print "#error, ERROR do not change",nameorg,"because of rlse array"
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
	#	    if n>0:
	#	    	for i in range(n):
	#			self.indentmap.append([linenumber,['enddo',fnlabel],len(self.indentstack),fnlabel,tok])
	if "goto" in tok:
		label=self.get_gotolabel(tok)
		if self.printit>9:
			print "goto"
		for nn in label:
			self.indentmap.append([linenumber,['goto',nn],len(self.indentstack),fnlabel,tok])
			if self.printsource==1:
				#takao# print >>sys.stderr , "ERROR, possible spagettie code at line",linenumber,"label=", nn
				pass

	if "go" in tok:
		#search go and to
		for n in range(len(tok)-1):
			if tok[n]=="go" and tok[n+1]=="to":
				label=self.get_gotolabel(tok)
				for nn in label:
					self.indentmap.append([linenumber,['goto',nn],len(self.indentstack),fnlabel,tok])
					if self.printsource==1:
						#takao# print >>sys.stderr , "ERROR, possible spagettie code at",linenumber,"label=", nn
						pass
 				break
					
				
	if "return" in tok:
		for name in tok:
			if self.w_varlistall_contains(name,self.execfuncname[0]):
				print "#error tok contains return and",name
				print >>sys.stderr, "#error tok contains return and",name
				sys.exit(10)
		if self.printit>9:
			print "return"
		self.indentmap.append([linenumber,['return'],len(self.indentstack),fnlabel,tok])
		if self.printsource and len(self.w_varlist)>0:
			print "#error, have return with len(w_varlist)>0 at line",linenumber
			print >>sys.stderr ,"ERROR, have return with len(w_varlist)>0",linenumber
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
		id= self.w_varlistundel_contains(struc,self.execfuncname[0])
		if id>=0:
		    if self.printsource==1:
			print thisprogram,"do not change ",self.w_varlistundel[id]
		else:
		    for n2 in self.w_varlistall:
			name,level,nameorig,namesize,type,prevar,execfuncname0 = n2
			if (execfuncname0==self.execfuncname[0] or execfuncname0=='*') and  name==struc:
				#print thisfunc,"1>",var,var0,typestr,supp,self.execfuncname
				#print thisfunc,"2>",name,level,nameorig,namesize,type,prevar,execfuncname0
				aline = type+' ,allocatable :: '+prevar+name+'(:)'
				#print thisfunc,"3>",aline
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
	name=w_var[0]
	prefix=w_var[5]
	#print thisfunc,'first',w_var,tok
	tok0=tok
	
	# w(var) or var itself
	ichanged=1
	while ichanged==1:
	    ichanged=0
	    for i in range(len(tok)):
		if tok[i]==name:
			if i-2>=0 and tok[i-2]=='w' and tok[i-1]=='(':
				self.indentmap.append([linenumber,['wvar',w_var[0],w_var[1]],len(self.indentstack),fnlabel,tok])
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
				break
				
			else:
				self.indentmap.append([linenumber,['wref',w_var[0],w_var[1]],len(self.indentstack),fnlabel,tok])
				self.idisable_change=1
				self.w_varlistundel_append_uniq('wref',w_var[0],w_var[6],linenumber,tok,0)
				#print thisfunc,'wref',w_var,linenumber,tok
				return [tok0,0]

	#print thisfunc,'changed:',tok
	return [tok,1]


    def indentmap_deldup(self):
	thisfunc="indentmap_deldup"
	dellist=[]
	for i in range(len(self.indentmap)):
		if isinstance(self.indentmap[i],list)==1:
			linen,kind,level,fnlabel,tok= self.indentmap[i]
			if kind[0]=='wref':
				j=i-1
				if j>=0:
					linen2,kind2,level2,fnlabel2,tok2= self.indentmap[j]
					if kind2[0][0:3]=="def" and kind[1]==kind2[1]:
						if self.printit==1:
							print thisfunc,"delete ref",linen,tok
							print thisfunc,"delete this",linen2,tok2
						dellist.append(i)
	tok=[]
	for i in range(len(self.indentmap)):
		if i in dellist:
			continue
		tok.append(self.indentmap[i])

	self.indentmap=tok


    def indentmap_delunusedlabel(self):
	# search goto
	thisfunc="indentmap_delunusedlabel"
	gotolist=[]
	for i in range(len(self.indentmap)):
		if isinstance(self.indentmap[i],list)==1:
			linen,kind,level,fnlabel,tok= self.indentmap[i]
			if kind[0]=='goto':
				gotolabel=kind[1]
				gotolist.append(gotolabel)
	if self.printit==1:
		print thisfunc,"gotolist=",gotolist
	dellist=[]
	for aline in self.indentmap:
		linen,kind,level,fnlabel,tok= aline
		if kind[0]=="label":
		    if self.printit==1:
		    	print "label=",kind
		    ifound=0
		    for gotolabel in gotolist:
			label=kind[1]
			if gotolabel==label:
				if self.printit==1:
					print "found"
				ifound=1	
				break
		    if ifound==0:
			dellist.append(kind[1])
	if self.printit==1:
		print thisfunc,"dellabellist=",dellist

	# delete label
	newlist=[]
	for aline in self.indentmap:
                linen,kind,level,fnlabel,tok= aline
                if kind[0]=="label":
			if not (kind[1] in dellist):
				newlist.append(aline)
		else:
			newlist.append(aline)
	self.indentmap=newlist	
					

    def indentmap_delloop(self):
	thisfunc='indentmap_delloop'
	iflag=1
	while iflag==1:
		iflag=0
		for i in range(len(self.indentmap)):
			linen,kind,level,fnlabel,tok=self.indentmap[i]
			if kind[0]=='do':
				#search "enddo"
				j=i+1
				if j>=len(self.indentmap):
					print thisfunc, "index ERROR(do)",j,len(self.indentmap)
					print thisfunc,"indentmap=",self.indentmap
					print >>sys.stderr, thisfunc, "index ERROR(do)",j,len(self.indentmap)
                                        print >>sys.stderr, thisfunc,"indentmap=",self.indentmap
					sys.exit(10)
				linen2,kind2,level2,fnlabel2,tok2=self.indentmap[j]
				if kind2[0]=='enddo':
					#delete i and j
					if self.printit==1:
						print thisfunc,"delete1", self.indentmap[i]
						print thisfunc,"delete2", self.indentmap[j]
					del self.indentmap[i:j+1]
					iflag=1
					break
			if iflag==1:
				break
		for i in range(len(self.indentmap)):
		    if self.printit==1:
		    	print i,self.indentmap[i]
		    if isinstance(self.indentmap[i],list)==1:
                        linen,kind,level,fnlabel,tok=self.indentmap[i]
			idelete=0
                        if kind[0]=='if':
				#search endif
				for j in range(i+1,len(self.indentmap)):
					if self.printit==1:
						print thisfunc,"search j",self.indentmap[j]
					linen2,kind2,level2,fnlabel2,tok2=self.indentmap[j]
					if not (kind2[0] in ['else','elseif','endif']):
						break
					if kind2[0]=='endif':
						ing=0
						for k in range(i+1,j):
							kind3=self.indentmap[k][1]
							if self.printit==1:
								print thisfunc,"elseorendif?",self.indentmap[k],kind3
							if not (kind3[0]=='else' or kind3[0]=='elseif'):	
								ing=1
								break
						if self.printit==1:
							print thisfunc,"ing=",ing
						if ing==0:
							if self.printit==1:
							    for k in range(i,j+1):
								print thisfunc,"deleteN",self.indentmap[k]
							del self.indentmap[i:j+1]
							iflag=1
							idelete=1
							break
			if idelete==1:
				break
				
			
    def indentmap_checkindent(self):
	thisfunc="indentmap_checkindent"
	if self.printit==1:
		print thisfunc,len(self.indentstack),self.indentstack
	for kind in self.indentstack:
		if kind[0]=='if' or kind[0]=='do' or kind[0]=='select':
			print thisfunc, "indent level ERROR, remains:", self.indentstack
			print >> sys.stderr, thisfunc, "indent level ERROR, remains:", self.indentstack
			sys.exit(10)
			

    def indentmap_delloopjumpoutside0(self):
	thisfunc='indentmap_delloopjumpoutside0'
	for i in range(len(self.indentmap)):	
		linen,kind,level,fnlabel,tok=self.indentmap[i]
		if kind[0]=='goto':
			if self.printit==1:
				print thisfunc,'found goto',i,kind
			gotolabel=kind[1]
			# search label=gotolabel
			for j in range(len(self.indentmap)):
				linen2,kind2,level2,fnlabel2,tok2=self.indentmap[j]
				if kind2[0]=='label' and kind2[1]==gotolabel:
					if self.printit==1:
						print thisfunc,'found label',j,kind2
					if i<j:
						# jump backward
						i1=i
						i2=j
						blocklabel=-1
						for k in range(i2,i1,-1):
							linen3,kind3,level3,fnlabel3,tok3=self.indentmap[k]
							if kind3[0] in ['enddo', 'endif', 'endselect']:
								blocklabel=kind3[1]
								break
						if blocklabel>0:
							for k in range(len(self.indentmap)):
								linen4,kind4,level4,fnlabel4,tok4=self.indentmap[k]
								if kind4[0] in ['if', 'do', 'select']:
									if kind4[1]==blocklabel:
										# found range [k:i2]
										ing=0
										for L in range(k,i2):
											linen5,kind5,level5,fnlabel5,tok5=self.indentmap[L]
											if not kind5[0] in ['if','else','elseif','endif','do','enddo','select','endselect']:
												ing=1
										if ing==0:
											# yes you can delete
											del self.indentmap[k:i2+1]
											return 1
										
									
					else:
						# jump forward
						i1=j
						i2=i
						if self.printit==1:
							print thisfunc,"jump forward",i1,i2
                                                blocklabel=-1
                                                for k in range(i1,i2):
                                                        linen3,kind3,level3,fnlabel3,tok3=self.indentmap[k]
                                                        if kind3[0] in ['do', 'if', 'select']:
                                                                blocklabel=kind3[1]
								if self.printit==1:
									print thisfunc,'found ',blocklabel,k,kind3
                                                                break
                                                if blocklabel>0:
							if self.printit==1:
								print thisfunc,'serach label=',blocklabel
                                                        for k in range(len(self.indentmap)):
                                                                linen4,kind4,level4,fnlabel4,tok4=self.indentmap[k]
                                                                if kind4[0] in ['endif', 'enddo', 'endselect']:
                                                                        if kind4[1]==blocklabel:
										if self.printit==1:
											print thisfunc,'found end',k,kind4
                                                                                # found range [i1:k]
                                                                                ing=0
                                                                                for L in range(i1,k):
                                                                                        linen5,kind5,level5,fnlabel5,tok5=self.indentmap[L]
                                                                                        if not kind5[0] in ['if','else','elseif','endif','do','enddo','select','endselect','label','goto']:
                                                                                                ing=1
                                                                                if ing==0:
											if self.printit==1:
												print thisfunc,'you can delete',self.indentmap[i1:k+1]
                                                                                        # yes you can delete
                                                                                        del self.indentmap[i1:k+1]
											return 1
	return 0


    def indentmap_delloopjumpoutside(self):
	thisfunc='indentmap_delloopjumpoutside'
	iaction=1
	while iaction==1:
	    iaction=self.indentmap_delloopjumpoutside0()
					

    def indentmap_relabel(self):
	thisfunc='indentmap_relabel'
	ilabel=0
	indentstack=[]
	for i in range(len(self.indentmap)):
                linen,kind,level,fnlabel,tok=self.indentmap[i]
		if self.printit==1:
			print thisfunc,linen,kind
		if kind[0] in ['do','if','select']:
			ilabel=ilabel+1
			indentstack.append([kind[0],ilabel])	
			if self.printit==1:
				print thisfunc,"stack=",indentstack
			if len(kind)==2:
				newkind=[kind[0],ilabel,kind[1]]
			elif len(kind)==1:
				newkind=[kind[0],ilabel]
			self.indentmap[i]=[linen,newkind,level,fnlabel,tok]
		if kind[0] in ['else','elseif','case']:
			a=indentstack.pop()
			labelnow=a[1]
			indentstack.append(a)
			if self.printit==1:
				print thisfunc,"stack=",indentstack
			if len(kind)==2:
				newkind=[kind[0],labelnow,kind[1]]
			elif len(kind)==1:
				newkind=[kind[0],labelnow]
			self.indentmap[i]=[linen,newkind,level,fnlabel,tok]
		if kind[0] in ['endif','enddo','endselect']:
			a=indentstack.pop()
			if self.printit==1:
				print thisfunc,"stack=",indentstack
                        labelnow=a[1]
			if len(kind)==2:
				newkind=[kind[0],labelnow,kind[1]]
			elif len(kind)==1:
				newkind=[kind[0],labelnow]
                        self.indentmap[i]=[linen,newkind,level,fnlabel,tok]
	
	if len(indentstack)>0:
		print thisfunc,"ERROR indentstack is not empty"
		print thisfunc,indentstack
		print >>sys.stderr, thisfunc,"ERROR indentstack is not empty"
		print >>sys.stderr, thisfunc,indentstack
		sys.exit(10)
    

    def  w_varlist_addall(self):
	nn=''
        for name in self.w_varlist:
                nn=nn+ name[0]+' '
	return nn


    def w_varlist_show(self):
	if len(self.w_varlist)>0:
	    for name in self.w_varlist:
		print name[0],
	    print ''
	else:
	    print '(None)'


    def w_varlistall_addall(self):
	nn=''
        for name in self.w_varlistall:
                nn=nn+name[0]+' '
	return nn


    def w_varlistall_show(self):
	for name in self.w_varlistall:
		print name[0],
	print ''

    def w_varlistundel_addall(self):
        nn=''
        for name in self.w_varlistundel:
	    if name[2]==self.execfuncname[0]:
                nn=nn+name[1]+'('+list2string(name[0],'','')+') '

	for name in self.w_varlistall:
		if name[0] in self.callvarlist and name[2]==self.execfuncname[0]:
			nn=nn+name[0]+'(call) '
	return nn

    def w_varlistnotundel_addall(self):
	listundel=[]
	for name in self.w_varlistundel:
	    if name[2]==self.execfuncname[0]:
		listundel.append(name[1])
	errorlist=[]
	for name in self.w_varlist:
		if not name[0] in listundel and name[2]==self.execfuncname[0]:
			errorlist.append(name[0])
	nn=''
	if len(errorlist)==0:
		return nn
	for name in errorlist:
		nn=nn+name+' '
	return nn	

    def w_varlist_contains(self,name,funcname): 
	thisfunc='w_varlist_contains'
	for w_var in self.w_varlist:
		if name==w_var[0] and (w_var[6]==funcname or funcname=='*'):
			return 1
	return 0

    def w_varlistall_contains(self,name,funcname): 
	thisfunc='w_varlistall_contains'
	for w_var in self.w_varlistall:
		if name==w_var[0] and (w_var[6]==funcname or funcname=='*'):
			return 1
	return 0

    def indentstring(self,level):
	spc='    '
	s=''
        for i in range(level):
		s=s+spc
	return s	


    def indentmap_print(self,iop=1):
        thisfunc='indentmap_print'
	if self.printmap==0:
		return

	try:
		fileout=open(self.mapfile,'a')
	except:
		print 'ERROR failed to open file', self.mapfile
		print >>sys.stderr , 'ERROR failed to open file', self.mapfile,' to append'
		sys.exit(10)

	if iop==0:
	    thisfunc='indentmap_print'
            #print "indentmap"
	    if len(self.indentmap)>0:
	    	fileout.write("indentmap\n")
            if len(self.indentmap)==0:
                #print "None"
		fileout.write("None\n")
                return
            spc="   "
            for name in self.indentmap:
                linen,kind,level,fnlabel,tok = name
                #print spc,linen,kind,level,fnlabel
		fileout.write(spc+" "+linen+" "+kind+" "+str(level)+" "+fnlabel)

	elif iop==1:
	    thisfunc='indentmap_print_w_var'
	    if len(self.indentmap)>0:
	    	fileout.write("indentmap\n")
	    spc='   '
	    self.w_varlist=[]
	    self.w_varlistall=[]
	    for map in self.indentmap:
		linen,kind,level,fnlabel,tok= map
		#print spc,linen,kind,'indent=',level,
		if kind[0] in ['if','select','do','else']:
			spcindent=self.indentstring(level-1)
		else:
			spcindent=self.indentstring(level)
		fileout.write(spc+" "+str(linen).zfill(5)+" "+spcindent+list2string(kind)+" ")
		if kind[0] in ['defrr','defdr','defr','defi','defdc','defcc','redfi','redfrr']:
			type,prevar=self.get_kind(kind[0])

			w_var=kind[1]
			idx=tok.index(w_var)
			nameorig,namesize=self.get_defname(tok,idx)
			#print 'orig=',nameorig,'size=',namesize
			fileout.write('name='+nameorig+' size='+namesize+"\n")
			self.w_varlist_append_uniq([w_var,level,nameorig,namesize,type,prevar,self.execfuncname[0]],linen,tok,0)
			self.w_varlistall_append_uniq([w_var,level,nameorig,namesize,type,prevar,self.execfuncname[0]],linen,tok)
		elif kind[0]=='rlse':
			#print tok
			#print 'old w_varlist=',
			fileout.write('\nold w_varlist=')
			listo=self.w_varlist_addall()
			fileout.write(listo+"\n")
			w_var=kind[1]
			#print "release",w_var
			fileout.write("release "+w_var+"\n")
			list_rel=self.w_varlist_release(w_var,level,linen,0)
			#print 'new w_varlist=',
			fileout.write('new w_varlist=')
			listo=self.w_varlist_addall()
			fileout.write(listo+"\n")
		elif kind[0] in ['wref','wvar']:
			#print thisfunc,tok
			fileout.write('name='+kind[1]+'\n')
			if not self.w_varlist_contains(kind[1],self.execfuncname[0]):
				fileout.write('ERROR, probably '+kind[1]+' is not defined yet.\n')
		elif kind[0] in ['return']:
			fileout.write('\n')
			if len(self.w_varlist)>0:
                		#print "ERROR",
                		#print "w_varlist remains:",
                		fileout.write('ERROR, return but w_varlist remains:')
                		listo=self.w_varlist_addall()
                		fileout.write(listo+'\n')
                elif kind[0] in ['label']:
                        if len(self.w_varlist)>0:
                                fileout.write('w_varlist=')
                                listo=self.w_varlist_addall()
                                fileout.write(listo)
                        fileout.write('\n')
		elif kind[0] in ['goto']:
			if len(self.w_varlist)>0:
				fileout.write('w_varlist=')
				listo=self.w_varlist_addall()
                                fileout.write(listo)
			fileout.write('\n')
		elif kind[0] in ['loadwlist']:
			self.load_wvar(self.execfuncname[0]+'.'+kind[1])
                        if len(self.w_varlist)>0:
                                fileout.write('w_varlist=')
                                listo=self.w_varlist_addall()
                                fileout.write(listo)
			fileout.write('\n')
		elif kind[0] in ['savewlist']:
                        if len(self.w_varlist)>0:
                                fileout.write('w_varlist=')
                                listo=self.w_varlist_addall()
                                fileout.write(listo)
			fileout.write('\n')
		elif kind[0] in ['subroutine','function']:
			self.execfuncname=[kind[1],kind[0]]
			fileout.write('\n')
		else:
			#print
			fileout.write('\n')

	    if len(self.w_varlist)>0:
		#print "ERROR",
		#print "w_varlist remains:",
		fileout.write('end but w_varlist remains:')
		listo=self.w_varlist_addall()
		fileout.write(listo+'\n')

            listo=self.w_varlistundel_addall()
	    if len(listo)>0:
                fileout.write('w_varlistundel(reason) :')
                fileout.write(listo)
	        fileout.write('\n')

	    listo=self.w_varlistnotundel_addall()
	    if len(listo)>0:
		fileout.write('ERROR, w_varlist(except undel) : ')
		fileout.write(listo)
		fileout.write('\n')

	    #print ''
            fileout.write('\n')	
	    if len(self.w_varlistall)>0:
		#print "w_varlist used (all):",
		fileout.write('w_varlist used (all): ')
		listo=self.w_varlistall_addall()
		fileout.write(listo+'\n')
		fileout.write('to deallcate all\n')
		ltok=self.make_statement_dealloc(self.w_varlistall)
		for nn in ltok:
			fileout.write(spc6+nn[0]+'\n')

        fileout.write('\n')
	fileout.close()

    def indentmap_clearifnow(self):
	thisfunc='indentmap_clearifnow'
	iw=0
	for map in self.indentmap:
		linen,kind,level,fnlabel,tok= map
		if kind[0] in ['wval','wref','rlese']:
			iw=1
			break
		if kind[0][0:3]=='def':
			iw=1
			break
		if kind[0][0:4]=='redf':
			iw=1
			break
	if iw==0:
		self.indentmap=[]

    def indentmap_simplify(self):
	thisfunc="indentmap_simplify"
	#self.indentmap_print()
	self.indentmap_deldup()
	self.indentmap_delunusedlabel()
	self.indentmap_delloop()
	self.indentmap_checkindent()
	self.indentmap_relabel()
	self.indentmap_delloopjumpoutside()
	self.indentmap_checkindent()
	self.indentmap_clearifnow()

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
	    if n2[2] == self.execfuncname[0]:
                nn.append(n2[1])
        if len(nn)>0 and self.printsource==1:
                print thisprogram,"w_varlistundel:",list2string(nn)
        # make list in self.w_varlist ,but not nn
        listshoulddel=[]
	listshouldnotdel=[]
        for n2 in self.w_varlist:
                if not n2[0] in nn:
                        listshoulddel.append(n2[0])
		else:
			listshouldnotdel.append(n2[0])
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
        if len(self.indentmap)>0 and self.printmap==1:
                self.indentmap_simplify()
                self.indentmap_print()

	return l_tok

    def  do_at_beginning_of_subroutine(self,tok):
		thisfunc="do_at_beginning_of_subroutine"
                self.var_dic.clear()
                self.indentstack=[]
                self.w_varlist=[]
                #self.w_varlistall=[]
                #self.w_varlistundel=[]
                self.indentmap=[]
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
	#print thisfunc,'search', name,funcname	
	for i in range(len(self.w_varlistundel)):
		nn=self.w_varlistundel[i]
		#print thisfunc,'?',name,funcname,nn
		if nn[1]==name and nn[2]==funcname:
			#print thisfunc,name,funcname,'ret=',i
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
		self.indentmap.append([linenumber,['subroutine',tok[1]],len(self.indentstack),fnlabel,tok])
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
		self.indentmap.append([linenumber,['function',tok[i+1]],len(self.indentstack),fnlabel,tok])
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

	self.analyze_wvar(tok,linenumber,fnlabel)

	#search w_varlist in tok
        iloop=1
	#print thisfunc,"w_varlist=",
	#self.w_varlist_show()


	while iloop>0:
            iloop=0
	    for w_var in self.w_varlistall:
		if w_var[0] in tok and self.execfuncname[0]==w_var[6] :
			if self.w_varlist_contains(w_var[0],self.execfuncname[0])==0 and self.printsource==1:
			    if self.printsource==1:
				print thisprogram,'warning(1) , probably ',w_var[0],' is not defined yet at linenumber=',linenumber
			    if tok[0]=="common":
				print
				print "#error w variable,",w_var[0],"is defined in common block"
				print
				print >> sys.stderr, "ERROR w variable,",w_var[0],"is defined in common block"
			# check w_var[0] is in self.w_varlistundel
			id= self.w_varlistundel_contains(w_var[0],self.execfuncname[0])
			if id>=0:
				self.indentmap.append([linenumber,['wref',w_var[0],w_var[6]],len(self.indentstack),fnlabel,tok])
				self.w_varlistundel_append_uniq('wref',w_var[0],w_var[6],linenumber,tok)
				if self.printsource==1:
					print thisprogram,'do not change',w_var[0],'because of',self.w_varlistundel[id][0]
				continue
			#else:
			#	self.w_varlistundel_append_uniq('wref',w_var,linenumber,tok)
			#print thisfunc,"change", w_var
			tok,iloop0=self.do_w_var(w_var,tok,linenumber,fnlabel)
			ichanged=ichanged+iloop0
			iloop=iloop+iloop0

	if ichanged>0:
		return [1,tok]
	else:
		return [0,'']


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
		else:
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


def analyze_option0(packline):
    thisfunc='analize_option0'
# CDELW savewlist name
# CDELW loadwlist name
    if len(packline)<=15:
	return []
    if packline[0:6]!='CDELW1':
	return []
    s=packline.split(' ')
    if len(s)>=3:
	opt,cmd,name=s[0:3]
	if opt=='CDELW1':
		if cmd=='savewlist' or cmd=='loadwlist':
			#print thisfunc,cmd,name
			return [cmd,name]
    return []


def analyze_option(packline):
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
	ret=analyze_option0(packline[0])
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



def do_main(contline,varfile0,mapfile0,printsource0,printmap0):
    thisfunc='do_main'
#parse constructor
    doparse=DoParse(varfile0,mapfile0,printsource0,printmap0)



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

	    wopt=analyze_option(packline[0])
	    if len(wopt)==2:
		if wopt[0]=='loadwlist':
			doparse.load_wvar(doparse.execfuncname[0]+'.'+wopt[1])
			doparse.indentmap.append([linenumber,[wopt[0],wopt[1]],len(doparse.indentstack),fnlabel,''])
		elif wopt[0]=='savewlist':
			doparse.save_wvar(doparse.execfuncname[0]+'.'+wopt[1])
			doparse.indentmap.append([linenumber,[wopt[0],wopt[1]],len(doparse.indentstack),fnlabel,''])

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

    #print "last"
    #doparse.indentmap_simplify()
    #doparse.indentmap_print()



#---main---
thisfunc="main"
#constructure
fline=mFLine.FLine()
contline=fline.get__in_ol()

try:
	os.unlink('newvar.dat')
except:
	idummy=0

do_main(contline,'newvar.dat','file.map',0,1)
do_main(contline,'newvar.dat','file.map',0,0)
do_main(contline,'newvar.dat','file.map',1,0)

