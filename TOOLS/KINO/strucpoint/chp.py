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


thisprogram='Cchp1'
strinfo='...info...'+spc10

def get_kind(callname):
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

	return type

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

class t_var3:
    vartype=''
    prefix=''
    name=''

    def __init__(self,vartype,prefix,name):
	self.prefix=prefix
	self.name=name
	self.vartype=vartype

    def __str__(self):
	return self.vartype+','+self.prefix+','+self.name

    def realname(self):
	return self.prefix+self.name

class t_var4:
	fullname=''
	name=''
	vartype=''
	attrib=''
	def __init__(self,name,name_full,vartype,attrib):
		self.name=name
		self.fullname=name_full
		self.vartype=vartype
		self.attrib=attrib
	def __str__(self):
		return "["+self.name+','+self.fullname+','+self.vartype+','+self.attrib+"]"

class DoParse:
    printit=0
    varlist=[]

    def __init__(self):
	printit=0

    def change_var_in_tok(self,var,tok0):
	thisfunc='change_var_in_tok'

        ichanged=0
        tok=tok0
	iloop=1
	while iloop>0:
	    iloop=0
	    for i in range(len(tok)):
		if tok[i]==var.name:
			if i-2>=0 and i+1<len(tok) and tok[i-2]=='w' and tok[i-1]=='(':
			    if tok[i+1]==')':
				tok.pop(i+1)
				tok.pop(i)
				tok.pop(i-1)
				tok.pop(i-2)
				tok.insert(i-2,var.realname())	
				iloop=1
				ichanged=1
				break
			    else:
				print "#error",thisfunc,"type is not w(foo)",tok
				sys.exit(10)

	return [ichanged,tok]

    def typesize(self,tok,istart0):
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


    def do_typevar(self,tok):
        thisfunc='do_typevar'
        typestr0=''
        typesize='0'
        istart=-1
        ERROR=0
        if tok[0]=='type':
                istart=1
                typesize=-1     
                typestr0=tok[0]
                l_out=self.typesize(tok,istart)
                if l_out[0]:
                        istart=l_out[1]
                        typesize=l_out[2]
                else:
                        ERROR=1
        elif tok[0]=='integer' or tok[0]=='real':
                istart=1
                typesize='0'
                typestr0=tok[0]
                l_out=self.typesize(tok,istart)
                if l_out[0]:
                        istart=l_out[1]
                        typesize=l_out[2]
                else:
                        ERROR=1
        elif tok[0]=='character':
                istart=1
                typesize='0'
                typestr0=tok[0]
                l_out=self.typesize(tok,istart)
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
                l_out=self.typesize(tok,istart)
                if l_out[0]:
                        istart=l_out[1]
                        typesize=l_out[2]
                else:
                        ERROR=1
        elif tok[0]=='complex':
                istart=1
                typesize='4'
                typestr0='complex'
                l_out=self.typesize(tok,istart)
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
                alist=t_var4(ll[0],ll[1],typedef,supp)
                tok2.append(alist)

        return tok2


    def manip_defvarlist(self,defvarlist):
	thisfunc='manip_defvarlist'
        defnamelist=[]
	for nn in defvarlist:
		defnamelist.append(nn.name)


	newvarlist=[]
	ichanged=0
	for var in self.varlist:
		if var.name in defnamelist:
			#delete var in defnamelist
			idx=defnamelist.index(var.name)
			defvarlist.pop(idx)
			newvarlist.append(t_var4(var.prefix+var.name,var.prefix+var.name,var.vartype,",pointer"))
			ichanged=1
			

	if len(newvarlist)==0:
		return [0,[]]

	toklist=[]
	tok=[]

	if ichanged==1:
		if len(defvarlist)>0:
			defstr=defvarlist[0].vartype+defvarlist[0].attrib+"::"
			tok.append(defstr)

		for i in range(len(defvarlist)):
			nn=defvarlist[i]
			tok.append(nn.fullname)
			if i!=len(defvarlist)-1:
				tok.append(",")
		toklist=tok
		for nn in newvarlist:
			tok=nn.vartype+nn.attrib+" :: "+nn.fullname+"(:)"
			toklist.append([tok])
		
	return [ichanged,toklist]
				
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


    def manip_eq(self,tok):
	thisfunc='manip_eq'
	pos_eq=-1
        for  i in range(len(tok)):
        	if tok[i]=='=':
                	pos_eq=i
                        break

	if pos_eq==-1:
		# no eq , return
		return [0,tok]

	tok_perc=do_split0__ils_ol(tok,'%')
        ichanged=0
	for v in self.varlist:
		if v.name in tok_perc:
			prefix=v.prefix
			l_tok=self.element_connect(tok)
			if l_tok[0].find('%')>=0:
			 # source =l_tok[0]
                                # target= l_tok[2]
				var = separate_var([l_tok[2]])
				struc = separate_var([l_tok[0]])
				eqtype="sv"
			elif l_tok[2].find('%')>=0:
                                var = separate_var([l_tok[0]])
                                struc = separate_var([l_tok[2]])
                                eqtype="vs"
			else:
				# not supported
				return [0,tok]
			varstr=prefix+list2str(var.tok,"")
			if var.name==v.name:
				if var.array[0]==0 and var.array[1]==0:
					varstr=v.realname()
				else:
					print "#error",var.name, "is array."
					return [0,tok]
			tok2=[]
			for nn in struc.tok:
				if nn==var.name:
					tok2.append(prefix)
					tok2.append(nn)
				else:
					tok2.append(nn)
			strucstr=list2str(tok2,"")
			if v.name==struc.element:
				if struc.array[0]==0 and struc.array[1]==0:
					strucstr=struc.name+"%"+v.realname()
				else:
					print "#error",struc.name,"%",struc.element,"is array." 
					return [0,tok]
			if eqtype=="vs":
				ichanged=1
				tok=[varstr,"=>",strucstr]
			elif eqtype=="sv":
				ichanged=1
				tok=[strucstr,"=>",varstr]

	return [ichanged,tok]


    def manip_redf(self,var,tok):
	thisfunc="manip_redf"
        l=len(tok)
	toklist=[]
	if tok[0]=='call' and (tok[1]=='redfi' or tok[1]=='redfrr') and tok[2]=='(' and tok[l-1]==')':
		type=tok[1]
		ok=0
		type=get_kind(tok[1])
		if type==var.vartype:
			ok=1
			if type=='real(8)':
				vartmp="rv_a_tmp"
			elif type=='integer':
				vartmp="iv_a_tmp"
			else:
				print "#error",thisfunc, "unknown type,",vartmp, "no change"
		if ok==0:
			print "#error",thisfunc, "type mismatch, no change"
			return [0,tok]
		toksize=tok[5:l-1]
		sizestr=list2str(toksize,'')
		tok2="i_data_size=size("+var.realname()+"); allocate("+vartmp+"(i_data_size))"
		toklist.append([tok2])
		tok2=vartmp+"="+var.realname()+"; deallocate("+var.realname()+")"
		toklist.append([tok2])
		tok2="i_data_size=min(i_data_size,"+sizestr+"); allocate("+var.realname()+"("+sizestr+"))"
		toklist.append([tok2])
		tok2=var.realname()+"(:i_data_size)="+vartmp+"(:i_data_size); deallocate("+vartmp+")"
		toklist.append([tok2])
		return [1,toklist]
		
	return [0,tok]
		
    def manip_def(self,var,calldeflist,tok):
	thisfunc='manip_def'
	toklist=[]
	ipos=-1
	for i in range(len(tok)):
		if tok[i] in calldeflist:
			ipos=i
			break

	if ipos==-1:
		print "#error, failed to find",calldeflist,"in",tok
		return [0,tok]

	# size
	istart=ipos
	ic=0
        for i in range(istart,len(tok)):
                if tok[i]=='(':
                        ic=ic+1
                elif tok[i]==')':
                        ic=ic-1
                        if ic==0:
                                iend=i
                                break

	#print "tok[istart:iend]=",tok[istart-1:iend+1]

	toklist=tok[:istart-1]
	if tok[0]=='if':
		toklist.append("then")
        if tok[istart+2]==var.name and tok[istart+3]=="," and tok[len(tok)-1]==')':
		sizetok=tok[istart+4:len(tok)-1]
		sizestr=list2str(sizetok,'')
		tok2="allocate("+var.realname()+"(abs("+sizestr+")))"
		toklist.append([tok2])
		if var.vartype=='integer':
			zerostr="0"
		else:
			zerostr="0.0d0"
		tok2="if ("+sizestr+"<0) "+ var.realname()+"(:)="+zerostr
		toklist.append([tok2])
		if tok[0]=='if':
			toklist.append("endif")
		return [1,toklist]
	else:
		print "#error",thisfunc,"unknown structure",var.name
	return [0,tok]

    def manip_rlse(self,var,calldeflist,tok):
	thisfunc='manip_rlse'
	toklist=[]
	ipos=-1
	for i in range(len(tok)):
		if tok[i] in calldeflist:
			ipos=i
			break

	if ipos==-1:
		print "#error, failed to find",calldeflist,"in",tok
		return [0,tok]

	# size
	istart=ipos
	ic=0
        for i in range(istart,len(tok)):
                if tok[i]=='(':
                        ic=ic+1
                elif tok[i]==')':
                        ic=ic-1
                        if ic==0:
                                iend=i
                                break


	toklist=tok[:istart-1]
	if tok[0]=='if':
		toklist.append("then")
        if tok[istart+2]==var.name and tok[istart+3]==")" and tok[len(tok)-1]==')':
		tok2="if (associated("+var.realname()+")) deallocate("+var.realname()+")"
		toklist.append([tok2])
		if tok[0]=='if':
			toklist.append("endif")
		print "#error",thisfunc,"deallocate",var.realname(),", check whether this is OK or not."
		return [1,toklist]
	else:
		print "#error",thisfunc,"unknown structure",var.name
	return [0,tok]



	

    def do_line(self,line,linenumber,fnlabel):
	thisfunc='do_line'
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

	# var definition
	if tok[0] in ['integer' ,'double','logical','real','type' ,'character','complex'] :
		defvarlist=self.do_typevar(tok)
#
		ichanged,l_tok2=self.manip_defvarlist(defvarlist)
		if ichanged==1:
			# changed
			return [1,l_tok2]
		else:
			return [0,tok]

	if tok[0] in ['parameter','format','use','subroutine','end','enddo','continue']:
		return [0,tok]
		
	# foo = struc%goo 
	ichanged,tok2=self.manip_eq(tok)
	if ichanged==1:
		return [ichanged,tok2]

	# call defps2(...)
	if len(tok)>2 and tok[0]=='call' and tok[1]=='defps2':
		for var in self.varlist:
			if var.name in tok:
				print "#error found defps2, no change"
				return [0,tok]

	# call def*
	ichange=0
        calldeflist=['defrr','defdr','defr','defi','defcc','defdc']
	for var in self.varlist:
	    ifound=0
	    if var.name in tok:
	      for nn in tok:
		if nn in calldeflist:
			ifound=1
			break
	    if ifound==1:
		ichange2,tok=self.manip_def(var,calldeflist,tok)
		if ichange2==1:
			ichange=1
	if ichange==1:
		return [ichange,tok]

        ichange=0
        calldeflist=['rlse']
        for var in self.varlist:
            ifound=0
            if var.name in tok:
              for nn in tok:
                if nn in calldeflist:
                        ifound=1
                        break
            if ifound==1:
                ichange2,tok=self.manip_rlse(var,calldeflist,tok)
                if ichange2==1:
                        ichange=1
        if ichange==1:
                return [ichange,tok]

	    

	# call redf*
	if len(tok)>2 and tok[0]=='call' and (tok[1]=='redfrr' or tok[1]=='redfi'):
	    ichange=0
	    for var in self.varlist:
		if var.name==tok[3]:
			ichange2,tok=self.manip_redf(var,tok)
			if ichange2==1:
				ichange=1
	    if ichange==1:
		return [ichange,tok]

	# change w(foo) -> prefix+foo
        ichanged=0
	for var in self.varlist:
		if var.name in tok:
			# try to change it
			ichanged2,tok=self.change_var_in_tok(var,tok)
			if ichanged2>0:
				ichanged=1
			break
		# check var in self.varlist end 

	# not w(foo) -> error
	for var in self.varlist:
		for i in range(len(tok)):
			if tok[i]==var.name:
				if i-2>=0 and i+1<len(tok) and tok[i-2]=='w' and tok[i-1]=='(' and tok[i+1]==')':
					donothing=1
				else:
					print "#error",thisfunc,"unknown structure",tok
					sys.exit(10)

	return [ichanged,tok]

def do_main(contline,cmdvar):
    thisfunc='do_main'
    doparse=DoParse()

    doparse.varlist.append(cmdvar)
    for packline in contline:
            fnlabel=get_fnlabel(packline[0])
            if len(packline[1])==0:
                fnlabel=''
            linenumber=get_linenumber(packline[0])
	    ichanged,tok=doparse.do_line(packline[1],linenumber,fnlabel)
	    if ichanged==0:
	    	print_original__isl_on('',packline)
	    elif ichanged==1:
		print_original__isl_on(thisprogram,packline)
		nspc=count_spc__il_oi(packline[0])
		cmd=make_list_recur__ili_ol(nspc,tok,fnlabel)
                for nn in cmd:
                	print nn


#---main---
thisfunc="main"
argv=sys.argv
argc=len(argv)
#print argc,argv
if argc!=4:
	print "error arg must be 3"
	sys.exit(10)
cmdvar=t_var3(argv[1],argv[2],argv[3])
print thisprogram,"cmdvar=",cmdvar
#constructure
fline=mFLine.FLine()
contline=fline.get__in_ol()


do_main(contline,cmdvar)

