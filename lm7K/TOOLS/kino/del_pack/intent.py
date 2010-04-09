import sys
from types import *
import string

import mKeyword
import mFLine
import mStrucList

from mToken import *

spc10='          '
spc30=spc10+spc10+spc10
#spc6='123456'
spc6='      '


thisprogram='Cgetarg'
strinfo='...info...'+spc10

#---------------------------------------------------------

class DoParse:
    execfuncname=['','']

    var_dic={};
    var_toadd=[]
    newvar=[]
    newvarfile='newvar.dat'

    keyword=[]
    struclist=[]
    uselist=[]
    use_toadd=[]

    reflist=[]

    printit=0

    def __init__(self):
	self.keyword=mKeyword.Keyword()
        self.struclist=mStrucList.StrucList()


    def read_newvar(self):
# store in var_toadd
	thisfunc='read_newvar'
	try:
		fr=open(self.newvarfile,'r')
	except:
		return 
	for line in fr:
		z,a,b,c=line[:-1].split(' ')
		if z=='var':
			aline=[a,b,c]
			self.var_toadd.append(aline)
			if self.printit==1:
				print thisprogram+strinfo+" structure",aline
		elif z=='use':
			aline=[a,b]
			self.use_toadd.append(aline)
			if self.printit==1:
				print thisprogram+strinfo+" use_to_add",aline
	fr.close()


    def write_newvar(self):
# write newvar
	thisfunc='write_newvar'
	fw=open(self.newvarfile,'w')
	for line in self.newvar:
		aline='var '+line[0]+' '+line[1]+' '+line[2]
		fw.write(aline+'\n')

	flist=[]
	for line in self.newvar:
		module=line[0]
		flist.append(module)
	flist=list(set(flist))

	# flist must use 'm_struc_def' 
	# 
	mstrucdeflist=[]
	for nn in self.uselist:
		usemodule=nn[1]
		funcname=nn[0][0]
		if usemodule=='m_struc_def':
			mstrucdeflist.append(funcname)

	for f in flist:
		if not (f in mstrucdeflist):
			aline='use '+f+' m_struc_def dummy'
			fw.write(aline+'\n')
        fw.close()


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
					print 'error',thisfunc,'delm=; ', elem
				i=i+1
				alias=elem[i]
			if elem[i]==',':
				if a[1]==0:
					print 'error',thisfunc,'delm=, ', elem
				i=i+1
				mask=elem[i]
			if elem[i]==':':
				if a[1]==0:
					print 'error',thisfunc,'delm=: ', elem
				i=i+1
				range1=elem[i]
			if elem[i]=='.':
				if [1]==0:
					print 'error',thisfunc,'delm=. ', elem
				i=i+1
				range2=elem[i]
			i=i+1
		if range1!='' and range2=='':
			range2=range1
		tok2.append([name,alias,mask,range1,range2])
	return tok2

# 'a b c,d' -> [name,alias,mask,range1,range2] of a,b,c
#input string : str0
    def do_string__is_ol(self,str0):
	thisfunc="do_string__is_ol"
	n=len(str0)
	if str0[0]=="'" and str0[n-1]=="'":
		str=str0[1:n-1]
	else:
		print "error in",thisfunc
		print 'str0=',str0
		exit()	
		
        tok=self.do_stringoption__is_ol(str)
	return tok

#	tok=str.split(",")
#	have_mask=0
#	mask=""
#        if len(tok)==2:
#		have_mask=1
#		mask=tok[1]
#	cmd0=tok[0]
#	cmd=cmd0.split(" ")	 # cmd = list
#	n=len(cmd) 
#	print "string>",n,cmd,'mask=',mask
#
#	tok2=[]
#	tok2.append(have_mask)
#	tok2.append(mask)
#	tok2=tok2+cmd 
#	return tok2


    def add_struc(self,strstruc,struc):
        alist=[self.execfuncname[0],strstruc,struc]
	for nn in self.newvar:
		if alist==nn:
			return
	self.newvar.append(alist)

# manipulate [ dgets ( 'lat vol' , slat )  ]
# return = string
#input, list: tok
    def do_arg1__il_ol(self,tok):
	thisfunc='do_arg1__il_ol'
        n=len(tok)
#        if n!=6:
#                print "error in do_arg1__il_ol"
#                print 'n, tok=',n,tok
#                exit()
        funcname=tok[0]
        str=tok[2]
        struc=tok[4]
        tok1=self.do_string__is_ol(str)

	member_l=[]
	alias_l=[]
	mask_l=[]
	range1_l=[]
	range2_l=[]
	for elem in tok1:
		member_l.append(elem[0])
		alias_l.append(elem[1])
		mask_l.append(elem[2])
		range1_l.append(elem[3])
		range2_l.append(elem[4])

	if range1_l[1]!='':
		if range1_l[1]!=range2_l[1]:
			print 'error in range ', range1_l,range2_l
			sys.exit(10)

        mask=mask_l[1]
	have_mask=1
	if mask=='':
        	have_mask=0
        member=member_l
						#    struc,member 
        member_info = self.struclist.find__iss_os(member[0],member[1])
        member_size=member_info[2]
        member_type=member_info[0]
	if member_info[3]=='string':
		print thisfunc,"type error"
		sys.exit(10)


	self.add_struc(member[0],struc)

        tok2=struc+ "%"+ member[1]
	if range1_l[1]!='':
		tok2=tok2+'('+range1_l[1]+')'
	else:
		if member_size!='1':
			tok2=tok2+'(1)'
        if funcname=="igets" or funcname=="dgets":          
                tok2=tok2
	if funcname=='igets':
		if member_type[0]!='i':
			print "type mismatch", funcname,member_info
			sys.exit(10)
		tok2='int('+tok2+')'
		if have_mask==1:
			tok2='iand('+mask+','+tok2+')'

	elif funcname=='dgets':
		if member_type[0]=='i':
			print "type mismatch", funcname,member_info
			sys.exit(10)

        elif funcname=="lgors":
		if member_type[0]=='r':
			print "type mismatch", funcname,member_info
			sys.exit(10)
		tok2='int('+tok2+')'
	
		if have_mask==0:
			mask='0'
                if mask<0:
                        tok2=tok2+" .ne.0"
                else:
                        tok2="iand("+mask+","+tok2+") .ne.0"
        return tok2


# manipulate [ igetss ( 'spec lmxa' , is , sspec ) ] 
# return = string
#input,list: tok
    def do_arg2__il_ol(self,tok):
	thisfunc='do_arg2__il_ol'
        n=len(tok)
        if n!=8:
                print thisfunc,"error in do_arg2__il_ol",n,tok
                sys.exit(10)
        funcname=tok[0]
        str=tok[2]
        index=tok[4]
        struc=tok[6]
        tok1=self.do_string__is_ol(str)
        member_l=[]
        alias_l=[]
        mask_l=[]
        range1_l=[]
        range2_l=[]
        for elem in tok1:
                member_l.append(elem[0])
                alias_l.append(elem[1])
                mask_l.append(elem[2])
                range1_l.append(elem[3])
                range2_l.append(elem[4])

	for s in mask_l:
		if s!='':
			print 'mask is not null',mask_l,tok
			sys.exit(10)
	for s in range1_l:
		if s!='':
			print 'range1 is not null',range1_l,tok
			sys.exit(10)
	for s in range2_l:
		if s!='':
			print 'range1 is not null',range2_l,tok
			sys.exit(10)

        member=member_l
						# structure, member 
        member_info = self.struclist.find__iss_os(member[0],member[1])
        member_size=member_info[2]
        member_type=member_info[0]
	if member_info[3]=='string':
		print thisfunc,"error type"
		sys.exit(10)

	self.add_struc(member[0],struc)

        tok2=struc+'('+index+')%'+ member[1]
        if range1_l[1]!='':
                tok2=tok2+'('+range1_l[1]+')'
        else:
                if member_size!='1':
                        tok2=tok2+'(1)'

	if funcname=="igetss":
		if member_type=='real(8)':
			print 'type mismatch, ',funcname,member_type,member
			sys.exit(10)
	elif funcname=='dgetss':
		if member_type=='integer(8)':
			print 'type mismatch, ',funcname,member_type,member
			sys.exit(10)

        tok3=""
        if member[0]=="spec":
                if funcname=="igetss":
			if member_type[0]!='i':
                        	print "error, mismatch",funcname,member_type
                elif funcname=="dgetss":
			if member_type[0]!='r':
                        	print "error, mismatch",funcname,member_type
                tok3="int("+struc+"("+index+")%"+member[1]+")"
        elif member[0]=="site":
                if funcname=="igetss":
			if member_type[0]!='i':
                        	print "error, mismatch",funcname,member_type
                elif funcname=="dgetss":
			if member_type[0]!='r':
                        	print "error, mismatch",funcname,member_type
                tok3="("+struc+"("+index+")%"+member[1]+")"
        else:
                print "error in do_arg2__il_ol"
                exit()


        return tok3

	

# return = list
#input, list: tok
    def do_pack__il_ol(self,tok):
	thisfunc='do_pack__il_ol'
        n=len(tok)
        funcname=tok[0]
        str=tok[2]
	tok1=self.do_string__is_ol(str)
        member_l=[]
        alias_l=[]
        mask_l=[]
        range1_l=[]
        range2_l=[]
        for elem in tok1:
                member_l.append(elem[0])
                alias_l.append(elem[1])
                mask_l.append(elem[2])
                range1_l.append(elem[3])
                range2_l.append(elem[4])

        for s in mask_l:
                if s!='':
                        print 'mask is not null',mask_l,tok
                        sys.exit(10)
        for s in range1_l:
                if s!='':
                        print 'range1 is not null',range1_l,tok
                        sys.exit(10)
        for s in range2_l:
                if s!='':
                        print 'range1 is not null',range2_l,tok
                        sys.exit(10)

	str_struc=member_l[0]
        str_member=member_l[1:]
	for str in str_member:
		member_info=self.struclist.find__iss_os(str_struc,str)	

        struc=tok[4]
	var=["0","0","0","0","0"]

	self.add_struc(str_struc,struc)

        n=len(tok)

	i=0
	for ind in range(6,n,2):
		var[i]=tok[ind]
		i=i+1	
	ith=-1
	if str_struc=="site" or str_struc=="spec":
		ith=var[0]
		var=var[1:]
	cmd=[]
	cmd2=''
	cmd.append(cmd2)
	if funcname=='lsets':
		member_name=str_member[0]
		member_info=self.struclist.find__iss_os(str_struc,member_name)
		member_size=member_info[2]
		member_type=member_info[0]
		lval=var[0]
		lmask=var[1]
		if member_type[0]!='i':
			print thisfunc,"type mismatch"
			print str_struc,str_member,struc,var
			sys.exit(10)
		if str_struc=='site' or  str_struc=='spec':
			print thisfunc,'error , site or spec'
			print str_struc,str_member,struc,var
                        sys.exit(10)
		strucvar=struc+'%'+member_name
		cmd2="call lsets_bitop_i8("+strucvar+","+member_size+","+lval+","+lmask+")"
		cmd.append(cmd2)
			
		
	elif funcname[0]=="p": 
		i=0
		for name in str_member:
			member_info=self.struclist.find__iss_os(str_struc,name)
			member_size=member_info[2]
			member_type=member_info[0]
			if member_info[3]=='string':
				print thisfunc,"error type",member_info
				sys.exit(10)
			if member_type[0]=='i':
				callname='ii8copy'
			else:
				callname='dcopy'
			if is_variable(var[i])==1:
				var_type=[]
				try:
					var_type=self.var_dic[separate_array(var[i])]
				except:
					if self.printit==1:
						print thisprogram+strinfo,"NG::",var[i],",",member_type,"::",struc+"%"+name,member_size

				if len(var_type)>1:
				    if sametype(var_type[1],member_type)==0:
					print thisfunc,"type mismatch",
               	                 	print "type of var[i],",var[i], var_type
                                	print "type of string",member_type
					sys.exit(10)
				    if self.printit==1:
				    	print thisprogram+strinfo,var_type[1],"::",var[i],",",member_type,"::",struc+"%"+name,member_size

			strucvar=''
			if str_struc=="site" or str_struc=="spec":
				strucvar=struc+"("+ith+")%"+name
			else:
				strucvar=struc+"%"+name
			if member_size=='1':
                            if str_struc=="site" or str_struc=="spec":
			    	cmd2=strucvar+"="+var[i]+" "
			    else:
			    	cmd2=strucvar+"="+var[i]+" "

			    cmd.append(cmd2)
			else:
			    cmd2='i_copy_size=size('+strucvar+') '
			    cmd.append(cmd2)
			    cmd2='call '+callname+'(i_copy_size,'+var[i]+',1,'+strucvar+',1) '
			    cmd.append(cmd2)
			i=i+1
	elif funcname[0]=="u" :
		i=0
		for name in str_member:
			member_info=self.struclist.find__iss_os(str_struc,name)
			member_size=member_info[2]
			member_type=member_info[0]
			if member_info[3]=='string':
				print thisfunc,"error type", member_info
				sys.exit(10)
			if member_type[0]=='i':
				callname='i8icopy'
			else:
				callname='dcopy'
			if is_variable(var[i])==1:
				var_type=[]
				try:
					var_type=self.var_dic[separate_array(var[i])]
				except:
					if self.printit==1:
						print thisprogram+strinfo,"NG::",var[i],",",member_type,"::",struc+"%"+name,member_size
				if len(var_type)>1:
				    if sametype(var_type[1],member_type)==0:
					print thisfunc,"type mismatch",
               	                 	print "type of var[i],",var[i], var_type
                                	print "type of string",member_type
					sys.exit(10)
				    if self.printit==1:
				    	print thisprogram+strinfo,var_type[1],"::",var[i],",",member_type,"::",struc+"%"+name,member_size

                        strucvar=''
                        if str_struc=="site" or str_struc=="spec":
                                strucvar=struc+"("+ith+")%"+name
                        else:
                                strucvar=struc+"%"+name

			if member_size=='1':
			    cmd2=var[i]+"="+strucvar
			    cmd.append(cmd2)
			else:
			    if str_struc=="site" or str_struc=="spec":
				cmd2='i_copy_size=size('+strucvar+') '
			    else:
			    	cmd2='i_copy_size=size('+strucvar+') '
                            cmd.append(cmd2)
                            cmd2='call '+callname+'(i_copy_size,'+strucvar+',1,'+var[i]+',1) '

                            cmd.append(cmd2)
			i=i+1
	return cmd


# return = string
#input,list: tok
    def do_spackv__il_ol(self,tok):
	thisfunc='do_spackv__il_ol'
        n=len(tok)
        funcname=tok[0]
	job=int(tok[2])  # a new part
        str=tok[4]
        tok1=self.do_string__is_ol(str)
        member_l=[]
        alias_l=[]
        mask_l=[]
        range1_l=[]
        range2_l=[]
        for elem in tok1:
                member_l.append(elem[0])
                alias_l.append(elem[1])
                mask_l.append(elem[2])
                range1_l.append(elem[3])
                range2_l.append(elem[4])

        for s in mask_l:
                if s!='':
                        print 'mask is not null',mask_l,tok
                        sys.exit(10)
        for s in range1_l:
                if s!='':
                        print 'range1 is not null',range1_l,tok
                        sys.exit(10)
        for s in range2_l:
                if s!='':
                        print 'range1 is not null',range2_l,tok
                        sys.exit(10)

        str_struc=member_l[0]
        str_member=member_l[1:]
	for str in str_member:
		member_info=self.struclist.find__iss_os(str_struc,str)	

        struc=tok[6]
        var=["0","0","0","0","0"]

	self.add_struc(str_struc,struc)

        n=len(tok)

        i=0
        for ind in range(8,n,2):
                var[i]=tok[ind]
                i=i+1
        cmd=[]
	if funcname=='spackv':
		if len(str_member)!=1:
			print "error in spackv, number of member"
			print str_member
			sys.exit(10)
		i1=var[0]
		i2=var[1]
		par=var[2]
		par_info=self.struclist.find__iss_os(str_struc,str_member[0])
		typestr=par_info[0]
		member_type=par_info[0]
		if is_variable(par)==1:
		    try:
		    	var_type=self.var_dic[separate_array(par)]
		    except:
			var_type=[]
		    if len(var_type)>0:
                    	if sametype(var_type[1],member_type)==0:
	                	print thisfunc,"type mismatch",
                        	print "type of var[i],",var[i], var_type
                        	print "type of string",member_type
                        	sys.exit(10)
		    else:
			if self.printit==1:
				print thisprogram+strinfo+ "unknown type",par

		if par_info[3]=='string':
			print thisfunc,"type error"
			sys.exit(10)

		if typestr[0]=='i':
			callname="spackv_array_copy_i8_i"
		else:
			callname="spackv_array_copy_r8_r8"
		if par_info[2]=='1':
			cmd2="i_data_size=1; "
		else:
			cmd2="i_data_size=size("+struc+"(1)%"+str_member[0]+")"
		cmd.append(cmd2)
		if str_struc!="site" and str_struc!="spec":
			print "error in spackv, str_struc=",str_struc
			sys.exit(10)
		if job%10==0:
			# unpack
			pu='u'
		else:
			pu='p'
		if job/10==0:
			cmd2="do i_spackv="+i1+","+i2+" "
			cmd.append(cmd2)
			cmd2="call "+callname+"('"+pu+"',"
			cmd2=cmd2+struc+"(i_spackv)%"+str_member[0]+",i_data_size,1,"+par+")"
			cmd.append(cmd2)
			cmd.append("enddo")
		else:
			cmd2="do i_spackv="+i1+","+i2+" "
			cmd.append(cmd2)
			cmd2="call "+callname+"('"+pu+"',"
			cmd2=cmd2+struc+"(i_spackv)%"+str_member[0]+",i_data_size,i_spackv+1-"+i1+","+par+")"
			cmd.append(cmd2)
			cmd.append("enddo")

	elif funcname=="spacks":
		par=var[0]
		i1=var[1]
		i2=var[2]
		par_info=self.struclist.find__iss_os(str_struc,str_member[0])
		typestr=par_info[3]

                membertype=par_info[0]
		var_type=[]
		try:
                	var_type=self.var_dic[separate_array(par)]
		except:
			if self.printit==1:
				print thisprogram+strinfo+"check_variable", par ,"not found"

		if typestr!='string':	
			print 'error typestr,',typestr,str_struc,str_member[0],par,i1,i2
			sys.exit(10)
		callname='spacks_copy'
		if job==0:
			pu='u'
		else:
			pu='p'
		cmd2="do i_spacks="+i1+","+i2+" "
		cmd.append(cmd2)
                cmd2="call "+callname+"('"+pu+"',"
                cmd2=cmd2+struc+"(i_spacks)%"+str_member[0]+','+i1+','+i2+','+\
                    par+',i_spacks+1-'+i1+")"
		cmd.append(cmd2)
		cmd.append("enddo")

        return cmd

	
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

# manipulate [ dgets ( 'lat vol' , slat )  ]
# and so on.
#input, string: tok
    def do_functions__il_ol(self,tok):
	n =len(tok)
	for i in range(n):
		if type(tok[i])==ListType:
			#print "recursive do_func"
			atok=self.do_functions__il_ol(tok[i])
			tok[i]="("+atok[2]+")"

        tok=self.merge_variableinline__il_ol(tok)
	
	return [1,tok]
	



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
	error=0
	if tok[0]=='type':
		istart=1
		typesize=-1	
		typestr0=tok[0]
		l_out=self.typesize__is_ol(tok,istart)
                if l_out[0]:
                        istart=l_out[1]
                        typesize=l_out[2]
                else:
                        error=1
	elif tok[0]=='integer' or tok[0]=='real':
		istart=1
		typesize='0'
		typestr0=tok[0]
		l_out=self.typesize__is_ol(tok,istart)
		if l_out[0]:
			istart=l_out[1]
			typesize=l_out[2]
		else:
			error=1
	elif tok[0]=='character':
		istart=1
		typesize='0'
		typestr0=tok[0]
		l_out=self.typesize__is_ol(tok,istart)
                if l_out[0]:
                        istart=l_out[1]
                        typesize=l_out[2]
                else:
                        error=1
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
                        error=1
	elif tok[0]=='complex':
		istart=1
		typesize='4'
		typestr0='complex'
		l_out=self.typesize__is_ol(tok,istart)
                if l_out[0]:
                        istart=l_out[1]
                        typesize=l_out[2]
                else:
                        error=1
	elif tok[0]=='double' and tok[1]=='complex':
		istart=2
		typesize='8'
		typestr0='complex'
	elif tok[0]=='external' or tok[0]=='parameter' or tok[0]=='data':
		tok2=[]
		return tok2	
	else:
		error=1
	if error==1:
		print thisfunc,"error",tok
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
	

    def do_use__il_on(self,tok):
	thisfunc="do_use__il_on"
	nn=tok[1]
	if nn[0].isalpha()==1:
		self.uselist.append([self.execfuncname,nn])

    def make_reflist(self,line):
	self.reflist=[]
        thisfunc="make_reflist",line
        for name in line:
		if name=='(' or name==')' or name==',':
			continue
		else:
			self.reflist.append([name,0])

	print "reflist",self.reflist


    def line2tok(self,line):
        tok=line.split(",");
        tok=add_c__ils_ol(tok,',')
        tok = do_split0__ils_ol(tok,"'")
        tok = merge_string__il_ol(tok)
        tok = do_split1__ils_ol(tok,"(")
        tok = do_split1__ils_ol(tok,")")
        tok = do_split1__ils_ol(tok,"=")
        tok = do_split1__ils_ol(tok,"//")
        tok = do_split1__ils_ol(tok," ")
        tok = do_split1__ils_ol(tok,":")
        tok = do_split1__ils_ol(tok,"*")
        tok= del_spc__il_ol(tok)
	tok2=[]
	i=0
	while i< len(tok):
		if tok[i]=="=":
			if tok[i+1]=="=":
				tok2.append("==")
				i=i+1
			else:
				tok2.append(tok[i])
		elif tok[i]==":":
			if tok[i+1]==":":
				tok2.append("::")
				i=i+1
			else:
				tok2.append(tok[i])
		else:
			tok2.append(tok[i])
		i=i+1

        return tok2

    def indent_bracket0(self,tok,ilast):
	thisfunc="indent_bracket0"
        tok2=[]
	i=ilast
	tok2.append(tok[i])
	i=i+1
	tok2.append(tok[i])
	i=i+1
        while i < len(tok):
		name=tok[i]
                if name=="(": 
                        a=tok2.pop()
                        if a=="(":
                                tok2.append(a)
                                atmp=self.indent_bracket0(tok,i)
                        else:
                                atmp=self.indent_bracket0(tok,i-1)
                        tok2.append(atmp[0])
                        i=atmp[1]+1
		elif name==")":
			tok2.append(tok[i])
			return [tok2,i]
                else:
                        tok2.append(tok[i])
		i=i+1
	return [tok2,i]

    def indent_bracket(self,tok):
	tok2=[]
        i=0
	while i < len(tok):
		name=tok[i]
		if name=="(":
			a=tok2.pop()
			if a=="(":
				tok2.append(a)
				atmp=self.indent_bracket0(tok,i)
			else:
				atmp=self.indent_bracket0(tok,i-1)
			tok2.append(atmp[0])
			i=atmp[1]+1
		else:
			tok2.append(tok[i])
			i=i+1
	return tok2

    def do_nothing(self):
	i=1

    def flag_refval(self,tok,val):
	thisfunc="flag_refval"
	print thisfunc,tok
	var=''
	if isinstance(tok[-1],basestring)==1:
		var=tok[-1]
	elif isinstance(tok,list)==1:
		var=tok[0][0]
	else:
		print thisfunc,"error"
		sys.exit(10)
		
	print thisfunc,"var>",var
	for i in range(len(self.reflist)):
		name=self.reflist[i]
		if name[0]==var:
			self.reflist[i][1] = self.reflist[i][1] | val
			print thisfunc,"reflist, ",self.reflist[i]
			break
		

    def do_intentout(self,tok):
	thisfunc="do_intentout"
	idx=0
	haveeq=0
	for i in range(len(tok)):
		name=tok[i]
		if isinstance(name,basestring)==1 and name=="=":
			idx= i
			haveeq=1
			break
	if haveeq==1:
		print thisfunc,"left>",tok[:idx]
		self.flag_refval(tok[:idx],1)		

#input,string: line
    def do_line__is_ol(self,line):
	thisfunc='do_line__is_ol'
	tok=self.line2tok(line)

	if len(tok)==0:
		return [0,'']

	if len(tok)==1 and tok[0]=='end':
		if len(self.execfuncname)==0:
			print "error,",thisfunc,", found 'end', but function name is unknown"
			print tok
			sys.exit(10)
		l_tok2=[1,['end',self.execfuncname[1],self.execfuncname[0]]]
		return l_tok2

	if tok[0]=='subroutine':
		self.execfuncname=[tok[1],'subroutine']
		self.var_dic.clear()
                print tok
		self.make_reflist(tok[2:])
		return [0,[]]

	def_func=0
	for i in range(len(tok)):
		if tok[i]=='function':
			j=i-1
			if tok[j]=='integer' or tok[j]=='precision' or tok[j]=='logical':
				self.execfuncname=[tok[i+1],'function']
				self.make_reflist(tok[i+2:])
				def_func=1
	if def_func==1:
		self.var_dic.clear()
		return [0,[]]


	if tok[0]=='use':
		self.do_use__il_on(tok)
			
	if tok[0]=='integer' or tok[0]=='double' or tok[0]=='logical' or tok[0]=='real' or  tok[0]=='type' or tok[0]=='character' :
		varlist=self.do_typevar__il_on(tok)
		print varlist
		# make dictionary
		for var1 in varlist:
			self.var_dic[var1[0]]=[var1[1],var1[2]]

		return [0,[]]
	if tok[0]=='external':
		return [0,[]]

	tok=self.indent_bracket(tok)
	print  "tok2>",tok

        self.do_intentout(tok)

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


def make_list_recur__ili_ol(nspc,list):
	thisfunc='make_list_recur__ili_ol'
	tok=[]
	pad=''
	for i in range(nspc):
		pad=pad+' '
		padcont=''
	for i in range(nspc):
		if i==5:
			padcont=padcont+'.'
		else:
			padcont=padcont+' '
	cmd=make_line_recur__ili_ol(list)
        for name in cmd:
	    if len(name)>0:
		if name[0]==0:
                	#print pad+name[1]
			tok.append(pad+name[1])
		else:
			#print padcont+name[1]
			tok.append(padcont+name[1])
	#print
	tok.append('')
	return tok



def do_main(contline,printit0=1):
    thisfunc='do_main'
#parse constructor
    doparse=DoParse()

    doparse.printit=printit0
    doparse.read_newvar()

# default rule
    doparse.var_toadd.append(['*','bz','sbz'])
    doparse.var_toadd.append(['*','array','sarray'])
    doparse.var_toadd.append(['*','ctrl','sctrl'])
    doparse.var_toadd.append(['*','gw','sgw'])
    doparse.var_toadd.append(['*','ham','sham'])
    doparse.var_toadd.append(['*','lat','slat'] )
    doparse.var_toadd.append(['*','mix','smix'])
    doparse.var_toadd.append(['*','move','smove'])
    doparse.var_toadd.append(['*','optic','soptic'])
    doparse.var_toadd.append(['*','ordn','sordn'])
    doparse.var_toadd.append(['*','pot','spot'])
    doparse.var_toadd.append(['*','spec','sspec'])
    doparse.var_toadd.append(['*','site','ssite'])
    doparse.var_toadd.append(['*','str','sstr'])
    doparse.var_toadd.append(['*','tb','tsb'] )

#print
    printcont=1

    have_module=0
    #for packline in contline:
#	tok=line2tok__is_ol(packline[1])
#	if len(tok)>0 and tok[0]=='module':
#		have_module=1
#		break
#
#    if len(doparse.var_toadd)>0 and have_module==0 and doparse.printit==1:
#	print spc6+"module func_"+doparse.var_toadd[0][0]
#	print spc6+"contains"

    for packline in contline:
	if printcont==1:
		# parse mode
	    listchanged=doparse.do_line__is_ol(packline[1])
	    if listchanged[0]==1 and doparse.printit==1:
	        print_original__isl_on(thisprogram,packline)
	    elif listchanged[0]==2 and doparse.printit==1:
		print_original__isl_on('',packline)
	    if listchanged[0]>0:
		nspc=count_spc__il_oi(packline[0])
		cmd=make_list_recur__ili_ol(nspc,listchanged[1])
		if doparse.printit==1:
		    for nn in cmd:
			print nn
	    else:
		if doparse.printit==1:
	        	print_original__isl_on('',packline)

	elif printcont==2 and doparse.printit==1:
		# print original line mode
	    for list in packline[0]:
		print "3>", list
	elif printcont==3 and doparse.printit==1:
		print packline

#    if len(doparse.var_toadd)>0 and have_module==0 and doparse.printit==1:
#	print spc6+"end module func_"+doparse.var_toadd[0][0]
			
#print doparse.var_dic

    doparse.write_newvar()


#---main---
#constructure
fline=mFLine.FLine()
contline=fline.get__in_ol()

do_main(contline,0)
print "again"
do_main(contline,1)

