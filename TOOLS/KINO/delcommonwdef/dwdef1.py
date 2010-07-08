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

thisprogram='Cdwdef1'
strinfo='...info...'+spc10

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

class DoParse:
    printit=0
    varlist=[]
    wdeffunc=[]
    wdic={}
    def __init__(self):
        printit=0
	varlist=[]
        wdeffunc=[]
    def find_subname(self,tok):
        if tok[0]=='subroutine':
                self.execfuncname=[tok[1],'subroutine']
                return [0,'']
	if tok[0]=='program':
		self.execfuncname=[tok[1],'program']
		return [0,'']
        if 'function' in tok:
            def_func=0
            for i in range(len(tok)):
                if tok[i]=='function':
                        j=i-1
                        if tok[j]=='integer' or tok[j]=='precision' or tok[j]=='logical':
                                self.execfuncname=[tok[i+1],'function']
                                def_func=1
                                break
            if def_func==1:
                if i+1 >= len(tok):
                        print thisfunc,"ERROR try to access tok[i+1], i+1=",i+1, "len=",len(tok)
                        print thisfunc, linenumber,tok
                        print >>sys.stderr, thisfunc,"ERROR try to access tok[i+1], i+1=",i+1, "len=",len(tok)
                        print  >>sys.stderr,thisfunc, linenumber,tok
                        # don't exit to know linenumber 
                        #sys.exit(10) 
                return [0,'']
	return [1,tok]

    def do_line_findwdef(self,line,linenumber,fnlabel):
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

	iret,newtok=self.find_subname(tok)

	if "w" in tok:
		if tok[0]=="common":
			self.wdeffunc.append(self.execfuncname[0])
			self.wdic[self.execfuncname[0]]=0
		elif tok[0]=="integer" or tok[0]=="real":
                        idummy=0
		else:
		    if self.execfuncname[0] in self.wdic:
			self.wdic[self.execfuncname[0]]=self.wdic[self.execfuncname[0]]+1
	
	return [0,'']	

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

        iret,newtok=self.find_subname(tok)
	if iret==0:
		return [0,'']

	if len(self.execfuncname)==0:
		return [0,'']
	if not self.execfuncname[0] in self.wdeffunc:
		return [0,'']

        if "w" in tok:
                if tok[0]=="common":
			if self.wdic[self.execfuncname[0]]>0:
                        	return [1,['integer','::' 'iwdummy']]
			else:
                        	return [2,'']

		elif tok[0]=="integer" or tok[0]=="real":
			if tok[1]=="w":
				return [2,'']
			else:
				print "#error unknown wdef",tok
				sys.exit(10)
		else:
			for i in range(len(tok)):
				if tok[i]=="w":
				    if tok[i+1]=="(":
					print "#error have w(...)", tok
					sys.exit(10)
				    else:
					tok[i]="iwdummy"
			return [1,tok]
				
        return [0,'']




def do_main(contline):
    thisfunc='do_main'
    doparse=DoParse()

    for packline in contline:
            fnlabel=get_fnlabel(packline[0])
            if len(packline[1])==0:
                fnlabel=''
            linenumber=get_linenumber(packline[0])
	    ichanged,tok=doparse.do_line_findwdef(packline[1],linenumber,fnlabel)

    if len(doparse.wdeffunc)>0:
    	for nn in  doparse.wdeffunc:
		print thisprogram,"have_common_w_in",nn,", #_of_w_access=",doparse.wdic[nn]

    if 1:
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
	    elif ichanged==2:
		print_original__isl_on(thisprogram,packline)



#---main---
thisfunc="main"
argv=sys.argv
argc=len(argv)
#print argc,argv
#print thisprogram,"cmdvar=",cmdvar
#constructure
fline=mFLine.FLine()
contline=fline.get__in_ol()


do_main(contline)

