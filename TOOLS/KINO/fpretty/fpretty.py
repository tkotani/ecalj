#
#
#  usage:
#   python thisfile < file.F > x
#
#
import os
import re
import sys
from types import *
import string

import mFLine
import mStrucList

from mToken import *

thisprogram="findent"

#spc6='123456'
spc6='      '

# change indentstr if you want to change indentlevel
indentstr='  '
indentstrcont='  '

class t_indentstack:
	word=''
	fnlabel=''
	linenumber=''
	def __init__(self,word,fnlabel,linenumber):
		self.word=word
		self.fnlabel=fnlabel
		self.linenumber=linenumber
	def __str__(self):
		p=' '
		return word+p+fnlabel+p+linenumber
	

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


class Parse:
    indentstack=[]
    fline=[]
    def __init__(self):
	self.fline=mFLine.FLine()
	contline=self.fline.get__in_ol()
	self.indentstack=[]
        self.do_main(contline)

    def indentstack_release(self,fnlabel,searchw):
        thisfunc="indentstack_release"
        if len(fnlabel)>0:      
                nfirst=len(self.indentstack)
                idel=-1
                for i in range(len(self.indentstack)):
                        #word,label=self.indentstack[i]
			word=self.indentstack[i].word
                        label=self.indentstack[i].fnlabel
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
                nfirst=len(self.indentstack)
                idel=-1
                for i in range(len(self.indentstack)):
                        #word,label=self.indentstack[i]
			word=self.indentstack[i].word
			label=self.indentstack[i].fnlabel
                        if word==searchw:
                                # find the last one
                                idel=i
                if idel>=0:
                        tok2=[]
                        for i in range(idel):
                                tok2.append(self.indentstack[i])
                        self.indentstack=tok2
                nlast=len(self.indentstack)

        return nfirst-nlast



    def checkindent(self,level,tok):
	thisfunc="checkindent"
	if len(self.indentstack)!=level:
		print  thisfunc, "#error ",thisprogram,"indent inconsistent",len(self.indentstack),level,tok
		print >>sys.stderr, thisfunc, "#error ",thisprogram,"indent inconsistent",len(self.indentstack),level,tok
		sys.exit(10)


    def do_line(self,fnlabel,linenumber,tok):
		if len(tok)==0:
			return
		if (len(tok)==1 and tok[0]=='end'):
			# end of subroutine/function
                        self.checkindent(0,tok)
			self.indentstack=[]
		if len(tok)>=2 and (tok[0]=='end' and (tok[1]=='subroutine' or tok[1]=='function')):
			# anther end of subrotine/function
                        self.checkindent(0,tok)
			self.indentstack=[]
		if tok[0]=='subroutine':
			# start subroutine 
                        self.checkindent(0,tok)
			self.indentstack=[]
		if 'function' in tok:
                        self.checkindent(0,tok)
			def_func=0
			for i in range(len(tok)):
			    if tok[i]=='function':
                        	j=i-1
                        	if tok[j]=='integer' or tok[j]=='precision' or tok[j]=='logical':
					def_func=1
					break
			if def_func==1:
				#start function
				self.indentstack=[]
		if tok[0]=='use'and tok[1]!='=':
			# use sentense
                        self.checkindent(0,tok)
			self.indentstack=[]
		if tok[0]=='integer' or tok[0]=='double' or tok[0]=='logical' or tok[0]=='real' or  tok[0]=='type ' or tok[0]=='character' :
			# type
                        self.checkindent(0,tok)
			self.indentstack=[]
		if tok[0]=='external' or tok[0]=='data' or tok[0]=='parameter':
			# external/data/parameter
                        self.checkindent(0,tok)
			self.indentstack=[]
		if len(tok)>=2 and tok[0]=='end' and (tok[1]=='module' or tok[1]=='subroutine' or tok[1]=="function"):
			# end/subroutine/module
                        self.checkindent(0,tok)
			self.indentstack=[]

		#--------------
		if tok[0]=="if":
			if "then" in tok:
				self.indentstack.append(t_indentstack('if','',linenumber))
		if tok[0]=="elseif" or (len(tok)>=2 and tok[0]=="else" and tok[1]=="if"):
                	if "then" in tok:
				self.indentstack.append(t_indentstack('elseif','',linenumber))
		if tok[0]=="else":
			self.indentstack.append(t_indentstack('else','',linenumber))
		if 'endif' in tok or (len(tok)>=2 and tok[0]=="end" and tok[1]=="if"):
			self.indentstack_release(fnlabel,'if')
		if tok[0]=="do":
			if len(tok)>1:
				a=re.match("[0-9]+",tok[1])
                        	if a != None:
					self.indentstack.append(t_indentstack("do",tok[1],linenumber))
				else:
					self.indentstack.append(t_indentstack("do",'',linenumber))
			else:
				self.indentstack.append(t_indentstack("do","",linenumber))
		if 'enddo' in tok or 'end' in tok:
			endlist=[]
                	for i in range(len(tok)):
                            if tok[i]=='enddo' or (i+1<len(tok) and  tok[i]=='end' and tok[i+1]=='do'):
                                endlist.append(i)
                	if len(endlist)>0:
                    	    for i in range(len(endlist)):
                        	n=self.indentstack_release('','do')
		if len(fnlabel)>0:
			n=self.indentstack_release(fnlabel,'do')
		if tok[0]=="select":
                	self.indentstack.append(t_indentstack("select",fnlabel,linenumber))
		if tok[0]=="endselect" or (len(tok)>=2 and tok[0]=="end" and tok[1]=="select"):
                	self.indentstack_release(fnlabel,"select")
		if tok[0]=='case':
			self.indentstack.append(t_indentstack("case",fnlabel,linenumber))

    def flush_indent(self,tok,level=0):
	thisfunc="flush_indent"
	ret=[]
	if isinstance(tok[0],basestring)==1:
	    return self.flush_indent([tok],level)
	if isinstance(tok[0],list)==0:
		print thisfunc,"ERROR",tok[0]
		sys.exit(10)
	#for line in tok:
	for id in range(len(tok)):
	    line=tok[id]
            if isinstance(line[0],basestring)==1:
                outline=line[0]
                name=line[0]
                if self.fline.is_comment(name):
                    outline=line[0]
#234567
#    .do  i=
                elif len(name)>6:
                        #separate and merge again to delete comments
		    str1=name[:6]
		    str2=name[6:]
		    str2=str2.lstrip()
		    outline=str1
		    for i in range(level):
			outline=outline+indentstr
		    if id!=0:
			outline=outline+indentstrcont
		    outline=outline+str2
		    	
                ret.append(outline)
            else:
		print thisfunc,"ERROR",line[0]

        return  ret

    def indentstack_level(self):
	level=0
	for nn in self.indentstack:
		if nn.word in ['if','select','do']:
			level=level+1
	return level
	
		
    def do_main(self,contline):
	for packline in contline:
		fnlabel=get_fnlabel(packline[0])
		if len(packline[1])==0:
                	fnlabel=''
		linenumber=get_linenumber(packline[0])
		tok=line2tok__is_ol(packline[1])
		self.do_line(fnlabel,linenumber,tok)
		indentlevel=self.indentstack_level()
		level=len(self.indentstack)
		flag1=0
		flag2=0
		if level>0 :
			flag1= linenumber== self.indentstack[level-1].linenumber
			flag2= self.indentstack[level-1].word in ['if','select','do','else','elseif','case']
			if flag1 and flag2:
				indentlevel=indentlevel-1
		ret=self.flush_indent(packline[0],indentlevel)
		if flag1:
			sflag1='T'
		else:
			sflag1='F'
		if flag2:
			sflag2='T'
		else:
			sflag2='F'

		for nn in ret:
			#print "L=",fnlabel,"I=",level,indentlevel,sflag1,sflag2,"#=",linenumber,nn
			print nn

#-------------------------------------

parse=Parse()

