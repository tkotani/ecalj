import sys
from types import *

def list2string(tok):
	if isinstance(tok,list)==0:
		return ''
	if len(tok)==0:
		return ''
	s='['
	for i in range(len(tok)):
		nn=tok[i]
		if i!=0:
			s=s+' '
		if isinstance(nn,basestring)==0:
			s=s+str(nn)
		else:
			if len(nn)==0:
				s=s+"''"
			else:
				s=s+nn
	s=s+']'
	return s

# 1 2 3 4 -> 1 c 2 c 3 c 4
#input list: tok
#input string: c
def add_c__ils_ol(tok,c):
        n=len(tok)
        i=0
        tok2=[]
        while (i<n):
                tok2.append( tok[i] )
                if (i != n-1 ):
                        tok2.append(c)
                i=i+1
        return tok2


# w, (, sss, ), a  -> 'w(sss)', a
#input list: tok
def merge_w__il_ol(tok):
        n=len(tok)
        i=0
        tok2=[]
        while i<n:
                s=tok[i]
                if  s=="w":
                        wstr=s
                        j=i+1
                        indent=0
                        while j<n:
                                t=tok[j]
                                wstr=wstr+t
                                if t=="(":
                                        indent=indent+1
                                if t==")":
                                        indent=indent-1
                                if (indent==0):
                                        i=j
                                        s=wstr
                                        break
                                j=j+1
                tok2.append(s)
                i=i+1
        return tok2


# 'var' '*' 'i' -> 'var*i'
def merge_ast__il_ol(tok):
        #print "merge_ast1>",tok
        tok2=[]
        i=0
        nmax=len(tok)
        while i<nmax:
                nn=tok[i]
                if  nn=='*':
                        str=tok2.pop()
                        i=i+1
                        str=str+nn+tok[i]
                        tok2.append(str)
                elif nn[0]=='*':
                        str=tok2.pop()
                        str=str+nn
                        tok2.append(str)
                else:
                        tok2.append(nn)
                i=i+1
        #print "merge_ast2>",tok2
        return tok2



# ', sss, ' -> 'sss'
#input, string: tok
def merge_string__il_ol(tok):
        n=len(tok)
        i=0
        tok2=[]
        while i<n:
                s=tok[i]
                if s=="'":
                        str=s
                        j=i+1
                        while j<n:
                                str=str+tok[j]
                                if tok[j]=="'":
                                        i=j
                                        break;
                                j=j+1
                        tok2.append(str)
                        i=i+1
                else:
                        tok2.append(s)
                        i=i+1
        return tok2

# 123c345c5 -> 123 345 5
#input, list: tok
#input,string: c
def do_split0__ils_ol(tok,c):
        tok2=[]
        for s in tok:
             s=s.split(c)
             s=add_c__ils_ol(s,c)
             tok2=tok2+s
        return tok2


# 123c345c5 -> 123 345 5
# but don't separate 'sss  bbb'
#input, list: tok
#input, string: c
def do_split1__ils_ol(tok,c):
        tok2=[]
        for s in tok:
                if s[0:1] == "'":
                        tok2.append(s)
                else:
                        s=s.split(c)
                        s=add_c__ils_ol(s,c)
                        tok2=tok2+s
        return tok2


# a b c ' ' d '' e -> a b c d e
#input, list: tok
def del_spc__il_ol(tok):
        tok2=[]
        for s in tok:
                if len(s)==0:
                        continue
                if s==" ":
                        continue
                tok2.append(s)
        return tok2

#   a, b, c, //, d, e
# -> a, b, c//d, e
#input, list: tok
def merge_stringslash__il_ol(tok):
        tok2=[]
        i=0
        n=len(tok)
        while i<n:
                s=tok[i]
                if i+1<n and tok[i+1]=="//":
                        s=tok[i]+tok[i+1]+tok[i+2]
                        i=i+2
                tok2.append(s)
                i=i+1
        return tok2


# spacks(0,'spec name',sspec,alabl, ips ( ib1 ) , ips ( ib1 ) )
# ->
# spacks(0,'spec name',sspec,alabl, ips(ib1) , ips(ib1))
#input,string: tok
def merge_variable__il_pi_ol(tok,indentlevel=2):
	thisfunc='merge_variable__il_pi_ol'
        tok2=[]
        n=len(tok)
        indent=0
        i=0
        while i<n:
                s=tok[i]
                tok2.append(s)
                if s=="(":
                        indent=indent+1
                elif s==")":
                        indent=indent-1
                if s=="(":
                    if indent==indentlevel:
                        str2=tok2.pop()
			str1=''
			try:
                        	str1=tok2.pop()
			except:
				str1=''
                        str=str1+str2
                        j=i+1
                        while j<n:
                                s=tok[j]
                                str=str+s
                                if s=="(":
                                        indent=indent+1
                                if s==")":
                                        indent=indent-1
                                if s==")" and indent==indentlevel-1:
                                        tok2.append(str)
                                        break
                                j=j+1
                        i=j


                i=i+1
        return tok2


# i = igets ( 'site pos' , ssite ) -> i = [ igets ( 'site pos' , ssite ) ]
#input, list: tok
#input,string:  name
def separate_function__ils_ol(tok,name):
        tok2=[]
        i=0
        n=len(tok)
        while i<n:
                s=tok[i]
                if s==name:
                        tok3=[]
                        tok3.append(s)
                        indent=0
                        j=i+1
                        while j<n:
                                t=tok[j]
                                if t=="(":
                                        indent=indent+1
                                if t==")":
                                        indent=indent-1
                                tok3.append(t)
                                if indent==0:
                                        s=tok3
                                        i=j
                                        break
                                j=j+1

                tok2.append(s)
                i=i+1
        return tok2


#input, string: line
def del_ln__is_os(line):
        str=''
        for c in line:
                if c!='\n':
                        str=str+c
        return str


def separate_array__il_ol(nn):
        list1=nn.split('(')
        list1 = do_split1__ils_ol(list1,"*")

        return [list1[0],nn]


def separate_array(s):
        if isinstance(s,basestring)==0:
                print "separate_array error type"
                sys.exit(10)
        if len(s)==0:
                return s
        a=''
        j=len(s)
        for i in range(len(s)):
                if s[i]=='(':
                        j=i
                        break
        if j>=0:
                return s[:j]
        return ''


def fix_list__il_ol(lis):
        thisfunc="fix_list__il_ol"
        tok=[]
        for n1 in lis:
                if isinstance(n1,list):
                        for n2 in n1:
                                tok.append(n2)
                else:
                        tok.append(n1)
        return tok


def delete_lastcall__il_ol(tok):
        name=tok.pop()
        if name!='call':
                print  thisfunc, "error last is not call"
                print tok
                sys.exit(10)
	#search if ... then
	add_then=0
	if len(tok)>0:
	    if tok[0]=='if':
		# search then
		have_then=0
		for name in tok:
			if name=='then':
				have_then=1
				break
		if have_then==0:
			tok.append('then')
			add_then=1
	
        return [add_then,tok]


def is_variable(s):
# return value
# 0: not variable
# 1: variable except w(...)
# 2: w variable
        if isinstance(s,basestring)==0:
                print "is_variable error type"
                sys.exit(10)
        if s[0]=='w' and len(s)>1:
                if s[1]=='(':
                        return 2
        c=s[0]
        if c.isalpha()==1:
                return 1
        return 0


def count_spc__is_oi(c):
# count # of spaces at the beginning of the line
# input c is string
        if isinstance(c,basestring)==0:
                print "count_spc__is_oi, error type"
                sys.exit(0)
        n=len(c)
        for i in range(len(c)):
                if c[i]!=' ':
                        n=i+1
                        break
        return n

def count_spc__il_oi(list,level=0):
# count # of spaces at the beginning of the line
# input list is a list
        if len(list)==0:
                print "error listsize=0"
                sys.exit(10)
        if type(list[0])==ListType:
                return  count_spc__il_oi(list[0],level+1)
        return  count_spc__is_oi(list[0])


def sametype(a,b):
        if isinstance(a,basestring)==0:
                print "is_variable error type, a"
                sys.exit(10)
        if isinstance(b,basestring)==0:
                print "is_variable error type, b"
                sys.exit(10)

        if a==b:
                return 1
        # integer(4) integer(8) integer
        if a[0]==b[0] and a[0]=='i':
                return 1
        return 0


def line2tok__is_ol(line):
        tok=line.split(",");
        tok=add_c__ils_ol(tok,',')
        tok = do_split0__ils_ol(tok,"'")
        tok = merge_string__il_ol(tok)
        tok = do_split1__ils_ol(tok,"(")
        tok = do_split1__ils_ol(tok,")")
        tok = do_split1__ils_ol(tok,"=")
        tok = do_split1__ils_ol(tok,"//")
        tok = do_split1__ils_ol(tok," ")
        tok = do_split1__ils_ol(tok,"::")
        tok = do_split1__ils_ol(tok,"*")
        tok = do_split1__ils_ol(tok,";")
        tok= del_spc__il_ol(tok)
        return tok

