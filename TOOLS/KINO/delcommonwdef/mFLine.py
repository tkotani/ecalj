import sys
from types import *
from mToken import *
import re

class FLine:

  linework=[]
	
  def del_ln__is_os(self,line):
#          isinstance(line,string)
        str=''
        for c in line:
                if c!='\n':
                        str=str+c
        return str

  def correctlist0__il_on(self,l_list):
	if isinstance(l_list,list)==0:
		print 'correctlist0__il_on>,error type'
		sys.exit(10)
	for nn in l_list:
		if len(nn)==0:
			continue
		elif isinstance(nn[0],list):
			self.linework.append(self.correctlist0__il_on(nn))	
		elif isinstance(nn[0],basestring):
			self.linework.append(nn)
		else:
			print "error,nn=", nn
			sys.exit(10)

  def correctlist__il_ol(self,l_list):
#   l_list=[[[A,B],[C,D],None],[E,F]]
#  -> return = [[A,B],[C,D],[E,F] 
	if isinstance(l_list,list)==0:
		print "correctlist__il_ol> error type"
		sys.exit(10)
	self.linework=[]
	self.correctlist0__il_on(l_list)
	tok=[]
	for nn in self.linework:
		if isinstance(nn,list)>0:
			tok.append(nn)
	return tok



  def flush__il_pi_os(self,line,level=0):
	thisfunc='flush__il_pi_os'
#    line=['    program foo",1]
# or 
#    line=[['    subroutine foo(i,',1],
#          ['   . j)',2]] 
# return= [ [['    subroutine foo(i,',1],['   . j)',2]],'subroutine foo(i,i)']

	if isinstance(line,list)==0:
		print "arg error, flush__il_os"
		sys.exit(10)
	if isinstance(line[0],basestring)==1:
		name=line[0]
		if self.is_comment(name):
		    outline=''
#234567
#    .do  i=
		elif len(name)>6:
			#separate and merge again to delete comments
		    ll=line2tok__is_ol(name[6:])
		    str=''
		    for nn in ll:
			# delete comment
			if nn[0]=='!':
				break
			str=str+nn+' '
		    outline=str
		else:
		    outline=''
		longline= outline
	else:
		longline=''
		for ll in line:
			longline=longline+ self.flush__il_pi_os(ll,level+1)
        return longline

  def is_comment(self,aline):
      thisfunc="is_comment"
      if len(aline)>0 and aline[0].isdigit():
	  ret=0
	  return ret
      if len(aline)>0 and aline[0]!=' ':
          ret= 1
	  return ret
      return 0

  def is_cont(self,aline):
      	if len(aline)>5 and aline[5]!=' ':
		return 1
	return 0

  def get__in_ol(self):
	thisfunc='get__in_ol'
#return = [ [[line1,linenumber1],[line2,linenumber2],...] , continued_line]
	linelist=[]
	iline=0
    	for aline in sys.stdin:
		iline=iline+1
        	aline=self.del_ln__is_os(aline)
		#aline=aline.replace('\t','        ')
		aline=re.sub("^\t",'        ',aline)
#1234c6
		newtok=[aline,iline]
		linelist.append(newtok)
		if self.is_comment(aline)==1:
			continue
		elif self.is_cont(aline)==1:
			# pop and repack
			tmplist=[]
			while 1==1:
				atmp=[]
				try:
					atmp=linelist.pop()
				except:
					print "atmp.pop error"
					break
				tmplist.append(atmp)
					
				if len(atmp)==0:
					break
				if isinstance(atmp[0],list): # list
					tmplist.reverse()
					tmplist=self.correctlist__il_ol(tmplist)
					linelist.append(tmplist)
					break
				if isinstance(atmp[0],basestring):
					str=atmp[0]
					if self.is_comment(str)==1:
						continue
					elif self.is_cont(str)==1:
						continue
					else:
						tmplist.reverse()
						tmplist=self.correctlist__il_ol(tmplist)
						linelist.append(tmplist)
						break
		
	linelist2=[]
	for nn in linelist:
		ll=self.flush__il_pi_os(nn)
		linelist2.append([nn,ll])
	return linelist2

#---main---
if __name__ == "__main__":
#class
    line=FLine()

#list
    contline=line.get__in_ol()

    job=2
    for nn in contline:
    	if job==1:
	    if isinstance(nn[0],list):
		for mm in nn:
			print "--->",mm
	    else:
		print nn
    	elif job==2:
	    print nn

	
