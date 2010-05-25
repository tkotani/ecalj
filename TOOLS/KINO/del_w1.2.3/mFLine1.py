import sys
from types import *

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
		if len(name)>1 and name[0]!=' ':
		    outline=''
#234567
#    .do  i=
		elif len(name)>6:
		    outline=name[6:]
		else:
		    outline=''
		longline= outline
	else:
		longline=''
		for ll in line:
			longline=longline+ self.flush__il_pi_os(ll,level+1)
        return longline



  def get__in_ol(self):
#return = [ [[line1,linenumber1],[line2,linenumber2],...] , continued_line]
	linelist=[]
	iline=0
    	for aline in sys.stdin:
		iline=iline+1
        	aline= self.del_ln__is_os(aline)
#1234c6
		newtok=[aline,iline]
		linelist.append(newtok)
		if len(aline)>0 and aline[0]!=' ':
			continue
		elif len(aline)>6 and aline[5]!=' ':
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
					if len(str)>0 and str[0]!=' ': # comment
						continue
					elif len(str)>6 and str[5]!=' ': # continue
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

	
