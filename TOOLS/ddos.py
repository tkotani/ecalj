#!/usr/bin/python
import sys
import string
import os
print
print "========= plotg  ========="
try:
	print "Readin two files =     ", sys.argv[1], sys.argv[2]
except:	
	print "usage: dqpu QPU1 QPU2"
	sys.exit()

fr = open(sys.argv[1],'rt')
oxx=fr.read()
oxx=string.split(oxx,'\n')
fr.close()

fr = open(sys.argv[2],'rt')
oyy=fr.read()
oyy=string.split(oyy,'\n')
fr.close()

#for ix in range(len(oxx)):
dmax=0e0
for ix in range(0,len(oxx)):
	iline=oxx[ix]
	ilin2=oyy[ix]

	if ix>=0:
		l1= string.split(iline)
		l2= string.split(ilin2)
		#print l1, l2,len(l1)
		for iw in range(0,len(l1)):
			try:
				w1=string.atof(l1[iw])
				w2=string.atof(l2[iw])
				aaa = '%6.6f ' % w1
				aaa = aaa+ '%6.6f ' % w2
				ddd= w1-w2
				aaa = aaa+ '%6.6f ' % ddd
				if(abs(ddd)>dmax):
					dmax=abs(ddd)
					print l1,l2, aaa
			except:
				print 'xxx'
print ' dmax=',dmax
