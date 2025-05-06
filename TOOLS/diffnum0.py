#!/usr/bin/env python3
import sys
#import string
import os
import re

def linewithkey(comparekeys,lines):
	aaa=[]
	for iline in lines:
		ix=0
		for key in comparekeys:
			if( re.match(key,iline) ): ix=1
		if(ix==1): aaa=aaa+ [iline]
	#print aaa
	return aaa

def comparenum(tol,file1,file2,comparekeys,printsw):
	pr = (printsw==1)
	oxx= file1.split('\n') 
	oyy= file2.split('\n')
	errmax=0.0
	ix1=0
	ix2=0
	if(len(comparekeys)>0):
		oxx=linewithkey(comparekeys+['written','comparing'],oxx)
		oyy=linewithkey(comparekeys+['written','comparing'],oyy)
#	for ix in range(len(oxx)):
#		print ix,oxx[ix]
#	print oxx
#	print len(oyy),max(len(oxx),len(oyy))
#	sys.exit()

	ierrl=0
	for ix in range(max(len(oxx),len(oyy))):
		try: #check end of oxx or not
			izz=ix1
			while izz < len(oxx):
				izz = izz+1
				iline = oxx[izz]
				iii1 = iline.split()
				if(not iii1 ==[]): break
		except:
			break

		while ix1 < len(oxx):
			ix1 = ix1+1
			iline= re.sub(r'D\+','e+',oxx[ix1])
			iline= re.sub(r'D\-','e-',iline)
			iii1= iline.split()
			if(not iii1 ==[]): break
		#print	
		#print '1=',ix1,iii1

		while ix2 < len(oyy):
			ix2 = ix2+1
			ilin2= re.sub(r'D\+','e+',oyy[ix2])
			ilin2= re.sub(r'D\-','e-',ilin2)
			iii2= ilin2.split()
			if(not iii2 ==[]): break
		#print '2=',ix2,iii2

		#skip lines including these keywords.
		#if 'written' in set(iii1):
		#	continue
		#if 'comparing' in set(iii1):
		#	continue


		for i in range(len(iii1)):
			ixx=0
			try:
				w1=float(iii1[i])
				ixx=1
			except:
				continue

			if(ixx==1):
				w2=float(iii2[i])
				#print 
				#print ' w1=',w1
				#print ' w2=',w2
				#print ' abs 1-w2= ',abs(w1-w2)
				if (  abs(w1-w2) > errmax ): errmax = abs(w1-w2)

				if (  abs(w1-w2) > tol ):
					errmax=abs(w1-w2)
					print()
					print( ix,iline)
					print( ix,ilin2)
					ierrl=ierrl+1
					if(ierrl>10):
						print ('Error lines (less than ten shown) with > tol=',tol)
						return errmax

	return errmax
