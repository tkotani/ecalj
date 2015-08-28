#!/usr/bin/python
import os, sys,re
f1 = open("PROCAR.UP", 'r')
lines = f1.readlines()
nline= len(lines)
#print nline
ionline=0
for iline in lines:
	if 'k-point' in iline: 
		x = (iline.split('x =')[1]).split('\n')[0]
		iband=0
		print
		print '# ' + iline
		ionline=0
	if 'band' in iline:
		if iband>0:
			www= '%12.5f' % weight
			print x +' '+y + www   # NOTE: this is for the previous band, where
		                               # we already 'weight' calculated.
		ionline=0
		weight=0.0
		iband=iband+1
		y= iline.split('energy')[1].split('#')[0]

	if ionline ==1:
		ilines=re.split('\s+',iline)
		if 'tot' in iline : continue #continue goto next iteration. Different from continue of 
		if len(ilines)<3  : continue
		#print ilines,len(ilines)
		if(int(ilines[1])==2) : #oxygen site
			weight = weight + float(ilines[3])+float(ilines[4])+float(ilines[5])  #py + pz +px
	if 'ion' in iline:
		ionline=1

sys.exit()
