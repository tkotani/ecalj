#!/usr/bin/env python
VESTA='~/VESTA-x86_64/VESTA'
import sys
import os
import re
rootpath=os.path.dirname(os.path.abspath(sys.argv[0])) 

argvs = sys.argv
argc = len(argvs)

print 'args= ',argvs,argc

if (argc != 2 or '--help' in argvs):
	print ' == Image Display VASP(input file of lmf) to ctrl =='
	print '    usage: viewvesta POSCAR_foo.vasp    '
	sys.exit(-1)

fname=argvs[1]
if argvs[1][0:5]=='ctrls':
	fname = re.sub('ctrls.','POSCAR_',argvs[1])+'.vasp'
elif argvs[1][0:4]=='ctrl':
	fname = re.sub('ctrl.','POSCAR_',argvs[1])+'.vasp'
os.system(rootpath+'/ctrl2vasp.py '+argvs[1])

print fname
vestaopen= VESTA + ' '+os.getcwd()+'/'+fname
print vestaopen
os.system(vestaopen) 

