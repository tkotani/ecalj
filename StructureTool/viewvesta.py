#!/usr/bin/env python3
import subprocess,sys
import sys
import os
import re

out = subprocess.getoutput("which VESTA")
if(len(out)==0): out = subprocess.getoutput("which vesta")
VESTA=str(out)
if(len(VESTA)==0):
        print('No vesta found!')
        sys.exit()
#print(VESTA)
rootpath=os.path.dirname(os.path.abspath(sys.argv[0])) 

argvs = sys.argv
argc = len(argvs)

print( 'args= ',argvs,argc)

if (argc != 2 or '--help' in argvs):
	print( ' == Image Display VASP(input file of lmf) to ctrl ==')
	print( '    usage: viewvesta POSCAR_foo.vasp    ')
	sys.exit(-1)

fname=argvs[1]
ix=0
if argvs[1][0:5]=='ctrls':
	fname = re.sub('ctrls.','POSCAR_',argvs[1])+'.vasp'
	ix=1
elif argvs[1][0:4]=='ctrl':
	fname = re.sub('ctrl.','POSCAR_',argvs[1])+'.vasp'
	ix=1
if ix==1:
	if os.path.isfile(rootpath+'/ctrl2vasp.py'):
		os.system(rootpath+'/ctrl2vasp.py '+argvs[1])
	else:
		os.system(rootpath+'/ctrl2vasp '+argvs[1])
print('cwd=',os.getcwd())
print(fname)
vestaopen= VESTA + ' '+os.getcwd()+'/'+fname
print( vestaopen)
os.system(vestaopen) 

