#!/usr/bin/env python
#os.path.abspath(sys.argv[0]) 
from math import *
import re,sys,decimal,os,string
#rootpath=os.path.dirname(os.path.abspath(sys.argv[0])) 
#print rootpath
#sys.path.append(rootpath+'/convert') #add path to /convert
#import lmf2ctrl
import convctrl

argvs = sys.argv
argc  = len(argvs)
print 'readin args= ',argvs,argc
if ('--help' in argvs): 
	print ' == Convert POSCAR file(input file of lmf) to ctrl =='
	print '    usage: vasp2ctrl POSCAR_foobar.vasp   [option] '
	print '          Then we have ctrls.foobar.vasp2ctrl        '
	print '     option: --alat=10.66 (ALAT, you can use eqs such as 10.66*1.88)'
	sys.exit(-1)
for ix in argvs:
	if 'POSCAR_' in ix: 
		print ' Read ', ix, '. It is converted to ctrl file now.'
		ext1=re.sub('POSCAR_','',ix)
		ext2=re.sub('.vasp','',ext1)
		break
	# elif '.cif' in ix: 
	# 	ext1 = re.sub('.cif','',ix)
	# 	ext2 = re.sub('.vasp','',ext1)
	# 	break
for ix in argvs:
	if '--alat=' in ix: 
		alatin=string.atof(eval(re.sub('--alat=','',ix) ))
	else:
		alatin=None
		ratioa=1.00000000000
titleinput = 'ctrls.'+ext2

vaspread = open(argvs[1]).read().split('\n') 
plat1 = vaspread[2].split()
plat2 = vaspread[3].split()
plat3 = vaspread[4].split()
print ' '+vaspread[7]

alat_val = convctrl.vasp2ctrl_alat(vaspread) #unit
if(alatin): ratioa = alatin/alat_val #conversion by given --alat=alatin
all_atom,NBAS_val = convctrl.vasp2ctrl_atomcount(vaspread)
#atom_list   = convctrl.vasp2ctrl_atom(vaspread,alat_val,NBAS_val)
atom_list = convctrl.vasp2ctrl_atom(vaspread,alat_val,NBAS_val,plat1,plat2,plat3)
convctrl.vasp2ctrl_write(vaspread,alat_val,NBAS_val,atom_list,titleinput,all_atom,ratioa)

print ' OK! we have ',titleinput+'.vasp2ctrl'


