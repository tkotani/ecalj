#!/usr/bin/env python
#os.path.abspath(sys.argv[0]) 
from math import *
import re,sys,decimal,os
rootpath=os.path.dirname(os.path.abspath(sys.argv[0])) 
print rootpath
sys.path.append(rootpath+'/convert') #add path to /convert
import lmf2ctrl

argvs = sys.argv
argc  = len(argvs)
print 'args= ',argvs,argc

if (argc != 2 or '--help' in argvs): 
	print ' == Convert POSCAR file(input file of lmf) to ctrl =='
	print '    usage: vasp2ctrl POSCAR_foofar    '
	print '          Then we have ctrls.foo.vasp2ctrl        '
	sys.exit(-1)
for ix in argvs:
	if 'POSCAR_' in ix: 
		print ' Read ', ix, '. It is converted to ctrl file now.'
		ext1=re.sub('POSCAR_','',ix)
		ext2=re.sub('.vasp','',ext1)
		break
	elif '.cif' in ix: 
		ext1 = re.sub('.cif','',ix)
		ext2 = re.sub('.vasp','',ext1)
		break
titleinput = 'ctrls.'+ext2

vaspread = open(argvs[1]).read().split('\n') 
plat1 = vaspread[2].split()
plat2 = vaspread[3].split()
plat3 = vaspread[4].split()
# print plat1
# print plat2
# print plat3
print ' '+vaspread[7]
#sys.exit()

alat_val = lmf2ctrl.vasp2ctrl_alat(vaspread) #unit
all_atom,NBAS_val = lmf2ctrl.vasp2ctrl_atomcount(vaspread)
#atom_list   = lmf2ctrl.vasp2ctrl_atom(vaspread,alat_val,NBAS_val)
atom_list = lmf2ctrl.vasp2ctrl_atom(vaspread,alat_val,NBAS_val,plat1,plat2,plat3)
lmf2ctrl.vasp2ctrl_write(vaspread,alat_val,NBAS_val,atom_list,titleinput,all_atom)

print ' OK! we have ',titleinput+'.vasp2ctrl'


