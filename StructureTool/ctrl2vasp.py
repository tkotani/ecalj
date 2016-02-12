#!/usr/bin/env python
import re,os,sys,decimal
from math import *
#rootpath = os.path.abspath(os.path.dirname(sys.argv[0]))
#sys.path.append(rootpath+'/convert') # Add path rootpath/convert to read lmf2vasp.py
#import lmf2vasp
import convctrl
argvs = sys.argv
argc = len(argvs)
#print
#print 'args= ',argvs,argc

if (argc != 2 or '--help' in argvs): # help
	print ' == Convert ctrl file(input file of lmf) to VASP POSCAR =='
	print '    usage: ctrl2vasp ctrl.foo '
	print '           (ctrls.foo instead)' 
	print '           Then we have POSCAR_foo.vasp  (now for Cartesian only)   '
	sys.exit(-1)
#if ('--Cartesian' in argvs): 
#       coordinates='Cartesian'   
#elif ('--Direct' in argvs): 
#	coordinates='Direct'   
#else:
#	print ' You need to set --Cartesian or --Direct.' 
#	sys.exit(-1)
coordinates='Cartesian'   

for ix in argvs:  #Get extensions for ctrl.ext or ctrls.ext
	if 'ctrl.' in ix:
		ext = ix.split('ctrl.')[1]
		break
	elif 'ctrls.' in ix:
		ext = ix.split('ctrls.')[1]
		break

titleinput = 'POSCAR_'+ext

openfile = open(argvs[1]).read().split('\n') 
perfectopen = convctrl.fileopen(openfile)
variable_val = convctrl.constlist(openfile)
print 'defined variables are:',variable_val

for const in range(len(variable_val)):
	exec variable_val[const] 
variable_name = convctrl.keywordname(variable_val)

ALATone = convctrl.alat(perfectopen,variable_name)
ALATone[0][1] = eval(ALATone[0][1])
angstrom = 0.529177
ALAT = ALATone[0][1]*angstrom # ALAT in angstrom
PLAT_list = convctrl.plat(perfectopen,variable_name)
#print 'plat=',PLAT_list

for line_ing in range(len(PLAT_list)):
	PLAT_list[line_ing] = re.sub('/','/1.0/',PLAT_list[line_ing]) # this trick replace 5/3 with 5/1.0/3. In python 2 this gives difference (5/3=1 in python2).
	PLAT_list[line_ing] = eval(PLAT_list[line_ing])

atomlist = convctrl.atom(perfectopen,variable_name)
#print atomlist
for line_all in range(len(atomlist)):
	for line_each in range(1,4): # / is replaced by /1.0/ to avoid "integer division" in python2.x
		atomlist[line_all][line_each] = re.sub('/','/1.0/',atomlist[line_all][line_each])
		atomlist[line_all][line_each] = eval(atomlist[line_all][line_each])
		#print atomlist[line_all][line_each]
		#print "coordinates=",coordinates
		if coordinates=='Cartesian' :
			atomlist[line_all][line_each] = atomlist[line_all][line_each] # *ALAT new VASP
#		elif coordinates=='Direct' :
#			atomlist[line_all][line_each] = atomlist[line_all][line_each]
#		print atomlist[line_all][line_each]


savefile = convctrl.savefile(ALAT,PLAT_list,atomlist,titleinput,coordinates)

print ' OK! we have ' ,titleinput




