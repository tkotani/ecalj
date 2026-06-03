#!/usr/bin/env python3
import re,os,sys,decimal
from math import *
import convctrl
argvs = sys.argv
argc = len(argvs)

if (argc != 2 or '--help' in argvs): # help
	print( ' == Convert ctrl* / ctrlg.*.toml file (input file of lmf) to VASP POSCAR ==')
	print( '    usage: ctrl2vasp ctrlg.foo.toml          (canonical, post-2026-05)')
	print( '           ctrl2vasp ctrls.foo               (legacy structure-only)')
	print( '           ctrl2vasp ctrl.foo                (legacy text ctrl)')
	print( '           Then we have POSCAR_foo.vasp  (Cartesian only)')
	sys.exit(-1)
coordinates='Cartesian'
angstrom = 0.529177

inp = argvs[1]

# canonical TOML: ctrlg.<sname>.toml  -- parse with tomllib and feed the
# same (ALAT, PLAT_list, atomlist) tuple shape convctrl.savefile() wants.
if inp.startswith('ctrlg.') and inp.endswith('.toml'):
	ext = inp[len('ctrlg.'):-len('.toml')]
	titleinput = 'POSCAR_' + ext
	import tomllib
	with open(inp, 'rb') as f:
		d = tomllib.load(f)
	# struc.alat is in Bohr; convert to angstrom to match savefile output.
	ALAT = d['struc']['alat'] * angstrom
	# struc.plat is a 3x3 nested list -> flatten to 9 floats.
	PLAT_list = [v for row in d['struc']['plat'] for v in row]
	# [[site]] array -> [[name, x, y, z], ...].  savefile() groups
	# consecutive same-name entries into one POSCAR atom-type block,
	# matching the original ctrl-text behaviour.
	atomlist = [[s['atom'], *s['pos']] for s in d['site']]
	convctrl.savefile(ALAT, PLAT_list, atomlist, titleinput, coordinates)
	print(' OK! we have ', titleinput)
	sys.exit(0)

# legacy text ctrl / ctrls path
for ix in argvs:  #Get extensions for ctrl.ext or ctrls.ext
        print(ix)
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
print ('defined variables are:',variable_val)
T=1
F=0
for const in range(len(variable_val)):
        print(variable_val[const])
        exec( variable_val[const] )
variable_name = convctrl.keywordname(variable_val)

ALATone = convctrl.alat(perfectopen,variable_name)
ALATone[0][1] = eval(ALATone[0][1])
ALAT = ALATone[0][1]*angstrom # ALAT in angstrom
PLAT_list = convctrl.plat(perfectopen,variable_name)

for line_ing in range(len(PLAT_list)):
	PLAT_list[line_ing] = re.sub('/','/1.0/',PLAT_list[line_ing]) # this trick replace 5/3 with 5/1.0/3. In python 2 this gives difference (5/3=1 in python2).
	PLAT_list[line_ing] = eval(PLAT_list[line_ing])

atomlist = convctrl.atom(perfectopen,variable_name)
for line_all in range(len(atomlist)):
        for line_each in range(1,4): # / is replaced by /1.0/ to avoid "integer division" in python2.x
                atomlist[line_all][line_each] = re.sub('/','/1.0/',atomlist[line_all][line_each])
                atomlist[line_all][line_each] = eval(atomlist[line_all][line_each])
                if coordinates=='Cartesian' :
                        atomlist[line_all][line_each] = atomlist[line_all][line_each] # *ALAT new VASP


savefile = convctrl.savefile(ALAT,PLAT_list,atomlist,titleinput,coordinates)

print (' OK! we have ' ,titleinput)
