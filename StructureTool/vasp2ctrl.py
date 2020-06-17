#!/usr/bin/env python3
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
#print 'readin args= ',argvs,argc
if ('--help' in argvs or argc == 1): 
        print()
        print( ' == Convert POSCAR file (vasp format) to ctrls (ecalj format) ==')
        print( '  Usage: vasp2ctrl POSCARfoobar   [option] ')
        print( '         POSCARfoobar is the file name in POSCAR format (vasp5)')
        print( '         Then we have ctrls.POSCARfoobar.vasp2ctrl     ')
        print( '  Option: --alat=10.66 (Set ALAT in a.u.; Use eqs such as 10.66*1.88)')
        print( '         To get clean ctrl file, set --alat= as a lattice constant in a.u.')
        print( '')
        sys.exit(-1)
#for ix in argvs:
#	if '.vasp' in ix: 
#		print( ' Read ', ix, '. It is converted to ctrl file now.'
ext2=argvs[1]
#		ext1=re.sub('POSCAR_','',ix)
#		ext2=re.sub('.vasp','',ext1)
#		break
	# elif '.cif' in ix: 
	# 	ext1 = re.sub('.cif','',ix)
	# 	ext2 = re.sub('.vasp','',ext1)
	# 	break
for ix in argvs:
        if '--alat=' in ix: 
#                alatin=string.atof(eval(re.sub('--alat=','',ix) ))
                alatin=float(eval(re.sub('--alat=','',ix) ))
        else:
                alatin=None
                ratioa=1.00000000000
titleinput = 'ctrls.'+ext2
print(ext2,' -> ', titleinput)
vaspread = open(argvs[1]).read().split('\n') 
plat1 = vaspread[2].split()
plat2 = vaspread[3].split()
plat3 = vaspread[4].split()
print( ' '+vaspread[7])

alat_val = convctrl.vasp2ctrl_alat(vaspread) #unit
if(alatin): ratioa = alatin/alat_val #conversion by given --alat=alatin
all_atom,NBAS_val = convctrl.vasp2ctrl_atomcount(vaspread)
#atom_list   = convctrl.vasp2ctrl_atom(vaspread,alat_val,NBAS_val)
atom_list = convctrl.vasp2ctrl_atom(vaspread,alat_val,NBAS_val,plat1,plat2,plat3)
convctrl.vasp2ctrl_write(vaspread,alat_val,NBAS_val,atom_list,titleinput,all_atom,ratioa)

print( ' OK! we have ',titleinput+'.vasp2ctrl')


