#!/usr/bin/env python
# -*- coding: utf-8 -*-

from math import *
import re
import sys,decimal
import os
os.path.abspath(sys.argv[0]) #---引数１(このプログラム)のpathを通す
rootpath=os.path.dirname(os.path.abspath(sys.argv[0])) #---引数１のディレクトリ名を特定し、pathを通す
print rootpath
sys.path.append(rootpath+'/convert') #---pathをサブディレクトリまで通して関数の引用可能にする
import lmf2ctrl #---関数記述のファイル名

argvs = sys.argv
argc = len(argvs)

print 'args= ',argvs,argc

if (argc != 2 or '--help' in argvs): #---help文
	print ' == Convert POSCAR file(input file of lmf) to ctrl =='
	print '    usage: vasp2ctrl POSCAR_foo.vasp    '
	print '          Then we have ctrl.foo        '
	sys.exit(-1)
for ix in argvs:
	if 'POSCAR_' in ix: #---条件分岐:ctrl2vaspで作成したPOSCAR形式への対応
		print ' Read ', ix, '. It is converted to ctrl file now.'
		ext1=re.sub('POSCAR_','',ix)
		ext2=re.sub('.vasp','',ext1)
		break
	elif '.cif' in ix: #---条件分岐：VESTAからExportしたcifファイルの対応
		ext1 = re.sub('.cif','',ix)
		ext2 = re.sub('.vasp','',ext1)
		break
titleinput = 'ctrls.'+ext2

vaspread = open(argvs[1]).read().split('\n') #---読み込むPOSCAR形式fileを段落ごとにリスト化

alat_val = lmf2ctrl.vasp2ctrl_alat(vaspread)

all_atom = lmf2ctrl.vasp2ctrl_atomcount(vaspread)

NBAS_val = lmf2ctrl.vasp2ctrl_NBAS(all_atom)

atom_list = lmf2ctrl.vasp2ctrl_atom(vaspread,alat_val,NBAS_val)

lmf2ctrl.vasp2ctrl_write(vaspread,alat_val,NBAS_val,atom_list,titleinput,all_atom)

print ' OK! we have ',titleinput #---書き込み完了メッセージ


