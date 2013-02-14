#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
from math import *
import os
import sys,decimal
os.path.abspath(sys.argv[0]) #---引数１(このプログラム)のpathを通す
rootpath = os.path.abspath(os.path.dirname(sys.argv[0])) #---引数１のディレクトリ名を特定し、pathを通す
print rootpath
sys.path.append(rootpath+'/convert') #---pathをサブディレクトリまで通して関数の引用可能にする
import lmf2vasp

argvs = sys.argv
argc = len(argvs)
print

print 'args= ',argvs,argc


if (argc != 3 or '--help' in argvs): #---helpを設置
	print ' == Convert ctrl file(input file of lmf) to Vasp POSCAR =='
	print '    usage: ctrl2vasp ctrl.foo [--Cartesian or --Direct]   '
	print '           Then we have POSCAR_foo.vasp        '
	sys.exit(-1)
if ('--Cartesian' in argvs): #---引数にcartesianをとるときに文字を関連付ける
	coordinates='Cartesian'   
elif ('--Direct' in argvs): #---引数にdirectをとるときに文字を関連付ける
	coordinates='Direct'   
else:
	print ' You need to set --Cartesian or --Direct.' #---errorメッセージ   
	sys.exit(-1)

for ix in argvs: #---入力ファイルが2通りあるので対応可能
	if 'ctrl.' in ix:
		ext = ix.split('ctrl.')[1] #---ctrl.OOのOO部分を取り出す
		break
	elif 'ctrls.' in ix:
		ext = ix.split('ctrls.')[1] #---ctrls.OOへの対応
		break

titleinput = 'POSCAR_'+ext

openfile = open(argvs[1]).read().split('\n') #---入力fileを改行でsplitして読み込む

perfectopen = lmf2vasp.fileopen(openfile)

keyword_val = lmf2vasp.constlist(openfile)
for const in range(len(keyword_val)):
	exec keyword_val[const] #---evalを使うので関数外で実行可能にする

keyword_name = lmf2vasp.keywordname(keyword_val)

ALATone = lmf2vasp.alat(perfectopen,keyword_name)
ALATone[0][1] = eval(ALATone[0][1])
angstrom = 0.529177
ALAT = ALATone[0][1]*angstrom #---定数をかけてオングストローム単位にする

PLAT_list = lmf2vasp.plat(perfectopen,keyword_name)
for line_ing in range(len(PLAT_list)):
	PLAT_list[line_ing] = eval(PLAT_list[line_ing])

atomlist = lmf2vasp.atom(perfectopen,keyword_name)
for line_all in range(len(atomlist)):
		for line_each in range(1,4): #---整数で割ると結果が0になるため/を/1.0/に変更
			atomlist[line_all][line_each] = re.sub('/','/1.0/',atomlist[line_all][line_each])
			atomlist[line_all][line_each] = eval(atomlist[line_all][line_each])
			if coordinates=='Cartesian' :
				atomlist[line_all][line_each] = atomlist[line_all][line_each]*ALAT
			elif coordinates=='Direct' :
				atomlist[line_all][line_each] = atomlist[line_all][line_each]


savefile = lmf2vasp.savefile(ALAT,PLAT_list,atomlist,titleinput,coordinates)

print ' OK! we have ' ,titleinput #---書き込み完了メッセージ




