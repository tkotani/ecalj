#!/usr/bin/env python
# -*- coding: utf-8 -*-
VESTA='~/VESTA-x86_64/VESTA '
import sys
import os
import re
os.path.abspath(sys.argv[0]) #---pathの取得

argvs = sys.argv
argc = len(argvs)

print 'args= ',argvs,argc

if (argc != 2 or '--help' in argvs): #---help文
	print ' == Image Display VASP(input file of lmf) to ctrl =='
	print '    usage: viewvesta POSCAR_foo.vasp    '
	sys.exit(-1)

vestaopen=VESTA +os.getcwd()+'/'+argvs[1] #---vestaの起動path 保存先により要変更
print vestaopen

os.system(vestaopen) #---vesta起動

