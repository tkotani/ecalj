# -*- coding: utf-8 -*-

import re
import sys

def vasp2ctrl_alat(vaspread): #---vasp形式ファイルを読み込む
	alat_list = vaspread[1].split()
	alat_val = float( alat_list[0] )/0.529177 #---Å単位から戻す
	return alat_val	

def vasp2ctrl_atomcount(vaspread): #---原子名をリストに格納する
	all_atom = []
	atom_pat = vaspread[5].split()
	atom_num = vaspread[6].split()
	for number in range(len(atom_pat)):
		for quant in range( int(atom_num[number]) ):
			all_atom.append(atom_pat[number])
	return all_atom

def vasp2ctrl_NBAS(all_atom): #---原子の数からNBASをとる
	NBAS_val = len(all_atom)
	return NBAS_val


def vasp2ctrl_atom(vaspread,alat_val,NBAS_val): #---ATOMのPOSをとる
	atom_list = []
	angs = alat_val*0.529177
	for i in range(NBAS_val):
		atom_val = vaspread[8+i].split()
		atom_list.append(atom_val)
	for atom1 in range(len(atom_list)):
		for atom2 in range(len(atom_list[atom1])):
			atom_list[atom1][atom2] = float(atom_list[atom1][atom2])/angs #---alatで割って元に戻す
	return atom_list

def vasp2ctrl_write(vaspread,alat_val,NBAS_val,atom_list,titleinput,all_atom): #---ファイルへの書き込み
	outputname = titleinput+'.vasp2ctrl' #---拡張子の追加
	ctrlwrite = open(outputname,'w')
	tes = all_atom
	outputlist1 = 'STRUC'+'\n'+' '*5+'ALAT='+str(alat_val)+'\n' \
	+' '*5+'PLAT='+vaspread[2]+vaspread[3]+vaspread[4]+'\n' \
	+' '*2+'NBAS='+str(NBAS_val)+'\n'+'SITE'+'\n' #---write数回記述せずに一つにまとめて書く
	ctrlwrite.write(outputlist1)
	for sort in range(len(tes)):
		outputlist2 = ' '*5+'ATOM='+tes[sort]+' '+'POS='+str(atom_list[sort][0]) \
		+' '+str(atom_list[sort][1])+' '+str(atom_list[sort][2])+'\n'
		ctrlwrite.write(str(outputlist2))
	ctrlwrite.close()




