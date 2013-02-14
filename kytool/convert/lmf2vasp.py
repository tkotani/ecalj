# -*- coding: utf-8 -*-

import re
import sys

def fileopen(openfile): #---ctrlfileからコメントアウト文を取り除く関数
	newctrl = []
	newctrl2 = []
	for readctrl in range(len(openfile)):
		sharp_remove = openfile[readctrl].split('#')
		newctrl.append(sharp_remove[0]) #---コメント文を取り除く
	for i in range(len(newctrl)):
		eee = re.sub(r'\s*=\s* ', '=', newctrl[i]) #---=前後のスペースを取り除く
		newctrl2.append(eee)
	return newctrl2

def constlist(openfile): #---ctrlfileから変数取得の関数
	keywordlist = []
	for open_read in openfile:
		hx = re.match('^',open_read)
		open_read = re.sub('\^','**',open_read) #---^はevalで使えないため**へ置換
		rx = re.match('/', open_read)
		open_read = re.sub('/', '/1.0/', open_read) #---/が整数の場合0をとるので/1.0/へ置換
		fx = re.match(r'\s*=\s* ',  open_read) 
		open_read = re.sub(r'\s*=\s* ', '=', open_read)
		ex = re.match(r'^\%\s*const ',  open_read) 
		open_read = re.sub(r'^\%\s*const ', '', open_read) #---変数の取得
		const_temp = open_read.split('#')
		const_ex = const_temp[0]
		if ex: #---上記でヒットしたときの条件分岐
			const_val = const_ex.split(' ')
			for const_form in range(len(const_val)):
				if const_val[const_form] != "" :
					keywordlist.append(const_val[const_form])
	return keywordlist

def keywordname(keyword_val): #---変数から変数名をとる関数
	keyword_temp = []
	for const_name in range(len(keyword_val)):
		splitword = keyword_val[const_name].split('=')
		keyword_temp.append(splitword[0])
	return keyword_temp

def alat(perfectopen,keyword_name): #---ALATを取得する
	for i in range(len(perfectopen)):
		if perfectopen[i].find("ALAT") != -1: #---ALATを見つける
			break
	struc_rem = re.sub('STRUC' , '' ,perfectopen[i]) #---STRUCの除去
	alat_temp = struc_rem.split()
	ALATone = []
	for follow_struc in range(len(alat_temp)):
		alat_val = alat_temp[follow_struc].split('=') #---[ALAT＝OO]リストを[ALAT,OO]に
		ALATone.append(alat_val)
	for alat_first in range(len(ALATone)):
		for alat_second in range(len(ALATone[alat_first])):
			for const_name in range(len(keyword_name)):
				ALATone[alat_first][alat_second] = re.sub('\{'+keyword_name[const_name]+'\}',keyword_name[const_name],ALATone[alat_first][alat_second]) #---ALATに含まれる{}をとる
	return ALATone

def plat(perfectopen,keyword_name): #---PLATを取得する
	for i in range(len(perfectopen)):
		if perfectopen[i].find("PLAT") != -1: #---PLATを見つける
			break
	plat_temp = perfectopen[i].split('=')
	PLATone = []
	for follow_struc in range(len(plat_temp)):
		plat_val = plat_temp[follow_struc].split()
		PLATone.append(plat_val)
	PLAT_list = PLATone[-1] #---作るリストの構造上、最後の成分がPLATの値をとる
	for line_ing in range(len(PLAT_list)):
		for const_name in range(len(keyword_name)):
			PLAT_list[line_ing] = re.sub('\{'+keyword_name[const_name]+'\}',keyword_name[const_name],PLAT_list[line_ing]) #---PLATに含まれる{}をとる
	return PLAT_list

def atom(perfectopen,keyword_name): #---ATOM,POSの取得
	temp1 = []
	temp3 = []
	temp5 = []
	for i in range(len(perfectopen)):
		temp2 = perfectopen[i].strip().split()
		temp1.append(temp2)
	for i in range(len(temp1)):
		for ii in range(len(temp1[i])):
			if re.match("POS", temp1[i][ii]) != None : #---POSが含まれるリストの成分の検索
				temp3.append(temp1[i])
	if re.match('SITE', temp3[0][0]) != None : #---POSにSITEが含まれるときとそうでないときがあるため条件分岐
		temp3[0].pop(0) #---POSと同じ行にSITEがある場合は削除
	for i in range(len(temp3)):	
		for ii in range(len(temp3[i])):
			vvv1 = re.sub('ATOM=','',temp3[i][ii])
			vvv2 = re.sub('POS=','',vvv1) #---ATOM=とPOS=の削除
			temp5.append(vvv2)
	for line_each in range(len(temp5)):
		for const_name in range(len(keyword_name)):
			temp5[line_each] = re.sub('\{'+keyword_name[const_name]+'\}',keyword_name[const_name],temp5[line_each]) #---POSに含まれる{}をとる
	sep_atom = []
	kaisuu = int(len(temp5)/4) #---以下の作業：ATOM,POS共に１つのリストに入るため４つずつリストに分ける
	j1 = -1
	for j in range(kaisuu):
		j1 = j1+1
		j2 = 4*j1
		sep_list = temp5[j2:j2+4]
		sep_atom.append(sep_list) #---[ATOMname,POS1,POS2,POS3]のリストがATOMの数だけ出来る
	return sep_atom

def savefile(ALAT,PLAT_list,atomlist,titleinput,coordinates): #---ファイルへの書き出し
	title = titleinput
	output = title+'.vasp' #---vesta対応の拡張子追加
	writefile = open(output,'w')
	outputlist1 = title+'\n'+str(ALAT)+'\n'
	writefile.write(outputlist1)
	for j in range(len(PLAT_list)):
		outputlist2 = str(PLAT_list[j])+' '
		writefile.write(outputlist2)
		if (j+1)%3 == 0: #---PLATは３つずつ改行挟んでの書き込み
			writefile.write('\n')
	for i in range(len(atomlist)):
		outputlist3 = atomlist[i][0]+' ' #---ATOMの書き込み
		writefile.write(outputlist3)
	writefile.write('\n')
	for k in range(len(atomlist)):
		outputlist4 = '1'+' '
		writefile.write(outputlist4)
	outputlist5 = '\n'+coordinates+'\n' #---direct,cartesianの表記
	writefile.write(outputlist5)
	for j in range(len(atomlist)):
		outputlist6 = str(atomlist[j][1])+' '+str(atomlist[j][2])+' '+str(atomlist[j][3])+'\n'
		writefile.write(outputlist6) #---上記でPOSを３つずつ改行挟んで書き込み
	writefile.close()


