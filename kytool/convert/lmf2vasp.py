#lmf2vasp functions. for ctrl2vasp.
import re
import sys

def fileopen(openfile): #return cleanuped contents of openfile.
	newctrl = []
	newctrl2 = []
	for readctrl in range(len(openfile)):
		sharp_remove = openfile[readctrl].split('#')
		newctrl.append(sharp_remove[0]) # remove comment lines
	for i in range(len(newctrl)):
		eee = re.sub(r'\s*=\s* ', '=', newctrl[i]) # remove spaces sandwitch =.
		newctrl2.append(eee)
	return newctrl2

def constlist(openfile): #---get input values from ctrls (openfile)
	keywordlist = []
	for open_read in openfile:
		hx = re.match('^',open_read)
		open_read = re.sub('\^','**',open_read) # x^y means x**y
		rx = re.match('/', open_read)
		open_read = re.sub('/', '/1.0/', open_read) # python /1.0/
		fx = re.match(r'\s*=\s* ',  open_read) 
		open_read = re.sub(r'\s*=\s* ', '=', open_read)
		ex = re.match(r'^\%\s*const ',  open_read) 
		open_read = re.sub(r'^\%\s*const ', '', open_read)
		const_temp = open_read.split('#') # After '#', neglected.
		const_ex = const_temp[0]
		if ex: 
			const_val = const_ex.split(' ')
			for const_form in range(len(const_val)):
				if const_val[const_form] != "" :
					keywordlist.append(const_val[const_form])
	return keywordlist

def keywordname(keyword_val): #geve valus from keyword
	keyword_temp = []
	for const_name in range(len(keyword_val)):
		splitword = keyword_val[const_name].split('=')
		keyword_temp.append(splitword[0])
	return keyword_temp

def alat(perfectopen,keyword_name): #get alat (unit)
	for i in range(len(perfectopen)):
		if perfectopen[i].find("ALAT") != -1: #find alat
			break
	struc_rem = re.sub('STRUC' , '' ,perfectopen[i]) # remove "STRUC"
	alat_temp = struc_rem.split()
	ALATone = []
	for follow_struc in range(len(alat_temp)):
		alat_val = alat_temp[follow_struc].split('=') # ALAT=foobar --> alat,foobar
		ALATone.append(alat_val)
	for alat_first in range(len(ALATone)):
		for alat_second in range(len(ALATone[alat_first])):
			for const_name in range(len(keyword_name)):
				ALATone[alat_first][alat_second] = re.sub('\{'+keyword_name[const_name]+'\}',keyword_name[const_name],ALATone[alat_first][alat_second]) #replace {foobar} from ALAT
	return ALATone

def plat(perfectopen,keyword_name): #---get PLAT
	fileall=''
	for ix in perfectopen:
		fileall +=ix
#	print fileall	
#	print perfectopen.connect()
#	plat_temp = perfectopen[i].split('PLAT=')
	plat_temp = fileall.split('PLAT=')
	PLAT_list = plat_temp[1].split()[0:9]
#	print PLAT_list
#	sys.exit()

#	for i in range(len(perfectopen)):
#		if perfectopen[i].find("PLAT") != -1: #---find PLAT
#			break
#	plat_temp = perfectopen[i].split('=')
#	PLATone = []
#	for follow_struc in range(len(plat_temp)):
#		plat_val = plat_temp[follow_struc].split()
#		PLATone.append(plat_val)
#	PLAT_list = PLATone[-1] # last element constains PLAT.
	for line_ing in range(len(PLAT_list)):
		for const_name in range(len(keyword_name)):
			PLAT_list[line_ing] = re.sub('\{'+keyword_name[const_name]+'\}',keyword_name[const_name],PLAT_list[line_ing]) # replace {foobar} with its values
	return PLAT_list

def atom(perfectopen,keyword_name): #---get ATOM,POS
	temp1 = []
	temp3 = []
	temp5 = []
	for i in range(len(perfectopen)):
		temp2 = perfectopen[i].strip().split()
		temp1.append(temp2)
	for i in range(len(temp1)):
		for ii in range(len(temp1[i])):
			if re.match("POS", temp1[i][ii]) != None : # look for POS
				temp3.append(temp1[i])
	if re.match('SITE', temp3[0][0]) != None : # wether POS is in SITE or not.
		temp3[0].pop(0) 
	for i in range(len(temp3)):	
		for ii in range(len(temp3[i])):
			vvv1 = re.sub('ATOM=','',temp3[i][ii])
			vvv2 = re.sub('POS=','',vvv1) #---remove ATOM and POS
			temp5.append(vvv2)
	for line_each in range(len(temp5)):
		for const_name in range(len(keyword_name)):
			temp5[line_each] = re.sub('\{'+keyword_name[const_name]+'\}',keyword_name[const_name],temp5[line_each]) 
	sep_atom = []
	kaisuu = int(len(temp5)/4) 
	j1 = -1
	for j in range(kaisuu):
		j1 = j1+1
		j2 = 4*j1
		sep_list = temp5[j2:j2+4]
		sep_atom.append(sep_list) 
	return sep_atom

def savefile(ALAT,PLAT_list,atomlist,titleinput,coordinates): 
	title = titleinput
	output = title+'.vasp' 
	writefile = open(output,'w')
	outputlist1 = title+'\n'+str(ALAT)+'\n'
	writefile.write(outputlist1)
	for j in range(len(PLAT_list)):
		outputlist2 = str(PLAT_list[j])+' '
		writefile.write(outputlist2)
		if (j+1)%3 == 0: 
			writefile.write('\n')
	for i in range(len(atomlist)):
		outputlist3 = atomlist[i][0]+' ' 
		writefile.write(outputlist3)
	writefile.write('\n')
	for k in range(len(atomlist)):
		outputlist4 = '1'+' '
		writefile.write(outputlist4)
	outputlist5 = '\n'+coordinates+'\n' #---direct,cartesian
	writefile.write(outputlist5)
	for j in range(len(atomlist)):
		outputlist6 = str(atomlist[j][1])+' '+str(atomlist[j][2])+' '+str(atomlist[j][3])+'\n'
		writefile.write(outputlist6) 
	writefile.close()


