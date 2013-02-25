import re
import sys
def vasp2ctrl_alat(vaspread): #read alat from VASP file
	alat_list = vaspread[1].split()
	alat_val = float( alat_list[0] )/0.529177 # convert it to a.u. from angstrome.
	return alat_val	

def vasp2ctrl_atomcount(vaspread): # set atom namet in all_atom[]
	all_atom = []
	atom_pat = vaspread[5].split()
	atom_num = vaspread[6].split()
	for number in range(len(atom_pat)):
		for quant in range( int(atom_num[number]) ):
			all_atom.append(atom_pat[number])
	NBAS_val = len(all_atom)
	return all_atom,NBAS_val

def vasp2ctrl_write(vaspread,alat_val,NBAS_val,atom_list,titleinput,all_atom):
	outputname = titleinput+'.vasp2ctrl'
	ctrlwrite = open(outputname,'w')
	tes = all_atom
	outputlist1 = 'STRUC'+'\n'+' '*5+'ALAT='+str(alat_val)+'\n' \
	+' '*5+'PLAT='+vaspread[2]+'\n'+' '*10+vaspread[3]+'\n'+' '*10+vaspread[4]+'\n' \
	+' '*2+'NBAS='+str(NBAS_val)+'\n'+'SITE'+'\n'
	ctrlwrite.write(outputlist1)
	for sort in range(len(tes)):
		outputlist2 = ' '*5+'ATOM='+tes[sort]+' '+'POS='+str(atom_list[sort][0]) \
		+' '+str(atom_list[sort][1])+' '+str(atom_list[sort][2])+'\n'
		ctrlwrite.write(str(outputlist2))
	ctrlwrite.close()

def vasp2ctrl_atom(vaspread,alat_val,NBAS_val,plat1,plat2,plat3):
	atom_list = []
	angs = alat_val*0.529177
	for i in range(NBAS_val):
		atom_val = vaspread[8+i].split()
		atom_list.append(atom_val)
	plat1= [float(plat1[i]) for i in range(len(plat1)) ] 
	plat2= [float(plat2[i]) for i in range(len(plat2)) ] 
	plat3= [float(plat3[i]) for i in range(len(plat3)) ] 
	for atom1 in range(len(atom_list)):
		atom_list[atom1]= [float(atom_list[atom1][i]) for i in range(3)]
		#print 'xxxxxxxxxxxxx',vaspread[7],atom_list[atom1]
		if vaspread[7]=='Direct':
			for ix in range(3):
				atom_list[atom1][ix] = atom_list[atom1][0]*plat1[ix] \
				+ atom_list[atom1][1]*plat2[ix] + atom_list[atom1][2]*plat3[ix]
	return atom_list
