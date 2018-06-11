#!/usr/bin/python
import numpy as np
import sys as sys
import os as os
import shutil as shutil
import glob as glob

### extracting the Hubbard interaction by cRPA ###
def read_cRPA_U(nwf,lines):
	print "making cRPA Interactions(final data)....."
	ilines = 0	
	for line in lines:		
		line_sp = line.split()		
		isp = int(line_sp[6])
		i = int(line_sp[7])
		j = int(line_sp[8])
		k = int(line_sp[9])
		l = int(line_sp[10])		
		omega = float(line_sp[12]) #[11] maybe eV, [12] is hartree=2*13.6eV
	
	        if omega == 0.0:
			v = v1[ilines]

			if i == j == k == l:
				U[i-1][i-1] =  float(line_sp[13]) + v
			elif i == j and k == l:
				U[i-1][k-1] = float(line_sp[13]) + v
			elif i == k and j == l:
				J[i-1][j-1] = float(line_sp[13]) + v
			elif i == l and k == j:
				J[i-1][k-1] = float(line_sp[13]) + v
			ilines += 1
		else :
			break		
		#print line_sp

### extracting the effective interaction for mRPA ###
def read_mRPA_W_UP(numwan,lines,v1): #
	print "making mRPA Inputs(seed data)....."
	fw1 = open('W0.tens.mRPA','w')
	fw2 = open('W0.mat.mRPA','w')
	ilines = 0 # lines distinguishing the orbital 
	for line in lines:
		dumline = line.split()
		i = int(dumline[7])
		j = int(dumline[8])
		k = int(dumline[9])
		l = int(dumline[10])
		m1 = j+numwan*(i-1)
		m2 = l+numwan*(k-1)

		if float(dumline[12]) == 0.0 : # omega=0(static)
			v= v1[ilines]
			Wout=float(dumline[13]) +v
			fw1.write("%5i %5i %5i %5i %13.7f\n" % (i,j,k,l,Wout))
			fw2.write("%5i %5i %13.7f\n" % (m1,m2,Wout))
			ilines += 1
	fw1.close()			   
	fw2.close()


########## make flags
mRPA="off"
cRPA="off"

if os.path.exists("hrotr.up") and os.path.exists("LMTO") :
	print "serching input files...."
else :
	print "error!!! hrotr.up or LMTO does not exist!!!"
	sys.exit

if os.path.exists("Screening_W-v.UP") and os.path.exists("Coulomb_v.UP") :
	mRPA="on"
	print "mRPA is set up"

if os.path.exists("Screening_W-v_crpa.UP") and os.path.exists("Coulomb_v.UP") :
        cRPA="on"
	print "cRPA is set up"

######### Reading nwf from 'hrotr.up' ########   

#f = open("GWinput",'r')
f = open("hrotr.up",'r')
lines = f.readlines()
#lread = False
for line in lines:
	line_sp = line.split()
	nwf = int(line_sp[1]) 	
	break
#	if "<MLWF>" in line:
#		lread = True
#	elif lread == True:
#		line_sp = line.split()
#		nwf = int(line_sp[0])
#		break
f.close()

######### Reading spin from 'LMTO' ##########

f = open("LMTO",'r')
lines = f.readlines()
lread = False
for line in lines:
	if "spin" in line:
		lread = True
	elif lread == True:
		line_sp = line.split()
		nsp = int(line_sp[0])
		break
f.close()

print "nwf = ", nwf
print "nsp = ", nsp
	

########## Opening Files########

if os.path.exists("Coulomb_v.UP") :
	f1 = open("Coulomb_v.UP",'r')
	lines1 = f1.readlines()
	f1.close()

if os.path.exists("Screening_W-v.UP") :
	f2 = open("Screening_W-v.UP",'r')
	lines2 = f2.readlines()
	f2.close()

if os.path.exists("Screening_W-v_crpa.UP") :
	f3 = open("Screening_W-v_crpa.UP",'r')
	lines3 = f3.readlines() 
	f3.close()

############### Reading Coulomb V values
print "Reading Coulomb interaction V values...\n"
# spin up
v1 = np.zeros((nwf*nwf*nwf*nwf))
# spin down
v2 = np.zeros((nwf*nwf*nwf*nwf))

U = np.zeros((nwf,nwf))
J = np.zeros((nwf,nwf))

ilines1 = 0
ilines2 = 0

for line1 in lines1:
	line1_sp = line1.split() 
	if line1_sp[1] == '1':
		v1[ilines1] = float(line1_sp[11])
		ilines1 += 1
	elif line1_sp[1] == '2':
		v2[ilines2] = float(line1_sp[11])
		ilines2 += 1

############### Writing the mRPA interaction value for FLEX
if mRPA == "on" :
	read_mRPA_W_UP(nwf,lines2,v1)
	if os.path.exists("W0.mat.mRPA") :
		print "=> OK. 'W0.mat.mRPA' is generated!"

############### Writing the cRPA interaction input for FLEX code

if cRPA == "on" :
	if os.path.exists("Interaction_cRPA.dat") :
		if os.path.exists("Interaction_cRPA.dat-old") :
			shutil.copyfile("Interaction_cRPA.dat","Interaction_cRPA.dat-old2")
			print "[waring] Interaction_cRPA.dat is copied to Interaction_cRPA.dat-old2!"
		else :
			shutil.copyfile("Interaction_cRPA.dat","Interaction_cRPA.dat-old")
                        print "[waring] Interaction_cRPA.dat is copied to Interaction_cRPA.dat-old1!"


	fout = open("Interaction_cRPA.dat",'w')
	fout.write("System_name\n")
	fout.write("%.10s %6i\n" % ("num_wan", nwf))
	fout.write("%.10s %6i\n" % ("spin", nsp))

	print "Reading cRPA U values..."

	read_cRPA_U(nwf,lines3)
	fout.write("%.10s" % ("Coulomb\n"))

	for iwf1 in range(0,nwf):
		for iwf2 in range(0,nwf):
			fout.write("%10.4f" % (U[iwf1][iwf2]))
		fout.write("\n")
	fout.write("\n")
			
	fout.write("%.10s\n" % ("Exchange"))

	for iwf1 in range(0,nwf):
		for iwf2 in range(0,nwf):
			fout.write("%10.4f" % (J[iwf1][iwf2]))
		fout.write("\n")
	fout.write("\n")
	if os.path.exists("Interaction_cRPA.dat") :
		print "=> OK. 'Interaction_cRPA.dat' is generated!"

