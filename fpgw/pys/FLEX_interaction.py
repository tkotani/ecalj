#!/usr/bin/python
import numpy as np

def read_cRPA_U(nwf,lines):
	ilines = 0
	for line in lines:
		line_sp = line.split()

		isp = int(line_sp[1])
		i = int(line_sp[2])
		j = int(line_sp[3])
		k = int(line_sp[4])
		l = int(line_sp[5])
		omega = float(line_sp[7])
	        if omega == 0.0:
			v = v1[ilines]

			if i == j == k == l:
				U[i-1][i-1] =  float(line_sp[8]) + v
			elif i == j and k == l:
				U[i-1][k-1] = float(line_sp[8]) + v
			elif i == k and j == l:
				J[i-1][j-1] = float(line_sp[8]) + v
			elif i == l and k == j:
				J[i-1][k-1] = float(line_sp[8]) + v
			ilines += 1
		else :
			break
		#print line_sp


########## Reading spin from LMTO, nwf from GWinput

f = open("GWinput",'r')
lines = f.readlines()
lread = False
for line in lines:
	if "<MLWF>" in line:
		lread = True
	elif lread == True:
		line_sp = line.split()
		nwf = int(line_sp[0])
		break
f.close()

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

	

#######################

f1 = open("Coulomb_v",'r')
f3 = open("Screening_W-v_crpa",'r')
lines1 = f1.readlines() 
lines3 = f3.readlines() 
f1.close()
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
		v1[ilines1] = float(line1_sp[6])
		ilines1 += 1
	elif line1_sp[1] == '2':
		v2[ilines2] = float(line1_sp[6])
		ilines2 += 1



################# Writing constrained RPA W : W_cRPA values


############### Writing the interaction input for FLEX code
fout = open("Interaction.dat",'w')

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

print "=> OK. 'Interaction.dat' is generated!"
