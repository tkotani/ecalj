#!/usr/bin/python

import numpy as np

print "Dviding the W interaction.\n"

################################################################
##### Input

nwf = 2
ns = 1				

#nwf : number of Wannier functions
#ns = 2 for spin polar, 1 for non-spin polar

#################################################################

f1 = open("Coulomb_v", 'r')
f2 = open("Screening_W-v", 'r')

lines1 = f1.readlines()
lines2 = f2.readlines()

f1.close()
f2.close()

v = np.empty((ns,nwf,nwf,nwf,nwf), dtype=object)
w = np.empty((ns,nwf,nwf,nwf,nwf), dtype=object)

for ins in range(0,ns):
	for i in range(0,nwf):
	      for j in range(0,nwf):
		    for k in range(0,nwf):
			  for l in range(0,nwf):
				v[ins][i][j][k][l] = '' 
				w[ins][i][j][k][l] = '' 


for line1 in lines1:
	v[int(line1.split()[1])-1][int(line1.split()[2])-1][int(line1.split()[3])-1][int(line1.split()[4])-1][int(line1.split()[5])-1] += line1 

for line2 in lines2:
	line_temp = line2.split() 
	line_temp[8] = str(float(line_temp[8]) + float(v[int(line2.split()[1])-1][int(line2.split()[2])-1][int(line2.split()[3])-1][int(line2.split()[4])-1][int(line2.split()[5])-1].split()[6]))
	w[int(line2.split()[1])-1][int(line2.split()[2])-1][int(line2.split()[3])-1][int(line2.split()[4])-1][int(line2.split()[5])-1] += " \t".join(line_temp) + "\n"

for ins in range(0,ns):
	for i in range(0,nwf):
	      for j in range(0,nwf):
		    for k in range(0,nwf):
			  for l in range(0,nwf):
				f3 = open("W"+str(ins+1)+str(i+1)+str(j+1)+str(k+1)+str(l+1)+".dat",'w')
				f3.write(w[ins][i][j][k][l])
				f3.close()

print "Dviding the W interaction was finished. W****.dat was generated.\n"
