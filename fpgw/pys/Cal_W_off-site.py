#!/usr/bin/python
import numpy as np

def w_name(irws, isp,i,j,k,l,lcrpa):
	if not lcrpa:
        	filename ="W"+str(irws)+"_"+str(isp)+str(i)+str(j)+str(k)+str(l)+".dat"
	elif lcrpa:
        	filename ="W_cRPA"+str(irws)+"_"+str(isp)+str(i)+str(j)+str(k)+str(l)+".dat"
        return filename

def write_W(nrws,prws,nsp,nwf,lines,lcrpa):
	for irws in range(1,nrws+1):
		for isp in range(1,nsp+1):
			for i in range(1,nwf+1):
				for j in range(1,nwf+1):
					for k in range(1,nwf+1):
						for l in range(1,nwf+1):
							if i == j == k == l:
								f = open(w_name(irws,isp,i,j,k,l,lcrpa),'w')
								f.write("")
								f.close()

	for line in lines:
		line_sp = line.split()
		
		irws = int(line_sp[prws])
		isp = int(line_sp[6])
		i = int(line_sp[7])
		j = int(line_sp[8])
		k = int(line_sp[9])
		l = int(line_sp[10])

		if i == j == k == l:

			line_sp[13] = str(float(line_sp[13]) + v[irws-1][isp-1][i-1][j-1][k-1][l-1])
			if i == 1 and j == 1 and k == 1 and l == 1:
				ilines = 0

			f = open(w_name(irws,isp,i,j,k,l,lcrpa),'a')
			line_join = " \t".join(line_sp) + "\n"
			f.write(line_join)
			f.close()
		ilines += 1
		#print line_sp

def write_static_W(nrws,nsp,nwf,lcrpa):
	Contents = ""
	for irws in range(1,nrws+1):
		for isp in range(1,nsp+1):
			for i in range(1,nwf+1):
		    		for j in range(1,nwf+1):
					if i==j:
		 				f1 = open(w_name(irws,isp,i,i,i,i,lcrpa),'r')
						Contents += "U : <nn|W|nn>\t" + f1.readlines()[0]
		Contents += "\n"
	#		    for j in range(1,nwf+1):
	#		      if i!=j:
	#			  f2 = open(w_name(isp,i,i,j,j,lcrpa),'r')
	#			  Contents += "U' : <nn|W|mm>\t" + f2.readlines()[0]
	#			  f3 = open(w_name(isp,i,j,i,j,lcrpa),'r')
	#			  Contents += "J : <nm|W|mn>\t" + f3.readlines()[0]
	#			  f4 = open(w_name(isp,i,j,j,i,lcrpa),'r')
	#			  Contents += "J' : <nm|W|nm>\t" + f4.readlines()[0]
	#		    Contents += "\n"

	if not lcrpa:
		f5 = open("Static_W.dat",'w')
	if lcrpa:
		f5 = open("Static_U.dat",'w')
	f5.write(Contents)
	f5.close()

	f1.close()
#	if nwf >= 2:
#		f2.close()
#		f3.close()
#		f4.close()

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
f2 = open("Screening_W-v",'r')
f3 = open("Screening_W-v_crpa",'r')
lines1 = f1.readlines() 
lines2 = f2.readlines() 
lines3 = f3.readlines() 
f1.close()
f2.close()
f3.close()

line_sp = lines1[-1].split()

nrws1 = int(line_sp[1])
nrws2 = int(line_sp[2])

print "nrws1 = ", nrws1
print "nrws2 = ", nrws2

if nrws1 != 1:
	nrws = nrws1
	prws = 1
elif nrws2 != 1:
	nrws = nrws2
	prws = 2

############### Reading Coulomb V values
print "Reading Coulomb interaction V values...\n"

# spin up
v = np.zeros((nrws,nsp,nwf,nwf,nwf,nwf))

#print v

for line1 in lines1:
      line1_sp = line1.split() 
      irws = int(line1_sp[prws])
      isp = int(line1_sp[6])
      iwf1 = int(line1_sp[7])
      iwf2 = int(line1_sp[8])
      iwf3 = int(line1_sp[9])
      iwf4 = int(line1_sp[10])
      v[irws-1][isp-1][iwf1-1][iwf2-1][iwf3-1][iwf4-1] = float(line1_sp[11])

#
################ Writing screened interaction W values 
#print "Writing screened Coulomb interaction W values..."
#lcrpa = False
#write_W(nrws,prws,nsp,nwf,lines2,lcrpa)
#print "=> W****.dat files were generated!\n"
#
#
#
################## Writing constrained RPA W : W_cRPA values
#
#print "Writing cRPA W values..."
#lcrpa = True
#write_W(nrws,prws,nsp,nwf,lines3,lcrpa)
#print "=> W_cRPA****.dat files were generated!\n"


################# Extracting static U, U', J and J' values

print "Extracting static U, U', J and J' values..."

lcrpa = False
write_static_W(nrws,nsp,nwf,lcrpa)
print "=> Static_W.dat file was generated!"

lcrpa = True
write_static_W(nrws,nsp,nwf,lcrpa)
print "=> Static_U.dat file was generated!"

