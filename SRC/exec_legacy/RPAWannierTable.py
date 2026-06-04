#!/usr/bin/env python
#  S.W.Jang Apr2015. (need python2.x)
#  Run this after you did genMLWF. See Wannier/REAMDE and/or ecaljmanual.
#  You may need to modify this if necessary for your purpose.

def w_name(isp,i,j,k,l,lcrpa):
	if not lcrpa:
        	filename ="W"+str(isp)+str(i)+str(j)+str(k)+str(l)+".dat"
	elif lcrpa:
        	filename ="W_cRPA"+str(isp)+str(i)+str(j)+str(k)+str(l)+".dat"
        return filename

def write_W(nsp,nwf,lines,lcrpa):
	for isp in range(1,nsp+1):
		for i in range(1,nwf+1):
			for j in range(1,nwf+1):
				for k in range(1,nwf+1):
					for l in range(1,nwf+1):
						f = open(w_name(isp,i,j,k,l,lcrpa),'w')
						f.write("")
						f.close()
	for line in lines:
		line_sp = line.split()

		isp = line_sp[1]
		i = line_sp[2]
		j = line_sp[3]
		k = line_sp[4]
		l = line_sp[5]

		if i == '1' and j == '1' and k == '1' and l == '1':
			ilines = 0
		line_sp[8] = str(float(line_sp[8]) + v[ilines])

		f = open(w_name(isp,i,j,k,l,lcrpa),'a')
		line_join = " \t".join(line_sp) + "\n"
		f.write(line_join)
		f.close()
		ilines += 1
		#print line_sp

def write_static_W(nsp,nwf,lcrpa):
	Contents = ""
	for isp in range(1,nsp+1):
		for i in range(1,nwf+1):
		    for j in range(1,nwf+1):
		      if i==j:
			  f1 = open(w_name(isp,i,i,i,i,lcrpa),'r')
			  Contents += "U : <nn|W|nn>\t" + f1.readlines()[0]
		    for j in range(1,nwf+1):
		      if i!=j:
			  f2 = open(w_name(isp,i,i,j,j,lcrpa),'r')
			  Contents += "U' : <nn|W|mm>\t" + f2.readlines()[0]
			  f3 = open(w_name(isp,i,j,i,j,lcrpa),'r')
			  Contents += "J : <nm|W|mn>\t" + f3.readlines()[0]
			  f4 = open(w_name(isp,i,j,j,i,lcrpa),'r')
			  Contents += "J' : <nm|W|nm>\t" + f4.readlines()[0]
		    Contents += "\n"

	if not lcrpa:
		f5 = open("Static_W.dat",'w')
	if lcrpa:
		f5 = open("Static_U.dat",'w')
	f5.write(Contents)
	f5.close()
	f1.close()
	f2.close()
	f3.close()
	f4.close()

### MAIN #######################3
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
print " --- Summarize data on RPA and cRPA --- "
print " NOTE:nwf is from GWinput. nsp is from LMTO "
print " nwf (# of Wannier fun,) = ", nwf
print " nsp (# of spin )        = ", nsp

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

############### Reading Coulomb V values
print "Reading Coulomb interaction V values...\n"
v = [0]*nwf*nwf*nwf*nwf
#print v
ilines1 = 0
for line1 in lines1:
      line1_sp = line1.split() 
      v[ilines1] = float(line1_sp[6])
      ilines1 += 1

############### Writing screened interaction W values 
print "Writing screened Coulomb interaction W values..."
lcrpa = False
write_W(nsp,nwf,lines2,lcrpa)
print "=> W****.dat files were generated!\n"

################# Writing constrained RPA W : W_cRPA values
print "Writing cRPA W values..."
lcrpa = True
write_W(nsp,nwf,lines3,lcrpa)
print "=> W_cRPA****.dat files were generated!\n"


################# Extracting static U, U', J and J' values
print "Extracting static U, U', J and J' values..."
lcrpa = False
write_static_W(nsp,nwf,lcrpa)
print "=> Static_W.dat file was generated!"
lcrpa = True
write_static_W(nsp,nwf,lcrpa)
print "=> Static_U.dat file was generated!"

