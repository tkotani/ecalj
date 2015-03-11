#!/usr/bin/env python
import re,os,commands,subprocess
aaa=\
"""
set terminal postscript enhanced color eps
set output "massfit.00A00EPS.eps"
set grid
set ylabel "mass/mass(electron)"
# Fitting region
set xrange [00A00EMIN: 00A00EMAX] 
 f(x) = a+ x**2*(b + x**2*c) 
 fit f(x) "temp" u ($5):(($8+$16)/2.0) via a,b, c
# Plotting region
set xrange [0.0: 0.12]
 plot \
 "00A00F1" u ($5):($8) lt 5 pt  5 w lp ti "band= 17",\
 "00A00F2" u ($5):($8) lt 6 pt  6 w lp ti "band= 18",\
 f(x) 
print 'mass=',f(0)
"""

############ block 1 ################# mhh L line 
file1= "Band005Syml004Spin1.mass"
file2= "Band006Syml004Spin1.mass"
dat='mhh_L'

EMIN="0.03"
EMAX="0.08"
temp='temp'+file1+file2
aaa=re.sub('00A00EPS',dat,aaa)
aaa=re.sub('00A00F1' ,file1 ,aaa)
aaa=re.sub('00A00F2' ,file2 ,aaa)
aaa=re.sub('00A00EMIN' ,EMIN,aaa)
aaa=re.sub('00A00EMAX' ,EMAX,aaa)
aaa=re.sub('00A00temp' ,temp,aaa)

os.system('paste '+ file1 +' '+file2 +' > '+ temp)
fff= 'mfit.'+dat+'.glt'
f = open(fff,'w')
f.write(aaa)
f.close
print fff+ 'for gnuplot'
#ccc='gnuplot '+ fff 
#os.system('which gnuplot')
#rrr=commands.getoutput(ccc)
#subprocess.call(['gnuplot', fff])
