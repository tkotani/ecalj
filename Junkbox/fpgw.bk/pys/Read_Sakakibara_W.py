#!/usr/bin/python
import os, sys

############### Input ########

system = "la2cuo4"
nwf = 2
ns = 1

#nwf : number of Wannier functions
#ns = 2 for spin polar, 1 for non-spin polar

###############################

Contents = ""

for i in range(0,ns):
        for j in range(0,nwf):
            for k in range(0,nwf):
              if j==k:
                  f1 = open("W"+str(i+1)+str(j+1)+str(j+1)+str(j+1)+str(j+1)+".dat",'r')
                  Contents += "U : <nn|W|nn>\t" + f1.readlines()[0]
              else:
                  f4 = open("W"+str(i+1)+str(j+1)+str(j+1)+str(k+1)+str(k+1)+".dat",'r')
                  Contents += "U' : <nn|W|mm>\t" + f4.readlines()[0]
                  f3 = open("W"+str(i+1)+str(j+1)+str(k+1)+str(j+1)+str(k+1)+".dat",'r')
                  Contents += "J : <nm|W|mn>\t" + f3.readlines()[0]
                  f2 = open("W"+str(i+1)+str(j+1)+str(k+1)+str(k+1)+str(j+1)+".dat",'r')
                  Contents += "J' : <nm|W|nm>\t" + f2.readlines()[0]
            Contents += "\n"

f5 = open("Sakakibara_W.dat",'w')
f5.write(Contents)
f5.close()


f6 = open("plot_UJ_"+system+".eps.gnu",'w')
contents_eps = """set term postscript eps font "Arial, 30" size 4,3 enhanced color
set output "UJ_%s.eps"
set key font "Arial, 20" #at "left top"
set key samplen 0.7
set xlabel "Frequency (eV)"
set ylabel "W (eV)"
set xzeroaxis 3
set mxtics 
set mytics
set xrange [0:50]
#set ytics 5
set yrange [-5:25]
 plot \
 "W11111.dat" u 8:9 lt 1 lc rgb "black" lw 3 pt 1 w l t "<d_{z^2}d_{z^2}|W|d_{z^2}d_{z^2}>",\\
 "W12222.dat" u 8:9 lt 1 lc rgb "red" lw 3 pt 1 w l t "<d_{x^2-y^2}d_{x^2-y^2}|W|d_{x^2-y^2}d_{x^2-y^2}>"\\
# "W11111.dat" u 8:($10*13.605) lt 1 lc rgb "black" lw 3 pt 1 w l t "<d_{z^2}d_{z^2}|W|d_{z^2}d_{z^2}>"\\
# "W12222.dat" u 8:($10*13.605) lt 1 lc rgb "red" lw 3 pt 1 w l t "<d_{x^2-y^2}d_{x^2-y^2}|W|d_{x^2-y^2}d_{x^2-y^2}>"
""" % (system)
f6.write(contents_eps)
f6.close()

os.system("gnuplot plot_UJ_"+system+".eps.gnu")
