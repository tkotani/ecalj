#!/usr/bin/env python
import os,re,sys
#echo $argv[1] #$argv[2] #$argv[3] $argv[4]
#echo $argv[2]
xx='7'  
yy='13' 
zz='1'
#dirname = `pwd`
#$dirname
#set xmax = 1
#set zz = argv[1]
#set xmax = $argv[4]
#ETOTeV.fc_${id}_222_rmt.800 dat1
#ETOTeV.fc_${id}_222_rmt.850 dat2
#set sft1 = 0
#1.4
#set id2 = 11
#sed id0 = $argv[1]
#pldimer1 > datg
gggtitle1 = open("gggtitle1",'r').read().split('\n')
eref=open('gggeref','r').read().split('\n')[0]
arglines=open('gggargs','r').read().split('\n')

argso=set([])
datg=sys.argv[1]
#print args

gggheader="""#!/usr/bin/gnuplot -persist
#set xrange [0: 1]
set encoding iso
set xlabel "Distance/\305"
set ylabel "Energy/Ry. +"""
gggheader=gggheader + str(eref) +'.0"\n'

gggtail="""
set term postscript enhanced color
set output 'pl.eps."""+datg+"'"+ \
"""
replot
set term png 
set output 'pl.png."""+datg+"'\n"+'replot\n'


i=0
datl='set title "'+datg+'\\n \\\n '
gggd='plot \\\n'
dsplito=set([])
for datlabel in gggtitle1:
   if len(datlabel)==0: continue
   dsplit=datlabel.split('/')
   dadd=''
   if ' '.join(dsplit[0:2]) in dsplito:
       pass
   else:
       dadd = dadd + '\\n'+'/'.join(dsplit[0:2])+'\\\n'
       #dadd = dadd + '  '+args+' : '

   args= arglines[i].split(' ')[3]#' '.join( arglines[i].split(' ')[1:4] )
   #dadd=dadd+'     '+dsplit[2]
   datl = datl+dadd+'('+str(i+1)+')'
   hhh=''
   if i>0 : hhh=', '
   gggd= gggd+ hhh+'"datg.'+datg+'"  index ' \
   + str(i) + " using (($"+xx+"+$"+zz+")*.529177):($"+yy+"+"+str(eref)+") with lp ti '(" \
   + str(i+1) + ') '+ dsplit[2] +args+"'"
   i=i+1 
   dsplito.add(' '.join(dsplit[0:2]))
   #print dsplito
datl=datl+'"\n'   
gggout = gggheader+datl+gggd+gggtail+' \n'
print gggout

