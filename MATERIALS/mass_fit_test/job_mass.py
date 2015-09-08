#!/usr/bin/python2
import os,sys
import commands
ggg1="""
### massfit by gnuplot ####################
outpng=ttt.'.png'
min(x,y) = (x < y) ? x : y
max(x,y) = (x > y) ? x : y

set terminal png size 1280,960
set out outpng
set grid
set ylabel " eV "
set xlabel " q/(2pi alat) (1/a.u)"

ry=13.605

# f1(x) = x**2*(a1 + x**2*b1 ) 
# fit [x=efitmin:efitmax] f1(x) 'data1' u ($5):($7) via a1,b1
# qmass1=ry/a1
# qe0=ry**2/(qmass1**2*b1)
# f2(x) = x**2*(a2 + x**2*b2 )
# fit [x=efitmin:efitmax] f2(x) 'data2' u ($5):($7) via a2,b2
# qmass2=ry/a2
# qe02=ry**2/(qmass2**2*b2)


#f2(x) = b2*(abs(x)*(1+c2*abs(x)))
#fit f2(x) 'data2' u ($7):($5**2) via b2,c2
#f1(x) = a*( sqrt(1.0 + abs(bb)*x**2) -1.0)
#fit f1(x) 'data1' u ($5):($7) via a,bb
#b=a*abs(bb)
#ddf1(x)=x
# f2(x) = a2*( sqrt(1.0 + abs(bb2)*x**2) -1.0)
# fit f2(x) 'data2' u ($5):($7) via a2,bb2
# b2=b*abs(bb2)
# ddf2(x)=x

qq1(x) = mass1ry*abs(x) + mass1ryE01m*x**2
fit [x=efitmin:efitmax] qq1(x) 'data1' u ($7):($5**2) via mass1ry,mass1ryE01m
qq2(x) = mass2ry*abs(x) + mass2ryE02m*x**2
fit [x=efitmin:efitmax] qq2(x) 'data2' u ($7):($5**2) via mass2ry,mass2ryE02m

E01m=mass1ryE01m/mass1ry
E02m=mass2ryE02m/mass2ry
E01=1/E01m
E02=1/E02m

ee1(x)=E01/2*(sqrt(1+4*x**2/(mass1ry*E01))-1)
ee2(x)=E02/2*(sqrt(1+4*x**2/(mass2ry*E02))-1)


"""

ggg2="""

#mass= ry/((b+b2)/2)
#massinv=1.0/mass

# dd10= 2.*ry/ddf1(0.)
# dd20= 2.*ry/ddf2(0.)
# dd1a= 2.*ry/ddf1(0.02)
# dd2a= 2.*ry/ddf2(0.02)
# dd1b= 2.*ry/ddf1(0.04)
# dd2b= 2.*ry/ddf2(0.04)

#tt= ttt.sprintf("=  %6.3f %6.3f %6.3f ; %6.3f %6.3f %6.3f",dd10,dd1a,dd1b,dd20,dd2a,dd2b) 

mxx=sprintf(" mass= %6.3f %6.3f; E0(eV)= %6.2f %6.2f", mass1ry*ry,mass2ry*ry,E01,E02)
formula='   q**2/(2*mass)=(E*(1+E/E0))*ry/2, where E in eV ry=13.605'
set title ttt.mxx.formula

set xrange[xl:1.5*xr] 
set yrange[0:0.08]

#plot \
#'data1' u ($5):($7) lt 1 pt  1 w p,\
#'data2' u ($5):($7) lt 2 pt  2 w p,\
#f1(x),f2(x),b*x**2

plot \
'data1' u ($5):(abs($7)) lt 1 pt  1 w p,'data2' u ($5):(abs($7)) lt 2 pt  2 w p,\
ee1(x),ee2(x)


mxxf=sprintf(" mass,E0(eV)= %6.3f %6.2f ", (mass1ry*ry+mass2ry*ry)/2,(E01+E02)/2)

#print mxx
print 'mmm3: ',ttt.mxx
print 'mmm4: ',ttt.mxxf
print 'mmm5: ',outpng
print 'mmm6: '
exit
"""

def cutew(fin,cute,fout,qupper):
            finread=fin.read().split('\n')
            fouth = open(fout,'w')
            i=0
            qmin=1e10
            qmax=-1e10
            for iline in finread:
                i=i+1
                if(i<3): continue
                if(len(iline)==0): break
                ilines=iline.split()
#                print i,ilines[6]
                #print ilines[6]
#                print iline
                if(float(ilines[4])>qupper): continue
                if(abs(float(ilines[6]))<cute): 
                    if(float(qmin)>float(ilines[4])): 
                        qmin=float(ilines[4])
                    if(float(qmax)<float(ilines[4])): 
                        qmax=float(ilines[4])
                    ilines[6]=str(abs(float(ilines[6]))) #energy positive
                    fouth.write(' '.join(ilines)+'\n')
#                    print iline
            fouth.close()
            return qmin,qmax

def eset(aaa,bbb):
    return "-e '"+aaa+"="+bbb+"' "
print "=== Band00A (A,B)=(5,6) for mhh, (A,B)=(3,4) for mlh, (A,B)=(1,2) for mso ZB case ==="
print " Syml00C C=4 for [111], C=5 for [100], C=6 for [110] "

### get directory which ends .mass
matlist=[]
for i in os.listdir('.'):
    if os.path.isdir(i) and '.mass' in i:
        matlist=matlist+ [i]

#matlist= commands.getoutput('ls -1 *|grep mass').split("\n")
matlistxxx="""
AlAs_so.qsgw.mass
AlSb_so.qsgw.mass
CdSe_so.qsgw.mass
CdSe_so.lda.mass
CdS_so.qsgw.mass
CdTe_so.qsgw.mass
GaAs_so.qsgw.mass
GaSb_so.qsgw.mass 
GaP_so.qsgw.mass
InP_so.qsgw.mass
InAs_so.qsgw.mass
InSb_so.qsgw.mass
InSb_so.qsgw.scaled.mass
MgS_so.qsgw.mass
MgS_so.qsgw.mass
ZnSe_so.qsgw.mass
ZnS_so.qsgw.mass
ZnTe_so.qsgw.mass
"""
#for mat in matlist.split():

for mat in matlist:
####
    print mat
    if '3csic' in mat.lower(): continue

    print 'mmm0: '+mat
    for di  in ['111', '100', '110']:
#    for di  in ['100']:
        if di=='111': syml='4'
        if di=='100': syml='5'
        if di=='110': syml='6'
        for mxx  in ['mso', 'mlh', 'mhh', 'mee']:
            yrange=''
            xrange=''
            if mxx=='mso':
                band1='1'
                band2='2'
            if mxx=='mlh':
                band1='3'
                band2='4'
            if mxx=='mhh':
                band1='5'
                band2='6'
            if mxx=='mee':
                band1='7'
                band2='8'
####### fitting region #############################
            efitmin='0.01'
            efitmax='0.05' 
            qupper=10000
########## special cases ###########
            if 'AlSb' in mat and mxx=='mee':
                qupper=0.05 #upper limit
            if 'HgSe' in mat and mxx=='mlh':
                qupper=0.02 #upper limit
                efitmax='0.03' 

            fname=mat+'/Band00'+band1+'Syml00'+syml+'Spin1.mass'
            print 'mmm1:',fname
            f1=open(fname,'r')
            qmin,qmax=cutew(f1,2*float(efitmax),'data1',qupper)

            f1.close()
            fname=mat+'/Band00'+band2+'Syml00'+syml+'Spin1.mass'
            print 'mmm2:',fname
            f2=open(fname,'r')
            qmin,qmax=cutew(f2,2*float(efitmax),'data2',qupper)
#            print qmin,qmax


            xl=str(qmin)
            xr=str(qmax)

#            xl='0.0'
#            xr='0.01'
#            if mxx=='mhh':
#                xl='0.0'
#                xr='0.1'

#             if 'GaNzb' in mat and mxx=='mhh':
#                 xl='0.0'
#                 xr='0.04'
#             if 'GaAs' in mat and mxx=='mhh':
#                 xl='0.0' 
#                 xr='0.04'
#             if 'GaSb' in mat and mxx=='mhh':
#                 xl='0.0' 
#                 xr='0.04'

#            print mat,di,mxx
#            if 'Cd' in mat and mxx=='mhh':
#                 xl='0.0' 
#                 xr='0.1'

            
#           f2.close()
#           sys.exit()

            f=open('JOBmassfit.glt','w')
            f.write(ggg1+xrange+yrange+ggg2)
            f.close()
              
#            data1= '"'+mat+'/Band00'+band1+'Syml00'+syml+'Spin1.mass"'
#            data2= '"'+mat+'/Band00'+band2+'Syml00'+syml+'Spin1.mass"'
#
	    ttt  = '"'+mat+'_'+mxx+di+'"'
#            aaa='gnuplot '+ eset('data1',data1)+eset('data2',data2)+eset('ttt',ttt) \
#                +eset('xl',xl)+ eset('xr',xr) +' JOBmassfit.glt'
            aaa='gnuplot '+eset('ttt',ttt) +eset('efitmin',efitmin)+eset('efitmax',efitmax) \
                +eset('xl',xl)+ eset('xr',xr) +' JOBmassfit.glt'
            out=commands.getoutput(aaa+' 2>&1 |grep mmm')
            print aaa
            print out
