#!/usr/bin/env python
import string,sys,os
def keydata(comm,key,sep):
    ppp=string.split(comm,key)
    ppp2=string.split(ppp[1],'@')[0]
    return string.split(ppp2,sep)

def enddata(comm):
    ppp2=string.split(comm,'@')[-1]
    return ppp2

def rangereal(ini,end,step):
    dat=[]
    d=ini
    while 1:
        if d<end+0.0001:
            aaa= '%4.2f' % d
            #print aaa
            dat.append(aaa)
            d=d+step
        else:
            break
    return dat

### main from here ###    
jobtemp = open(sys.argv[2],'rt').read() #open('jobtemplate','rt').read()
patt    = string.split( open(sys.argv[1],'rt').read(),'\n') 
for comm in patt:
    if len(comm) <5 : continue
    if comm[0]=='#': continue
# distance
    dat = keydata(comm, 'distance=',';')
    distance = rangereal(float(dat[0]),float(dat[1]),0.1)
    fsmoms = keydata(comm, 'fsmom=',',')
    rstars = keydata(comm, 'rstar=',',')
    enddat = enddata(comm)
    print 
    print '=== ',distance,' ',fsmoms,' ',rstars,' ',enddat
    print comm
    
    for dis in distance:
        for fsmom in fsmoms:
            for rstar in rstars:
                command='jobmoldimer1 '+dis+', fsmom='+fsmom+'@ rstar='+rstar+'@ '+enddat
                fff=string.replace(jobtemp,'REPLACEHERE___',command)
                jobname=string.replace(command,' ','_')
                jobname=string.replace(jobname,'$','')
                jobfile = open(jobname,'w')
                jobfile.write(fff)
                jobfile.close()
                #print 'xxxxxxxxxxxxxxxxxxxxxxx',jobname
                os.system("echo qsub -f fai "+jobname)
                os.system("qsub -q fai "+jobname)
