#!/usr/bin/env python
import string,sys,os,commands
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

### template file for job script ###
jobtempbg = """\
#!/bin/bash
#@$-lP 1
#@$-lp 1
#@$-lm 3.5gb
#@$-eo
#@$-oi
#cd ${QSUB_WORKDIR}
#LD_LIBRARY_PATH=/home/usr5/f70205a/local/lib:/home/usr5/f70205a/local/lib64:$LD_LIBRARY_PATH
#LD_RUN_PATH=/home/usr5/f70205a/local/lib:/home/usr5/f70205a/local/lib64:$LD_RUN_PATH
#PATH=${HOME}/bin:${HOME}/local/bin:$PATH
source atomlist.bash
source homodimerdistance.bash
source extra.bash
#source extra_nrel.bash
REPLACEHERE___
"""

jobtemppjsub = """\
#!/bin/bash
#PJM -L "node=1"
#PJM -L "rscgrp=fx-single"
#PJM -L "elapse=48:00"
#PJM -X
#PJM --no-stging
PATH=${HOME}/bin:${HOME}/local/bin:$PATH
source atomlist.bash
source homodimerdistance.bash
source extra.bash
REPLACEHERE___
"""

jobtemppjsubnrel = """\
#!/bin/bash
#PJM -L "node=1"
#PJM -L "rscgrp=fx-single"
#PJM -L "elapse=48:00"
#PJM -X
#PJM --no-stging
PATH=${HOME}/bin:${HOME}/local/bin:$PATH
source atomlist.bash
source homodimerdistance.bash
#source extra.bash
source extra_nrel.bash
REPLACEHERE___
"""

jobtempbgnrel = """\
#!/bin/bash
#@$-lP 1
#@$-lp 1
#@$-lm 3.5gb
#@$-eo
#@$-oi
#cd ${QSUB_WORKDIR}
#LD_LIBRARY_PATH=/home/usr5/f70205a/local/lib:/home/usr5/f70205a/local/lib64:$LD_LIBRARY_PATH
#LD_RUN_PATH=/home/usr5/f70205a/local/lib:/home/usr5/f70205a/local/lib64:$LD_RUN_PATH
#PATH=${HOME}/bin:${HOME}/local/bin:$PATH
source atomlist.bash
source homodimerdistance.bash
#source extra.bash
source extra_nrel.bash
REPLACEHERE___
"""

jobtempqsub = """\
#!/bin/bash
@$-lP 1
@$-lp 1
@$-lm 3.5gb
@$-eo
@$-oi
cd ${QSUB_WORKDIR}
LD_LIBRARY_PATH=/home/usr5/f70205a/local/lib:/home/usr5/f70205a/local/lib64:$LD_LIBRARY_PATH
LD_RUN_PATH=/home/usr5/f70205a/local/lib:/home/usr5/f70205a/local/lib64:$LD_RUN_PATH
PATH=${HOME}/bin:${HOME}/local/bin:$PATH
source atomlist.bash
source homodimerdistance.bash
source extra.bash
#source extra_nrel.bash
REPLACEHERE___
"""

### main from here ###    
jobtemp = """#!/bin/bash
source atomlist.bash
source homodimerdistance.bash
source extra.bash
#source extra_nrel.bash
REPLACEHERE___
"""

#open(sys.argv[2],'rt').read() #open('jobtemplate','rt').read()
print sys.argv
if '--bg' in sys.argv:
    mode='bg'
    print ##### Generate jobque files for backgound mode #####'
elif '--pjsub' in sys.argv:
    print '##### Genarate jobque files for pjsub mode #####'
    mode='pjsub'
else:
    print 'no option for jobquescript'
    sys.exit()

if '--atom' in sys.argv:
    print 'Atom mode'
    atom=1
else :
    print 'Dimer mode'
    atom=0

if '--nrel' in sys.argv:
    print 'nrel mode'
    nrel=1
    if(mode=='bg'): mode='bgnrel'
    if(mode=='pjsub'): mode='pjsubnrel'
else :
    print 'rel mode'
    nrel=0
if '--continue' in sys.argv:
    cont=1
else:
    cont=0


patt    = string.split( open(sys.argv[1],'rt').read(),'\n') 
for comm in patt:
    if len(comm) <5 : continue
    if comm[0]=='#': continue
# distance
    dat = keydata(comm, 'distance=',';')
    distance = rangereal(float(dat[0]),float(dat[1]),0.1)
#    print distance
    if(atom==1): distance=['0.0']
    fsmoms = keydata(comm, 'fsmom=',',')
    rstars = keydata(comm, 'rstar=',',')
    enddat = enddata(comm)
    print 
    print '=== ',distance,' ',fsmoms,' ',rstars,' ',enddat
    print comm
    jjj='jobque'
    if cont==1: jjj=jjj+'.continue'
    jobque     = open(jjj,'a')
    #jobquebgnrel = open('jobque.bgnrel','a')
    #jobqueqsub   = open('jobque.qsub','a')
    for dis in distance:
        for fsmom in fsmoms:
            for rstar in rstars:
                if(atom==1): command='jobatom1.py '+dis+', fsmom='+fsmom+'@ rstar='+rstar+'@ '+enddat
                if(atom==0): command='jobmoldimer1.py '+dis+', fsmom='+fsmom+'@ rstar='+rstar+'@ '+enddat

                # Generate ctrl files
                if cont==0:
                    jobfile = open('jobtempexec','w')
                    fff=string.replace(jobtemp,'REPLACEHERE___',command)
                    jobfile.write(fff)
                    jobfile.close()
                    os.system("bash jobtempexec")
                    os.system("rm jobtempexec")
                # directory name for working (where we have ctrl files).
                f=open('ctrldir','r')
                ddd=f.read()
                #print 'ctrldir= ',ddd
                #print 'cwd    = ',os.getcwd()
                f.close()

                # Check previous calculations. tail save.dimer, and findout where we restart.
                if cont==1:
                    initic=-1
                    tailsave=commands.getstatusoutput('tail -1 '+ddd+'/save.dimer')
                    if tailsave[0]==0: 
                        print 'tail savedimer=',tailsave[1]
                        pwe=tailsave[1].split('pwemax=')[1].split(' ')[0]
                        bzw=tailsave[1].split('bzw=')[1].split(' ')[0]
                        print 'pwe bzw from the end of save.dimer=',pwe,bzw
                        if(pwe==2 and bzw==.01):
                            initic=1
                        else:
                            initic=int(pwe)
                else:
                    initic=-9999

                # all ctrl files sorted.
                clines=''
                ctrlfiles=string.split(commands.getoutput('ls -1 '+ddd+'/ctrl.dimer.*'),'\n')
                for ic in sorted(ctrlfiles):
                    pwexx=string.split(ic,'ctrl.dimer.')[1]
                    cname= 'ctrl.dimer.'+pwexx
                    if int(pwexx) >= initic:
                        clines= clines+'cp '+cname+' ctrl.dimer;bash ctrl.dimer;rm mixm.dimer\n'
                #print clines

                # job file
                jobname=string.replace(string.replace(command,' ','_'),'$','')
                if cont==1: jobname=jobname+'.continue'
                print '  jobfile=',jobname
                jobfile = open(jobname,'w')
                if(mode=='bg'):    jobtempmmm=jobtempbg
                if(mode=='pjsub'): jobtempmmm=jobtemppjsub
                if(mode=='bgnrel'):jobtempmmm=jobtempbgnrel
                if(mode=='pjsubnrel'):jobtempmmm=jobtemppjsubnrel
                fff=string.replace(jobtempmmm,'REPLACEHERE___','cd '+ddd+'\n'+clines)
                jobfile.write(fff)
                jobfile.close()

                #job que file
                jobque.write("./"+jobname+"\n")

                #jobquebgnrel.write("bash "+jobname+" &\n")
                #jobqueqsub.write("qsub -f fai "+jobname+" &\n")

                #os.system("echo qsub -f fai "+jobname)
                #os.system("qsub -q fai "+jobname)
                #os.system("echo bash "+jobname)
                #os.system("bash "+jobname+" &")
                #os.system("echo qsub -f fai "+jobname)
                #os.system("qsub -q fai "+jobname)
jobque.close()
