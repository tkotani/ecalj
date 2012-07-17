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
if sys.argv[2]=='--bg':
    mode='bg'
    print 'jobque file for backgound mode'
elif sys.argv[2]=='--pjsub':
    print 'jobque file for pjsub mode'
    mode='pjsub'
else:
    print 'no option for jobquescript'
    sys.exit()

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
    
    jobque     = open('jobque','a')
    #jobquebgnrel = open('jobque.bgnrel','a')
    #jobqueqsub   = open('jobque.qsub','a')
    for dis in distance:
        for fsmom in fsmoms:
            for rstar in rstars:
                command='jobmoldimer1.py '+dis+', fsmom='+fsmom+'@ rstar='+rstar+'@ '+enddat

                # Generate ctrl files
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

                # all ctrl files sorted.
                clines=''
                ctrlfiles=string.split(commands.getoutput('ls -1 '+ddd+'/ctrl.dimer.*'),'\n')
                for ic in sorted(ctrlfiles):
                    cname= 'ctrl.dimer'+string.split(ic,'ctrl.dimer')[1]
                    clines= clines+'cp '+cname+' ctrl.dimer;bash ctrl.dimer\n'
                #print clines
                   
                # job file
                jobname=string.replace(string.replace(command,' ','_'),'$','')
                print '  jobfile=',jobname
                jobfile = open(jobname,'w')
                if(mode=='bg'): jobtempmmm=jobtempbg
                if(mode=='pjsub'): jobtempmmm=jobtemppjsub
                fff=string.replace(jobtempmmm,'REPLACEHERE___','cd '+ddd+'\n'+clines)
                #fff=string.replace(jobtempbgnrel,'REPLACEHERE___','cd '+ddd+'\n'+clines)
                jobfile.write(fff)
                jobfile.close()

                #job que file
                if(mode=='bg') : jobque.write("bash "+jobname+" &\n")
                if(mode=='pjsub') : jobque.write("pjsub "+jobname+" \n")

                #jobquebgnrel.write("bash "+jobname+" &\n")
                #jobqueqsub.write("qsub -f fai "+jobname+" &\n")

                #os.system("echo qsub -f fai "+jobname)
                #os.system("qsub -q fai "+jobname)
                #os.system("echo bash "+jobname)
                #os.system("bash "+jobname+" &")
                #os.system("echo qsub -f fai "+jobname)
                #os.system("qsub -q fai "+jobname)
jobque.close()
