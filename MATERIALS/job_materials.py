#!/usr/bin/env python
import re,sys,string,os,signal,time,commands
scfile = open("Materials.ctrls.database",'rt').read().split('\n')
#print scfile
sec={}
aaa=''
findsection=False

### find sections ###
for line in scfile:
    if findsection:
        if line[0:3]=='---' :
            sec[tag]=aaa
            aaa=''
            findsection= False
        else:    
            aaa = aaa+ line+'\n'

    if line=='' : continue
    if line[0]=='#': continue
    if line[0]=='=': continue

    m=re.match(r'^\w+:',line)
    if m: 
        findsection= True
        tag=m.group()
        #print 'find a section=', tag
        aaa=''
#print sec.keys()
print
material={}
for line in sec['DATASECTION:'].split('\n'):
    if line=='' : continue
    if line[0]=='#': continue
    if line[0]==' ': continue
    if line[0]=='\t': continue
    print line
    m=re.match(r'\w+\s',line)
    mater=m.group().strip()
    content=line[len(mater):]
    material[mater]= content.lstrip()
    #print mater,'--->',material[mater]
AllMaterials= material.keys()
AllMaterials=sorted(AllMaterials)
print
print  '=== Materials in Materials.ctrls.database ==='
aaa=''
ix=0
for x in AllMaterials:
    aaa= aaa + '  '+ x
    ix=ix+1
    if(ix==10): 
        print aaa
        ix=0
        aaa=''

###############
nargv = len(sys.argv) -1
args=sys.argv[1:]
argset= set(sys.argv[1:])
#print argset


#####################################################
print 
if not '--noexec' in  argset and '--all' in argset:
    for i in ['LaGaO3','4hSiC','Bi2Te3']:
        print ' To reduce CPUtime, --all skip ',i
        AllMaterials.remove(i)
#print AllMaterials
#sys.exit()
#AllMaterials=AllMaterials.split('\s*')
#####################################################


if ('--all' in  argset):
	choosedmaterials=AllMaterials
elif(nargv==0 or '--help' in argset):
      print
      print " PURPOSE: perform GGA/LDA calculations for sysmtems in Materials.ctrl.database"
      print " USAGE:   job_materials.py [options] "
      print "   [options] are material names above separeted by space."
      print "      --all : all materials are specified (in an hour now.)"
      print "      --noexec: not do lmf. Just generates directories and ctrl files"
      print " NOTE: to set number of cores for lmf-MPIK, please modify this code."
      print "       Search lmf-MPIK in this file"
      sys.exit(-1)
else:
    choosedmaterials=args

###########################################
print 
print ' We are going to do LDA for ',choosedmaterials
dummy=raw_input('Push return to continue!')
if dummy != '': sys.exit()

################################################
for matr in choosedmaterials:
    lll = material[matr]
    m = re.search(r'\w+:',lll)
    STRUCTURE= m.group()
    #print STRUCTURE
    aaa= re.sub(STRUCTURE,'',material[matr]).lstrip()
    matdat= re.split(r'\s+',aaa)
    #print  matdat
    #sys.exit()
    constdat=''
    option=''
    optionlmf=''
    for ii in matdat:
        if re.match(r'@[0-9]+=',ii):
            pass
        elif re.match(r'--\w*',ii):
            option=option +' ' + ii
        elif re.match(r'lmf-+\w*',ii):
            optionlmf=optionlmf +' ' + ii.lstrip('lmf')
        else:
            constdat=constdat+' '+ii
        
    aaa= '#id  = '+ matr +'\n'
    #print STRUCTURE
    #print sec
    structemp = sec[STRUCTURE]
    #print '---- -----'
    #m = re.findall(r'@[0-9]+',sec[STRUCTURE])
    #print m
    #print '*************************'
    for ix in matdat:
        try:
            m=re.match(r'@[0-9]+=',ix)
            mid= m.group()
            mat=re.split(mid,ix)
            #print 'vvvvvvvv',mid[0:-1],mat[1]
            structemp=re.sub(mid[0:-1]+' ',mat[1]+' ',structemp)
            structemp=re.sub(mid[0:-1]+'\n',mat[1]+'\n',structemp)
        except:
            pass
    stot=''    
    for iss in structemp.split('\n'):
        if re.match(r'\%\s*const',iss):
            iss= iss + constdat
        stot=stot+iss+'\n'
    aaa = aaa+  stot
    ext=string.lower(matr)
    ctrlsnm = "ctrls."+ext
    ctrlgenc= 'ctrlgenM1.py '+ext+' ' + option 

    print '==='+matr+' /  ',ctrlsnm,' '+option+ ' ==='
    print ' command=',ctrlgenc
    print aaa

    ctrls = open(ctrlsnm ,'w')
    ctrls.write(aaa)
    ctrls.close()

    os.system(ctrlgenc)
    os.system('mkdir '+matr)
    os.system('cp ctrlgenM1.ctrl.'+ext+' ctrl.'+ext)
    os.system('cp -f *.'+ext +' '+matr)
    if os.path.exists("./RSEQ_ERROR"): os.system('cp -f RSEQ* '+matr)
    os.system('rm -f *tmp* RSEQ*')
    rdir=os.path.dirname(os.path.abspath(sys.argv[0]))
    wdir= os.path.dirname(os.path.abspath(sys.argv[0]))+'/'+matr

    os.chdir(wdir)
    os.system('pwd')
    os.system('lmfa '+ext+optionlmf+' >llmfa')
    #joblmf='lmf  '+ext+optionlmf+' >llmf'   
    joblmf='mpirun -np 2 lmf-MPIK  '+ext+optionlmf+' >llmf'
    print 'See joblmf file for lmf-MPIK command arguments.'
    os.system('echo '+joblmf +'>joblmf')
    if ('--noexec' in  argset):
        pass
    else:
        try:
            commands.getoutput(joblmf)
        except KeyboardInterrupt:
            'Keyboad interrrupt by user'
            sys.exit()

    os.chdir(rdir)
    os.system('pwd')

