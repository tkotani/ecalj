#!/usr/bin/env python
import re,sys,string,os,signal,time
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
if ('--all' in  argset):
	choosedmaterials=AllMaterials
elif(nargv==0 or '--help' in argset):
      print
      print " PURPOSE: perform GGA/LDA calculations for sysmtems in Materials.ctrl.database"
      print " USAGE:   job_materials.py [options] "
      print "   [options] are material names above separeted by space"
      print "   or --all for all calculations. (in an hour now.)"
      print
      sys.exit(-1)
else:
    choosedmaterials=args
################
print 
print ' !!! Go Materials for',choosedmaterials
print 

for matr in choosedmaterials:
    lll = material[matr]
    m = re.search(r'\w+:',lll)
    STRUCTURE= m.group()
    aaa= re.sub(STRUCTURE,'',material[matr]).lstrip()
    matdat= re.split(r'\s+',aaa)
    #print  matdat
    constdat=''
    option=''
    for ii in matdat:
        if re.match(r'@[0-9]+=',ii):
            pass
        elif re.match(r'--\w*',ii):
            option=option +' ' + ii
        else:
            constdat=constdat+' '+ii
        
    aaa= '#id  = '+ matr +'\n' 
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
    os.system('mv *.'+ext +' '+matr)
    os.system('rm *tmp*')
    rdir=os.path.dirname(os.path.abspath(sys.argv[0]))
    wdir= os.path.dirname(os.path.abspath(sys.argv[0]))+'/'+matr

    os.chdir(wdir)
    os.system('pwd')
    os.system('lmfa '+ext+' >llmfa')
    os.system('lmf  '+ext+' >llmf')
    os.chdir(rdir)
    os.system('pwd')

 
