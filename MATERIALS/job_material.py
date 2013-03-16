#!/usr/bin/env python
import re,sys,string,os
scfile = open("Materials.ctrls",'rt').read().split('\n')
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
        print 'find a section=', tag
        aaa=''
print sec.keys()

print '------ Material list --------'
material={}
for line in sec['DATASECTION:'].split('\n'):
    if line=='' : continue
    if line[0]=='#': continue
    m=re.match(r'\w+\s',line)
    mater=m.group().strip()
    content=line[len(mater):]
    material[mater]= content.lstrip()
    #print mater,'--->',material[mater]
print material.keys()

mmmkeys= material.keys()

lll=''
for matr in ['Cu']:#mmmkeys:
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
    ctrlgenc= './ctrlgenM1.py '+ext+' ' + option 

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
    os.system('rm *.tmp*')
    rdir=os.path.dirname(os.path.abspath(sys.argv[0]))
    wdir= os.path.dirname(os.path.abspath(sys.argv[0]))+'/'+matr

    os.chdir(wdir)
    os.system('pwd')
    os.system('lmfa '+ext+' >llmfa')
    os.system('lmf  '+ext+' >llmf')
    os.chdir(rdir)
    os.system('pwd')
