#!/usr/bin/env python
import string,sys,re,copy,os

def getcaller(oxx,routine,lineall):
    routinex =set([])
    for i in range(len(oxx)-1):
        line = oxx[i]
        for r in routine:
            #print line,'xxxx',i,r
            ppp= r+'->(?P<caller>\w+)'
            #print ppp
            ff=re.match(ppp,line)
            if(ff):  
                lineall=lineall+line+'\n'
                rrr= ff.group('caller')
                routinex.add(rrr)
                break
    return routinex,lineall

###############################
src = open("callcaller.dotdata",'rt').read()
oxx= string.split(src,'\n') 
argset= sys.argv[1:]

lineall=''
routine=set(argset)
routinex=set([])
l=1
while l>0:
    #print 'iiiiiiii ',routine
    routinex,lineall = getcaller(oxx,routine,lineall)
    #print 'eeeeeeee ',routinex,len(routinex)
    #print
    l=len(routinex)
    routine = copy.deepcopy(routinex)
title='_'.join([x for x in argset])
aaa='digraph '+ title + ' {\n'+ lineall+'}\n'

######################
open(title+".dotdata",'wt').write(aaa)
ccc ='dot -Tgif '+ title + ".dotdata -o " + title+ ".gif" 
os.system(ccc)
print 'Wrote '+ title+".dotdata and "+ title+".gif"


#print 'rrrrrrrrrr'
#print lineall
#print 'rrrrrrrrrr'
#print routine

sys.exit()

