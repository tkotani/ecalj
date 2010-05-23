#!/usr/bin/env python
# A general program templete for find and replace for files.
#
import sys,os,string,re
thisdir= os.path.dirname(os.path.abspath(__file__))
#pathname=thisdir+ ' '+thisdir+'/delw0419'
#delw=thisdir+"/KINO/del_w2.2/delw.py"
#print pathname
#sys.path.append(pathname)

#pf2=re.compile("allocated",re.I)
#pf4=re.compile(p4,re.I)
#targetf = "("+p1+"|"+p2+"|"+p3+")"
#targetf = "("+p1+"|"+p3+")"
#targetf = "("+p4+")"
#targetf = "("+p1+A")"

p1 = "\w+_w_\w+"
pf1 =re.compile(p1,re.I)

def rep_kinos(line):
    itry= True
    ix=0
    while itry:
        a=pf1.search(line)
        #b=pf2.search(line)
        #c=pf4.search(line)
        if(a): # and (not b) and (not c)):
            #ifind=1
            #zzz= zzz + '%d' % i + line +'\n'
            #linenew = line
            word=a.group()
            word2=re.split('_w_o',word)
            wx=word2[1]+'_'+word2[0]
            line=re.sub(a.group(),wx,line,1)
            ix=1
        else:
            itry=False
    #if(ix==1): print 'xxxxxx: ',line
    return line



nargv = len(sys.argv) -1
#if(nargv==0 or sys.argv[1]=='--help'):
#    print ' === Convert W(oxxx) to regular f90 allocation ==='
#    print ' usage: deawall.py  file1.F file2.F ... '
#    sys.exit()
argset= sys.argv[1:]

#os.system("ls")
if(not os.path.exists('converted')): os.mkdir('./converted')








##############################################3
p1 = "^\s.*if\s*\(.*\s*\)\s*\w*_w_o\w"
p2 = "^\s.*\Ww\W"
p3 = "^\s.*common.*\/.*w.*\/"
p4 = "^\s.*if\s*\(\s*-.*\s*\)\s*\w*_w_o\w"
p5 = "^\s.*\Ww\s*\(\s*"
#targetf = "("+p1+"|"+p2+"|"+p3+")"
#targetf = "("+p1+"|"+p3+")"
#targetf = "("+p4+")"
#targetf = "("+p1+")"
p2 = "^Cdelw"
pf2=re.compile(p2,re.I)
p3 = "^Cgetarg"
pf3=re.compile(p3,re.I)
for fn in argset:
#    print '------- file name:', i, '-------------'
    zzz= '------- file name:'+ fn + '-------------\n'
    src = open(fn,'rt').read()
    oxx= string.split(src,'\n') 
    i=0
    ifind=0
    output=''
    for ix in range( len(oxx)-1):
        i=i+1
        line=oxx[ix]
        a=pf2.search(line)
        b=pf3.search(line)
#        c=pf4.search(line)
        if(a or b):
            ifind=1
            #zzz= zzz + '%d' % i + line +'\n'
            #linenew='Ctakao_ZeroClear_NotRequiered '+ line
        else:
            pass
            linenew=line    
            output=output+linenew+'\n'
###
    oxx= string.split(output,'\n') 
    i=0
    ifind=0
    output=''
    for ix in range( len(oxx)-1):
        i=i+1
        line=oxx[ix]
        output = output+ rep_kinos(line)+'\n'
    out = open('./converted/'+fn,'w')
    out.write(output)
    out.close()


sys.exit()











#### A case to comment out zeroclear in Kino's
#--- patterns ---
for fn in argset:
    zzz= '------- file name:'+ fn + '-------------\n'
    src = open(fn,'rt').read()
    oxx= string.split(src,'\n') 
    i=0
    ifind=0
    output=''
    for ix in range( len(oxx)-1):
        i=i+1
        line=oxx[ix]
        output = output+ rep_kinos(line)+'\n'
#    print output
#     out = open('./converted/'+fn,'w')
#     out.write(output)



##############################################3
p1 = "^\s.*if\s*\(.*\s*\)\s*\w*_w_o\w"
p2 = "^\s.*\Ww\W"
p3 = "^\s.*common.*\/.*w.*\/"
p4 = "^\s.*if\s*\(\s*-.*\s*\)\s*\w*_w_o\w"
p5 = "^\s.*\Ww\s*\(\s*"
#targetf = "("+p1+"|"+p2+"|"+p3+")"
#targetf = "("+p1+"|"+p3+")"
#targetf = "("+p4+")"
#targetf = "("+p1+")"
p2 = "^Cdelw"
pf2=re.compile(p2,re.I)
p3 = "^Cgetarg"
pf3=re.compile(p3,re.I)
for fn in argset:
#    print '------- file name:', i, '-------------'
    zzz= '------- file name:'+ fn + '-------------\n'
    src = open(fn,'rt').read()
    oxx= string.split(src,'\n') 
    i=0
    ifind=0
    output=''
    for ix in range( len(oxx)-1):
        i=i+1
        line=oxx[ix]
        a=pf2.search(line)
        b=pf3.search(line)
#        c=pf4.search(line)
        if(a or b):
            ifind=1
            #zzz= zzz + '%d' % i + line +'\n'
            #linenew='Ctakao_ZeroClear_NotRequiered '+ line
        else:
            pass
            linenew=line    
            output=output+linenew+'\n'
    #if(ifind==1): 
        #print zzz
    out = open('./converted/'+fn,'w')
    out.write(output)
sys.exit()
