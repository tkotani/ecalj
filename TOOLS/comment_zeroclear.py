#!/usr/bin/env python
# comment out intiallization 
import sys,os,string,re
thisdir= os.path.dirname(os.path.abspath(__file__))
#pathname=thisdir+ ' '+thisdir+'/delw0419'
delw=thisdir+"/KINO/del_w2.2/delw.py"
#print pathname
#sys.path.append(pathname)

nargv = len(sys.argv) -1
if(nargv==0 or sys.argv[1]=='--help'):
    print ' === Convert W(oxxx) to regular f90 allocation ==='
    print ' usage: deawall.py  file1.F file2.F ... '
    sys.exit()
argset= sys.argv[1:]
#os.system("ls")
if(not os.path.exists('converted')): os.mkdir('./converted')

p1 = "^\s.*if\s*\(.*\s*\)\s*\w*_w_o\w"
p2 = "^\s.*\Ww\W"
p3 = "^\s.*common.*\/.*w.*\/"

p4 = "^\s.*if\s*\(\s*-.*\s*\)\s*\w*_w_o\w"

#targetf = "("+p1+"|"+p2+"|"+p3+")"
#targetf = "("+p1+"|"+p3+")"
targetf = "("+p1+")"

#targetf = "("+p4+")"

pf =re.compile(targetf,re.I)
pf2=re.compile("allocated",re.I)

pf4=re.compile(p4,re.I)

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
        a=pf.search(line)
        b=pf2.search(line)
        c=pf4.search(line)
        if(a and (not b) and (not c)):
            ifind=1
            zzz= zzz + '%d' % i + line +'\n'
            linenew='Ctakao_ZeroClear_NotRequiered '+ line
        else:
            linenew=line    
        output=output+linenew+'\n'
    if(ifind==1): 
        print zzz
        out = open('./converted/'+fn,'w')
        out.write(output)

sys.exit()
