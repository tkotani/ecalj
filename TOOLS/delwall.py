#!/usr/bin/env python
# Program unit analyzer.
import sys,os
thisdir= os.path.dirname(os.path.abspath(__file__))
#pathname=thisdir+ ' '+thisdir+'/delw0419'
delw=thisdir+"/KINO/del_w2.3/delw.py"
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
for i in argset:
    print '------- file name:', i, '-------------'
    os.system("python "+delw+" < "+ i + " > ./converted/"+i)
    os.system("cp file.map ./converted/"+i+'.file.map')
sys.exit()
