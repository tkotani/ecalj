#!/usr/bin/env python
import sys,os
thisdir= os.path.dirname(os.path.abspath(__file__))
argset= sys.argv[1:]
for ffile in argset:
    print('=========== '+ ffile + ' start ===========')
    os.system('jobindent '+ ffile)
    os.system('mv '+ ffile+' OLD')
