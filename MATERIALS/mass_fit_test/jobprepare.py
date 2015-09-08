#!/usr/bin/python2
# usage: run this command and save it to 'jobrun' file
# then run 'bash jobrun'
#
syml='''\
51  0 0 0   .5 .5  .5    Gamma  L
51  0 0 0    1.  0  0    Gamma  X
51  0 0 0   .75 .75 0    Gamma  K
-888 !note -888 start Mass line. Here is a ZB case
1025  0 0 0   .5 .5  .5   1  512     0.1 0.01    Gamma  L
1025  0 0 0    1.  0  0   1  512     0.1 0.01    Gamma  X
1025  0 0 0   .75 .75 0   1  512     0.1 0.01    Gamma  K
0 !terminator  
'''
open('syml.init','w').write(syml)

import commands
exe='~/bin'
for i in commands.getoutput('ls -1|grep "\.mass"').split("\n"):
    head= i.lower().split('_')[0]
    print 'echo ----------- ',i
    print 'pushd .'
    print 'cd '+ i
    print 'cp ../syml.init syml.'+head
    print exe+'/job_band '+ head +' -np 24 -vnspin=2 -vso=1 NoGnuplot'
    print 'popd'
