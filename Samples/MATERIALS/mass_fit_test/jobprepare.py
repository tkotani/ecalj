#!/usr/bin/python2
# usage: run this command and save it to 'jobrun' file
# then run 'bash jobrun'
#
syml='''\
### ndiv2, ninit2 nend2 etolv etolc are for mass mode.
#ndiv qleft(1:3) qright(1:3) llabel rlabel  ndiv2 ninit2 nend2 etolv(Ry) etolc(Ry)
51    0 0 0      .5 .5  .5   Gamma  L       513     1    81    0.1       0.01  
51    0 0 0      1.  0  0    Gamma  X       513     1    81    0.1       0.01  
51    0 0 0      .75 .75 0   Gamma  K       513     1    81    0.1       0.01  
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
    print exe+'/job_band '+ head +' -np 4 -vnspin=2 -vso=1 NoGnuplot'
    print 'popd'
