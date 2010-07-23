#!/usr/bin/env python
# python 2.4.3
# This routine checks module-dependency in fortran90 and compile them in right order.
#
import os
import sys
import string
import re
import glob
import shutil
import time
import resource
import os.path
import filecmp

def cpu_time_child():
	return resource.getrusage(resource.RUSAGE_CHILDREN)[0]


def testrun(testname, commanddir,datadir,workdir,commands,start,testc,enforce,deletetemp):
	if(enforce==1): shutil.rmtree(workdir,'ignore_errors')

# set up work directory and copy start files into it.
	try:
		os.mkdir(workdir)
	except:
		print workdir +' already exist. You can use --enfoece, which clear the directory first.'
		sys.exit(-1)
	print '--> start run for '+ dir + ' ==='
	print ' commanddir =',commanddir
	print ' datadir    =',datadir
	print ' workdir    =',workdir
	for ff in string.split(start):
		print '   cp '+ ff + ' to workdir'
		shutil.copy(datadir+ff,workdir)

# Run commands successively
 	os.chdir(workdir)
	print "--> We will perform commands: "
	for cm in commands:
		print '   '+commanddir+cm
	print "--> now we statrt them successively."
	for cm in commands:
		print '   Runnning '+commanddir+ cm
		cpu_init= cpu_time_child()
		ret=os.system(commanddir+cm)
		if(ret==0):
			print '     ' + cm + ' --> OK! '+' UsedCPUtime=',cpu_time_child()-cpu_init
		else:
			print 'Error exit! Failed for "'+cm +'" in ' + workdir
			print '  (or killed by yourself)'
			sys.exit(-1)
# Compare results.
	print
	print '--> Compare your results with previous results'
	rett=0
	for test in testc:
		print '   Compare results by ' + commanddir+ test
		ret=os.system(commanddir+test)
		rett=rett+abs(ret)
		if(ret==0):
			print "   Compare --> OK!"
			print
			pass
		else:
			print '=== Failed for test ' + testname,' ==='
			sys.exit(-1)
			
	if(rett==0 & deletetemp==1):
		shutil.rmtree(workdir,'ignore_errors')
		print '--> we remove workdir=', workdir
	else:
		print '--> we did calculations at workdir=', workdir
	print
	print


##### main routine ##########################
testdir = os.getcwd()
ecaldir = os.path.dirname(os.path.dirname(testdir)) #+'/ecal'
exe= ecaldir+'/fpgw/exec/'
print 'exe=',exe
print 'testdir=',testdir
print 'ecal directory=',ecaldir
print ' --- Is ecal diretory OK? (return or name of directory. or edit testgw.py '
#dat=raw_input(' --- Is ecal diretory OK? (return or name of directory. or edit testgw.py) --- ')
#if(dat==''):
#	pass
#else:
#	ecaldir=dat
#	print 'ecal=',ecaldir

### Definition of tests
testname=[]
commands=[]
startfile=[]
testcommand=[]


###########
comparekey=" 'fp pot' 'fp evl'"



#
testn='gas_eps_lmfh'
testname.append(testn)
startfile.append('ctrl.gas GWinput')
commands.append(['lmfa gas > llmfa','lmf gas>llmf','eps_lmfh gas'])
datadir  = testdir+ '/' + testn+'/'
ef1='EPS0001.nlfc.dat '
ef3='EPS0003.nlfc.dat '
ef4='EPS0004.nlfc.dat '
comp1= 'diffnum '+ef1+ datadir+ef1+comparekey
comp3= 'diffnum '+ef3+ datadir+ef3+comparekey
comp4= 'diffnum '+ef4+ datadir+ef4+comparekey
ef1a='EPS0001.dat '
ef3a='EPS0003.dat '
ef4a='EPS0004.dat '
comp1a= 'diffnum '+ef1a+ datadir+ef1a+comparekey
comp3a= 'diffnum '+ef3a+ datadir+ef3a+comparekey
comp4a= 'diffnum '+ef4a+ datadir+ef4a+comparekey
testcommand.append([comp1,comp3,comp4,comp1a,comp3a,comp4a ])

#
testn='gas_epsPP_lmfh'
testname.append(testn)
startfile.append('ctrl.gas GWinput')
commands.append(['lmfa gas > llmfa','lmf gas>llmf','epsPP_lmfh gas'])
datadir  = testdir+ '/' + testn+'/'
ef1='EPS0001.nlfc.dat '
ef3='EPS0003.nlfc.dat '
ef5='EPS0005.nlfc.dat '
comp1= 'diffnum '+ef1+ datadir+ef1+comparekey
comp3= 'diffnum '+ef3+ datadir+ef3+comparekey
comp5= 'diffnum '+ef5+ datadir+ef5+comparekey
testcommand.append([comp1,comp3,comp5 ])

#
testn='fe_epsPP_lmfh_chipm'
testname.append(testn)
startfile.append('ctrl.fe GWinput')
commands.append(['lmfa fe > llmfa','lmf fe>llmf','epsPP_lmfh_chipm fe'])
datadir  = testdir+ '/' + testn+'/'
ef1='ChiPM0001.nlfc.mat '
ef2='ChiPM0002.nlfc.mat '
ef3='ChiPM0003.nlfc.mat '
ef4='ChiPM0004.nlfc.mat '
ef5='ChiPM0005.nlfc.mat '
comp1= 'diffnum '+ef1+ datadir+ef1+comparekey
comp2= 'diffnum '+ef2+ datadir+ef2+comparekey
comp3= 'diffnum '+ef3+ datadir+ef3+comparekey
comp4= 'diffnum '+ef4+ datadir+ef4+comparekey
comp5= 'diffnum '+ef5+ datadir+ef5+comparekey
testcommand.append([comp1,comp2,comp3,comp4,comp5 ])

#
testn='si_gw_lmfh'
testname.append(testn)
startfile.append('ctrl.si GWinput')
commands.append(['lmfa si > llmfa','lmf si > llmf_lda','gw_lmfh si'])
datadir  = testdir+ '/' + testn+'/'
testcommand.append(['dqpu QPU '+ datadir+ 'QPU'])

#
testn='gas_pw_gw_lmfh'
testname.append(testn)
startfile.append('ctrl.gas GWinput')
commands.append(['lmfa gas > llmfa','lmf gas > llmf_lda','gw_lmfh gas'])
datadir  = testdir+ '/' + testn+'/'
testcommand.append(['dqpu QPU '+ datadir+ 'QPU'])

#
testn='si_gwsc'
testname.append(testn)
startfile.append('ctrl.si GWinput')
commands.append(['lmfa si > llmfa','gwsc1shot si'])
datadir  = testdir+ '/' + testn+'/'
testcommand.append(['dqpu QPU '+datadir+'QPU','diffnum log.si '+datadir+'log.si'+comparekey])
#testcommand.append(['dqpu QPU '+datadir+'QPU','diff log.si '+datadir+'log.si'])

#
testn='gas_gwsc'
testname.append(testn)
startfile.append('ctrl.gas GWinput')
commands.append(['lmfa gas > llmfa','gwsc1shot gas'])
datadir  = testdir+ '/' + testn+'/'
#testcommand.append(['dqpu QPU '+datadir+'QPU','diff log.gas '+datadir+'log.gas'])
testcommand.append(['dqpu QPU '+datadir+'QPU','diffnum log.gas '+datadir+'log.gas'+comparekey])


# This is still too simplified --> too bad answer. But it is a test for instalation for NSPIN=2.
#
testn='nio_gwsc'
testname.append(testn)
startfile.append('ctrl.nio GWinput')
commands.append(['lmfa nio > llmfa','gwsc1shot nio'])
datadir  = testdir+ '/' + testn+'/'
testcommand.append(['dqpu QPU '+datadir+'QPU','diffnum log.nio '+datadir+'log.nio'+comparekey])
#testcommand.append(['dqpu QPU '+datadir+'QPU','diff log.nio '+datadir+'log.nio'])


### Readin flags ####################################
nargv = len(sys.argv) -1
argset= set(sys.argv[1:])
if (nargv ==0 or '--help' in argset):
	print ' --- Install test for GW ---'
	print '   usage :'
	print '     testgw.py [options] testname testname ... '
	print '   options:'
	print '     --help   :  this help'
	print '     --enforce:  remove temp_* directories at first even if they exists'
	print '     --all    :  test all cases (no testname is needed)'
	print '     --show   :  show all test name defined in this routine'
	print '     --deletetemp: delete temp_* directory if tests are succeeded'
	print '   testname :  testname, e.g. si_gw_lmfh as shown below'
	for i in range(len(testname)):
		print '            :', testname[i]
	sys.exit()
if ('--enforce' in  argset):
	enforce=1
else:
	enforce=0
if ('--all' in  argset):
	all=1
else:
	all=0
if ('--deletetemp' in  argset):
	deletetemp=1
else:
	deletetemp=0

### Initialization ##################################
print '### TEST run ###'
print '  ecaldir =',ecaldir
commanddir = ecaldir+'/fpgw/exec/'
print '  initialize: copy lm7K/lmfa lmf lmfgw lmfgw fpgw/exec/'
if(not ecaldir+'/lm7K/lmfa',commanddir+'lmfa'):
	shutil.copy(ecaldir+'/lm7K/lmfa',  commanddir)
if(not ecaldir+'/lm7K/lmf',commanddir+'lmf'):
	shutil.copy(ecaldir+'/lm7K/lmf',  commanddir)
if(not ecaldir+'/lm7K/lmfgw',commanddir+'lmfgw'):
	shutil.copy(ecaldir+'/lm7K/lmfgw',  commanddir)
if(not ecaldir+'/lm7K/lmf2gw',commanddir+'lmf2gw'):
	shutil.copy(ecaldir+'/lm7K/lmf2gw',  commanddir)
if(not ecaldir+'/TOOLS/diffnum',commanddir+'diffnum'):
	shutil.copy(ecaldir+'/TOOLS/diffnum',  commanddir)
print


### Run tests
for i in range(len(testname)):
	if all==0:
		if testname[i] in  argset:
			pass
		else:
			continue
			
	print '=== '+ testname[i]+' ========================================'
	dir        = testname[i]+'/'
	datadir    = testdir+ '/' + dir
	workdir    = testdir+ '/temp_' + dir
	testrun(testname[i],commanddir,datadir,workdir,commands[i],startfile[i],testcommand[i],enforce,deletetemp)

sys.exit()


