#!/usr/bin/env python
import os, sys, string, re

#---------------------------------------------------
def lineReadfile(filename):
	""" read file and make list of the content of the file, and \n->'' """
#	"input:filename output=readlines() "
	f = open(filename)
	list1 =[]
	while 1:
		s = f.readline()
		if s=="":
			break
		s=string.replace(s,"\n","")
		if s=="":
			continue
		list1.append(s)
	f.close()
	return list1

### help sections ###
if(len(sys.argv) -1 !=2):
    print \
""" 
rescalectrl.py
---------------
 Purpose : 
   Rescale numbers in ctrls file.
 Example (how to use): 
   >rescalectrl.py 1.88972687777*8.085 ctrls.nd2fe14b
   the first argument is a number in python format.
"""
    sys.exit()
alatinput= sys.argv[1]
fname    = sys.argv[2]
print "# ... convert scale of ctrls file... Keep the same structure"
fff = lineReadfile(fname)
outfile=''

platon=False
nplat=0
for line in fff:
    linen=line
#    print
#    print 'LINEIN='+line

## ALAT token
    if( 'ALAT' in line): 
        alatold = float(line.split('ALAT')[1].split('=')[1])
        alatnew = eval(alatinput)
        ratio=alatnew/alatold
        linen= re.sub(r'ALAT\s*=\s*[\.\d]*','ALAT='+alatinput,line)+' '
        print '# conversion by resclalectrl.py'
        print '#  alatold=',alatold, ' ==> alatnew=',alatnew
        print '#  ratio=alatnew/alatold=',str(ratio)

## PLAT token
    if( 'PLAT' in line): 
        platon=True
    if platon and nplat <9 : #get nine numbers after PLAT=
        mplat = re.search(r'PLAT\s*=\s*',line)
        if mplat:
            lll='      PLAT= '
            lsp=re.split(r'PLAT\s*=\s*',line)[1]
        if not mplat:  #if data of PLAT is in next lines
            lll='            '
            lsp=line
#        print lsp.split()
#        pdat=''
#        pdat= (pdat+str(float(x)/ratio) for x in lsp.split())
#        print 'ppp ',pdat  #(float(x)/ratio for x in nnn)
        for nnn in lsp.split(): # lsp.split() contains numbers
            pdat = "%20.15f" % (float(nnn)/ratio)
            nplat= nplat+1
            lll=lll+ str(pdat)+' '
        linen=lll

### ATOM POS=
    if 'ATOM' in line and 'POS' in line:
        linep= re.search(r'POS\s*=\s*[\.\d]* \s*[\.\d]* \s*[\.\d]*',line)
        numline=' '
        for aaa in linep.group(0).split():
            aaanum= re.search(r'[\d\.]+',aaa)
            if aaanum : numline=numline +' '+ "%20.15f" % (float(aaanum.group(0))/ratio )
        linen= re.sub(r'POS\s*=\s*[\.\d]* \s*[\.\d]* \s*[\.\d]*','POS='+numline,line)
#       print 'posout=', linen
    
### output lines are connected 
    outfile = outfile + linen +'\n'

#    print 'LINOUT:'+linen
#print '-------outfile-----------------'
print outfile

#f=open(fname,'rt')
#fin = f.read() 
#f.close()
