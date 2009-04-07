#!/usr/bin/python2
# This routine checks module-dependency in fortran90 and compile them in right order.
#
import os
import sys
import string
import re
import glob

#---------------------
def rsp(dat):
	dat= string.lower(re.sub(' '  ,'',dat))
	return dat

def connect(alist):
	result=''
	for lll in alist:
		result = result + ' ' + lll
	return result

def uniqlist(input): # uniq
	output = []
	for i in input:
		if not i in output:
			output.append(i)
	return output

def skipc(flines,fin):  #  remove comment, treat continuous lines,
	            #  and make list [number 
	re_nocomm = re.compile('^\ *subroutine', re.IGNORECASE)
	re_char   = re.compile('("[^"]*")|(\'[^\']*\')')
	il=0
	fff =[]
	init=1
	for fl in flines:
		il=il+1
		if(not re.search('( |#)',fl[0:1])): continue  # skip comment line
		flx = re.split('!',fl)[0]
		if(re.search('#',fl[0:1]) or re.search(' ',fl[5:6]) ):
			if(init==0):
				flo =string.lower(flo)
				flo = re_char.sub('@@@',flo) #
				flo  = re.sub('^ *','',flo) 
				flo  = re.sub(' +',' ',flo) 
				flo  = re.sub('end *do','enddo',flo, re.IGNORECASE) 
				flo  = re.sub('end *if','endif',flo, re.IGNORECASE) 
				flo  = re.sub('end *interface','endinterface',flo, re.IGNORECASE) 
				fff.append([fin,ilo,flo,floo])
			else:
				init=0
			flo =  flx    #new sentence start
			floo=  fl +'\n'   #new sentence start
			ilo= il
		else:   #continuous line            
			flo= flo + flx[6:]
			floo= floo + fl + '\n'
	flo  = re.sub('^ *','',flo) 
	flo  = re.sub(' +',' ',flo) 
	fff.append([fin,ilo,string.lower(flo),floo])
#	for x in fff:
#		print x[0:3]
#	sys.exit()
	return fff
	    
re_mod      =  re.compile( '^\ *module(?!( *procedure))', re.IGNORECASE)
rcall       =  re.compile( 'call', re.IGNORECASE)
re_function =  re.compile( \
	'((^\ *(complex|real|double|character|integer|logical)[^!]* *function))', re.IGNORECASE)
re_function2 =  re.compile( \
	'((^\ *function))', re.IGNORECASE)
re_end       = re.compile( '((^\ *end( |\Z)))', re.IGNORECASE)
re_end_mod   = re.compile( '((^\ *end *module))', re.IGNORECASE)
re_contains  = re.compile( '((^\ *contains))', re.IGNORECASE)
re_subroutine= re.compile( '^\ *subroutine', re.IGNORECASE)
re_blockdata = re.compile( '^\ *blockdata', re.IGNORECASE)
re_call      = re.compile( '^\ *call', re.IGNORECASE)
re_use       = re.compile( '^\ *use', re.IGNORECASE)
re_program   = re.compile( '^\ *program', re.IGNORECASE)
re_entry     = re.compile( '^\ *entry', re.IGNORECASE)
reff         = re.compile( 'function', re.IGNORECASE)
rcal         = re.compile( 'call', re.IGNORECASE)
ruse         = re.compile( 'use *', re.IGNORECASE)
rsub         = re.compile( 'subroutine', re.IGNORECASE)
rmod         = re.compile( 'module', re.IGNORECASE)
rprog        = re.compile( 'program', re.IGNORECASE)
rent         = re.compile( 'entry', re.IGNORECASE)
rblock       = re.compile( 'blockdata', re.IGNORECASE)



#tmp ='chekckmodule.tmp'
#tmp = 'tmp'
#print '######### check module ############# ',sys.argv[1:]
src  = connect(sys.argv[1:])
#print '######### check module ############# src =', src
#fff = os.system(zzz+ ' > ' + tmp)
#sys.exit()
ffiles=src

# tmp = 'makedat'
# oxx = string.split(open(tmp,'rt').read(),'\n')
# ### Find all fortran files.
# ffiles=[]
# #ix=0
# for iline in oxx:
# #	ix=ix+1
# #	if(iline[0:3]!='ifc'): continue
# 	offf = re.split(' +',iline)
# 	for iff in offf[1:]:
# 		if(re.search('\.o',iff)):
# 			ffile = re.sub('\.o','.F',iff)
# 			ffiles.append(ffile)
# #			print ffile
#ffiles = uniqlist(ffiles)
ffiles = re.split(' *',ffiles)[1:]
#print ffiles


####################################################
# test
#ffiles=[]
#ffiles= ['../gwsrc/genallcf_mod.F']
#ffiles= ['../gwsrc/tetwt4.F','../main/h_uumatrix.m.F','../gwsrc/pointops.F','../gwsrc/genallcf_mod.F','../gwsrc/readeigen.F', '../main/hx0fp0.m.F']
#ffiles= ['../gwsrc/tetwt4.F']#,'../gwsrc/pointops.F','../main/h_uumatrix.m.F']#,'../gwsrc/genallcf_mod.F','../gwsrc/readeigen.F', '../main/hx0fp0.m.F']
#ffiles = re.split(' +',ffiles)
#print ffiles
####################################################



### skip comments ##############################
#fall=[]
sdef={}
mdef={}
for ffilein in ffiles:
	print 
	print ffilein ," -------------------------"
	flines = string.split(open(ffilein,'rt').read(),'\n')
	flines = skipc(flines,ffilein) #flines are list

### find subroutines and range
#for ff in fall: #ff contains line data for each file
	level=0
	calllist=[]
	slines=[]
	sstack=[]
	fdat= flines[0]
	fff = ['main',fdat[0],fdat[1],fdat[3]]
	sstack.append(fff)

# fdat[0]   FileNM
# fdat[1]   line number
# fdat[2]   line (comment removed)
# fdat[3]   line data  
	for fdat in flines:
		fline =fdat[2]
		num =  fdat[1]
		slines.append([num,fline])
#		slines[num]=fline
#		print fline
		if(  re_mod.search(fline)):
			level=level+1
			subname= re.sub('^ *','',rmod.split(fline)[1])
			subname= re.split("( |^Z)",subname)[0]
			fff= ['mod: '+subname,fdat[0],fdat[1],fdat[3]]
			sstack.append(fff)
		elif(re_end_mod.search(fline)): #module
			level=level-1
			fff = sstack.pop()
			subname = string.lower(fff[0])  #[5:]
			mdef[subname]= fff[1:3] + [fdat[1]] + [fff[3]] +[''] +['']#add range end
		elif(re_program.search(fline) ):
			level=level+1
			subname= string.split( rprog.split(fline)[1]," " )[1]
			subname= rsp(subname)
			fff= ['mai: '+subname,fdat[0],fdat[1],fdat[3]]
			sstack.append(fff)
		elif(re_subroutine.search(fline) ):
			level=level+1
			subname= string.split( rsub.split(fline)[1],"(" )[0]
			subname= rsp(subname)
			fff= ['sub: '+subname,fdat[0],fdat[1],fdat[3]]
			sstack.append(fff)
		elif(re_entry.search(fline) ):
			level=level
			subname= string.split( rent.split(fline)[1],"(" )[0]
			subname= rsp(subname)
			fff= ['sub: '+subname,fdat[0],fdat[1],fdat[3]]
			sstack.append(fff)
		elif(re_function.search(fline) or re_function2.search(fline) ):
			level=level+1
			subname= string.split( reff.split(fline)[1],"(" )[0]
			subname= rsp(subname)
			fff= ['fun: '+subname,fdat[0],fdat[1],fdat[3]]
			sstack.append(fff)
		elif(re_end.search(fline)):
			level=level-1
			fff = sstack.pop()
#			print ' fline =',fline
#			print ' sstack =',fff
 			if(sdef.has_key(fff[0])):
				fff[0]=fff[0]+'+'
				while(1) :
					if(sdef.has_key(fff[0])):
						fff[0]=fff[0]+'+'
						continue
					break
# 				print ' subroutine:: ',fff[0],sdef[fff[0]]
# 				print ' subroutine:: ',[fdat[0],fdat[1],fdat[3]]
# 				print " ERROR: duplicated subroutine name ..."
# 				for x in sdef:
# 					print x
# 				sys.exit()
			subname = string.lower(fff[0]) #[5:]
			sdef[subname]= fff[1:3] + [fdat[1]] + [fff[3]] + [calllist] + [slines] #add range end
			slines=[]
			calllist=[]
		elif(re_call.search(fline)):
			subname= rsp( string.split( rcall.split(fline)[1],"(" )[0] )
			fff= [subname,fdat[0],fdat[1],fdat[3]]
			calllist.append(fff)
		elif(re_use.search(fline)):
			subname= re.sub('^ *','',ruse.split(fline)[1])
			subname= rsp( re.split("(,|^Z| )", subname)[0] )
			fff= [subname,fdat[0],fdat[1],fdat[3]]
			calllist.append(fff)
		elif(re_blockdata.search(fline) ):
			level=level+1
			subname= rsp( re.split("(\  |\n)", rblock.split(fline)[1])[0] )
			fff= [subname,fdat[0],fdat[1],fdat[3]]
			sstack.append(fff)
		else:
			continue
#		print  (' %5i ' % fdat[1]) + ( '%2i ' %  level),fdat[2]
		continue

		if( re_contains.search(fline)):	level=level+1

##########################################
# print '################## module ##########################'
# for i in mdef:
#  	print i+' : ',mdef[i][0:3]
#  	print mdef[i][3]
#  	print 'def@  ', i, mdef[i][0], mdef[i][1], mdef[i][2]
sdef.update(mdef)
skeys=sdef.keys()
skeys.sort()


### find functions #######################
fun=[]
rfun='('
init=1
for k in sdef.keys():
	if(k[0:5]=='fun: '):
		fun = fun + [k[5:]]
		if(init==1):
			aa=''
			init=0
		else:
			aa='|'
		rrr  = re.sub( '\+','',k[5:] )
		rfun = rfun + aa + rrr
rfun = '(?<=\W)'+ rfun +')(?=\W)'
print 'rfun=',rfun
rfff = re.compile(rfun)
#print ' fun=',fun
#sys.exit()


###########find function calls ###############################
for i in skeys:
	filen= sdef[i][0] 
	slines=sdef[i][5]
#	init = sdef[i][1]
#	end  = sdef[i][2]
#	print 'def@  ', i,filen,init,end 
#	for ifun in fun:
#		slines = re.sub(ifun,'~'+ifun+'~',slines)
	for ix in slines:
		num= ix[0]
#		aaa=re.findall(rfun,ix[1])
		aaa= rfff.findall(ix[1])
		output = []
#		print 'aaa',ix
		for aaax in aaa:
			if(not (aaax==i[5:])):
				output.append(aaax)
		if(not (output==[])):
			for iout in output:
#				print num,iout,ix[1]
				fff= [iout,filen,num]
				sdef[i][4].append(fff)
# This is full data
#	for line in sdef[i][5]:
#		print line # iline =  [line number, line text]

#	print sdef[i][3],
#	for ic in sdef[i][4]:
#		print 'cal@  ', i,' ',ic[0],ic[1],ic[2]


### Write definition and call. ####
print '#############################################'
for i in skeys:
	filen= sdef[i][0] 
	init = sdef[i][1]
	end  = sdef[i][2]
	print
	print 'def@  ', i,filen,init,end 
# This is full data
#	for line in sdef[i][5]:
#		print line # iline =  [line number, line text]

#	print sdef[i][3],
	for ic in sdef[i][4]:
		print 'cal@  ', i,' ',ic[0],ic[1],ic[2]


### remove header from sdef ############
for i in skeys:
	ii  = i[5:]
	sss = sdef[i]
	del sdef[i]
	sdef[ii]= sss

### tree generation ####################
i='lmfp'
once = 1
print '###### call tree for ', i,' ##########'

calling1 = sdef[i][4]
for ic1 in calling1:
	print 'tree0 ', ic1[0],ic1[1],ic1[2]
	ir1 = ic1[0]
	try:
		calling2= sdef[ir1][4]
	except:
		continue
	
	icalled2=[]
	for ic2 in calling2:
		ir2 = ic2[0]
		if once==1 and ir2 in icalled2:	continue
		icalled2.append(ir2)
		print 'tree1  |  ',ic2[0],ic2[1],ic2[2]
		try:
			calling3= sdef[ir2][4]
		except:
			continue

		icalled3=[]
		for ic3 in calling3:
			ir3 = ic3[0]
			if once==1 and ir3 in icalled3:	continue
			icalled3.append(ir3)
			print 'tree2  |   | ',ic3[0],ic3[1],ic3[2]
			try:
				calling4= sdef[ir3][4]
			except:
				continue

			icalled4 = []
			for ic4 in calling4:
				ir4 = ic4[0]
				if once==1 and ir4 in icalled4:	continue
				icalled4.append(ir4)
				
				print 'tree3  |   |  | ',ic4[0],ic4[1],ic4[2]
				try:
					calling5= sdef[ir4][4]
				except:
					continue


sys.exit()





