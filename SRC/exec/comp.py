#!/usr/bin/env python3
import os,sys,re
def runprogs(runlist):
    for irun in runlist:
        print(irun)
        err = os.system(irun)
        if(err): print('Error exit!')
        if(err): sys.exit(-1)
        
def comp(f1,f2,label,tol,key,key2=None):
    ix=0
    if(key2): ix=1
    val=[-9999,-9999]
    for ifnum,ifi in enumerate([f1,f2]):
        ifile=open(ifi,'rt').read().split('\n')
        for i,line in enumerate(ifile): #print(line.split(' '))
            if(re.search(key,line)):
                if(ix==0 or (ix!=0 and re.search(key2,line))):
                    line=re.sub(r'=\s+','=',line)
                    line=re.sub(r':\s+',':',line)
                    #print(key,line)
                    val[ifnum]= float(re.split(key,line)[1].split(' ')[0])
    #sys.exit()
    if(val[0]==-9999 or val[1]==-9999): return 'skip! '
    diff=val[0]-val[1]
    out='ERR! '
    if(abs(diff)<float(tol)): out='OK! '
    keys=label+(26-len(label))*' '+out
    keys=keys+' '+str(val[0])+' '+str(val[1])
    keys=keys+(60-len(keys))*' '+' tol='+str(tol)+' diff='+str(diff)+' Key1:'+key
    if(ix!=0): keys= keys+ ' Key2:'+key2
    print(keys)
    return out

def compall(f1in,f2in,tol):
    f1= open(f1in,'rt').read().split('\n')
    f2= open(f2in,'rt').read().split('\n')
    diff=0
    for ifnum,ifi in enumerate(f1):
        #print(f1[ifnum])
        iline1=re.split(r'\s+',f1[ifnum])
        iline2=re.split(r'\s+',f2[ifnum])
        #print(ifnum, iline1, iline2)
        if re.match(r'^\s*#', f1[ifnum]): continue
        #if(iline1[0][0:1]=='#'): continue
        iline1=[float(i) for i in iline1 if i!=''] #and float(i)!=0]
        iline2=[float(i) for i in iline2 if i!=''] #and float(i)!=0]
        #print(ifnum, iline1, iline2)
        for i,idat1 in enumerate(iline1):
            diff=max(idat1-iline2[i],diff)
    if(diff>float(tol)) :
        print('ERROR: max deviation =',diff,' tolerance =',tol)
        return 'ERR! '
    if(diff<=float(tol)):
        print('max deviation =',diff,' tolerance =',tol)
        return 'OK! '
    
def compeval(f1in,f2in,key,lineeval,evalso,tol):
    f1= open(f1in,'rt').read().split('\n')
    f2= open(f2in,'rt').read().split('\n')
    for ifi in [1,2]:
        if(ifi==1): fdat=f1
        if(ifi==2): fdat=f2
        ifg=0
        for ifnum,idat in enumerate(fdat):
            if(re.search(key,idat)): ifg=1
            if(ifg>0): ifg=ifg+1
            if(ifg==int(lineeval)+1):
                ival=[i for i in re.split(r'\s+',idat) if i!='']
                if(ifi==1):
                    evhomo1=float(ival[int(evalso)-1])
                    evlumo1=float(ival[int(evalso)])
                    diff1=(evlumo1-evhomo1)*13.605
                if(ifi==2):
                    evhomo2=float(ival[int(evalso)-1])
                    evlumo2=float(ival[int(evalso)])
                    diff2=(evlumo2-evhomo2)*13.605
    print("Spin-Orbit splitting =",diff2,' eV ( states ',evhomo2,' and ',evlumo2,')')
    print("                      ",diff1,' eV ( states ',evhomo1,' and ',evlumo1,')')
    if(abs(diff1-diff2)<float(tol)):
        out='OK!'
        aaa='PASSED! '+f2in
    else:    
        out='ERR!'
        aaa='FAILED at '+f2in
    with open("summary.txt", "a") as aout: print(aaa, file=aout)
        
    return out

def diffnum(f1,f2,tol,comparekeys):
    import diffnum0
    file1 = open(f1,'rt').read()
    file2 = open(f2,'rt').read()
    errmax = diffnum0.comparenum(tol,file1,file2, comparekeys, printsw=0)
    if(errmax<tol):
        print(' Comparison OK!  MaxDiff=',errmax,'< tol=',tol,' for ', os.path.basename(f2))
        out='ok! '
        aaa='PASSED! '+f2 #.split('work/')[1]
    else:
        print()
        print() 
        print(' Error! MaxDiff=',errmax,'> tol=',tol,' for ', os.path.basename(f2))
        out='err! '
        aaa='FAILED at '+f2 #.split('work/')[1]
    with open("summary.txt", "a") as aout: print(aaa, file=aout)
    return out

def compareqpu(qpu1,qpu2,printsw):
	pr = (printsw==1)
	oxx= qpu1.split('\n') 
	oyy= qpu2.split('\n')
	errmax=0.0
	for ix in range( max(len(oxx),len(oyy))):
		iline=oxx[ix]
		ilin2=oyy[ix]
		if ix==5:
			if pr: print (iline)
		if ix>=6:
			if pr: print (iline[0:32],end='')
			try:
				for iw in range(4,17,1):
					w1=float(iline.split()[iw])
					w2=float(ilin2.split()[iw])
					if (iw >=15) & pr : print( '%9.5f' % (w1-w2),end='')
					if (iw <=14) & pr : print( '%6.2f' % (w1-w2),end='')
					if( abs(w1-w2)>errmax ): errmax=abs(w1-w2)
					#errmax=errmax + abs(w1-w2)
				if pr: print()
			except:
				if pr: print(iline)
	return errmax

def dqpu(f1,f2):
        qpu1 = open(f1,'rt').read()
        qpu2 = open(f2,'rt').read()
        errmax = compareqpu(qpu1,qpu2, printsw=0)
        #errmax = compareqpu(qpu1,qpu2, printsw=0)
        etol=1.1e-2
        if(errmax<etol):
                print(' Comparison OK! max(abs(QPU-QPU))=',errmax,' <etol=',etol,'for',os.path.basename(f2))
                out='ok! '
                aaa='PASSED! '+ f2 #f2.split('work/')[1]
        else:
                errmax = compareqpu(qpu1,qpu2, printsw=1)
                print() 
                print(' Error! max(abs(QPU-QPU))=',errmax,'>etol=',etol,'for',os.path.basename(f2))
                out='err! '
                aaa='FAILED at '+f2 #.split('work/')[1]
        with open("summary.txt", "a") as aout: print(aaa, file=aout)
        return out

defatol1=1e-5
dehf1tol1=1e-5
dfmax1tol1=0.1
dmom1tol1=1e-4
dehf1toln=1e-5
drmsqtol1=1e-4
dosclstol=0.001 
dorbmtol=1e-5
    
def test1_check(f1,f2):
    print('compare '+f1 +' and '+f2)
    test=  comp(f1,f2,'FA etot (last species)  ',defatol1, 'etot=' )
    test+= comp(f1,f2,'1st  iter ehf.eV        ',dehf1tol1,r'ehf\(eV\)=','^h')
    test+= comp(f1,f2,'1st  iter ehk.eV        ',dehf1tol1,r'ehk\(eV\)=','^h')
    test+= comp(f1,f2,'1st  iter mmom          ',dmom1tol1,'mmom=','^h')
    test+= comp(f1,f2,'2nd  iter ehf.Ry        ',dehf1toln,'ehf=','it  2')
    test+= comp(f1,f2,'2nd  iter ehk.Ry        ',dehf1toln,'ehk=','it  2')
    test+= comp(f1,f2,'9th  iter ehf.Ry        ',dehf1toln,'ehf=','it  9')
    test+= comp(f1,f2,'9th  iter ehk.Ry        ',dehf1toln,'ehk=','it  9')
    test+= comp(f1,f2,'last iter x ehf.eV      ',dehf1toln,r'ehf\(eV\)=','^x')
    test+= comp(f1,f2,'last iter x ehk.eV      ',dehf1toln,r'ehk\(eV\)=','^x')
    test+= comp(f1,f2,'last iter c ehf.eV      ',dehf1toln,r'ehf\(eV\)=','^c')
    test+= comp(f1,f2,'last iter c ehk.eV      ',dehf1toln,r'ehk\(eV\)=','^c')
    test+= comp(f1,f2,'last iter E(LDA+U).eV   ',dehf1toln,r'Etot\(LDA+U\)=')
    test+= comp(f1,f2,'last iter max force     ',dfmax1tol1,'Maximum Harris force =')
    test+= comp(f1,f2,'last iter mmom          ',dmom1tol1,'mmom=')
    test+= comp(f1,f2,'chk1ch last iter RMS dq ',drmsqtol1,'RMS DQ=')
    test+= comp(f1,f2,'Orbital moment          ',dorbmtol,'total orbital moment   1:')
    test+= comp(f1,f2,'last iter ehf.eV E(MTO+PW)',dehf1toln,r'pwmode=[^0].*ehf\(eV\)=')
    test+= comp(f1,f2,'last iter ehk E(MTO+PW)  ',dehf1toln, r'pwmode=[^0].*ehk\(eV\)=')
    #print(test)
    if('ERR!' in test) :
        aaa=' FAILED! TEST 1 compare files:'+f1+' and '+f2
        out='err! '
    else: 
        aaa='PASSED! TEST 1 '+os.path.basename(f2)
        out='ok! '
    with open("summary.txt", "a") as aout: print(aaa, file=aout)
    return out

def test2_check(f1,f2,tol=dosclstol):
    test=compall(f1,f2,tol)
    if('ERR!' in test) :
        aaa='FAILED! TEST 2 comparison comparison files:'+f1+' and '+f2
        out='err! '
    else: 
        aaa='PASSED! TEST 2 '+os.path.basename(f2)
        out='ok! '
    with open("summary.txt", "a") as aout: print(aaa, file=aout)
    return out
