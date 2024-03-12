#!/usr/bin/env python3
import os,sys,re
#print(sys.argv[1:])
#f1= sys.argv[1]
#f2= sys.argv[2]
#label=sys.argv[3]
#tol=sys.argv[4]
#key=sys.argv[5]
def comp(f1,f2,label,tol,key,key2=None):
#    ix=0
#    for i,dat in enumerate(sys.argv[:]):
#        if(dat=='-v'):
#            ix=i+1
#            break
#    if(ix!=0): key2=sys.argv[ix]
    #print(len(f1))
    ix=0
    if(key2): ix=1
    val=[-9999,-9999]
    for ifnum,ifi in enumerate([f1,f2]):
        ifile=open(ifi,'rt').read().split('\n')
        for i,line in enumerate(ifile): #print(line.split(' '))
            if(re.search(key,line)):
                if(ix==0 or (ix!=0 and re.search(key2,line))):
                    line=re.sub('=\s+','=',line)
                    line=re.sub(':\s+',':',line)
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

#import os,sys,re
#f1= open(sys.argv[1],'rt').read().split('\n')
#f2= open(sys.argv[2],'rt').read().split('\n')
#tol= sys.argv[3]
def compall(f1in,f2in,tol):
    f1= open(f1in,'rt').read().split('\n')
    f2= open(f2in,'rt').read().split('\n')
    diff=0
    for ifnum,ifi in enumerate(f1):
        iline1=re.split('\s+',f1[ifnum])
        iline2=re.split('\s+',f2[ifnum])
        #print(ifnum, iline1, iline2)
        if(iline1[0][0:1]=='#'): continue
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
    #key = sys.argv[3]
    #lineeval=sys.argv[4]
    #evalso=sys.argv[5]
    #tol=sys.argv[6]
    for ifi in [1,2]:
        if(ifi==1): fdat=f1
        if(ifi==2): fdat=f2
        ifg=0
        for ifnum,idat in enumerate(fdat):
            if(re.search(key,idat)): ifg=1
            if(ifg>0): ifg=ifg+1
            if(ifg==int(lineeval)+1):
                ival=[i for i in re.split('\s+',idat) if i!='']
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
    else:    
        out='ERR!'
    return out

######## diffnum
def diffnum(f1,f2,tol,comparekeys):
    import diffnum0
    # ############### main ###################################
    # try:
    # 	print( "   Readin two files =     ", sys.argv[1], sys.argv[2])
    # except:	
    # 	print( "   take difference of numbers in two files")
    # 	print( "   usage: diffnum FILE1 FILE2 'comparekeys' ")
    # 	sys.exit()
    file1 = open(f1,'rt').read()
    file2 = open(f2,'rt').read()
    #comparekeys=sys.argv[3:]
    #print( ' Comparekeys=',comparekeys)
    errmax = diffnum0.comparenum(tol,file1,file2, comparekeys, printsw=0)
    #print( 'end of comparenum')
    if(errmax<tol):
        print(' Comparison OK!  MaxDiff=',errmax,'< tol=',tol,' for ', os.path.basename(f2))
        out='ok! '
        aaa='ok! for '+f2.split('work')[1]
    else:
        print()
        print() 
        print(' Error! MaxDiff=',errmax,'> tol=',tol,' for ', os.path.basename(f2))
        out='err! '
        aaa='failed at '+f2.split('/work')[1]
    with open("../summary.txt", "a") as aout: print(aaa, file=aout)
    return out
