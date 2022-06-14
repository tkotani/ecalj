#!/usr/bin/env python3
import re,sys #,string
def vasp2ctrl_alat(vaspread): #read alat from VASP file
    alat_list = vaspread[1].split()
    alat_val = float( alat_list[0] )/0.529177 # convert it to a.u. from angstrome.
    return alat_val    

def vasp2ctrl_atomcount(vaspread): # set atom namet in all_atom[]
    all_atom = []
    atom_pat = vaspread[5].split()
    atom_num = vaspread[6].split()
    for number in range(len(atom_pat)):
        for quant in range( int(atom_num[number]) ):
            all_atom.append(atom_pat[number])
    NBAS_val = len(all_atom)
    return all_atom,NBAS_val

def vasp2ctrl_write(vaspread,alat_val,NBAS_val,atom_list,titleinput,all_atom,ratioa):
    outputname = titleinput+'.vasp2ctrl'
    ctrlwrite = open(outputname,'w')
    tes = all_atom
    avec= [float(x)/ratioa for x in re.split(r'[\s+]',vaspread[2]) if x!='']
    bvec= [float(x)/ratioa for x in re.split(r'[\s+]',vaspread[3]) if x!='']
    cvec= [float(x)/ratioa for x in re.split(r'[\s+]',vaspread[4]) if x!='']
    avecs = ' '.join(['%18.11f ' % x for x in avec])
    bvecs = ' '.join(['%18.11f ' % x for x in bvec])
    cvecs = ' '.join(['%18.11f ' % x for x in cvec])
    plats = 'PLAT=  '+ avecs +' \n' +' '*12+ bvecs +'\n'+' '*12+ cvecs +'\n'
    outputlist1 = 'STRUC'+'\n'+' '*5+'ALAT=' + str(alat_val*ratioa) +'\n' \
    +' '*5+ plats \
    +' '*2+'NBAS='+str(NBAS_val)+'\n'+'SITE'+'\n'
    ctrlwrite.write(outputlist1)
    for sort in range(len(tes)):
        position  = [float(x)/ratioa for x in atom_list[sort][0:3] ]
        positions = ' '.join(['%18.11f ' % x for x in position])
        atoms = '%s' % tes[sort]
        outputlist2 = ' '*5+'ATOM='+atoms+' '+'POS='+positions +'\n'
        ctrlwrite.write(str(outputlist2))
    ctrlwrite.close()

def vasp2ctrl_atom(vaspread,alat_val,NBAS_val,plat1,plat2,plat3):
    atom_list = []
    angs = alat_val*.529177
    for i in range(NBAS_val):
        atom_val = vaspread[8+i].split()
        atom_list.append(atom_val)
    plat1= [float(plat1[i]) for i in range(len(plat1)) ] 
    plat2= [float(plat2[i]) for i in range(len(plat2)) ] 
    plat3= [float(plat3[i]) for i in range(len(plat3)) ] 

    for atom1 in range(len(atom_list)):
        x = float(atom_list[atom1][0])
        y = float(atom_list[atom1][1])
        z = float(atom_list[atom1][2])

        if vaspread[7].split()[0].lower()=='cartesian':
            atom_list[atom1][0] = x #/alat_val/.529177 #for new VASP
            atom_list[atom1][1] = y #/alat_val/.529177
            atom_list[atom1][2] = z #/alat_val/.529177
            print('read cartesian poscar',x,y,z)
        if vaspread[7].split()[0].lower()=='direct':
            atom_list[atom1][0] = x*plat1[0] + y*plat2[0]  + z*plat3[0]
            atom_list[atom1][1] = x*plat1[1] + y*plat2[1]  + z*plat3[1]
            atom_list[atom1][2] = x*plat1[2] + y*plat2[2]  + z*plat3[2]
            print('read direct(fractional) poscar',x,y,z)

    return atom_list

# for ctrl2vasp.
def fileopen(openfile): #return cleanuped contents of openfile.
    newctrl = []
    newctrl2 = []
    for readctrl in range(len(openfile)):
        sharp_remove = openfile[readctrl].split('#')
        newctrl.append(sharp_remove[0]) # remove comment lines
    for i in range(len(newctrl)):
        eee = re.sub(r'\s*=\s* ', '=', newctrl[i]) # remove spaces sandwitch =.
        newctrl2.append(eee)
    return newctrl2

def constlist(openfile): #---get input values from ctrls (openfile)
    keywordlist = []
    for open_read in openfile:
        hx = re.match('^',open_read)
        open_read = re.sub('\^','**',open_read) # x^y means x**y
        rx = re.match('/', open_read)
        open_read = re.sub('/', '/1.0/', open_read) # python /1.0/
        fx = re.match(r'\s*=\s* ',  open_read) 
        open_read = re.sub(r'\s*=\s* ', '=', open_read)
        ex = re.match(r'^\%\s*const ',  open_read) 
        open_read = re.sub(r'^\%\s*const ', '', open_read)
        const_temp = open_read.split('#') # After '#', neglected.
        const_ex = const_temp[0]
        if ex: 
            const_val = const_ex.split(' ')
            for const_form in range(len(const_val)):
                if const_val[const_form] != "" :
                    keywordlist.append(const_val[const_form])
    return keywordlist

def keywordname(keyword_val): #geve valus from keyword
    keyword_temp = []
    for const_name in range(len(keyword_val)):
        splitword = keyword_val[const_name].split('=')
        keyword_temp.append(splitword[0])
    return keyword_temp

def alat(perfectopen,keyword_name): #get alat (unit)
    for i in range(len(perfectopen)):
        if perfectopen[i].find("ALAT") != -1: #find alat
            break
    struc_rem = re.sub('STRUC' , '' ,perfectopen[i]) # remove "STRUC"
    alat_temp = struc_rem.split()
    ALATone = []
    for follow_struc in range(len(alat_temp)):
        alat_val = alat_temp[follow_struc].split('=') # ALAT=foobar --> alat,foobar
        ALATone.append(alat_val)
    for alat_first in range(len(ALATone)):
        for alat_second in range(len(ALATone[alat_first])):
            for const_name in range(len(keyword_name)):
                ALATone[alat_first][alat_second] = re.sub('\{'+keyword_name[const_name]+'\}',keyword_name[const_name],ALATone[alat_first][alat_second]) #replace {foobar} from ALAT
    return ALATone

def plat(perfectopen,keyword_name): #---get PLAT
    fileall=''
    for ix in perfectopen:
        fileall +=ix+' '
#    print fileall    
#    print perfectopen.connect()
#    plat_temp = perfectopen[i].split('PLAT=')
    plat_temp = fileall.split('PLAT=')
#    print plat_temp[1]
    PLAT_list = plat_temp[1].split()[0:9]
#    print PLAT_list
#    sys.exit()

#    for i in range(len(perfectopen)):
#        if perfectopen[i].find("PLAT") != -1: #---find PLAT
#            break
#    plat_temp = perfectopen[i].split('=')
#    PLATone = []
#    for follow_struc in range(len(plat_temp)):
#        plat_val = plat_temp[follow_struc].split()
#        PLATone.append(plat_val)
#    PLAT_list = PLATone[-1] # last element constains PLAT.
    for line_ing in range(len(PLAT_list)):
        for const_name in range(len(keyword_name)):
            PLAT_list[line_ing] = re.sub('\{'+keyword_name[const_name]+'\}',keyword_name[const_name],PLAT_list[line_ing]) # replace {foobar} with its values
    return PLAT_list

def atom(perfectopen,keyword_name): #---get ATOM,POS
    temp1 = []
    temp3 = []
    temp5 = []
    for i in range(len(perfectopen)):
        temp2 = perfectopen[i].strip().split()
        temp1.append(temp2)
    for i in range(len(temp1)):
        for ii in range(len(temp1[i])):
            if re.match("POS", temp1[i][ii]) != None : # look for POS
                temp3.append(temp1[i])
    if re.match('SITE', temp3[0][0]) != None : # wether POS is in SITE or not.
        temp3[0].pop(0) 
    for i in range(len(temp3)):    
#        print 'zzzzz ',temp3[i]
        vvv1 = re.sub('ATOM=','',temp3[i][0])
#        print 'atom=',vvv1
        temp5.append(vvv1)
        for ii in range(len(temp3[i])):
            dat= temp3[i][ii]
            if re.match("POS=",dat):
                temp5.append(re.sub('POS=','',dat))
                temp5.append(temp3[i][ii+1])
                temp5.append(temp3[i][ii+2])
    for line_each in range(len(temp5)):
#        print 'pppp',temp5[line_each]
        for const_name in range(len(keyword_name)):
            temp5[line_each] = re.sub('\{'+keyword_name[const_name]+'\}',keyword_name[const_name],temp5[line_each]) 
    sep_atom = []
    count = int(len(temp5)/4) 
    j1 = -1
    for j in range(count):
        j1 = j1+1
        j2 = 4*j1
        sep_list = temp5[j2:j2+4]
        sep_atom.append(sep_list) 
    return sep_atom

def savefile(ALAT,PLAT_list,atomlist,titleinput,coordinates): 
    title = titleinput
    output = title+'.vasp' 
    writefile = open(output,'w')
    outputlist1 = title+'\n'+str(ALAT)+'\n'
    writefile.write(outputlist1)
    for j in range(len(PLAT_list)):
        outputlist2 = str(PLAT_list[j])+' '
        writefile.write(outputlist2)
        if (j+1)%3 == 0: 
            writefile.write('\n')
    atomo=''
    for i in range(len(atomlist)):
        if atomo != atomlist[i][0]:        
            outputlist3 = atomlist[i][0]+' ' 
            writefile.write(outputlist3)
        atomo=atomlist[i][0]
    writefile.write('\n')
    natom=0
    atomo=atomlist[0][0]
    for k in range(len(atomlist)):
        natom=natom+1
        if atomo != atomlist[k][0]:
            natomx= natom-1
            outputlist4 = '%i ' % natomx
            writefile.write(outputlist4)
            natom=1
        atomo=atomlist[k][0]
    if natom >0:
        outputlist4 = '%i ' % natom
        writefile.write(outputlist4)

    outputlist5 = '\n'+coordinates+'\n' #---direct,cartesian
    writefile.write(outputlist5)
    for j in range(len(atomlist)):
        outputlist6 = str(atomlist[j][1])+' '+str(atomlist[j][2])+' '+str(atomlist[j][3])+'\n'
        writefile.write(outputlist6) 
    writefile.close()


