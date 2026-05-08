### Get superlattice structure ###########
# How to build a superlattice
# 1. Choose two lattices, e.g. InAs and GaSb
# 2. Choose the number of layers in the superlattice, e.g. h1=3 and h2=3



#data from Van de Walle Phys. Rev. B 39, 1871 (1989) Table I.  The microsoft copilot did mistakes -->corrected by hand.
data = {
    "SiSi": {"a": 5.43, "c11": 1.675, "c12": 0.650, "c44": 0.801, "D001": 0.776, "G001": 3.641, "D110": 0.515, "G110": 4.417, "D111": 0.444, "G111": 4.628},
    "GeGe": {"a": 5.65, "c11": 1.315, "c12": 0.494, "c44": 0.684, "D001": 0.751, "G001": 2.876, "D110": 0.450, "G110": 3.570, "D111": 0.371, "G111": 3.751},
    "GaAs": {"a": 5.65, "c11": 1.223, "c12": 0.571, "c44": 0.600, "D001": 0.934, "G001": 2.522, "D110": 0.580, "G110": 3.359, "D111": 0.489, "G111": 3.574},
    "AlAs": {"a": 5.65, "c11": 1.250, "c12": 0.534, "c44": 0.542, "D001": 0.854, "G001": 2.656, "D110": 0.616, "G110": 3.207, "D111": 0.550, "G111": 3.361},
    "InAs": {"a": 6.08, "c11": 0.833, "c12": 0.453, "c44": 0.396, "D001": 1.088, "G001": 1.587, "D110": 0.674, "G110": 2.306, "D111": 0.570, "G111": 2.487},
    "GaPp": {"a": 5.43, "c11": 1.439, "c12": 0.652, "c44": 0.714, "D001": 0.906, "G001": 3.000, "D110": 0.559, "G110": 3.953, "D111": 0.470, "G111": 4.198},
    "AlPp": {"a": 5.43, "c11": 1.320, "c12": 0.630, "c44": 0.615, "D001": 0.955, "G001": 2.697, "D110": 0.623, "G110": 3.554, "D111": 0.536, "G111": 3.778},
    "InPp": {"a": 5.87, "c11": 1.022, "c12": 0.576, "c44": 0.460, "D001": 1.127, "G001": 1.897, "D110": 0.727, "G110": 2.768, "D111": 0.625, "G111": 2.990},
    "GaSb": {"a": 6.08, "c11": 0.908, "c12": 0.413, "c44": 0.445, "D001": 0.910, "G001": 1.891, "D110": 0.569, "G110": 2.482, "D111": 0.480, "G111": 2.635},
    "AlSb": {"a": 6.08, "c11": 0.877, "c12": 0.434, "c44": 0.408, "D001": 0.990, "G001": 1.763, "D110": 0.641, "G110": 2.372, "D111": 0.550, "G111": 2.530},
    "InSb": {"a": 6.48, "c11": 0.659, "c12": 0.356, "c44": 0.300, "D001": 1.080, "G001": 1.261, "D110": 0.698, "G110": 1.785, "D111": 0.600, "G111": 1.920},
    "ZnSe": {"a": 5.65, "c11": 0.826, "c12": 0.498, "c44": 0.400, "D001": 1.206, "G001": 1.447, "D110": 0.716, "G110": 2.340, "D111": 0.597, "G111": 2.556},
    "ZnSs": {"a": 5.40, "c11": 1.067, "c12": 0.666, "c44": 0.456, "D001": 1.248, "G001": 1.803, "D110": 0.814, "G110": 2.845, "D111": 0.704, "G111": 3.109},
    "ZnTe": {"a": 6.08, "c11": 0.713, "c12": 0.407, "c44": 0.312, "D001": 1.142, "G001": 1.311, "D110": 0.751, "G110": 1.907, "D111": 0.651, "G111": 2.060},
    "CdTe": {"a": 6.48, "c11": 0.562, "c12": 0.394, "c44": 0.206, "D001": 1.402, "G001": 0.807, "D110": 0.974, "G110": 1.386, "D111": 0.863, "G111": 1.535},
    "HgTe": {"a": 6.48, "c11": 0.597, "c12": 0.415, "c44": 0.226, "D001": 1.390, "G001": 0.870, "D110": 0.949, "G110": 1.499, "D111": 0.837, "G111": 1.660}
}
data['AsGa'] = data['GaAs'] # data is added so that data works even for AsGa (gives the same result as GaAs), others as well.
data['AsAl'] = data['AlAs']
data['AsIn'] = data['InAs']
data['PpGa'] = data['GaPp']
data['PpAl'] = data['AlPp']
data['PpIn'] = data['InPp']
data['SbGa'] = data['GaSb']
data['SbAl'] = data['AlSb']
data['SbIn'] = data['InSb']
data['SeZn'] = data['ZnSe']
data['SsZn'] = data['ZnSs']
data['TeZn'] = data['ZnTe']
data['TeCd'] = data['CdTe']
data['TeHg'] = data['HgTe']    

def cspec(aa):
    aout=aa
    if aa=='Pp': aout='P'
    if aa=='Ss': aout='S'
    return aout

def superlattice(lattice1,lattice2, h1,h2):
    # We first build ZB no-deformed superlattice along 001,
    # we identify layers by z coordinate. 
    # Note that 001 contains layers, which contains one atom per layer,
    # For example, 4+4 model contains 16 layers for 001.
    # 110 as well, # while 110 contains two atoms per layer. Thus  8 layers for 110
     
    # Problem is how to determined the layer spacing.
    # We use the logic of van de Walle, Phys. Rev. B 39, 1871 (1989). Eq.(1a) and (2a)
    # We first determine ainplane (instead of two layer case at Eq.(1a), we consider 16 layers for 001 ).
    # Spacing between layers contains cation-anion bond, corresponding to the bulk material.

    # In the case of 111, we have two atoms per layer. Thus two types of bonds can be at the boundaries.
    # That is the bondary between GaSb layer and InAs layer contains GaAs and InSb bonds--> In this case,
    # we take an average of the two bulk materials GaAs and InSc.
    '''
    Generate superlattice structures based on the given parameters.
    lattice1: str, first lattice name (e.g. 'InAs')
    lattice2: str, second lattice name (e.g. 'GaSb')
    h1: int, number of layers in the first lattice
    h2: int, number of layers in the second lattice
    WARN!--->We name Si as SiSi. P-->Pp S-->Ss to simpify the code.
    '''
    import numpy as np
    if((h1+h2)%2!=0):
        print('h1+h2 should be even')
        exit(1)
    ### format
    def fmt(arr):
        return '[{}]'.format(', '.join(f'{x:.6f}' for x in arr))  
    ### ZB lattice
    avec=np.array([0,.5,.5])
    bvec=np.array([.5,0,.5])
    cvec=np.array([.5,.5,0])
    vol0=np.dot(np.cross(avec,bvec),cvec)

    a1=lattice1[0:2]
    a2=lattice1[2:4]
    a3=lattice2[0:2]
    a4=lattice2[2:4]
    for direction in ['001','110']:
        if(direction=='001'):   #otsuka orthogonal superlattice
            AA=cvec
            BB=avec-bvec
            CC=-cvec+bvec+avec
            anion= np.array([1,1,1])/4
        elif(direction=='110'): #hinuma orthogonal superlattice
            AA= avec-bvec
            BB=-avec-bvec+cvec
            CC=cvec
            anion= np.array([1,1,1])/4
        center=(AA+BB+CC)/2
        vol = np.dot(np.cross(AA,BB),CC)
        ea=AA/np.linalg.norm(AA)
        eb=BB/np.linalg.norm(BB)
        ec=CC/np.linalg.norm(CC)
        projec=np.array([[ea[0], ea[1], ea[2]], 
                        [eb[0], eb[1], eb[2]], 
                        [ec[0], ec[1], ec[2]]])
        AA=projec@AA #components along ea,eb,ec
        BB=projec@BB
        CC=projec@CC
        center= projec@center
        anion=  projec@anion
        aabbcc=np.array([
            [AA[0],AA[1],AA[2]],
            [BB[0],BB[1],BB[2]],
            [CC[0],CC[1],CC[2]]
        ])        
        reciprocal= np.linalg.inv(aabbcc).T #reciprocal lattice
        #print('reciprocal=', reciprocal) #OK but trivial for orthogonal case
        center = (reciprocal@center)%1  #fractional coordinates
        anion  = (reciprocal@anion)%1   #fractional coordinates
        fu=vol/vol0
        
        print(f'============================= superlattice {a1}{a2} {a3}{a4} {direction} ============================') 
        print('formula unit vol/vol0=',fu)
        print('total number of atoms=',(h1+h2)*2)
        print(f'Layer width: {h1} {h2}') 
        print(f'C-axis direction:',direction)    
        print('supercell AA=',AA)
        print('supercell BB=',BB)
        print('supcecell CC=',CC,' z axis')
        zero= np.array([0,0,0])
        print('atom positions:')
        if direction=='001':
            print(f' === Four layer per C periodic cell total layer= {h1+h2}*2')
            nlayer=(h1+h2)*2
            print(f' === 001 has one atom per layer. layer is specified by z coordinate')
            print(f'    {fmt(zero)} {a1} layer0')
            print(f'    {fmt(anion)} {a2} layer1')
            print(f'    {fmt(center)} {a1} layer2')
            print(f'    {fmt((center+anion))} {a2} layer3')
            print(f'    nlayer={nlayer}')
            print('-------------------------------------------------------')
            lbond = [None] * nlayer
            layers= [None] * nlayer
            lpos  = [None] * nlayer
            atoms = [None] * nlayer
            i=0
            aa0=a1
            aa1=a2
            ic=0 
            ix=0
            zz1=np.array([0,0,1])
            layerindex=-1
            aay=a4
            while True:
                i+=1
                if(i>h1):
                    aa0=a3
                    aa1=a4
                    ix=1

                layerindex+=1
                lbond [layerindex]=aay+aa0
                lpos  [layerindex]=zero  +ic*zz1
                atoms [layerindex]=aa0
                layers[layerindex]=f'{layerindex} {ix} layer0 {fmt(lpos[layerindex])} {atoms [layerindex]} : interlayerbond {layerindex-1}-{layerindex} {lbond[layerindex]}'
                aay=aa0
                
                layerindex+=1
                lbond [layerindex]=aay+aa1
                lpos  [layerindex]=anion +ic*zz1
                atoms [layerindex]=aa1
                layers[layerindex]=f'{layerindex} {ix} layer1 {fmt(anion +ic*zz1)} {aa1} : interlayerbond {layerindex-1}-{layerindex} {lbond[layerindex]}'
                aay=aa1
                
                i+=1 
                if(i>h1):
                    aa0=a3
                    aa1=a4
                    ix=1
                    
                layerindex+=1
                lbond [layerindex]=aay+aa0
                lpos  [layerindex]=center+ic*zz1
                atoms [layerindex]=aa0
                layers[layerindex] =f'{layerindex} {ix} layer2 {fmt(center+ic*zz1)} {aa0} : interlayerbond {layerindex-1}-{layerindex} {lbond[layerindex]}'
                aay=aa0

                layerindex+=1
                lbond [layerindex]=aay+aa1
                lpos  [layerindex]=center+anion+ic*zz1
                atoms [layerindex]=aa1
                layers[layerindex] =f'{layerindex} {ix} layer3 {fmt(center+anion+ic*zz1)} {aa1} : interlayerbond {layerindex-1}-{layerindex} {lbond[layerindex]}'
                aay=aa1
                
                if(i>h1+h2-2):
                    break
                ic+=1
                
        #WandeWalle method.  We assume layer spacing is 0.25a van de Walle
            aden = sum(                      data[lbond[j]]['G001'] * (data[lbond[j]]['a'] * 0.25) for j in range(nlayer))
            anum = sum( data[lbond[j]]['a'] *data[lbond[j]]['G001'] * (data[lbond[j]]['a'] * 0.25) for j in range(nlayer))
            ain  = anum/aden # Eq.(1a) van de Walle
            #ain = 6.08112
            print('ainplane=',ain) #note direction of deposition of thin file. c axis is perpendicular to the plane of the film 
            aw=0
            aww= np.zeros(nlayer)
            for j in range(nlayer):
                aav=data[lbond[j]]['a']
                Dav=data[lbond[j]]['D001']
                aspace=0.25*aav*( 1  -  Dav * (ain/aav  -1) ) #Eq.(2a) van de Walle
                if(j==0): 
                    aspace0=aspace
                else:
                    aw=aw+aspace
                aww[j]=aw  
                #print(layers[j],f'{lpos[j]} {aw}')
            sumspace=aw+aspace0
            print('c-axis length=',sumspace)
            for j in range(nlayer):
                print(f'{lpos[j][0]:.6f} {lpos[j][1]:.6f} {aww[j]/sumspace:.6f} {atoms[j]}   !',layers[j])
            #POSCAR        
            csa1,csa2,csa3,csa4=[cspec(aa) for aa in [a1,a2,a3,a4]]
            poscar_filename = f'POSCAR.{csa1}{csa2}{h1}{csa3}{csa4}{h2}_{direction}'
            with open(poscar_filename, 'w') as fpos:
                fpos.write(f'{csa1}{csa2}{h1}{csa3}{csa4}{h2}_{direction}\n')
                fpos.write('1.0\n')
                fpos.write(' '.join(f'{x:.6f}' for x in (AA*ain)) + '\n')
                fpos.write(' '.join(f'{x:.6f}' for x in (BB*ain)) + '\n')
                fpos.write(' '.join(f'{x:.6f}' for x in (CC*sumspace)) + '\n')
                atomline={a1:[],a2:[],a3:[],a4:[]}
                for j in range(nlayer):
                    atomline[atoms[j]].append(f'{lpos[j][0]:.6f} {lpos[j][1]:.6f} {aww[j]/sumspace:.6f} {cspec(atoms[j])}')
                atomlist=''
                atomnum=''
                for key,value in atomline.items():
                    if(len(value)!=0):
                        atomlist=atomlist + f'{cspec(key)} '
                        atomnum =atomnum  + f'{len(value)} '
                fpos.write(atomlist+'\n')
                fpos.write(atomnum+'\n')
                fpos.write('direct\n')
                for atom in atomline:
                    for line in atomline[atom]:
                        fpos.write(line+'\n')    
            print(f'--- POSCAR written to {poscar_filename} ---\n')
            
        if direction=='110':
            centeranion=(center+anion)%1
            print(f' === Two layer per C periodic cell total layer= {h1+h2}')
            nlayer=h1+h2
            print(f' === two atoms per layer. layer is specified by z coordinate')
            print(f'    {fmt(zero)} {a1} layer0')
            print(f'    {fmt(centeranion)} {a2} layer0')
            print(f'    {fmt(center)} {a1} layer1')
            print(f'    {fmt(anion)} {a2} layer1')
            print(f'    nlayer={nlayer}')
            print('-------------------------------------------------------')

            lbond = [None] * nlayer
            layers= [None] * nlayer
            lpos  = [None] * nlayer
            atoms = [None] * nlayer

            aa0=a1
            aa1=a2
            ic=0 
            ix=0
            i=0
            layerindex=-1
            zz1=np.array([0,0,1])
            app0=a3
            app1=a4
            while True: #note two atoms per layer. For h1=h2=4 (16 atoms), we have 8 layers \times a atoms. layer is specified by z coordinate.
                i+=1
                if(i>h1):
                    aa0=a3
                    aa1=a4
                    ix=1

                layerindex+=1
                lbond [layerindex]=[app0+aa1,aa0+app1] # we have two bonds between two layers
                lpos  [layerindex]=[zero+ic*zz1,centeranion+ic*zz1]
                atoms [layerindex]=[aa0,aa1] #two atoms per layer
                # we have two atoms per layer.
                layers[layerindex]= [f' {layerindex} {ix} layer0 {fmt(lpos[layerindex][0])} {aa0}',
                                    f' {layerindex} {ix} layer0 {fmt(lpos[layerindex][1])} {aa1} : interlayerbond {layerindex-1}-{layerindex} {lbond[layerindex]}']
                app0=aa0
                app1=aa1
                
                i+=1 
                if(i>h1):
                    aa0=a3
                    aa1=a4
                    ix=1
                    
                layerindex+=1
                lbond [layerindex]=[app0+aa1,     aa0+app1]     #we have two bonds between two layers
                lpos  [layerindex]=[center+ic*zz1,anion+ic*zz1]
                atoms [layerindex]=[aa0,aa1]  #two atoms per layer
                # we have two atoms per layer.
                layers[layerindex]= [f' {layerindex} {ix} layer1 {fmt(lpos[layerindex][0])} {aa0}',
                                    f' {layerindex} {ix} layer1 {fmt(lpos[layerindex][1])} {aa1}: interlayerbond {layerindex-1}-{layerindex} {lbond[layerindex]}']
                app0=aa0
                app1=aa1
                
                if(i>h1+h2-2):
                    break
                ic+=1
        
        # We assume layer spacing is 0.25a van de Walle 1989
            aden = (sum(                         data[lbond[j][0]]['G110'] * (data[lbond[j][0]]['a'] * 0.25) for j in range(nlayer))
                +   sum(                         data[lbond[j][1]]['G110'] * (data[lbond[j][1]]['a'] * 0.25) for j in range(nlayer)) )/2
            anum = (sum( data[lbond[j][0]]['a']* data[lbond[j][0]]['G110'] * (data[lbond[j][0]]['a'] * 0.25) for j in range(nlayer))
                +   sum( data[lbond[j][1]]['a']* data[lbond[j][1]]['G110'] * (data[lbond[j][1]]['a'] * 0.25) for j in range(nlayer)) )/2
            ain  = anum/aden  # Eq.(1a) van de Walle
            print('ainplane=',ain)
            aw=0
            aww= np.zeros(nlayer)
            for j in range(nlayer):
                aav= (data[lbond[j][0]]['a']   +data[lbond[j][1]]['a'])/2
                Dav= (data[lbond[j][0]]['D110']+data[lbond[j][1]]['D110'])/2
                aspace=(0.5*aav)*( 1  -  Dav  * (ain/aav  -1) ) #Eq.(2a) van de Walle
                if(j==0): 
                    aspace0=aspace
                else:
                    aw=aw+aspace
                aww[j]=aw  
            sumspace=aw+aspace0
            print('c-axis length=',sumspace)
            for j in range(nlayer):
                print(f'{lpos[j][0][0]:.6f} {lpos[j][0][1]:.6f} {aww[j]/sumspace:.6f} {atoms[j][0]} ! {layers[j][0]}' )
                print(f'{lpos[j][1][0]:.6f} {lpos[j][1][1]:.6f} {aww[j]/sumspace:.6f} {atoms[j][1]} ! {layers[j][1]}' )
        #POSCAR
            csa1,csa2,csa3,csa4=[cspec(aa) for aa in [a1,a2,a3,a4]]
            poscar_filename = f'POSCAR.{csa1}{csa2}{h1}{csa3}{csa4}{h2}_{direction}'
            with open(poscar_filename, 'w') as fpos:
                fpos.write(f'{csa1}{csa2}{h1}{csa3}{csa4}{h2}_{direction}\n')
                fpos.write('1.0\n')
                fpos.write(' '.join(f'{x:.6f}' for x in (AA*ain)) + '\n')
                fpos.write(' '.join(f'{x:.6f}' for x in (BB*ain)) + '\n')
                fpos.write(' '.join(f'{x:.6f}' for x in (CC*sumspace)) + '\n')
                atomline={a1:[],a2:[],a3:[],a4:[]}
                #print(nlayer,atomline)
                for j in range(nlayer):
                    #print(atoms[j][0],atoms[j][1])
                    atomline[atoms[j][0]].append(f'{lpos[j][0][0]:.6f} {lpos[j][0][1]:.6f} {aww[j]/sumspace:.6f} {cspec(atoms[j][0])}')
                    atomline[atoms[j][1]].append(f'{lpos[j][1][0]:.6f} {lpos[j][1][1]:.6f} {aww[j]/sumspace:.6f} {cspec(atoms[j][1])}')
                atomlist=''
                atomnum=''
                #print('aaaaaaa',atomline)
                for key,value in atomline.items():
                    #print(key,value)
                    if(len(value)!=0):
                        atomlist=atomlist + f'{cspec(key)} '
                        atomnum =atomnum  + f'{len(value)} '
                fpos.write(atomlist+'\n')
                fpos.write(atomnum+'\n')
                fpos.write('direct\n')
                for atom in atomline:
                    #print('atom=',atom)
                    for line in atomline[atom]:
                        fpos.write(line+'\n')    
            print(f'--- POSCAR written to {poscar_filename} ---\n')    

if __name__ == "__main__":
    import sys
    h1=3
    h2=3
    ccaa=['GaAs','AlAs','InAs','GaPp','AlPp','InPp','GaSb','AlSb','InSb',  'ZnSe','ZnSs','ZnTe','CdTe','HgTe']
    i=0
    for ca1 in ccaa:
        for ca2 in ccaa:
            if(ca1==ca2): continue
            c1=ca1[0:2]
            c2=ca2[0:2]
            a1=ca1[2:4]
            a2=ca2[2:4]
            try:
                temp=data[ca1]['a']
                temp=data[ca2]['a']
                temp=data[c1+a2]['a']
                temp=data[c2+a1]['a']
                i+=1
                print(f'Superlattice: {i} {ca1} {ca2}')
                superlattice(ca1,ca2, h1,h2)
            except:
                continue
    import sys
    sys.exit()

            
