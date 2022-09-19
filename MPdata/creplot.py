import pymatgen.core as mg
from pymatgen.ext.matproj import MPRester
import os,sys,subprocess
import editglt

#from mp_api import MPRester
#from mp_api.matproj import MPRester

def createplot(mpid,key,fout):
    #num=['mp',mpid]
    #aa = mpid.split("-")
    #print(aa)
    num = 'mp-'+ mpid
    print()
    print('num=',num)
#    print('',file=fout)
#    print('num=',num,file=fout)
    API_KEY = key
    #print(API_KEY)
    with MPRester(API_KEY) as mpr:
        material = mpr.get_data(num)
    if len(material) == 0:
        print("APIdate is empty")
        print("APIdate is empty",file=fout)
        return

    #print(material) #many date of material, dic type
    #sys.exit()
    else:
 
        for item in material:
            mate = item["pretty_formula"]
            os.makedirs(num,exist_ok=True)
            os.chdir(num)
            #print(item['material_id'])    # result is 'mp-101'
            struc = mpr.get_structure_by_material_id(item['material_id'])
            struc.to(fmt='poscar', filename='POSCAR')

        subprocess.run("pwd")
        f=open('llmf','w')
        subprocess.run(["vasp2ctrl POSCAR"],shell=True,stdout=f)
        subprocess.run(["cp ctrls.POSCAR.vasp2ctrl ctrls."+num],shell=True,stdout=f)
        subprocess.run(["ctrlgenM1.py "+num],shell=True,stdout=f,stderr=f)
        subprocess.run(["cp ctrlgenM1.ctrl."+num+" ctrl."+num],shell=True,stdout=f)

        subprocess.run(["lmchk "+num+"|grep conf"],shell=True)        
        subprocess.run(["lmfa",num],stdout=f)
        print('lmfa finished')
        subprocess.run(["mpirun","-np","4","lmf-MPIK",num],stdout=f)
        print('lmf-MPIK finished')
        #subprocess.run(["tail -n1 save."+num],shell=True)
        outc=subprocess.check_output(["tail -n1 save."+num],shell=True)
        print(outc.decode('utf-8'),end='')
        subprocess.run(["getsyml",num,'-nobzview'],stdout=f,stderr=f)
        subprocess.run(["job_band "+num+" -np 4 NoGnuplot"],stdout=f,shell=True)
        subprocess.run(["job_tdos "+num+" -np 4 NoGnuplot"],stdout=f,shell=True)
        print('getsyml job_band job_tdos finished')

        #os.chdir("../")
        #subprocess.run("pwd")
        #subprocess.run(["editglt.py",mate,num])
        #editglt.edit(mate,num)
        #os.chdir(mate)
        for iglt in "bandplot.isp1.glt","bandplot.isp2.glt","tdos."+num+".glt":
            if os.path.isfile(iglt) == True:
                subprocess.run(["gnuplot",iglt])
        os.chdir("../")
        return outc
