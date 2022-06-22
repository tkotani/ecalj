import pymatgen.core as mg
from pymatgen.ext.matproj import MPRester
import os,sys,subprocess
import editglt

#from mp_api import MPRester
#from mp_api.matproj import MPRester

def createplot(mpid,key):
   
    #num=['mp',mpid]

    #aa = mpid.split("-")
    #print(aa)
    num = 'mp-'+ mpid
    print('num=',num)
    
    API_KEY = key
    print(API_KEY)

    with MPRester(API_KEY) as mpr:
        material = mpr.get_data(num)

    if len(material) == 0:
        print("APIdate is empty")
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

        subprocess.run(["vasp2ctrl POSCAR"],shell=True)
        subprocess.run(["cp ctrls.POSCAR.vasp2ctrl ctrls."+num],shell=True)
        subprocess.run(["ctrlgenM1.py "+num],shell=True)
        subprocess.run(["cp ctrlgenM1.ctrl."+num+" ctrl."+num],shell=True)
        subprocess.run(["lmfa",num])
        subprocess.run(["mpirun","-np","4","lmf-MPIK",num])
        subprocess.run(["getsyml",num,'-nobzview'],stdout=f)
        subprocess.run(["job_band "+num+" -np 4 NoGnuplot"],shell=True)
        subprocess.run(["job_tdos "+num+" -np 4 NoGnuplot"],shell=True)

        #os.chdir("../")
        #subprocess.run("pwd")
        #subprocess.run(["editglt.py",mate,num])
        #editglt.edit(mate,num)
        #os.chdir(mate)
        for iglt in "bandplot.isp1.glt","bandplot.isp2.glt","tdos."+num+".glt":
            if os.path.isfile(iglt) == True:
                subprocess.run(["gnuplot",iglt])
        os.chdir("../")
        return 
