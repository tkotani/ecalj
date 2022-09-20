import pymatgen.core as mg
from pymatgen.ext.matproj import MPRester
import os,sys,subprocess
import editglt

def createplot(mpid,key,ncore):
    num = 'mp-'+ mpid
    print('num=',num)
    API_KEY = key
    with MPRester(API_KEY) as mpr:
        material = mpr.get_data(num)
    if len(material) == 0:
        print("APIdate is empty")
        print()
        return ''
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

        subprocess.run(["mpirun","-np",ncore,"lmf-MPIK",num],stdout=f)
        print('lmf-MPIK finished')

        outc=subprocess.check_output(["tail -n1 save."+num],shell=True)
        print(outc.decode('utf-8'),end='')

        aaa="getsyml "+num+' -nobzview > lgetsyml '
        print(aaa)
        os.system(aaa)

        aaa="job_band "+num+" -np "+ncore+" NoGnuplot > ljobband "
        print(aaa)
        os.system(aaa)
        aaa="job_tdos "+num+" -np "+ncore+" NoGnuplot > ljobtdos"
        print(aaa)
        os.system(aaa)

        for iglt in "bandplot.isp1.glt","bandplot.isp2.glt","tdos."+num+".glt":
            if os.path.isfile(iglt) == True:
                os.system('grep -v x11    '+ iglt+' >temp0.glt')
                os.system('grep -v replot temp0.glt >temp.glt')
                os.system("gnuplot temp.glt")
        os.chdir("../")
        return outc
