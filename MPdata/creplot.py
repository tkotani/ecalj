import pymatgen.core as mg
from pymatgen.ext.matproj import MPRester
import os,sys,subprocess

def createplot(mpid,key,ncore,lmxa6):
    num = 'mp-'+ mpid
    print('num=',num)
    if True:
        os.system('cp POSCARALL/POSCAR.'+num+' POSCAR')
        subprocess.run("pwd")
        f=open('llmf','w')
        subprocess.run(["vasp2ctrl POSCAR"],shell=True,stdout=f)
        subprocess.run(["cp ctrls.POSCAR.vasp2ctrl ctrls."+num],shell=True,stdout=f)
        subprocess.run(["ctrlgenM1.py "+num],shell=True,stdout=f,stderr=f)
        subprocess.run(["cp ctrlgenM1.ctrl."+num+" ctrl."+num],shell=True,stdout=f)
        if lmxa6:
            os.system("sed -e 's/LMXA=4/LMXA=6/g' ctrl."+num+'> temp')
            os.system('cp temp ctrl.'+num)
        subprocess.run(["lmchk "+num+"|grep conf"],shell=True)        
        subprocess.run(["lmfa",num],stdout=f)
        print('lmfa finished')

        subprocess.run(["mpirun","-np",ncore,"lmf-MPIK",num],stdout=f)
        print('lmf-MPIK finished')

        aaa="tail -n1 save."+num+" >savet"
        os.system(aaa)
        with open('savet') as f:
            outc=f.read()
            print(outc,end='')
        #outc=subprocess.check_output(["tail -n1 save."+num],shell=True)
        #print(outc.decode('utf-8'),end='')

        aaa="getsyml "+num+' -nobzview > lgetsyml '
        print(aaa)
        os.system(aaa)
        aaa="job_band "+num+" -np "+ncore+" NoGnuplot > ljobband "
        print(aaa)
        os.system(aaa)
        aaa="job_tdos "+num+" -np "+ncore+" NoGnuplot > ljobtdos"
        print(aaa)
        os.system(aaa)
        aaa="job_pdos "+num+" -np "+ncore+" NoGnuplot > ljobpdos"
        print(aaa)
        os.system(aaa)
        
        for iglt in os.listdir('.'):
            files = os.listdir('.')
            if not ('.glt' in iglt): continue
            os.system('grep -v x11    '+ iglt+' >temp0.gxx')
            os.system('grep -v replot temp0.gxx >temp.gxx')
            os.system("gnuplot temp.gxx")
        os.chdir("../")
        return outc
