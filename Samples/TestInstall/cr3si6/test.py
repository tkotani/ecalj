import os
from comp import test1_check,test2_check,runprogs
def test(args,bindir,testdir,workdir):
        lmfa= f'mpirun -np 1 {bindir}/lmfa '
        lmf = f'mpirun -np {args.np} {bindir}/lmf '
        message1='''
        # Case cr3si6: a hexagonal environment with several atoms, two classes
        # --- Test 1.  Basic check of programs lmfa,lmf ---
         Case cr3si6: a hexagonal environment with several atoms, two classes
           Other checks:       verbose output, insulator, forces(mode 1)
         --- Test 1.  Basic check of programs lmfa,lmf ---
         Checks that program lmfa produces a sensible atm file
         and that program lmf iterates to the proper energy.
         The cr3si6 test also checks and illustrates the following:
         1.  insulator
         2.  forces with correction mode 1 (FORCES=1)
         3.  verbose output
         4.  A case with low kmxa (kmxa=2)
         '''       
        print(message1)
        outfile='out.lmf.cr3si6'
        runprogs([
                 lmfa+" cr3si6 --pr51 -vnit=2 --time=6 >  "+outfile,
                 lmf +" cr3si6 --pr51 -vnit=2 --time=6 >> "+outfile
        ])
        tall=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
        return tall
      