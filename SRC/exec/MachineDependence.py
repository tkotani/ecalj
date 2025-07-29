### Machine dependence ###################
import platform
machine_info = platform.uname()
#print(machine_info)

def mpiRUN(ncore):
    if("kugui" in machine_info):
        mpirun= f'mpirun --bind-to none --map-by node -np {ncore} ' if ncore!=0 else ''
    elif("ohtaka" in machine_info):
        mpirun= f'srun -n {ncore} ' if ncore!=0 else ''
    else:    
        mpirun= f'mpirun --bind-to core --map-by core -np {ncore} ' if ncore != 0 else ''
    return mpirun
def qsubCOMMAND(jobx):
    if("ohtaka" in machine_info):
        qsubc= f'sbatch {jobx}'
    else:    
        qsubc= f'qsub {jobx}'
    return qsubc
