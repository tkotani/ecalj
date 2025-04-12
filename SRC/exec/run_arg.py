# run_arg.py
import subprocess
import sys

def run_arg(argin, mpi_size, nfpgw, command, output, *target):
    echo_run = True  # standard
    mpi_run = f"srun -n {mpi_size}"  # standard

    target_str = ' '.join(target)
    command_str = f"{nfpgw}{command} {target_str}"

    if echo_run:
        #改行しない
        print(f"OK! --> Start", end=' ')
        #print(f"{argin} > _IN_")

    with open('_IN_', 'w') as f:
        f.write(argin)

    if mpi_size == '0':
        if echo_run:
            print(f"echo {argin} | {command_str} > {output}")
        result = subprocess.run(f"{command_str} < _IN_ > {output}", shell=True)
    else:
        if echo_run:
            print(f"echo {argin} | {mpi_run} {command_str} > {output}")
        result = subprocess.run(f"{mpi_run} {command_str} < _IN_ > {output}", shell=True)

    if result.returncode != 0:
        if echo_run:
            print(f"Error in {command} input_arg={argin}. See OutputFile={output}")
        sys.exit(10)

    #if echo_run:
    #    print(f"NOTE: Use run_arg defined in {nfpgw}/run_arg")
