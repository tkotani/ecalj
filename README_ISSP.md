# For the users of the ISSP supercomputer

## Common Information for System B and C
The default version of Python is outdated. We will prepare the latest Python in local.

### Installing `python` using `mise`
> [!TIP]
> mise is a type of package management software.

1. Add the following settings to `~/.bashrc` for the automatic installation and activation of `mise`:
```bash ~/.bashrc
type mise > /dev/null 2>&1 || curl https://mise.run | sh
eval "$(~/.local/bin/mise activate bash)"
```

2. Update `~/.bashrc` to install `mise`:
```bash
source ~/.bashrc
```

3. Install `python` using `mise`:
```bash
mise use python@latest -g
```

4. Install the required `python` libraries:
```bash
pip install numpy pandas seekpath spglib pymatgen mp-api scipy plotly
```

## Server-Specific Information

## System B: Ohtaka

1. Append the following to the end of `~/.bashrc`:
```bash
ulimit -s unlimited
export PATH=$HOME/bin:$PATH
module purge
module load openmpi/4.1.5-oneapi-2023.0.0-classic  
```

## System C: Kugui

1. Append the following to the end of `~/.bashrc`:
```bash
ulimit -s unlimited
export PATH=$HOME/bin:$PATH
module purge
module load nvhpc-nompi/24.7 openmpi_nvhpc compiler-rt tbb mkl
if which nvidia-cuda-mps-control > /dev/null 2>&1 ; then
  export CUDA_MPS_PIPE_DIRECTORY=$(pwd)/nvidia-mps-$(hostname)
  export CUDA_MPS_LOG_DIRECTORY=$(pwd)/nvidia-log-$(hostname)
  echo "start nvidia-cuda-mps-control at" $(hostname)
  nvidia-cuda-mps-control -d
fi
```
> [!TIP]
> Script for loading the required `modules` and starting MPS.

### about warning output

> [!INFO]
>The following messages may appear on the log file of the calculation, but you can ignore them.
>```
>[cpu121:54969] 7 more processes have sent help message help-mpi-common-cuda.txt / dlopen failed
>[cpu121:54969] Set MCA parameter "orte_base_help_aggregate" to 0 to see all help / error messages
>--------------------------------------------------------------------------
>The library attempted to open the following supporting CUDA libraries,
>but each of them failed.  CUDA-aware support is disabled.
>libcuda.so.1: cannot open shared object file: No such file or directory
>libcuda.dylib: cannot open shared object file: No such file or directory
>/usr/lib64/libcuda.so.1: cannot open shared object file: No such file or directory
>/usr/lib64/libcuda.dylib: cannot open shared object file: No such file or directory
>If you are not interested in CUDA-aware support, then run with
>--mca opal_warn_on_missing_libcuda 0 to suppress this message.  If you are interested
>in CUDA-aware support, then try setting LD_LIBRARY_PATH to the location
>of libcuda.so.1 to get passed this issue.
>```

## Further Details
ecalj's Ohttaka & Kugui branches are already specialized for them.
See [ecaljdoc](https://ecalj.github.io/ecaljdoc/guide/server_config) for the modification.
