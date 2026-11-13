
# Installation on kugui and ohtaka (GPU) at ISSP

> English version of the page [install/installISSP](../../install/installISSP.md). The pages at the top of ecaljdoc are the reference; this one was made from install/installISSP.md of the tree of c9f9a3003 (source blob 6ad7f826aed2).

> 📝 **Note** — This is a memo to install ecalj to the machines in
> ISSP@u-tokyo as of 2025. Module names, software versions, and
> filesystem layouts on `kugui` / `ohtaka` change over time; if any
> step fails, check the latest ISSP documentation
> ([SCC page](https://www.issp.u-tokyo.ac.jp/supercom/)) and adjust
> module versions / paths accordingly. The configuration below has
> been verified to work for the author's typical workflow but is not
> automatically updated when ISSP's environment changes.

## Installation on kugui

### 1. Log in

### 2. [Install python and the tools with pyenv or similar](./install.md)

### 3. Edit `~/.bashrc`

Insert the following at the end of `~/.bashrc`.

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
> An example script that loads the required `module`s and starts MPS.
> OneAPI_MPI is slightly faster but slightly less stable.
> Intel MKL can be used with both nvfortran and gfortran.


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
>

### 4. Install and test ecalj

```bash
git clone https://github.com/tkotani/ecalj.git
cd ecalj
```

Then run the install + test command.

```bash
qsub jobinstall_kugui.sh
```

This compiles and then runs the test calculations. `jobinstall_kugui.sh` installs both the CPU and the GPU versions. The
InstallAll.py called in it creates the directory `$HOME/bin` and copies the ecalj binaries and scripts there.
If this is a problem, back up the current copy and let ecalj overwrite it, or search for --bidir bin in InstallAll.py and rewrite it. bin must be in your PATH.

The test proceeds printing OK. Error-like messages such as the following can be ignored (normal operation).

```text
[cpu121:54969] 7 more processes have sent help message help-mpi-common-cuda.txt / dlopen failed
[cpu121:54969] Set MCA parameter "orte_base_help_aggregate" to 0 to see all help / error messages
--------------------------------------------------------------------------
The library attempted to open the following supporting CUDA libraries,
but each of them failed.  CUDA-aware support is disabled.
libcuda.so.1: cannot open shared object file: No such file or directory
libcuda.dylib: cannot open shared object file: No such file or directory
/usr/lib64/libcuda.so.1: cannot open shared object file: No such file or directory
/usr/lib64/libcuda.dylib: cannot open shared object file: No such file or directory
If you are not interested in CUDA-aware support, then run with
--mca opal_warn_on_missing_libcuda 0 to suppress this message.  If you are interested
in CUDA-aware support, then try setting LD_LIBRARY_PATH to the location
of libcuda.so.1 to get passed this issue.
```

The test ends like

```text
PASSED! ni_crpa/Screening_W-v_crpa.h
PASSED! srvo3_crpa/Screening_W-v.h
PASSED! srvo3_crpa/Screening_W-v_crpa.h
OK! ALL PASSED ===
  See work/summary.txt
Elapsed time for make        : 51 seconds
Elapsed time for testecalj.py: 409 seconds
```

This completes the installation and the test (about 10 minutes in total).

During the test,

```bash
tail -f ecalj/Samples/TestInstall/summary.txt
```

shows the progress. It is OK if it ends with `PASSED! srvo3_crpa/Screening_W-v_crpa.h` and there is no FAIL or error.
Even with nvfortran the test runs on the CPU.


### 5. GPU calculation test

Copy `ecalj/Samples/BenchmarkTest/inas2gasb2` to a working directory (the working files are large), put
`ecalj/Samples/BenchmarkTest/ISSP/job_kugui.sh` there, and run

```bash
qsub job_kugui.sh
```

. Check with `cat lgwsc`; if it ends with

```
===== QSGW iteration end   iter 2 ===
OK! ==== All calclation finished for  gwsc ====
```

it is OK.


### 6. Usage memo

[Usage memo](../manual/UsageISSP.md)


## Installation on ohtaka

Almost the same. Add the following at the end of `~/.bashrc`.

```bash
ulimit -s unlimited
export PATH=$HOME/bin:$PATH
module purge
module load openmpi/4.1.5-oneapi-2023.0.0-classic
```

The installation is

```bash
sbatch jobinstall_ohtaka.sh
```

. If you do not see `OK! ALL PASSED!` at the end at this stage, the installation has failed.

After confirming it, for the calculation test copy `ecalj/Samples/BenchmarkTest/inas2gasb2` to a working directory, put
`ecalj/Samples/BenchmarkTest/ISSP/job_ohtaka.sh` there, and

```bash
sbatch job_ohtaka.sh
```

run it. The 16-atom `inas4gasb4` follows the same procedure (rewrite `id` in the script).
The scripts were adapted on 2026-09-30 to the current input (`ctrlg.<sname>.toml`) and options, and have not been run on the ISSP machines since.


On kugui, srun is used instead of mpirun, among other things, so this is needed when running gwsc and others.

## Customization for the ISSP systems

The machine dependence is described in `SRC/exec/MachineDependence.py` of the ecalj source code.
On kugui, srun is used instead of mpirun, among other things; this is needed when running gwsc and others.
