
# Installation on kugui (GPU) and ohtaka (CPU) at ISSP

> English version of the page [install/installISSP](../../install/installISSP.md). The pages at the top of ecaljdoc are the reference; this one was made from install/installISSP.md of commit 296542c3a (source blob ba9bd6f61d8f).

> 📝 **Note** — This is a memo, as of 2025, for installing ecalj on the computers of the Institute for Solid State Physics (ISSP), the University of Tokyo.
> The module names, software versions, and file system layout of `kugui` / `ohtaka` change over time.
> If a step does not work, see the latest ISSP documentation ([SCC page](https://www.issp.u-tokyo.ac.jp/supercom/)) and
> adjust the module versions and paths. The settings below were confirmed to work in the author's usual use,
> but they are not updated automatically when the ISSP environment changes.

## Installation on kugui

### Log in

### [Install Python and tools with pyenv etc.](./install.md)

### Editing `~/.bashrc`

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
> This is an example script that loads the required `module`s and starts MPS.
> OneAPI_MPI is slightly faster but slightly less stable.
> Intel MKL can be used with both nvfortran and gfortran.

The following message may appear at login or at mpirun, but it is harmless (§1.4).

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

### Installing and testing ecalj

```bash
git clone https://github.com/tkotani/ecalj.git
cd ecalj
```

Then run the install + test command.

`jobinstall_kugui.sh` was removed from the repository (2026-10-01), so create a job script with the following content (e.g. `jobinstall_kugui.sh`) in `ecalj` and submit it. This is because the MPI tests cannot be run on the front end.

```bash
#!/bin/sh
#PBS -q i2cpu
#PBS -l select=1:ncpus=16:mpiprocs=8:ompthreads=1
#PBS -N install
ulimit -s unlimited
module purge
module load nvhpc-nompi/24.7 openmpi_nvhpc compiler-rt tbb mkl
export PATH=$HOME/bin:$PATH
./InstallAll.py --fc nvfortran --gpu
```

```bash
qsub jobinstall_kugui.sh
```

When you run this, the test calculations proceed after compilation. `jobinstall_kugui.sh` installs both the CPU version and the GPU version. The
InstallAll.py called in it creates the directory `$HOME/bin` and copies the ecalj binaries and scripts there.
If this is a problem, back up the current copy and let ecalj overwrite it, or pass `--bindir <directory>` to InstallAll.py (default `$HOME/bin`). The directory must be in your path.

The tests proceed, printing OK. You can ignore error-like messages such as the following (this is normal operation).

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

The tests end with something like

```text
PASSED! ni_crpa/Screening_W-v_crpa.h
PASSED! srvo3_crpa/Screening_W-v.h
PASSED! srvo3_crpa/Screening_W-v_crpa.h
OK! ALL PASSED ===
  See work/summary.txt
Elapsed time for make        : 51 seconds
Elapsed time for testecalj.py: 409 seconds
```

This completes the installation and tests (about 10 minutes in total).

During the tests,

you can check the progress in the standard output of the job (with PBS, `install.o<job number>`). The result of each target is written to `ecalj/Samples/TestInstall/<target>_work/summary.txt` (`See work/summary.txt` in the output example above is what an old version printed). It is OK if `PASSED! srvo3_crpa/Screening_W-v_crpa.h` is shown at the end and there is no FAIL or error.
Even with nvfortran, the tests run on the CPU.

### GPU calculation test

Copy `ecalj/Samples/BenchmarkTest/inas2gasb2` to a working directory (the work files are large), put
`ecalj/Samples/BenchmarkTest/ISSP/job_kugui.sh` there, and

```bash
qsub job_kugui.sh
```

run it. Do `cat lgwsc`; it is OK if it ends with

```
===== QSGW iteration end   iter 2 ===
OK! ==== All calclation finished for  gwsc ====
```

Since the sample contains the previous result `QPU.1run`, the single iteration of `gwsc 1` is counted as number 2, and the result is written to `QPU.2run`. Compare it with `QPU.1run` of the sample. The script passes `-c kugui_gpu` (the GPU node setting in `SRC/exec/clusters.toml`).

### Production calculations

The queues for production calculations are `F1accs` and others (`i1accs` is for tests of up to 30 minutes). Example:

```bash job.sh
#!/bin/sh
#PBS -q F1accs
#PBS -l select=1:ncpus=64:mpiprocs=64:ompthreads=1
#
ulimit -s unlimited
id=inas6gasb6
gwsc 5 -np 64 -np2 4 --gpu -c kugui_gpu $id > lgwsc
```

- `-np` is the number of MPI processes for the parts running on the CPU, and `-np2` is the number of processes for the parts running on the GPU (usually the number of GPUs). The first number (5) is the number of iterations. For details see [GPU version](../manual/ecaljgpu).
- The GPU version is not very effective for small systems, but with 4 GPUs it can handle up to about 40 atoms (depending on the calculation conditions; the cost grows as the square of the number of k points).

## Installation on ohtaka

It is almost the same. Add the following at the end of `~/.bashrc`.

```bash
ulimit -s unlimited
export PATH=$HOME/bin:$PATH
module purge
module load openmpi/4.1.5-oneapi-2023.0.0-classic
```

`jobinstall_ohtaka.sh` was also removed from the repository, so create a script with the following content in `ecalj` and submit it.

```bash
#!/bin/sh
#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -c 1
#SBATCH --job-name=install
#SBATCH --ntasks-per-node=8
./InstallAll.py --fc ifort
```

```bash
sbatch jobinstall_ohtaka.sh
```

If you do not see `OK! ALL PASSED!` at the end at this stage, the installation has failed.

After confirming this, for the calculation test copy `ecalj/Samples/BenchmarkTest/inas2gasb2` to a working directory, put
`ecalj/Samples/BenchmarkTest/ISSP/job_ohtaka.sh` there, and

```bash
sbatch job_ohtaka.sh
```

run it. Here too the result is written to `QPU.2run`, so compare it with `QPU.1run` of the sample. The script passes `-c ohtaka` (launch with srun). The 16-atom `inas4gasb4` uses the same procedure (edit `id` in the script).
The scripts were adapted to the current input (`ctrlg.<sname>.toml`) and options on 2026-09-30, and have not been run on the ISSP machines since then.

- `InstallAll.py` uses `mpirun` internally, but it works because it is a single node.
- An example of the parallel settings for a large calculation (16 nodes, 512 processes). If memory is insufficient, reduce the cores per node to increase the memory per core.

```bash
#SBATCH -p F16cpu
#SBATCH -N 16
#SBATCH -n 512
#SBATCH -c 4
#SBATCH --ntasks-per-node=32
#SBATCH -t 24:00:00
```

- If an MKL error appears during the run, use `-c 1`.

## Customization for the ISSP systems

The machine dependence is described in `SRC/exec/clusters.toml` (selected with `-c <name>` of `gwsc` etc.) and `SRC/exec/pylib/MachineDependence.py` of the ecalj source code.
On ohtaka srun is used instead of mpirun, and on the GPU nodes of kugui `mpirun --bind-to none` is used; this is needed when running gwsc etc.
