
# Installation on <del class="clfix">kugui and ohtaka (GPU)</del><span class="clfix">kugui (GPU) and ohtaka (CPU)</span> at ISSP

> English version of the page [install/installISSP](../../install/installISSP.md). The pages at the top of ecaljdoc are the reference; this one was made from install/installISSP.md of the tree of c9f9a3003 (source blob 42fac9413b5c).

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


<span class="clfix">The following message can appear at login or with mpirun; it is harmless (see §4).</span>

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

<div class="clfix-old">

```bash
qsub jobinstall_kugui.sh
```

</div>

<div class="clfix">

`jobinstall_kugui.sh` was removed from the repository (2026-10-01), so write a job script with the following content (e.g. `jobinstall_kugui.sh`) in `ecalj` and submit it. The MPI tests cannot be run on the front end.

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

</div>

This compiles and then runs the test calculations. `jobinstall_kugui.sh` installs both the CPU and the GPU versions. The
InstallAll.py called in it creates the directory `$HOME/bin` and copies the ecalj binaries and scripts there.
If this is a problem, back up the current copy and let ecalj overwrite it, or <del class="clfix">search for --bidir bin in InstallAll.py and rewrite it. bin must be in your PATH.</del><span class="clfix">pass `--bindir <dir>` to InstallAll.py (default `$HOME/bin`). That directory must be in your PATH.</span>

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

<del class="clfix">`tail -f ecalj/Samples/TestInstall/summary.txt`</del>

<span class="clfix">the standard output of the job (with PBS, `install.o<job id>`)</span> shows the progress. <span class="clfix">The result of each target is written to `ecalj/Samples/TestInstall/<target>_work/summary.txt` (the `See work/summary.txt` in the output example above is from an older version).</span> It is OK if it ends with `PASSED! srvo3_crpa/Screening_W-v_crpa.h` and there is no FAIL or error.
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

it is OK.<span class="clfix"> The sample contains the earlier result `QPU.1run`, so the one iteration of `gwsc 1` is counted as number 2 and its result is written to `QPU.2run`. Compare it with `QPU.1run` of the sample. The script passes `-c kugui_gpu` (the GPU-node settings in `SRC/exec/clusters.toml`; added 2026-11-14).</span>


### 6. Production runs

<div class="clfix">

(Moved from the UsageISSP page, 2026-11-14.)

The queue for production runs is `F1accs` etc. (`i1accs` is the test queue, up to 30 minutes). Example:

```bash job.sh
#!/bin/sh
#PBS -q F1accs
#PBS -l select=1:ncpus=64:mpiprocs=64:ompthreads=1
#
ulimit -s unlimited
id=inas6gasb6
gwsc 5 -np 64 -np2 4 --gpu -c kugui_gpu $id > lgwsc
```

- `-np` is the number of MPI processes of the parts run on the CPU, `-np2` that of the parts run on the GPU (usually the number of GPUs). The first number (5) is the number of iterations. See [GPU version](../manual/ecaljgpu).
- The GPU version is not so effective for small systems, but four GPUs can handle about 40 atoms (it depends on the conditions; the cost grows with the square of the number of k points).

</div>


## Installation on ohtaka

Almost the same. Add the following at the end of `~/.bashrc`.

```bash
ulimit -s unlimited
export PATH=$HOME/bin:$PATH
module purge
module load openmpi/4.1.5-oneapi-2023.0.0-classic
```

The installation is

<div class="clfix-old">

```bash
sbatch jobinstall_ohtaka.sh
```

</div>

<div class="clfix">

`jobinstall_ohtaka.sh` was also removed from the repository, so write a script with the following content in `ecalj` and submit it.

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

</div>

. If you do not see `OK! ALL PASSED!` at the end at this stage, the installation has failed.

After confirming it, for the calculation test copy `ecalj/Samples/BenchmarkTest/inas2gasb2` to a working directory, put
`ecalj/Samples/BenchmarkTest/ISSP/job_ohtaka.sh` there, and

```bash
sbatch job_ohtaka.sh
```

run it.<span class="clfix"> Here too the result is written to `QPU.2run`; compare it with `QPU.1run` of the sample. The script passes `-c ohtaka` (launch with srun; added 2026-11-14).</span> The 16-atom `inas4gasb4` follows the same procedure (rewrite `id` in the script).
The scripts were adapted on 2026-09-30 to the current input (`ctrlg.<sname>.toml`) and options, and have not been run on the ISSP machines since.


<del class="clfix">On kugui, srun is used instead of mpirun, among other things, so this is needed when running gwsc and others.</del><span class="clfix">(removed: it repeats the next section)</span>

<div class="clfix">

(Moved from the UsageISSP page, 2026-11-14.)

- `InstallAll.py` uses `mpirun` inside; it works because it runs on one node.
- An example of the parallel settings of a large calculation (16 nodes, 512 processes). If memory runs short, reduce the cores per node to get more memory per core.

```bash
#SBATCH -p F16cpu
#SBATCH -N 16
#SBATCH -n 512
#SBATCH -c 4
#SBATCH --ntasks-per-node=32
#SBATCH -t 24:00:00
```

- If MKL errors appear during the run, set `-c 1`.

</div>

## Customization for the ISSP systems

The machine dependence is described in <del class="clfix">`SRC/exec/MachineDependence.py`</del><span class="clfix">`SRC/exec/clusters.toml` (selected with `-c <name>` of `gwsc` and others) and `SRC/exec/pylib/MachineDependence.py`</span> of the ecalj source code.
<del class="clfix">On kugui, srun is used instead of mpirun, among other things</del><span class="clfix">On ohtaka srun is used instead of mpirun, and on the kugui GPU nodes `mpirun --bind-to none`, among other things</span>; this is needed when running gwsc and others.
