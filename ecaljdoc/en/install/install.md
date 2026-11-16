# Installing ecalj {#ecalj-install}

> English version of the page [install/install](../../install/install.md). The pages at the top of ecaljdoc are the reference; this one was made from install/install.md of commit 4f74d9132 (source blob 5ceb8f1c83f8).

> ⚠️ **Migration to TOML (2026-05)** — The Fortran binaries read only `ctrlg.<sname>.toml`. Convert an old `ctrl.<sname>` / `GWinput` with `Legacy2toml.py <sname>`. For the whole guide and the list of migrated samples, see [Input files](../manual/toml_migration).

To install ecalj, you need to install the software and the Python modules that ecalj uses.
In addition, we make some soft links for convenience. Most of the steps are almost automatic. Follow the steps below.



## Required software and native tools {#_1-required-software-and-native-tools}

The core of ecalj is written in modern Fortran (Fortran 90/2003 style). The Fortran source files are in the following places.

```text
ecalj/SRC/main/*.f90
ecalj/SRC/subroutines/*.f90
```

The cpp switches for GPU support (OpenACC/CuBLAS) can be turned on/off automatically at compile time, as described below.

Recommended or required native tools and libraries (an Ubuntu/Debian example is given later):

* `git` (to clone and update the ecalj repository). `gitk` is useful for browsing the history visually.
* Fortran compiler: choose one available on your machine among `gfortran` (GNU), `ifort` (Intel), and `nvfortran` (NVIDIA). `nvfortran` is for GPU.
* BLAS/LAPACK/FFT libraries. Intel MKL is recommended if available. OpenBLAS also works.
* MPI implementation: Open MPI or MPICH. (On Ubuntu 24.04 it worked with Open MPI 4.x.)
* Build tools: `cmake`, `make`, `bash`.
* Plotting and utilities: `gnuplot` (used by plotting scripts and tests; not mandatory).

On Debian/Ubuntu, install these tools with the following commands.

```bash
sudo apt update
sudo apt install -y git build-essential cmake make gfortran libopenmpi-dev openmpi-bin \
  libblas-dev liblapack-dev libfftw3-dev libopenblas-dev gnuplot
```


Notes and hints:

* Intel MKL works well both with the Intel compilers and with gfortran.
<!-- * For GPU builds use `nvfortran` and pass `--gpu` to the installer; CUDA toolkit and drivers must be present. (See 3.Install below) -->
* Nowadays you may need to install such tools under your home directory.


Example: in my case (ThinkPad T14, Ubuntu 24.04) I use the following.

* Open MPI 4.1.6
* CMake 3.28.x
* make 4.3
* gfortran 13.x
* Intel MKL 2020.x (if available)


`git` is useful for updating and for tracking your local changes (e.g. `git diff` in `ecalj/`).



## Python and Python modules (using pyenv) {#_2-python-and-python-modules-using-pyenv}

> [!TIP]
> We recommend managing the Python version used for ecalj with `pyenv`. It leaves the system Python untouched and makes the environment easier to reproduce. You do not need `pyenv-virtualenv` unless you want to manage several environments. If you have only one environment for ecalj, `pyenv` alone is enough.
>
> **Note:** On recent systems such as Ubuntu, it works best to keep the project's Python environment separate from the system Python.

Python 3.9 or later is required. Below is the recommended procedure to install `pyenv`, set up the shell, install a Python version, and install the Python packages required by ecalj.


1. Install the prerequisites (Debian/Ubuntu example).

```bash
sudo apt update
sudo apt install -y build-essential curl libssl-dev zlib1g-dev libbz2-dev \
  libreadline-dev libsqlite3-dev wget llvm libncurses5-dev libncursesw5-dev \
  xz-utils tk-dev libffi-dev liblzma-dev python-openssl git
```

1. Install `pyenv` (recommended installation method).

```bash
curl https://pyenv.run | bash
```

1. Add the `pyenv` initialization to your shell startup file (`~/.bashrc` or `~/.profile`).

```bash
# pyenv
export PATH="$HOME/.pyenv/bin:$PATH"
eval "$(pyenv init --path)"
eval "$(pyenv init -)"
```

Then reload the shell.

```bash
source ~/.bashrc
```

1. Install the Python version used for ecalj. Example: install the latest 3.11 and set it as the local version.

```bash
pyenv install 3.11.6    # choose a concrete, available patch version
pyenv local 3.11.6      # or pyenv global 3.11.6 for system-wide in that account
```

1. Upgrade pip and install the required packages.

```bash
pip install --upgrade pip
pip install numpy pandas seekpath spglib pymatgen mp-api scipy plotly gnuplot cif2cell
```


Notes:
* If you prefer conda/anaconda, you may substitute its procedure. However, `pyenv` leaves the system Python untouched and coexists well with the environments of several projects.
* On some systems it may be better to install the build requirements of Python with the package manager (the example above is for Debian/Ubuntu).


If needed, you can create a `requirements.txt` listing the packages above in your project, to reproduce the same installation.

```text
numpy
pandas
seekpath
spglib
pymatgen
mp-api
scipy
plotly
gnuplot
cif2cell
```


Then install them as follows.

```bash
pip install -r requirements.txt
```



## Install and install test {#_3-install-and-installtest}

**On ohtaka and kugui (GPU machine) at ISSP, skip this section and see [installISSP.md](./installISSP.md).**


### Getting the ecalj package {#get-ecalj-package}

The ecalj package is at [https://github.com/tkotani/ecalj/](https://github.com/tkotani/ecalj/).
Check out the master branch.

For the ISSP machines, see [https://github.com/msobt/ecaljdoc](https://github.com/msobt/ecaljdoc).

To clone ecalj, do the following.

```bash
git clone https://github.com/tkotani/ecalj.git
```


Cloning creates a directory `ecalj/`.
You can view the history of ecalj as follows.

```bash
cd ecalj
gitk --all
```


### InstallAll.py {#installall-py}

Use `InstallAll.py [Options]` to install. The help is shown by the following.

```bash
./InstallAll.py --help
```


Example of the help output:

```text
usage: InstallAll [-h] [-np NP] [--clean] [--gpu] [--bindir BINDIR] --fc FC [--notest] [--verbose]

Install ecalj and run tests.

options:
  -h, --help       show this help message and exit
  -np NP           number of mpi cores for install test
                   default: 8
                   specify the number of MPI parallelization in test calculation
  --clean          Clean CMakeCache CMakeFiles before make
                   default: none
                   delete the cache files before compiling
  --gpu            nvfortran for GPU
                   default: none
                   compile the GPU and GPU-MP version
  --bindir BINDIR  ecalj binaries and scripts
  --fc FC          fortran compiler gfortran/ifort/ifx/nvfortran
  --notest         no test. only compile
  --verbose        verbose on for debug
```


`InstallAll.py` copies all the ecalj binaries and scripts to the directory given by `--bindir <path>/` (default `$HOME/bin/`). Choose `<path>/` as you like. `--fc` is mandatory.
For GPU, use `--fc nvfortran` together with `--gpu`.



The tests are in `ecalj/Samples/TestInstall/`. They are run by the
`testecalj` script under `ecalj/SRC/exec/` (it is installed into
`BINDIR` together with the binaries).


## Installing ecalj {#install-ecalj}

`InstallAll.py` compiles and installs ecalj, and then runs the install tests in
`ecalj/Samples/TestInstall/` (run by `testecalj`, with the comparison helper
`SRC/exec/comp.py`).
On success it prints 'OK! ALL PASSED!' at the end. It takes 3 to 10 minutes.

To install ecalj, first look at the help.

```bash
./InstallAll.py --help
```

A typical install command:

```bash
./InstallAll.py --fc gfortran
```

For a GPU build, do the following.

```bash
./InstallAll.py --fc nvfortran --gpu
```

To start from scratch, add `--clean`.

This compiles and links the Fortran sources, copies all programs and scripts to bin, and runs the minimal install tests in `ecalj/Samples/TestInstall/` (3 to 10 minutes).

*Hint:* For details of where the binaries are installed, see `./InstallAll.py` (default `BINDIR=~/bin`). Tools such as `getsyml`, `viewvesta`, and `vasp2ctrl` become soft links.

When all tests pass, the output looks like the following.

```text
PASSED! copt TEST 1 out.lmf.copt
PASSED! te TEST 1 out.lmf.te
...
PASSED! nio_gwsc /home/<you>/ecalj/Samples/TestInstall/nio_gwsc_work/QPU
PASSED! nio_gwsc /home/<you>/ecalj/Samples/TestInstall/nio_gwsc_work/log.nio
PASSED! fe_gwsc  /home/<you>/ecalj/Samples/TestInstall/fe_gwsc_work/QPU
PASSED! fe_gwsc  /home/<you>/ecalj/Samples/TestInstall/fe_gwsc_work/log.fe
PASSED! ni_crpa  /home/<you>/ecalj/Samples/TestInstall/ni_crpa_work/Screening_W-v.h
PASSED! srvo3_crpa /home/<you>/ecalj/Samples/TestInstall/srvo3_crpa_work/Screening_W-v.h
OK! ALL PASSED ===

real    4m49.712s
user    19m38.952s
sys     3m38.583s
```

In `ecalj/Samples/TestInstall/`, you can also run an individual test as follows.

```bash
cd ecalj/Samples/TestInstall
testecalj si_gwsc -np 8
```

* Each test has its own directory (e.g. `si_gwsc`).
* `testecalj` creates `si_gwsc_work/` next to it, runs the calculation there,
  and compares the results with the reference data in `ecalj/Samples/TestInstall/si_gwsc/`.
* The procedure of each test is in `test.py` (a small file) in each directory.
  The comparison helper is in `ecalj/SRC/exec/comp.py`. It is easy to add
  your own tests.

If the installation fails, you may need to start over from a clean state. Delete `CMakeFiles` and `CMakeCache.txt` in `ecalj/SRC/exec/`. Old `*.mod` files left under `SRC/` can also cause problems. Usually running `./InstallAll.py --clean` fixes it.

If the parallel make (`make -j`) fails, run make by hand.

```bash
cd ecalj/SRC/exec/
rm -rf CMakeFiles CMakeCache.txt
FC=nvfortran cmake .
make
```


After installation, set the command path (BINDIR). For example, if you moved all ecalj binaries to `~/bin`, add the following to `.bashrc`.

```bash
PATH="~/bin/:$PATH"
```




## Uninstall {#_4-uninstall}

See [Uninstall](./uninstall). There are three commands: `python3 uninstall.py`,
open a new shell, and `rm -rf <ecalj-checkout>`.



## Installing VESTA and getsyml {#install-vesta-and-getsyml}

VESTA is useful for viewing structures.
(Example: I installed `VESTA-gtk3.tar.bz2` (ver. 3.5.8, built on August 11, 2022, 23.8MB) on Ubuntu 24.)

There is a `viewvesta` command in `ecalj/StructureTool/`. Try the following.

```bash
viewvesta ctrlg.si.toml      # canonical (post-2026-05)
viewvesta ctrls.si           # legacy structure-only seed file
```

This lets you check the structure in the viewer. `/StructureTool` contains the conversion tools
`vasp2ctrl` and `ctrl2vasp`. With them you can convert structures to and from POSCAR.
`ctrl2vasp` accepts the same three forms of input as `viewvesta`.

`getsyml`, which makes the symmetry lines for band plots, is installed into bin by `InstallAll.py` (as a soft link).
The generated `syml.*` files are used for band plots in ecalj.
If `spglib` and `seekpath` are available, usually no extra steps are needed.
For details see the upstream memo [ecalj/GetSyml/README.md](https://github.com/tkotani/ecalj/blob/main/GetSyml/README.md).



## Additional memo {#additional-memo}

* When `InstallAll.py` finishes, all binaries and shell scripts have been copied to `--bindir <path>` (default `~/bin/`).

* Cleanup:

  If something goes wrong, run `InstallAll.py` again with `--clean`.

  You may need to do the following.

  ```bash
  rm -f CMakefiles CMakeCache.txt
  ```

  Do this in `ecalj/SRC/exec/` (and also delete all `*.mod` files under `SRC/`). Usually this is not necessary.

* Compiling only the Fortran part:

  To compile only the Fortran sources, go to `ecalj/SRC/exec/` and run the following.

  ```bash
  FC=gfortran cmake . -D CMAKE_BUILD_TYPE=Debug  # or FC=ifort / FC=nvfortran
  make
  ```

  You may want to look at `CMakeLists.txt`. For a clean build, delete `CMakeCache.txt` and `CMakeFiles/` first.

* Compiler bugs: a compiler can fail with `-O2` optimization. As written in `CMakeLists.txt`, using a lower optimization option (`-O1` or `-O0`) depending on the source file usually avoids such bugs. We usually do not use `-O3`. You may set up such conditional compilation settings; see `CMakeLists.txt`.

* The source code and the make system are under `SRC/`, and the test examples have been
  moved to `Samples/`.

  ```text
  ecalj/
  ├── SRC/
  │   ├── exec        : CMakeLists.txt, build scripts, testecalj driver, comp.py
  │   ├── main        : All main *.f90
  │   ├── subroutines : All subroutines *.f90
  │   └── external    : bundled third-party (toml-f, ...)
  └── Samples/
      ├── TestInstall : Root of install tests (testecalj --all)
      └── EPS / PROCAR / MLOsamples / ... : topical samples (the old Legacy/ was rebuilt and removed on 2026-09-30)
  ```

  All Fortran code is in `SRC/main/` and `SRC/subroutines/`.
  `SRC/exec/CMakeLists.txt` generates the Makefile with cmake.
  (A remnant of the past: `SRC/TestInstall_is_moved_to_under_ecaljSamples`
  points to the old location of TestInstall.)

* The test system: the install tests are in `ecalj/Samples/TestInstall/`.

  Each subdirectory has its own `test.py` (the procedure of that test).
  The common driver is `ecalj/SRC/exec/testecalj` (installed into `BINDIR` as
  `testecalj`), and the common comparison helper is in
  `ecalj/SRC/exec/comp.py`. Adding your own test is not difficult.

  * First, compute something.
  * Put the input and minimal results in a new directory under `Samples/TestInstall/`
    (or somewhere under `Samples/`).
  * Write how to compare them in `test.py` in that directory.

  To run the install tests, do the following.

  ```bash
  cd ecalj/Samples/TestInstall
  testecalj --all -np 8                # all install tests
  testecalj --gwall -np 8               # GW-only subset
  testecalj si_gwsc nio_gwsc -np 8      # specific tests
  ```

* Note: on Ubuntu 22, I saw gfortran + openmpi fail occasionally. If you have problems, use mpich. (The current status is unknown.)
