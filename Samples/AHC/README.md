# Test for BCC Fe for Anomalous Hall conductivity: 

We check the case  with SO=2 (SOC with LzSz. This means without off-diagonal). 
We and 
the case SO=1 

# 0. About your python environment
Recent OS requires your own python environments independent from the python used by the system. 
Thus you need your own environments even when you do not need multiple environments. Many environment tools. venv is one of the standard included in python.

** venv case **:  
```
cd ~/
python3 -m venv .venv
source .venv/bin/activate    # This activete .venv. We use python and modules under .venv
```
This is a case you make .venv under your root directory. 
Then we have to install required packages as
```
pip install numpy mpi4py matplotlib
```

--- 

** Conda case **:  
If you use conda instead of venv
```bash
conda create -n mpipy python=3.10
```
You can confirm the new environment as follows.
```bash
conda info -e
```
Then, let us install mpi4py in your new environment.
```bash
conda activate mpipy
pip install mpi4py
```


# How to test?
We have two srcipts job1 and job2.
```
./job1
./job2
```
job1 is for preparation. job2 runs for different fermi energies.
`lmf` may take a little long because of many k points setting in ctrl.

Band plot poped up during job1.


# 1. AHC SO=2

Here, we calculate AHC with so=2 (lzsz mode).
We will use 'phispinsym' mode, which assumes the same atomic orbitals for up and down spin.
```bash
lmfa Fe > llmfa
mpirun -np 4 lmf Fe --phispinsym > llmf
job_band Fe -np 4 --phispinsym  > job_band.out
job_AHC Fe -np 4 --phispinsym > job_AHC.out
mpirun -np 4 python hx0ahc.py -4. 4. 51 --phispinsym > lhx0ahc
```
- ahc_tet.isp**.dat : tetrahedron
- ahc_sp.isp**.dat : direct sum

The spin index will be 'uu', 'ud', 'du', and 'dd'.
For example, 'ud' means the contributions from excitations from up to down state. 

# 2. AHC SO=1 without off-diagonal

Next, let's try AHC calculation with so=1 (l.s mode).
We omit the off-diagonal component of soc hamiltonian intentionally with '--testso' option.
Therefore, we should obtain the same results with so=2 case.
```bash
lmfa Fe > llmfa
mpirun -np 4 lmf Fe --phispinsym --testso -vso=1 > llmf
job_band Fe -np 4 --phispinsym --testso -vso=1 > job_band.out
job_AHC Fe -np 4 --phispinsym --testso -vso=1 > job_AHC.out
mpirun -np 4 python hx0ahc.py -4. 4. 51 --phispinsym --testso -vso=1 > lhx0ahc
```

3. plot the figures

```bash
gnuplot -persist plot_AHC_tet.glt
gnuplot -persist plot_AHC_direct.glt
```

### Tetrahedron method
![](figs/AHC_tet.png)

### Direct summation
![](figs/AHC_direct.png)