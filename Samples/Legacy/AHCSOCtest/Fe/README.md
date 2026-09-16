# Test for BCC Fe

## Anomalous Hall conductivity: check difference of SO=2 and SO=1 without off-diagonal

0. About mpi4py

Before calculation, please install mpi4py in your environment.
I recommend you to create a new environment for mpi4py (for example mpipy) by the following command.
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

1. first step: lmf calculation and construct UU matrix

You can run the test as follows. You should modify "exe" in job1 before the calculation. 
```bash
./job1
```

Or you can also manually perform the calculation.

For AHC with so=2 (lzsz mode):
```bash
mkdir k4_so2; cd k4_so2
cp ../*.Fe .; cp ../GWinput .
lmfa Fe > llmfa
mpirun -np 4 lmf Fe --phispinsym > llmf
job_band Fe -np 4 --phispinsym  > job_band.out
job_AHC Fe -np 4 --phispinsym > job_AHC.out
cd ..
```
Here, we will use 'phispinsym' mode, which assumes the same atomic orbitals for up and down spin.

For AHC with so=1 (l.s mode):
```bash
cp -r k4_so2 k4_so1; cd k4_so1
mpirun -np 4 lmf Fe --phispinsym --testso -vso=1 > llmf
job_band Fe -np 4 --phispinsym --testso -vso=1 > job_band.out
job_AHC Fe -np 4 --phispinsym --testso -vso=1 > job_AHC.out
cd ..
```
We omit the off-diagonal component of soc hamiltonian intentionally with '--testso' option.
Therefore, we should obtain the same results with so=2 case.

2. second step: AHC as a function of Fermi energy

You can run the test as follows.
```bash
./job2
```

Or you can also manually perform the calculation.

Before the calculation, you need to activate the environment, which you can use mpi4py
```bash
conda activate mpipy
```

Then, perform the following command:
```bash
cd k4_so2
mpirun -np 4 python hx0ahc.py -4. 4. 51 --phispinsym > lhx0ahc
cd ..
cd k4_so1
mpirun -np 4 python hx0ahc.py -4. 4. 51 --phispinsym --testso -vso=1 > lhx0ahc
cd ..
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
    

4. for finer k-point mesh

You should use cluster computers to perform the calculation.
```bash
qsub job.sh
```
job.sh conducts job_k4-k32_1st and job_k4-k32_2nd

```bash
gnuplot -persist plot_AHC_direct_k4-k32_so2.glt
gnuplot -persist plot_AHC_direct_k4-k32_so1.glt
```