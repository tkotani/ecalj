---
project: calj: MAE calculation based on Force theorem. See Liqin PhysRevB.99.054418,
author: takao kotani
email: takaokotani@gmail.com
---
LICENCE: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.html)

## Calculate MagneticAnisotropyEnergy by the dfference of band energy.
We calculation band energy difference with SOC. We calculate one-shot band energy with SOC,
on top of the result without SOC.

See appendix A of [K.Liqin](https://link.aps.org/doi/10.1103/PhysRevB.99.054418). 
We compare contributions by the SOC Hamiltonian $H_{\rm SOC(100)}$ or by $H_{\rm SOC(110)}$.

hammhso in m_bandcal.f90 contains the SOC Hamiltonian,
which is generated in aughsoc.f90 defined in m_augmbl.f90.
I implemente the SOC Hamiltonian following the appendix A.

## Example for FePt: 
We explain how to calculate MAE from the difference of band energy.
#### Step0. Settings in ctrl file.
In ctrl, we set
```
nspin=2 
MMOM=0,0,2.2   !for Fe site initial condition.
LMXA=6         !for better results (because we enlarge l max to expande eigenfunction within MT)
nk1,nk2,nk3 8 8 6   !this is for test. We may need more, such as nk 16 16 12.
xcfun=103      !GGA-PBE (or =1 for LDA. I am not so sure how it effects to.)
so=0           !no SOC at first
PWMODE=11      !|q+G| cutoff. q-dependent APW basis. This is (probably better) to keep crystal symmetry 
               !PWMODW=1 is q-independent APW basis.
%const soa1=0 soa2=0 soa3=1 
    SOCAXIS={soa1} {soa2} {soa3} 	!Set spin axis. This is in HAM category
```
See ctrl.fept. In addition, we may(or maynot) need to supply accurate atomic positions.

#### Step1. Self-consistent calculation.
```
mpirun -np 4 ~/ecalj/SRC/exec/lmf fept --phispinsym -vso=0 |tee llmf`
```

Note that --phispinsym is to use the same basis function independent of spins.
Then the end of save.fept file shows the totale energy and band energy as
```
Start LMF fept --phispinsym -vso=0
...
c mmom= 3.3033 ehf(eV)=-535698.889261 ehk(eV)=-535698.889267 sev(eV)=-367.926757
```
ehf(Harris-Folkner energy) should be in agreement with ehk (Hohenberg-Kohn energy) when
we have perfect convergence except numerical error
(thus the agreement is a check to confirm calculations are performed correctly).
Note that sev shows the band energy.

* NOTE:
 If we do not use spin-symmetric basis without --phispinsym, we got
```
c mmom= 3.3218 ehf(eV)=-535698.910109 ehk(eV)=-535698.910140 sev(eV)=-367.556242
```

The difference from --phispinsym may indicate slight poorness of basis set. (We may need 4d for PZ?---not checked yet).

#### Step2.(only for QSGW). Extend sigm to sigm_fbz. 
```
mpirun -np 1 lmf fept --phispinsym --wsig_fbz 
```

This generates sigm_fbz.fept (Sigma-Vxc), which is for all k points.
Original sigm_fbz contains Sigma-Vxc only at irreducibel k points.

#### Step 3. Get band energy for SO=0, SO=1(001) and SO=1(110) 
Run three calculations. These are one-shot band calculations with no space-group symmetry.
```
mpirun -np 4 ~/ecalj/SRC/exec/lmf fept --phispinsym --nosym --quit=band -vso=0 >& llmfso0 
mpirun -np 4 ~/ecalj/SRC/exec/lmf fept --phispinsym --nosym --quit=band -vso=1 -vsoa1=0 -vsoa2=0 -vsoa3=1 >& llmso001
mpirun -np 4 ~/ecalj/SRC/exec/lmf fept --phispinsym --nosym --quit=band -vso=1 -vsoa1=1 -vsoa2=1 -vsoa3=0 >& llmso110
```

Then we get save.copt showing that 
```
...
Start LMF fept --phispinsym --nosym --quit=band -vso=0
i mmom= 3.3031 ehf(eV)=-535698.889270 ehk(eV)=-535698.889280 sev(eV)=-367.932517
Start LMF fept --phispinsym --nosym --quit=band -vso=1 -vsoa1=0 -vsoa2=0 -vsoa3=1
i mmom= 3.2702 ehf(eV)=-535699.203668 ehk(eV)=-535699.193569 sev(eV)=-368.246915
Start LMF fept --phispinsym --nosym --quit=band -vso=1 -vsoa1=1 -vsoa2=1 -vsoa3=0
i mmom= 3.2645 ehf(eV)=-535699.201112 ehk(eV)=-535699.190232 sev(eV)=-368.244360
```

Thus our result shows that 001 is stable with MAE= -368.246915(001) - (-368.244360(100)) = 2555 micro eV.
You may need to enlarge number of k points for better results.

- The first calculation is to reproduce the original at Step1 without space-group symmetry.
We see slight difference -367.932517 from -367.926757 for sev 
because of numerical error(good enough?). We calculate one-shot band energy with SOC,
on top of the result without SOC.

- Option --quit=band do not alter rst file. 
 sev is the band energy for the potential in rst file.
 This sev is the reference band energy without SO. 

- ` HAM_SOCAXIS       opt    r8v      3,  3,   *,  3       0 0 1`
is shown at the beginig of console output 
We allow only 110 and 001 for SOAXIS now because of the setting at the begging of m_augmbl.f90.
It should be not so difficult to modify this code.

- We use --phispinsym (the same orbital for up and dn) for lmf.
We need good convergence (at least for lmf with --phispinsym).
Good congergence of lmf should assure numerically accurate space-group 
symmetry of rst file, as well as sigm file.


## Example for MnGa: 
We do similar procedure (we use nk=12 12 12 in the next example. 
I saw better congergence rather than nk=888 for xcfunc=103).
```
mpirun -np 4 lmf mnga -vso=0 --phispinsym >llmfso0
mpirun -np 4 ~/ecalj/SRC/exec/lmf mnga --phispinsym --nosym --quit=band -vso=0 >& llmfso0 
mpirun -np 4 ~/ecalj/SRC/exec/lmf mnga --phispinsym --nosym --quit=band -vso=1 -vsoa1=0 -vsoa2=0 -vsoa3=1 >& llmso001
mpirun -np 4 ~/ecalj/SRC/exec/lmf mnga --phispinsym --nosym --quit=band -vso=1 -vsoa1=1 -vsoa2=1 -vsoa3=0 >& llmso110
```

Then we have save.mnga
```
Start LMF mnga --phispinsym --nosym --quit=band -vso=0
i mmom= 2.4771 ehf(eV)=-84429.149985 ehk(eV)=-84429.150004 sev(eV)=-440.821571
Start LMF mnga --phispinsym --nosym --quit=band -vso=1 -vsoa1=0 -vsoa2=0 -vsoa3=1
i mmom= 2.4782 ehf(eV)=-84429.169923 ehk(eV)=-84429.169937 sev(eV)=-440.841509
Start LMF mnga --phispinsym --nosym --quit=band -vso=1 -vsoa1=1 -vsoa2=1 -vsoa3=0
i mmom= 2.4783 ehf(eV)=-84429.169524 ehk(eV)=-84429.169539 sev(eV)=-440.841110
These shows -440.841509-(-440.841110)=-399 microeV (xcfun=103 GGAPBE)
```

> ---old result (2021?.just memo)----\
> c so=0 mmom=2.4082439 ehf=-6195.6056894 ehk=-6195.6056905\
> i so=0 mmom=2.4082689                      ehf=-6195.605 6905 ehk=-6195.6056917 sumev=  -32.994 029\ 
> i so=1 soa1=0 soa2=0 soa3=1 mmom=2.4094348 ehf=-6195.607 1509 ehk=-6195.607152  sumev=  -32.995 490\ 
> i so=1 soa1=1 soa2=1 soa3=0 mmom=2.4091936 ehf=-6195.607 1208 ehk=-6195.607122  sumev=  -32.995 459\
> ===> gave -421 micro eV for xcfun=1.


----
#### NOTE for implementation and improvements.
I recommend you to see ovarall calling flow lmfp->bndfp->m_band_cal and m_subze_bzintegration.
Note that sev is the band energy.
At m_bandcal_init, we set up Hamiltonian at m_bandcal.f90@L129
          hamm(:,:,1:2)= hamm(:,:,1:2)+hammhso(:,:,1:2) !spin-diag SOC elements (1,1), (2,2) added
          hamm(:,:, 3) = hammhso(:,:,3)                 !spin-offdiagonal SOC elements (1,2) added
SO part is included in hammhso(:,:,isp) where isp=3 is right-upper block of spin matrix.
hammhso is generated at call aughsoc in m_bandcal.
Look into hammhso if necessary 
(T.Kotani implemented only the case 001 and 110, but it should be easy to make it general case).
We set the SOC axis by HAM_SOCAXIS.


	


