# Test for dielectric function for homogeneious electron gas (Lindhard function) 

The main rouitne is main/hhomogas.f90-->subroutines/main_hhomogas.f90

Look into job file before performing it.
Run
```bash
job
```
When job have finished, you will get x0homo.dat 
The original results is saved as x0homo.result.dat.

---
Let me review the job. 

This is for the empty sphere. (now Z=0.0001 to detour stop(I will fix)).
The lattice is taken from bulk Li case.
job does band plot at first.
Then qg4gw gives the q points and tetrahedron division information. 

Then we finally run the key routine hhomogas. 
Please look into main_hhomogas.f90
Some parameter settings are fixed in the code.
So you may need to change them for your purpose (e.g. change fermi energy).

Repeat only hhomogas to learn how to use tetrahedron method.
When you change number of k points use job2 after setting n1n2n3 in GWinput.

---

Our final results is the plot by
>gnuplit -p egaschi0.glt
for the final results x0homo.dat,
which contains the real and imag part of the non-interacting Lindhard polarization function x0.

See Fetter-Walecka https://drive.google.com/file/d/0B4EhZYLsWtz-VGFoRlRNcU1qVzA/view?usp=drive_link&resourcekey=0-POrF9ZW73o-IalpgR9cXWA
Fig.12.9 is reproduced. 

Set any q points in main_hhomogas.f90:L279

Note volume factor vol.
This is contained in the definition of rcxq (=imag part of zxq).