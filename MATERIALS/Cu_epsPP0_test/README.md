# Get epsilon(omega) at q=0

In ecalj, we have to take q->0 limit to obtain epsilon(omega,q=0).
No local field correction.

1.
Let cleanup directory with ctrl.cu and GWinput only
Get converged results.

2.
Make sure, GWinput contains
-------------------
QforEPSunita on 
<QforEPS>
 0d0 0d0 0.001
 0d0 0d0 0.0014142
 0d0 0d0 0.002
 0d0 0d0 0.0028284
 0d0 0d0 0.004
</QforEPS>
-------------------
. 'QforEPSunita on' means QforEPS is in unit of 2pi/alat.
T.K think this setting is good enough for all materials.
But not checked well.
Curent version askes you to set like this (x,y direction as well).
To get good results, we have to use large enough number of k points (n1n2n3 in GWinput).
Roughly speaking, 16x16x16 or more for Si.

3.
Run
>epsPP0 -np 4 cu
This gives
a.epsilon interband contribution
  epsinter.dat  epsinter.glt
b.eplison intraband contribution
  epsintra.dat  epsintra.glt
c. sum of them
  epsall.dat  epsall.glt  
------------
The bare data files EPS* contains numerical errors.
So I implemented a new scheme to cancel out numerical error to get epsilon at q=0.

4.
Run
>gnuplot -p epsinter.glt
gnuplot files epsinter.glt and epsintra.glt are for check.
Final results are summarized in espall.dat

Imaginary part of interband contribution is only at omega=0.
This is the Fermi surface contribution.
Because of our numerical technique, we use finite q (but very small).
We can zoom up to see the Im part of Imaginary part.
It scales well, resulting almost the same (well-superposed) real parts. 

2023-05-08
