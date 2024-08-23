# For dielectric function for homogeneious electron gas (Lindhard function)

The main rouitne is main/hhomogas.f90-->subroutines/main_hhomogas.f90

Look into job file before performing it. When job finishes, you will get x0homo.dat 
The original results is saved as x0homo.result.dat.

(in principle, you can skip --jobgw=1 , but some small files may be needed).

---
Let me review the job. 

This is for the empty sphere. (now Z=0.0001 to detour stop(I will fix)).
The lattice is taken from bulk Li case. 
We do a band plot at first.
Then qg4gw gives the q points and tetrahedron division information. 

Then we finally run the key routine hhomogas. 
Please look into main_hhomogas.f90
Some parameter settings are fixed in the code. So you may need to
change them for your purpose (e.g. change fermi energy).


Repeat hhomogas to learn how to use tetrahedron method.

---

Our final results is the plot by egaschi0.glt for the final results x0homo.dat,
which contains the real and imag part of the non-interacting Lindhard polarization function x0.

You should get a smooth good result for enlarged n1n2n3 in GWinput
(it gets closer to the analytic result of Lindhard. See Fetter & Walecka, for example).
As egaschi0.glt show only a part of x0homo.dat,
you may need to edit egaschi0.glt to show more data in x0homo.dat.

For enlarged n1n2n3, you may run only a few cycles of do 1001 loop to
reduce computational time (do 1001 iq = iqxini,iqxend ).
The loop index iq corresponds to the irreducible BZ q points for which we calculate x0.

We use single band case for any q (we take only the lowest
band in this example). Thus results are correct only for small q case,
in the current setting of main_hhomogas.f90.

