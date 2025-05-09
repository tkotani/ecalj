## ecalj: Mmaximally localized Wannier with effective interaction in CRPA

NOTE: Maxloc will be replaced with MLO

* SRC/wanniergw/ contains the Wannier function generator on top of ecalj resouces.

* Srcipt genMLWF
>genMLWF 
is the script to generate the Wanneir functions.
In addtion, it gives effective interaction <ij|W|kl> in CRPA.
Required setting is written in the GWinput file (mkGWIN_lmf2 generate template in GWinput.tmp).
 

* How to run samples?
** We have several samples for genMLWF
Samples are: at ecaj/Samples/MLWF_samples
├── CuMLWFs n1n2n3 4 4 4, rough but small sample to test it at firts. Within ~five minutes.
├── Cu
├── Fe
├── La2CuO4
├── NiOMLWF
├── SrVO3MLWF

** Run a sample at Samples/MLWF_samples/CuMLWFs
Run self-consistent calculation as
>lmfa cu
>mpirun -np 4 lmf-MPIK cu

Then we make band plot by
>job_band cu -np 4
, where 4 means number of cores.
This gives correct Fermi energy, and stored into bnds.* file

Then we run main script of maxloc wannier with effective interaction W as
>genMLWF -np 12 cu

(1) Look into ecalj/Samples/CuMLWFs/bandplot.MLWF.cu.glt
    This is for interpolated band.
    A line "bnds.maxloc.up" u ($5):($6+de) lt 3 w l ti "Wannier" is added to usual output of
    bandplot.cu.isp* given by job_band.
    
(2) Plot psi.xsf file. 
    This contains MaxLoc Wannier function. (not necessary; you can
    comment out "wanplot" line if you don't need to plot Wannier.)

(3) grep Wan lwmat*
    This gives the matrix element of the Coulomb interaction and W-v.
    (off-diagonal elments are also calculated. Sorry, I have not yet
    document how to read this...)
    
(4) We have to specify inner window (<worb> section in GWinput)
    and so on related to wannier in GWinput.
    The template is given by mkGWIN_lmf2 from ctrl file.

** OUTPUT

*** Wannier funciton plot
  Use wan.xsf by Xcrysden

*** Band plot including Wannier band
  See bandplot.MLWF.isp1.glt

*** Plot U, J parameters.

1) Run the srcipt genMLWF; please look into this script.
   (originally at SRC/exec/genMLWF ).

   This is made from two stages. 
   >argin=2; run_arg $argin $NO_MPI $nfpgw /hmaxloc lmaxloc2
   (equilavent to echo 2| hmaxloc >lmaxloc2) is the end of Wannier function generation.

   At this point, you can make band plot to check whether your setting for
   Wannier work well or not; the model-Hilbert space by band plot.
   (we need syml.* file and run job_band to get original energy bands)
   Then plot wannier band on top of it. See
   Samples/MLWF_samples/CuMLWFs/bandplot.MLWF.isp1.glt as an example.
   
`   If the plot is strange, you need to choose windows for Wannier.
   (Repeat echo 2| hmaxloc >lmaxloc2 until you have satisfactry
    fitting with changing the setting wannier part in GWinput ).

2) We get three files (see genMLWF)
      three files containing v and W-v information.
          grep "Wannier" lwmatK1 > Coulomb_v
          grep "Wannier" lwmatK2 > Screening_W-v
          grep "Wannier" lwmatK3 > Screening_W-v_crpa

These are text files <ab|W|cd> element
a,b,c,d are index of Wannier functions (sorry, ask us if necessary).

  Then we have Static_W.dat (RPA) and Static_U.dat (cRPA)
  These contains static U, U', J, and J' (\omega = 0).

grep '    1    1    1    1    1'             Coulmb_v
grep '    1    1    1    1    1    0.000000' Screening_W-v.UP 
grep '    1    1    1    1    1    0.000000' Screening_W-v.crpa
shows
 Coulomb_v.UP:          Wannier     1    1    0.000000    0.000000    0.000000    1    1    1    1    1   23.499183   -0.000000
 Screening_W-v.UP:      Wannier     1    1    0.000000    0.000000    0.000000    1    1    1    1    1    0.000000    0.000000  -20.317956   -0.000000
 Screening_W-v_crpa.UP: Wannier     1    1    0.000000    0.000000    0.000000    1    1    1    1    1    0.000000    0.000000  -20.188076   -0.000000
This means 
<11|W|11>     =23.499183-20.317956
<11|U_CRPA|11>=23.499183-20.188076
Note that this is a test example.

* NOTE: 
For practical calculations, we need to stop right after a line
>argin=2; run_arg $argin $NO_MPI $nfpgw /hmaxloc lmaxloc2  
in the genMLWF. 
(Need to modify genMLWF by hand; note that "run_arg" is a subrouitne of bash. 
 This line is equivalent to  >echo 2| hmaxloc >& lmaxloc2)
Then we need to check the band plot, with superposing Wannier band
"bnds.maxloc.up" u ($5):($6).
If the Wannier band is not reasonable, change settings of inner and outer windows,
and run >echo 2| hmaxloc |tee lmaxloc2
again and again until you have a reasonable Wannier band plot.
Then go ahead to the next steps (you may use genMLWF2).
-----------
genMLWF  
      I added the "grep" commands in the last part of genMLWF to make

      hwmatK_MPI.F
         print all matrix elements of W and frequency in eV. called in genMLWF

* CAUTION:
you must run job_band in advance to genMLWF!
If not you may need to follow instruction of "Efermi shift" as
follows;
------
NOTE: Efermi shift:
 genMLWF requies bnds.${target} to read the Fermi energy.
 To generate it, we need to run job_band in advance.
 Or run, 
 >echo 2 | hmaxloc  > lmaxloc2
 (need syml*); this can be runned after genMLWF.
----
(Or need to shift Ef by hand as follows in gnuplot script.)
 ----------------------------------------------------
 de = ((ef shown in "lmaxloc2") - (ef in llmf_ef(bnds.${target}))*13.605
 plot \
 "bnd1.dat" u 2:3 lt 1 pt 1 not w l,\
 "bnd2.dat" u 2:3 lt 1 pt 1 not w l,\
 "bnd3.dat" u 2:3 lt 1 pt 1 not w l,\
 "bnd4.dat" u 2:3 lt 1 pt 1 not w l,\
 "bnd5.dat" u 2:3 lt 1 pt 1 not w l,\
 "bnds.maxloc.up" u ($5):($6+de) lt 3 w l ti "Wannier"

* Known problems:
** Range of plot looks not good;
   Especially, vis_wan_ubound, vis_wan_lbound should be not integer.
   Probably, need to improve/(bug fix) wanplot.F.
** xsf files (phi.xsf, wan.xsf) are special for Xcrysden,
   thus it might be inconvenient. 
   You need to use GUI of Xcrysden. But I have not tested no other format.

* development history:
Started from Apr2015, T.kotani take developments of MaxLoc Wannier,
which is implemented in the ecalj by T.Miyake R.Sakakma and H.Kino.
2009; maxloc090910 (T.Miyake)
2009; Its documentation and Visualizer at Kino's https://github.com/nim-hrkn/visualize
2014Aug; T.Kotani modified it so as to fit to latest ecalj.
2015jan: Sengwoo Jang's modification.
2015apr: Sengwoo's RPAWannierTable.py.

2020aug: reorganize and checke README_wannier.org
