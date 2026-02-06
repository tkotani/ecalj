# MTO-based localized orbital (MLO) generation and W for MLO.

1. Run lmf and band calculations.
2. Set model bands in GWinput. See
3. Set
      mlo_method 0
      mlo_emax 0 for seimconductor, or 7 for Al2O3_Cr. (above energy eV for localized orbitals relative to Efermi).
4. Then Run job_mlo. and See gnuplot -p bandplot_MLO.isp1.glt

>job_mlo gaas -np 8
>gnuplot -p bandplot_MLO.isp1.glt 

-----

5. To get W in RPA, run
>job_mloW fe -np 8

And then see the sum of Coulomb_v.* and Screening_W-v.*
(we need minor fix for cRPA mode if necessary)

