# MTO-based localized orbital (MLO) generation

1. Run lmf and band calculations.
2. Set model bands in GWinput
3. Set
       mlo_method 1 (for GaAs) or
       mlo_method 2 (for Ni3d and O2p. To extract such narrow bands).
4. Then Run job_mlo. and See gnuplot -p bandplot_MLO.isp1.glt

>job_mlo gaas -np 8
>gnuplot -p bandplot_MLO.isp1.glt 



