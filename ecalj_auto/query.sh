#!/bin/bash
# I use pymatgen 2024.10.3 with Python 3.13.3
#       pymatgen 2023.5.10 with Python 3.9.13

time=$(date +"%Y%m%d-%H%M%S")
dir=$time
### for test (Si Ga Al)
auto_mpquery.py --nsites 1 2 --include Si Ga Al --dir=testSGAtest

exit

#auto_mpquery.py --nsites 1 2 --include Si Ga Al --dir=$dir/SimpleTest
#auto_mpquery.py --nsites 2 3 --fname joblist_nsites3 --dir=testSGA_nsite23

### gw1000 ###
auto_mpquery.py --dir=gw1500 --nsites 1 8 --fname joblist

exit
#auto_mpquery.py --dir=$dir/default
auto_mpquery.py --dir=gw1000
auto_mpquery.py --dir=gw1500 --nsites 1 8 

### a sample including Tl only but with the gw1000 condition. No LN AC, only expt
#auto_mpquery.py --include Tl --dir testTl

### for compounds with 4f-elements (LN, AC)
auto_mpquery.py --include LN AC --exclude NG --dir=w4f

### for magnetic compounds
auto_mpquery.py --mag true --dir=mag

### for magnetic compounds with not excluding NG, LN, AC
auto_mpquery.py --mag true --exclude none --dir=mag_all

### for a dataset to compare with experimental bandgaps
auto_mpquery.py --nsites 1 100 --mpid expt --mag false --theoretical both --metal both --dir=expt
