#!/bin/bash
echo ##############
echo  This script do not require arguments. 
echo    
echo    This generate a call-caller data for *.F files in lmf.
echo    
echo  --- Now generating a file 'callcaller.dat' ... Wait 1 minutes or so! 
echo ##############
#../TOOLS/ANALYZE/analyze.py freeat fp/*.F subs/*.F slatsm/*.F gwd/*.F >callcaller.dat ; grep tree callcaller.dat > lmfa_tree
../TOOLS/f_calltree.py fp/*.F subs/*.F slatsm/*.F gwd/*.F >callcaller.dat 2>callcaller.err

