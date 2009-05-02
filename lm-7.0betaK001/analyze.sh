#!/bin/bash
echo ##############
echo  This script do not require arguments. 
echo    
echo    This generate a call tree file 'lmfp_tree' for fp/lmfp.f, which is called from the main routine lmv7.f.
echo    If you like, change keyword 'lmfp' in TOOLS/ANALYZE/analyze.py
echo    This routine generate a tree up to 3rd lower level.
echo    
echo  --- Now generating a file 'lmfp_tree' ... Wait 10~30sec or so! 
echo ##############
../TOOLS/ANALYZE/analyze.py lmfp fp/*.F subs/*.F slatsm/*.F gwd/*.F >callcaller.dat ; grep tree callcaller.dat > lmfp_tree

