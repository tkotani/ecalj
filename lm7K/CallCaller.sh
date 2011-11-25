#!/bin/bash
echo ##############
echo  'This generate a call-caller data set for fp/\*.F and subs/\*.F .'
echo  'No arguments required. This drives ../TOOLS/f_calltree.py '
echo  'HELP --> ../TOOLS/FparserTools/f_calltree.py --help, and read CallCaller.sh'
echo
echo  '--- Now generating a file 'callcaller.dat' ... It takes 1 minute or so!'
echo
echo ##############
#../TOOLS/ANALYZE/analyze.py freeat fp/*.F subs/*.F slatsm/*.F gwd/*.F >callcaller.dat ; grep tree callcaller.dat > lmfa_tree
../TOOLS/FparserTools/f_calltree.py lmv7.F lmfav7.F lmv7util.F fp/*.F subs/*.F subs/*.F90 slatsm/*.F >callcaller.dat 2>callcaller.err
#../TOOLS/FparserTools/f_calltree.py lmv7.F lmfav7.F lmv7util.F fp/*.F subs/*.F subs/*.F90 >callcaller.dat 2>callcaller.err
#../TOOLS/f_calltree.py lmv7.F fp/*.F >callcaller.dat 2>callcaller.err
egrep -e '^(ERROR|Error)' callcaller.err
echo
echo '------------------------------------------------------------------------------'
echo '--- If no ERROR is shown above, it is succeeded. To check, tail callcaller.err.'
echo '--- If ERROR is shown above, look into callcaller.err. Something wrong.'
echo 
echo ' If you want to make a callcaller tree, try'
echo ' >GenCCtree.sh callcaller.dotdata'
echo ' --> Then you get cc.ps.' 
echo ' Note that you need graphviz for GenCCtree.sh.'
