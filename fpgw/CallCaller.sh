#!/bin/bash
echo ##############
echo  'This generate a call-caller data set for lmv7.F fp/*.F and subs/\*.F subs/*.F90 slatsm/*.F .'
echo  'No arguments required. This drives ../TOOLS/f_calltree.py '
echo  'HELP --> ../TOOLS/FparserTools/f_calltree.py --help, and read CallCaller.sh'
echo
echo  '--- Now generating a file 'callcaller.dat' ... Wait!!! It takes 1 minute or so!'
echo '        If you like to apply this to other programs, modify this script'
echo  ' NOTE: T.Kotani is not sure whether this is relaiable enough or not... let me know something wrong...'
echo ##############

../TOOLS/FparserTools/f_calltree.py main/hsfp0.m.F gwsrc/*.F >callcaller.dat 2>callcaller.err
#../TOOLS/FparserTools/f_calltree.py fp/*.F subs/*.F subs/*.F90 >callcaller.dat 2>callcaller.err

#../TOOLS/FparserTools/f_calltree.py lmv7.F lmfav7.F lmv7util.F fp/*.F subs/*.F subs/*.F9#0 slatsm/*.F >callcaller.dat 2>callcaller.err

#../TOOLS/FparserTools/f_calltree.py lmv7.F lmfav7.F lmv7util.F fp/*.F subs/*.F subs/*.F90 >callcaller.dat 2>callcaller.err
#../TOOLS/f_calltree.py lmv7.F fp/*.F >callcaller.dat 2>callcaller.err
egrep -e '^(ERROR|Error)' callcaller.err
echo
echo '------------------------------------------------------------------------------'
echo '--- If no ERROR is shown above (if ERROR is not in callcaller.err), it is succeeded. ---'
echo '       Note that Unsed files might be used by other mainprogram.'
echo '--- If ERROR is shown above, look into callcaller.err. Something wrong.'
echo 
echo ' If you want to make a callcaller-tree picture, try'
echo ' >GenCCtree.sh callcaller.dotdata'
echo ' --> Then you get ccmap.ps.; it is better to use smaller callcaller.dotdata(need to modify this script to make it).' 
echo ' Note that you need graphviz for GenCCtree.sh. as apt-get install graphviz'
