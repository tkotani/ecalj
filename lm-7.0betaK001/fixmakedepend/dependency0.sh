#!/bin/bash
#
# execute at the directory where Makefile exists ( or you can change path if you change this script.) 
#
# usage: bash fixmakedepend/thisscript targetmakefile_or_include
# e.g.   bash fixmakedepend/thisscript MAKEINC/Make.inc.ifort 
#
#
topdir=fixmakedepend
target=$1
# write modulelist and uselist
grep -i '^ *module' *.F  */*.F | grep -v -i 'module *procedure'  > ${topdir}/modulelist
grep -i '^ *use [A-Za-z]' *.F */*.F   > ${topdir}/uselist
# write dependency of makefile
(cd ${topdir}; gawk -f makedepend.awk > depend0 )
# fix duplicated dependency 
gawk -f ${topdir}/fixmakedepend.awk  ${topdir}/depend0  ${target}
(cd ${topdir}; rm -f depend0   modulelist  uselist)

