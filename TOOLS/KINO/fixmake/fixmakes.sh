#!/bin/sh
path=/home/kino2/GW/ecalj_2011_0520-0602/TOOLS/KINO/fixmake
src=MAKEINC/Make.inc.ifort
 gawk -f $path/fromfiles.awk subs/*.F fp/*.F 
 cat $src __makefile__ | gawk -f $path/frommake.awk > x
 mv x $src

src=MAKEINC/Make.inc.ifort_mpi
 gawk -f $path/fromfiles.awk subs/*.F fp/*.F 
 cat $src __makefile__ | gawk -f $path/frommake.awk > x
 mv x $src


src=MAKEINC/Make.inc.ifort_mpik
 gawk -f $path/fromfiles.awk subs/*.F fp/*.F 
 cat $src __makefile__ | gawk -f $path/frommake.awk > x
 mv x $src

src=MAKEINC/Make.inc.gfortran
 gawk -f $path/fromfiles.awk subs/*.F fp/*.F 
 cat $src __makefile__ | gawk -f $path/frommake.awk > x
 mv x $src

src=MAKEINC/Make.inc.gfortran_mpi
 gawk -f $path/fromfiles.awk subs/*.F fp/*.F 
 cat $src __makefile__ | gawk -f $path/frommake.awk > x
 mv x $src

src=MAKEINC/Make.inc.gfortran_mpik
 gawk -f $path/fromfiles.awk subs/*.F fp/*.F 
 cat $src __makefile__ | gawk -f $path/frommake.awk > x
 mv x $src

