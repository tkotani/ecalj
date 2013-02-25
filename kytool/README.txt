kytool is a converter between ctrl(or ctrls) and POSCAR.
Just treat crystal structure part of ctrl file.
If you see something wrong or inconvenient, let t.kotnai know it.

========
INSTALL:
 You may need to add path to here, kytool/ to invoke following commands such as vasp2ctrl.py.

 As a viewer, we  can use VESTA (http://jp-minerals.org/vesta/jp/download.html).
 This is required to invoke viewer by a command viewvesta.py in the followings
 At the begining of viewvesta.py. Set VESTA path.

========
Usage:
We now have three commands below. These command uses functions in convert/, so relative path is important.

*vasp2ctrl.py
　Convert POSCAR to ctrl file. For help, type this command without
  arguments.
  E.g: 
    mkdir TEST
    cd TEST
    cp ../sample/ctrl.nio .
    ../vasp2ctrl.py ctrl.nio 
　
*ctrl2vasp.py
  ctrl file to POSCAR file. (current version is for Cartesian, but not so difficult if you like Direct).
  E.g: 
    mkdir TEST
    cd TEST
    cp ../sample/10-Opal.cif.vasp POSCAR_opal
    ../vasp2ctrl.py POSCAR_opal

*viewvesta.py
  this is to see POSCAR by VESTA. VESTA is poor for external control.
  probably, it is better to use another kinds of viewer in future.
  Currently, this is a very simple program just to invole VESTA for
  given file.

You can test it in your new directory with some files in sample/ directory.
Crystal structure files *.cif are taken from http://staff.aist.go.jp/nomura-k/japanese/itscgallary.htm




