kytool is a converter between ctrl(or ctrls) and POSCAR.
Just treat crystal structure part of ctrl file.

========
INSTALL:
 You may need to add path to here, kytool/ to invoke following commands such as vasp2ctrl.py.

 As a viewer, we  can use VESTA (http://jp-minerals.org/vesta/jp/download.html).
 This is required to invoke viewer by a command viewvesta.py in the followings
 At the begining of viewvesta.py. Set VESTA path.

========
Samples:
 There are samples in sample/ directories.
 Crystal structure files *.cif are taken from http://staff.aist.go.jp/nomura-k/japanese/itscgallary.htm

========
Usage:
We now have three commands below. These command uses functions in
convert/, so relative path is important.

*vasp2ctrl.py
　Convert POSCAR to ctrl file. 
  For help, type this command without arguments.
   E.g: ../vasp2ctrl.py ctrl.nio --Cartesian
　
*ctrl2vasp.py
  ctrl file to POSCAR file. 
  Note that you have to supply --Direct or Cartesian.
   E.g: ../ctrl2vasp.py ctrl.nio --Cartesian

*viewvesta.py
  this is to see POSCAR by VESTA. VESTA is poor for external control.
  probably, it is better to use another kinds of viewer in future.
  Currently, this is a very simple program just to invole VESTA for
  given file.


You can test it in your new directory with some files in sample/
directory.

====
this kytool/ is by Kazuyoshi Yamamoto, and modified by t.kotani.
Let me know if you find something wrong.




