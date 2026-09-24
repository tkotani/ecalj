from setcomm import *
'''An example to run MPI-fortran codes successively.
   To get ecaljF.so, run >mpif90 -j -shared -fPIC -o ecaljfortran.so *.f90'''

ecaljF = setcomm("ecaljF.so")

callF( ecaljF.hello )
#callF( ecaljF.hello )
#callF( ecaljF.hello3, [True, 34, 8.2, 9.3] )
