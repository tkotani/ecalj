      subroutine defwsr(wsr,z)

C- Returns default value of radii for given nuclear charge
C ----------------------------------------------------------------------
Ci Inputs:
Ci   z     :nuclear charge
Co Outputs:
Co   wsr   :Wigner-Seitz sphere radius (in atomic units)
Cr Remarks:
Cr These are equilibrium average Wigner-Seitz radii for closed-packed
Cr structure of pure element.
Cu Updates
Cu   19 Apr 02 Adapted from Stuttgart LMTO56
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      double precision wsr,z
C Local variables:
      integer iprint,iz,lgunit
      double precision cpwsr(0:100)
C External calls:
      external isanrg
C Intrinsic functions:
      intrinsic idnint
C Data statements:
      data cpwsr/2.00d0,
     .1.39d0,2.55d0,3.04d0,2.27d0,1.96d0,1.66d0,1.90d0,1.90d0,
     .2.17d0,2.89d0,3.76d0,3.25d0,2.95d0,2.63d0,2.56d0,2.70d0,
     .2.85d0,3.71d0,4.66d0,3.88d0,3.31d0,2.99d0,2.76d0,2.64d0,
     .2.57d0,2.52d0,2.52d0,2.55d0,2.62d0,2.78d0,2.75d0,2.79d0,
     .2.83d0,2.94d0,3.13d0,4.32d0,4.95d0,4.22d0,3.61d0,3.28d0,
     .3.03d0,2.91d0,2.82d0,2.77d0,2.78d0,
     .2.84d0,2.95d0,3.14d0,3.30d0,3.45d0,
     .3.30d0,3.31d0,3.50d0,4.31d0,5.30d0,
     .4.20d0,3.91d0,3.80d0,3.75d0,3.70d0,
     .3.65d0,3.60d0,3.55d0,3.52d0,3.61d0,
     .3.67d0,3.70d0,3.73d0,3.75d0,3.56d0,
     .3.44d0,3.23d0,3.04d0,2.93d0,2.86d0,
     .2.82d0,2.83d0,2.88d0,2.98d0,3.27d0,
     .3.57d0,3.62d0,3.37d0,3.46d0,3.63d0,
     .4.44d0,5.81d0,4.30d0,3.84d0,3.52d0,
     .3.32d0,3.13d0,3.02d0,2.96d0,2.93d0,
     .2.93d0,2.95d0,2.99d0,3.05d0,3.17d0/

      iz = idnint(z)
      call isanrg(iz,0,100,'DEFWSR:','z',.true.)
      wsr = cpwsr(iz)
      if (iprint() .ge. 120) write(lgunit(1),300) iz,wsr

  300 format(/' DEFWSR: default radius for z=',i2,': ',f8.2)
      end subroutine defwsr


