      subroutine getmom(nsp,nl,qnu,nclass,nrclas,idxdn,zval)
C- Get total number of valence electrons from moments
C ----------------------------------------------------------------------
Ci Inputs
Ci   nsp,nl,nclass,nrclas
Ci   qnu(1,l,isp,ic) - zeroth moments for each L, spin, and class
Ci   idxdn(l,ic) - downfolding switches (if > 1, exclude orbitals)
Co Outputs
Co   zval: total number of valence electrons
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nsp,nl,nclass,nrclas(1),idxdn(nl,nclass)
      double precision zval,qnu(3,nl,nsp,1)
C Local variables
      integer l,ic,isp
      double precision d1mach

      do  3  ic = 1, nclass
        do  2  isp = 1, nsp
          do  1  l = 1, nl
            if (idxdn(l,ic) .gt. 1) goto 1
            zval = zval + nrclas(ic)*qnu(1,l,isp,ic)
    1     continue
    2   continue
    3 continue

C      print 10, zval
C      10 format('GETMOM: zval = ',f10.6)

      if (zval .lt. 10*d1mach(3)) zval = 0d0

      end

