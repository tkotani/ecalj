      subroutine atsev(ic,nrc,ves,nl,nsp,qnu,sevat,sevs)
C- Accumulates atom contribution to sumev and shift q*V(rmax)
C ----------------------------------------------------------------------
Ci Inputs
Ci   ic    :class index
Ci   nrc   :number of classes of this class
Ci   ves   :potential at MT boundary for this class
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   qnu   :energy-weighted moments of the sphere charges
Ci   sevat :atom sum-of-eigenvalues
Co Outputs
Co   sevs  :sevs(1) = sum_sites sum-of-eigenvalues (KKR, atom)
Co         :sevs(2) = sum_sites (charge * ves)
Cl Local variables
Cl         :
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
      integer ic,nrc(ic),nl,nsp,i
      double precision ves(ic),qnu(3,nl*nsp,ic),sevat,sevs(2),q

      q = 0
      do  10  i = 1, nl*nsp
   10 q = q + qnu(1,i,ic)
      sevs(1) = sevs(1) + sevat*nrc(ic)
      sevs(2) = sevs(2) + q*ves(ic)*nrc(ic)

      end




