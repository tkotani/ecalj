      subroutine shortn(p,p1,dlat,nkd)
C- Get p1 = shortest vector such that p1-p is a lattice vector.
C ----------------------------------------------------------------
Ci Inputs
Ci   p     :vector to shorten
Ci   dlat  :lattice vectors, sorted by increasing length
Ci   nkd   :number of dlat
Co Outputs
Co   p1    :shortened vector
Cr Remarks
Cr   A slightly skewed norm is used to make result unique.
Cr   The first vector in the list must be the zero vector.
c ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nkd
      double precision p(3),p1(3),dlat(3,nkd)
C ... Local parameters
      integer irep,k0,k
      double precision anrm2,x,y,z,p2,critk0,dd,crit
      anrm2(x,y,z) = x*x*1.00001d0 + y*y*1.00002d0 + z*z*1.00003d0
     .                -x*0.000004d0 - y*0.000003d0 - z*0.000002d0

      p1(1) = p(1)
      p1(2) = p(2)
      p1(3) = p(3)
      do  88  irep = 1, 20
        p2 = anrm2(p1(1),p1(2),p1(3))
        k0 = 1
        critk0 = p2
        do  52  k = 1, nkd
          dd = dlat(1,k)**2 + dlat(2,k)**2 + dlat(3,k)**2
C         Only valid if dlat is sorted:
          if (dd .gt. p2*4d0) goto 53
           crit = anrm2(p1(1)+dlat(1,k),p1(2)+dlat(2,k),p1(3)+dlat(3,k))
          if (crit .lt. critk0) then
            k0 = k
            critk0 = crit
          endif
  52    continue
  53    if (k0 .eq. 1) return
        p1(1) = p1(1) + dlat(1,k0)
        p1(2) = p1(2) + dlat(2,k0)
        p1(3) = p1(3) + dlat(3,k0)
   88 continue
      call rx('shortn: shortest vector not found')
      end

