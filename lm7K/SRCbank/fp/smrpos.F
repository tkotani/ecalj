      subroutine smrpos(smrho,k1,k2,k3,n1,n2,n3)
C- Make smrho real and positive
C ----------------------------------------------------------------------
Ci Inputs
Ci   k1,k2,k3 dimensions smooth crystal densities, potentials on a mesh
Ci   n1,n2,n3 uniform mesh on which smrho,smcor,cmvxc defined
Co Outputs
Co   smrho :negative elements
Cr Remarks
Cu Updates
Cu   1 May 00 Adapted from nfp rhopos.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer k1,k2,k3,n1,n2,n3
      double complex smrho(k1,k2,k3)
C ... Local parameters
      integer stdo,lgunit,nneg,i3,i2,i1,iprint
      double precision rmin,xx

      stdo = lgunit(1)
      nneg = 0
      rmin = 999
      do  10  i3 = 1, n3
      do  10  i2 = 1, n2
      do  10  i1 = 1, n1
        xx = dble(smrho(i1,i2,i3))
        rmin = min(rmin,xx)
        if (xx .lt. 0) then
          nneg = nneg+1
          xx = 1d-8
        endif
        smrho(i1,i2,i3) = xx
   10 continue

      if (nneg .gt. 0 .and. iprint() .ge. 20) write(stdo,333) nneg,rmin
  333 format(' rhopos (warning): mesh density negative at',i6,
     .  ' points.  min=',f13.8)

      end

