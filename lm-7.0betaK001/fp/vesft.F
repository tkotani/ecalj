      subroutine vesft(slat,ng,gv,kv,cv,k1,k2,k3,smrho,smpot,sum)
C- Make electrostatic potential of density given in recip space
C ----------------------------------------------------------------------
Ci Inputs
Ci   slat  :struct containing information about the lattice
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   kv    :indices for gather/scatter operations (gvlist.f)
Ci   cv    :work array holding smrho and smpot in glist form
Ci   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
Ci   smrho :FT of smooth density on uniform mesh
Co Outputs
Co   smpot :FT of smooth electrostatic potential
Co   sum   :integral pot*density
Cr Remarks
Cu Updates
Cu   22 Apr 00 Adapted from nfp ves_ft.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ng
      integer kv(ng,3),k1,k2,k3
      double precision slat(1),gv(ng,3),sum
      double complex smrho(k1,k2,k3),smpot(k1,k2,k3),cv(ng)
C ... Local parameters
      integer i
      double precision pi,pi8,alat,vol,tpiba,g2
      double complex ccc

      call tcn('vesft')
      call upack2('lat alat vol',slat,alat,vol)
      pi   = 4d0*datan(1d0)
      pi8  = 8*pi
      tpiba=2*pi/alat

C ... Gather density coefficients
      call gvgetf(ng,1,kv,k1,k2,k3,smrho,cv)

C ... smpot(G) = 8 pi /G**2 smrho(G)
C     call info2(30,1,0,' vesft:  smooth density coeff to (0,0,0) '//
C    .  'is %;12,6D',cv(1),0)
cccccccccccccc
c      print *,'ng=',ng
c      stop 'xxxxxxxxxx vesft'
cccccccccccccc
      sum = 0d0
      cv(1) = 0
      do  i = 2, ng
        g2 = tpiba*tpiba*(gv(i,1)**2+gv(i,2)**2+gv(i,3)**2)
        ccc = (pi8/g2)*cv(i)
        sum = sum + dconjg(cv(i))*ccc
        cv(i) = ccc
      enddo
      sum = vol*sum

C ... Scatter smooth potential into array
      call gvputf(ng,1,kv,k1,k2,k3,cv,smpot)
      call tcx('vesft')
      end

