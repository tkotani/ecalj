      module m_mkrou3
      contains
      subroutine mkrou3(mode,lmxa,nlml,nsp,pnz,dmatl,hab,sab,qsum,hsum)
C- l-decomposed charges and eigenvalue sum
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing
Ci         :1 make qsum only
Ci         :2 make qsum and hsum both
Ci   lmxa  :augmentation l-cutoff
Ci   dmatl :dmatl(l1,l2,mlm,i,j,isp) holds coefficients to a y_lm
Ci         :expansion of the function products f_i(l1) f_j(l2)
Ci         :where f_i,f_j form this table.
Ci         :    (uu  us  uz)
Ci         :    (su  ss  sz)
Ci         :    (zu  zs  zz)
Ci   hab   :<u,s | h | u,s> for each pair uu, us, su, ss; see remarks
Ci   sab   :<u,s | 1 | u,s>
Co Outputs
Co   qsum  :l-decomposed sphere charge
Co   hsum  :l-decomposed one-electron energy
Cr Remarks
Cr   qsum and hsum are used to find the band centers of gravity
Cr   u and s ae linear combinations radial wave functions defined as:
Cr   u has val=1, slo=1 at rmax, s has val=0, slo=1
Cu Updates
Cu   28 Aug 01 Extended to local orbitals.
C ----------------------------------------------------------------------
      implicit none
      integer mode,lmxa,n0,nab,nlml,nsp
      parameter (n0=10,nab=9)
      real(8),optional:: hab(nab,n0,nsp), hsum(n0,nsp)
      real(8):: qsum(n0,nsp) ,sab(nab,n0,nsp),
     .dmatl(0:lmxa,0:lmxa,nlml,3,3,nsp),pnz(n0,2)
      integer l,isp,m
      double precision pi,srfpi,qz,hz
      if (mode .eq. 0) return
C     call tcn('mkrou3')
      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)
      do  isp = 1, nsp
        do  l = 0, lmxa
          if (mode .ge. 1) qsum(l+1,isp) = 0d0
          if (mode .ge. 2) hsum(l+1,isp) = 0d0
        enddo
      enddo
      do  isp = 1, nsp
        do  l = 0, lmxa
          m = l+1
          qsum(m,isp) =
     .    + dmatl(l,l,1,1,1,isp)*sab(1,m,isp)*srfpi
     .    + dmatl(l,l,1,1,2,isp)*sab(2,m,isp)*srfpi
     .    + dmatl(l,l,1,2,1,isp)*sab(3,m,isp)*srfpi
     .    + dmatl(l,l,1,2,2,isp)*sab(4,m,isp)*srfpi
          if (mode .ge. 2) then
            hsum(m,isp) =
     .      + dmatl(l,l,1,1,1,isp)*hab(1,m,isp)*srfpi
     .      + dmatl(l,l,1,1,2,isp)*hab(2,m,isp)*srfpi
     .      + dmatl(l,l,1,2,1,isp)*hab(3,m,isp)*srfpi
     .      + dmatl(l,l,1,2,2,isp)*hab(4,m,isp)*srfpi
          endif
C         ... uz, sz, zu, zs, zz terms
          if (pnz(m,1) .ne. 0) then
            qz =
     .      + dmatl(l,l,1,1,3,isp)*sab(5,m,isp)*srfpi
     .      + dmatl(l,l,1,2,3,isp)*sab(6,m,isp)*srfpi
     .      + dmatl(l,l,1,3,1,isp)*sab(5,m,isp)*srfpi
     .      + dmatl(l,l,1,3,2,isp)*sab(6,m,isp)*srfpi
     .      + dmatl(l,l,1,3,3,isp)*sab(7,m,isp)*srfpi
! 13sep2012takao add qz and hz
            qsum(m,isp) = qsum(m,isp) + qz  !qsum is including local orbital 13sep2012takao
            if (mode .ge. 2) then
              hz =
     .        + dmatl(l,l,1,1,3,isp)*hab(5,m,isp)*srfpi
     .        + dmatl(l,l,1,2,3,isp)*hab(6,m,isp)*srfpi
     .        + dmatl(l,l,1,3,1,isp)*hab(5,m,isp)*srfpi
     .        + dmatl(l,l,1,3,2,isp)*hab(6,m,isp)*srfpi
     .        + dmatl(l,l,1,3,3,isp)*hab(7,m,isp)*srfpi
! 13sep2012takao add qz and hz
              hsum(m,isp)=hsum(m,isp)+hz  !hsum is including local orbital 13sep2012takao     
            endif
          endif

        enddo
      enddo
C     call tcx('mkrou3')
      end subroutine mkrou3
      end module m_mkrou3
