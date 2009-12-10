      subroutine rcnsl0(tau,a,lmax,alat,rlat,nkr,dlat,nkd,vol,cy,dl)
C  reduced structure constants on lattice for e=0 and q=0.
C  result is periodic sum of (2*l-1)!!*ylm/r**(l+1). additive
C  constant for l=0 is chosen so that function averages to zero.
C     implicit none
C passed parameters 
      integer lmax,nkd,nkr
      double precision a,alat,vol
      double precision tau(3),dlat(3,nkd),rlat(3,nkr),cy(1),dl(1)
C local parameters 
      integer ilm,nlm
      external rcsl01,rcsl02

      call tcn('k-space Ewald')
      call rcsl01(tau,a,lmax,alat,rlat,nkr,vol,dl)
      call tcx('k-space Ewald')
      call tcn('r-space Ewald')
      call rcsl02(tau,a,lmax,alat,dlat,nkd,dl)
      call tcx('r-space Ewald')
      nlm = (lmax+1)**2
      do  ilm = 1, nlm
        dl(ilm) = dl(ilm)*cy(ilm)
      enddo
      return
      end
      subroutine rcsl01(tau,a,lmax,alat,rlat,nkr,vol,dl)
C  k-space part of reduced structure constants for e=0 and q=0
C     implicit none
C passed parameters 
      integer lmax,nkr
      double precision a,alat,r2,vol
      double precision tau(3),rlat(3,nkr),dl(1)
C local parameters 
      integer ilm,ir,l,m,nlm,lmaxx
      parameter(lmaxx=17)
      double precision fpibv,gamma,scalp,tpi,tpiba,yyy,eiphi(0:1),xxx
      double precision r(3),yl(lmaxx**2)
      external sylm

      if (lmax .gt. lmaxx) call rx('rcnsl0: increase lmaxx')
      tpi=8d0*datan(1d0)
      gamma = .25d0/(a*a)
      fpibv = 2d0*tpi/vol
      tpiba = tpi/alat
      nlm = (lmax+1)**2
      do  ilm = 1, nlm
        dl(ilm) = 0d0
      enddo
      dl(1) = -fpibv*gamma
      do  ir = 2, nkr
        r(1) = tpiba*rlat(1,ir)
        r(2) = tpiba*rlat(2,ir)
        r(3) = tpiba*rlat(3,ir)
        scalp = alat*(r(1)*tau(1) + r(2)*tau(2) + r(3)*tau(3))
        eiphi(0) = dcos(scalp)
        eiphi(1) = dsin(scalp)
        call sylm(r,yl,lmax,r2)
        yyy = fpibv*dexp(-gamma*r2)/r2
        ilm = 0
        do  l = 0, lmax
          do  m = 1, 2*l+1
            ilm = ilm+1
            dl(ilm) = dl(ilm) + yl(ilm)*yyy*eiphi(0)
          enddo
C eiphi *= (0,-1):
          xxx = eiphi(0)
          eiphi(0) = eiphi(1)
          eiphi(1) = -xxx
        enddo
      enddo

      end
      subroutine rcsl02(tau,a,lmax,alat,dlat,nkd,dl)
C  real space summation
C     implicit none
C passed parameters 
      integer lmax,nkd
      double precision a,alat,r2
      double precision tau(3),dlat(3,nkd),dl(1)
C local parameters 
      integer ilm,ir,ir1,l,m,lmaxx
      parameter(lmaxx=17)
      double precision a2,cc,gl,r1,srpi,ta2,derfc
      double precision r(3),yl(lmaxx**2),chi(0:lmaxx)
      external sylm

      srpi=dsqrt(4d0*datan(1d0))
      a2 = a*a
      ta2 = 2d0*a2
      cc = 4d0*a2*a/srpi
      ir1 = 2
      if (tau(1)**2+tau(2)**2+tau(3)**2 .gt. 1d-6) ir1=1
      do  ir = ir1, nkd
        r(1) = alat*(tau(1)-dlat(1,ir))
        r(2) = alat*(tau(2)-dlat(2,ir))
        r(3) = alat*(tau(3)-dlat(3,ir))
        call sylm(r,yl,lmax,r2)
        r1 = dsqrt(r2)
        chi(0) = derfc(a*r1)/r1
        gl = -cc*dexp(-a2*r2)/ta2
        do  l = 1, lmax
          chi(l) = ((2*l-1)*chi(l-1) - gl)/r2
          gl = ta2*gl
        enddo
        ilm = 0
        do  l = 0, lmax
          do  m = 1, 2*l+1
            ilm = ilm+1
            dl(ilm) = dl(ilm) + yl(ilm)*chi(l)
          enddo
        enddo
      enddo
C --- add dl3 for diagonal sructure constants ------
      if (ir1 .eq. 2) dl(1) = dl(1) - 2d0*a/srpi
      return
      end
