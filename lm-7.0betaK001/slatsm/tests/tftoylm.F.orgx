C Extracts coefficients a_L of function f_L(cost,phi) = sum_L a_L Y_L
C using Legendre gaussian quadrature

C Because any Y_L Y_L' can be integrated exactly provided
C nx, np are at least L+L'+1, all coefficents a_L of a function
C f = sum_L<=L' a_L Y_L can be obtained exactly.
      subroutine fmain
C     program test
      integer l,m,nx,np,ip,ix,nlm,lmax
      double precision cy(1000),YL(100),xyz(3),aL(100)
      double precision xtab(100),w(100),wphi
      double precision pi2,f,phi,sint,cosphi,sinphi,wgt

      call sylmnc(cy,20)
      pi2 = 8*datan(1d0)
      print *, 'Coeffs determined exactly if nx and np > L+L'' '
      print *,
     .  'See subroutine fx for function def and L'' truncation ...' 
      print * 

    5 print *, 'nx, np are no of points in cost and phi;'
      print *, 'L is max l for coefficient a_L.  Enter nx, np, L: '
      read(*,*) nx, np, lmax
      call mklegw(nx,xtab,w,0)

      do  7  nlm = 1, (lmax+1)**2
    7 aL(nlm) = 0

      do  10  ip = 1, np
        phi = pi2*dble(ip)/np
        wphi = pi2/np
        cosphi = dcos(phi)
        sinphi = dsin(phi)
        do  10  ix = 1, nx
          sint = dsqrt(1-xtab(ix)**2)
          xyz(1) = cosphi*sint
          xyz(2) = sinphi*sint
          xyz(3) = xtab(ix)
          call makylm(xyz,cy,lmax,YL)
          call fx(xyz,cy,f)
          wgt = wphi*w(ix)*f
          do  20  l = 0, lmax
            do  20  m = -l, l
              nlm = l**2 + l + m + 1
              aL(nlm) = aL(nlm) + wgt*YL(nlm)
   20     continue
   10 continue

      do  30  l = 0, lmax
        do  30  m = -l, l
        nlm = l**2 + l + m + 1
        print 333, l, m, aL(nlm), 1/aL(nlm)
  333   format(' l=',i2,', m=',i2,', aL,1/aL=',2f15.10)
   30 continue

      goto 5
      end

      subroutine fx(xyz,cy,fL)
C- Returns f(xyz)
      double precision xyz,cy(1),fL
C Local variables
      integer ilm,lmax,nlm
      double precision aL2(100),YL(200)

      lmax = 6
      nlm = (lmax+1)**2
      call makylm(xyz,cy,lmax,yl)

c aL = 1 / (2 + nlm)

      do  10  ilm = 1, nlm
   10 aL2(ilm) = 1/(2+dble(ilm))

      fL = 0
      do  20  ilm = 1, nlm
   20 fL = fL + aL2(ilm)*YL(ilm)
        
      end
      subroutine makylm(p,cy,lmax,yl)
C- Make proper ylm for point p
      implicit none
      integer lmax,ilm,nlm
      double precision r2s,p(3),cy(1),yl(1)
      call sylm(p,yl,lmax,r2s)
      nlm = (lmax+1)**2
      do 2 ilm=1,nlm
    2 yl(ilm)=yl(ilm)*cy(ilm)
      end
