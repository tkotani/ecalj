C Tests numerical integral of Ylm and Ylm**2
      subroutine fmain
      implicit none
      integer i,l,m,nx,np,nxp,ip,ix,nlm,lmax,lm,ilm
      double precision cy(1000),yl(100),xyz(3,10000),w(10000)
      double precision fx,fz,pi2,f,f1,f2,phi,flm(100)
      common /static/ xyz,w,cy

      call sylmnc(cy,30)
      print *, 'nx, np are no of points in cost and phi;'
      print *, 'lmax is max l for Ylm.  Enter nx, np, lmax: '
      nx = -122
      np = 0
      lmax = 6
      read(*,*) nx, np, lmax

      call fpiint(nx,np,nxp,xyz,w)
      print *, 'found ',nxp,' points'
C      do  5  i = 1, nxp
C    5 print 345, i, xyz(1,i), xyz(2,i), xyz(3,i), w(i)
C  345 format(i4, 4f15.10)
      do  20  l = 0, lmax
        do  20  m = -l, l
        nlm = l**2 + l + m + 1
        f1 = 0
        f2 = 0
        do  10  i = 1, nxp
          f = fz(xyz(1,i), l, nlm, cy)
          f1 = f1 + w(i)*f
          f2 = f2 + w(i)*f**2
   10   continue
        print 333, l, m, f1, f2
C       print 333, l, m, f1
  333   format(' l=',i2,', m=',i2,', num intgl f, f**2=',2f18.10)
   20 continue

      print *, 'numerically decompose Y(lm) into Y_L ... ... enter nlm'
      nlm = 16
      read(*,*) nlm
      do  40  lm = 1, nlm
        call dpzero(flm,nlm)
        do  42  i = 1, nxp
        do  42  ilm = 1, nlm
          flm(ilm) = flm(ilm) + 
     .      fz(xyz(1,i),l,ilm,cy)*fz(xyz(1,i),l,lm,cy)*w(i)
   42   continue
        write(*,928) lm, (flm(ilm), ilm = 1, nlm)
  928   format(i4,1x,25F8.5)
   40 continue


      end
      double precision function fz(xyz,l,nlm,cy)
      implicit none
      integer l,nlm
      double precision cost,phi,cy
      double precision sint,xyz(3),yl(200)

      call makylm(xyz,cy,l,yl)
      fz = yl(nlm)
      end
      double precision function fx(cost,phi,l,nlm,cy)
      implicit none
      integer l,nlm
      double precision cost,phi,cy
      double precision sint,xyz(3),yl(200)

      sint = dsqrt(1-cost**2)
      xyz(1)=dcos(phi)*sint
      xyz(2)=dsin(phi)*sint
      xyz(3)=cost
      call makylm(xyz,cy,l,yl)
      fx = yl(nlm)
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
