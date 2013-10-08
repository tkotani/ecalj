c-Taken from Ferdi's GW  -------------------------------------------------------------------
      subroutine genqbz (icase,qbas,n1,n2,n3,
     o qbz,wbz, nstbz)
c 92.02.16
c generates the k-points in the 1BZ
c the 1BZ is a parallepiped formed by G1,G2,G3 (qbas(3,3))
c this is divided into microcells defined by G1/n1,G2/n2,G3/n3
c the k-points may be thought of as being centred at each microcell
c the sampling weight for each k-point is the same (1/n1*n2*n3)

c qbas = base reciprocal vectors G1,G2,G3
c n1,n2,n3 = divisions along G1,G2,G3

c qbz  = k-points in the 1BZ
c wbz  = sampling weight for qbz

      implicit real*8 (a-h,o-z)
      dimension qbas(3,3)
      dimension qbz(3,n1*n2*n3),wbz(n1*n2*n3)
      dimension qmic(3,3),w1(3),w2(3),w3(3)
      integer(4):: icase,nstbz(n1*n2*n3)

c vectors forming microcells
c      call cv      (1.d0/dble(n1),qbas(1,1),3,qmic(1,1))
c      call cv      (1.d0/dble(n2),qbas(1,2),3,qmic(1,2))
c      call cv      (1.d0/dble(n3),qbas(1,3),3,qmic(1,3))
      qmic(:,1)= qbas(:,1)/dble(n1)
      qmic(:,2)= qbas(:,2)/dble(n2)
      qmic(:,3)= qbas(:,3)/dble(n3)
      nstbz=0
      half=0d0
      if(icase==2) half=.5d0
c sampling weight
      weight     = 1.d0/dble(n1*n2*n3)

      kount      = 0
      do      i1 = 1,n1
c      call cv      (dble(i1-1),qmic(1,1),3,w1)
        w1= (i1-1+half)*qmic(:,1)
        do      i2 = 1,n2
c      call cv      (dble(i2-1),qmic(1,2),3,w2)
          w2= (i2-1+half)*qmic(:,2)
          do      i3 = 1,n3
c      call cv      (dble(i3-1),qmic(1,3),3,w3)
            w3= (i3-1+half)*qmic(:,3)
            kount      = kount + 1
            qbz(1,kount) = w1(1) + w2(1) + w3(1)
            qbz(2,kount) = w1(2) + w2(2) + w3(2)
            qbz(3,kount) = w1(3) + w2(3) + w3(3)
            wbz(kount) = weight
            if(icase==2.and.(i1==1.or.i1==n1).and.(i2==1.or.i2==n2).and.(i3==1.or.i3==n3)) then
              nstbz(kount) = 2*2*2
            endif
          end do
        end do
      end do
Cstop2rx 2013.08.09 kino      if (kount .ne. n1*n2*n3)stop 'genqbz: wrong no. k-points'
      if (kount .ne. n1*n2*n3)call rx( 'genqbz: wrong no. k-points')
      return
      end
c--------------------
c      subroutine cv (c,v,n,
c     o w )
c forms w(i) = c * v(i)
c
c      implicit real*8(a-h,o-z)
c      dimension v(n)
c      dimension w(n)
c
c      do       i = 1,n
c      w(i)       = c*v(i)
c      end do
c
c      return
c      end
c----------------

c$$$c--- from lm6.11--------------------
c$$$      subroutine shorbz(p,pout,plat,qlat)
c$$$C- Shortens vector to equivalent in first Brillouin zone.
c$$$C ----------------------------------------------------------------
c$$$Ci Inputs:
c$$$Ci   plat,qlat lattice vectors and inverse
c$$$Ci   p         vector to shorten
c$$$Co Outputs:
c$$$Co   pout      shortened p
c$$$Cr Remarks
c$$$Cr   Switch around plat,qlat to shorten reciprocal space vectors.
c$$$Cr   Jan 1997 Adapted from shorps to fix bug:  Example:
c$$$Cr   plat=  -0.5  0.5  1.7517  0.5  -0.5  1.7517  0.5  0.5  -1.7517
c$$$Cr   p= 0.0 -0.5 -1.26384
c$$$Cr   Should get pout -> 0.5 0.0 0.48786, not -0.5 1.0 0.48786.
c$$$C ----------------------------------------------------------------
c$$$C     implicit none
c$$$      double precision p(3),pout(3),plat(3,3),qlat(3,3),x(3),x0,xx,a2,ap
c$$$      double precision tol
c$$$      parameter (tol=-1d-10)
c$$$      integer i,j,m,j2min,j3min,j1,j2,j3
c$$$
c$$$C --- Reduce to unit cell centered at origin ---
c$$$      do  1  i = 1, 3
c$$$C ... x is projection of pin along plat(i), with multiples of p removed
c$$$        x0 = p(1)*qlat(1,i)+p(2)*qlat(2,i)+p(3)*qlat(3,i)
c$$$        xx = idnint(x0)
c$$$        x(i) = x0-xx
c$$$    1 continue
c$$$C ... pout is x rotated back to Cartesian coordinates
c$$$      do  2  m = 1, 3
c$$$        pout(m) = x(1)*plat(m,1)+x(2)*plat(m,2)+x(3)*plat(m,3)
c$$$    2 continue
c$$$
c$$$C --- Try shortening by adding +/- basis vectors ---
c$$$   15 continue
c$$$      do  10  j1 =  0, 1
c$$$        j2min = -1
c$$$        if (j1 .eq. 0) j2min = 0
c$$$      do  10  j2 = j2min, 1
c$$$        j3min = -1
c$$$        if (j1 .eq. 0 .and. j2 .eq. 0) j3min = 0
c$$$      do  10  j3 = j3min, 1
c$$$
c$$$C     ... (-1,0,1) (plat(1) + (-1,0,1) plat(2)) + (-1,0,1) plat(3))
c$$$        do  17  i = 1, 3
c$$$          x(i) = plat(i,1)*j1 + plat(i,2)*j2 + plat(i,3)*j3
c$$$   17   continue
c$$$        a2 = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
c$$$        ap = pout(1)*x(1) + pout(2)*x(2) + pout(3)*x(3)
c$$$        j = 0
c$$$        if (a2 + 2*ap .lt. tol) j = 1
c$$$        if (a2 - 2*ap .lt. tol) j = -1
c$$$        if (j .ne. 0) then
c$$$          pout(1) = pout(1) + j*x(1)
c$$$          pout(2) = pout(2) + j*x(2)
c$$$          pout(3) = pout(3) + j*x(3)
c$$$          goto 15
c$$$        endif
c$$$   10 continue
c$$$
c$$$      end
c$$$C      subroutine fmain
c$$$C      implicit none
c$$$C      double precision plat(9),qlat(9),p(3),p1(3),xx
c$$$C      integer mode(3)
c$$$C
c$$$C      data plat /-0.5d0,0.5d0,1.7517d0,
c$$$C     .            0.5d0,-.5d0,1.7517d0,
c$$$C     .            0.5d0,0.5d0,-1.7517d0/
c$$$C
c$$$C      integer w(10000)
c$$$C      common /w/ w
c$$$C
c$$$C      data p /0.0d0,-0.5d0,-1.2638400000000001d0/
c$$$C
c$$$C      call wkinit(10000)
c$$$C
c$$$CC ... qlat = (plat^-1)^T so that qlat^T . plat = 1
c$$$C      call mkqlat(plat,qlat,xx)
c$$$C
c$$$C      call shorbz(p,p1,plat,qlat)
c$$$C      call prmx('p1 from shorbz',p1,1,1,3)
c$$$C
c$$$C      mode(1) = 2
c$$$C      mode(2) = 2
c$$$C      mode(3) = 3
c$$$C      call shorps(1,plat,mode,p,p1)
c$$$C      call prmx('p1 from shorps',p1,1,1,3)
c$$$C      end
