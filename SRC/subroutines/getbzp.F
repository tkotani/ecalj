      subroutine getbzp(nfilqp,npln,nwmx,nqp,nw,e,nq,vx,vy,x,y,z,q)
C- Get BZ planes and BZ sampling parameters for making BZ maps
C ----------------------------------------------------------------------
Ci Inputs
Ci   nfilqp: file-handle for BZPL
Ci   npln: total number of BZ planes
Ci   nwmx: maximum number of energy windows
Ci   nqp: total number of q-points
Co Outputs
Co   nw: number of energy windows for each BZ plane
Co   e(1), e(2): lower and upper limits of energy windows
Co   nq,vx,vy,x,y,z: see Remarks
Co   q: q-points
Cr Remarks
Cr   The form of the free-format file BZPL is:
Cr      nw1  np1
Cr      e(1,1)  e(2,1)
Cr      ...
Cr      e(1,nw1) e(2,nw1)
Cr      (vx(i),i=1,3) x(1) x(2) nq(1)  (vy(i),i=1,3) y(1) y(2) nq(2)  z
Cr      ... (np1 BZ planes) ...
Cr      nw2  np2
Cr      ...
Cr
Cr   The vector vy is first made orthogonal to vx (i.e. a rectangular BZ
Cr   slice) and then both are normalized.  The integers nq(1) and nq(2)
Cr   are the number of divisions along the orthonormalized vectors
Cr   vx and vy, respectively.  The numbers x(1), x(2), y(1), y(2), and z
Cr   are in units of (2 pi / alat) (as in file SYML).  The nq(1)*nq(2)
Cr   q-points are constructed as follows:
Cr      q(i,1)           = x(1)*vx(i) + y(1)*vy(i) + z*vn(i)
Cr      ...
Cr      q(i,nq(1))       = x(2)*vx(i) + y(1)*vy(i) + z*vn(i)
Cr      ...
Cr      q(i,nq(1)*nq(2)) = x(2)*vx(i) + y(2)*vy(i) + z*vn(i)
Cr
Cr   where vn is the unit normal vector obtained from the cross product
Cr   of the orthonormalized vx and vy vectors.  The file BZPL can
Cr   contain input lines corresponding to multiple sets of BZ
Cr   planes/energy windows; the resulting q-points are numbered
Cr   sequentially.
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nfilqp,npln,nwmx,nqp
      integer nw(npln),nq(2,npln)
      double precision e(2,nwmx,npln),vx(3,npln),vy(3,npln),x(2,npln),
     .y(2,npln),z(npln),q(3,nqp)
C Local variables
      integer i,iprint,ip,nptot,nqtot,nwin,np0,iw,nq0,iq,ix,iy
      double precision rx,ry,dp,ddot,dqx,dqy,qx,qy
      double precision vx0(3),vy0(3),vn(3)

      if (iprint() .gt. 30) then
        print *
        print *, ' GETBZP: Reading BZ planes from disc ...'
      endif

      nptot = 0
      nqtot = 0
      rewind nfilqp

C --- Begin loop over sets of BZ planes/energy windows ---
   10 read (nfilqp,*,err=999,end=999) nwin,np0
      call rxx(nwin .lt. 1,'GETBZP: bad nwin')
      call rxx(nwin .gt. nwmx,'GETBZP: nwin gt nwmx')
      call rxx(nptot+np0 .gt. npln,'GETBZP: nptot gt npln')
      do  20  iw = 1, nwin
        read (nfilqp,*,err=999,end=999) (e(i,iw,nptot+1),i=1,2)
        call rxx(e(1,iw,nptot+1).gt.e(2,iw,nptot+1),'GETBZP: e1 gt e2')
   20 continue

C --- Begin loop over BZ planes in this set ---
      do  80  ip = nptot+1, nptot+np0
        nw(ip) = nwin
        if (ip .gt. nptot+1)
     .  call dcopy(2*nwin,e(1,1,nptot+1),1,e(1,1,ip),1)
        read (nfilqp,*,err=999,end=999)
     .  (vx(i,ip),i=1,3),(x(i,ip),i=1,2),nq(1,ip),
     .  (vy(i,ip),i=1,3),(y(i,ip),i=1,2),nq(2,ip),z(ip)
        call rxx(nq(1,ip) .lt. 2,'GETBZP: bad n1')
        call rxx(nq(2,ip) .lt. 2,'GETBZP: bad n2')
        nq0 = nq(1,ip)*nq(2,ip)
        nqtot = nqtot + nq0
        call rxx(nqtot .gt. nqp,'GETBZP: nqtot gt nqp')

C --- Orthonormalize axes, construct normal vector ---
        rx = ddot(3,vx(1,ip),1,vx(1,ip),1)
        call rxx(dabs(rx) .lt. 1.d-06,'GETBZP: bad vx')
        dp = dabs(ddot(3,vx(1,ip),1,vy(1,ip),1)) / rx

        do  30  i = 1, 3
          vy0(i) = vy(i,ip) - dp*vx(i,ip)
   30   continue

        rx = dsqrt(rx)
        ry = dsqrt(ddot(3,vy0,1,vy0,1))
        call rxx(dabs(ry) .lt. 1.d-06,'GETBZP: vx and vy are colinear')

        do  40  i = 1, 3
          vx0(i) = vx(i,ip) / rx
          vy0(i) = vy0(i) / ry
   40   continue

        call cross(vx0,vy0,vn)

        if (iprint() .ge. 40)
     .  write (*,500) ip,nq0,vx0,(x(i,ip),i=1,2),nq(1,ip),
     .  vy0,(y(i,ip),i=1,2),nq(2,ip),vn,z(ip)
        if (iprint() .gt. 50) write(*,510)
  500   format(/' Plane ',i2,'.',i7,' k-points:'/
     .  ' VX= (',3f10.5,' )  X1,X2=',2f10.5,'  NX=',i5/
     .  ' VY= (',3f10.5,' )  Y1,Y2=',2f10.5,'  NY=',i5/
     .  ' VN= (',3f10.5,' )      Z=',f10.5)
  510   format(/19x,'K')

C --- Construct q-points ---
        dqx = (x(2,ip) - x(1,ip)) / (nq(1,ip) - 1d0)
        dqy = (y(2,ip) - y(1,ip)) / (nq(2,ip) - 1d0)
        iq = nqtot - nq0
        do  70  iy = 1, nq(2,ip)
          qy = y(1,ip) + (iy-1)*dqy
          do  60  ix = 1, nq(1,ip)
            qx = x(1,ip) + (ix-1)*dqx
            iq = iq + 1
            do  50  i = 1, 3
              q(i,iq) = qx*vx0(i) + qy*vy0(i) + z(ip)*vn(i)
              if (dabs(q(i,iq)) .lt. 1d-12) q(i,iq) = 0d0
   50       continue
            if (iprint() .gt. 50) write (*,520) iq, (q(i,iq),i=1,3)
  520       format(i5,3f10.6)
   60     continue
   70   continue
   80 continue

C --- Get next set ---
      nptot = nptot + np0
      goto 10
  999 call rxx(nqtot .eq. 0,'GETBZP: empty or badly formed BZPL file')
      call rxx(nptot .ne. npln,'GETBZP: file mismatch, nptot ne npln')
      call rxx(nqtot .ne. nqp,'GETBZP: file mismatch, nqtot ne nqp')

      end

