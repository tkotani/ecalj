#define VANILLA 1
#define NOQUAD 1
CSFCPP#define F90 1
#define FFTW 1
#define IN_PLACE 1
      subroutine ropylg(lp,lmax,ndim,nrx,nr,x,y,z,r2,yl,gyl)
C- Gradients of YL's (polynomials) for a set of points, with YL as input
C ----------------------------------------------------------------------
Ci Inputs
Ci   lp    :if nonzero, adds term  r^l grad (r^-l Yl).
Ci   lmax  :maximum l for a given site
Ci   ndim  :dimensions gyl.  Must be at least (lmax+1)**2
Ci   nrx   :leading dimension of yl,gyl
Ci   nr    :number of points
Ci   x,y,z :cartesian coordinates of points
Ci   r2    :x^2+y^2+z^2
Ci   yl    :Spherical harmonic polynomials YL.  YL's must be normalized
Ci         :and must be made through lmax+1 (i.e. nlm=1..(lmax+2)**2)
Co Outputs
Co   gyl   :gradient of yl
Cl Local variables
Cl         :
Cr Remarks
Cr   
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lp,lmax,ndim,nrx,nr
      double precision x(nr),y(nr),z(nr),yl(nrx,*),gyl(nrx,ndim,3),r2(*)
C ... Local parameters
      integer ilm,kx1,kx2,ky1,ky2,kz,l,m,i
      double precision cx1,cx2,cy1,cy2,cz,f

      if ((lmax+1)**2 .gt. ndim) call rx('ropylg: ndim too small')

C --- Gradients of yl's ---
      ilm = 0
      do  10  l = 0, lmax
      do  10  m = -l, l
      ilm = ilm+1
      call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
      do  10  i = 1, nr
      f = (2*l+1)/r2(i)
      gyl(i,ilm,1) = (yl(i,ilm)*x(i) - cx1*yl(i,kx1) - cx2*yl(i,kx2))*f
      gyl(i,ilm,2) = (yl(i,ilm)*y(i) - cy1*yl(i,ky1) - cy2*yl(i,ky2))*f
      gyl(i,ilm,3) = (yl(i,ilm)*z(i) - cz*yl(i,kz))*f
   10 continue

      if (lp .eq. 0) return

C --- Add r**l (grad r**-l) yl ---
      ilm = 0
      do  20  l = 0, lmax
      do  20  m = -l, l
      ilm = ilm+1
      do  20  i = 1, nr
      gyl(i,ilm,1) = gyl(i,ilm,1) - l*x(i)/r2(i)*yl(i,ilm)
      gyl(i,ilm,2) = gyl(i,ilm,2) - l*y(i)/r2(i)*yl(i,ilm)
      gyl(i,ilm,3) = gyl(i,ilm,3) - l*z(i)/r2(i)*yl(i,ilm)
   20 continue

      end
C$$$#if TEST
C$$$      subroutine tl(np,ilmx,nlm,rp,grp,ggrp,frl,yl,gyl,wp,xp,yp,zp,r2)
C$$$      implicit none
C$$$      integer np,ilmx,nlm
C$$$      double precision yl(np,nlm),gyl(np,nlm,3),xp(np),yp(np),zp(np),
C$$$     .  wp(np),rp(np),grp(np,3),ggrp(np,3),frl(nlm),r2(np)
C$$$      integer ll,ip,lmax,j,ilm,l

C$$$      lmax = ll(nlm)
C$$$      call ropyln(np,xp,yp,zp,lmax+1,np,yl,r2)
C$$$      call xyl(nlm,np,wp,yl(1,ilmx),yl,frl)
C$$$      call prmr('yl(ilmx)',frl,nlm)

C$$$C ... Show laplacian yl is 0:
C$$$      call ropylg(0,lmax,nlm,np,np,xp,yp,zp,r2,yl,gyl)
C$$$      call dpzero(ggrp,np*3)
C$$$      do  120  j = 1, 3
C$$$        call xyl(nlm,np,wp,gyl(1,ilmx,j),yl,frl)
C$$$        call prmr('component of grad',frl,nlm)
C$$$        do  118  ilm = 1, nlm
C$$$        do  118  ip = 1, np
C$$$  118   ggrp(ip,j) = ggrp(ip,j) + frl(ilm)*gyl(ip,ilm,j)
C$$$        call xyl(nlm,np,wp,ggrp(1,j),yl,frl)
C$$$        call prmr('component of nabla',frl,nlm)
C$$$  120 continue
C$$$      call dpadd(ggrp,ggrp(1,2),1,np,1d0)
C$$$      call dpadd(ggrp,ggrp(1,3),1,np,1d0)
C$$$      call xyl(nlm,np,wp,ggrp,yl,frl)
C$$$      call prmr('laplacian',frl,nlm)
C$$$C     call cexit(1,1)

C$$$C ... Show laplacian r^-l yl is -l(l+1) yl; 
C$$$C     Use grad(r^-1 yl) = 1/r sum_L (a_L yl)
C$$$      call ropylg(1,lmax,nlm,np,np,xp,yp,zp,r2,yl,gyl)
C$$$C ... Term 1/r grad r (grad r^-l yl)
C$$$      call dpzero(ggrp,np*3)
C$$$      do  20  j = 1, 3
C$$$        call xyl(nlm,np,wp,gyl(1,ilmx,j),yl,frl)
C$$$        call prmr('component of grad',frl,nlm)
C$$$        do  18  ilm = 1, nlm
C$$$        do  18  ip = 1, np
C$$$   18   ggrp(ip,j) = ggrp(ip,j) + frl(ilm)*gyl(ip,ilm,j)
C$$$        call xyl(nlm,np,wp,ggrp(1,j),yl,frl)
C$$$C        call prmr('1st term of nabla',frl,nlm)
C$$$C ... Term grad (1/r) . (grad r^-l yl)
C$$$        do  22  ip = 1, np
C$$$        if (j.eq.1) ggrp(ip,j) = ggrp(ip,j) - xp(ip)*gyl(ip,ilmx,j)
C$$$        if (j.eq.2) ggrp(ip,j) = ggrp(ip,j) - yp(ip)*gyl(ip,ilmx,j)
C$$$        if (j.eq.3) ggrp(ip,j) = ggrp(ip,j) - zp(ip)*gyl(ip,ilmx,j)
C$$$   22   continue
C$$$        call xyl(nlm,np,wp,ggrp(1,j),yl,frl)
C$$$C        call prmr('component of nabla',frl,nlm)
C$$$   20 continue
C$$$      call dpadd(ggrp,ggrp(1,2),1,np,1d0)
C$$$      call dpadd(ggrp,ggrp(1,3),1,np,1d0)
C$$$      call xyl(nlm,np,wp,ggrp,yl,frl)
C$$$      call prmr('laplacian',frl,nlm)
C$$$      call cexit(1,1)
C$$$      end
C$$$      subroutine xyl(nlm,np,wp,fp,yl,fl)
C$$$C- Yl-projection of function tabulated on an angular mesh
C$$$      implicit none
C$$$      integer nlm,np,ip,ilm
C$$$      double precision fl(nlm),fp(np),yl(np,nlm),wp(np)

C$$$      call dpzero(fl,nlm)
C$$$      do  20  ip = 1, np
C$$$      do  20  ilm = 1, nlm
C$$$   20 fl(ilm) = fl(ilm) + fp(ip)*wp(ip)*yl(ip,ilm)
C$$$      end
C$$$      subroutine prmr(strn,f,nl)
C$$$      implicit none
C$$$      integer nl,j,fopna,ifi
C$$$      double precision f(nl)
C$$$      character*(10) fmt, strn*(*)
C$$$      ifi = 19
C$$$      open(ifi,file='out')
C$$$      write(ifi,*) nl, 2
C$$$      do  10  j = 1, nl
C$$$   10 write(ifi,333) j, f(j)
C$$$  333 format(i4, f15.10)
C$$$      close(ifi)
C$$$      print *, strn
C$$$      pause
C$$$      end

C$$$C Test program to check ropylg
C$$$      subroutine fmain
C$$$      implicit none
C$$$      integer nrx,lmx,nlmx,nlm2,nr,lmax,nlm1,ir,ii,i,l,ilm,i1,i2,nsize,
C$$$     .  nnn
C$$$      parameter (nrx=20,lmx=6,nlmx=(lmx+1)**2,nlm2=(lmx+2)**2,
C$$$     .  nsize=100000, nnn=300)
C$$$      double precision cy(16**2),x(nrx),y(nrx),z(nrx),r2(nrx),tops,
C$$$     .  ylv(nrx,nlm2),gylv(nrx,nlmx,3),yl(nlm2),gyl(nlm2,3),dr(3),
C$$$     .  p(3,nnn),wp(nnn)
C$$$      integer oxp,oyp,ozp,or2,oyl,ogyl,orp,ogrp,oggrp,ofrl,ll,lp,np,
C$$$     .  nlmf,nph,nth,w(nsize)
C$$$      common /w/ w
C$$$      common /static/ cy

C$$$      call wkinit(nsize)
C$$$      call sylmnc(cy,15)

C$$$C --- Laplacian of Yl ---
C$$$      print *, 'ilm:'
C$$$      ilm = 8
C$$$      read(*,*) ilm
C$$$      lmax = ll(ilm)+1
C$$$      nth=lmax+2
C$$$      nph=2*nth
C$$$      nlmf = (lmax+2)**2
C$$$      call fpiint(nth,nph,np,p,wp)
C$$$      print *, np, ' angular points'

C$$$      call defrr(oxp,     np)
C$$$      call defrr(oyp,     np)
C$$$      call defrr(ozp,     np)
C$$$      call defrr(or2,     np)
C$$$      call defrr(oyl,     (lmax+3)**2*np)
C$$$      call defrr(ofrl,    nlmf)
C$$$      call dcopy(np,p(1,1),3,w(oxp),1)
C$$$      call dcopy(np,p(2,1),3,w(oyp),1)
C$$$      call dcopy(np,p(3,1),3,w(ozp),1)
C$$$      call defrr(ogyl,    nlmf*np*3)
C$$$      call defrr(orp,     np)
C$$$      call defrr(ogrp,    np*3)
C$$$      call defrr(oggrp,   np*3)
C$$$      call tl(np,ilm,nlmf,w(orp),w(ogrp),w(oggrp),w(ofrl),
C$$$     .  w(oyl),w(ogyl),wp,w(oxp),w(oyp),w(ozp),w(or2))

C$$$C --- Compare ropylg against ylg ---
C$$$      nr = 5
C$$$   99 print *, 'lmax='
C$$$      read(*,*) lmax
C$$$      if (lmax .gt. lmx) stop 'increase lmx in main'

C$$$      call makr(0d0,nr,x,y,z)

C$$$C ... Make grad ylm's
C$$$      nlm1 = (lmax+2)**2
C$$$      call ropyln(nr,x,y,z,lmax+1,nrx,ylv,r2)
C$$$      call ropylg(0,lmax,nlmx,nrx,nr,x,y,z,r2,ylv,gylv)

C$$$C ... Check against ylg
C$$$      tops = 0d0
C$$$      do  10  ir = 1, nr
C$$$        dr(1) = x(ir)
C$$$        dr(2) = y(ir)
C$$$        dr(3) = z(ir)
C$$$        call ylg(dr,lmax,nlm2,cy,yl,gyl)

C$$$        do  11  i = 1, 3

C$$$        do  12  l = 0, lmax
C$$$        i1 = l*l+1
C$$$        i2 = (l+1)**2
C$$$   12   print 333, (gyl(ii,i),ii=i1,i2)
C$$$  333   format(9f8.5)
C$$$        print *
C$$$        do  14  l = 0, lmax
C$$$        i1 = l*l+1
C$$$        i2 = (l+1)**2
C$$$        do  16  ii = i1, i2
C$$$   16   tops = max(tops,dabs(gyl(ii,i)-gylv(ir,ii,i)))
C$$$   14   print 333, (gylv(ir,ii,i)-gyl(ii,i),ii=i1,i2)
C$$$        print *, '----------- end of ir,i=', ir,i
C$$$   11   continue
C$$$   10 continue

C$$$      print 335, tops
C$$$  335 format(' max errors for grad h:',f12.6)

C$$$      end
C$$$      subroutine makr(rsm,nr,x,y,z)
C$$$      implicit none
C$$$      integer nr,i,ir
C$$$      double precision rs,rsm,x(1),y(1),z(1)
C$$$      real ran1
C$$$      rs = rsm
C$$$      if (rsm .lt. 1d-9) rs = .5d0
C$$$      call ran1in(1)
C$$$      do  10  i = 1, nr
C$$$        ir = i+1
C$$$        x(i) = abs((ran1()-.5d0)*5*rs)
C$$$        y(i) = (ran1()-.5d0)*5*rs
C$$$        z(i) = (ran1()-.5d0)*5*rs
C$$$   10 continue

C$$$      x(1) = .3d0*dsqrt(2d0)
C$$$      y(1) = .4d0*dsqrt(2d0)
C$$$      z(1) = .5d0*dsqrt(2d0)
C$$$      end
C$$$#endif

