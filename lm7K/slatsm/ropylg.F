#define VANILLA 1
#define NOQUAD 1
#define F90 1
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
      double precision x(nr),y(nr),z(nr),yl(nrx,1),gyl(nrx,ndim,3),r2(1)
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
#if TEST
      subroutine tl(np,ilmx,nlm,rp,grp,ggrp,frl,yl,gyl,wp,xp,yp,zp,r2)
      implicit none
      integer np,ilmx,nlm
      double precision yl(np,nlm),gyl(np,nlm,3),xp(np),yp(np),zp(np),
     .  wp(np),rp(np),grp(np,3),ggrp(np,3),frl(nlm),r2(np)
      integer ll,ip,lmax,j,ilm,l

      lmax = ll(nlm)
      call ropyln(np,xp,yp,zp,lmax+1,np,yl,r2)
      call xyl(nlm,np,wp,yl(1,ilmx),yl,frl)
      call prmr('yl(ilmx)',frl,nlm)

C ... Show laplacian yl is 0:
      call ropylg(0,lmax,nlm,np,np,xp,yp,zp,r2,yl,gyl)
      call dpzero(ggrp,np*3)
      do  120  j = 1, 3
        call xyl(nlm,np,wp,gyl(1,ilmx,j),yl,frl)
        call prmr('component of grad',frl,nlm)
        do  118  ilm = 1, nlm
        do  118  ip = 1, np
  118   ggrp(ip,j) = ggrp(ip,j) + frl(ilm)*gyl(ip,ilm,j)
        call xyl(nlm,np,wp,ggrp(1,j),yl,frl)
        call prmr('component of nabla',frl,nlm)
  120 continue
      call dpadd(ggrp,ggrp(1,2),1,np,1d0)
      call dpadd(ggrp,ggrp(1,3),1,np,1d0)
      call xyl(nlm,np,wp,ggrp,yl,frl)
      call prmr('laplacian',frl,nlm)
C     call cexit(1,1)

C ... Show laplacian r^-l yl is -l(l+1) yl; 
C     Use grad(r^-1 yl) = 1/r sum_L (a_L yl)
      call ropylg(1,lmax,nlm,np,np,xp,yp,zp,r2,yl,gyl)
C ... Term 1/r grad r (grad r^-l yl)
      call dpzero(ggrp,np*3)
      do  20  j = 1, 3
        call xyl(nlm,np,wp,gyl(1,ilmx,j),yl,frl)
        call prmr('component of grad',frl,nlm)
        do  18  ilm = 1, nlm
        do  18  ip = 1, np
   18   ggrp(ip,j) = ggrp(ip,j) + frl(ilm)*gyl(ip,ilm,j)
        call xyl(nlm,np,wp,ggrp(1,j),yl,frl)
C        call prmr('1st term of nabla',frl,nlm)
C ... Term grad (1/r) . (grad r^-l yl)
        do  22  ip = 1, np
        if (j.eq.1) ggrp(ip,j) = ggrp(ip,j) - xp(ip)*gyl(ip,ilmx,j)
        if (j.eq.2) ggrp(ip,j) = ggrp(ip,j) - yp(ip)*gyl(ip,ilmx,j)
        if (j.eq.3) ggrp(ip,j) = ggrp(ip,j) - zp(ip)*gyl(ip,ilmx,j)
   22   continue
        call xyl(nlm,np,wp,ggrp(1,j),yl,frl)
C        call prmr('component of nabla',frl,nlm)
   20 continue
      call dpadd(ggrp,ggrp(1,2),1,np,1d0)
      call dpadd(ggrp,ggrp(1,3),1,np,1d0)
      call xyl(nlm,np,wp,ggrp,yl,frl)
      call prmr('laplacian',frl,nlm)
      call cexit(1,1)
      end
      subroutine xyl(nlm,np,wp,fp,yl,fl)
C- Yl-projection of function tabulated on an angular mesh
      implicit none
      integer nlm,np,ip,ilm
      double precision fl(nlm),fp(np),yl(np,nlm),wp(np)

      call dpzero(fl,nlm)
      do  20  ip = 1, np
      do  20  ilm = 1, nlm
   20 fl(ilm) = fl(ilm) + fp(ip)*wp(ip)*yl(ip,ilm)
      end
      subroutine prmr(strn,f,nl)
      implicit none
      integer nl,j,fopna,ifi
      double precision f(nl)
      character*(10) fmt, strn*(*)
      ifi = 19
      open(ifi,file='out')
      write(ifi,*) nl, 2
      do  10  j = 1, nl
   10 write(ifi,333) j, f(j)
  333 format(i4, f15.10)
      close(ifi)
      print *, strn
      pause
      end

C Test program to check ropylg
      subroutine fmain
      implicit none
      integer nrx,lmx,nlmx,nlm2,nr,lmax,nlm1,ir,ii,i,l,ilm,i1,i2,nsize,
     .  nnn
      parameter (nrx=20,lmx=6,nlmx=(lmx+1)**2,nlm2=(lmx+2)**2,
     .  nsize=100000, nnn=300)
      double precision cy(16**2),x(nrx),y(nrx),z(nrx),r2(nrx),tops,
     .  ylv(nrx,nlm2),gylv(nrx,nlmx,3),yl(nlm2),gyl(nlm2,3),dr(3),
     .  p(3,nnn),wp(nnn)
      integer oxp,oyp,ozp,or2,oyl,ogyl,orp,ogrp,oggrp,ofrl,ll,lp,np,
     .  nlmf,nph,nth,w(nsize)
      common /w/ w
      common /static/ cy

      call wkinit(nsize)
      call sylmnc(cy,15)

C --- Laplacian of Yl ---
      print *, 'ilm:'
      ilm = 8
      read(*,*) ilm
      lmax = ll(ilm)+1
      nth=lmax+2
      nph=2*nth
      nlmf = (lmax+2)**2
      call fpiint(nth,nph,np,p,wp)
      print *, np, ' angular points'

      call defrr(oxp,     np)
      call defrr(oyp,     np)
      call defrr(ozp,     np)
      call defrr(or2,     np)
      call defrr(oyl,     (lmax+3)**2*np)
      call defrr(ofrl,    nlmf)
      call dcopy(np,p(1,1),3,w(oxp),1)
      call dcopy(np,p(2,1),3,w(oyp),1)
      call dcopy(np,p(3,1),3,w(ozp),1)
      call defrr(ogyl,    nlmf*np*3)
      call defrr(orp,     np)
      call defrr(ogrp,    np*3)
      call defrr(oggrp,   np*3)
      call tl(np,ilm,nlmf,w(orp),w(ogrp),w(oggrp),w(ofrl),
     .  w(oyl),w(ogyl),wp,w(oxp),w(oyp),w(ozp),w(or2))

C --- Compare ropylg against ylg ---
      nr = 5
   99 print *, 'lmax='
      read(*,*) lmax
      if (lmax .gt. lmx) stop 'increase lmx in main'

      call makr(0d0,nr,x,y,z)

C ... Make grad ylm's
      nlm1 = (lmax+2)**2
      call ropyln(nr,x,y,z,lmax+1,nrx,ylv,r2)
      call ropylg(0,lmax,nlmx,nrx,nr,x,y,z,r2,ylv,gylv)

C ... Check against ylg
      tops = 0d0
      do  10  ir = 1, nr
        dr(1) = x(ir)
        dr(2) = y(ir)
        dr(3) = z(ir)
        call ylg(dr,lmax,nlm2,cy,yl,gyl)

        do  11  i = 1, 3

        do  12  l = 0, lmax
        i1 = l*l+1
        i2 = (l+1)**2
   12   print 333, (gyl(ii,i),ii=i1,i2)
  333   format(9f8.5)
        print *
        do  14  l = 0, lmax
        i1 = l*l+1
        i2 = (l+1)**2
        do  16  ii = i1, i2
   16   tops = max(tops,dabs(gyl(ii,i)-gylv(ir,ii,i)))
   14   print 333, (gylv(ir,ii,i)-gyl(ii,i),ii=i1,i2)
        print *, '----------- end of ir,i=', ir,i
   11   continue
   10 continue

      print 335, tops
  335 format(' max errors for grad h:',f12.6)

      end
      subroutine makr(rsm,nr,x,y,z)
      implicit none
      integer nr,i,ir
      double precision rs,rsm,x(1),y(1),z(1)
      real ran1
      rs = rsm
      if (rsm .lt. 1d-9) rs = .5d0
      call ran1in(1)
      do  10  i = 1, nr
        ir = i+1
        x(i) = abs((ran1()-.5d0)*5*rs)
        y(i) = (ran1()-.5d0)*5*rs
        z(i) = (ran1()-.5d0)*5*rs
   10 continue

      x(1) = .3d0*dsqrt(2d0)
      y(1) = .4d0*dsqrt(2d0)
      z(1) = .5d0*dsqrt(2d0)
      end
#endif

