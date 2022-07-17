module m_ropyln
  public ropyln,ropylg
  private
  contains
subroutine ropyln(n,x,y,z,lmax,nd,yl,rsq)
  !- Normalized spheric harmonic polynomials (vectorizes).
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   n     :number of points for which to calculate yl
  !i   x     :x component of Cartesian coordinate
  !i   y     :y component of Cartesian coordinate
  !i   z     :z component of Cartesian coordinate
  !i   lmax  :compute yl for l=0..lmax
  !i   nd    :leading dimension of yl; nd must be >= n
  !o Outputs
  !o   yl    :Ylm(i,ilm) are the (real) spherical harmonics
  !o         :for l=0..lmax and point i.
  !o   rsq   :rsq(i) square of length of point i
  !r Remarks
  !r   yl = real harmonics (see Takao's GW note) * r^l
  !u Updates
  !u  25 Jun 03 (Kino) initialize cx to zero
  ! ----------------------------------------------------------------------
  implicit none
  integer:: nd , n , i , m , lmax , l , kk=-999
  real(8) ,allocatable :: cm_rv(:)
  real(8) ,allocatable :: sm_rv(:)
  real(8) ,allocatable :: q_rv(:)
  real(8) ,allocatable :: h_rv(:)
  double precision :: x(*),y(*),z(*),yl(nd,*),rsq(*),cx(3)
  double precision :: fpi,f2m
  !     call tcn('ropyln')
  if (n > nd) call rx('ropyln: nd too small')
  fpi = 16*datan(1d0)
  allocate(cm_rv(n),sm_rv(n),q_rv(n*2),h_rv(n))
  do  2  i = 1, n
     rsq(i) = x(i)*x(i)+y(i)*y(i)+z(i)*z(i)
2 enddo
  cx = 0d0
  ! --- Loop over m: cx are normalizations for l, l-1, l-2 ---
  f2m = 1d0
  do  10  m = 0, lmax
     call ropcsm ( m , n , x , y , h_rv , cm_rv , sm_rv )
     if (m == 0) then
        cx(1) = dsqrt(1/fpi)
     else
        f2m = f2m*2*m*(2*m-1)
        cx(1) = dsqrt((2*m+1)*2/fpi/f2m)
     endif
     do  11  l = m, lmax
        call ropqln ( m , l , n , rsq , z , cx , q_rv , kk )
        call ropynx ( m , l , kk , n , q_rv , cm_rv , sm_rv, nd , yl )
        cx(3) = cx(2)
        cx(2) = cx(1)
        cx(1) = cx(1)*dsqrt(dble((l+1-m)*(2*l+3))/dble((l+1+m)*(2*l+1)))
11   enddo
10 enddo
  if (allocated(h_rv)) deallocate(h_rv)
  if (allocated(q_rv)) deallocate(q_rv)
  if (allocated(sm_rv)) deallocate(sm_rv)
  if (allocated(cm_rv)) deallocate(cm_rv)
  !     call tcx('ropyln')
end subroutine ropyln

subroutine ropqln(m,l,n,r2,z,cx,q,kk)
  !- Makes qml for m,l. Must be called in sequence l=m,m+1... for fixed m
  !  Returns kk, which points to the current component of q.
  !  These subroutines are utility routines called by ropyln.f.
  !  Kept separate from ropyln because some optimizing compilers have bugs
  !  (e.g. intel ifort version 11).
  !  These routines are the time-critical steps.
  !     implicit none
  integer :: mm,n,i,l,m,kk,k2,k1
  double precision :: q(n,2),r2(n),z(n),cx(3)
  double precision :: a,b,xx,yy
  ! --- Case l=m ---
  if (l == m) then
     a = 1d0
     do  1  mm = 0, m-1
        a = a*(2*mm+1)
1    enddo
     kk = 1
     a = a*cx(1)
     do  2  i = 1, n
        q(i,kk) = a
2    enddo
     return
  endif
  ! --- Case l=m+1 ---
  if (l == m+1) then
     b = 1d0
     do  3  mm = 0, m
        b = b*(2*mm+1)
3    enddo
     b = b*cx(1)
     kk = 2
     do  4  i = 1, n
        q(i,kk) = b*z(i)
4    enddo
     return
  endif
  ! --- Case l=m+2 and higher by recursion ---
  if (l >= m+2) then
     k2 = kk
     k1 = kk+1
     if (k1 == 3) k1 = 1
     xx = -(l+m-1d0)/(l-m)*cx(1)/cx(3)
     yy = (2*l-1d0)/(l-m)*cx(1)/cx(2)
     do  6  i = 1, n
        q(i,k1) = xx*r2(i)*q(i,k1)+yy*z(i)*q(i,k2)
6    enddo
     kk = k1
     return
  endif
end subroutine ropqln

subroutine ropynx(m,l,kk,n,q,cm,sm,nd,yl)
  !     implicit none
  integer :: lav,n,nd,l,i,m,kk
  double precision :: q(n,2),cm(n),sm(n),yl(nd,1)
  lav = l*(l+1)+1
  do  1  i = 1, n
     yl(i,lav+m) = cm(i)*q(i,kk)
1 enddo
  if (m == 0) return
  do  2  i = 1, n
     yl(i,lav-m) = sm(i)*q(i,kk)
2 enddo
end subroutine ropynx

subroutine ropcsm(m,n,x,y,w,cm,sm)
  !- Makes cm and sm. Must be called in sequence m=0,1,2...
  implicit none
  integer :: m,n,i
  double precision :: x(n),y(n),w(n),cm(n),sm(n)
  ! --- Case m=0 ---
  if (m == 0) then
     do  1  i = 1, n
        cm(i) = 1d0
1    enddo
     do  2  i = 1, n
        sm(i) = 0d0
2    enddo
     return
  endif
  ! --- Case m=1 ---
  if (m == 1) then
     do  3  i = 1, n
        cm(i) = x(i)
3    enddo
     do  4  i = 1, n
        sm(i) = y(i)
4    enddo
     return
  endif
  ! --- Case m ge 2 ---
  if (m >= 2) then
     do  5  i = 1, n
        w(i) = cm(i)
5    enddo
     do  6  i = 1, n
        cm(i) = x(i)*cm(i) - y(i)*sm(i)
6    enddo
     do  7  i = 1, n
        sm(i) = y(i)*w(i) + x(i)*sm(i)
7    enddo
     return
  endif
end subroutine ropcsm

subroutine ropylg(lp,lmax,ndim,nrx,nr,x,y,z,r2,yl,gyl)
  !- Gradients of YL's (polynomials) for a set of points, with YL as input
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   lp    :if nonzero, adds term  r^l grad (r^-l Yl).
  !i   lmax  :maximum l for a given site
  !i   ndim  :dimensions gyl.  Must be at least (lmax+1)**2
  !i   nrx   :leading dimension of yl,gyl
  !i   nr    :number of points
  !i   x,y,z :cartesian coordinates of points
  !i   r2    :x^2+y^2+z^2
  !i   yl    :Spherical harmonic polynomials YL.  YL's must be normalized
  !i         :and must be made through lmax+1 (i.e. nlm=1..(lmax+2)**2)
  !o Outputs
  !o   gyl   :gradient of yl
  !l Local variables
  !l         :
  !r Remarks
  !r
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: lp,lmax,ndim,nrx,nr
  double precision :: x(nr),y(nr),z(nr),yl(nrx,*),gyl(nrx,ndim,3),r2(*)
  ! ... Local parameters
  integer :: ilm,kx1,kx2,ky1,ky2,kz,l,m,i
  double precision :: cx1,cx2,cy1,cy2,cz,f
  if ((lmax+1)**2 > ndim) call rx('ropylg: ndim too small')
  ! --- Gradients of yl's ---
  ilm = 0
  do    l = 0, lmax
     do    m = -l, l
        ilm = ilm+1
        call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
        do    i = 1, nr
           f = (2*l+1)/r2(i)
           gyl(i,ilm,1) = (yl(i,ilm)*x(i) - cx1*yl(i,kx1) - cx2*yl(i,kx2))*f
           gyl(i,ilm,2) = (yl(i,ilm)*y(i) - cy1*yl(i,ky1) - cy2*yl(i,ky2))*f
           gyl(i,ilm,3) = (yl(i,ilm)*z(i) - cz*yl(i,kz))*f
        enddo
     enddo
  enddo
  if (lp == 0) return
  ! --- Add r**l (grad r**-l) yl ---
  ilm = 0
  do    l = 0, lmax
     do    m = -l, l
        ilm = ilm+1
        do    i = 1, nr
           gyl(i,ilm,1) = gyl(i,ilm,1) - l*x(i)/r2(i)*yl(i,ilm)
           gyl(i,ilm,2) = gyl(i,ilm,2) - l*y(i)/r2(i)*yl(i,ilm)
           gyl(i,ilm,3) = gyl(i,ilm,3) - l*z(i)/r2(i)*yl(i,ilm)
        enddo
     enddo
  enddo
end subroutine ropylg
end module m_ropyln

!$$$#if TEST
!$$$C Test program to check ropyln
!$$$      subroutine fmain
!$$$      implicit none
!$$$      integer nrx,nlmx,nr,lmax,nlm1,ir,ii,l,ilm,i1,i2,nsize
!$$$      parameter (nrx=20,nlmx=49,nsize=100000)
!$$$      double precision cy(16**2),x(nrx),y(nrx),z(nrx),r2(nrx),
!$$$     .ylv(nrx,nlmx),ylok(nrx,nlmx),dr(3),tops,ylm(nlmx)
!$$$
!$$$      call wkinit(nsize)
!$$$      call wkfast(.true.)
!$$$      call sylmnc(cy,15)
!$$$
!$$$      lmax = 2
!$$$   99 print *, 'lmax='
!$$$      read(*,*) lmax
!$$$
!$$$      call makr(0d0,nr,x,y,z)
!$$$
!$$$C ... Make nonvectorized ylm's
!$$$      nlm1 = (lmax+1)**2
!$$$      do  ir = 1, nr
!$$$        dr(1) = x(ir)
!$$$        dr(2) = y(ir)
!$$$        dr(3) = z(ir)
!$$$        call sylm(dr,ylm,lmax,r2)
!$$$        do  ilm = 1, nlm1
!$$$          ylok(ir,ilm) = cy(ilm)*ylm(ilm)
!$$$        enddo
!$$$      enddo
!$$$C     Test: Y_1-1 = sqrt(3/4/pi) y
!$$$C     print *, y(1) * dsqrt(0.75/4/atan(1d0))
!$$$C     print *, ylok(1,2)
!$$$
!$$$      call ropyln(nr,x,y,z,lmax,nrx,ylv,r2)
!$$$
!$$$      tops = 0
!$$$      do  10  ir = 1, nr
!$$$        do  12  l = 0, lmax
!$$$          i1 = l*l+1
!$$$          i2 = (l+1)**2
!$$$          print 333, (ylok(ir,ii),ii=i1,i2)
!$$$   12   continue
!$$$  333   format(9f8.5)
!$$$        print *
!$$$        do  14  l = 0, lmax
!$$$          i1 = l*l+1
!$$$          i2 = (l+1)**2
!$$$          do  16  ii = i1, i2
!$$$            tops = max(tops,dabs(ylok(ir,ii)-ylv(ir,ii)))
!$$$   16     continue
!$$$          print 333, (ylok(ir,ii)-ylv(ir,ii),ii=i1,i2)
!$$$   14   continue
!$$$        print *, '----------------'
!$$$   10 continue
!$$$
!$$$      print 335, tops
!$$$  335 format(' max error for ylm:',1pe12.2)
!$$$      end
!$$$      subroutine makr(rsm,nr,x,y,z)
!$$$      implicit none
!$$$      integer nr,i,ir
!$$$      double precision rs,rsm,x(1),y(1),z(1)
!$$$      real ran1
!$$$      rs = rsm
!$$$      if (rsm .lt. 1d-9) rs = .5d0
!$$$      call ran1in(1)
!$$$      nr = 5
!$$$      do  10  i = 1, nr
!$$$        ir = i+1
!$$$        x(i) = abs((ran1()-.5d0)*5*rs)
!$$$        y(i) = (ran1()-.5d0)*5*rs
!$$$        z(i) = (ran1()-.5d0)*5*rs
!$$$   10 continue
!$$$
!$$$      end
!$$$#endif


!$$$#if TEST !======================================================================
!$$$      subroutine tl(np,ilmx,nlm,rp,grp,ggrp,frl,yl,gyl,wp,xp,yp,zp,r2)
!$$$      implicit none
!$$$      integer np,ilmx,nlm
!$$$      double precision yl(np,nlm),gyl(np,nlm,3),xp(np),yp(np),zp(np),
!$$$     .  wp(np),rp(np),grp(np,3),ggrp(np,3),frl(nlm),r2(np)
!$$$      integer ll,ip,lmax,j,ilm,l
!$$$      lmax = ll(nlm)
!$$$      call ropyln(np,xp,yp,zp,lmax+1,np,yl,r2)
!$$$      call xyl(nlm,np,wp,yl(1,ilmx),yl,frl)
!$$$      call prmr('yl(ilmx)',frl,nlm)
!$$$C ... Show laplacian yl is 0:
!$$$      call ropylg(0,lmax,nlm,np,np,xp,yp,zp,r2,yl,gyl)
!$$$      call dpzero(ggrp,np*3)
!$$$      do  120  j = 1, 3
!$$$        call xyl(nlm,np,wp,gyl(1,ilmx,j),yl,frl)
!$$$        call prmr('component of grad',frl,nlm)
!$$$        do  118  ilm = 1, nlm
!$$$        do  118  ip = 1, np
!$$$  118   ggrp(ip,j) = ggrp(ip,j) + frl(ilm)*gyl(ip,ilm,j)
!$$$        call xyl(nlm,np,wp,ggrp(1,j),yl,frl)
!$$$        call prmr('component of nabla',frl,nlm)
!$$$  120 continue
!$$$      call dpadd(ggrp,ggrp(1,2),1,np,1d0)
!$$$      call dpadd(ggrp,ggrp(1,3),1,np,1d0)
!$$$      call xyl(nlm,np,wp,ggrp,yl,frl)
!$$$      call prmr('laplacian',frl,nlm)
!$$$C     call cexit(1,1)
!$$$C ... Show laplacian r^-l yl is -l(l+1) yl;
!$$$C     Use grad(r^-1 yl) = 1/r sum_L (a_L yl)
!$$$      call ropylg(1,lmax,nlm,np,np,xp,yp,zp,r2,yl,gyl)
!$$$C ... Term 1/r grad r (grad r^-l yl)
!$$$      call dpzero(ggrp,np*3)
!$$$      do  20  j = 1, 3
!$$$        call xyl(nlm,np,wp,gyl(1,ilmx,j),yl,frl)
!$$$        call prmr('component of grad',frl,nlm)
!$$$        do  18  ilm = 1, nlm
!$$$        do  18  ip = 1, np
!$$$   18   ggrp(ip,j) = ggrp(ip,j) + frl(ilm)*gyl(ip,ilm,j)
!$$$        call xyl(nlm,np,wp,ggrp(1,j),yl,frl)
!$$$C        call prmr('1st term of nabla',frl,nlm)
!$$$C ... Term grad (1/r) . (grad r^-l yl)
!$$$        do  22  ip = 1, np
!$$$        if (j.eq.1) ggrp(ip,j) = ggrp(ip,j) - xp(ip)*gyl(ip,ilmx,j)
!$$$        if (j.eq.2) ggrp(ip,j) = ggrp(ip,j) - yp(ip)*gyl(ip,ilmx,j)
!$$$        if (j.eq.3) ggrp(ip,j) = ggrp(ip,j) - zp(ip)*gyl(ip,ilmx,j)
!$$$   22   continue
!$$$        call xyl(nlm,np,wp,ggrp(1,j),yl,frl)
!$$$C        call prmr('component of nabla',frl,nlm)
!$$$   20 continue
!$$$      call dpadd(ggrp,ggrp(1,2),1,np,1d0)
!$$$      call dpadd(ggrp,ggrp(1,3),1,np,1d0)
!$$$      call xyl(nlm,np,wp,ggrp,yl,frl)
!$$$      call prmr('laplacian',frl,nlm)
!$$$      call cexit(1,1)
!$$$      end
!$$$      subroutine xyl(nlm,np,wp,fp,yl,fl)
!$$$C- Yl-projection of function tabulated on an angular mesh
!$$$      implicit none
!$$$      integer nlm,np,ip,ilm
!$$$      double precision fl(nlm),fp(np),yl(np,nlm),wp(np)
!$$$      call dpzero(fl,nlm)
!$$$      do  20  ip = 1, np
!$$$      do  20  ilm = 1, nlm
!$$$   20 fl(ilm) = fl(ilm) + fp(ip)*wp(ip)*yl(ip,ilm)
!$$$      end
!$$$      subroutine prmr(strn,f,nl)
!$$$      implicit none
!$$$      integer nl,j,fopna,ifi
!$$$      double precision f(nl)
!$$$      character*(10) fmt, strn*(*)
!$$$      ifi = 19
!$$$      open(ifi,file='out')
!$$$      write(ifi,*) nl, 2
!$$$      do  10  j = 1, nl
!$$$   10 write(ifi,333) j, f(j)
!$$$  333 format(i4, f15.10)
!$$$      close(ifi)
!$$$      print *, strn
!$$$      pause
!$$$      end


!$$$C Test program to check ropylg !===================================================
!$$$      subroutine fmain
!$$$      implicit none
!$$$      integer nrx,lmx,nlmx,nlm2,nr,lmax,nlm1,ir,ii,i,l,ilm,i1,i2,nsize,
!$$$     .  nnn
!$$$      parameter (nrx=20,lmx=6,nlmx=(lmx+1)**2,nlm2=(lmx+2)**2,
!$$$     .  nsize=100000, nnn=300)
!$$$      double precision cy(16**2),x(nrx),y(nrx),z(nrx),r2(nrx),tops,
!$$$     .  ylv(nrx,nlm2),gylv(nrx,nlmx,3),yl(nlm2),gyl(nlm2,3),dr(3),
!$$$     .  p(3,nnn),wp(nnn)
!$$$      integer oxp,oyp,ozp,or2,oyl,ogyl,orp,ogrp,oggrp,ofrl,ll,lp,np,
!$$$     .  nlmf,nph,nth,w(nsize)
!$$$      common /w/ w
!$$$      common /static/ cy
!$$$      call wkinit(nsize)
!$$$      call sylmnc(cy,15)
!$$$C --- Laplacian of Yl ---
!$$$      print *, 'ilm:'
!$$$      ilm = 8
!$$$      read(*,*) ilm
!$$$      lmax = ll(ilm)+1
!$$$      nth=lmax+2
!$$$      nph=2*nth
!$$$      nlmf = (lmax+2)**2
!$$$      call fpiint(nth,nph,np,p,wp)
!$$$      print *, np, ' angular points'
!$$$      call defrr(oxp,     np)
!$$$      call defrr(oyp,     np)
!$$$      call defrr(ozp,     np)
!$$$      call defrr(or2,     np)
!$$$      call defrr(oyl,     (lmax+3)**2*np)
!$$$      call defrr(ofrl,    nlmf)
!$$$      call dcopy(np,p(1,1),3,w(oxp),1)
!$$$      call dcopy(np,p(2,1),3,w(oyp),1)
!$$$      call dcopy(np,p(3,1),3,w(ozp),1)
!$$$      call defrr(ogyl,    nlmf*np*3)
!$$$      call defrr(orp,     np)
!$$$      call defrr(ogrp,    np*3)
!$$$      call defrr(oggrp,   np*3)
!$$$      call tl(np,ilm,nlmf,w(orp),w(ogrp),w(oggrp),w(ofrl),
!$$$     .  w(oyl),w(ogyl),wp,w(oxp),w(oyp),w(ozp),w(or2))
!$$$C --- Compare ropylg against ylg ---
!$$$      nr = 5
!$$$   99 print *, 'lmax='
!$$$      read(*,*) lmax
!$$$      if (lmax .gt. lmx) stop 'increase lmx in main'
!$$$      call makr(0d0,nr,x,y,z)
!$$$C ... Make grad ylm's
!$$$      nlm1 = (lmax+2)**2
!$$$      call ropyln(nr,x,y,z,lmax+1,nrx,ylv,r2)
!$$$      call ropylg(0,lmax,nlmx,nrx,nr,x,y,z,r2,ylv,gylv)
!$$$C ... Check against ylg
!$$$      tops = 0d0
!$$$      do  10  ir = 1, nr
!$$$        dr(1) = x(ir)
!$$$        dr(2) = y(ir)
!$$$        dr(3) = z(ir)
!$$$        call ylg(dr,lmax,nlm2,cy,yl,gyl)
!$$$        do  11  i = 1, 3
!$$$        do  12  l = 0, lmax
!$$$        i1 = l*l+1
!$$$        i2 = (l+1)**2
!$$$   12   print 333, (gyl(ii,i),ii=i1,i2)
!$$$  333   format(9f8.5)
!$$$        print *
!$$$        do  14  l = 0, lmax
!$$$        i1 = l*l+1
!$$$        i2 = (l+1)**2
!$$$        do  16  ii = i1, i2
!$$$   16   tops = max(tops,dabs(gyl(ii,i)-gylv(ir,ii,i)))
!$$$   14   print 333, (gylv(ir,ii,i)-gyl(ii,i),ii=i1,i2)
!$$$        print *, '----------- end of ir,i=', ir,i
!$$$   11   continue
!$$$   10 continue
!$$$      print 335, tops
!$$$  335 format(' max errors for grad h:',f12.6)
!$$$      end
!$$$      subroutine makr(rsm,nr,x,y,z)
!$$$      implicit none
!$$$      integer nr,i,ir
!$$$      double precision rs,rsm,x(1),y(1),z(1)
!$$$      real ran1
!$$$      rs = rsm
!$$$      if (rsm .lt. 1d-9) rs = .5d0
!$$$      call ran1in(1)
!$$$      do  10  i = 1, nr
!$$$        ir = i+1
!$$$        x(i) = abs((ran1()-.5d0)*5*rs)
!$$$        y(i) = (ran1()-.5d0)*5*rs
!$$$        z(i) = (ran1()-.5d0)*5*rs
!$$$   10 continue
!$$$      x(1) = .3d0*dsqrt(2d0)
!$$$      y(1) = .4d0*dsqrt(2d0)
!$$$      z(1) = .5d0*dsqrt(2d0)
!$$$      end
!$$$#endif

