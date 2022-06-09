module m_ropyln
  public ropyln
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
  !     implicit none
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
  cx(1) = 0
  cx(2) = 0
  cx(3) = 0

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
        call ropynx ( m , l , kk , n , q_rv , cm_rv , sm_rv &
             , nd , yl )
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

end module m_ropyln

!     Separate subroutines below to avoid problems with some
!     optimizing compilers.
!      subroutine ropqln(m,l,n,r2,z,cx,q,kk)
!C- Makes qml for m,l. Must be called in sequence l=m,m+1... for fixed m
!c  Returns kk, which points to the current component of q.
!      implicit none
!      integer mm,n,i,l,m,kk,k2,k1
!      double precision q(n,2),r2(n),z(n),cx(3)
!      double precision a,b,xx,yy

!C --- Case l=m ---
!      if (l .eq. m) then
!        a = 1d0
!        do  1  mm = 0, m-1
!    1   a = a*(2*mm+1)
!        kk = 1
!        a = a*cx(1)
!        do  2  i = 1, n
!    2   q(i,kk) = a
!        return
!      endif

!C --- Case l=m+1 ---
!      if (l .eq. m+1) then
!        b = 1d0
!        do  3  mm = 0, m
!    3   b = b*(2*mm+1)
!        b = b*cx(1)
!        kk = 2
!        do  4  i = 1, n
!    4   q(i,kk) = b*z(i)
!        return
!      endif

!C --- Case l=m+2 and higher by recursion ---
!      if (l .ge. m+2) then
!        k2 = kk
!        k1 = kk+1
!        if (k1 .eq. 3) k1 = 1
!        xx = -(l+m-1d0)/(l-m)*cx(1)/cx(3)
!        yy = (2*l-1d0)/(l-m)*cx(1)/cx(2)
!        do  6  i = 1, n
!    6   q(i,k1) = xx*r2(i)*q(i,k1)+yy*z(i)*q(i,k2)
!        kk = k1
!        return
!      endif
!      end
!      subroutine ropynx(m,l,kk,n,q,cm,sm,nd,yl)
!      implicit none
!      integer lav,n,nd,l,i,m,kk
!      double precision q(n,2),cm(n),sm(n),yl(nd,1)
!      lav = l*(l+1)+1
!      do  1  i = 1, n
!    1 yl(i,lav+m) = cm(i)*q(i,kk)
!      if (m .eq. 0) return
!      do  2  i = 1, n
!    2 yl(i,lav-m) = sm(i)*q(i,kk)
!      end
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

