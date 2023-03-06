module m_ropyln
  public ropyln,ropylg
  private
contains
  subroutine ropyln(n,x,y,z,lmax,nd,yl,rsq) !- Normalized spheric harmonic polynomials (vectorizes).
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
2   enddo
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
11     enddo
10  enddo
    if (allocated(h_rv)) deallocate(h_rv)
    if (allocated(q_rv)) deallocate(q_rv)
    if (allocated(sm_rv)) deallocate(sm_rv)
    if (allocated(cm_rv)) deallocate(cm_rv)
  end subroutine ropyln
  subroutine ropqln(m,l,n,r2,z,cx,q,kk) !- Makes qml for m,l. Must be called in sequence l=m,m+1... for fixed m
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
1      enddo
       kk = 1
       a = a*cx(1)
       do  2  i = 1, n
          q(i,kk) = a
2      enddo
       return
    endif
    ! --- Case l=m+1 ---
    if (l == m+1) then
       b = 1d0
       do  3  mm = 0, m
          b = b*(2*mm+1)
3      enddo
       b = b*cx(1)
       kk = 2
       do  4  i = 1, n
          q(i,kk) = b*z(i)
4      enddo
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
6      enddo
       kk = k1
       return
    endif
  end subroutine ropqln
  subroutine ropynx(m,l,kk,n,q,cm,sm,nd,yl)
         implicit none
    integer :: lav,n,nd,l,i,m,kk
    double precision :: q(n,2),cm(n),sm(n),yl(nd,1)
    lav = l*(l+1)+1
    do  1  i = 1, n
       yl(i,lav+m) = cm(i)*q(i,kk)
1   enddo
    if (m == 0) return
    do  2  i = 1, n
       yl(i,lav-m) = sm(i)*q(i,kk)
2   enddo
  end subroutine ropynx
  subroutine ropcsm(m,n,x,y,w,cm,sm) !- Makes cm and sm. Must be called in sequence m=0,1,2...
    implicit none
    integer :: m,n,i
    double precision :: x(n),y(n),w(n),cm(n),sm(n)
    ! --- Case m=0 ---
    if (m == 0) then
       do  1  i = 1, n
          cm(i) = 1d0
1      enddo
       do  2  i = 1, n
          sm(i) = 0d0
2      enddo
       return
    endif
    ! --- Case m=1 ---
    if (m == 1) then
       do  3  i = 1, n
          cm(i) = x(i)
3      enddo
       do  4  i = 1, n
          sm(i) = y(i)
4      enddo
       return
    endif
    ! --- Case m ge 2 ---
    if (m >= 2) then
       do  5  i = 1, n
          w(i) = cm(i)
5      enddo
       do  6  i = 1, n
          cm(i) = x(i)*cm(i) - y(i)*sm(i)
6      enddo
       do  7  i = 1, n
          sm(i) = y(i)*w(i) + x(i)*sm(i)
7      enddo
       return
    endif
  end subroutine ropcsm
  subroutine ropylg(lp,lmax,ndim,nrx,nr,x,y,z,r2,yl,gyl)!- Gradients of YL's (polynomials) for a set of points, with YL as input
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
    implicit none
    integer :: lp,lmax,ndim,nrx,nr
    double precision :: x(nr),y(nr),z(nr),yl(nrx,*),gyl(nrx,ndim,3),r2(*)
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
