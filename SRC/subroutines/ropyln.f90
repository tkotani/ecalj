module m_ropyln ! Normalized spheric harmonic polynomials (vectorizes).
  public ropyln,ropylg, ropyln_d
  private
contains
  subroutine ropyln(n,x,y,z,lmax,nd, yl,rsq) ! Normalized spheric harmonic polynomials (vectorizes).
    !i   n     :number of points for which to calculate yl
    !i   x,y,z:f Cartesian coordinate
    !i   lmax  :compute yl for l=0..lmax
    !i   nd    :leading dimension of yl; nd must be >= n
    !o Outputs
    !o   yl  :Ylm(i,ilm) are the real spherical harmonics for l=0..lmax and point i.
    !         yl = real harmonics (see Takao's GW note) * r^l
    !o   rsq   :rsq(i) square of length of point i
    implicit none
    integer:: nd , n , i , m , lmax , l , kk=-999,lav,k1,mm,k2
    real(8):: yl(nd,*), x(n),y(n),z(n),rsq(n),cx(3),f2m,a,b,xx,yy,fff(0:lmax+1)
    real(8),parameter:: fpi = 16*datan(1d0)
    real(8),allocatable :: cm_rv(:),sm_rv(:), q_rv(:,:), h_rv(:)
    if (n > nd) call rx('ropyln: nd too small')
    allocate(cm_rv(n),sm_rv(n),q_rv(n,2),h_rv(n))
    rsq = x**2+y**2+z**2
    cx = 0d0 
    mloop: do  10  m = 0, lmax     ! --- Loop over m: cx are normalizations for l, l-1, l-2 ---
       if (m == 0) then
          cm_rv=1d0
          sm_rv=0d0
       elseif (m == 1) then
          cm_rv=x
          sm_rv=y
       elseif (m >= 2) then
          h_rv=cm_rv
          cm_rv = x*cm_rv - y*sm_rv
          sm_rv = y*h_rv +  x*sm_rv
       endif
       if (m == 0) then !call ropcsm ( m , n , x , y , h_rv , cm_rv , sm_rv )
          cx(1) = dsqrt(1/fpi)
       else
          fff(1:m)=[(dble(2*mm*(2*mm-1)),mm=1,m)]
          f2m = product(fff(1:m))
          cx(1) = dsqrt((2*m+1)*2/fpi/f2m)
       endif
       lloop: do  11  l = m, lmax   !  These routines are the time-critical steps.
          if (l == m) then
             kk = 1    !  Returns kk, which points to the current component of q_rv.
             fff(0:m)=[(dble(2*mm+1),mm=0,m-1),cx(1)]
             q_rv(:,kk) = product(fff(0:m))
          elseif(l == m+1) then
             kk = 2
             fff(0:m+1)=[(dble(2*mm+1),mm=0,m),cx(1)]
             q_rv(:,kk) = product(fff(0:m+1))*z(:)
          elseif (l >= m+2) then
             k1 = kk+1
             if (k1 == 3) k1 = 1
             k2 = kk
             q_rv(:,k1) = -(l+m-1d0)/(l-m)*cx(1)/cx(3)*rsq(:)*q_rv(:,k1) + (2*l-1d0)/(l-m)*cx(1)/cx(2)*z(:)*q_rv(:,k2)
             kk = k1
          endif
          lav = l*(l+1)+1
          yl(1:n,lav+m)          = cm_rv(:)*q_rv(:,kk)
          if(m/=0) yl(1:n,lav-m) = sm_rv(:)*q_rv(:,kk)
          cx(1:3) = [cx(1)*dsqrt(dble((l+1-m)*(2*l+3))/dble((l+1+m)*(2*l+1))), cx(1), cx(2)]
11     enddo lloop
10  enddo mloop
  end subroutine ropyln
  subroutine ropyln_d(n,x,y,z,lmax,nd, yl,rsq) ! Normalized spheric harmonic polynomials (vectorizes).
    !i   n     :number of points for which to calculate yl
    !i   x,y,z:f Cartesian coordinate
    !i   lmax  :compute yl for l=0..lmax
    !i   nd    :leading dimension of yl; nd must be >= n
    !o Outputs
    !o   yl  :Ylm(i,ilm) are the real spherical harmonics for l=0..lmax and point i.
    !         yl = real harmonics (see Takao's GW note) * r^l
    !o   rsq   :rsq(i) square of length of point i
    !$acc routine seq
    implicit none
    integer:: nd , n , i , m , lmax , l , kk=-999,lav,k1,mm,k2
    real(8):: yl(nd,*), x(n),y(n),z(n),rsq(n),cx(3),f2m,a,b,xx,yy,fff(0:lmax+1)
    real(8),parameter:: fpi = 16*datan(1d0)
    real(8),allocatable :: cm_rv(:),sm_rv(:), q_rv(:,:), h_rv(:)
    ! if (n > nd) call rx('ropyln: nd too small')
    allocate(cm_rv(n),sm_rv(n),q_rv(n,2),h_rv(n))
    rsq = x**2+y**2+z**2
    cx = 0d0 
    mloop: do  10  m = 0, lmax     ! --- Loop over m: cx are normalizations for l, l-1, l-2 ---
       if (m == 0) then
          cm_rv=1d0
          sm_rv=0d0
       elseif (m == 1) then
          cm_rv=x
          sm_rv=y
       elseif (m >= 2) then
          h_rv=cm_rv
          cm_rv = x*cm_rv - y*sm_rv
          sm_rv = y*h_rv +  x*sm_rv
       endif
       if (m == 0) then !call ropcsm ( m , n , x , y , h_rv , cm_rv , sm_rv )
          cx(1) = dsqrt(1/fpi)
       else
          fff(1:m)=[(dble(2*mm*(2*mm-1)),mm=1,m)]
          f2m = product(fff(1:m))
          cx(1) = dsqrt((2*m+1)*2/fpi/f2m)
       endif
       lloop: do  11  l = m, lmax   !  These routines are the time-critical steps.
          if (l == m) then
             kk = 1    !  Returns kk, which points to the current component of q_rv.
             fff(0:m)=[(dble(2*mm+1),mm=0,m-1),cx(1)]
             q_rv(:,kk) = product(fff(0:m))
          elseif(l == m+1) then
             kk = 2
             fff(0:m+1)=[(dble(2*mm+1),mm=0,m),cx(1)]
             q_rv(:,kk) = product(fff(0:m+1))*z(:)
          elseif (l >= m+2) then
             k1 = kk+1
             if (k1 == 3) k1 = 1
             k2 = kk
             q_rv(:,k1) = -(l+m-1d0)/(l-m)*cx(1)/cx(3)*rsq(:)*q_rv(:,k1) + (2*l-1d0)/(l-m)*cx(1)/cx(2)*z(:)*q_rv(:,k2)
             kk = k1
          endif
          lav = l*(l+1)+1
          yl(1:n,lav+m)          = cm_rv(:)*q_rv(:,kk)
          if(m/=0) yl(1:n,lav-m) = sm_rv(:)*q_rv(:,kk)
          cx(1:3) = [cx(1)*dsqrt(dble((l+1-m)*(2*l+3))/dble((l+1+m)*(2*l+1))), cx(1), cx(2)]
11     enddo lloop
10  enddo mloop
  end subroutine ropyln_d
  subroutine ropylg(lp,lmax,ndim,nrx,nr,x,y,z,r2,yl,gyl)!- Gradients of YL's (polynomials) for a set of points, with YL as input
    use m_ll,only:ll
    use m_scg,only:scglp1
    !i   lp    :if nonzero, adds term  r^l grad (r^-l Yl).
    !i   lmax  :maximum l for a given site
    !i   ndim  :dimensions gyl.  Must be at least (lmax+1)**2
    !i   nr    :number of points
    !i   x,y,z :cartesian coordinates of points
    !i   r2    :x^2+y^2+z^2
    !i   yl    :Spherical harmonic polynomials YL.  YL's must be normalized
    !i         :and must be made through lmax+1 (i.e. nlm=1..(lmax+2)**2)
    !o Outputs
    !o   gyl   :gradient of yl
    implicit none
    integer :: lp,lmax,ndim,nrx,nr, ilm,kx1,kx2,ky1,ky2,kz,l,m,i
    real(8) :: x(nr),y(nr),z(nr),yl(nr,*),gyl(nr,ndim,3),r2(*),cx1,cx2,cy1,cy2,cz,f
    if ((lmax+1)**2 > ndim) call rx('ropylg: ndim too small')
    do ilm = 1,(lmax+1)**2
       l=ll(ilm)
       call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
       do i = 1, nr
          gyl(i,ilm,1:3) = (2*l+1)/r2(i)* [& !! --- Gradients of yl's ---
               yl(i,ilm)*x(i) - cx1*yl(i,kx1) - cx2*yl(i,kx2),&
               yl(i,ilm)*y(i) - cy1*yl(i,ky1) - cy2*yl(i,ky2),&
               yl(i,ilm)*z(i)                                 - cz*yl(i,kz)] &
               - merge(1,0,lp/= 0) * l/r2(i)*yl(i,ilm)*[x(i), y(i), z(i)] !Add r**l (grad r**-l) yl if lp/=0
       enddo
    enddo
  end subroutine ropylg
end module m_ropyln
