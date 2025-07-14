module m_ropyln ! Normalized spheric harmonic polynomials (vectorizes).
  public ropyln,ropylg
  private
  integer,parameter::lmax=20
contains
  real(8) function f2m(m)
    integer:: mm,m
    real(8),save:: fff(lmax)
    logical,save:: init=.true.
    if(lmax<m) call rx('f2m')
    if(init) then
      fff(1:lmax)=[(dble(2*mm*(2*mm-1)),mm=1,lmax)]
      init=.false.
    endif
    f2m = product(fff(1:m))
  end function f2m
  real(8) function pfff(m)
    integer:: mm,m
    real(8),save:: fff(0:lmax)
    logical,save::init=.true.
    if(lmax<m) call rx('pfff')
    if(init) then
      fff(0:lmax)=[(dble(2*mm+1),mm=0,lmax)]
      init=.false.
    endif
    pfff=product(fff(0:m))
  end function pfff

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
    integer:: nd , n , i , m , lmax , l , kk=-999,lav,mm,kp
    real(8):: yl(nd,*), x(n),y(n),z(n),rsq(n),cx(3),a,b,xx,yy,fff(0:lmax+1),f2mx,cx1(0:lmax+1,0:lmax),cx2,cx3,cx0
    real(8),parameter:: fpi = 16*datan(1d0)
    real(8) :: cm_rv(n,0:lmax),sm_rv(n,0:lmax), q_rv(n,0:lmax)
    mloop0: do m = 0, lmax
      if(m==0) then
        cm_rv(:,0)=1d0
        sm_rv(:,0)=0d0
      elseif(m==1) then  
        cm_rv(:,1)=x
        sm_rv(:,1)=y
      else
        cm_rv(:,m) = x*cm_rv(:,m-1) -  y*sm_rv(:,m-1)
        sm_rv(:,m) = y*cm_rv(:,m-1) +  x*sm_rv(:,m-1)
      endif  
    enddo mloop0
    
    rsq = x**2+y**2+z**2
    cx = 0d0 
    cx0 = dsqrt(1/fpi)
    mloop: do  10  m = 0, lmax     ! --- Loop over m: cx are normalizations for l, l-1, l-2 ---
       if (m >0) cx0 = dsqrt((2*m+1)*2/fpi/f2m(m))
       cx1(m,m)=cx0 !for l=m
       do l = m, lmax   !  These routines are the time-critical steps.
         cx1(l+1,m) = cx1(l,m)*dsqrt(dble((l+1-m)*(2*l+3))/dble((l+1+m)*(2*l+1)))
       enddo
       lloop0: do l = m, lmax   !  These routines are the time-critical steps.
         kk = l-m !mod(l-m,2)+1    !  Returns kk, which points to the current component of q_rv.     !         fff(0:m-1)=[(dble(2*mm+1),mm=0,m-1)]
          if (kk==0) then
             q_rv(:,kk) = pfff(m-1)*cx1(l,m) 
          elseif(kk==1) then
             q_rv(:,kk) = pfff(m)*cx1(l,m)*z(:)
          else 
            q_rv(:,kk) = -(l+m-1d0)/(l-m)*cx1(l,m)/cx1(l-2,m)*rsq(:) *q_rv(:,kk-2) &
                 + (2*l-1d0)/(l-m)*cx1(l,m)/cx1(l-1,m)*z(:) *q_rv(:,kk-1)
          endif
       enddo lloop0
       lloop: do l = m, lmax   
         yl         (1:n, l*(l+1)+1+m) = cm_rv(:,m)*q_rv(:,l-m)
         if(m/=0) yl(1:n, l*(l+1)+1-m) = sm_rv(:,m)*q_rv(:,l-m)
       enddo lloop
10  enddo mloop
  end subroutine ropyln
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
