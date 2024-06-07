module m_purewintz !Gaussian integral along im omg axis. pure, we need it in module for interface.
  public wintzsg_npm_wgtim
contains
  pure subroutine wintzsg_npm_wgtim(npm,ua_,expa_,we,esmr, wgtim) !Gaussian integral along im omg axis
    !     loop over w' = (1-x)/x, frequencies in Wc(k,w')
    !     {x} are gaussian-integration points between (0,1)
    use m_genallcf_v3,only: niw
    use m_readfreq_r,only: wt=>wwx,x=>freqx
    implicit none
    integer,intent(in)::            npm
    real(8),intent(in)::                ua_,expa_(niw),we,esmr
    real(8),intent(out)::                              wgtim(0:npm*niw)
    real(8):: we2,weh,wel,weh2,wel2,cons(1:niw),omd(niw),omd2,rup,rdn,sss,sig2,omd2w(niw), wintsf
    integer :: ie,i,ix,verbose
    real(8),parameter :: pi=3.1415926535897932d0, rmax=2d0 !rmax =2 is by ferdi. Is it optimum? See wintz
    real(8)::sig, smxowx, ee, omg, ww,cons1,cons2,xx,aw,eee,aw2
    logical :: test
    if(esmr==0d0) then ! if w = e the integral = -v(0)/2 ! frequency integral
       wgtim(0) = 0d0
       cons = 1d0/(we**2*x**2 + (1d0-x)**2)
       wgtim(1:niw)= we*cons(:)*wt(:)*(-1d0/pi)
       if(npm==2) wgtim(niw+1:2*niw)= -cons(1:niw)* (1d0/x(1:niw)-1d0)*wt(1:niw)*(-1d0/pi)
       if(dabs(we)<rmax/ua_) wgtim(0)= -sum(wgtim(1:niw)*expa_(1:niw)) -0.5d0*dsign(1.d0,we)*dexp(we**2*ua_**2)*erfc(ua_*dabs(we))
    else
       sig  = .5d0*esmr
       sig2 = 2d0*sig**2
       wgtim=0d0
       omd   = 1d0/x - 1d0
       omd2w = omd**2 + we**2
       where(omd2w/sig2  > 5d-3) cons = (1d0 - exp (- omd2w/sig2))/omd2w
       where(omd2w/sig2 <= 5d-3) cons = (1d0/sig2 - omd2w/sig2**2/2d0 + omd2w**2/sig2**3/6d0 - omd2w**3/sig2**4/24d0 &
                  + omd2w**4/sig2**5/120d0- omd2w**5/sig2**6/720d0 )
       wgtim(1:niw) = -we*cons*wt/(x**2)/pi
       wgtim(0)     =  we*sum(cons*expa_*wt/(x**2))/pi
       if(npm==2) wgtim(niw+1:2*niw) = cons*omd*wt(1:niw)/(x**2)/pi !Asymmetric contribution for
       !! --- Gaussian part. We use intrincic fortran function now. ------------
       aw = abs(ua_*we)
       aw2 = aw**2
       eee = we**2/sig2 
       wgtim(0)=wgtim(0)+ dsign(1d0,we)*.5d0*exp(aw2)*( erfc(sqrt(aw2+eee)) -erfc(aw) ) !2023feb
    endif
  end subroutine wintzsg_npm_wgtim
endmodule m_purewintz

subroutine wintzsg_npm_wgtim(npm,a,expa,we,esmr, wgtim) !Gaussian integral along im omg axis
!     loop over w' = (1-x)/x, frequencies in Wc(k,w')
!     {x} are gaussian-integration points between (0,1)
!---------------------------------------------------------------------
  use m_genallcf_v3,only: nx=>niw
  use m_readfreq_r,only: wt=>wwx,x=>freqx
  implicit none
  integer,intent(in):: npm!,nx
  real(8),intent(in):: expa(nx),we,esmr,a !,wt(nx) x(nx),
  real(8):: we2,weh,wel,weh2,wel2,cons,omd,omd2,rup,rdn,sss,sig2,omd2w
  real(8)    :: wintsf,errsum=0d0 !,derfcx,derfc
  integer :: ie,nav = 2000,i,ix,verbose
  real(8),parameter :: pi=3.1415926535897932d0, rmax=2d0 !rmax =2 is by ferdi. Is it optimum? See wintz
  real(8)::sig, smxowx, ee, omg, ww,cons1,cons2,xx,aw,eee,aw2
  logical :: test
  real(8):: wgtim(0:npm*nx)
  if(esmr==0d0) then
     call wintz_npm_wgtim(npm,x,wt,a,expa,we,nx, wgtim)
     return
  endif
  sig  = .5d0*esmr
  sig2 = 2d0*sig**2
  we2  = we**2
  wgtim=0d0
  do  i = 1,nx
     omd   = 1d0/x(i) - 1d0
     omd2w = omd**2 + we2
     ! pole weight is given by   1- exp (-R^2/sig2) = 1/N \int_0^R exp(- x^2/sig2) 2 pi r dr
     ! Gauss theorem---> but is this correct? Not three dimentional...
     if(omd2w/sig2 > 5d-3) then
        cons = (1d0 - exp (- omd2w/sig2))/omd2w
     else
        cons = ( 1d0/sig2 - omd2w/sig2**2/2d0 &
             + omd2w**2/sig2**3/6d0  - omd2w**3/sig2**4/24d0 &
             + omd2w**4/sig2**5/120d0- omd2w**5/sig2**6/720d0 )
     endif
     wgtim(i)= wgtim(i)+ we*cons*wt(i)/(x(i)**2)*(-1d0/pi)
     wgtim(0)= wgtim(0)- we*cons*expa(i)*wt(i)/(x(i)**2)*(-1d0/pi)
     if(npm==2) then !Asymmetric contribution for
        wgtim(i+nx) = wgtim(i+nx) - cons*omd*wt(i)/(x(i)**2)*(-1d0/pi)
     endif
  enddo
  !! --- Gaussian part -------------------------------------
  aw = abs(a*we)
  aw2 = aw**2
!  if(aw<6d0) then
     eee = we**2/sig2
     wgtim(0)=wgtim(0)+ dsign(1d0,we)*.5d0*exp(aw2)*( erfc(sqrt(aw2+eee)) -erfc(aw) )
!  else !April2004
!     wgtim(0)= wgtim(0) -.5d0* (2d0/(a*we)/pi)*sqrt(pi)/2d0* &
!          (1d0 - 1d0/2d0/aw2* (1d0 - 3d0/2d0/aw2* (1d0 - 5d0/2d0/aw2* &
!          (1d0 - 7d0/2d0/aw2*  (1d0 - 9d0/2d0/aw2* (1d0 - 11d0/2d0/aw2 ))))))
!  endif
end subroutine wintzsg_npm_wgtim
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine wintz_npm_wgtim(npm,x,wt,a,expa,we,nx, wgtim)
  ! takao complex version of wint by ferdi 92.02.24
  ! wintz = < [w'=-inf,inf] (i/2pi) v(w')/(w+w'-e) >
  ! the integration is done along the imaginary axis
  ! w' ==> iw', w' is now real
  ! wintz = - < [w'=0,inf] v(iw') (1/pi) (w-e)/{(w-e)^2 + w'^2} >
  ! transform: x = 1/(1+w')
  ! wintz = - < [x=0,1] v(iw') (1/pi) (w-e)/{(w-e)^2 + w'^2}x^2 >
  ! the integrand is peak around w'=0 or x=1 when we=w-e=0
  ! to handel the problem, add and substract the singular part as follows:
  ! wintz = - < [x=0,1] {v(iw') - v(0)exp(-a^2 w'^2)}
  !                   *(1/pi) (w-e)/{(w-e)^2 + w'^2}x^2 >
  !        - (1/2) v(0) sgn(w-e) exp(a^2 (w-e)^2) erfc(a|w-e|)
  ! the second term of the integral can be done analytically, which
  ! results in the last term
  ! when we=w-e ==> 0, (1/pi) (w-e)/{(w-e)^2 + w'^2} ==> delta(w')
  ! the integral becomes -v(0)/2
  ! v       = v(iw')
  ! v0      = v(0)
  ! x       = s.o.
  ! wt      = weights for integration
  ! a       = a constant determining the range of exp(-a^2 w'^2)
  ! expa(x) = exp(-a^2 w'^2)
  ! we      = w - e
  ! nx      = no. x points
  implicit real*8 (a-h,o-z)
  integer:: nx,i,npm
  real(8):: x(nx),wt(nx),expa(nx)
  real(8) :: rmax=2d0 !rmax =2 is by ferdi. Is it optimum? See wintz
  real(8),parameter:: pi=3.1415926535897932d0,tol=1d-8
  real(8):: wgtim(0:npm*nx)
  ! if w = e the integral = -v(0)/2 ! frequency integral
  if (dabs(we)< tol) call rx1( 'wintz: |w-e| < tol',we)
  we2  = we*we
  sum  = 0d0
  wgtim = 0d0
  if (dabs(we) < rmax/a) then
     do       i = 1,nx
        omd   = 1d0/x(i) - 1d0
        onemx      = 1.d0 - x(i)
        cons       = 1d0/(we2*x(i)*x(i) + onemx*onemx)
        wgtim(i)= wgtim(i)+ we*cons*wt(i)*(-1d0/pi)
        wgtim(0)= wgtim(0)+ we*cons*(-expa(i))*wt(i)*(-1d0/pi)
        if(npm==2) then     !Asymmetric contribution for
           wgtim(i+nx) = wgtim(i+nx) - cons*omd*wt(i)*(-1d0/pi)
        endif
     enddo
     wgtim(0)= wgtim(0)-0.5d0*dsign(1.d0,we)*dexp(we2*a*a)*erfc(a*dabs(we))
  else
     do       i = 1,nx
        omd   = 1d0/x(i) - 1d0 !this was missing. I added this at 25June2008
        onemx      = 1.d0 - x(i)
        cons       = 1d0/(we2*x(i)*x(i) + onemx*onemx)
        wgtim(i)= wgtim(i)+ we*cons*wt(i)*(-1d0/pi)
        if(npm==2) then     !Asymmetric contribution for
           wgtim(i+nx)= wgtim(i+nx) - cons* omd*wt(i)*(-1d0/pi)
        endif
     enddo
  endif
end subroutine wintz_npm_wgtim

!> gaussian smeared-pole version of wintz. ! takao developed from wintz by ferdi.
!! Assume that each eigenvale, and w-e as a result, has the width of esmr( Ry).
complex(8) function wintzsg_npm(npm,v,v0,x,wt,a,expa,we, nx,esmr)
  implicit none
  integer(4),intent(in)::npm,nx
  complex(8),intent(in) ::v(npm*nx),v0
  real(8),intent(in):: x(nx),wt(nx),expa(nx),we,esmr
  integer(4) :: i
  real(8):: a,we2,weh,wel,weh2,wel2,cons, &
       omd,omd2,rup,rdn,sss,sig2,omd2w
  complex(8) :: wwz,wintz_npm, img=(0d0,1d0),sum,   wintzsg1
  real(8)    :: wintsf,errsum=0d0!,derfcx,derfc
  integer(4) :: ie,nav = 2000
  real(8)    :: pi=3.1415926535897932d0, rmax=2d0
  ! max =2 is by ferdi. Is it optimum? See wintz
  real(8)::sig, smxowx, ee, omg, ww,cons1,cons2,xx,aw,eee,aw2
  integer(4) ix,verbose
  logical :: test
  complex(8)::sumgauss,sumgauss1,sumgauss2
  if(esmr==0d0) then
     wintzsg_npm = wintz_npm(npm,v,v0,x,wt,a,expa,we,nx)
     return
  endif
  sig  = .5d0*esmr
  sig2 = 2d0*sig**2
  we2  = we**2
  sum  = 0d0
  !! simple integration scheme.
  do  i = 1,nx
     omd   = 1d0/x(i) - 1d0
     omd2w = omd**2 + we2
     ! pole weight is given by   1- exp (-R^2/sig2) = 1/N \int_0^R exp(- x^2/sig2) 2 pi r dr
     ! Gauss theorem---> but is this correct? Not three dimentional...
     if(omd2w/sig2 > 5d-3) then
        cons = (1d0 - exp (- omd2w/sig2))/omd2w
     else
        cons = ( 1d0/sig2 - omd2w/sig2**2/2d0 &
             + omd2w**2/sig2**3/6d0  - omd2w**3/sig2**4/24d0 &
             + omd2w**4/sig2**5/120d0- omd2w**5/sig2**6/720d0 )
     endif
     sum  = sum + we*cons*( v(i)-v0*expa(i) ) *wt(i)/(x(i)**2)
     if(npm==2) then !Asymmetric contribution for
        sum  = sum - cons* v(i+nx)*omd*wt(i)/(x(i)**2)
     endif
  enddo
  !! --- Gaussian part -------------------------------------
  aw = abs(a*we)
  aw2 = aw**2
  if(aw<6d0) then
     eee = we**2/sig2
     sumgauss = dsign(1d0,we)*.5d0*v0*exp(aw2)*( erfc(sqrt(aw2+eee)) -erfc(aw) )
  else !April2004
     sumgauss =  -.5d0*v0* (2d0/(a*we)/pi)* &
          sqrt(pi)/2d0* &
          (1d0 - 1d0/2d0/aw2* &
          (1d0 - 3d0/2d0/aw2* &
          (1d0 - 5d0/2d0/aw2* &
          (1d0 - 7d0/2d0/aw2* &
          (1d0 - 9d0/2d0/aw2* &
          (1d0 - 11d0/2d0/aw2 &
          ))))))
  endif
  wintzsg_npm = -sum/pi  + sumgauss
  if(verbose()>90) &
       write(1116, "(' we sig wintzsg_npm =',2f8.3,4f14.6)" )we,sig,wintzsg_npm
END function wintzsg_npm
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
complex(8) function wintz_npm(npm,v,v0,x,wt,a,expa,we, nx) ! takao complex version of wint by ferdi
  implicit real*8 (a-h,o-z)
  integer:: nx,i,npm
  real(8):: x(nx),wt(nx),expa(nx)
  complex(8) ::v(nx),v0,sum,img=(0d0,1d0),wintz
  real(8) :: rmax=2d0 !rmax =2 is by ferdi. Is it optimum? See wintz
  data pi/3.1415926535897932d0/, tol/1.d-8/
  if (dabs(we) < tol) then
     print *, ' we=',we
     call rx( 'wintz: |w-e| < tol')
  endif
  we2        = we*we
  sum        = 0.d0
  if (dabs(we) < rmax/a) then
     do       i = 1,nx
        omd   = 1d0/x(i) - 1d0
        onemx      = 1.d0 - x(i)
        cons       = 1d0/(we2*x(i)*x(i) + onemx*onemx)
        sum        = sum + we*cons*(v(i) - v0*expa(i))*wt(i)
        if(npm==2) then !Asymmetric contribution for
           sum  = sum - cons* v(i+nx)*omd*wt(i)
        endif
     enddo
     wintz = -sum/pi - 0.5d0*v0*dsign(1.d0,we)*dexp(we2*a*a)*erfc(a*dabs(we))
  else
     do       i = 1,nx
        omd   = 1d0/x(i) - 1d0
        onemx      = 1.d0 - x(i)
        cons       = 1d0/(we2*x(i)*x(i) + onemx*onemx)
        sum        = sum + we*cons*v(i)*wt(i)
        if(npm==2) then !Asymmetric contribution for
           sum  = sum - cons* v(i+nx)*omd*wt(i)
        endif
     enddo
     wintz      = -sum/pi
  endif
  wintz_npm= wintz
  return
END function wintz_npm
! ssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine matzwz(zw,zmel, ntp0,nstate,ngb, zwz)
  implicit none
  integer(4) :: nstate,ntp0,itp,it,ngb
  complex(8) :: zw(ngb,ngb),zmel(ngb,nstate,ntp0),zwz(nstate,ntp0)
  complex(8), allocatable :: CC(:,:,:)
  allocate(CC(ngb,nstate,ntp0) )
  call matm(zw,zmel,cc, ngb, ngb, nstate*ntp0)
  do itp = 1,ntp0
     do  it = 1,nstate
        zwz(it,itp) = sum( dconjg(zmel(1:ngb,it,itp))*CC(1:ngb,it,itp))
     enddo
  enddo
  deallocate(CC)
end subroutine matzwz
subroutine alagr3z2wgt(x,xi, wgt)
  implicit none
  intent(in) ::          x,xi
  logical:: ieqj
  real(8) :: amatinv(3,3),amat(3,3),ratio,detxx,wgt(3),x,xi(3)
  if(x<0d0) call rx( ' alagr3z2: x<0d0')
  amat(1:3,1) = 1d0
  amat(1:3,2) = xi(1:3)**2
  amat(1:3,3) = xi(1:3)**4
  call minv33(amat,amatinv)
  wgt=matmul([1d0,x**2,x**4], amatinv)
end subroutine alagr3z2wgt
! sssssssssssssssssssssssssssssssssssssssssss
complex(8) function alagr3zz(x,xi,fi)
  ! even function version of alagr3z ! return the interpolated value on x for fi(xi).
  ! Imag part is corrected to be >0
  implicit none
  real(8)::  xi(3), amatinv(3,3),amat(3,3),x,detxx
  complex(8) :: fi(3)
  amat(1:3,1) = 1d0
  amat(1:3,2) = xi(1:3)**2
  amat(1:3,3) = xi(1:3)**4
  call minv33(amat,amatinv)
  alagr3zz=dcmplx( &
       sum (matmul(amatinv,dreal(fi)) * (/1d0,x**2,x**4/) ), &
       sum (matmul(amatinv,dimag(fi)) * (/1d0,x**2,x**4/) ) )
  if(dimag(alagr3zz)>0d0) alagr3zz = dcmplx( dreal(alagr3zz),0d0)
ENDfunction alagr3zz
! sssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine timeshow(info)
#ifdef __GPU
    use cudafor
#endif 
  character*(*) :: info
#ifdef __GPU
    integer :: ierr
    ierr = cudadevicesynchronize() 
#endif 
  write(6,'(a,$)')info
  call cputid(0)
  call flush(6)
end subroutine timeshow
