!     loop over w' = (1-x)/x, frequencies in Wc(k,w')
!     {x} are gaussian-integration points between (0,1)
!---------------------------------------------------------------------
subroutine wintzsg_npm_wgtim(npm,a,expa,we,esmr, wgtim)
  use m_genallcf_v3,only: nx=>niw
  use m_readfreq_r,only: wt=>wwx,x=>freqx
  implicit none
  integer,intent(in):: npm!,nx
  real(8),intent(in):: expa(nx),we,esmr,a !,wt(nx) x(nx),
  real(8):: we2,weh,wel,weh2,wel2,cons,omd,omd2,rup,rdn,sss,sig2,omd2w
  real(8)    :: wintsf,errsum=0d0,derfcx,derfc
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
  if(aw<6d0) then
     eee = we**2/sig2
     wgtim(0)=wgtim(0)+ dsign(1d0,we)*.5d0*exp(aw2)*( derfcx(sqrt(aw2+eee)) -derfcx(aw) )
  else !April2004
     wgtim(0)= wgtim(0) -.5d0* (2d0/(a*we)/pi)*sqrt(pi)/2d0* &
          (1d0 - 1d0/2d0/aw2* (1d0 - 3d0/2d0/aw2* (1d0 - 5d0/2d0/aw2* &
          (1d0 - 7d0/2d0/aw2*  (1d0 - 9d0/2d0/aw2* (1d0 - 11d0/2d0/aw2 ))))))
  endif
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
     wgtim(0)= wgtim(0)-0.5d0*dsign(1.d0,we)*dexp(we2*a*a)*derfc(a*dabs(we))
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
  real(8)    :: wintsf,errsum=0d0,derfcx,derfc
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
     sumgauss = dsign(1d0,we)*.5d0*v0*exp(aw2)*( derfcx(sqrt(aw2+eee)) -derfcx(aw) )
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
  !$$$      return
  !$$$c---test code
  !$$$      print *,' sum= ', sum
  !$$$      print *,' awe1= ', a*we
  !$$$      print *,' erf1= ', derfcx(a*we)
  !$$$      print *,' awe2= ', sqrt(a**2+1d0/sig2)*we
  !$$$      print *,' erf2= ', derfcx(sqrt(a**2+1d0/sig2)*we)
  !$$$      print *,' dexp= ', dexp(we2*a**2)
  !$$$      print *,' erf3= ', derfcx(aw)
  !$$$      print *,' erf4= ', derfcx(sqrt(aw2+eee))
  !$$$      print *,' aw2= ', aw2
  !$$$      aw = a*we
  !$$$      write(6, "(' wintzsg_npm chk1=',2f8.3, 4d18.10)" )we,sig,
  !$$$     &  - .5d0*v0*  dexp(aw**2)*( derfcx(aw) - derfcx(sqrt(a**2+1d0/sig2)*we)  )
  !$$$      write(6, "(' wintzsg_npm chk2=',2f8.3, 4d18.10)" )we,sig,sumgauss1
  !$$$      write(6, "(' wintzsg_npm chk3=',2f8.3, 4d18.10)" )we,sig,sumgauss2
  !$$$      return
END function wintzsg_npm
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
complex(8) function wintz_npm(npm,v,v0,x,wt,a,expa,we, nx)
  ! takao complex version of wint by ferdi
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
     wintz = -sum/pi - 0.5d0*v0*dsign(1.d0,we)*dexp(we2*a*a)*derfc(a*dabs(we))
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
!$$$!sssssssssssssssssssssssssssssssssssssssssssssssssssssss
!$$$      subroutine matzwz2(iSigma_en ,zw,zmel, ntq,nstate,ngb, zwz)
!$$$      use m_mpi, only: mpi__rank
!$$$      implicit none
!$$$      integer(4), intent(in) :: iSigma_en
!$$$      complex(8), intent(in) :: zw(ngb,ngb),zmel(ngb,nstate,ntq)
!$$$      integer(4), intent(in) :: nstate,ntq,ngb
!$$$      complex(8), intent(out) :: zwz(nstate,ntq,ntq)
!$$$      integer(4) :: itp,itpp,it
!$$$      complex(8) :: zdotc
!$$$      complex(8), allocatable :: CC(:,:,:)
!$$$      complex(8),allocatable:: z1r(:,:),z2r(:,:), zwzi(:,:)
!$$$      integer :: ivc,verbose
!$$$      allocate(CC(ngb,nstate,ntq) )
!$$$      call matm(zw,zmel,cc, ngb, ngb, nstate*ntq) !most time consuming part
!$$$      if (iSigma_en==1.or.iSigma_en==5) then
!$$$        zwz=0d0
!$$$        do itp = 1,ntq
!$$$          do  it = 1,nstate
!$$$            zwz(it,itp,itp)=
!$$$     &    zdotc(ngb,zmel(1,it,itp),1,CC(1,it,itp),1 )
!$$$          enddo
!$$$        enddo
!$$$      elseif (iSigma_en==2 .or. iSigma_en==3 ) then
!$$$        if(verbose()>39)write(*,*)'info: USE GEMM FOR SUM (zwz=zmel*cc)'
!$$$        allocate(z1r(ntq,ngb),z2r(ngb,ntq),zwzi(ntq,ntq))
!$$$        do  it = 1, nstate
!$$$          do  itp = 1, ntq
!$$$            do  ivc = 1, ngb
!$$$              z1r(itp,ivc) = dconjg(zmel(ivc,it,itp))
!$$$              z2r(ivc,itp) = CC(ivc,it,itp)
!$$$            enddo
!$$$          enddo
!$$$          call zgemm('N','N',ntq,ntq,ngb,(1d0,0d0),z1r,ntq,
!$$$     .      z2r,ngb,(0d0,0d0),zwzi,ntq)
!$$$          do  itp = 1, ntq
!$$$            do itpp = 1, ntq
!$$$              zwz(it,itp,itpp) = zwzi(itp,itpp)
!$$$            enddo
!$$$          enddo
!$$$        enddo
!$$$        deallocate(z1r,z2r,zwzi)
!$$$      else
!$$$        call rx( "sxcf, matzwz2: iSigma_en /= 0,1,2,3")
!$$$      endif
!$$$      deallocate(CC)
!$$$      end subroutine matzwz2
!$$$!sssssssssssssssssssssssssssssssssssssssssssssssssssssss
complex(8) function alagr3zwgt (x,xi,wgt)
  implicit real*8 (a-h,o-z)
  ! three-point interpolation with arbitrary mesh
  ! f(x) = [ { (x-x2)(x-x3) } / { (x1-x2)(x1-x3) } ] f1
  !      + [ { (x-x1)(x-x3) } / { (x2-x1)(x2-x3) } ] f2
  !      + [ { (x-x1)(x-x2) } / { (x3-x1)(x3-x2) } ] f3
  !      = wg1(1)* f1+ wgt(2)*f2 +wgt(3)*f3
  ! x     = the point at which the function is to be interpolated
  ! xi(3) = points where the function is given
  real(8):: xi(3),wgt(3)
  xx1        = x-xi(1)
  xx2        = x-xi(2)
  xx3        = x-xi(3)
  x12        = xi(1)-xi(2)
  x13        = xi(1)-xi(3)
  x23        = xi(2)-xi(3)
  wgt(1)     =   xx2*xx3/(x12*x13)
  wgt(2)     =   xx1*xx3/(-x12*x23)
  wgt(3)     =   xx1*xx2/(x13*x23)
END function alagr3zwgt
! sssssssssssssssssssssssssssssssssssssssssss
!> even function version of alagr3z2wgt
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
  alagr3zz = dcmplx ( &
       sum (matmul(amatinv,dreal(fi)) * (/1d0,x**2,x**4/) ), &
       sum (matmul(amatinv,dimag(fi)) * (/1d0,x**2,x**4/) ) )
  if(dimag(alagr3zz)>0d0) alagr3zz = dcmplx( dreal(alagr3zz),0d0)
END function alagr3zz
! sssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine timeshow(info)
  character*(*) :: info
  write(*,'(a,$)')info
  call cputid(0)
end subroutine timeshow
