!> Fermi-Dirac smeared-pole version of wintz (2026-09-20; Gaussian 2004-2026). takao, from wintz by ferdi.
complex(8) function wintzsg_npm(npm,v,v0,x,wt,a,expa,we, nx,esmr)
  !! Imaginary-axis integral of one level for Sigma_c (one-shot GW, sxcf_fal2), Eq. 57 of PRB 76, 165106,
  !! for a level smeared by the Fermi-Dirac kernel of width esmr (= kBT [Ry] of [gw] t_sigmaw; the same
  !! kernel as the pole-term weights wfacx2 / weavx2).  we = (omega - e)/2 (Hartree) of the unsmeared level.
  !! The sharp-level integral wintz_npm is averaged over the kernel, e -> e + xk, xk = kBT ln(u/(1-u)),
  !! u = (j-1/2)/nqfd (the FD cumulative is the measure); its discontinuous half residue
  !! -v0 sign(we)/2 is averaged analytically, -v0 (2 Phi_FD(omega-e) - 1)/2, the smooth rest numerically.
  !! Levels farther than 30 kBT from omega see a smooth integrand and use the sharp formula at we.
  !! (Until 2026-09-20 this was the Gaussian regularization 'sig = esmr/2' of Eq. 57, which is the Gaussian
  !! level smearing of the esmr era; the same change was made in m_sxcf_sc.)
  use m_wfac, only: fd_cdf
  implicit none
  integer(4),intent(in)::npm,nx
  complex(8),intent(in) ::v(npm*nx),v0
  real(8),intent(in):: x(nx),wt(nx),expa(nx),we,esmr,a
  integer(4), parameter :: nqfd = 40
  integer(4) :: jq
  real(8) :: u, xk, wej, kbt_ha
  complex(8) :: wintz_npm, wintz_npm_smooth, sum
  if(esmr<=0d0 .or. abs(2d0*we) >= 30d0*esmr) then
     wintzsg_npm = wintz_npm(npm,v,v0,x,wt,a,expa,we,nx)
     return
  endif
  kbt_ha = .5d0*esmr                        ! we is in Hartree, esmr in Ry
  sum = 0d0
  do jq = 1, nqfd
     u   = (jq - .5d0)/nqfd
     xk  = kbt_ha*log(u/(1d0-u))            ! level shift (Hartree)
     wej = we - xk
     sum = sum + wintz_npm_smooth(npm,v,v0,x,wt,a,expa,wej,nx)
  enddo
  wintzsg_npm = sum/nqfd - .5d0*v0*(2d0*fd_cdf(2d0*we, esmr) - 1d0)   ! + the FD-averaged half residue
END function wintzsg_npm
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
complex(8) function wintz_npm_smooth(npm,v,v0,x,wt,a,expa,we, nx)
  !! wintz_npm without its half residue -v0 sign(we)/2: the remainder -v0 sign(we)/2 (exp(a^2 we^2) erfc(a|we|) - 1)
  !! is continuous at we = 0, so the level average of wintzsg_npm can quadrature it.
  implicit none
  integer(4),intent(in)::npm,nx
  complex(8),intent(in) ::v(npm*nx),v0
  real(8),intent(in):: x(nx),wt(nx),expa(nx),we
  real(8), intent(in) :: a
  complex(8) :: wintz_npm
  real(8) :: aw, rmax=2d0
  aw = a*abs(we)
  if (abs(we) < rmax/a) then
     wintz_npm_smooth = wintz_npm_nostep(npm,v,v0,x,wt,a,expa,we,nx) - 0.5d0*v0*dsign(1d0,we)*(dexp(aw*aw)*erfc(aw) - 1d0)
  else
     wintz_npm_smooth = wintz_npm(npm,v,v0,x,wt,a,expa,we,nx)
  endif
contains
  complex(8) function wintz_npm_nostep(npm,v,v0,x,wt,a,expa,we, nx)   ! numeric part of wintz_npm (|we| < rmax/a branch)
    integer(4),intent(in)::npm,nx
    complex(8),intent(in) ::v(npm*nx),v0
    real(8),intent(in):: x(nx),wt(nx),expa(nx),we,a
    real(8) :: pi=3.1415926535897932d0, we2, omd, onemx, cons
    complex(8) :: sum
    integer :: i
    we2 = we*we; sum = 0d0
    do i = 1,nx
       omd   = 1d0/x(i) - 1d0
       onemx = 1d0 - x(i)
       cons  = 1d0/(we2*x(i)*x(i) + onemx*onemx)
       sum   = sum + we*cons*(v(i) - v0*expa(i))*wt(i)
       if(npm==2) sum = sum - cons*v(i+nx)*omd*wt(i)
    enddo
    wintz_npm_nostep = -sum/pi
  end function wintz_npm_nostep
END function wintz_npm_smooth
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
!------------------------------------------------------------------------
complex(8) function alagr2zz(x,xi,fi)
  ! 2-point LINEAR interpolation in u=x^2 (even function), drop-in replacement for
  ! alagr3zz using only the bracketing pair xi(1)<=x<xi(2). Weights are the convex
  ! pair (1-t),t with t in [0,1] -> the result stays between fi(1) and fi(2):
  ! NO negative weights, NO overshoot, NO sign flip (unlike the 3-point Lagrange,
  ! whose end cardinals go negative and overshoot a curved/cusped Wc near omega=0).
  ! The imag clamp below is therefore redundant (kept for parity with alagr3zz).
  implicit none
  real(8):: xi(2), x, t, d
  complex(8):: fi(2)
  d = xi(2)**2 - xi(1)**2
  t = 0d0
  if(d/=0d0) t = (x**2 - xi(1)**2)/d
  if(t<0d0) t = 0d0
  if(t>1d0) t = 1d0
  alagr2zz = (1d0-t)*fi(1) + t*fi(2)
  if(dimag(alagr2zz)>0d0) alagr2zz = dcmplx( dreal(alagr2zz),0d0)
ENDfunction alagr2zz
! sssssssssssssssssssssssssssssssssssssssssssssssssssss
! subroutine timeshow(info)
! #ifdef __GPU
!     use cudafor
! #endif 
!   character*(*) :: info
! #ifdef __GPU
!     integer :: ierr
!     ierr = cudadevicesynchronize() 
! #endif 
!   write(6,'(a,$)')info
!   call cputid(0)
!   call flush(6)
! end subroutine timeshow
