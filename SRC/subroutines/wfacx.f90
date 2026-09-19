!> Smearing weights for the self-energy occupation.
!! Gaussian smearing (width esmr): wfacx(x) = \int_el^eh 1/sqrt(2)/esmr exp(-(x-ek/esmr)**2)).
!!
!! 2026-06-15 (t_sigmakbt): optional Fermi-Dirac smearing for the self-energy occupation.
!! When the Sigma-side finite-T is enabled (set_sigma_fd(.true.,kbt) called once before the
!! sxcf loops), wfacx/wfacx2/weavx2 use a Fermi-Dirac thermal kernel at temperature kbt[Ry]
!! instead of the Gaussian(esmr) kernel. This makes the FS occupation in Gv (exchange) and
!! G(W-v) (correlation) a true f(eps;T) at the finite-T chemical potential (EFERMI_kbt),
!! parallel to tetrakbt on the chi0 side. Default (sig_fd=.false.) reproduces the exact
!! Gaussian behaviour. T->0 (kbt->0) reduces to the sharp step.
module  m_wfac
  implicit none
  public :: wfacx2, weavx2, set_sigma_fd, fd_cdf, wcdf, wcut
  ! --- Sigma-side finite-T (Fermi-Dirac) state, set once via set_sigma_fd ---
  logical, save, public :: sig_fd  = .false.  ! .true. -> use Fermi-Dirac kernel at sig_kbt
  real(8), save, public :: sig_kbt = 0d0      ! electronic temperature in Ry (kBT)
contains
  !> Enable/disable the Fermi-Dirac self-energy smearing and set its temperature [Ry].
  subroutine set_sigma_fd(on, kbt)
    logical, intent(in) :: on
    real(8), intent(in) :: kbt
    sig_fd  = on
    sig_kbt = kbt
  end subroutine set_sigma_fd
  !> Fermi-Dirac cumulative (logistic) CDF: Phi(x)=1/(1+exp(-x/kbt)).
  !! Occupation of a level at ek for a boundary at e is fd_cdf(e-ek,kbt).
  pure real(8) function fd_cdf(x, kbt)
    real(8), intent(in) :: x, kbt
    real(8) :: t
    if(kbt<=0d0) then          ! T->0 : sharp step
       if(x>0d0) then; fd_cdf=1d0; elseif(x<0d0) then; fd_cdf=0d0; else; fd_cdf=0.5d0; endif
       return
    endif
    t = x/kbt
    if(t >  40d0) then; fd_cdf = 1d0
    elseif(t < -40d0) then; fd_cdf = 0d0
    else; fd_cdf = 1d0/(1d0+exp(-t)); endif
  end function fd_cdf
  !> stable ln(cosh(z))
  pure real(8) function lncosh(z)
    real(8), intent(in) :: z
    real(8) :: za
    za = abs(z)
    lncosh = za + log(0.5d0*(1d0+exp(-2d0*za)))
  end function lncosh
  !> antiderivative of x*g(x) where g=Phi' is the FD thermal kernel of width kbt:
  !! I(x) = 0.5*( x*tanh(x/(2kbt)) - 2*kbt*lncosh(x/(2kbt)) ),  dI/dx = x*g(x).
  pure real(8) function fd_iav(x, kbt)
    real(8), intent(in) :: x, kbt
    real(8) :: a
    a = 2d0*kbt
    fd_iav = 0.5d0*( x*tanh(x/a) - a*lncosh(x/a) )
  end function fd_iav

  !> CDF of the occupation kernel used by wfacx2: boundary at x = e - ek (Gaussian esmr, or FD when sig_fd).
  pure real(8) function wcdf(x, esmr)
    real(8), intent(in) :: x, esmr
    if (sig_fd) then
       wcdf = fd_cdf(x, sig_kbt)
    elseif (esmr == 0d0) then
       wcdf = merge(1d0, 0d0, x >= 0d0)
    else
       wcdf = 0.5d0*erfc(-x/sqrt(2d0)/esmr)
    endif
  end function wcdf
  !> half-width (Ry) beyond which the kernel of wcdf is negligible (0 for the sharp step).
  pure real(8) function wcut(esmr)
    real(8), intent(in) :: esmr
    if (sig_fd) then
       wcut = 30d0*sig_kbt
    else
       wcut = 6d0*esmr
    endif
  end function wcut
  pure real(8) function wfacx2(e1,e2, ek,esmr)
    real(8), intent(in) :: e1, e2, ek, esmr
    real(8) ::el,eh
    real(8),parameter :: ewidthcut=1d-6
    el=min(e1,e2)  !May2006
    eh=max(e1,e2)
    if(eh-el< ewidthcut) then !July2006
       wfacx2=0d0
       return
    endif
    if(sig_fd) then
       wfacx2 = fd_cdf(eh-ek,sig_kbt) - fd_cdf(el-ek,sig_kbt)
       return
    endif
    if(esmr==0d0) then
       wfacx2=0d0
       if(el <= ek .AND. ek <eh ) wfacx2=1d0
       return
    endif
    wfacx2 = .5d0*erfc(-(eh-ek)/sqrt(2d0)/esmr) - .5d0*erfc(-(el-ek)/sqrt(2d0)/esmr)
  END function wfacx2
  !! Averaged energy in window[el, eh] for the smearing kernel centred on ek.
  pure real(8) function weavx2(e1,e2, ek,esmr)
    real(8), intent(in) :: e1, e2, ek, esmr
    real(8) ::el,eh,wtt,sig2,xl,xh
    real(8),parameter:: pi=3.1415926535897932d0
    el=min(e1,e2)
    eh=max(e1,e2)
    if(sig_fd) then
       xl=el-ek; xh=eh-ek
       wtt = fd_cdf(xh,sig_kbt) - fd_cdf(xl,sig_kbt)
       if(wtt < 1d-12) then        ! degenerate window: fall back to midpoint
          weavx2 = 0.5d0*(el+eh)
          return
       endif
       weavx2 = ek + ( fd_iav(xh,sig_kbt) - fd_iav(xl,sig_kbt) )/wtt
       return
    endif
    wtt=    0.5d0*erfc(-(eh-ek)/sqrt(2d0)/esmr) -0.5d0*erfc(-(el-ek)/sqrt(2d0)/esmr)
    sig2= 2d0*esmr**2
    weavx2 = ek+ esmr/sqrt(2d0*pi) *( -exp(-(eh-ek)**2/sig2) + exp(-(el-ek)**2/sig2) )/wtt
  END function weavx2
end module m_wfac

!> Standalone exchange occupation weight (kept as external function for existing callers).
pure real(8) function wfacx(el,eh, ek,esmr)
  use m_wfac, only: sig_fd, sig_kbt, fd_cdf
  implicit none
  real(8), intent(in) :: el,eh,ek,esmr
  if(sig_fd) then
     wfacx = fd_cdf(eh-ek,sig_kbt) - fd_cdf(el-ek,sig_kbt)
     return
  endif
  if(esmr==0d0) then
     wfacx=0d0
     if(el <= ek .AND. ek <eh ) wfacx=1d0
     return
  endif
  wfacx = .5d0*erfc(-(eh-ek)/sqrt(2d0)/esmr) - .5d0*erfc(-(el-ek)/sqrt(2d0)/esmr)
end function wfacx
