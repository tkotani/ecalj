!> Occupation kernel of the self-energy: Fermi-Dirac of width kbt (Ry), [gw] t_sigmaw (K).
!!
!! Every intermediate level ek of Sigma_x = Gv and of the pole term of Sigma_c = G(W-v) is treated
!! as a level smeared by g(x) = -df/dx = f(1-f)/kbt (f the Fermi function of width kbt); the weight
!! of the level inside an energy window [el,eh] is Phi(eh-ek)-Phi(el-ek) with Phi the cumulative of g.
!! kbt = 0 is the sharp step (core levels).  Until 2026-09-20 this was a Gaussian of width esmr;
!! t_sigmaw = 1000 K is close to the old esmr = 0.01 Ry (kbt = 0.0063 Ry vs sigma = 0.01 Ry).
module  m_wfac
  implicit none
  public :: wfacx2, weavx2, fd_cdf, wcdf, wcut, sig_window, pole_weights
contains
  !> Fermi-Dirac cumulative (logistic) CDF: Phi(x)=1/(1+exp(-x/kbt)); occupation of a level at ek
  !! for a boundary at e is fd_cdf(e-ek,kbt).  kbt<=0: sharp step.
  pure real(8) function fd_cdf(x, kbt)
    real(8), intent(in) :: x, kbt
    real(8) :: t
    if(kbt<=0d0) then
       if(x>0d0) then; fd_cdf=1d0; elseif(x<0d0) then; fd_cdf=0d0; else; fd_cdf=0.5d0; endif
       return
    endif
    t = x/kbt
    if(t >  40d0) then; fd_cdf = 1d0
    elseif(t < -40d0) then; fd_cdf = 0d0
    else; fd_cdf = 1d0/(1d0+exp(-t)); endif
  end function fd_cdf
  !> same as fd_cdf (name kept for the callers that think of it as "the CDF of the kernel")
  pure real(8) function wcdf(x, kbt)
    real(8), intent(in) :: x, kbt
    wcdf = fd_cdf(x, kbt)
  end function wcdf
  !> half-width (Ry) beyond which the kernel is negligible: g(15 kbt)/g(0) ~ 3e-7.  0 for the sharp step.
  pure real(8) function wcut(kbt)
    real(8), intent(in) :: kbt
    wcut = 15d0*kbt
  end function wcut
  !> Half-width (Ry) of the window around E_F inside which intermediate levels can be partially
  !! occupied (batch sizing in sxcf_scz_count, the state window of the pole term, the real-axis
  !! omega range in getwemax): the kernel tail.
  pure real(8) function sig_window(kbt)
    real(8), intent(in) :: kbt
    sig_window = wcut(kbt)
  end function sig_window
  !> stable ln(cosh(z))
  pure real(8) function lncosh(z)
    real(8), intent(in) :: z
    real(8) :: za
    za = abs(z)
    lncosh = za + log(0.5d0*(1d0+exp(-2d0*za)))
  end function lncosh
  !> antiderivative of x*g(x): I(x) = 0.5*( x*tanh(x/(2kbt)) - 2*kbt*lncosh(x/(2kbt)) ), dI/dx = x*g(x).
  pure real(8) function fd_iav(x, kbt)
    real(8), intent(in) :: x, kbt
    real(8) :: a
    a = 2d0*kbt
    fd_iav = 0.5d0*( x*tanh(x/a) - a*lncosh(x/a) )
  end function fd_iav
  !> Weight of the level ek inside the window [e1,e2].
  pure real(8) function wfacx2(e1,e2, ek,kbt)
    real(8), intent(in) :: e1, e2, ek, kbt
    real(8), parameter :: ewidthcut=1d-6
    real(8) :: el, eh
    el=min(e1,e2); eh=max(e1,e2)
    wfacx2 = 0d0
    if(eh-el < ewidthcut) return
    wfacx2 = fd_cdf(eh-ek,kbt) - fd_cdf(el-ek,kbt)
  END function wfacx2
  !> Mean energy of the smeared level ek inside the window [e1,e2] (used when wcsmear=false).
  pure real(8) function weavx2(e1,e2, ek,kbt)
    real(8), intent(in) :: e1, e2, ek, kbt
    real(8) ::el,eh,wtt,xl,xh
    el=min(e1,e2); eh=max(e1,e2)
    if(kbt<=0d0) then; weavx2 = ek; return; endif
    xl=el-ek; xh=eh-ek
    wtt = fd_cdf(xh,kbt) - fd_cdf(xl,kbt)
    if(wtt < 1d-12) then        ! degenerate window: fall back to midpoint
       weavx2 = 0.5d0*(el+eh)
       return
    endif
    weavx2 = ek + ( fd_iav(xh,kbt) - fd_iav(xl,kbt) )/wtt
  END function weavx2
  !> Shared by the QSGW pole term (m_sxcf_sc) and the one-shot one (sxcf_fal2).
  subroutine pole_weights(omg, ef, ek, esmr, smear, nw, freq_r, wfac, scale, iw1, iw2, wts)
    ! Weights of one (it,itp) pair of the real-axis pole term of Sigma_c on the W mesh points
    ! iw1:iw2 (wts(iw), sum = wfac*scale), or iw2 < iw1 when the pair contributes nothing.
    ! The intermediate level ek is smeared by the Fermi-Dirac kernel of width esmr (= kBT of
    ! [gw] t_sigmaw, m_wfac); wfac is its weight inside the window between ef and omg (PRB 76, 165106)
    ! and scale carries the k-point weight and the sign of (omg - ef).
    !  smear=.false. (default): W_c is taken at the mean energy omega_bar = |omg - weavx2|/2 (Eq. 58),
    !     interpolated on the 3 bracketing mesh points (omega^2 Lagrange; 2-point linear with --WVR2ptRaxis).
    !  smear=.true.  ([gw] wcsmear, default): the kernel is applied to W_c(omega) itself.  The level e' runs over
    !     the window with the kernel g(e'-ek), omega = |omg - e'|/2 (Hartree); mesh point iw owns
    !     omega in [(freq_r(iw-1)+freq_r(iw))/2, (freq_r(iw)+freq_r(iw+1))/2] and gets
    !     Phi(b-ek)-Phi(a-ek) over the e' interval mapped from that cell, clipped to the window.
    !  esmr=0 (core levels): sharp step, always the mean-energy path.
    use m_cmdopt_registry, only: c0_WVR2ptRaxis
    real(8), intent(in) :: omg, ef, ek, esmr, wfac, scale
    logical, intent(in) :: smear
    integer, intent(in) :: nw
    real(8), intent(in) :: freq_r(0:nw)
    integer, intent(out) :: iw1, iw2
    real(8), intent(out) :: wts(0:nw)
    real(8) :: el, eh, sgn, wa, wb, a, b, t, cut, x, xi(3), amat(3,3), tt2p, w3(3)
    integer :: iw, ia, ib, ixs
    wts = 0d0; iw1 = 0; iw2 = -1
    if (.not. smear .or. wcut(esmr) <= 0d0) then          ! ---- W_c at the mean energy (Eq. 58)
      x   = .5d0*abs(omg - weavx2(omg, ef, ek, esmr))     ! \bar{omega_epsilon}
      ixs = findloc(freq_r(1:nw) > x, value=.true., dim=1)
      if (ixs < 1 .or. ixs > nw-1) return                 ! findloc miss (0) or beyond the mesh; ixs=1 (bin 0 = static W) is valid
      xi = freq_r(ixs-1:ixs+1)
      if (c0_WVR2ptRaxis) then                            ! 2-point linear in omega^2: convex, no overshoot
        tt2p = max(0d0, min(1d0, (x**2 - xi(1)**2)/(xi(2)**2 - xi(1)**2)))
        w3 = [1d0-tt2p, tt2p, 0d0]
      else                                                ! 3-point Lagrange in omega^2 (alagr3zz)
        amat(1:3,1) = 1d0; amat(1:3,2) = xi**2; amat(1:3,3) = xi**4
        w3 = matmul([1d0, x**2, x**4], inverse33(amat))
      endif
      iw1 = ixs-1; iw2 = ixs+1
      wts(iw1:iw2) = wfac*scale*w3
      return
    endif
    el = min(omg, ef); eh = max(omg, ef)                  ! ---- kernel-integrated weights ([gw] wcsmear)
    cut = wcut(esmr)
    el = max(el, ek - cut); eh = min(eh, ek + cut)        ! where the kernel is non-negligible
    if (eh <= el) return
    sgn = merge(1d0, -1d0, omg < ef)                      ! e' = omg + sgn*2*omega lies between omg and ef
    wa = .5d0*min(abs(omg-el), abs(omg-eh)); wb = .5d0*max(abs(omg-el), abs(omg-eh))
    ia = 0; ib = nw
    do iw = 1, nw                                         ! cells reaching [wa, wb]
      if (.5d0*(freq_r(iw-1)+freq_r(iw)) <= wa) ia = iw
      if (.5d0*(freq_r(iw-1)+freq_r(iw)) >= wb) then; ib = iw - 1; exit; endif
    enddo
    do iw = ia, ib
      a = omg + sgn*2d0*merge(0d0, .5d0*(freq_r(iw-1)+freq_r(iw)), iw == 0)
      b = omg + sgn*2d0*merge(freq_r(nw), .5d0*(freq_r(iw)+freq_r(iw+1)), iw == nw)
      if (a > b) then; t = a; a = b; b = t; endif
      a = max(a, el); b = min(b, eh)
      if (b <= a) cycle
      t = wcdf(b - ek, esmr) - wcdf(a - ek, esmr)   ! kernel weight of this cell; these sum to wfac
      if (t < 1d-12) cycle
      wts(iw) = scale*t
      if (iw2 < iw1) iw1 = iw
      iw2 = iw
    enddo
  end subroutine pole_weights
  pure function crossf(a,b) result(c)
    implicit none
    intent(in):: a,b
    real(8):: a(3),b(3),c(3)
    c(1)=a(2)*b(3)-a(3)*b(2)
    c(2)=a(3)*b(1)-a(1)*b(3)
    c(3)=a(1)*b(2)-a(2)*b(1)
  end function crossf
  pure function inverse33(matrix) result(inverse) !Inverse of 3X3 matrix
    implicit none
    real(8),intent(in) :: matrix(3,3)
    real(8) :: inverse(3,3), det
    inverse(:,1)= crossf(matrix(:,2),matrix(:,3))
    inverse(:,2)= crossf(matrix(:,3),matrix(:,1))
    inverse(:,3)= crossf(matrix(:,1),matrix(:,2))
    det = sum(matrix(:,1)*inverse(:,1))
    inverse = transpose(inverse)
    inverse = inverse/det
  end function inverse33
end module m_wfac

!> Standalone occupation weight (kept as external function for existing callers).
pure real(8) function wfacx(el,eh, ek,kbt)
  use m_wfac, only: fd_cdf
  implicit none
  real(8), intent(in) :: el,eh,ek,kbt
  wfacx = fd_cdf(eh-ek,kbt) - fd_cdf(el-ek,kbt)
end function wfacx
