!> Contour-decomposition test of Sigma_c for ONE intermediate level and a model W_c (matrix elements = 1).
!!
!! Model W_c(z) = -v w_p^2 / (w_p^2 - z^2 - i gam z)  (Drude-Lorentz, analytic in the upper half plane;
!! W_c(i w') = -v w_p^2/(w_p^2 + w'^2 + gam w') is real and has a LINEAR term in w' like a metal).
!!
!! Code path under test (same routines / same formulae as hgw):
!!   P  : m_wfac::pole_weights (window weight Phi, [gw] wcsmear cell weights) x Re W_c on the real-axis mesh
!!   I  : the imaginary-axis weights of m_sxcf_sc (CorrelationSelfEnergyImagAxis, copied verbatim below) on
!!        niw Gauss-Legendre points, W_c(i w') exact.
!! Reference: I_ex(w_e) = -(1/pi) int_0^inf dw' w_e W_c(i w')/(w_e^2+w'^2) by adaptive quadrature and
!!            P_ex = s int_{el}^{eh} de g(e-e') Re W_c(|omg-e|/2) by fine quadrature, both FD-averaged
!!            (the sharp-level contour identity is exact for a W_c with a spectral representation).
!! Sweep omg across e' and print I, P, I+P for code and reference.
!! Units: energies (omg, e', ef, kbt) in Ry as in hgw; W_c arguments (w_e, w', nu) in Hartree.
!! usage: contour_test [niw] [kbt_K] [gam_eV] [wcsmear 0/1]
program contour_test
  use m_wfac, only: pole_weights, fd_cdf, wfacx2
  implicit none
  integer :: niw, nw, nwhis, iw, iw1, iw2, jq, nq, i, io, nqfd
  real(8) :: kbt, gam, wp, v, ef, ep, omg, esmr, we, aw, aw2, u, xk, estep, wfac, scale, ua_, pi, rmax
  real(8) :: aa, bb, ratio, dw, omg2max, rydberg, hartree, tK, ge, isum, psum, iex, pex, x1, s
  real(8), allocatable :: freqx(:), wwx(:), freqw(:), expa_(:), cons(:), wgtim(:), wts(:), freq_r(:), frhis(:)
  real(8), allocatable :: wci(:), wcr(:)
  logical :: smear
  character(64) :: arg
  pi = 4d0*atan(1d0); rydberg = 13.605693d0; hartree = 27.211386d0; rmax = 2d0; ua_ = 1d0
  niw = 10; tK = 1000d0; gam = 0.1d0; smear = .true.
  if (command_argument_count() >= 1) then; call get_command_argument(1, arg); read(arg,*) niw; endif
  if (command_argument_count() >= 2) then; call get_command_argument(2, arg); read(arg,*) tK; endif
  if (command_argument_count() >= 3) then; call get_command_argument(3, arg); read(arg,*) gam; endif
  if (command_argument_count() >= 4) then; call get_command_argument(4, arg); read(arg,*) i; smear = i /= 0; endif
  kbt = 8.6171d-5*tK/rydberg          ! Ry
  esmr = kbt
  wp = 1.8d0/hartree; gam = gam/hartree; v = 1d0   ! W_c model in Hartree
  ef = 0d0                             ! Fermi level (Ry)
  ep = 0.02d0                          ! intermediate level e' = +0.27 eV above E_F (unoccupied, like t2g)
  ! ---- imaginary-axis GL mesh as m_readfreq_r
  allocate(freqx(niw), wwx(niw), freqw(niw), expa_(niw), cons(niw), wgtim(0:niw), wci(0:niw))
  call gauss(niw, 0d0, 1d0, freqx, wwx)
  do iw = 1, niw
    freqw(iw) = (1d0 - freqx(iw))/freqx(iw)
    expa_(iw) = exp(-ua_**2*freqw(iw)**2)
    wci(iw) = wc_imag(freqw(iw))
  enddo
  wci(0) = wc_imag(0d0)
  ! ---- real-axis histogram mesh as m_freq (HistBin_ratio 1.03, HistBin_dw 1e-5 Ry -> Hartree/2?)  m_freq works in Ry? frhis is
  !      in Hartree ("omg is in unit of Hartree"); dw_in = HistBin_dw (Ry) is used as is.  We copy: dw=1e-5, ratio=1.03.
  dw = 1d-5; ratio = 1.03d0; aa = ratio - 1d0; bb = dw/aa; omg2max = 1.0d0   ! up to 27 eV
  iw = 0
  do; iw = iw + 1; if (bb*(exp(aa*(iw-1)) - 1d0) > omg2max) exit; enddo
  nwhis = iw + 2
  allocate(frhis(nwhis+1)); do iw = 1, nwhis+1; frhis(iw) = bb*(exp(aa*(iw-1)) - 1d0); enddo
  nw = nwhis - 1
  allocate(freq_r(0:nw), wcr(0:nw), wts(0:nw))
  freq_r(0) = 0d0; freq_r(1:nw) = (frhis(1:nw) + frhis(2:nw+1))/2d0
  do iw = 0, nw; wcr(iw) = wc_real(freq_r(iw)); enddo
  write(*,'("# contour_test: niw=",i3,"  T=",f7.1," K (kbt=",es10.3," Ry)  gam=",f6.3," eV  wp=1.8 eV  wcsmear=",l1,"  nw=",i4)') &
       niw, tK, kbt, gam*hartree, smear, nw
  write(*,'(a,f7.3,a,es9.2,a,es9.2,a)') "# e'-E_F = ", ep*rydberg, " eV; W mesh spacing near 0: ", &
       freq_r(2)-freq_r(1), " Ha; smallest imag. w' = ", minval(freqw), " Ha"
  write(*,'(a)') "# omg-e'(eV)   I_code      P_code      Sum_code   |  I_ref       P_ref       Sum_ref   | dI   dP   dSum"

  nqfd = 40
  do io = -60, 60
    omg = ep + io*0.0025d0        ! sweep +-0.15 Ry = +-2 eV around e' in 34 meV steps
    ! ================= code path: imaginary axis (m_sxcf_sc block, esmr>0, valence level)
    nq = 1; if (abs(omg - ep) < 30d0*esmr) nq = nqfd
    wgtim = 0d0
    do jq = 1, nq
      xk = 0d0
      if (nq > 1) then; u = (jq - .5d0)/nq; xk = esmr*log(u/(1d0-u)); endif
      we = .5d0*(omg - ep - xk)
      aw = abs(ua_*we); aw2 = aw*aw
      cons = 1d0/(we**2*freqx**2 + (1d0-freqx)**2)
      block
        real(8) :: w1(0:niw)
        w1(1:niw) = we*cons*wwx*(-1d0/pi)
        estep = dsign(1d0,we)*dexp(aw2)*erfc(aw)
        if (nq > 1) estep = estep - dsign(1d0,we)
        w1(0) = merge(-sum(w1(1:niw)*expa_) - 0.5d0*estep, 0d0, mask=dabs(we) < rmax/ua_)
        wgtim = wgtim + w1/nq
      end block
    enddo
    if (nq > 1) wgtim(0) = wgtim(0) - 0.5d0*(2d0*fd_cdf(omg - ep, esmr) - 1d0)
    isum = sum(wgtim*wci)
    ! ================= code path: pole term (pole_weights, scale = sign(omg-ef), wkkr = 1)
    wfac = wfacx2(omg, ef, ep, esmr)
    scale = dsign(1d0, omg - ef)
    psum = 0d0
    if (wfac >= 1d-10) then
      call pole_weights(omg, ef, ep, esmr, smear, nw, freq_r(0:nw), wfac, scale, iw1, iw2, wts)
      if (iw2 >= iw1) psum = sum(wts(iw1:iw2)*wcr(iw1:iw2))
    endif
    ! ================= reference (exact quadratures, FD-averaged)
    iex = 0d0; pex = 0d0
    do jq = 1, 400                       ! FD average over e = e' + x, x = kbt ln(u/(1-u)), u midpoint
      u = (jq - .5d0)/400; xk = esmr*log(u/(1d0-u))
      iex = iex + i_exact(.5d0*(omg - ep - xk))/400
    enddo
    s = dsign(1d0, omg - ef)
    x1 = min(omg, ef); ge = max(omg, ef)         ! window [el, eh]
    do jq = 1, 4000
      u = (jq - .5d0)/4000; xk = esmr*log(u/(1d0-u))   ! level at e = e' + xk with measure g
      if (ep + xk > x1 .and. ep + xk < ge) pex = pex + s*wc_real(.5d0*abs(omg - (ep + xk)))/4000
    enddo
    write(*,'(f10.3,3f12.5,2x,3f12.5,2x,3f9.5)') (omg-ep)*rydberg, isum, psum, isum+psum, iex, pex, iex+pex, &
         isum-iex, psum-pex, isum+psum-iex-pex
  enddo
contains
  real(8) function wc_imag(w)     ! W_c(i w'), w' in Hartree
    real(8), intent(in) :: w
    wc_imag = -v*wp**2/(wp**2 + w**2 + gam*w)
  end function wc_imag
  real(8) function wc_real(nu)    ! Re W_c(nu), nu >= 0 in Hartree
    real(8), intent(in) :: nu
    complex(8) :: z
    z = -v*wp**2/(wp**2 - nu**2 - (0d0,1d0)*gam*nu)
    wc_real = dble(z)
  end function wc_real
  real(8) function i_exact(we)    ! -(1/pi) int_0^inf dw' we W_c(i w')/(we^2+w'^2), substitution w' = |we| tan(t)
    real(8), intent(in) :: we
    integer :: k, n
    real(8) :: t, dt, wpr, sgn, acc
    n = 20000; dt = (pi/2)/n; acc = 0d0; sgn = dsign(1d0, we)
    if (abs(we) < 1d-12) then; i_exact = 0d0; return; endif
    do k = 1, n
      t = (k - .5d0)*dt; wpr = abs(we)*tan(t)
      acc = acc + wc_imag(wpr)*dt        ! we/(we^2+w'^2) dw' = sgn dt
    enddo
    i_exact = -sgn*acc/pi
  end function i_exact
end program contour_test
