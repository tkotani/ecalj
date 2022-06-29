subroutine corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc, rfoc,z)
  use m_struc_def
  use m_lmfinit,only: pnux=>pnu,pzx=>pz
  use m_hansr,only: hansmr
  !- Returns parameters for smooth core+nucleus representation
  ! ----------------------------------------------------------------------
  !i Inputs
  !i     Elts read: lfoca rfoca qc z ctail etail stc orhoc lmxb pz rmt rg
  !i     Stored:    *
  !i     Passed to: *
  !i   is    :species index
  !o Outputs
  !o   cofg  :coefficient to Gaussian part of pseudocore density
  !o         :assigned so that pseudocore charge = true core charge
  !o   cofh  :coefficient to Hankel part of pseudocore density
  !o         :Hankel contribution is determined by inputs
  !o         :(qcorh,ceh,rfoc) and should accurately represent the
  !o         :true core core density for r>rmt
  !o   qcorg :charge in the gaussian part; see Remarks
  !o   qcorh :charge in the Hankel part; see Remarks
  !o   qsc   :number of electrons in semicore treated by local orbitals
  !o   lfoc  :switch specifying treatment of core density.
  !o          0 => val,slo = 0 at sphere boundary
  !o          1 => core tails included explicitly with valence
  !o          2 => tails included perturbatively
  !o
  !o   rfoc :smoothing radius for hankel head fitted to core tail
  !o   z     :nuclear charge
  !r Remarks
  !r   qcorg and qcorh are the charges in the gaussian and hankel parts.
  !r   The hankel part is used when the core is allowed to spill out of
  !r   the augmentation sphere.
  !r
  !r   cofg and cofh are the coefficients in front of the standard
  !r   gaussian and smoothed hankel functions for l=0.
  !r   That is: the pseudocore density is
  !r      cofg*g0(rg;r)*Y0 + cofh*h0(rfoca;r)*Y0        (1)
  !r   ceh and rfoc are the energy and sm.-radius for the hankel part.
  !r   cofg is set so that qc = integral of eq. 1 above.
  !r
  !r   For lfoc=0 there is no Hankel part; qc carried entirely by Gausian
  !r   For lfoc>0 there is no Hankel part; Gaussian carries difference
  !r              between qc and charge in Hankel part.
  !r
  !r   To add to the radial density 4*pi*r**2*rho_true, multiply
  !r   cofg,cofh by srfpi.
  !l Local variables
  !l    ccof :coefficient for core tail, for a smoothed Hankel.
  !l          ccof is differs from spec->ctail because ctail is
  !l          constructed for an unsmoothed Hankel.
  !u Updates
  !u   15 Sep 01 Generates qsc.  Argument list changed.
  !u   24 Apr 00 Adapted from nfp corpars.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: is,i_copy_size
  real(8):: qcorg , qcorh , qsc , cofg , cofh , ceh , rfoc , z
  type(s_spec)::sspec(is)
  integer:: n0 , lfoc , lmxb , l,isp
  parameter (n0=10)
  real(8):: pnu(n0),pz(n0),ccof,q0,q1,qc,rmt,rsm,stc,x0(0:n0), xi(0:n0),dgetss
  real(8),parameter:: fpi = 16d0*datan(1d0), srfpi = dsqrt(fpi), y0 = 1d0/srfpi
  lfoc=sspec(is)%lfoca
  rfoc=sspec(is)%rfoca
  qc=sspec(is)%qc
  z=sspec(is)%z
  ccof=sspec(is)%ctail
  ceh=sspec(is)%etail
  stc=sspec(is)%stc
  lmxb=sspec(is)%lmxb
  pnu= pnux(1:n0,1,is) !sspec(is)%p
  pz = pzx(1:n0,1,is)  !sspec(is)%pz
  rmt = sspec(is)%rmt
  if ( rfoc <= 1d-5 ) rfoc = (sspec(is)%rg)
  qsc = 0
  isp=1 !we assme int pz(:,1)=pz(:,2) int pnu as well
  do  l = 0, lmxb
     if (int(pz(l+1)) /= 0) then
        if (int(mod(pz(l+1),10d0)) < int(pnu(l+1))) qsc = qsc + 4*l+2
     endif
  enddo
  ! ... Scale smoothed hankel coeff for exact spillout charge
  !     q1 = spillout charge in sm. Hankel
  !     q0 = spillout charge in ordinary Hankel
  if (ccof /= 0) then
     call hansmr(rmt,0d0,1/rfoc,x0,1)
     call hansmr(rmt,ceh,1/rfoc,xi,1)
     q1 = srfpi/ceh*(-dexp(rfoc**2/4*ceh) &
          - rmt**3*(xi(1)-dexp(rfoc**2/4*ceh)*x0(1)))
     rsm = 0.05d0
     call hansmr(rmt,0d0,1/rsm,x0,1)
     call hansmr(rmt,ceh,1/rsm,xi,1)
     q0 = srfpi/ceh*(-dexp(rsm**2/4*ceh) &
          - rmt**3*(xi(1)-dexp(rsm**2/4*ceh)*x0(1)))
     q0 = q0*y0
     q1 = q1*y0
     ccof = ccof*q0/q1
  endif
  ! ... Set gaussian and hankel charges
  qcorg = qc
  qcorh = 0d0
  if (lfoc > 0) then
     qcorh = -ccof*dexp(ceh*rfoc*rfoc/4d0)/ceh
     qcorg = qc-qcorh
  endif
  ! ... Coeffients to the the gaussian and hankel terms
  cofh = -y0*qcorh*ceh*dexp(-ceh*rfoc*rfoc/4d0)
  cofg = y0*qcorg
  !      write (6,352) is,qcorg,qcorh,cofg,cofh
  !  352 format(' spec',i3,'  qcorg,qcorh=',2f10.6,'  cofg,cofh=',2f12.4)
end subroutine corprm


