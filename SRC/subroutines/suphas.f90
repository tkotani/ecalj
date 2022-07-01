subroutine suphas(q,p,ng,iv,n1,n2,n3,qlat,cosgp,singp)
  !- Makes exp(-i p* (q+G)) for a list of reciprocal lattice vectors
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   q     :Bloch wave number to be added to G
  !i   p     :position
  !i   ng    :number of G-vectors
  !i   iv    :G-vectors, as integer multiples of qlat
  !i   n1    :maximum first component of iv
  !i   n2    :maximum second component of iv
  !i   n3    :maximum third component of iv
  !i   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
  !o Outputs
  !i   cosgp :cos(-i p (q+G))
  !i   singp :sin(-i p (q+G))
  !r Remarks
  !r   iv may be calculated from setup suphs0.
  !u Updates
  !u   19 May 00 Adapted from nfp su_phs.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: ng,iv(ng,3),n1,n2,n3
  double precision :: p(3),cosgp(ng),singp(ng),qlat(3,3),q(3),gg(3,ng)
  complex(8):: cc(ng),img=(0d0,1d0)
  real(8):: tpi = 8d0*datan(1d0)
  call tcn('suphas')
  cc = exp(-img*tpi*sum(p*q)) * exp(-img*tpi*matmul(p, matmul(qlat, transpose(iv))))
  cosgp=dreal(cc)
  singp=dimag(cc)
  call tcx('suphas')
end subroutine suphas
subroutine suphs0(plat,ng,gv,iv)
  integer :: ng,iv(ng,3)
  double precision :: gv(ng,3),plat(3,3)
  iv = nint(matmul(gv,plat))
end subroutine suphs0

