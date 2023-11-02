module m_madmat
  contains
!      subroutine madmat(nbas,bas,awald,alat,vol,dlat,nkd,glat,nkg,dmad)
!        call madmat ( nbas , rv_a_opos , lat_awald , lat_alat , lat_vol , rv_a_odlv
!     .      , lat_nkd , rv_a_oqlv , lat_nkq , rv_a_omad )
subroutine madmat(dmad)
  use m_lgunit,only:stdo
  use m_lmfinit,only:nbas,alat=>lat_alat
  use m_lattic,only: awald=>lat_awald,vol=>lat_vol,dlat=>rv_a_odlv,nkd=>lat_nkd,glat=>rv_a_oqlv, nkg=>lat_nkq,bas=>rv_a_opos
  use m_shortn3_plat,only: shortn3_plat,nout,nlatout
  use m_lattic,only: qlat=>lat_qlat,plat=>lat_plat
  !- Coefficients to Madelung matrix
  ! ----------------------------------------------------------------
  !i Inputs
  !i   nbas,bas,awald,alat,vol,dlat,nkd,glat,nkg
  !o Outputs
  !o   dmad
  !r Remarks
  !r   The potential from a basis of more than one atom is obtained by
  !r   superposing the potential from each sublattice.  Matrix element
  !r   dmad(i,j) is the (1/2)potential at position tau(i) from unit
  !r   charges on sublattice j, compensated by a uniform background of
  !r   density 1/vol.  (if tau(i)=0, the potential from the singularity
  !r   0 is not included.)
  !r
  !r   Call lattc to generate inputs awald,vol,dlat,nkd,glat,nkg.
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  !      integer nbas,nkd,nkg
  !      double precision awald,alat,vol
  !      double precision bas(3,nbas),dlat(3,nkd),glat(3,nkg)
  double precision :: dmad(nbas,nbas)
  ! local variables
  integer :: ibas,jbas,iprint,i
  double precision :: tau(3)
  real(8):: pp(3)
  dmad=0.0d0
  ! --- Generate Madelung matrix ---
  do  101  ibas = 1, nbas
     do  10  jbas = ibas, nbas
        do  15  i = 1, 3
           tau(i) = bas(i,jbas)-bas(i,ibas)
15      enddo
        !call shortn(tau,tau,dlat,nkd)
        pp = matmul(transpose(qlat),tau)
        call shortn3_plat(pp)
        tau = matmul(plat,pp+nlatout(:,1))
        call strx00(tau,awald,alat,vol,glat,nkg,dlat,nkd, &
             dmad(ibas,jbas))
        dmad(jbas,ibas) = dmad(ibas,jbas)
10   enddo
101 enddo
  ! --- Printout ---
  if (iprint() < 90) return
  write(stdo,100) nbas
100 format(/' Madelung matrix: nbas=',i5)
  do  20  ibas = 1, nbas
     print 110, (dmad(ibas,jbas),jbas=1,nbas)
110  format(7f11.7)
20 enddo
end subroutine madmat
subroutine strx00(tau,awald,alat,vol,glat,nkg,dlat,nkd,dl)
  !- Widget to make structure constant DL for L=0,E=0,K=0.
  ! ----------------------------------------------------------------
  !i Inputs
  !i   TAU,awald,ALAT,
  !i   VOL:     cell volume
  !i   glat,nkg:reciprocal lattice vectors and number
  !i   DLAT,NKD:real (direct) space lattice vectors and number
  !o Outputs
  !o   dl: 'potential' Phi at tau (see remarks below)
  !r Remarks
  !r   dl is in inverse atomic units.
  !r
  !r  Formally, the kappa=0 structure constants H_L
  !r  may be evaluated as the limit as e-> 0 of
  !r  F_L = H_L(e,\r) - C(e)\delta(L,0), with C(e) = -4\pi/\omega Y0 1/e
  !r  This latter term subtracts of a constant background which
  !r  makes int_cell ( ) cell vanish.  As e->0, both terms
  !r  diverge, but the difference remains finite.  In that limit,
  !r  \nabla F_0 = -4\pi\delta_{L,0}(r) + 4\pi/\omega Y_0; thus
  !r  the F_0 structure constant corresponds \nabla^{-1} of a
  !r  superposition of Y0 delta(r-R), compensated by a
  !r  uniform background of density Y0/omega.  This routine
  !r  returns dl=F_0/Y_0, which is half the potential of a superposition
  !r  of unit charges at lattice points R.
  !r
  !r  Evaluation of the Ewald sum:
  !r  The potential is broken into a "damped" and a "smooth" part as
  !r  (1) h(r) = 1/r = erfc(a*r)/r + erf(a*r)/r = h^d(r) + h^s(r)
  !r  which is independent of the Ewald parameter a, and the sum
  !r  (2) Phi(r) = \sum_R h(r-R) - uniform background
  !r  is calculated by adding contributions from h^d and h^s separately;
  !r  by subtracting the uniform background total potential is finite.
  !r  In the case r -> 0, The onsite term h(r) is subracted making the
  !r  potential finite at the origin.
  !r  h^s is calculated in Fourier space as
  !r  (3) \sum_R h^s(r-R) =
  !r     4\pi / \omega \sum_G G^{-2} \exp{-G^2/4 a^2 + i G \dot r}
  !r  The term G=0 is infinite in consequence of the infinite range of
  !r  the 1/r potential, but is canceled by a uniform background, up to
  !r  a constant of integration.  The constant is chosen to make Phi
  !r  independent of a.  Equation (2) becomes:
  !r  (4) Phi(r) = \sum_R h^d (r-R) +
  !r               \sum_{G \ne 0} 4\pi / \omega G^{-2}
  !r                              \exp{-G^2/4 a^2 + i G \dot r}
  !r  and in the special case r->0 it is noted that
  !r  (5) \sum_R h^d (r-R) - 1/r = \sum_{R \ne 0} h^d (r-R) - h^s(r),
  !r  and the contribution h^s(0) = 2*a/\sqrt(\pi) is subtracted.
  !r
  !r  First lattice vector must be the zero vector.
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: nkg,nkd
  double precision :: tau(3),glat(3,nkg),dlat(3,nkd)
  double precision :: awald,alat,vol,dl
  ! local parameters
  integer :: ir,ir1
  double precision :: pi,tpi,gamma,tpiba2,r1,r2,d1mach

  pi = 4*datan(1d0)
  tpi = 2*pi
  gamma = 0.25d0/awald**2
  tpiba2 = (tpi/alat)**2
  dl = -gamma
  ! --- Reciprocal space sums, backwards to sum small numbers first ---
  do  26  ir = nkg, 2, -1
     r2 = tpiba2*(glat(1,ir)**2 + glat(2,ir)**2 + glat(3,ir)**2)
     dl = dl + dexp(-gamma*r2)/r2*dcos(tpi* &
          (glat(1,ir)*tau(1)+ glat(2,ir)*tau(2)+ glat(3,ir)*tau(3)))
26 enddo
  dl = dl*4*pi/vol

  ! --- Real space sums, backwards to sum small numbers first  ---
  if (tau(1)**2 + tau(2)**2 + tau(3)**2 > d1mach(3)) then
     ir1 = 1
  else
     ir1 = 2
     dl = dl - 2d0*awald/dsqrt(pi)
  endif
  do  20  ir = nkd, ir1, -1
     r1 = alat*dsqrt((tau(1)-dlat(1,ir))**2 + &
          (tau(2)-dlat(2,ir))**2 + &
          (tau(3)-dlat(3,ir))**2)
     dl = dl + erfc(awald*r1)/r1
20 enddo
end subroutine strx00
endmodule m_madmat
