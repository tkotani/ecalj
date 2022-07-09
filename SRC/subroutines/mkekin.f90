subroutine mkekin(sv_p_osig,sv_p_otau,sv_p_oppi,sv_p_oqkkl,vconst,smpot,smrho,sumev, sumtv)
  use m_struc_def
  use m_lmfinit,only:lso,nkaph,nsp,nspc,stdo,nbas,ispec,sspec=>v_sspec,nlmto
  use m_lattic,only: lat_vol
  use m_supot,only: lat_nabc,k1,k2,k3
  use m_orbl,only: Orblib,ktab,ltab,offl,norb

  !- Evaluate the valence kinetic energy
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   nlmto :dimension of hamiltonian (MTO)
  !i   lcplxp=1 only now: 0 if ppi is real; 1 if ppi is complex
  !i   osig  :augmentation overlap integrals
  !i   otau  :augmentation kinetic energy integrals
  !i   oppi  :augmentation kinetic + potential integrals
  !i   oqkkl :local density-matrix; see rlocbl.f
  !i   k1..3 :dimensions smpot,smrho
  !i   vconst:constant potential added to hamiltonian
  !i   smpot :smooth input potential on uniform mesh (mkpot.f)
  !i         :smpot = Ves~ + vxc = Ves(n0 + compensating gaussians) + vxc
  !i   smrho :smooth output density n0 (no gaussians n0~-n0)
  !i   sumev :sum of eigenvalues
  !o Outputs
  !o   sumtv :kinetic energy
  !l Local variables
  !l   sraugm:sum_ib q * (tau+ppi-tau) : corresponds to valftr in mkpot.f
  !l   smresh:sm rho * sm V ; corresponds to valfsm in mkpot.f
  !l   lso   :1 include L.S coupling; 2 include LzSz part only
  !r Remarks
  !r   The valence kinetic energy is evaluated in the usual way as
  !r        sumtv = sumev - srhov
  !r   where sumev is the sum-of-eigenvalues and srhov is the integral
  !r   of the output density and input potential.
  !r   Integrals of the density with the xc potential are folded into the
  !r   electrostatic parts, that is:
  !r     V0 = V0es + V0xc  V1 = V1es + V1xc  V2 = V2es + V2xc
  !r   and are not discussed here.
  !r
  !r   mkekin make the electrostatic integral
  !r     int (n0~ Ves~ + n1 Ves1 - n2 Ves2~)                         (40)
  !r   as described in
  !r      M. Methfessel, M. van Schilfgaarde, and R. A. Casali,
  !r      Lecture Notes in Physics, {\bf 535}. H. Dreysse,
  !r      ed. (Springer-Verlag, Berlin) 2000.
  !r   for an output density defined through (smrho,qkkl) and an input
  !r   potential defined through smpot and the matrix elements (ppi-tau).

  !r   Consider only one atom/cell for simplicity. The output density is
  !r     n = sum_ij Dij F~i F~j with Dij = {sum_n w_n z*_in z_jn}   (32)
  !r   i.e. the contraction of the density matrix with partial densities
  !r     F~i F~j = Fi Fj +
  !r               sum_kLk'L' C^i_kL (P~kLP~k'L' - PkLPk'L') C^j_k'L' (33)
  !r             = n0_ij + n1_ij - n2_ij
  !r   Note that local parts of the partial densities have two `levels'
  !r   of decomposition, namely at the `ij' level as in Eq. 33, or
  !r   at a still finer level in which the (kLk'L') indices are not
  !r   summed over.  Thus
  !r     n{1,2} = sum_ij D_ij n{1,2}_ij
  !r     n{1,2}_ij = sum_kLk'L' C^i_kL n{1,2}_kL,k'L' C^j_k'L'
  !r     n{1,2}_kL,k'L' = PkL Pk'L'
  !r   Note also that the 'k index' may be a sum over polynomials, or when
  !r   function `heads' are dealt with, the function itself, as described
  !r   in augmat.f.  As in making the matrix elements, we have to deal
  !r   with three cases, HH; HP; PP, but this is a inessential detail
  !r   needed only because representing H with sums of polynomials tends
  !r   to be ill-conditioned, and in the description below we ignore it.
  !r
  !r   Densities n0 and n2 have corresponding n0~ and n2~ which include
  !r   the additional multipole terms that guarantee n1 and n2~ have
  !r   the same multipole moments.  Thus:
  !r     n0~ = n0 + sum_M q_M G_M
  !r     n2~ = n2 + sum_M q_M G_M
  !r   where q_M are the difference in multipole moments between n1 and n2
  !r     q_M = int dr Y_M r^m (n1 - n2)
  !r   We can define partial densities for multipole contributions as well
  !r     n2~-n2 = sum_ij D_ij (n2~-n2)_ij
  !r     (n2~-n2)_ij = sum_M Q_ijM G_M
  !r                 = sum_kLk'L'M C^i_kL Q_kkLL'M G_M C^j_k'L'
  !r   with the two forms decomposing q_M into two levels:
  !r     q_M = sum_ij D_ij Q_ijM
  !r     Q_ijM = sum_kLk'L' C^i_kL Q_kkLL'M C^j_k'L'
  !r     Q_kkLL'M = int dr Y_M r^m (P~kL P~k'L' - PkL Pk'L')         (27)
  !r
  !r   Using the identity
  !r     n2~ - n2 = n0~ - n0 = sum_M q_M G_M
  !r   Eq. 40 is evaluated as
  !r     int n0 Ves0~ + n1 Ves1 - n2 Ves2~ + sum_M q_M G_M (Ves0~-Ves2~)
  !r   The first term is evaluated on the mesh and stored in srmesh
  !r   The remaining terms amount to products of the density-matrix
  !r   and the ppi matrix elements.  Thus:
  !r     int n1 Ves1 = sum_ij D_ij int n1_ij Ves1
  !r     int n1_ij Ves1 = sum_kLk'L' C^i_kL int n1_kL,k'L' Ves1 C^j_k'L'
  !r                    = sum_kLk'L' C^i_kL pi1_kk'LL' C^j_k'L'
  !r   where pi1 is the first term of the pi matrix element, Eq. 29:
  !r     pi1_kk'LL' = P~kL V1 P~k'L'
  !r   Similarly for the second term, substituting n2 for n1 and
  !r   Ves2~ for Ves1.
  !r     int n2 Ves2~ = sum_ij D_ij int n2_ij Ves2~
  !r     int n2_ij Ves2~ = sum_kLk'L' C^i_kL int n2_kL,k'L' Ves2~ C^j_k'L'
  !r                    = sum_kLk'L' C^i_kL pi2_kk'LL' C^j_k'L'
  !r     pi2_kk'LL' = P~kL V1 P~k'L'
  !r   The last term just amounts to products of the density-matrix and
  !r   the remaining parts of the ppi matrix element:
  !r     pi_kk'LL'  = pi1_kk'LL' - pi2_kk'LL' + pi3_kk'LL'
  !r     pi3_kk'LL' = sum_M Q_kkLL'M int G_M (Ves0~ - Ves2~)
  !r   Evaluating the last term in the electrostatic integral we have
  !r     rhoV_MP = int sum_M q_M G_M (Ves0~ - Ves2~)
  !r             = int sum_ij D_ij sum_M Q_ijM G_M (Ves0~ - Ves2~)
  !r             = sum_ij D_ij sum_kLk'L'M C^i_kL pi3_kk'LL' C^j_k'L'
  !r   which follows using the relationship between Q_kkLL'M and Q_ijM
  !r   Using the definition of the local density-matrix (see rlocbl.f)
  !r      qpp_kLk'L' = sum_ij D_ij C^i_kL C^j_k'L'
  !r   the electrostatic integral then becomes
  !r     int rhoVes = int n0 Ves0~ + n1 Ves1 - n2 Ves2~ + rhoV_MP
  !r                = int n0 Ves0~
  !r                + sum_ij D_ij sum_kLk'L' C^i_kL pi_kk'LL' C^j_k'L'
  !r                = int n0 Ves0~ + sum_kLk'L' qpp'LL' pi_kk'LL'
  !r
  !u Updates
  !u   29 Jun 05 (MvS) SO hamiltonian not included in d.c. terms when
  !u             evaluating kinetic energy.
  !u   27 Aug 01 Extended to local orbitals.
  !u   18 Jun 00 spin polarized
  !u   20 Jun 00 adapted from nfp get_ekin
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: i_copy_size
  type(s_cv1),target :: sv_p_oppi(3,nbas)
  type(s_rv1),target :: sv_p_otau(3,1)
  type(s_rv1),target :: sv_p_osig(3,1)
  type(s_rv1),target :: sv_p_oqkkl(3,nbas)
  real(8),pointer:: OQPP(:), OQHP(:),OQHH(:), &
       OSIGHH(:),OSIGHP(:),OSIGPP(:), &
       OTAUHH(:),OTAUHP(:),OTAUPP(:)
  complex(8),pointer:: OPPIHH(:),OPPIHP(:),OPPIPP(:)
  real(8):: sumev , sumtv , vconst
  double complex smpot(k1,k2,k3,2),smrho(k1,k2,k3,2)
  ! ... Local parameters
  integer :: ib,igetss,ipr,is,kmax,lgunit,lmxa,lmxh,n0,n1,n2,n3, &
       ngabc(3),nglob,nkap0,nlma,nlmh!,norb
  ! lso,isw
  logical :: lgors
  parameter (n0=10,nkap0=3)
  !      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
  integer:: ntab(n0*nkap0),blks(n0*nkap0)
  double precision :: qum1,qum2,sraugm,srhov,srmesh,sum1,sum2,sumh, &
       sumq,sumt,vol,xx
  equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
  include 'mpif.h'
  integer:: procid=0,ier=0
  integer,parameter::master=0
  logical:: iprx
  call mpi_comm_rank(mpi_comm_world,procid,ier)
  iprx=.false.
  if(procid==master) iprx= .TRUE. 
  call tcn('mkekin')
  !      ndimh=ham_ldham(1)
  call getpr(ipr)
  ngabc=lat_nabc
  vol=lat_vol

  ! --- Integral n0(out) (Ves0~ + Vxc0), contribution from mesh ---
  !     Note that it does not include the term (n0~-n0) Ves0~
  !      print *,' mkekin: procid sumsmpot sumsmrho sum1 sum2=',procid,sum(abs(smpot)),sum(abs(smrho)),sum1,sum2
  call mshdot(vol,nsp,n1,n2,n3,k1,k2,k3,smpot,smrho,sum1,sum2)
  call mshint(vol,nsp,n1,n2,n3,k1,k2,k3,smrho,qum1,qum2)
  srmesh = sum1 + vconst*qum1
  ! --- Integral rhout*veff, part from augmentation ---
  sraugm = 0d0
  do  ib = 1, nbas
     is = ispec(ib) !int(ssite(ib)%spec)
     lmxa=sspec(is)%lmxa
     kmax=sspec(is)%kmxt
     lmxh=sspec(is)%lmxb
     if (lmxa == -1) goto 10
     call orblib(ib) ! norb , ltab , ktab , offl 
     !       Block into groups of consecutive l
     call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)
     nlma = (lmxa+1)**2
     nlmh = (lmxh+1)**2
     OQHH => sv_p_oqkkl(3,ib)%v
     OQHP => sv_p_oqkkl(2,ib)%v
     OQPP => sv_p_oqkkl(1,ib)%v
     OTAUHH =>sv_p_otau(3,ib)%v
     OTAUHP =>sv_p_otau(2,ib)%v
     OTAUPP =>sv_p_otau(1,ib)%v
     OSIGHH =>sv_p_osig(3,ib)%v
     OSIGHP =>sv_p_osig(2,ib)%v
     OSIGPP =>sv_p_osig(1,ib)%v
     OPPIHH =>sv_p_oppi(3,ib)%cv
     OPPIHP =>sv_p_oppi(2,ib)%cv
     OPPIPP =>sv_p_oppi(1,ib)%cv
     call pvgtkn ( kmax , lmxa , nlma , nkaph , norb , ltab , ktab &
          , blks , lmxh , nlmh , OTAUHH , OSIGHH , OPPIHH &
          , OTAUHP , OSIGHP , OPPIHP &
          , OTAUPP , OSIGPP , OPPIPP  &
          , lso , OQHH , OQHP , OQPP , nsp , nspc , sumt &
          , sumq , sumh )
     !       Add site augmentation contribution to rhout * (ham - ke)
     sraugm = sraugm + sumh - sumt
10   continue
  enddo
  srhov = srmesh + sraugm
  sumtv = sumev - srhov
  if (ipr >= 30) write(stdo,340) srmesh,sraugm,srhov,sumev,sumtv
340 format(/' srhov:',3f14.6,' sumev=',f12.6,'   sumtv=',f12.6)
  call tcx('mkekin')
end subroutine mkekin


subroutine pvgtkn(kmax,lmxa,nlma,nkaph,norb,ltab,ktab,blks,lmxh, &
     nlmh,tauhh,sighh,ppihhz,tauhp,sighp,ppihpz, &
     taupp,sigpp,ppippz,lso,qhh,qhp,qpp,nsp,nspc,&
     sumt,sumq,sumh)
  !- Local contribution to kinetic energy for one site
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   kmax  :cutoff in PkL expansion
  !i   lmxa  :dimensions sigpp, taupp
  !i   nlma  :L cutoff in PkL expansion
  !i   nkaph :dimensions augmentation matrices
  !i   norb  :number of orbitals for this site
  !i   ltab  :table of l quantum numbers for the orbitals
  !i   ktab  :table of k numbers (orbital type) for the orbitals
  !i   blks  :block size for grouping orbitals into blocks (gtbls1)
  !i   lmxh  :dimensions sighh, sighp, tauhh, tauhp
  !i   nlmh  :dimensions heads ppi and qhh and qhp
  !i   tauhh :head-head kinetic energy integrals (augmat.f)
  !i   sighh :head-head overlap integrals (augmat.f)
  !i   ppihh :head-head kinetic + potential integrals (augmat.f)
  !i   tauhp :head-tail kinetic energy integrals (augmat.f)
  !i   sighp :head-tail overlap integrals (augmat.f)
  !i   ppihp :head-tail kinetic + potential integrals (augmat.f)
  !i   taupp :tail-tail kinetic energy integrals (augmat.f)
  !i   sigpp :tail-tail overlap integrals (augmat.f)
  !i   ppipp :tail-tail potential integrals (augmat.f)
  !i   lcplxp=1 only:0 if ppi is real; 1 if ppi is complex
  !i   lso   :1 include L.S coupling; 2 include LzSz part only
  !i   qhh   :head-head density matrix for this site
  !i   qhp   :head-tail density matrix for this site
  !i   qpp   :tail-tail density matrix for this site
  !i   nsp   :number of spin channels
  !i   nspc  :2 for coupled spins; otherwise 1
  !o Outputs
  !o   sumt  :site contribution to kinetic energy
  !o   sumq  :site contribution to overlap (charge ?)
  !o   sumh  :site contribution to kinetic energy + potential
  !r Remarks
  !u Updates
  !u    1 Sep 04 Adapted to handle complex ppi
  !u   28 Aug 01 Extended to local orbitals.
  ! ----------------------------------------------------------------------
  implicit none
  integer :: kmax,lmxa,nlma,lmxh,nlmh,nsp,nspc,lso
  integer :: nkaph,norb,ltab(norb),ktab(norb),blks(norb)
  double precision :: &
       tauhh(nkaph,nkaph,0:lmxh,1),sighh(nkaph,nkaph,0:lmxh,1), &
       tauhp(nkaph,0:kmax,0:lmxh,1),sighp(nkaph,0:kmax,0:lmxh,1), &
       taupp(0:kmax,0:kmax,0:lmxa,1),sigpp(0:kmax,0:kmax,0:lmxa,1), &
       qhh(nkaph,nkaph,nlmh,nlmh,1), &
       qhp(nkaph,0:kmax,nlmh,nlma,1), &
       qpp(0:kmax,0:kmax,nlma,nlma,1)
  double complex &
       ppihhz(nkaph,nkaph,nlmh,nlmh,1), &
       ppihpz(nkaph,0:kmax,nlmh,nlma,1), &
       ppippz(0:kmax,0:kmax,nlma,nlma,1)
  ! ... Local parameters
  integer :: ilm1,ilm2,k1,k2,ll,nlm11,nlm12,nlm21,nlm22,i
  double precision :: sumt,sumq,sumh,xx
  integer :: io1,io2,l1,l2

  !   ... Remove SO part from potential.
  ! 26gau2021 takao: I don't exactly know why I can skip this.
  ! From the beginnig definition of ppi only includes z-axis contribution (25aug2021).
  ! Since ppi is the pi integral plus LzSz part), do we need to remove LzSz contribution from ppihhz?
  ! But, anyway, it seems that skiping 'call ppi2z2' passs old tests'. I don't know why...

  !      if ( lso .ne. 0) then
  !        call ppi2z2(6,nsp,nspc,nkaph,nkaph,nlmh,nlmh,  ppihhz)
  !        call ppi2z2(6,nsp,nspc,nkaph,1+kmax,nlmh,nlma, ppihpz)
  !        call ppi2z2(6,nsp,nspc,1+kmax,1+kmax,nlma,nlma,ppippz)
  !      endif
  ! ccccccccccccccccccccccccccccccccccccc

  sumt = 0d0
  sumq = 0d0
  sumh = 0d0
  ! ... Pkl*Pkl
  do  i = 1, nsp
     do  k1 = 0, kmax
        do  k2 = 0, kmax
           do  ilm1 = 1, nlma
              l1 = ll(ilm1)
              sumt = sumt + qpp(k1,k2,ilm1,ilm1,i)*taupp(k1,k2,l1,i)
              sumq = sumq + qpp(k1,k2,ilm1,ilm1,i)*sigpp(k1,k2,l1,i)
              do  ilm2 = 1, nlma
                 sumh = sumh + qpp(k1,k2,ilm1,ilm2,i)*ppippz(k1,k2,ilm1,ilm2,i)
              enddo
           enddo
        enddo
     enddo
     ! ... Hsm*Hsm
     do  io2 = 1, norb
        if (blks(io2) /= 0) then
           !       k2,l2 = k and starting l index for this block
           l2 = ltab(io2)
           k2 = ktab(io2)
           nlm21 = l2**2+1
           nlm22 = nlm21 + blks(io2)-1
           do  io1 = 1, norb
              if (blks(io1) /= 0) then
                 !         k1,l1 = k and starting l index for this block
                 l1 = ltab(io1)
                 k1 = ktab(io1)
                 nlm11 = l1**2+1
                 nlm12 = nlm11 + blks(io1)-1
                 do  ilm1 = nlm11, nlm12
                    l1 = ll(ilm1)
                    do  ilm2 = nlm21, nlm22
                       if (ilm1 == ilm2) then
                          sumt = sumt + qhh(k1,k2,ilm1,ilm2,i)*tauhh(k1,k2,l1,i)
                          sumq = sumq + qhh(k1,k2,ilm1,ilm2,i)*sighh(k1,k2,l1,i)
                       endif
                       sumh = sumh+qhh(k1,k2,ilm1,ilm2,i)*ppihhz(k1,k2,ilm1,ilm2,i)
                    enddo
                 enddo
              endif
           enddo
        endif
     enddo
     ! ... Hsm*Pkl
     do  io1 = 1, norb
        if (blks(io1) /= 0) then
           !       k1,l1 = k and starting l index for this block
           l1 = ltab(io1)
           k1 = ktab(io1)
           nlm11 = l1**2+1
           nlm12 = nlm11 + blks(io1)-1
           do  k2 = 0, kmax
              do  ilm1 = nlm11, nlm12
                 l1 = ll(ilm1)
                 do  ilm2 = 1, nlma
                    if (ilm1 == ilm2) then
                       sumt = sumt + qhp(k1,k2,ilm1,ilm2,i)*tauhp(k1,k2,l1,i)
                       sumq = sumq + qhp(k1,k2,ilm1,ilm2,i)*sighp(k1,k2,l1,i)
                    endif
                    sumh = sumh+qhp(k1,k2,ilm1,ilm2,i)*ppihpz(k1,k2,ilm1,ilm2,i)
                 enddo
              enddo
           enddo
        endif
     enddo
  enddo
end subroutine pvgtkn

