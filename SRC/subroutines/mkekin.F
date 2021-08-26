      subroutine mkekin(sv_p_osig,sv_p_otau,sv_p_oppi,sv_p_oqkkl,vconst,smpot,smrho,sumev, sumtv)
      use m_struc_def  
      use m_lmfinit,only:lso,nkaph,nsp,nspc,stdo,nbas,ssite=>v_ssite,sspec=>v_sspec,ham_ldham,iprmb !lcplxp,
      use m_lattic,only: lat_vol
      use m_supot,only: lat_nabc,k1,k2,k3
C- Evaluate the valence kinetic energy
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ndimh :dimension of hamiltonian (ham_ldham)
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa kmxt lmxb
Ci     Stored:    *
Ci     Passed to: *
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: nabc vol
Ci     Stored:    *
Ci     Passed to: *
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: oindxo
Ci     Stored:    *
Ci     Passed to: lgors
Ci   lcplxp=1 only now: 0 if ppi is real; 1 if ppi is complex
Ci   osig  :augmentation overlap integrals
Ci   otau  :augmentation kinetic energy integrals
Ci   oppi  :augmentation kinetic + potential integrals
Ci   oqkkl :local density-matrix; see rlocbl.f
Ci   k1..3 :dimensions smpot,smrho
Ci   vconst:constant potential added to hamiltonian
Ci   smpot :smooth input potential on uniform mesh (mkpot.f)
Ci         :smpot = Ves~ + vxc = Ves(n0 + compensating gaussians) + vxc
Ci   smrho :smooth output density n0 (no gaussians n0~-n0)
Ci   sumev :sum of eigenvalues
Co Outputs
Co   sumtv :kinetic energy
Cl Local variables
Cl   sraugm:sum_ib q * (tau+ppi-tau) : corresponds to valftr in mkpot.f
Cl   smresh:sm rho * sm V ; corresponds to valfsm in mkpot.f
Cl   lso   :1 include L.S coupling; 2 include LzSz part only
Cr Remarks
Cr   The valence kinetic energy is evaluated in the usual way as
Cr        sumtv = sumev - srhov
Cr   where sumev is the sum-of-eigenvalues and srhov is the integral
Cr   of the output density and input potential.
Cr   Integrals of the density with the xc potential are folded into the
Cr   electrostatic parts, that is:
Cr     V0 = V0es + V0xc  V1 = V1es + V1xc  V2 = V2es + V2xc
Cr   and are not discussed here.
Cr
Cr   mkekin make the electrostatic integral
Cr     int (n0~ Ves~ + n1 Ves1 - n2 Ves2~)                         (40)
Cr   as described in
Cr      M. Methfessel, M. van Schilfgaarde, and R. A. Casali,
Cr      Lecture Notes in Physics, {\bf 535}. H. Dreysse,
Cr      ed. (Springer-Verlag, Berlin) 2000.
Cr   for an output density defined through (smrho,qkkl) and an input
Cr   potential defined through smpot and the matrix elements (ppi-tau).
C
Cr   Consider only one atom/cell for simplicity. The output density is
Cr     n = sum_ij Dij F~i F~j with Dij = {sum_n w_n z*_in z_jn}   (32)
Cr   i.e. the contraction of the density matrix with partial densities
Cr     F~i F~j = Fi Fj +
Cr               sum_kLk'L' C^i_kL (P~kLP~k'L' - PkLPk'L') C^j_k'L' (33)
Cr             = n0_ij + n1_ij - n2_ij
Cr   Note that local parts of the partial densities have two `levels'
Cr   of decomposition, namely at the `ij' level as in Eq. 33, or
Cr   at a still finer level in which the (kLk'L') indices are not
Cr   summed over.  Thus
Cr     n{1,2} = sum_ij D_ij n{1,2}_ij
Cr     n{1,2}_ij = sum_kLk'L' C^i_kL n{1,2}_kL,k'L' C^j_k'L'
Cr     n{1,2}_kL,k'L' = PkL Pk'L'
Cr   Note also that the 'k index' may be a sum over polynomials, or when
Cr   function `heads' are dealt with, the function itself, as described
Cr   in augmat.f.  As in making the matrix elements, we have to deal
Cr   with three cases, HH; HP; PP, but this is a inessential detail
Cr   needed only because representing H with sums of polynomials tends
Cr   to be ill-conditioned, and in the description below we ignore it.
Cr
Cr   Densities n0 and n2 have corresponding n0~ and n2~ which include
Cr   the additional multipole terms that guarantee n1 and n2~ have
Cr   the same multipole moments.  Thus:
Cr     n0~ = n0 + sum_M q_M G_M
Cr     n2~ = n2 + sum_M q_M G_M
Cr   where q_M are the difference in multipole moments between n1 and n2
Cr     q_M = int dr Y_M r^m (n1 - n2)
Cr   We can define partial densities for multipole contributions as well
Cr     n2~-n2 = sum_ij D_ij (n2~-n2)_ij
Cr     (n2~-n2)_ij = sum_M Q_ijM G_M
Cr                 = sum_kLk'L'M C^i_kL Q_kkLL'M G_M C^j_k'L'
Cr   with the two forms decomposing q_M into two levels:
Cr     q_M = sum_ij D_ij Q_ijM
Cr     Q_ijM = sum_kLk'L' C^i_kL Q_kkLL'M C^j_k'L'
Cr     Q_kkLL'M = int dr Y_M r^m (P~kL P~k'L' - PkL Pk'L')         (27)
Cr
Cr   Using the identity
Cr     n2~ - n2 = n0~ - n0 = sum_M q_M G_M
Cr   Eq. 40 is evaluated as
Cr     int n0 Ves0~ + n1 Ves1 - n2 Ves2~ + sum_M q_M G_M (Ves0~-Ves2~)
Cr   The first term is evaluated on the mesh and stored in srmesh
Cr   The remaining terms amount to products of the density-matrix
Cr   and the ppi matrix elements.  Thus:
Cr     int n1 Ves1 = sum_ij D_ij int n1_ij Ves1
Cr     int n1_ij Ves1 = sum_kLk'L' C^i_kL int n1_kL,k'L' Ves1 C^j_k'L'
Cr                    = sum_kLk'L' C^i_kL pi1_kk'LL' C^j_k'L'
Cr   where pi1 is the first term of the pi matrix element, Eq. 29:
Cr     pi1_kk'LL' = P~kL V1 P~k'L'
Cr   Similarly for the second term, substituting n2 for n1 and
Cr   Ves2~ for Ves1.
Cr     int n2 Ves2~ = sum_ij D_ij int n2_ij Ves2~
Cr     int n2_ij Ves2~ = sum_kLk'L' C^i_kL int n2_kL,k'L' Ves2~ C^j_k'L'
Cr                    = sum_kLk'L' C^i_kL pi2_kk'LL' C^j_k'L'
Cr     pi2_kk'LL' = P~kL V1 P~k'L'
Cr   The last term just amounts to products of the density-matrix and
Cr   the remaining parts of the ppi matrix element:
Cr     pi_kk'LL'  = pi1_kk'LL' - pi2_kk'LL' + pi3_kk'LL'
Cr     pi3_kk'LL' = sum_M Q_kkLL'M int G_M (Ves0~ - Ves2~)
Cr   Evaluating the last term in the electrostatic integral we have
Cr     rhoV_MP = int sum_M q_M G_M (Ves0~ - Ves2~)
Cr             = int sum_ij D_ij sum_M Q_ijM G_M (Ves0~ - Ves2~)
Cr             = sum_ij D_ij sum_kLk'L'M C^i_kL pi3_kk'LL' C^j_k'L'
Cr   which follows using the relationship between Q_kkLL'M and Q_ijM
Cr   Using the definition of the local density-matrix (see rlocbl.f)
Cr      qpp_kLk'L' = sum_ij D_ij C^i_kL C^j_k'L'
Cr   the electrostatic integral then becomes
Cr     int rhoVes = int n0 Ves0~ + n1 Ves1 - n2 Ves2~ + rhoV_MP
Cr                = int n0 Ves0~
Cr                + sum_ij D_ij sum_kLk'L' C^i_kL pi_kk'LL' C^j_k'L'
Cr                = int n0 Ves0~ + sum_kLk'L' qpp'LL' pi_kk'LL'
Cr
Cu Updates
Cu   29 Jun 05 (MvS) SO hamiltonian not included in d.c. terms when
Cu             evaluating kinetic energy.
Cu   27 Aug 01 Extended to local orbitals.
Cu   18 Jun 00 spin polarized
Cu   20 Jun 00 adapted from nfp get_ekin
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ndimh,i_copy_size
      type(s_rv1),target :: sv_p_oppi(3,nbas)
      type(s_rv1),target :: sv_p_otau(3,1)
      type(s_rv1),target :: sv_p_osig(3,1)
      type(s_rv1),target :: sv_p_oqkkl(3,nbas)
      real(8),pointer:: OQPP(:), OQHP(:),OQHH(:),
     &       OSIGHH(:),OSIGHP(:),OSIGPP(:),
     &       OTAUHH(:),OTAUHP(:),OTAUPP(:),
     &       OPPIHH(:),OPPIHP(:),OPPIPP(:)
      real(8):: sumev , sumtv , vconst
      double complex smpot(k1,k2,k3,2),smrho(k1,k2,k3,2)
C ... Local parameters
      integer ib,igetss,ipr,is,kmax,lgunit,lmxa,lmxh,n0,n1,n2,n3,
     .ngabc(3),nglob,nkap0,nlma,nlmh,norb
c lso,isw
      logical lgors
      parameter (n0=10,nkap0=3)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
     .ntab(n0*nkap0),blks(n0*nkap0)
      double precision qum1,qum2,sraugm,srhov,srmesh,sum1,sum2,sumh,
     .sumq,sumt,vol,xx
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
#if (MPI|MPIK)
       include 'mpif.h'
       integer:: procid=0,ier=0
       integer,parameter::master=0
       logical:: iprx
       call mpi_comm_rank(mpi_comm_world,procid,ier)
       iprx=.false.
       if(procid==master) iprx=.true.
#endif
      ndimh=ham_ldham(1) 
      call getpr(ipr)
      ngabc=lat_nabc
      vol=lat_vol

C --- Integral n0(out) (Ves0~ + Vxc0), contribution from mesh ---
C     Note that it does not include the term (n0~-n0) Ves0~
c      print *,' mkekin: procid sumsmpot sumsmrho sum1 sum2=',procid,sum(abs(smpot)),sum(abs(smrho)),sum1,sum2
      call mshdot(vol,nsp,n1,n2,n3,k1,k2,k3,smpot,smrho,sum1,sum2)
      call mshint(vol,nsp,n1,n2,n3,k1,k2,k3,smrho,qum1,qum2)
      srmesh = sum1 + vconst*qum1
C --- Integral rhout*veff, part from augmentation ---
      sraugm = 0d0
      do  ib = 1, nbas
        is = int(ssite(ib)%spec)
        lmxa=sspec(is)%lmxa
        kmax=sspec(is)%kmxt
        lmxh=sspec(is)%lmxb
        if (lmxa .eq. -1) goto 10
        call orbl ( ib , 0 , ndimh , iprmb , norb , ltab , ktab 
     .   , xx , offl , xx )
C       Block into groups of consecutive l
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
        OPPIHH =>sv_p_oppi(3,ib)%v
        OPPIHP =>sv_p_oppi(2,ib)%v
        OPPIPP =>sv_p_oppi(1,ib)%v
        call pvgtkn ( kmax , lmxa , nlma , nkaph , norb , ltab , ktab
     .              , blks , lmxh , nlmh , OTAUHH , OSIGHH , OPPIHH
     .              , OTAUHP , OSIGHP , OPPIHP
     .              , OTAUPP , OSIGPP , OPPIPP !, lcplxp
     .              , lso , OQHH , OQHP , OQPP , nsp , nspc , sumt 
     .              , sumq , sumh )
C       Add site augmentation contribution to rhout * (ham - ke)
        sraugm = sraugm + sumh - sumt
   10   continue
      enddo
      srhov = srmesh + sraugm
      sumtv = sumev - srhov
      if (ipr .ge. 30) write(stdo,340) srmesh,sraugm,srhov,sumev,sumtv
  340 format(/' srhov:',3f14.6,' sumev=',f12.6,'   sumtv=',f12.6)
      end subroutine mkekin


      subroutine pvgtkn(kmax,lmxa,nlma,nkaph,norb,ltab,ktab,blks,lmxh,
     .nlmh,tauhh,sighh,ppihhz,tauhp,sighp,ppihpz,
     .taupp,sigpp,ppippz,lso,qhh,qhp,qpp,nsp,nspc, !lcplxp,
     .sumt,sumq,sumh)

C- Local contribution to kinetic energy for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   kmax  :cutoff in PkL expansion
Ci   lmxa  :dimensions sigpp, taupp
Ci   nlma  :L cutoff in PkL expansion
Ci   nkaph :dimensions augmentation matrices
Ci   norb  :number of orbitals for this site
Ci   ltab  :table of l quantum numbers for the orbitals
Ci   ktab  :table of k numbers (orbital type) for the orbitals
Ci   blks  :block size for grouping orbitals into blocks (gtbls1)
Ci   lmxh  :dimensions sighh, sighp, tauhh, tauhp
Ci   nlmh  :dimensions heads ppi and qhh and qhp
Ci   tauhh :head-head kinetic energy integrals (augmat.f)
Ci   sighh :head-head overlap integrals (augmat.f)
Ci   ppihh :head-head kinetic + potential integrals (augmat.f)
Ci   tauhp :head-tail kinetic energy integrals (augmat.f)
Ci   sighp :head-tail overlap integrals (augmat.f)
Ci   ppihp :head-tail kinetic + potential integrals (augmat.f)
Ci   taupp :tail-tail kinetic energy integrals (augmat.f)
Ci   sigpp :tail-tail overlap integrals (augmat.f)
Ci   ppipp :tail-tail potential integrals (augmat.f)
Ci   lcplxp=1 only:0 if ppi is real; 1 if ppi is complex
Ci   lso   :1 include L.S coupling; 2 include LzSz part only
Ci   qhh   :head-head density matrix for this site
Ci   qhp   :head-tail density matrix for this site
Ci   qpp   :tail-tail density matrix for this site
Ci   nsp   :number of spin channels
Ci   nspc  :2 for coupled spins; otherwise 1
Co Outputs
Co   sumt  :site contribution to kinetic energy
Co   sumq  :site contribution to overlap (charge ?)
Co   sumh  :site contribution to kinetic energy + potential
Cr Remarks
Cu Updates
Cu    1 Sep 04 Adapted to handle complex ppi
Cu   28 Aug 01 Extended to local orbitals.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer kmax,lmxa,nlma,lmxh,nlmh,nsp,nspc,lso !lcplxp,
      integer nkaph,norb,ltab(norb),ktab(norb),blks(norb)
      double precision
     .tauhh(nkaph,nkaph,0:lmxh,1),sighh(nkaph,nkaph,0:lmxh,1),
     .tauhp(nkaph,0:kmax,0:lmxh,1),sighp(nkaph,0:kmax,0:lmxh,1),
     .taupp(0:kmax,0:kmax,0:lmxa,1),sigpp(0:kmax,0:kmax,0:lmxa,1),
     .qhh(nkaph,nkaph,nlmh,nlmh,1),
     .qhp(nkaph,0:kmax,nlmh,nlma,1),
     .qpp(0:kmax,0:kmax,nlma,nlma,1)
      double complex
     .ppihhz(nkaph,nkaph,nlmh,nlmh,1),
     .ppihpz(nkaph,0:kmax,nlmh,nlma,1),
     .ppippz(0:kmax,0:kmax,nlma,nlma,1)
C ... Local parameters
      integer ilm1,ilm2,k1,k2,ll,nlm11,nlm12,nlm21,nlm22,i
      double precision sumt,sumq,sumh,xx
      integer io1,io2,l1,l2

!   ... Remove SO part from potential.
! 26gau2021 takao: I don't exactly know why I can skip this.
! From the beginnig definition of ppi only includes z-axis contribution (25aug2021).
! Since ppi is the pi integral plus LzSz part), do we need to remove LzSz contribution from ppihhz?
! But, anyway, it seems that skiping 'call ppi2z2' passs old tests'. I don't know why...      
c    
c      if ( lso .ne. 0) then
c        call ppi2z2(6,nsp,nspc,nkaph,nkaph,nlmh,nlmh,  ppihhz)
c        call ppi2z2(6,nsp,nspc,nkaph,1+kmax,nlmh,nlma, ppihpz)
c        call ppi2z2(6,nsp,nspc,1+kmax,1+kmax,nlma,nlma,ppippz)
c      endif
ccccccccccccccccccccccccccccccccccccccc
      
      sumt = 0d0
      sumq = 0d0
      sumh = 0d0

C ... Pkl*Pkl
      do  i = 1, nsp
        do  k1 = 0, kmax
          do  k2 = 0, kmax
c$$$            if (lcplxp .eq. 0) then
c$$$              do  ilm1 = 1, nlma
c$$$                l1 = ll(ilm1)
c$$$                sumt = sumt + qpp(k1,k2,ilm1,ilm1,i)*taupp(k1,k2,l1,i)
c$$$                sumq = sumq + qpp(k1,k2,ilm1,ilm1,i)*sigpp(k1,k2,l1,i)
c$$$                do  ilm2 = 1, nlma
c$$$                  sumh = sumh + qpp(k1,k2,ilm1,ilm2,i)*ppipp(k1,k2,ilm1,ilm2,i)
c$$$                enddo
c$$$              enddo
c$$$
c$$$            else
              do  ilm1 = 1, nlma
                l1 = ll(ilm1)
                sumt = sumt + qpp(k1,k2,ilm1,ilm1,i)*taupp(k1,k2,l1,i)
                sumq = sumq + qpp(k1,k2,ilm1,ilm1,i)*sigpp(k1,k2,l1,i)
                do  ilm2 = 1, nlma
                  sumh = sumh + qpp(k1,k2,ilm1,ilm2,i)*ppippz(k1,k2,ilm1,ilm2,i)
                enddo
              enddo
c$$$            endif
          enddo
        enddo

C ... Hsm*Hsm
        do  io2 = 1, norb
          if (blks(io2) .ne. 0) then
C       k2,l2 = k and starting l index for this block
            l2 = ltab(io2)
            k2 = ktab(io2)
            nlm21 = l2**2+1
            nlm22 = nlm21 + blks(io2)-1
            do  io1 = 1, norb
              if (blks(io1) .ne. 0) then
C         k1,l1 = k and starting l index for this block
                l1 = ltab(io1)
                k1 = ktab(io1)
                nlm11 = l1**2+1
                nlm12 = nlm11 + blks(io1)-1
c$$$                if (lcplxp .eq. 0) then
c$$$                  do  ilm1 = nlm11, nlm12
c$$$                    l1 = ll(ilm1)
c$$$                    do  ilm2 = nlm21, nlm22
c$$$                      if (ilm1 .eq. ilm2) then
c$$$                        sumt = sumt + qhh(k1,k2,ilm1,ilm2,i)*tauhh(k1,k2,l1,i)
c$$$                        sumq = sumq + qhh(k1,k2,ilm1,ilm2,i)*sighh(k1,k2,l1,i)
c$$$                      endif
c$$$                      sumh = sumh+ qhh(k1,k2,ilm1,ilm2,i)*ppihh(k1,k2,ilm1,ilm2,i)
c$$$                    enddo
c$$$                  enddo
c$$$
c$$$                else
                  do  ilm1 = nlm11, nlm12
                    l1 = ll(ilm1)
                    do  ilm2 = nlm21, nlm22
                      if (ilm1 .eq. ilm2) then
                        sumt = sumt + qhh(k1,k2,ilm1,ilm2,i)*tauhh(k1,k2,l1,i)
                        sumq = sumq + qhh(k1,k2,ilm1,ilm2,i)*sighh(k1,k2,l1,i)
                      endif
                      sumh = sumh+qhh(k1,k2,ilm1,ilm2,i)*ppihhz(k1,k2,ilm1,ilm2,i)
                    enddo
                  enddo
c$$$                endif
C          Needs to be checked ...
C          if (nspc .eq. 2) then
C            call pvgtk1(nkaph,nlmh,k1,k2,nlm11,nlm12,nlm21,nlm22,
C     .        qhh,ppihhz,sumh)
C          endif

              endif
            enddo

          endif
        enddo

C ... Hsm*Pkl
        do  io1 = 1, norb
          if (blks(io1) .ne. 0) then
C       k1,l1 = k and starting l index for this block
            l1 = ltab(io1)
            k1 = ktab(io1)
            nlm11 = l1**2+1
            nlm12 = nlm11 + blks(io1)-1
            do  k2 = 0, kmax

c$$$              if (lcplxp .eq. 0) then
c$$$                do  ilm1 = nlm11, nlm12
c$$$                  l1 = ll(ilm1)
c$$$                  do  ilm2 = 1, nlma
c$$$                    if (ilm1 .eq. ilm2) then
c$$$                      sumt = sumt + qhp(k1,k2,ilm1,ilm2,i)*tauhp(k1,k2,l1,i)
c$$$                      sumq = sumq + qhp(k1,k2,ilm1,ilm2,i)*sighp(k1,k2,l1,i)
c$$$                    endif
c$$$                    sumh = sumh+qhp(k1,k2,ilm1,ilm2,i)*ppihp(k1,k2,ilm1,ilm2,i)
c$$$                  enddo
c$$$                enddo
c$$$              else
                do  ilm1 = nlm11, nlm12
                  l1 = ll(ilm1)
                  do  ilm2 = 1, nlma
                    if (ilm1 .eq. ilm2) then
                      sumt = sumt + qhp(k1,k2,ilm1,ilm2,i)*tauhp(k1,k2,l1,i)
                      sumq = sumq + qhp(k1,k2,ilm1,ilm2,i)*sighp(k1,k2,l1,i)
                    endif
                    sumh = sumh+qhp(k1,k2,ilm1,ilm2,i)*ppihpz(k1,k2,ilm1,ilm2,i)
                  enddo
                enddo
c$$$              endif

            enddo
          endif
        enddo
      enddo
c      if (.false.) then !lso==1) then ! .ne. 0) then
c        call ppi2z2(5,nsp,nspc,nkaph,nkaph,nlmh,nlmh,ppihhz)
c        call ppi2z2(5,nsp,nspc,nkaph,1+kmax,nlmh,nlma,ppihpz)
c        call ppi2z2(5,nsp,nspc,1+kmax,1+kmax,nlma,nlma,ppippz)
c      endif
      end subroutine pvgtkn

c$$$      subroutine ppi2z2(job,nsp,nspc,nf1,nf2,nlx1,nlx2,ppiz)
c$$$      integer job,nsp,nspc,nf1,nf2,nlx1,nlx2
c$$$      double complex  ppiz(nf1,nf2,nlx1,nlx2,nsp,nspc,2)
c$$$      integer i1,i
c$$$      real(8):: sz
c$$$      if (job==5) sz =  1d0
c$$$      if (job==6) sz = -1d0
c$$$      call daxpy(2*nf1*nf2*nlx1*nlx2*nsp*nspc,sz,ppiz(1,1,1,1,1,1,2),1,ppiz,1)
c$$$      end
