          subroutine mkekin ( nbas , ndimh , ssite , sspec , slat , sham 
     .    , sv_p_osig , sv_p_otau , sv_p_oppi , lcplxp , sv_p_oqkkl , k1 
     .    , k2 , k3 , vconst , smpot , smrho , sumev , sumtv ) 





      use m_struc_def  !Cgetarg
      use m_globalvariables
C- Evaluate the valence kinetic energy
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ndimh :dimension of hamiltonian
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
Ci   lcplxp:0 if ppi is real; 1 if ppi is complex
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
C     implicit none
C ... Passed parameters
      integer k1,k2,k3,nbas,ndimh,lcplxp
          type(s_rv1) :: sv_p_oppi(3,nbas)

         type(s_rv1) :: sv_p_otau(3,1)

        type(s_rv1) :: sv_p_osig(3,1)

       type(s_rv1) :: sv_p_oqkkl(3,nbas)

      real(8):: sumev , sumtv , vconst
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat
      type(s_ham)::sham

      double complex smpot(k1,k2,k3,2),smrho(k1,k2,k3,2)
C ... Local parameters
      integer ib,igetss,ipr,is,kmax,lgunit,lmxa,lmxh,n0,n1,n2,n3,
     .ngabc(3),nglob,nkap0,nlma,nlmh,norb,nsp,nspc,stdo,lso,isw
                   real(8),pointer :: rv_p_oppipp(:) =>NULL()

                  real(8),pointer :: rv_p_oppihp(:) =>NULL()

                 real(8),pointer :: rv_p_oppihh(:) =>NULL()

                real(8),pointer :: rv_p_otaupp(:) =>NULL()

               real(8),pointer :: rv_p_otauhp(:) =>NULL()

              real(8),pointer :: rv_p_otauhh(:) =>NULL()

             real(8),pointer :: rv_p_osigpp(:) =>NULL()

            real(8),pointer :: rv_p_osighp(:) =>NULL()

           real(8),pointer :: rv_p_osighh(:) =>NULL()

          real(8),pointer :: rv_p_oqpp(:) =>NULL()

         real(8),pointer :: rv_p_oqhp(:) =>NULL()

        real(8),pointer :: rv_p_oqhh(:) =>NULL()

ckino Dec.14.2011:         integer,pointer :: iv_p_oiprmb(:) =>NULL()

      logical lgors
      parameter (n0=10,nkap0=3)
      integer nkaph,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
     .ntab(n0*nkap0),blks(n0*nkap0)
      double precision qum1,qum2,sraugm,srhov,srmesh,sum1,sum2,sumh,
     .sumq,sumt,vol,xx
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
C ... Heap
#if DEBUG
      integer fopna
#endif

      stdo = lgunit(1)
Changenglob      nsp  = nglob('nsp')
      nsp  = globalvariables%nsp
Changenglob      nspc = nglob('nspc')
      nspc = globalvariables%nspc
Changenglob      nkaph = nglob('nkaph')
      nkaph = globalvariables%nkaph
      call getpr(ipr)

      i_copy_size=size(slat%nabc)
      call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)
      vol=slat%vol


ckino Dec.14.2011:         iv_p_oiprmb => sham%iv_p_oindxo 


      lso = isw ( iand(4,int(sham%lncol)) .ne.0 ) + 2 * isw ( iand(32,int(sham%lncol)) .ne.0 
     .) + 3 * isw ( iand(64,int(sham%lncol)) .ne.0 )


C     call zprm3('sm rho-out',0,smrho,k1,k2,k3*nsp)
#if DEBUG
      call info0(0,0,0,' rewinding pi file ...')
      is = fopna('ppi',-1,4)
      rewind is
#endif

C --- Integral n0(out) (Ves0~ + Vxc0), contribution from mesh ---
C     Note that it does not include the term (n0~-n0) Ves0~
      call mshdot(vol,nsp,n1,n2,n3,k1,k2,k3,smpot,smrho,sum1,sum2)
      call mshint(vol,nsp,n1,n2,n3,k1,k2,k3,smrho,qum1,qum2)
      srmesh = sum1 + vconst*qum1
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(6,"('vvv: mkekin :sum1 vconst*qum1 ',10d13.6)")
c     &  sum1,vconst*qum1
c      write(6,"(a,12d13.5)")'rrrrs: smchk smrho 1'
c     &  , sum(smrho(:,:,:,1)), sum(abs(smrho(:,:,:,1)))
c      write(6,"(a,12d13.5)")'rrrrs: smchk smrho 2'
c     &  , sum(smrho(:,:,:,2)), sum(abs(smrho(:,:,:,2)))
cccccccccccccccccccccccccccccccccccccccccccccccccccc

C --- Integral rhout*veff, part from augmentation ---
      sraugm = 0d0
      do  ib = 1, nbas
        is = int(ssite(ib)%spec)


        lmxa=sspec(is)%lmxa
        kmax=sspec(is)%kmxt
        lmxh=sspec(is)%lmxb

        if (lmxa .eq. -1) goto 10

ckino Dec.14.2011:           call orbl ( ib , 0 , ndimh , iv_p_oiprmb , norb , ltab , ktab 
ckino 12月.27.2011:           call orbl ( ib , 0 , ndimh , sham%iv_p_oindxo , norb , ltab , ktab 
         call orbl ( ib , 0 , ndimh , sham%iv_a_oindxo , norb , ltab , ktab 
     .   , xx , offl , xx ) 

C       Block into groups of consecutive l
        call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)

        nlma = (lmxa+1)**2
        nlmh = (lmxh+1)**2

        rv_p_oqhh => sv_p_oqkkl(3,ib)%v
        rv_p_oqhp => sv_p_oqkkl(2,ib)%v
        rv_p_oqpp => sv_p_oqkkl(1,ib)%v
        rv_p_otauhh => sv_p_otau(3,ib)%v
        rv_p_otauhp => sv_p_otau(2,ib)%v
        rv_p_otaupp => sv_p_otau(1,ib)%v
        rv_p_osighh => sv_p_osig(3,ib)%v
        rv_p_osighp => sv_p_osig(2,ib)%v
        rv_p_osigpp => sv_p_osig(1,ib)%v
        rv_p_oppihh => sv_p_oppi(3,ib)%v
        rv_p_oppihp => sv_p_oppi(2,ib)%v
        rv_p_oppipp => sv_p_oppi(1,ib)%v

                    call pvgtkn ( kmax , lmxa , nlma , nkaph , norb , ltab , ktab 
     .              , blks , lmxh , nlmh , rv_p_otauhh , rv_p_osighh , rv_p_oppihh 
     .              , rv_p_oppihh , rv_p_otauhp , rv_p_osighp , rv_p_oppihp , rv_p_oppihp 
     .              , rv_p_otaupp , rv_p_osigpp , rv_p_oppipp , rv_p_oppipp , lcplxp 
     .              , lso , rv_p_oqhh , rv_p_oqhp , rv_p_oqpp , nsp , nspc , sumt 
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
     .nlmh,tauhh,sighh,ppihh,ppihhz,tauhp,sighp,ppihp,ppihpz,
     .taupp,sigpp,ppipp,ppippz,lcplxp,lso,qhh,qhp,qpp,nsp,nspc,
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
Ci   lcplxp:0 if ppi is real; 1 if ppi is complex
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
C     implicit none
C ... Passed parameters
      integer kmax,lmxa,nlma,lmxh,nlmh,nsp,nspc,lcplxp,lso
      integer nkaph,norb,ltab(norb),ktab(norb),blks(norb)
      double precision
     .tauhh(nkaph,nkaph,0:lmxh,1),sighh(nkaph,nkaph,0:lmxh,1),
     .tauhp(nkaph,0:kmax,0:lmxh,1),sighp(nkaph,0:kmax,0:lmxh,1),
     .taupp(0:kmax,0:kmax,0:lmxa,1),sigpp(0:kmax,0:kmax,0:lmxa,1),
     .ppihh(nkaph,nkaph,nlmh,nlmh,1),qhh(nkaph,nkaph,nlmh,nlmh,1),
     .ppihp(nkaph,0:kmax,nlmh,nlma,1),qhp(nkaph,0:kmax,nlmh,nlma,1),
     .ppipp(0:kmax,0:kmax,nlma,nlma,1),qpp(0:kmax,0:kmax,nlma,nlma,1)
      double complex
     .ppihhz(nkaph,nkaph,nlmh,nlmh,1),
     .ppihpz(nkaph,0:kmax,nlmh,nlma,1),
     .ppippz(0:kmax,0:kmax,nlma,nlma,1)
C ... Local parameters
      integer ilm1,ilm2,k1,k2,ll,nlm11,nlm12,nlm21,nlm22,i
      double precision sumt,sumq,sumh,xx
      integer io1,io2,l1,l2


C ... Remove SO part from potential
      if (lso .ne. 0) then
C     io1 = nkaph*nkaph*nlmh*nlmh
C     call zprm('ppiz+so',2,ppihhz,io1,io1,nsp*nspc)
        call ppi2z(6,nsp,nspc,nkaph,nkaph,nlmh,nlmh,xx,ppihhz,ppihhz)
        call ppi2z(6,nsp,nspc,nkaph,1+kmax,nlmh,nlma,xx,ppihpz,ppihpz)
        call ppi2z(6,nsp,nspc,1+kmax,1+kmax,nlma,nlma,xx,ppippz,ppippz)
      endif

#if DEBUG
      integer fopna
      i = fopna('ppi',-1,4)
      call info0(0,0,0,' PVGTKN reading ppi')
      read (i) ppihh
      read (i) ppipp
      read (i) ppihp
#endif

C      print *, '!!'
C      ppihh = 0
C      tauhh = 0
C      ppihp = 0
C      tauhp = 0
C      ppipp = 0
C      taupp = 0

      sumt = 0d0
      sumq = 0d0
      sumh = 0d0

C ... Pkl*Pkl
      do  i = 1, nsp
        do  k1 = 0, kmax
          do  k2 = 0, kmax
            if (lcplxp .eq. 0) then
              do  ilm1 = 1, nlma
                l1 = ll(ilm1)
                sumt = sumt + qpp(k1,k2,ilm1,ilm1,i)*taupp(k1,k2,l1,i)
                sumq = sumq + qpp(k1,k2,ilm1,ilm1,i)*sigpp(k1,k2,l1,i)
                do  ilm2 = 1, nlma
                  sumh = sumh + qpp(k1,k2,ilm1,ilm2,i)*ppipp(k1,k2,ilm1,ilm2,i)
                enddo
              enddo

            else
              do  ilm1 = 1, nlma
                l1 = ll(ilm1)
                sumt = sumt + qpp(k1,k2,ilm1,ilm1,i)*taupp(k1,k2,l1,i)
                sumq = sumq + qpp(k1,k2,ilm1,ilm1,i)*sigpp(k1,k2,l1,i)
                do  ilm2 = 1, nlma
                  sumh = sumh + qpp(k1,k2,ilm1,ilm2,i)*ppippz(k1,k2,ilm1,ilm2,i)
                enddo
              enddo
            endif
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
                if (lcplxp .eq. 0) then
                  do  ilm1 = nlm11, nlm12
                    l1 = ll(ilm1)
                    do  ilm2 = nlm21, nlm22
                      if (ilm1 .eq. ilm2) then
                        sumt = sumt + qhh(k1,k2,ilm1,ilm2,i)*tauhh(k1,k2,l1,i)
                        sumq = sumq + qhh(k1,k2,ilm1,ilm2,i)*sighh(k1,k2,l1,i)
                      endif
                      sumh = sumh+ qhh(k1,k2,ilm1,ilm2,i)*ppihh(k1,k2,ilm1,ilm2,i)
                    enddo
                  enddo

                else
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
                endif
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

              if (lcplxp .eq. 0) then
                do  ilm1 = nlm11, nlm12
                  l1 = ll(ilm1)
                  do  ilm2 = 1, nlma
                    if (ilm1 .eq. ilm2) then
                      sumt = sumt + qhp(k1,k2,ilm1,ilm2,i)*tauhp(k1,k2,l1,i)
                      sumq = sumq + qhp(k1,k2,ilm1,ilm2,i)*sighp(k1,k2,l1,i)
                    endif
                    sumh = sumh+qhp(k1,k2,ilm1,ilm2,i)*ppihp(k1,k2,ilm1,ilm2,i)
                  enddo
                enddo
              else
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
              endif

            enddo
          endif
        enddo
      enddo

      if (lso .ne. 0) then
C       io1 = nkaph*nkaph*nlmh*nlmh
C       call zprm('ppiz+so',2,ppihhz,io1,io1,nsp*nspc)
        call ppi2z(5,nsp,nspc,nkaph,nkaph,nlmh,nlmh,xx,ppihhz,ppihhz)
        call ppi2z(5,nsp,nspc,nkaph,1+kmax,nlmh,nlma,xx,ppihpz,ppihpz)
        call ppi2z(5,nsp,nspc,1+kmax,1+kmax,nlma,nlma,xx,ppippz,ppippz)
      endif


C     call info2(0,0,0,' pvgtkn sumh= %,6;6d sumt= %,6;6d ',sumh,sumt)

      end subroutine pvgtkn


C      subroutine pvgtk1(nkaph,nlmh,k1,k2,nlm11,nlm12,nlm21,nlm22,qhh,
C     .  ppihhz,sumh)
CC- Contribution to double-counting from off-diagonal spin blocks
CC ----------------------------------------------------------------------
CCi Inputs
CCi   nkaph :number of types of one l-quantum number in the basis
CCi   nkaph :number of NMTO energies
CCi   nlmh  :dimensions qhh and ppihhz
CCi   k1    :left orbital channel (block) for H*H
CCi   k2    :right orbital channel (block) for H*H
CCi   nlm11 :starting L for left channel
CCi   nlm12 :ending L for left channel
CCi   nlm21 :starting L for right channel
CCi   nlm22 :ending L for right channel
CCi   qhh   :qhh(..,2) head-head density matrix for this site, 12 block
CCi   ppihhz:ppihh(..,3..4) head-head matrix elements of LxSx+LySy
CCo Outputs
CCo   sumh
CCl Local variables
CCl         :
CCr Remarks
CCr
CCu Updates
CCu   17 Jun 05  First created
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer k1,k2,nkaph,nlmh,nlm11,nlm12,nlm21,nlm22
C      double complex qhh(nkaph,nkaph,nlmh,nlmh,2)
C      double complex ppihhz(nkaph,nkaph,nlmh,nlmh,3)
C      double precision sumh
CC ... Local parameters
C      integer ilm1,ilm2,l1,ll
C      double complex sumz
C
C      call zprm('qhh',2,qhh(1,1,1,1,2),nlmh,nlmh,nlmh)
C      call zprm('ppihh',2,ppihhz(1,1,1,1,3),nlmh,nlmh,nlmh)
C
C      sumz = 0
C      do  ilm1 = nlm11, nlm12
C      do  ilm2 = nlm21, nlm22
C        sumz = sumz +
C     .         qhh(k1,k2,ilm1,ilm2,2)*ppihhz(k1,k2,ilm1,ilm2,3)
C        enddo
C      enddo
C
C      print *, 'sumz=',sumz
C
CC      sumh = sumh + sumz
C
C      end

