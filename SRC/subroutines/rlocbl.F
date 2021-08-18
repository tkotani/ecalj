      module m_prlcb
      contains
      subroutine prlcb2(job,ia,nkaph,iprmb,nlmha,kmax,nlma,isp,cPkL,
     .     nlmto,evec,ewgt,evl,qhh,qhp)
C- Add one and two-center terms to density coeffs
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :0 accumulate local density-matrix
Ci         :1 accumulate local density-matrix weighted by energy
Ci   ia    :site of augmentation
Ci   nkaph :dimensions qhh,qhp
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nlmha :dimensions qhh,qhp
Ci   kmax  :polynomial cutoff
Ci   nlma  :augmentation L-cutoff
Ci   isp   :spin channel
Ci   cPkL  :PkL expansion eigenvector at site ia.
Ci   nlmto :dimension of lmto component of basis
Ci   evec  :eigenvector
Ci   ewgt  :eigenvector weight
Ci   evl   :energy weight (job=1)
Co Outputs
Co   qhh   :one-center density-matrix for PkL expansion (job=0)
Co         :energy-weighted matrix (job=1)
Co   qhp   :two-center density-matrix for PkL expansion (job=0)
Co         :energy-weighted matrix (job=1)
Cr Remarks
Cu Updates
Cu   05 Jul 08 (T. Kotani)
Cu             Option to accumulate energy-weighted output density
Cu   27 Aug 01 Extended to local orbitals.
C ----------------------------------------------------------------------
      implicit none
      integer job,ia,kmax,nkaph,isp,nlmto,nlma,nlmha,iprmb(1)
      double precision qhh(nkaph,nkaph,nlmha,nlmha,isp),
     .     qhp(nkaph,0:kmax,nlmha,nlma,isp),ewgt !(numq)
      real(8),optional::evl
      double complex evec(nlmto),cPkL(0:kmax,nlma)
C ... Local parameters
      integer i1,i2,ilm1,ilm2,ilma,io1,io2,iq,k,ik1,ik2,
     .     l1,l2,n0,nkap0,nlm11,nlm12,nlm21,nlm22
      parameter (n0=10,nkap0=3)
      integer norb,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
     .     blks(n0*nkap0),ntab(n0*nkap0)
      double precision xx
      if (nlmto .eq. 0) return
      call tcn('prlcb2')
C --- Loop over all orbitals centered at site ia, incl. local orbs ---
      call orbl(ia,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
C     Block into groups of consecutive l
      call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)
      do io1 = 1, norb
        l1  = ltab(io1)
        ik1 = ktab(io1)
        nlm11 = l1**2+1
        nlm12 = nlm11 + blks(io1)-1
        i1 = offl(io1) !  i1 = hamiltonian offset for first orbital in block
        do  ilm1 = nlm11, nlm12
          i1 = i1+1
C     ... Accumulate products H*Pkl
c          do  iq = 1, numq
            if (job .eq. 0) then
              do  k = 0, kmax
                do  ilma = 1, nlma
                  qhp(ik1,k,ilm1,ilma,isp)= qhp(ik1,k,ilm1,ilma,isp)
     .            + 2d0*dconjg(evec(i1))*cPkL(k,ilma)*ewgt !(iq)
                enddo
              enddo
            else
              do  k = 0, kmax
                qhp(ik1,k,ilm1,ilm1,isp)= qhp(ik1,k,ilm1,ilm1,isp)
     .          + evl*2d0*dconjg(evec(i1))*cPkL(k,ilm1)*ewgt !(iq)
              enddo
            endif
c          enddo
C     ... Accumulate products H*H
          do  io2 = 1, norb
            l2  = ltab(io2)
            ik2 = ktab(io2)
            nlm21 = l2**2+1
            nlm22 = nlm21 + blks(io2)-1
            i2 = offl(io2) !   i2 = orbital index in iprmb order
            do  ilm2 = nlm21, nlm22
              i2 = i2+1
              if (job .eq. 0) then
c                do  iq = 1, numq
                  qhh(ik1,ik2,ilm1,ilm2,isp) =
     .            qhh(ik1,ik2,ilm1,ilm2,isp) +
     .            dconjg(evec(i1))*evec(i2)*ewgt !(iq)
c                enddo
              elseif (job .eq. 1 .and. ilm1 .eq. ilm2) then
c                do  iq = 1, numq
                  qhh(ik1,ik2,ilm1,ilm2,isp) =
     .            qhh(ik1,ik2,ilm1,ilm2,isp) +
     .            evl*dconjg(evec(i1))*evec(i2)*ewgt !(iq)
c                enddo
              endif
            enddo
          enddo
        enddo
      enddo
      call tcx('prlcb2')
      end subroutine prlcb2
!!
      subroutine prlcb3(job,kmax,nlma,isp,cPkL,ewgt,evl,qpp)
C- Add to local density coefficients for one state
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :0 accumulate local density-matrix
Ci         :1 accumulate local density-matrix weighted by energy
Ci   kmax  :polynomial cutoff in PkL expansion
Ci   nlma  :L cutoff in PkL expansion
Ci   isp   :spin channel
Ci   cPkL  :coefficients to PkL expansion of evec
Ci   numq  :number of trial fermi levels
Ci   ewgt  :eigenvector weights
Ci   evl   :energy weight (job=1)
Co Outputs
Co   qpp   :local density matrix for PkL expansion (job=0)
Co         :energy-weighted local density matrix (job=1)
Cr Remarks
Cu Updates
Cu   05 Jul 08 (T. Kotani)
Cu             Option to accumulate energy-weighted output density
C ----------------------------------------------------------------------
      implicit none
      integer job,kmax,nlma,isp!,numq
      double complex cPkL(0:kmax,nlma)
      double precision qpp(0:kmax,0:kmax,nlma,nlma,isp),ewgt !(numq)
      real(8),optional::evl
      double precision fac
      integer iq,ilm2,ilm1,k1,k2
      call tcn('prlcb3')
c      do  iq = 1, numq
        fac = ewgt !(iq)
        if (job .eq. 1) fac = evl*ewgt !(iq)
        do  ilm2 = 1, nlma
          do  ilm1 = 1, nlma
            do  k1 = 0, kmax
              do  k2 = 0, kmax
                qpp(k1,k2,ilm1,ilm2,isp)= qpp(k1,k2,ilm1,ilm2,isp)
     .          + fac*dconjg(cPkL(k1,ilm1))*cPkL(k2,ilm2)
              enddo
            enddo
          enddo
        enddo
c      enddo
      call tcx('prlcb3')
      end subroutine prlcb3
      end module m_prlcb
!!
      subroutine rlocbl ( ssite , sspec , lfrce , nbas , isp  ! slat ,
     .     , q , ndham , ndimh , nspc , napw , igvapw , iprmb , nevec 
     .     , evec , ewgt , evl , sv_p_osig , sv_p_otau , sv_p_oppi , lcplxp 
     .     , lekkl , sv_p_oqkkl , sv_p_oeqkkl , f )
      use m_struc_def,only: s_site,s_spec,s_rv1
      use m_prlcb,only:   prlcb2,prlcb3
      use m_lmfinit,only: iv_a_oidxcg,iv_a_ojcg,rv_a_ocy,rv_a_ocg,lat_alat,nkaph
      use m_lattic,only: lat_qlat
      use m_bstrux,only: bstrux
C- Accumulates the local atomic densities.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: bstrux
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa kmxt lmxb rsma
Ci     Stored:    *
Ci     Passed to: bstrux
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: ocg ojcg oidxcg ocy alat qlat
Ci     Stored:    *
Ci     Passed to: bstrux
Ci   lfrce :if nonzero, accumulate contribution to force
Ci   nbas  :size of basis
Ci   isp   :spin channel
Ci   q     :Bloch wave number
Ci   ndham :leanding dimension of evl
Ci   ndimh :dimension of evec
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   napw  :number of G vectors in PW basis (gvlst2.f)
Ci   igvapw:G vectors in PW basis, units of qlat (gvlst2.f)
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nevec :number of occupied eigenvectors
Ci   evec  :eigenvectors
Ci   ewgt  :eigenvector weights
Ci   evl   :eigenvalues
Ci   osig  :overlap matrix of P_kL
Ci   otau  :kinetic energy matrix of P_kL (not used here)
Ci   oppi  :potential matrix of P_kL
Ci   lcplxp:0 if ppi is real; 1 if ppi is complex
Ci   lekkl :0 do not accumulate oeqkkl; 1 do accumulate oeqkkl
Co Outputs
Co   oqkkl :local density-matrix; see Remarks
Co   oeqkkl:local part of energy-weighted density matrix
Co   f     :local contribution to forces is added
Cl Local variables
Cl   ispc  :the current spin index in the coupled spins case.
Cl         :Some quantities have no separate address space for each
Cl         :spin in the indepedent-spins case (evec,evl,ewgt) but do
Cl         :in the coupled-spins case.  A separate loop ispc=1..nspc
Cl         :must be added for the latter case
Cl         :ispc is the appropriate index for objects which distinguish
Cl         :spins in the spin-coupled case only
Cl   isp   :isp  is the appropriate index for objects which distinguish
Cl         :spins in the spin-uncoupled case only
Cl   ksp   :the current spin index in both independent and coupled
Cl         :spins cases.
Cl         :ksp is appropriate spin index for quantities that have
Cl         :separate address space for each spin in every case
Cl         :(potential- and density-like objects).
Cr Remarks
Cr   The qkkl are contractions of the proper density-matrix
Cr      Dij = {sum_n w_n evec*_in evec_jn}
Cr   and the coefficients to the one-center expansion of the wave
Cr   function inside the augmentation sphere
Cr     F~i = Fi + sum_kL C^i_kL (P~kL - PkL)
Cr   As usual, we neglect cross terms when making function products.
Cr   Thus function products are of the form
Cr     F~i F~j = Fi Fj +
Cr             = sum_kLk'L' C^i_kL (P~kL P~k'L' - PkL Pk'L') C^j_k'L'
Cr             = sum_kLk'L' C^i_kL (n1kLk'L' - n2kLk'L') C^j_k'L'
Cr   the qkkl are defined as, e.g.
Cr      qpp_kLk'L' = sum_ij D_ij C^i_kL C^j_k'L'
Cr   so that the local part of the output density is
Cr      n1 - n2 = sum_kLk'L' qpp_kLk'L' (n1kLk'L' - n2kLk'L')
Cu Updates
Cu   05 Jul 08 (T. Kotani) output density for new PW part
Cu             Option to accumulate energy-weighted output density
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   16 Jun 05 Makes spin-off-diagonal density matrix, noncollinear case
Cu   23 Dec 04 Extended to spin-coupled case
Cu    1 Sep 04 Adapted to handle complex ppi
Cu   25 Aug 04 Adapted to extended local orbitals
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   27 Aug 01 Extended to local orbitals.
Cu   17 Jun 00 spin polarized
Cu   25 May 00 Adapted from nfp rloc_q.f
C ----------------------------------------------------------------------
      implicit none
      integer lfrce,nbas,isp,ndimh,nspc,nevec,lcplxp,lekkl,
     .iprmb(1),ndham,napw,igvapw(3,napw)
      type(s_rv1),target:: sv_p_oppi(3,1), sv_p_otau(3,1), sv_p_osig(3,1),
     &                     sv_p_oeqkkl(3,1), sv_p_oqkkl(3,1)
      real(8),pointer:: OQPP(:), OQHP(:),OQHH(:),OEQPP(:),OEQHP(:),OEQHH(:),OPPIPP(:)
     &     ,OSIGPP(:),OPPIHP(:), OSIGHP(:)
      real(8):: q(3), ewgt(nevec) , evl(ndham,isp) , f(3,nbas)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      double complex evec(ndimh,nspc,nevec)
      integer is,nlmbx,nlmx,ktop0,npmx,nkap0,n0
      parameter (nlmbx=25,  nkap0=3, n0=10) !npmx=32,
      integer kmaxx,nlmax,igetss,nglob,nlmto !mp,
      double precision alat,qlat(3,3)
      integer ia,isa,ivec,kmax,lmxa,nlma,oxf,lmxha,nlmha,ispc,ksp
      integer:: ob , odb
      double precision pa(3),rsma,pi,tpiba
      complex(8),allocatable:: cPkL(:),da(:),wk(:)
      integer ::iwdummy, iaini,iaend
      complex(8),allocatable::w_ob(:),w_odb(:),w_oxf(:)
      if (nevec .le. 0) return
      call tcn('rlocbl')
c      nkaph = globalvariables%nkaph
C ... Find maximum sizes needed to allocate strux; allocate them
      nlmax = 0
      kmaxx = 0
      do  ia = 1, nbas
        isa = int(ssite(ia)%spec)
        lmxa=sspec(isa)%lmxa
        kmax=sspec(isa)%kmxt
        nlma = (lmxa+1)**2
        kmaxx = max(kmaxx,kmax)
        nlmax = max(nlmax,nlma)
      enddo
      nlmto = ndimh - napw
      alat = lat_alat
      qlat = lat_qlat
      pi = 4d0*datan(1d0)
      tpiba = 2d0*pi/alat
      nlmx  = nlmax
      ktop0 = kmaxx
      allocate(cPkL((ktop0+1)*nlmx),da((ktop0+1)*nlmx*3),wk((ktop0+1)*nlmx))
      if (nlmax .gt. nlmx)  call rxi('rlocbl: nlmx.lt.nlma=',nlmax)
      if (kmaxx .gt. ktop0) call rxi('rlocbl: ktop0.lt.kmax=',kmax)
      allocate(w_ob(ndimh*nlmax*(kmaxx+1)), w_odb(ndimh*nlmax*(kmaxx+1)*3))
      if (lfrce .ne. 0) then
        allocate(w_oxf(3*nbas))
        w_oxf=0d0
      endif
!! Loop over augmentation sites ---
      do ia = 1,nbas ! iaini,iaend
        isa= ssite(ia)%spec
        pa = ssite(ia)%pos(1:3)
        lmxa= sspec(isa)%lmxa
        if (lmxa .eq. -1) cycle !goto 10
        lmxha=sspec(isa)%lmxb
        kmax= sspec(isa)%kmxt
        rsma= sspec(isa)%rsma
        nlmha = (lmxha+1)**2
        nlma  = (lmxa+1)**2
          OQPP => sv_p_oqkkl(1,ia)%v
          OQHP => sv_p_oqkkl(2,ia)%v   
          OQHH => sv_p_oqkkl(3,ia)%v  
         OEQPP => sv_p_oeqkkl(1,ia)%v
         OEQHP => sv_p_oeqkkl(2,ia)%v
         OEQHH => sv_p_oeqkkl(3,ia)%v
        OPPIPP => sv_p_oppi(1,ia)%v 
        OSIGPP => sv_p_osig(1,ia)%v 
        OPPIHP => sv_p_oppi(2,ia)%v 
        OSIGHP => sv_p_osig(2,ia)%v 
C   --- Strux to expand all orbitals and their gradients at site ia ---
        call bstrux ( 1 , ia , pa , rsma , q , 
     .      kmax , nlma , ndimh , napw , igvapw , w_ob , w_odb )
C   --- Loop over eigenstates ---
C       In noncollinear case, isp=1 always => need internal ispc=1..2
C       ksp is the current spin index in both cases:
C       ksp = isp  in the collinear case
C           = ispc in the noncollinear case
C       whereas ispc is spin index in the noncoll case, but 1 for coll.
        do ivec = 1, nevec
        do ispc = 1, nspc
            ksp = max(ispc,isp)
!!  ... Pkl expansion of eigenvector
            call rlocb1(ndimh,nlma,kmax, evec(1,ispc,ivec), w_ob, cPkL)
!!  ... Add to local density coefficients for one state
            call prlcb3 ( job=0 , kmax=kmax , nlma=nlma , isp=ksp ,
     .        cpkl=cpkl, ewgt=ewgt(ivec),                                qpp=OQPP)
            call prlcb2 ( job=0 , ia=ia , nkaph=nkaph , iprmb=iprmb ,
     .         nlmha=nlmha , kmax=kmax, nlma=nlma, isp=ksp , cpkl=cpkl, nlmto=nlmto , 
     .         evec=evec(1,ispc,ivec), ewgt=ewgt(ivec),    qhh=OQHH, qhp=OQHP)
            if (lekkl .eq. 1) then
               call prlcb3(1, kmax, nlma, ksp, cpkl, ewgt(ivec), evl(ivec,isp) , OEQPP)
               call prlcb2(1, ia, nkaph,iprmb,nlmha,kmax,nlma ,ksp , cpkl , nlmto ,  
     .              evec(1,ispc,ivec), ewgt(ivec),  evl ( ivec , isp ) , OEQHH , OEQHP)
            endif
!! ... Contribution to forces
            if (lfrce .ne. 0) then
              call rxx(nspc.ne.1,'forces not implemented in noncoll case')
              call flocbl ( nbas , ia , kmax , nkaph , lmxha , nlmha , nlma 
     .            , lmxa , nlmto , ndimh , iprmb , ksp , evl ( ivec , isp ) , evec 
     .            ( 1 , ispc , ivec ) , ewgt (ivec ) ,  cpkl , w_odb
     .            , da , wk , OPPIPP , OPPIPP , OSIGPP , OPPIHP
     .            , OPPIHP , OSIGHP , lcplxp , w_oxf )
            endif
        enddo
        enddo
c   10   continue
      enddo !end loop over ia
      if (lfrce .ne. 0) call daxpy(3*nbas,1d0,w_oxf,1,f,1)
      deallocate(cPkL,da,wk,w_ob,w_odb)
      if(lfrce /=0) deallocate(w_oxf)
      call tcx('rlocbl')
      end subroutine rlocbl

      subroutine rlocb1(ndimh,nlma,kmax,evec,b,  cPkL)
C- Pkl expansion of wave function at one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   ndimh :dimension of evec
Ci   nlma  :augmentation L-cutoff in PkL expansion
Ci   kmax  :k- cutoff in PkL expansion
Ci   evec  :eigenvector coefficients
Ci   b     :strux to expand of orbitals from other sites in PkL
Ci         :b = b(ndimh,nlma,0:kmax)
Co Outputs
Co   cPkL  :coefficients to PkL expansion of evec
      implicit none
      integer kmax,ndimh,nlma,i,k,ilma
      double complex b(ndimh,nlma,0:kmax),cPkL(0:kmax,nlma),evec(ndimh)
      call tcn('rlocb1')
      do  k = 0, kmax
      do  ilma = 1, nlma
         cPkL(k,ilma) =  sum(evec(1:ndimh)*b(1:ndimh,ilma,k))
      enddo
      enddo
      call tcx('rlocb1')
      end subroutine rlocb1



