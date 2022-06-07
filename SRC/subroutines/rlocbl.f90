module m_prlcb
contains
  subroutine prlcb2(job,ia,nkaph,nlmha,kmax,nlma,isp,cPkL, &
       nlmto,evec,ewgt,evl,qhh,qhp)
    use m_orbl,only: Orblib,ktab,ltab,offl,norb
    !- Add one and two-center terms to density coeffs
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   job   :0 accumulate local density-matrix
    !i         :1 accumulate local density-matrix weighted by energy
    !i   ia    :site of augmentation
    !i   nkaph :dimensions qhh,qhp
    !i   nlmha :dimensions qhh,qhp
    !i   kmax  :polynomial cutoff
    !i   nlma  :augmentation L-cutoff
    !i   isp   :spin channel
    !i   cPkL  :PkL expansion eigenvector at site ia.
    !i   nlmto :dimension of lmto component of basis
    !i   evec  :eigenvector
    !i   ewgt  :eigenvector weight
    !i   evl   :energy weight (job=1)
    !o Outputs
    !o   qhh   :one-center density-matrix for PkL expansion (job=0)
    !o         :energy-weighted matrix (job=1)
    !o   qhp   :two-center density-matrix for PkL expansion (job=0)
    !o         :energy-weighted matrix (job=1)
    !r Remarks
    !u Updates
    !u   05 Jul 08 (T. Kotani)
    !u             Option to accumulate energy-weighted output density
    !u   27 Aug 01 Extended to local orbitals.
    ! ----------------------------------------------------------------------
    implicit none
    integer :: job,ia,kmax,nkaph,isp,nlmto,nlma,nlmha
    double precision :: qhh(nkaph,nkaph,nlmha,nlmha,isp), &
         qhp(nkaph,0:kmax,nlmha,nlma,isp),ewgt !(numq)
    real(8),optional::evl
    double complex evec(nlmto),cPkL(0:kmax,nlma)
    ! ... Local parameters
    integer :: i1,i2,ilm1,ilm2,ilma,io1,io2,iq,k,ik1,ik2, &
         l1,l2,n0,nkap0,nlm11,nlm12,nlm21,nlm22
    parameter (n0=10,nkap0=3)
    integer ::  blks(n0*nkap0),ntab(n0*nkap0)
    ! orb,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
    double precision :: xx
    if (nlmto == 0) return
    call tcn('prlcb2')
    ! --- Loop over all orbitals centered at site ia, incl. local orbs ---
    call orblib(ia)! norb,ltab,ktab,offl
    !     Block into groups of consecutive l
    call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)
    do io1 = 1, norb
       l1  = ltab(io1)
       ik1 = ktab(io1)
       nlm11 = l1**2+1
       nlm12 = nlm11 + blks(io1)-1
       i1 = offl(io1) !  i1 = hamiltonian offset for first orbital in block
       do  ilm1 = nlm11, nlm12
          i1 = i1+1
          !     ... Accumulate products H*Pkl
          !          do  iq = 1, numq
          if (job == 0) then
             do  k = 0, kmax
                do  ilma = 1, nlma
                   qhp(ik1,k,ilm1,ilma,isp)= qhp(ik1,k,ilm1,ilma,isp) &
                        + 2d0*dconjg(evec(i1))*cPkL(k,ilma)*ewgt !(iq)
                enddo
             enddo
          else
             do  k = 0, kmax
                qhp(ik1,k,ilm1,ilm1,isp)= qhp(ik1,k,ilm1,ilm1,isp) &
                     + evl*2d0*dconjg(evec(i1))*cPkL(k,ilm1)*ewgt !(iq)
             enddo
          endif
          !          enddo
          !     ... Accumulate products H*H
          do  io2 = 1, norb
             l2  = ltab(io2)
             ik2 = ktab(io2)
             nlm21 = l2**2+1
             nlm22 = nlm21 + blks(io2)-1
             i2 = offl(io2) 
             do  ilm2 = nlm21, nlm22
                i2 = i2+1
                if (job == 0) then
                   !                do  iq = 1, numq
                   qhh(ik1,ik2,ilm1,ilm2,isp) = &
                        qhh(ik1,ik2,ilm1,ilm2,isp) + &
                        dconjg(evec(i1))*evec(i2)*ewgt !(iq)
                   !                enddo
                elseif (job == 1 .AND. ilm1 == ilm2) then
                   !                do  iq = 1, numq
                   qhh(ik1,ik2,ilm1,ilm2,isp) = &
                        qhh(ik1,ik2,ilm1,ilm2,isp) + &
                        evl*dconjg(evec(i1))*evec(i2)*ewgt !(iq)
                   !                enddo
                endif
             enddo
          enddo
       enddo
    enddo
    call tcx('prlcb2')
  end subroutine prlcb2
  !!
  subroutine prlcb3(job,kmax,nlma,isp,cPkL,ewgt,evl,qpp)
    !- Add to local density coefficients for one state
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   job   :0 accumulate local density-matrix
    !i         :1 accumulate local density-matrix weighted by energy
    !i   kmax  :polynomial cutoff in PkL expansion
    !i   nlma  :L cutoff in PkL expansion
    !i   isp   :spin channel
    !i   cPkL  :coefficients to PkL expansion of evec
    !i   numq  :number of trial fermi levels
    !i   ewgt  :eigenvector weights
    !i   evl   :energy weight (job=1)
    !o Outputs
    !o   qpp   :local density matrix for PkL expansion (job=0)
    !o         :energy-weighted local density matrix (job=1)
    !r Remarks
    !u Updates
    !u   05 Jul 08 (T. Kotani)
    !u             Option to accumulate energy-weighted output density
    ! ----------------------------------------------------------------------
    implicit none
    integer :: job,kmax,nlma,isp!,numq
    double complex cPkL(0:kmax,nlma)
    double precision :: qpp(0:kmax,0:kmax,nlma,nlma,isp),ewgt !(numq)
    real(8),optional::evl
    double precision :: fac
    integer :: iq,ilm2,ilm1,k1,k2
    call tcn('prlcb3')
    !      do  iq = 1, numq
    fac = ewgt !(iq)
    if (job == 1) fac = evl*ewgt !(iq)
    do  ilm2 = 1, nlma
       do  ilm1 = 1, nlma
          do  k1 = 0, kmax
             do  k2 = 0, kmax
                qpp(k1,k2,ilm1,ilm2,isp)= qpp(k1,k2,ilm1,ilm2,isp) &
                     + fac*dconjg(cPkL(k1,ilm1))*cPkL(k2,ilm2)
             enddo
          enddo
       enddo
    enddo
    !      enddo
    call tcx('prlcb3')
  end subroutine prlcb3
end module m_prlcb
!!
subroutine rlocbl ( ssite , sspec , lfrce , nbas , isp  &
  , q , ndham , ndimh , nspc , napw , igvapw ,  nevec &
       , evec , ewgt , evl , sv_p_osig , sv_p_otau , sv_p_oppi &
  , lekkl , sv_p_oqkkl , sv_p_oeqkkl , f )
  use m_struc_def,only: s_site,s_spec,s_rv1,s_cv1
  use m_prlcb,only:   prlcb2,prlcb3
  use m_lmfinit,only: iv_a_oidxcg,iv_a_ojcg,rv_a_ocy,rv_a_ocg,lat_alat,nkaph
  use m_lattic,only: lat_qlat
  use m_bstrux,only: bstrux
  !- Accumulates the local atomic densities.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec pos
  !i     Stored:    *
  !i     Passed to: bstrux
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: lmxa kmxt lmxb rsma
  !i     Stored:    *
  !i     Passed to: bstrux
  !i   slat  :struct for lattice information; see routine ulat
  !i     Elts read: ocg ojcg oidxcg ocy alat qlat
  !i     Stored:    *
  !i     Passed to: bstrux
  !i   lfrce :if nonzero, accumulate contribution to force
  !i   nbas  :size of basis
  !i   isp   :spin channel
  !i   q     :Bloch wave number
  !i   ndham :leanding dimension of evl
  !i   ndimh :dimension of evec
  !i   nspc  :2 for coupled spins; otherwise 1
  !i   napw  :number of G vectors in PW basis (gvlst2.f)
  !i   igvapw:G vectors in PW basis, units of qlat (gvlst2.f)
  !i   nevec :number of occupied eigenvectors
  !i   evec  :eigenvectors
  !i   ewgt  :eigenvector weights
  !i   evl   :eigenvalues
  !i   osig  :overlap matrix of P_kL
  !i   otau  :kinetic energy matrix of P_kL (not used here)
  !i   oppi  :potential matrix of P_kL
  !i   lcplxp=1 only now  :0 if ppi is real; 1 if ppi is complex
  !i   lekkl :0 do not accumulate oeqkkl; 1 do accumulate oeqkkl
  !o Outputs
  !o   oqkkl :local density-matrix; see Remarks
  !o   oeqkkl:local part of energy-weighted density matrix
  !o   f     :local contribution to forces is added
  !l Local variables
  !l   ispc  :the current spin index in the coupled spins case.
  !l         :Some quantities have no separate address space for each
  !l         :spin in the indepedent-spins case (evec,evl,ewgt) but do
  !l         :in the coupled-spins case.  A separate loop ispc=1..nspc
  !l         :must be added for the latter case
  !l         :ispc is the appropriate index for objects which distinguish
  !l         :spins in the spin-coupled case only
  !l   isp   :isp  is the appropriate index for objects which distinguish
  !l         :spins in the spin-uncoupled case only
  !l   ksp   :the current spin index in both independent and coupled
  !l         :spins cases.
  !l         :ksp is appropriate spin index for quantities that have
  !l         :separate address space for each spin in every case
  !l         :(potential- and density-like objects).
  !r Remarks
  !r   The qkkl are contractions of the proper density-matrix
  !r      Dij = {sum_n w_n evec*_in evec_jn}
  !r   and the coefficients to the one-center expansion of the wave
  !r   function inside the augmentation sphere
  !r     F~i = Fi + sum_kL C^i_kL (P~kL - PkL)
  !r   As usual, we neglect cross terms when making function products.
  !r   Thus function products are of the form
  !r     F~i F~j = Fi Fj +
  !r             = sum_kLk'L' C^i_kL (P~kL P~k'L' - PkL Pk'L') C^j_k'L'
  !r             = sum_kLk'L' C^i_kL (n1kLk'L' - n2kLk'L') C^j_k'L'
  !r   the qkkl are defined as, e.g.
  !r      qpp_kLk'L' = sum_ij D_ij C^i_kL C^j_k'L'
  !r   so that the local part of the output density is
  !r      n1 - n2 = sum_kLk'L' qpp_kLk'L' (n1kLk'L' - n2kLk'L')
  !u Updates
  !u   05 Jul 08 (T. Kotani) output density for new PW part
  !u             Option to accumulate energy-weighted output density
  !u   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
  !u   16 Jun 05 Makes spin-off-diagonal density matrix, noncollinear case
  !u   23 Dec 04 Extended to spin-coupled case
  !u    1 Sep 04 Adapted to handle complex ppi
  !u   25 Aug 04 Adapted to extended local orbitals
  !u   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
  !u   15 Feb 02 (ATP) Added MPI parallelization
  !u   27 Aug 01 Extended to local orbitals.
  !u   17 Jun 00 spin polarized
  !u   25 May 00 Adapted from nfp rloc_q.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: lfrce,nbas,isp,ndimh,nspc,nevec,lekkl, ndham,napw,igvapw(3,napw)
  type(s_cv1),target:: sv_p_oppi(3,1)
  type(s_rv1),target:: sv_p_otau(3,1), sv_p_osig(3,1), &
       sv_p_oeqkkl(3,1), sv_p_oqkkl(3,1)
  real(8),pointer:: OQPP(:), OQHP(:),OQHH(:),OEQPP(:),OEQHP(:),OEQHH(:),OSIGPP(:), OSIGHP(:)
  complex(8),pointer:: OPPIPP(:),OPPIHP(:)
  real(8):: q(3), ewgt(nevec) , evl(ndham,isp) , f(3,nbas)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  double complex evec(ndimh,nspc,nevec)
  integer :: is,nlmbx,nlmx,ktop0,npmx,nkap0,n0
  parameter (nlmbx=25,  nkap0=3, n0=10) !npmx=32,
  integer :: kmaxx,nlmax,igetss,nglob,nlmto !mp,
  double precision :: alat,qlat(3,3)
  integer :: ia,isa,ivec,kmax,lmxa,nlma,lmxha,nlmha,ispc,ksp
  integer:: ob , odb
  double precision :: pa(3),rsma,pi,tpiba
  complex(8),allocatable:: cPkL(:),da(:),wk(:)
  integer ::iwdummy, iaini,iaend
  complex(8),allocatable::w_ob(:),w_odb(:,:,:,:),force(:)
  if (nevec <= 0) return
  call tcn('rlocbl')
  !      nkaph = globalvariables%nkaph
  ! ... Find maximum sizes needed to allocate strux; allocate them
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
  if (nlmax > nlmx)  call rxi('rlocbl: nlmx < nlma=',nlmax)
  if (kmaxx > ktop0) call rxi('rlocbl: ktop0 < kmax=',kmax)
  allocate(w_ob(ndimh*nlmax*(kmaxx+1)), w_odb((kmaxx+1),nlmax,ndimh,3))
  if (lfrce /= 0) then
     allocate(force(3*nbas))
     force=0d0
  endif
  !! Loop over augmentation sites ---
  do ia = 1,nbas ! iaini,iaend
     isa= ssite(ia)%spec
     pa = ssite(ia)%pos(1:3)
     lmxa= sspec(isa)%lmxa
     if (lmxa == -1) cycle !goto 10
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
     OPPIPP => sv_p_oppi(1,ia)%cv
     OPPIHP => sv_p_oppi(2,ia)%cv
     OSIGPP => sv_p_osig(1,ia)%v
     OSIGHP => sv_p_osig(2,ia)%v
     !   --- Strux to expand all orbitals and their gradients at site ia ---
     call bstrux(1, ia, pa, rsma, q, kmax, nlma, ndimh, napw, igvapw, w_ob , w_odb )
     !   --- Loop over eigenstates ---
     !       In noncollinear case, isp=1 always => need internal ispc=1..2
     !       ksp is the current spin index in both cases:
     !       ksp = isp  in the collinear case
     !           = ispc in the noncollinear case
     !       whereas ispc is spin index in the noncoll case, but 1 for coll.
     do ivec = 1, nevec
        do ispc = 1, nspc
           ksp = max(ispc,isp)
           !!  ... Pkl expansion of eigenvector
           call rlocb1(ndimh,nlma,kmax, evec(1,ispc,ivec), w_ob, cPkL)
           !!  ... Add to local density coefficients for one state
           call prlcb3 ( job=0 , kmax=kmax , nlma=nlma , isp=ksp , cpkl=cpkl, ewgt=ewgt(ivec),   qpp=OQPP)
           call prlcb2 ( job=0 , ia=ia , nkaph=nkaph , &
                nlmha=nlmha , kmax=kmax, nlma=nlma, isp=ksp , cpkl=cpkl, nlmto=nlmto , &
                evec=evec(1,ispc,ivec), ewgt=ewgt(ivec),    qhh=OQHH, qhp=OQHP)
           if (lekkl == 1) then
              call prlcb3(1, kmax, nlma, ksp, cpkl, ewgt(ivec), evl(ivec,isp) , OEQPP)
              call prlcb2(1, ia, nkaph,nlmha,kmax,nlma ,ksp , cpkl , nlmto , &
                   evec(1,ispc,ivec), ewgt(ivec),  evl ( ivec , isp ) , OEQHH , OEQHP)
           endif
           !! ... Contribution to forces
           if (lfrce/=0) then
              call rxx(nspc.ne.1,'forces not implemented in noncoll case')
              call flocbl ( nbas , ia , kmax , nkaph , lmxha , nlmha , nlma, &
                   lmxa, nlmto, ndimh, ksp, evl(ivec,isp),evec(1,ispc,ivec),ewgt (ivec ),cpkl, w_odb,& 
                   OSIGPP, OSIGHP , OPPIPP, OPPIHP, force ) 
           endif
        enddo
     enddo
     !   10   continue
  enddo !end loop over ia
  if (lfrce /= 0) call daxpy(3*nbas,1d0,force,1,f,1)
  deallocate(cPkL,da,wk,w_ob,w_odb)
  if(lfrce /=0) deallocate(force)
  call tcx('rlocbl')
end subroutine rlocbl

subroutine rlocb1(ndimh,nlma,kmax,evec,b,  cPkL)
  !- Pkl expansion of wave function at one site
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ndimh :dimension of evec
  !i   nlma  :augmentation L-cutoff in PkL expansion
  !i   kmax  :k- cutoff in PkL expansion
  !i   evec  :eigenvector coefficients
  !i   b     :strux to expand of orbitals from other sites in PkL
  !i         :b = b(ndimh,nlma,0:kmax)
  !o Outputs
  !o   cPkL  :coefficients to PkL expansion of evec
  implicit none
  integer :: kmax,ndimh,nlma,i,k,ilma
  double complex b(ndimh,nlma,0:kmax),cPkL(0:kmax,nlma),evec(ndimh)
  call tcn('rlocb1')
  do  k = 0, kmax
     do  ilma = 1, nlma
        cPkL(k,ilma) =  sum(evec(1:ndimh)*b(1:ndimh,ilma,k))
     enddo
  enddo
  call tcx('rlocb1')
end subroutine rlocb1



