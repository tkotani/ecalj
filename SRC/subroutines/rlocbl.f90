module m_rlocbl
  public rlocbl 
  private
contains
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
  !o   qhh   :one-center density-matrix for PkL expansion
  !o   eqhh :energy-weighted matrix 
  !o   qhp   :two-center density-matrix for PkL expansion 
  !o   eqhh  :energy-weighted matrix
  !o   qpp   :local density matrix for PkL expansion 
  !o  eqpp   :energy-weighted local density matrix 
  subroutine rlocbl(lfrce,nbas,isp, q,ndham,ndimh,nspc,napw,igvapw, nevec &
       ,evec,ewgt,evl,sv_p_osig,sv_p_otau,sv_p_oppi,lekkl,sv_p_oqkkl,sv_p_oeqkkl,f )
    use m_struc_def,only: s_spec,s_rv1,s_cv1,s_rv5
    use m_lmfinit,only: lat_alat,nkaph,ispec,sspec=>v_sspec
    use m_lattic,only: lat_qlat,rv_a_opos
    use m_bstrux,only: bstrux_set,bstr,dbstr
    !- Accumulates the local atomic densities.
    ! ----------------------------------------------------------------------
    !i Inputs
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
    type(s_rv1),target:: sv_p_otau(3,1), sv_p_osig(3,1)
    type(s_rv5),target:: sv_p_oeqkkl(3,1), sv_p_oqkkl(3,1)
    real(8),pointer:: qpp(:,:,:,:,:),eqpp(:,:,:,:,:),&
         qhp(:,:,:,:,:),qhh(:,:,:,:,:),eqhp(:,:,:,:,:),eqhh(:,:,:,:,:),&
         OSIGPP(:), OSIGHP(:)
    complex(8),pointer:: OPPIPP(:),OPPIHP(:)
    real(8):: q(3),f(3,nbas)
    real(8),target:: ewgt(nevec),evl(ndham,isp)
    complex(8),target:: evec(ndimh,nspc,nevec)
    integer :: is,nlmbx,nlmx,ktop0,npmx,nkap0,n0
    parameter (nlmbx=25,  nkap0=3, n0=10) !npmx=32,
    integer :: kmaxx,nlmax,igetss,nglob,nlmto !mp,
    double precision :: alat,qlat(3,3)
    integer :: ia,isa,ivec,kmax,lmxa,nlma,lmxha,nlmha,ispc,ksp
    integer:: ob,odb
    double precision :: pa(3),pi,tpiba
    complex(8),allocatable:: cPkL(:),da(:),wk(:)
    integer ::iwdummy, iaini,iaend
    complex(8),allocatable::w_ob(:),w_odb(:,:,:,:)
    real(8),allocatable:: force(:)
    if (nevec <= 0) return
    call tcn('rlocbl')
    ! ... Find maximum sizes needed to allocate strux; allocate them
    nlmax = 0
    kmaxx = 0
    do  ia = 1, nbas
       isa = ispec(ia) !(ssite(ia)%spec)
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
    allocate(da((ktop0+1)*nlmx*3),wk((ktop0+1)*nlmx)) !cPkL((ktop0+1)*nlmx),
    if (nlmax > nlmx)  call rxi('rlocbl: nlmx < nlma=',nlmax)
    if (kmaxx > ktop0) call rxi('rlocbl: ktop0 < kmax=',kmax)
    allocate(w_ob(ndimh*nlmax*(kmaxx+1)), w_odb((kmaxx+1),nlmax,ndimh,3))
    if (lfrce /= 0) then
       allocate(force(3*nbas))
       force=0d0
    endif
    ialoop: do ia = 1,nbas !Loop over augmentation sites ---
       isa = ispec(ia) 
       pa  = rv_a_opos(:,ia)
       lmxa= sspec(isa)%lmxa
       if (lmxa == -1) cycle
       lmxha=sspec(isa)%lmxb
       kmax= sspec(isa)%kmxt
       nlmha = (lmxha+1)**2
       nlma  = (lmxa+1)**2
       call bstrux_set(ia,q) !Get bstr
       
       OPPIPP => sv_p_oppi(1,ia)%cv
       OPPIHP => sv_p_oppi(2,ia)%cv
       OSIGPP => sv_p_osig(1,ia)%v
       OSIGHP => sv_p_osig(2,ia)%v

       qpp  => sv_p_oqkkl(1,ia)%v
       qhp  => sv_p_oqkkl(2,ia)%v
       qhh  => sv_p_oqkkl(3,ia)%v
       eqpp => sv_p_oeqkkl(1,ia)%v
       eqhp => sv_p_oeqkkl(2,ia)%v
       eqhh => sv_p_oeqkkl(3,ia)%v
       
       !   --- Loop over eigenstates ---
       !       In noncollinear case, isp=1 always => need internal ispc=1..2
       !       ksp is the current spin index in both cases:
       !       ksp = isp  in the collinear case
       !           = ispc in the noncollinear case
       !       whereas ispc is spin index in the noncoll case, but 1 for coll.
       cPklblock: block
         use m_orbl,only: Orblib, norb,ltab,ktab,offl, ntab,blks
         complex(8):: cPkl(0:kmax,nlma),wtt
         complex(8),pointer::evecc(:)
         real(8),pointer:: ewgtt,evll
         integer:: k, ilma,ilm1,ilm2,k1,k2,io1,l1,ik1,i1,io2,l2,ik2,i2
         do ivec = 1, nevec
            do ispc = 1, nspc
               ksp = max(ispc,isp)
               evecc=> evec(1:ndimh,ispc,ivec)
               ewgtt=> ewgt(ivec)
               evll => evl(ivec,isp)
               do  k = 0, kmax
                  cPkL(k,:) =  matmul(evec(1:ndimh,ispc,ivec),bstr(1:ndimh,:,k)) ! Pkl expansion of eigenvector
               enddo
               do  ilm2 = 1, nlma !Add to local density coefficients for one state
                  do  k2 = 0, kmax
                     qpp (:,k2,:,ilm2,ksp)= qpp(:,k2,:,ilm2,ksp) + ewgtt*dconjg(cPkL(:,:))*cPkL(k2,ilm2)
                     eqpp(:,k2,:,ilm2,ksp)= eqpp(:,k2,:,ilm2,ksp)+ ewgtt*dconjg(cPkL(:,:))*cPkL(k2,ilm2)*evll
                  enddo
               enddo
               call orblib(ia)! norb,ltab,ktab,offl      !     Block into groups of consecutive l
               do io1 = 1, norb
                  l1  = ltab(io1)
                  ik1 = ktab(io1)
                  i1 = offl(io1) !  i1 = hamiltonian offset for first orbital in block
                  do ilm1 = l1**2+1, l1**2+1 + blks(io1)-1
                     i1 = i1+1
                     qhp (ik1,:,ilm1,:,ksp)    = qhp (ik1,:,ilm1,:,ksp)+ 2d0*dconjg(evecc(i1))*cPkL(:,:)*ewgtt 
                     eqhp(ik1,:,ilm1,ilm1,ksp)= eqhp(ik1,:,ilm1,ilm1,ksp) + evll*2d0*dconjg(evecc(i1))*cPkL(:,ilm1)*ewgtt 
                     do  io2 = 1, norb
                        l2  = ltab(io2)
                        ik2 = ktab(io2)
                        i2 = offl(io2) 
                        do  ilm2 = l2**2+1, l2**2+1 + blks(io2)-1
                           i2 = i2+1
                           qhh(ik1,ik2,ilm1,ilm2,ksp) =  qhh(ik1,ik2,ilm1,ilm2,ksp) + dconjg(evecc(i1))*evecc(i2)*ewgtt !(iq)
                           if(ilm1 == ilm2) then
                              eqhh(ik1,ik2,ilm1,ilm2,ksp) = eqhh(ik1,ik2,ilm1,ilm2,ksp) +evll*dconjg(evecc(i1))*evecc(i2)*ewgtt
                           endif
                        enddo
                     enddo
                  enddo
               enddo
               if (lfrce/=0) then! ... Contribution to forces
                  call rxx(nspc.ne.1,'forces not implemented in noncoll case')
                  call flocbl( nbas,ia,kmax,nkaph,lmxha,nlmha,nlma, &
                       lmxa,nlmto,ndimh,ksp,evll,evecc,ewgtt,cpkl,dbstr,&  
                       OSIGPP, OSIGHP,OPPIPP,OPPIHP, force ) 
               endif
            enddo
         enddo
       endblock cPklblock
    enddo ialoop
    if (lfrce /= 0) call daxpy(3*nbas,1d0,force,1,f,1)
    deallocate(da,wk,w_ob,w_odb) !cPkL,
    if(lfrce /=0) deallocate(force)
    call tcx('rlocbl')
  end subroutine rlocbl
  subroutine flocbl(nbas,ia,kmax,nkaph,lmxha,nlmha,nlma,lmxa,nlmto, &
       ndimh,isp,evll,evecc,ewgtt,cPkL,db, sigpp,sighp,ppippz,ppihpz,f)
    use m_orbl,only: Orblib, norb,ltab,ktab,offl, ntab,blks
    !- Force contribution from augmentation at site ia.
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nbas  :size of basis
    !i   ia    :site of augmentation
    !i   kmax  :polynomial cutoff
    !i   nkaph :leading dimension of sigph,ppihp at site ia
    !i   nlmha :dimensions ppihp
    !i   nlma  :augmentation L-cutoff
    !i   lmxa  :augmentation L-cutoff
    !i   ndimh :dimension of hamiltonian
    !i   evll   :eigenvalue
    !i   evecc  :eigenvector
    !i   ewgtt  :eigenvector weight
    !i   cPkL  :PkL expansion eigenvector at site ia.
    !i   b     :structure constants, needed for PW contribution
    !i   db    :gradient of structure constants, needed to make grad psi
    !i   da    :work array holding grad psi
    !i   wk    :work array holding (ppi-evll*sig)*evecc
    !i   ppipp :local tail-tail potential matrix
    !i   ppippz:local tail-tail potential matrix, complex form
    !i         :NB: only ppipp or ppippz is used, depending on lcplxp
    !i   sigpp :local tail-tail overlap matrix
    !i   ppihp :local head-tail potential matrix
    !i   ppihpz:local head-tail potential matrix, complex form
    !i         :NB: only ppihp or ppihpz is used, depending on lcplxp
    !i   sighp :local head-tail overlap matrix
    !i   lcplxp=1 only now: ppi is complex
    !o Outputs
    !o   f     :local contribution to forces is added
    !r Remarks
    !u Updates
    !u   05 Jul 08 (T. Kotani) Contribution from PW part of basis
    !u    1 Sep 04 Adapted to handle complex ppi
    !u   17 Jun 00 spin polarized
    !u   25 May 00 Adapted from nfp floc_q.f
    ! ----------------------------------------------------------------------
    implicit none
    integer :: ia,kmax,lmxa,nkaph,nbas,nlmto,ndimh,nlma,lmxha,nlmha,isp
    double precision :: evll,f(3,nbas),ewgtt, &
         sigpp(0:kmax,0:kmax,0:lmxa,isp),sighp(nkaph,0:kmax,0:lmxha,isp)
    double complex db(ndimh,nlma,0:kmax,3),da(0:kmax,nlma,3), &
         evecc(ndimh),cPkL(0:kmax,nlma),wk(0:kmax,nlma)
    double complex ppihpz(nkaph,0:kmax,nlmha,nlma,isp)
    double complex ppippz(0:kmax,0:kmax,nlma,nlma,isp)
    integer :: i,ib,ilm,ilmb,io,iq,k,l2,m,nlm1,n0,nkap0,ik,ilma,jlm,k1,k2,l,ll
    parameter (n0=10,nkap0=3)
    integer ::oi,iblk ! blks(n0*nkap0),ntab(n0*nkap0),
    double precision :: wt,xx,ssum(3)
    if (nlmto == 0) return
    call tcn('flocbl')
    ! ... Make wk=(ppi-evll*sig)*psi 
!    call flocb2(ia,nlmto,kmax,nkaph,nlmha,nlma,evll,evecc, &
!         ppippz(0,0,1,1,isp),sigpp(0,0,0,isp),ppihpz(1,0,1,1,isp),sighp(1,0,0,isp), cPkL,wk)
    ! ... Add Hsm*Pkl block of ppi-evll*sig times evecc
    !     Block evll*ppi contribution in groups of consecutive l
    call orblib(ia)!norb,ltab,ktab,offl
    wk=0d0
    do io = 1, norb  ! orbital index      
       l  = ltab(io) !  l index
       ik = ktab(io) !  radial index
       nlm1 = l**2+1
       do  ilmb = nlm1, (l+1)**2   ! evll*sig contribution requires explicit knowledge of l
          wk(:,ilmb) = wk(:,ilmb) - evll*sighp(ik,:,l,isp)*evecc(offl(io)+1+ilmb-nlm1)
       enddo
       if (blks(io) /= 0) then ! ppi contribution: loop over largest blocks possible
          do ilmb = nlm1, nlm1 + blks(io)-1
             wk(:,:) = wk(:,:) + ppihpz(ik,:,ilmb,:,isp)*evecc(offl(io)+1+ilmb-nlm1) 
          enddo
       endif
    enddo
    do  ilm = 1, nlma !Add Pkl*Pkl block of ppi-evll*sig times cPkL
       l = ll(ilm)
       wk(:,ilm) = wk(:,ilm) +[(sum(ppippz(k1,:,ilm,:,isp)*cPkL(:,:)) - evll*sum(sigpp(k1,:,l,isp)*cPkL(:,ilm)),k1=0,kmax)]
    enddo
    do ib = 1, nbas ! ... Loop over ib, virtual shift of wavefct part centered there
       if (ib == ia) cycle
       call orblib(ib)!norb,ltab,ktab,offl  !   ... Grad of psi expansion coeffs from a virtual shift at site ib
       da=0d0 !da can be written as {d cPkL}/{d R}.See rlocb1 to make cPkL. b is used instead of db
       do io = 1, norb
          oi = offl(io)  
          do iblk=1,blks(io)  !if blsk(io)=0, no loop cycle
             da(0:kmax,:,1:3)=da(0:kmax,:,1:3)+ evecc(oi+iblk)*reshape(db(oi+iblk,:,0:kmax,1:3),shape(da),order=[2,1,3])
          enddo
       enddo
       ! --- Force term is (grad psi_kL) * (ppi-evll*sig)*evecc ---
       ssum = ewgtt* 2d0* [(sum(dconjg(da(:,:,m))*wk(:,:)),m=1,3)]
       f(:,ib) = f(:,ib) - ssum
       f(:,ia) = f(:,ia) + ssum
    enddo
    ! --- Force at site ia from PWs ---
    da=0d0
    do  i = nlmto+1, ndimh
       da(:,:,:) = da(:,:,:) + evecc(i)*reshape(db(i,:,:,:),shape=shape(da),order=[2,1,3]) !db(ilm,:,:)
    enddo
    f(:,ia) = f(:,ia) + ewgtt*2d0*[(sum(dconjg(da(:,:,m))*wk(:,:)),m=1,3)] !Force term is (grad psi_kL) * (ppi-evll*sig)*evecc
    call tcx('flocbl')
  end subroutine flocbl

!   subroutine flocb2(ia,nlmto,kmax,nkaph,nlmha,nlma,evl,evec, &
!        ppippz,sigpp,ppihpz,sighp,cPkL,wk)!lcplxp,
!     use m_orbl,only: Orblib, norb,ltab,ktab,offl, ntab,blks
!     !- Make (ppi-evl*sig)*evec
!     ! ----------------------------------------------------------------------
!     !i Inputs
!     !i   ia    :site of augmentation
!     !l   nlmto :dimension of LMTO part of hamiltonian, for hh and ht blocks
!     !i   nlmha :dimensions ppihp
!     !i   nlma  :augmentation L-cutoff
!     !i   kmax  :polynomial cutoff
!     !i   evl   :eigenvalue
!     !i   evec  :eigenvector
!     !i   nkaph :leading dimension of sigph,ppihp
!     !i   ppipp :local tail-tail potential matrix
!     !i   ppippz:local tail-tail potential matrix, complex form
!     !i         :NB: only ppipp or ppippz is used, depending on lcplxp
!     !i   sigpp :local tail-tail overlap matrix
!     !i   ppihp :local head-tail potential matrix
!     !i   ppihpz:local head-tail potential matrix, complex form
!     !i         :NB: only ppihp or ppihpz is used, depending on lcplxp
!     !i   sighp :local head-tail overlap matrix
!     !i   lcplxp=1 only now: ppi is complex
!     !i   cPkL  :PkL expansion eigenvector at site ia.
!     !o Outputs
!     !o   wk    :(ppi-evl*sig)*evec
!     !r Remarks
!     !u Updates
!     ! ----------------------------------------------------------------------
!     !     implicit none
!     ! ... Passed parameters
!     integer :: ia,kmax,nkaph,nlma,nlmha,nlmto
!     double precision :: &
!          evl,sighp(nkaph,0:kmax,0:*),& ! & ppihp(nkaph,0:kmax,nlmha,nlma),
!          sigpp(0:kmax,0:kmax,0:*)!,ppipp(0:kmax,0:kmax,nlma,nlma)
!     double complex wk(0:kmax,nlma),evec(1),cPkL(0:kmax,nlma)
!     double complex ppihpz(nkaph,0:kmax,nlmha,nlma)
!     double complex ppippz(0:kmax,0:kmax,nlma,nlma)
!     ! ... Local parameters
!     integer :: i,ilm,ilma,ilmb,io,jlm,k,k1,k2,l,ll,nlm1,nlm2,n0,nkap0,ik
!     !     .norb
!     parameter (n0=10,nkap0=3)
! !    integer :: blks(n0*nkap0),ntab(n0*nkap0)!ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
!     double precision :: xx
!     wk=0d0
!     ! ... Add Hsm*Pkl block of ppi-evl*sig times evec
!     call orblib(ia)!norb,ltab,ktab,offl
!     !     Block evl*ppi contribution in groups of consecutive l
! !    call gtbsl4(ia) !1(4,norb,ltab,ktab,xx,xx,ntab,blks)
!     do  io = 1, norb
!        !       l,ik = l and kaph indices, needed for sigma
!        l  = ltab(io)
!        ik = ktab(io)
!        nlm1 = l**2+1
!        nlm2 = (l+1)**2
!        i  = offl(io)
!        !   ... evl*sig contribution requires explicit knowledge of l
!        do  ilmb = nlm1, nlm2
!           i = i+1
!           do  k = 0, kmax
!              wk(k,ilmb) = wk(k,ilmb) - evl*sighp(ik,k,l)*evec(i)
!           enddo
!        enddo
!        !   ... ppi contribution: loop over largest blocks possible
!        if (blks(io) /= 0) then ! .AND. lcplxp == 1) then
!           nlm2 = nlm1 + blks(io)-1
!           i  = offl(io)
!           do  ilmb = nlm1, nlm2
!              i = i+1
!              do  ilma = 1, nlma
!                 do  k = 0, kmax
!                    wk(k,ilma) = wk(k,ilma) + ppihpz(ik,k,ilmb,ilma)*evec(i)
!                 enddo
!              enddo
!           enddo
!        endif
!     enddo
!     ! ... Add Pkl*Pkl block of ppi-evl*sig times cPkL
!     do  ilm = 1, nlma
!        l = ll(ilm)
!        do  k1 = 0, kmax
!           do  jlm = 1, nlma
!              do  k2 = 0, kmax
!                 wk(k1,ilm) = wk(k1,ilm)+ppippz(k1,k2,ilm,jlm)*cPkL(k2,jlm)
!              enddo
!           enddo
!           do  k2 = 0, kmax
!              wk(k1,ilm) = wk(k1,ilm) - evl*sigpp(k1,k2,l)*cPkL(k2,ilm)
!           enddo
!        enddo
!     enddo
!   end subroutine flocb2
end module m_rlocbl
