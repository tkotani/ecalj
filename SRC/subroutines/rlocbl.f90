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
    real(8),pointer:: qpp(:,:,:,:,:),eqpp(:,:,:,:,:),qhp(:,:,:,:,:),qhh(:,:,:,:,:),eqhp(:,:,:,:,:),eqhh(:,:,:,:,:)
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
    real(8),allocatable:: force(:,:)
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
       allocate(force(3,nbas))
       force=0d0
    endif
    ialoop: do ia = 1,nbas !Loop over augmentation sites 
       isa = ispec(ia) 
       pa  = rv_a_opos(:,ia)
       lmxa= sspec(isa)%lmxa
       if (lmxa == -1) cycle
       lmxha=sspec(isa)%lmxb
       kmax= sspec(isa)%kmxt
       nlmha = (lmxha+1)**2
       nlma  = (lmxa+1)**2
       call bstrux_set(ia,q) !Get bstr
       qpp  => sv_p_oqkkl(1,ia)%v
       qhp  => sv_p_oqkkl(2,ia)%v
       qhh  => sv_p_oqkkl(3,ia)%v
       eqpp => sv_p_oeqkkl(1,ia)%v
       eqhp => sv_p_oeqkkl(2,ia)%v
       eqhh => sv_p_oeqkkl(3,ia)%v
       !       In noncollinear case, isp=1 always => need internal ispc=1..2
       !       ksp is the current spin index in both cases: ksp = isp in the collinear case;  ksp= ispc in the noncollinear case,
       !       whereas ispc is spin index in the noncoll case, but 1 for coll.
       cPklblock: block
         use m_orbl,only: Orblib, norb,ltab,ktab,offl, ntab,blks
         use m_lmfinit,only:nsp
         complex(8):: cPkl(0:kmax,nlma),wtt
         complex(8),pointer::evecc(:)
         real(8),pointer:: ewgtt,evll
         integer:: k, ilma,ilm1,ilm2,k1,k2,io1,l1,ik1,i1,io2,l2,ik2,i2
         real(8)::    sigpp(0:kmax,0:kmax,0:lmxa,nsp),sighp(nkaph,0:kmax,0:lmxha,nsp)
         complex(8):: ppihpz(nkaph,0:kmax,nlmha,nlma,nsp),  ppippz(0:kmax,0:kmax,nlma,nlma,nsp)
         sigpp= reshape(sv_p_osig(1,ia)%v,shape(sigpp))
         sighp= reshape(sv_p_osig(2,ia)%v,shape(sighp))
         ppippz=reshape(sv_p_oppi(1,ia)%cv,shape(ppippz))
         ppihpz=reshape(sv_p_oppi(2,ia)%cv,shape(ppihpz))
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
               lfrceT: if(lfrce/=0) then! ... Contribution to forces
                  call rxx(nspc.ne.1,'forces not implemented in noncoll case')
                  flocblblock: block !& !Force contribution from augmentation at site ia.
                    use m_orbl,only: Orblib, norb,ltab,ktab,offl, ntab,blks
                    integer :: isp,i,ib,ilm,ilmb,io,iq,k,l2,m,nlm1,nkap0,ik,ilma,jlm,k1,k2,l,ll,oi,iblk 
                    complex(8):: da(0:kmax,nlma,3), wk(0:kmax,nlma) !wk=(ppi-evll*sig)*psi 
                    real(8):: ssum(3)
                    isp=ksp
                    call orblib(ia) !norb,ltab,ktab,offl
                    wk=0d0
                    do io = 1, norb  ! orbital index      !Block evll*ppi contribution in groups of consecutive l
                       l  = ltab(io); ik = ktab(io)  !  l index, radial index
                       nlm1 = l**2+1
                       do  ilmb = nlm1, (l+1)**2   ! evll*sig contribution requires explicit knowledge of l
                          wk(:,ilmb) = wk(:,ilmb) - evll*sighp(ik,:,l,isp)*evecc(offl(io)+1+ilmb-nlm1) ! Add Hsm*Pkl block of ppi-evll*sig times evecc     
                       enddo
                       do ilmb = nlm1, nlm1 + blks(io)-1 ! ppi contribution: loop over largest blocks possible !blks(io) /= 0
                          wk(:,:) = wk(:,:) + ppihpz(ik,:,ilmb,:,isp)*evecc(offl(io)+1+ilmb-nlm1) 
                       enddo
                    enddo
                    do  ilm = 1, nlma !Add Pkl*Pkl block of ppi-evll*sig times cPkL
                       wk(:,ilm)= wk(:,ilm)+ &
                            [(sum(ppippz(k1,:,ilm,:,isp)*cPkL(:,:)) -evll*sum(sigpp(k1,:,ll(ilm),isp)*cPkL(:,ilm)), k1=0,kmax)]
                    enddo
                    do ib = 1, nbas ! ... Loop over ib, virtual shift of wavefct part centered there
                       if (ib == ia) cycle
                       call orblib(ib)!norb,ltab,ktab,offl  !   ... Grad of psi expansion coeffs from a virtual shift at site ib
                       da=0d0 !da can be written as {d cPkL}/{d R}.See rlocb1 to make cPkL. b is used instead of db
                       do io = 1, norb
                          do iblk=offl(io)+1, offl(io)+blks(io)  !if blsk(io)=0, no loop cycle
                             da(0:kmax,:,1:3)=da(0:kmax,:,1:3)&
                                  + evecc(iblk)* reshape(dbstr(iblk,:,0:kmax,1:3),shape(da),order=[2,1,3]) !transpose 1st and 2nd
                          enddo
                       enddo
                       ssum = ewgtt* 2d0* [(sum(dconjg(da(:,:,m))*wk(:,:)),m=1,3)] ! --- Force term is (grad psi_kL) * (ppi-evll*sig)*evecc ---
                       force(:,ib) = force(:,ib) - ssum
                       force(:,ia) = force(:,ia) + ssum
                    enddo
                    da=0d0 !Force at site ia from PWs ---
                    do  i = nlmto+1, ndimh
                       da(:,:,:) = da(:,:,:) + evecc(i)*reshape(dbstr(i,:,:,:),shape=shape(da),order=[2,1,3]) !db(ilm,:,:)
                    enddo
                    force(:,ia) = force(:,ia) + ewgtt*2d0*[(sum(dconjg(da(:,:,m))*wk(:,:)),m=1,3)] !Force term is (grad psi_kL) * (ppi-evll*sig)*evecc
                  endblock flocblblock
               endif lfrceT
             enddo
          enddo
        endblock cPklblock
    enddo ialoop
    if(lfrce/=0) f = force+f 
    call tcx('rlocbl')
  end subroutine rlocbl
end module m_rlocbl
!   subroutine flocbl(nbas,ia,kmax,nkaph,lmxha,nlmha,nlma,lmxa,nlmto, & !Force contribution from augmentation at site ia.
!        ndimh,isp,evll,evecc,ewgtt,cPkL,db, sigpp,sighp,ppippz,ppihpz,f)
!     use m_lmfinit,only: nsp
!     use m_orbl,only: Orblib, norb,ltab,ktab,offl, ntab,blks
!     !i   nbas  :size of basis
!     !i   ia    :site of augmentation
!     !i   kmax  :polynomial cutoff
!     !i   nkaph :leading dimension of sigph,ppihp at site ia
!     !i   nlmha :dimensions ppihp
!     !i   nlma  :augmentation L-cutoff
!     !i   lmxa  :augmentation L-cutoff
!     !i   ndimh :dimension of hamiltonian
!     !i   evll   :eigenvalue
!     !i   evecc  :eigenvector
!     !i   ewgtt  :eigenvector weight
!     !i   cPkL  :PkL expansion eigenvector at site ia.
!     !i   b     :structure constants, needed for PW contribution
!     !i   db    :gradient of structure constants, needed to make grad psi
!     !i   da    :work array holding grad psi
!     !i   wk    :work array holding (ppi-evll*sig)*evecc
!     !i   ppipp :local tail-tail potential matrix
!     !i   ppippz:local tail-tail potential matrix, complex form
!     !i         :NB: only ppipp or ppippz is used, depending on lcplxp
!     !i   sigpp :local tail-tail overlap matrix
!     !i   ppihp :local head-tail potential matrix
!     !i   ppihpz:local head-tail potential matrix, complex form
!     !i         :NB: only ppihp or ppihpz is used, depending on lcplxp
!     !i   sighp :local head-tail overlap matrix
!     !i   lcplxp=1 only now: ppi is complex
!     !o Outputs
!     !o   f     :local contribution to forces is added
!     implicit none
!     integer :: ia,kmax,lmxa,nkaph,nbas,nlmto,ndimh,nlma,lmxha,nlmha,isp
!     real(8):: evll,f(3,nbas),ewgtt, sigpp(0:kmax,0:kmax,0:lmxa,nsp),sighp(nkaph,0:kmax,0:lmxha,nsp)
!     complex(8):: db(ndimh,nlma,0:kmax,3),da(0:kmax,nlma,3),  evecc(ndimh),cPkL(0:kmax,nlma),wk(0:kmax,nlma),&
!          ppihpz(nkaph,0:kmax,nlmha,nlma,nsp),     ppippz(0:kmax,0:kmax,nlma,nlma,nsp)
!     integer :: i,ib,ilm,ilmb,io,iq,k,l2,m,nlm1,nkap0,ik,ilma,jlm,k1,k2,l,ll
!     integer ::oi,iblk 
!     double precision :: wt,xx,ssum(3)
! !    if (nlmto == 0) return
! !    call tcn('flocbl')
!     ! ... Make wk=(ppi-evll*sig)*psi  ! Add Hsm*Pkl block of ppi-evll*sig times evecc     Block evll*ppi contribution in groups of consecutive l
!     call orblib(ia)!norb,ltab,ktab,offl
!     wk=0d0
!     do io = 1, norb  ! orbital index      
!        l  = ltab(io) !  l index
!        ik = ktab(io) !  radial index
!        nlm1 = l**2+1
!        do  ilmb = nlm1, (l+1)**2   ! evll*sig contribution requires explicit knowledge of l
!           wk(:,ilmb) = wk(:,ilmb) - evll*sighp(ik,:,l,isp)*evecc(offl(io)+1+ilmb-nlm1)
!        enddo
!        if (blks(io) /= 0) then ! ppi contribution: loop over largest blocks possible
!           do ilmb = nlm1, nlm1 + blks(io)-1
!              wk(:,:) = wk(:,:) + ppihpz(ik,:,ilmb,:,isp)*evecc(offl(io)+1+ilmb-nlm1) 
!           enddo
!        endif
!     enddo
!     do  ilm = 1, nlma !Add Pkl*Pkl block of ppi-evll*sig times cPkL
!        wk(:,ilm)= wk(:,ilm)+ [(sum(ppippz(k1,:,ilm,:,isp)*cPkL(:,:)) -evll*sum(sigpp(k1,:,ll(ilm),isp)*cPkL(:,ilm)), k1=0,kmax)]
!     enddo
!     do ib = 1, nbas ! ... Loop over ib, virtual shift of wavefct part centered there
!        if (ib == ia) cycle
!        call orblib(ib)!norb,ltab,ktab,offl  !   ... Grad of psi expansion coeffs from a virtual shift at site ib
!        da=0d0 !da can be written as {d cPkL}/{d R}.See rlocb1 to make cPkL. b is used instead of db
!        do io = 1, norb
!           do iblk=offl(io)+1, offl(io)+blks(io)  !if blsk(io)=0, no loop cycle
!              da(0:kmax,:,1:3)=da(0:kmax,:,1:3)+ evecc(iblk)*reshape(db(iblk,:,0:kmax,1:3),shape(da),order=[2,1,3])
!           enddo
!        enddo
!        ssum = ewgtt* 2d0* [(sum(dconjg(da(:,:,m))*wk(:,:)),m=1,3)] ! --- Force term is (grad psi_kL) * (ppi-evll*sig)*evecc ---
!        f(:,ib) = f(:,ib) - ssum
!        f(:,ia) = f(:,ia) + ssum
!     enddo
!     da=0d0 !Force at site ia from PWs ---
!     do  i = nlmto+1, ndimh
!        da(:,:,:) = da(:,:,:) + evecc(i)*reshape(db(i,:,:,:),shape=shape(da),order=[2,1,3]) !db(ilm,:,:)
!     enddo
!     f(:,ia) = f(:,ia) + ewgtt*2d0*[(sum(dconjg(da(:,:,m))*wk(:,:)),m=1,3)] !Force term is (grad psi_kL) * (ppi-evll*sig)*evecc
! !    call tcx('flocbl')
!   end subroutine flocbl
!end module m_rlocbl
