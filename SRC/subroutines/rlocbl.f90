module m_rlocbl
  use m_ll,only:ll
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
    use m_struc_def,only: s_rv1,s_cv1,s_rv5,s_rv4,s_cv5
    use m_lmfinit,only: alat=>lat_alat,ispec
    use m_lmfinit,only: lmxa_i=>lmxa,lmxb_i=>lmxb,kmxt_i=>kmxt
    use m_lattic,only: qlat=>lat_qlat,rv_a_opos
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
    implicit none
    integer :: lfrce,nbas,isp,ndimh,nspc,nevec,lekkl, ndham,napw,igvapw(3,napw)
    integer :: is,nlmx,ktop0,npmx, kmaxx,nlmax,nlmto, ia,isa,ivec,kmax,lmxa,nlma,lmxha,nlmha,ispc,ksp
    type(s_cv5),target:: sv_p_oppi(3,nbas)
    type(s_rv4),target:: sv_p_otau(3,nbas)
    type(s_rv4),target:: sv_p_osig(3,nbas)
    type(s_rv5),target:: sv_p_oeqkkl(3,nbas), sv_p_oqkkl(3,nbas)
    real(8),pointer:: qpp(:,:,:,:,:),eqpp(:,:,:,:,:),qhp(:,:,:,:,:),qhh(:,:,:,:,:),eqhp(:,:,:,:,:),eqhh(:,:,:,:,:)
    real(8):: q(3),f(3,nbas)
    real(8),target:: ewgt(nevec),evl(ndham,isp)
    complex(8),target:: evec(ndimh,nspc,nevec)
    real(8):: force(3,nbas)
    real(8),   pointer:: sigpp(:,:,:,:),   sighp(:,:,:,:) 
    complex(8),pointer:: ppihpz(:,:,:,:,:),ppippz(:,:,:,:,:)
    if (nevec <= 0) return
    call tcn('rlocbl')
    nlmax = 0
    kmaxx = 0
    do  ia = 1, nbas
       isa = ispec(ia) 
       lmxa= lmxa_i(isa)
       kmax= kmxt_i(isa)
       nlma = (lmxa+1)**2
       kmaxx = max(kmaxx,kmax)
       nlmax = max(nlmax,nlma)
    enddo
    nlmto = ndimh - napw
    nlmx  = nlmax
    ktop0 = kmaxx
    force=0d0
    if(nlmax > nlmx)  call rxi('rlocbl: nlmx < nlma=',nlmax)
    if(kmaxx > ktop0) call rxi('rlocbl: ktop0 < kmax=',kmax)
    ialoop: do ia = 1,nbas !Loop over augmentation sites 
       isa = ispec(ia) 
       lmxa= lmxa_i(isa)
       if (lmxa == -1) cycle
       lmxha=lmxb_i(isa)
       kmax= kmxt_i(isa)
       nlmha = (lmxha+1)**2
       nlma  = (lmxa+1)**2
       qpp  => sv_p_oqkkl(1,ia)%v
       qhp  => sv_p_oqkkl(2,ia)%v
       qhh  => sv_p_oqkkl(3,ia)%v
       eqpp => sv_p_oeqkkl(1,ia)%v
       eqhp => sv_p_oeqkkl(2,ia)%v
       eqhh => sv_p_oeqkkl(3,ia)%v
       sigpp=> sv_p_osig(1,ia)%v  !sigpp(0:kmax,0:kmax,0:lmxa,nsp),sighp(nkaph,0:kmax,0:lmxha,nsp)
       sighp=> sv_p_osig(2,ia)%v 
       call bstrux_set(ia,q)  !Get bstr dbstr structure constant 
       ! In noncollinear case, isp=1 always => need internal ispc=1..2
       !   ksp is the current spin index in both cases: ksp = isp in the collinear case;  ksp= ispc in the noncollinear case,
       !   whereas ispc is spin index in the noncoll case, but 1 for coll.
       cPklblock: block
         use m_orbl,only: Orblib, norb,ltab,ktab,offl, ntab,blks
         use m_lmfinit,only:nsp
         complex(8):: cPkl(0:kmax,nlma),wtt
         complex(8),pointer::evecc(:)
         real(8),pointer:: ewgtt,evll
         integer:: k, ilma,ilm1,ilm2,k1,k2,io1,l1,ik1,i1,io2,l2,ik2,i2
         ppippz=>sv_p_oppi(1,ia)%cv !ppihpz(nkaph,0:kmax,nlmha,nlma,nsp),  ppippz(0:kmax,0:kmax,nlma,nlma,nsp)
         ppihpz=>sv_p_oppi(2,ia)%cv
         ivecloop: do ivec = 1, nevec
            ispcloop: do ispc = 1, nspc
               ksp = max(ispc,isp) !for lso=1, we use isp=1 only. thus ksp=ispc, for lso/=1, ksp=isp since nspc=1
               evecc=> evec(1:ndimh,ispc,ivec)
               ewgtt=> ewgt(ivec)
               evll => evl(ivec,isp)
               ! do  k = 0, kmax
               !    cPkL(k,:) =  matmul(evec(1:ndimh,ispc,ivec),bstr(1:ndimh,:,k)) ! Pkl expansion of eigenvector
               ! enddo
               ! MO replaced the above loop with the following 2024-11-07
               block
                 use m_blas, only: zmv => zmv_h, m_op_T
                 integer :: istat
                 complex(8) :: cPkLT(nlma, 0:kmax)
                 istat = zmv(bstr, evec(1,ispc,ivec), cPkLT, m=ndimh, n=nlma*(kmax+1), opA=m_op_T)
                 cPkL = transpose(cPkLT)
               endblock
               do  ilm2 = 1, nlma !Add to local density coefficients for one state
                  do  k2 = 0, kmax
                     qpp (:,k2,:,ilm2,ksp)= qpp(:,k2,:,ilm2,ksp) + ewgtt*dconjg(cPkL(:,:))*cPkL(k2,ilm2)
                     eqpp(:,k2,:,ilm2,ksp)= eqpp(:,k2,:,ilm2,ksp)+ ewgtt*dconjg(cPkL(:,:))*cPkL(k2,ilm2)*evll
                  enddo
               enddo
               call orblib(ia)! norb,ltab,ktab,offl      !     Block into groups of consecutive l
               iorbital1: do io1 = 1, norb 
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
               enddo iorbital1
               if(lfrce==0) cycle ! ... Contribution to forces
               call rxx(nspc.ne.1,'forces not implemented in noncoll case')
               Flocblblock: block ! Force contribution from augmentation at site ia.
                 use m_orbl,only: Orblib, norb,ltab,ktab,offl, ntab,blks
                 integer :: isp,i,ib,ilm,ilmb,io,iq,k,l2,m,nlm1,ik,ilma,jlm,k1,k2,l,oi,iblk 
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
               endblock Flocblblock
            enddo ispcloop
         enddo ivecloop
       endblock cPklblock
    enddo ialoop
    if(lfrce/=0) f = force+f 
    call tcx('rlocbl')
  end subroutine rlocbl
end module m_rlocbl
