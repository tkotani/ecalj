module m_hambl
  public hambl
contains
  subroutine hambl(isp,qin,smpot,vconst,osig,otau,oppi, h,s)! Make LDA/GGA Hamiltonian and overlap matrix for a k-point. No SOC added.
    use m_lmfinit,only: nbas , sspec=>v_sspec,nsp
    use m_igv2x,only: napw, igvapwin=>igv2x, ndimh
    use m_supot,only: n1,n2,n3
    use m_struc_def,only: s_rv1,s_cv1,s_rv4,s_cv5
    use m_lattic,only:plat=>lat_plat,qlat=>lat_qlat
    use m_augmbl,only: augmbl
    use m_hsibl,only:hsibl
    !i   isp   :spin index
    !i   qin    :Bloch vector (k-point).
    !i   k1,k2,k3 dimensions of smpot
    !i   smpot :smooth potential on uniform mesh (mkpot.f)
    !i   vconst: additional constant potential
    !i   osig,otau,oppi:  augmentation matrices
    !i   ndimh :dimension of hamiltonian and ovrlap h,s
    !o Outputs
    !o   h     :Hamiltonian matrix Bloch sum of T_ij + V_ij.
    !    See Eq.(C.2) and (C.3) in JPSJ84,034702(2015)
    !    ecalj/Document/PAPERandPRESENTATION/KotaniKinoAkai2015FormulationPMT.pdf [1]
    !o   s     :overlap matrix
    !    See Eq.(C.1)           in [1]
    !r  qpg(ig) = tpiba * ( qin + matmul(qlat,igapwin(1:3,ig))), ig=1,napw
    implicit none
    integer:: mode,isp,i
    type(s_cv5) :: oppi(3,nbas)
    type(s_rv4) :: otau(3,nbas)
    type(s_rv4) :: osig(3,nbas)
    real(8):: vconst, qin(3)
    complex(8):: smpot(n1,n2,n3,nsp), h(ndimh,ndimh),s(ndimh,ndimh)
    call tcn('hambl')
    h = 0d0 !Hamiltonian for the basis of MTO+APW
    s = 0d0 !Overlap matrix for the basis of MTO+APW
    call augmbl(isp,qin,osig,otau,oppi,ndimh, h,s)! Augmentation parts of h,s
    !    product sum f structure constant C_akL^i in Eq.(C.1)-(C.2) in Ref.[1].
    call smhsbl(vconst,qin,ndimh,napw,igvapwin,          h,s)!Smooth and Constant potential parts.
    call hsibl(n1,n2,n3,smpot,isp,qin,ndimh,napw,igvapwin, h)!Smooth potential part, 1st term of (C.3) in Ref.[1]
    do i=1,ndimh
       h(i+1:ndimh,i)=dconjg(h(i,i+1:ndimh)) !RU part = LD part*
       s(i+1:ndimh,i)=dconjg(s(i,i+1:ndimh))
    enddo
    call tcx('hambl')
  endsubroutine hambl
  subroutine smhsbl(vavg,q,ndimh, napw,igapw, h,s)!- Smoothed Bloch Hamiltonian (constant potential) and overlap matrix
    use m_lmfinit,only: alat=>lat_alat,nbas,nkaphh,lhh, ispec,sspec=>v_sspec,lmxa_i=>lmxa
    use m_lattic,only: lat_plat,rv_a_opos,qlat=>lat_qlat,vol=>lat_vol
    use m_uspecb,only:uspecb
    use m_orbl,only: Orblib1,Orblib2,ktab1,ltab1,offl1,norb1,ktab2,ltab2,offl2,norb2
    use m_ropyln,only: ropyln
    use m_smhankel,only:hhibl
    implicit none
    intent(in)::    vavg,q,ndimh, napw,igapw 
    !i Inputs
    !i   vavg  :constant potential (MT zero) to be added to h
    !i   q     : q wave vector
    !i   ndimh : dimension of hamiltonian
    !o Outputs
    !o   h     :smooth Bloch hamiltonian added to h (= s * vavg)
    !o   s     :smooth Bloch overlap added to s
    !b Bugs
    !b   Use of qpgv(1,1:napw) is poor design ... rewrite
    !r Remarks
    !r  *How orbital information is extracted and deployed.
    !r   Orbital specification requires the following information:
    !r     1. orbital type, specified by the triplet (l,rsmh,eh,pz)
    !r        (pz is relevant only for local orbitals)
    !r        Information is stored in spec->orbp, unpacked by uspecb
    !r        For local orbitals spec->pz is also required.
    !r        Note that each orbital type has 2*l+1 orbitals.
    !r     3. Location of orbitals (offsets) in hamiltonian
    !r
    !r   Routines needing this information unpack the data as follows:
    !r
    !r      call    Purpose, Data and format
    !r     (Input)
    !r     -------  ---------------------------------
    !r
    !r     uspecb    Extract orbital parameters for all orbital types.
    !r               rsmh(l,ik),eh(l,ik), l=0,lh(ik), ik=1..nkapi.
    !r               Entries for which rsmh(l,ik)>0 have envelopes.
    !r               Entries for which ik=nkapi and pz(l)>0 are local orbitals.
    !r               Note that these possiblities can simultaneously occur.
    !r
    !r     orbl      norb : (total number of orbital types)
    !r
    !r     gtbsl    Serves two purposes.
    !r     (norb,    1. Blocks orbitals with common (rsmh,eh) and
    !r     ltab,     consecutive l so that mesh tabulations are more
    !r     ktab,     efficiently generated.
    !r     rsmh,eh)  2. Marks each orbital type as being :
    !r                  a. valence function with envelope;
    !r                  b. local orbital with envelope;
    !r                  c. local orbital without envelope.
    !r               blks(io) = size of contiguous block of orbitals.  Note
    !r               that if io+1, io+2, ... are contiguous to io,
    !r               blks(io+1)=blks(io+2)=...=0
    !r
    !r
    !r  *The Bloch sum of s(q) is computed as s = sum_T s_T exp(iq.T)
    !r   where T are the set of lattice vectors.
    !r
    !r  *NB: The matrix subblock in s for sites (i,j) is computed with
    !r       the connecting vector pi-pj.  This convention is different
    !r       from the tight-binding convention, where s is computed for
    !r       pj-pi.  Thus this (i,j) matrix subblock is the hermitian
    !r       conjugate of the same subblock in the tight-binding case.
    !r       tight-binding case.  Note, however, that the entire s
    !r       would not correspond to the hermitian conjugate of the tight
    !r       binding case.
    !r
    !u Updates see github log after 2009
    !u   05 Jul 08 (T. Kotani) output density for new PW part
    !u   12 Aug 04 First implementation of extended local orbitals
    !u   14 Aug 02 Added overlap-only option
    !u   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
    !u   15 Feb 02 (ATP) Added MPI parallelization
    !u   27 Aug 01 Extended to local orbitals.
    !u   02 Mar 01 Bug fix for multiple-kappa case
    !u   19 May 00 Adapted from nfp smhs_q.f
    ! ----------------------------------------------------------------------
    integer :: nlms,kdim,n0,nkap0
    parameter (nlms=25, kdim=1, n0=10, nkap0=3)
    integer :: procid,master, mode,ndimh,napw,igapw(3,napw)
    real(8):: q(3) , vavg
    complex(8):: h(ndimh,ndimh),s(ndimh,ndimh)
    integer :: nlmto, i1,i2,ib1,ib2,ilm1,ilm2,io1,io2,is1,is2,nlm1,nlm2,l1,l2,ig,&
         lmxax,lmxa,nlmax, lh1(nkap0),lh2(nkap0),nkap1,nkap2, &
         ik1,blks1(n0*nkap0),ntab1(n0*nkap0), ik2,blks2(n0*nkap0),ntab2(n0*nkap0)
    integer:: iloop,iloopmx
    real(8) :: e1(n0,nkap0),rsm1(n0,nkap0),p1(3),e2(n0,nkap0),rsm2(n0,nkap0),p2(3),xx
    complex(8):: s0(nlms,nlms,0:kdim,nkap0,nkap0)
    real(8) :: qpg2,srvol,tpiba,denom,gam,ddot
    real(8),allocatable:: yl(:),ylv(:,:),qpgv(:,:),qpg2v(:)
    complex(8),allocatable:: srm1l(:)
    complex(8):: ovl,srm1,phase,fach
    real(8),parameter::pi = 4d0*datan(1d0),fpi = 4*pi
    parameter (srm1=(0d0,1d0))
    call tcn('smhsbl')
    procid = 0
    master = 0
    nlmto = ndimh-napw
    tpiba = 2d0*pi/alat
    srvol = dsqrt(vol)
    if(napw > 0) then
       lmxax = -1
       do  ib1 = 1, nbas 
          is1=ispec(ib1) 
          lmxa=lmxa_i(is1)
          lmxax = max(lmxax,lmxa)
       enddo
       nlmax=(lmxax+1)**2
       allocate(ylv(napw,nlmax),yl(nlmax),qpgv(3,napw),qpg2v(napw))
       do  ig = 1, napw
          qpgv(:,ig) = tpiba * (q + matmul(qlat, igapw(:,ig)))
       enddo
       call ropyln(napw,qpgv(1,1:napw),qpgv(2,1:napw),qpgv(3,1:napw), lmxax,napw,ylv,qpg2v)
       allocate(srm1l(0:lmxax))
       srm1l(0) = 1d0
       do  l1 = 1, lmxax
          srm1l(l1) = (srm1)**l1
       enddo
    endif
    if(nlmto >0 ) then
       ib1loop: do 1010 ib1=1,nbas
          is1=ispec(ib1) 
          p1 =rv_a_opos(:,ib1)   !site at ib1
          call uspecb(is1,rsm1,e1) 
          call orblib1(ib1) !norb1,ltab1,ktab1,offl1
          call gtbsl8(norb1,ltab1,ktab1,rsm1,e1,ntab1,blks1)
          ib2loop: do  ib2 = ib1, nbas !Hsm x Hsm part
             is2=ispec(ib2) 
             p2=rv_a_opos(:,ib2) !site at ib2
             call uspecb(is2,rsm2,e2)
             call orblib2(ib2) !norb2,ltab2,ktab2,xx,offl2,xx)
             call gtbsl8(norb2,ltab2,ktab2,rsm2,e2,ntab2,blks2)
             !     ... M.E. <1> and <T> between all envelopes connecting ib1 and ib2
             do  i1 = 1, nkaphh(is1) !nkap1 for ib1
                do  i2 = 1, nkaphh(is2) !nkap2 for ib2
                   nlm1 = (lhh(i1,is1)+1)**2
                   nlm2 = (lhh(i2,is2)+1)**2
                   if (nlm1 > nlms .OR. nlm2 > nlms) call rx('smhsbl: increase nlms')
                   call hhibl(p1,p2,q,rsm1(1,i1),rsm2(1,i2),e1(1,i1),e2(1,i2),nlm1,nlm2,1,nlms,nlms,s0(1,1,0,i1,i2)) 
                   !F0*F0 and F0*\laplacian F0 integrals (smH parts)
                enddo
             enddo
             do  io2 = 1, norb2 ! io2 for ib2 for orbital indices, poke block of integrals into s,h
                if (blks2(io2) == 0) cycle
                l2  = ltab2(io2) ! l2,ik2 = l and kaph indices, needed to locate block in s0
                ik2 = ktab2(io2)
                i2 = offl2(io2)
                do  ilm2 = l2**2+1, (l2+1)**2
                   i2 = i2+1
                   do  io1 = 1, norb1 !io1 for ib1
                      if (blks1(io1)==0) cycle
                      l1  = ltab1(io1)! l1,ik1 = l and kaph indices, needed to locate block in s0
                      ik1 = ktab1(io1)
                      i1 = offl1(io1)
                      do  ilm1 = l1**2+1, (l1+1)**2
                         i1 = i1+1
                         s(i1,i2)= s(i1,i2) + s0(ilm1,ilm2,0,ik1,ik2)
                         h(i1,i2)= h(i1,i2) - s0(ilm1,ilm2,1,ik1,ik2) +vavg*s0(ilm1,ilm2,0,ik1,ik2)
                         !                                 1:kinetic                    !0: constant
                      enddo
                   enddo
                enddo
             enddo
          enddo ib2loop
          ! Hsm x PW part  Hsm(i1) x PW (i2) Similar logic in fsmbl
          igloop: do  ig = 1, napw
             i2 = ig + nlmto
             qpg2 = qpg2v(ig)
             phase = exp(srm1*alat*ddot(3,qpgv(1,ig),1,p1,1))
             do  io1 = 1, norb1
                if (blks1(io1) == 0) cycle
                l1  = ltab1(io1)
                ik1 = ktab1(io1)
                i1  = offl1(io1)
                denom = e1(l1+1,ik1) - qpg2
                gam   = 1d0/4d0*rsm1(l1+1,ik1)**2
                fach  = -fpi/denom * phase * srm1l(l1) * exp(gam*denom)
                do  ilm1 = l1**2+1, (l1+1)**2
                   i1 = i1+1
                   !     <Hsm^bloch | exp(i q+G r)/srvol > and
                   !     <Hsm^bloch | nabla + vavg |  exp(i q+G r)/srvol > in a cell.
                   ovl = fach * ylv(ig,ilm1)/srvol ! JMP Eq.(9.4)
                   s(i1,i2) = s(i1,i2) + ovl
                   h(i1,i2) = h(i1,i2) + qpg2*ovl + vavg*ovl
                enddo
             enddo
          enddo igloop
1010   enddo ib1loop
    endif
    do  ig = 1, napw !! ... PW x PW part (diagonal matrix)
       i2 = ig + nlmto
       s(i2,i2) = s(i2,i2) + 1d0
       h(i2,i2) = h(i2,i2) + qpg2v(ig) + vavg
    enddo
    if (napw > 0)deallocate(yl,ylv,qpgv,qpg2v,srm1l)
    do  i1 = 1, ndimh !! ... Occupy second half of matrix
       do  i2 = i1, ndimh
          h(i2,i1) = dconjg(h(i1,i2))
          s(i2,i1) = dconjg(s(i1,i2))
       enddo
    enddo
    call tcx('smhsbl')
  end subroutine smhsbl
endmodule m_hambl
