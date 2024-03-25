!>Accumulating rxcq
subroutine x0gpu(rcxq,npr,nwhis,npm) 
  use m_readQG,only:      ngpmx
  use m_genallcf_v3,only: nlmto
  use m_readhbe,only:     nband
  implicit none
  integer:: nwhis,npm,npr                !dim of rcxq
  complex(8):: rcxq(npr,npr,nwhis,npm) !accumulating to rcxq
  ! SetByCPU
  integer:: ngp1, ngp2, ngvecpB1(3,ngpmx),ngvecpB2(3,ngpmx),nadd(3),ncc
  real(8):: symope(3,3),kvec(3)
  complex(8):: cphiq(nlmto,nband), cphim(nlmto,nband)
  complex(8),allocatable:: geigq(:,:),dgeigqk(:,:)
  integer,allocatable:: ngveccR(:,:)
  SetByCPU :block
    use m_x0kf,only: qrk,qq,ispm,ispq
    use m_hamindex,only: symgg=>symops
    use m_readeigen,only: readcphif,readgeigf
    use m_read_ppovl,only: getppx2, ngvecc,ngcread
    use m_readVcoud,only: ngc,ngb 
    use m_itq,only: itq,ntq
    use m_genallcf_v3,only: plat,nctot
    use m_read_bzdata,only: qlat
    use m_readQG,only: ngpmx,ngcmx,Readqg
    real(8)::qk(3),rkvec(3),qdiff(3),qkt(3),qrkt(3)
    kvec =  qq
    qk   =  qrk - qq ! qk = rk. rk is inside 1st BZ, not restricted to the irreducible BZ
    ncc=merge(0,nctot,npm==1)
    associate(cphitemp=> readcphif(qrk,ispq))    
      cphiq(1:nlmto,1:ntq) = cphitemp(1:nlmto,itq(1:ntq)) 
    endassociate
    cphim = readcphif(qk, ispm)
    symope= reshape([1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0],shape=[3,3]) !symgg(:,:,irot) !irot=1 should be E=[1,0,0, 0,1,0, 0,0,1]
    if(ngc/=0) then
      call readqg('QGpsi',qrk, qrkt, ngp1, ngvecpB1) ! qrk= q+rk is mapped to qrkt in BZ
      call readqg('QGpsi',qk,   qkt, ngp2, ngvecpB2) ! qk = rk
      qdiff = matmul(symope,kvec)-qrkt+qkt ! q - qrkt + qkt can be not zero. <M(q) Phi(qk) |Phi(qrk)>
      nadd = nint(matmul(transpose(plat), qdiff)) !nadd: difference in the unit of reciprocal lattice vectors.
      call getppx2(qlat,kvec) ! initialize m_read_ppovl
      if(ngc/=ngcread) call rx( 'melpln2t: ngc/= ngcx by getppx:PPOVLG')
      allocate(ngveccR(1:3,1:ngc),source=ngvecc) ! call rotgvec(symope, 1, ngc, [ngc], qlat, ngvecc, ngveccR)
    endif
    allocate(geigq(ngpmx,nband),dgeigqk(ngpmx,nband))
    geigq   = readgeigf(qrk, ispq) ! IPW part at qrk !G1 for ngp1
    dgeigqk = readgeigf(qk,ispm)   ! IPW part at qk  !G2 for ngp2
    dgeigqk = dconjg(dgeigqk)
  end block SetByCPU
  !We have to send all variables SetByCPU. rcxq should be kept in GPU during kloop:do 1500 in x0kf_v4hz
  GPUblock:block
    use m_x0kf,only: ns1,ns2,ispm,ispq,nqini,nqmax,icounkmink,icounkmaxk
    use m_x0kf,only: iwini,iwend,itc,itpc,jpmc,icouini,whwc
    use m_zmel,only: nbb,miat,tiat,shtvg,ppbir,ppovlz
    use m_genallcf_v3,only: nctot,nclass,natom, nlmto, iclass,nlnmv,nlnmc,icore,ncore,ecore,nlnmx
    use m_readVcoud,only: ngc,ngb 
    use m_rdpp,only: mdimx,nbloch,nblocha 
    use m_read_bzdata,only: nqbz,nqibz, qlat,ginv,qbz,qibz,wbz
    use m_read_ppovl,only: igggi,igcgp2i,nxi,nxe,nyi,nye,nzi,nze,nvgcgp2,ngcgp,ggg,ppovlinv
    use m_itq,only: itq,ntq
    complex(8),parameter:: img=(0d0,1d0),tpi= 8d0*datan(1d0)
    integer:: isp,ncnv,ncorec,nccc,mdim,it,ia,irot=1
    integer:: i,iap,ias,ib,ic,icp,nc,nc1,nv,ics,itp,iae,ims,ime
    integer:: icoun,nmini,nmmax,nmtot,nqtot,nt0,ntp0,invr,igb1,igb2,iw
    integer:: iasx(natom),icsx(natom),iatomp(natom),imdim(natom)
    real(8):: ppb(nlnmx,nlnmx,mdimx,nclass) ! ppb= <Phi(SLn,r-R)_q,isp1 |Phi(SL'n',r-R)_qk,isp2 B_k(S,i,rot^{-1}(r-R))>
    real(8):: tr(3,natom),shtv(3)
    complex(8),allocatable :: zmel(:,:,:) ! zmel(nbb,nmtot, nqtot), nbb:mixproductbasis, nmtot:middlestate, nq
    logical:: zmelconjg=.true.
    ! Former part for generating zmel is the copy of get_zmel_init in m_zmel.f90
    ppb = ppbir(:,:,:,:,irot,ispq)           !MPB has no spin dependence
    invr  = 1 !invgx(irot=1)       !invrot (irot,invgx,ngrp) ! Rotate atomic positions invrot*R = R' + T
    tr    = tiat(:,:,invr)
    iatomp= miat(:,invr)
    shtv  = shtvg(:,invr) !matmul(symope,shtvg(:,invr))
    imdim = [( sum(nblocha(iclass(1:ia-1)))+1  ,ia=1,natom)]
    iasx=[(sum(nlnmv(iclass(1:ia-1)))+1,ia=1,natom)]
    icsx=[(sum(ncore(iclass(1:ia-1))),ia=1,natom)]
    nmini = ns1          !starting index of middle state  (nctot+nvalence order)
    nmmax = ns2-nctot    !end      index of middle state  (nctot+nvalence order)
    nt0= nmmax-nmini+1
    ntp0=nqmax-nqini+1
    nmtot  = nctot + nt0     ! = phi_middle nmtot=ns2-ns1+1
    nqtot  = ncc   + ntp0    ! = phi_end
    ZmelBlock:block
      complex(8):: zmelt(1:nbloch+ngc,nmtot,nqtot)
      zmelt=0d0
      ZmelWithinMT: block !- Calculates <psi_q(itp) |psi_qk(it) B_k(rot(r-R))> 
        complex(8):: phasea(natom) 
        phasea = [(exp(-img *tpi* sum(kvec*tr(:,ia))),ia=1,natom)] 
        iatomloop: do concurrent(ia = 1:natom)
          ic    = iclass(ia)
          nc    = nlnmc(ic) !nlnmc      = number of l,n,m for core states
          nv    = nlnmv(ic) !nlnmv      = number of l,n,m for valence
          nc1   = nc + 1
          ncnv  = nc+nv  !ncore + nvalence
          iap   = iatomp(ia)  
          icp   = iclass(iap)
          ias   = iasx(ia)  !start of nlnmv for ia
          iae   = ias+nv-1  !end  
          ics   = icsx(ia)   
          ims   = imdim(iap) !start of PB for ia
          ime   = ims-1+nblocha(icp)
          mdim  = nblocha(icp)
          ncorec= merge(ncore(ic),0,nctot>0)
          nccc  = merge(ncore(ic),0,ncc>0)
          associate(      &
               dcphiqk => dconjg(cphim(ias:iae,nmini:nmmax)),&
               tdcphiqk=> transpose(dconjg(cphim(ias:iae,nmini:nmmax))),&
               cphiq   => cphiq(ias:iae,nqini:nqmax),&
               ppbc    => ppb(nc1:ncnv,icore(1:ncorec,ic),:,icp)) !all readin but we need only bands nmini: or nqini:
            do concurrent(i=1:mdim) !valence-valence  !=== this may be time-consuming block ==================
              zmelt(i-1+ims,nctot+1:nctot+nt0,ncc+1:ncc+ntp0)=phasea(ia) &
                   *matmul(tdcphiqk,matmul(transpose(ppb(nc1:ncnv,nc1:ncnv,i,icp)),cphiq))
            enddo
            do concurrent(it=1:ncorec) !core-valence
              zmelt(ims:ime, ics+it, ncc+1:ncc+ntp0)=phasea(ia)* matmul(transpose(ppbc(:,it,1:mdim)),cphiq(:,1:ntp0))
            enddo
            do concurrent(itp=1:nccc) !valence-core
              zmelt(ims:ime, nctot+1:nctot+nt0, ics+itp)=phasea(ia)*matmul(transpose(ppbc(:,itp,1:mdim)),dcphiqk(:,1:nt0))
            enddo                              !^^^^^^^^^phasea(ia) right? 2024-1-7 or phasea(ia) or dcongj(phasea(ia))?
          endassociate
        enddo iatomloop
      endblock ZmelWithinMT
      if(ngc/=0)then
        ZmelIPW:block  !> Mattrix elements <Plane psi |psi> from interstitial plane wave.
          integer:: igcgp2,nn(3),iggg,igp1,itp,igc,igp2,igcgp2i_(ngc,ngp2)
          complex(8):: zmelp0(ngc,nt0,ntp0),phase(ngc) , ggitp(ngcgp,ntp0),gggmat(ngcgp,ngp1)
          do concurrent(igcgp2=1:ngcgp,igp1=1:ngp1) !G synthesized 
            nn = ngvecpB1(:,igp1)- nvgcgp2(:,igcgp2) - nadd ! G1 -(Gc+G2) - Gadd.  Note that -Gadd= -rk + qt -qkt
            if(nn(1)<nxi .OR. nxe<nn(1) .OR. nn(2)<nyi .OR. nye<nn(2) .OR. nn(3)<nzi .OR. nze<nn(3)) cycle
            iggg = igggi(nn(1),nn(2),nn(3)) !inversion table
            if(iggg>=0) gggmat(igcgp2,igp1)=ggg(iggg)
          enddo
          do concurrent (igc=1:ngc,igp2=1:ngp2) !igc for B !G synthesized 
            nn = ngveccR(:,igc) + ngvecpB2(:,igp2)
            igcgp2i_(igc,igp2)=igcgp2i(nn(1),nn(2),nn(3))
          enddo
          phase(:)=[(exp( -img*tpi*sum((matmul(symope,kvec)+matmul(qlat,ngveccR(:,igc)))*shtv) ),igc=1,ngc)]
          !          associate(&
          !               geigq   => readgeigf(qrk, ispq),&         !read IPW part at q   !G1 for ngp1
          !               dgeigqk => dconjg(readgeigf(qk,ispqk))) !read IPW part at qk  !G2 for ngp2
          ggitp(:,:)= matmul(gggmat,geigq(1:ngp1,itq(nqini:nqmax)))
          do concurrent (itp= 1:ntp0) !=== this may be time-consuming block (or maynot)==================
            associate( ggitp_=>reshape([((ggitp(igcgp2i_(igc,igp2),itp),igc=1,ngc),igp2=1,ngp2)],shape=[ngc,ngp2]))
              zmelp0(:,:,itp)= matmul(ggitp_,dgeigqk(1:ngp2,nmini:nmmax))
            endassociate
          enddo
          forall(igc=1:ngc) zmelp0(igc,:,:)=phase(igc)*zmelp0(igc,:,:) 
          !          endassociate
          call matm(ppovlinv,zmelp0,zmelt(nbloch+1:nbloch+ngc,nctot+1:nctot+nt0,ncc+1:ncc+ntp0),ngc,ngc,ntp0*nt0)
        endblock ZmelIPW
        deallocate(geigq,dgeigqk)
      endif
      allocate(zmel(nbb,ns1:ns2, nqtot))
      call matm(dconjg(transpose(ppovlz)), zmelt, zmel,nbb,ngb,nmtot*nqtot) !MultiplePPOVLZ
      if(zmelconjg) zmel=dconjg(zmel)
    endblock ZmelBlock
    icounloop: do 1000 icoun=icounkmink,icounkmaxk
      TimeConsumingRcxq: block 
        complex(8):: zmelzmel(npr,npr) 
        associate( &
             jpm => jpmc(icoun),&  !\pm omega 
             it  => itc (icoun),&  !occ      at k
             itp => itpc(icoun))   !unocc    at q+k
          do concurrent(igb1=1:npr,igb2=1:npr) 
            zmelzmel(igb1,igb2)= dconjg(zmel(igb1,it,itp))*zmel(igb2,it,itp) 
          enddo
          do iw=iwini(icoun),iwend(icoun) !)& !rcxq is hermitian, thus, we can reduce computational time half.
            rcxq(:,:,iw,jpm)=rcxq(:,:,iw,jpm)+ whwc(iw-iwini(icoun)+icouini(icoun))* zmelzmel(:,:) ! Use zaxpy and symmetrize?
          enddo
          !forall(iw=iwini(icoun):iwend(icoun)) nwj(iw,jpm)=nwj(iw,jpm)+iwend(icoun)-iwini(icoun)+1 !counter check
        endassociate
      endblock TimeConsumingRcxq
1000 enddo icounloop
    deallocate(zmel)
  endblock GPUblock
end subroutine x0gpu
