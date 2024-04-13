!>Accumulating rxcq
subroutine x0gpu(rcxq,npr,nwhis,npm) 
  use m_readQG,only:      ngpmx
  use m_genallcf_v3,only: nlmto,natom
  use m_readhbe,only:     nband
  use m_rdpp,only:        mdimx
  use m_genallcf_v3,only: nlnmx,nclass
  use m_x0kf,only:        nqini,nqmax,ns1,ns2
  use m_rdpp,only:        nbloch,nblocha
  use m_kind,only: kindrcxq
  implicit none
  intent(in)::        npr,nwhis,npm
  intent(inout)::rcxq
  integer:: nwhis,npm,npr                !dim of rcxq
  complex(8),parameter:: img=(0d0,1d0),tpi= 8d0*datan(1d0)
  complex(kindrcxq):: rcxq(npr,npr,nwhis,npm) !accumulating to rcxq
  ! SetByCPU
  real(8):: symope(3,3),kvec(3)
  integer:: ncc,ngp1, ngp2, ngvecpB1(3,ngpmx),ngvecpB2(3,ngpmx),nadd(3)
  complex(8):: cphiq(nlmto,nband), cphim(nlmto,nband)
  complex(8),allocatable:: geigq(:,:),dgeigqk(:,:)
  integer,allocatable:: ngveccR(:,:)
  real(8):: ppb(nlnmx,nlnmx,mdimx,nclass) ! ppb= <Phi(SLn,r-R)_q,isp1 |Phi(SL'n',r-R)_qk,isp2 B_k(S,i,rot^{-1}(r-R))>
  real(8):: tr(3,natom),shtv(3)

  logical :: tgpu = .true.
#ifdef __GPU
  character(len=3) :: acway = 'gpu'
#else
  character(len=3) :: acway = 'cpu'
#endif
  real(8) :: t1, t2

  complex(8):: phasea(natom) 
  integer:: iasx(natom),icsx(natom),imdim(natom),iatomp(natom),nmini,nmmax,nmtot,nqtot,nt0,ntp0
  call cpu_time(t1)
  SetByCPU :block
    use m_zmel,only: miat,tiat,shtvg,ppbir
    use m_x0kf,only: qrk,qq,ispm,ispq
    use m_hamindex,only: symgg=>symops
    use m_readeigen,only: readcphif,readgeigf
    use m_read_ppovl,only: getppx2, ngvecc,ngcread
    use m_readVcoud,only: ngc,ngb 
    use m_itq,only: itq,ntq
    use m_genallcf_v3,only: plat,nctot,iclass,nlnmv,ncore
    use m_read_bzdata,only: qlat
    use m_readQG,only: ngpmx,ngcmx,Readqg
    integer::ia,invr
    real(8)::qk(3),rkvec(3),qdiff(3),qkt(3),qrkt(3)
    symope= reshape([1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0],shape=[3,3]) !symgg(:,:,irot) !irot=1 should be E=[1,0,0, 0,1,0, 0,0,1]
    kvec =  qq
    qk   =  qrk - qq ! qk = rk. rk is inside 1st BZ, not restricted to the irreducible BZ
    ncc=merge(0,nctot,npm==1)
    associate(cphitemp=> readcphif(qrk,ispq))    
      cphiq(1:nlmto,1:ntq) = cphitemp(1:nlmto,itq(1:ntq)) 
    endassociate
    cphim = readcphif(qk, ispm)
    if(ngc/=0) then
      call readqg('QGpsi',qrk, qrkt, ngp1, ngvecpB1) ! qrk= q+rk is mapped to qrkt in BZ
      call readqg('QGpsi',qk,   qkt, ngp2, ngvecpB2) ! qk = rk
      qdiff = kvec-qrkt+qkt ! q - qrkt + qkt can be not zero. <M(q) Phi(qk) |Phi(qrk)>
      nadd = nint(matmul(transpose(plat), qdiff)) !nadd: difference in the unit of reciprocal lattice vectors.
      call getppx2(qlat,kvec) ! initialize m_read_ppovl
      if(ngc/=ngcread) call rx( 'melpln2t: ngc/= ngcx by getppx:PPOVLG')
      allocate(ngveccR(1:3,1:ngc),source=ngvecc) ! call rotgvec(symope, 1, ngc, [ngc], qlat, ngvecc, ngveccR)
    endif
    allocate(geigq(ngpmx,nband),dgeigqk(ngpmx,nband))
    geigq   = readgeigf(qrk, ispq) ! IPW part at qrk !G1 for ngp1
    dgeigqk = readgeigf(qk,ispm)   ! IPW part at qk  !G2 for ngp2
    dgeigqk = dconjg(dgeigqk)
    ppb = ppbir(:,:,:,:,1,ispq)           !MPB has no spin dependence irot=1
    invr  = 1 !invgx(irot=1)       !invrot (irot,invgx,ngrp) ! Rotate atomic positions invrot*R = R' + T
    tr    = tiat(:,:,invr)
    iatomp= miat(:,invr)
    imdim = [( sum(nblocha(iclass(1:ia-1)))+1  ,ia=1,natom)]
    iasx=[(sum(nlnmv(iclass(1:ia-1)))+1,ia=1,natom)]
    icsx=[(sum(ncore(iclass(1:ia-1))),ia=1,natom)]
    shtv  = shtvg(:,invr) !matmul(symope,shtvg(:,invr))
    nmini = ns1          !starting index of middle state  (nctot+nvalence order)
    nmmax = ns2-nctot    !end      index of middle state  (nctot+nvalence order)
    nt0= nmmax-nmini+1
    ntp0=nqmax-nqini+1
    nmtot  = nctot + nt0     ! = phi_middle nmtot=ns2-ns1+1
    nqtot  = ncc   + ntp0    ! = phi_end
    phasea = [(exp(-img *tpi* sum(kvec*tr(:,ia))),ia=1,natom)] 
  end block SetByCPU
  call cpu_time(t2)
  print '(1x,A,2F10.6)','x0setcpu' ,(t2-t1)
  !We have to send all variables SetByCPU. rcxq should be kept in GPU during kloop:do 1500 in x0kf_v4hz
  GPUblock:block
    use m_zmel,only: ppovlz,nbb
    use m_x0kf,only: icounkmink,icounkmaxk, iwini,iwend,itc,itpc,jpmc,icouini,whwc
    use m_genallcf_v3,only: nctot,natom, nlmto, iclass,nlnmv,nlnmc,icore,ncore,ecore,nlnmx
    use m_read_ppovl,only: igggi,igcgp2i,nxi,nxe,nyi,nye,nzi,nze,nvgcgp2,ngcgp,ggg,ppovlinv
    use m_read_bzdata,only: qlat
    use m_readVcoud,only: ngc,ngb 
    use m_itq,only: itq,ntq
#ifdef __GPU
    use openacc
    use cudafor
    use m_math_gpu, only: mm_op_c, mm_op_n, zmm, zmm_sb
#else
    use m_math_cpu, only: mm_op_c, mm_op_n, zmm, zmm_sb
#endif
    integer :: ierr
    complex(8), allocatable :: zmel(:,:,:) ! zmel(nbb,nmtot,nqtot), nbb:mixproductbasis, nmtot:middlestate, nq
    complex(8), allocatable :: zmeltvv(:,:,:)
#ifdef __GPU
    attributes(device) :: zmel
    attributes(device) :: zmeltvv
#endif
    integer:: icoun,igb1,igb2,iw,ia
    logical,parameter   :: zmelconjg=.true.

    call cpu_time(t1)
    ZmelBlock: block ! ZmelBlock is a part of copy of get_zmel_init in m_zmel.f90
      complex(8):: zmelt(1:nbloch+ngc,nmtot,nqtot)
      complex(8) :: zmelt_d(1:nbloch+ngc,nmtot,nqtot)
#ifdef __GPU
      attributes(device) :: zmelt_d
#endif
      zmelt=0d0
      ZmelWithinMT: block !- Calculates <psi_q(itp) |psi_qk(it) B_k(rot(r-R))> 
        integer:: i,iap,ias,ib,ic,icp,nc,nc1,nv,ics,itp,iae,ims,ime, ncnv,ncorec,nccc,mdim,it
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

      ! print *, 'nbloch, ngc, nctot, ncc, ntp0, nt0:', nbloch, ngc, nctot, ncc, ntp0, nt0
      zmelt_d = zmelt
      if(ngc/=0)then
        ZmelIPW:block  !> Mattrix elements <Plane psi |psi> from interstitial plane wave.
          integer:: igcgp2,nn(3),iggg,igp1,itp,igc,igp2,igcgp2i_(ngc,ngp2), it
          complex(8):: zmelp0(ngc,nt0,ntp0),phase(ngc) , ggitp(ngcgp,ntp0),gggmat(ngcgp,ngp1)
          complex(8):: ggitp_all(ngc, ngp2, ntp0)
#ifdef __GPU
          attributes(device) :: ggitp_all
#endif
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
          ggitp(:,:)= matmul(gggmat,geigq(1:ngp1,itq(nqini:nqmax)))

          !$acc data copyin(ggitp, dgeigqk, igcgp2i_) create(zmelp0)
          !$acc kernels
          !$acc loop independent collapse(3)
          do itp = 1, ntp0
            do igp2 = 1, ngp2
              do igc = 1, ngc
                ggitp_all(igc, igp2, itp) = ggitp(igcgp2i_(igc,igp2), itp)
              enddo
            enddo
          enddo
          !$acc end kernels
          ierr = zmm_sb(mm_op_n, mm_op_n, ngc, nt0, ngp2, (1d0,0d0), &
                  &  ggitp_all, ngc, int(ngc*ngp2,8), dgeigqk(1,nmini), ngpmx, 0_8, (0d0,0d0), &
                  &  zmelp0, ngc, int(ngc*nt0,8), ntp0) 
          !$acc kernels
          !$acc loop independent collapse(3)
          do itp = 1, ntp0
            do it = 1, nt0
              do igc = 1, ngc
                zmelp0(igc,it,itp) = phase(igc)*zmelp0(igc,it,itp)
              enddo
            enddo
          enddo
          !$acc end kernels

          allocate (zmeltvv(ngc,nt0,ntp0))
          !$acc data copyin(ppovlinv(1:ngb,1:nbb))
          ierr = zmm(mm_op_n, mm_op_n, ngc, ntp0*nt0, ngc, (1d0,0d0), &
                  &  ppovlinv, ngc, zmelp0, ngc, (0d0,0d0), zmeltvv, ngc) 
          !$acc end data
          !$acc kernels
          zmelt_d(nbloch+1:nbloch+ngc,nctot+1:nctot+nt0,ncc+1:ncc+ntp0) = zmeltvv(1:ngc,1:nt0,1:ntp0)
          !$acc end kernels
          !$acc end data
          deallocate(zmeltvv)

        endblock ZmelIPW
        deallocate(geigq,dgeigqk)
      endif

      allocate(zmel(nbb,ns1:ns2, nqtot))
     !$acc data copyin(ppovlz(1:ngb,1:nbb)) 
     ierr = zmm(mm_op_c, mm_op_n, nbb, nmtot*nqtot, ngb, (1d0,0d0), ppovlz, ngb, zmelt_d, ngb, (0d0,0d0), zmel, ngb) 
     if(zmelconjg) then
       !$acc kernels
       zmel = dconjg(zmel)
       !$acc end kernels
     endif
     !$acc end data
    endblock ZmelBlock
    call cpu_time(t2)
    print '(1x,A,A,2F10.6)','zmel:', acway,(t2-t1)

    call cpu_time(t1)
    TimeConsumingRcxq: block 
      integer :: jpm, it, itp
      integer :: iwmin, iwmax, nprpr, iprpr
      real(kindrcxq), allocatable :: rcxqr(:,:,:), rcxqi(:,:,:)
      ! complex(kindrcxq) :: zwz, zz
      complex(8) :: zwz, zz

      ! print '(A,3I7)', 'iconn:', icounkmink,icounkmaxk, icounkmaxk-icounkmink
      ! print '(A,2I7)', 'npr, nbb:', npr, nbb
      ! print '(A,5I7)', 'ns1, ns2, nmtot, nqtot', ns1, ns2, nmtot, nqtot, nmtot*nqtot
      ! print '(A,3I7)', 'iwmin,iwmax:', iwmin, iwmax, iwmax-iwmin
      nprpr = (npr*(npr+1))/2
      allocate(rcxqr(nprpr,nwhis,npm))
      allocate(rcxqi(nprpr,nwhis,npm))
      iwmin = icouini(icounkmink)
      iwmax = icouini(icounkmaxk)+iwend(icounkmaxk)-iwini(icounkmaxk)
      !$acc data create(rcxqr, rcxqi), copyin(jpmc(icounkmink:icounkmaxk), &
      !$acc&     itc(icounkmink:icounkmaxk), itpc(icounkmink:icounkmaxk), &
      !$acc&     iwini(icounkmink:icounkmaxk), iwend(icounkmink:icounkmaxk), &
      !$acc&     icouini(icounkmink:icounkmaxk), whwc(iwmin:iwmax))

      !$acc kernels
      rcxqr(1:nprpr, 1:nwhis, 1:npm) = 0d0
      !$acc end kernels
      !$acc kernels
      rcxqi(1:nprpr, 1:nwhis, 1:npm) = 0d0
      !$acc end kernels

      !$acc kernels
      !$acc loop independent gang private(it, itp, jpm)
      do icoun = icounkmink, icounkmaxk
        jpm = jpmc(icoun)
        it  = itc (icoun)
        itp = itpc(icoun)

        !$acc loop independent collapse(2) vector private(iw, zwz, zz, iprpr)
        do igb2 = 1, npr
          do igb1 = 1, npr
            if(igb1 > igb2) cycle
            iprpr = ((igb2-1)*igb2)/2 + igb1
            zz = dconjg(zmel(igb1,it,itp))*zmel(igb2,it,itp)
            !$acc loop seq
            do iw=iwini(icoun), iwend(icoun)
              zwz = whwc(iw-iwini(icoun)+icouini(icoun))*zz
              !$acc atomic update
              rcxqr(iprpr, iw, jpm) = rcxqr(iprpr, iw, jpm) + real(zwz,kind=kindrcxq)
              !$acc end atomic
              !$acc atomic update
              rcxqi(iprpr, iw, jpm) = rcxqi(iprpr, iw, jpm) - real(zwz*(0d0,1d0),kind=kindrcxq)
              !$acc end atomic
            enddo
          enddo
        enddo
      enddo
      !$acc end kernels

      !$acc kernels
      !$acc loop independent collapse(4) private(iprpr)
      do jpm = 1, npm
        do iw = 1, nwhis
          do igb2 = 1, npr
            do igb1 = 1, npr
              if(igb1 > igb2) then
                iprpr = ((igb1-1)*igb1)/2 + igb2
                rcxq(igb1, igb2, iw, jpm) = rcxq(igb1, igb2, iw, jpm) &
                                        & + rcxqr(iprpr,iw,jpm) - rcxqi(iprpr,iw,jpm)*(0d0,1d0)
              else
                iprpr = ((igb2-1)*igb2)/2 + igb1
                rcxq(igb1, igb2, iw, jpm) = rcxq(igb1, igb2, iw, jpm) &
                                        & + rcxqr(iprpr,iw,jpm) + rcxqi(iprpr,iw,jpm)*(0d0,1d0)
              endif
            enddo
          enddo
        enddo
      enddo
      !$acc end kernels
      !$acc end data

      deallocate(rcxqr, rcxqi)
    endblock TimeConsumingRcxq
    call cpu_time(t2)
    print '(1x,A,A,2F10.6)','x0:', acway,(t2-t1)

    deallocate(zmel)
  endblock GPUblock
end subroutine x0gpu
