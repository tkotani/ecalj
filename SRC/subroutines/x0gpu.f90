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
  !$ use omp_lib
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

  complex(8):: phasea(natom) 
  integer:: iasx(natom),icsx(natom),imdim(natom),iatomp(natom),nmini,nmmax,nmtot,nqtot,nt0,ntp0
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
  !We have to send all variables SetByCPU. rcxq should be kept in GPU during kloop:do 1500 in x0kf_v4hz
  GPUblock:block
    use m_zmel,only: ppovlz,nbb
    use m_x0kf,only: icounkmink,icounkmaxk, iwini,iwend,itc,itpc,jpmc,icouini,whwc
    use m_genallcf_v3,only: nctot,natom, nlmto, iclass,nlnmv,nlnmc,icore,ncore,ecore,nlnmx
    use m_read_ppovl,only: igggi,igcgp2i,nxi,nxe,nyi,nye,nzi,nze,nvgcgp2,ngcgp,ggg,ppovlinv,nnxi,nnxe, nnyi,nnye,nnzi,nnze,nggg
    use m_read_bzdata,only: qlat
    use m_readVcoud,only: ngc,ngb 
    use m_itq,only: itq,ntq
    use m_stopwatch
#ifdef __GPU
    use openacc
    use cudafor
#endif
    use m_blas, only: m_op_c, m_op_n, m_op_t, zmm, zmm_batch !CPU or GPU versions specifed by macro __GPU
    integer :: ierr
    complex(8), allocatable :: zmel(:,:,:) ! zmel(nbb,nmtot,nqtot), nbb:mixproductbasis, nmtot:middlestate, nq
    complex(8), allocatable :: zmelt_d(:,:,:)
#ifdef __GPU
    attributes(device) :: zmel, zmelt_d
#endif
    integer:: icoun,igb1,igb2,iw,ia
    logical,parameter   :: zmelconjg=.true.

    call stopwatch_start(t_sw_zmel)
    ZmelBlock: block ! ZmelBlock is a part of copy of get_zmel_init in m_zmel.f90
      complex(8):: zmelt(1:nbloch+ngc,nmtot,nqtot)
#ifdef __GPU
      attributes(device) :: zmelt
#endif
      !$acc kernels
      zmelt=0d0
      !$acc end kernels
      ZmelWithinMT: block !- Calculates <psi_q(itp) |psi_qk(it) B_k(rot(r-R))> 
        integer:: i,iap,ias,ib,ic,icp,nc,nc1,nv,ics,itp,iae,ims,ime, ncnv,ncorec,nccc,mdim,it
        complex(8), allocatable :: ppbvphiq_d(:,:,:), cphim_d(:,:), cphiq_d(:,:), ppbc_d(:,:,:), ppbv_d(:,:,:)
#ifdef __GPU
        attributes(device) :: ppbvphiq_d, cphim_d, cphiq_d, ppbc_d, ppbv_d
#endif
        iatomloop: do ia = 1, natom
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

          allocate(ppbvphiq_d(nv,ntp0,mdim))
          !copy from CPU to GPU
          allocate(cphim_d(nv,nt0), source = cphim(ias:iae,nmini:nmmax))
          allocate(cphiq_d(nv,ntp0), source = cphiq(ias:iae,nqini:nqmax))
          allocate(ppbv_d(nv,nv,mdim))
          ppbv_d(1:nv,1:nv,1:mdim) = dcmplx(ppb(nc1:ncnv,nc1:ncnv,1:mdim,icp)) 


          ! valence-valence part
          allocate(zmelt_d(nt0,ntp0,mdim))
          ierr = zmm_batch(ppbv_d, cphiq_d, ppbvphiq_d, nv, ntp0, nv, mdim, opA = m_op_T, sameB = .true.)
          ierr = zmm_batch(cphim_d, ppbvphiq_d, zmelt_d, nt0, ntp0, nv, mdim, alpha = phasea(ia), opA = m_op_C, sameA = .true.)
          !$acc kernels
          do i = 1, mdim
            zmelt(i-1+ims,nctot+1:nctot+nt0,ncc+1:ncc+ntp0) = zmelt_d(1:nt0,1:ntp0,i)
          enddo
          !$acc end kernels
          deallocate(zmelt_d)

          ! core-valence part NOT TESTED
          allocate(zmelt_d(mdim,ntp0,max(1,ncorec)), ppbc_d(mdim,nv,max(1,ncorec)))
          !$acc kernels
          do i = 1, mdim
            ppbc_d(i,1:nv,1:max(1,ncorec)) = dcmplx(ppb(nc1:ncnv,1:(max(1,ncorec)),i,icp))
          enddo
          !$acc end kernels
          ierr = zmm_batch(ppbc_d, cphiq_d, zmelt_d, mdim, ntp0, nv, ncorec, alpha = phasea(ia), sameB = .true.)
          !$acc kernels
          do it = 1, ncorec
            zmelt(ims:ime,ics+it,ncc+1:ncc+ntp0) = zmelt_d(1:mdim,1:ntp0,it)
          enddo
          !$acc end kernels
          deallocate(zmelt_d)
          ! valence-core part NOT TESTED
          allocate(zmelt_d(mdim,nt0,max(nccc,1)))
          ierr =  zmm_batch(ppbc_d, cphim_d, zmelt_d, mdim, nt0, nv, nccc, alpha = phasea(ia), sameB = .true.)
          !$acc kernels
          do itp = 1, nccc
            zmelt(ims:ime,nctot+1:nctot+nt0,ics+itp) = zmelt_d(1:mdim,1:nt0,itp)
          enddo
          !$acc end kernels
          deallocate(zmelt_d, ppbvphiq_d, cphim_d, cphiq_d, ppbv_d, ppbc_d)
        enddo iatomloop
      endblock ZmelWithinMT

      if(ngc/=0)then
        ZmelIPW:block  !> Mattrix elements <Plane psi |psi> from interstitial plane wave.
          integer:: igcgp2,nn(3),iggg,igp1,itp,igc,igp2,igcgp2i_(ngc,ngp2), it
          complex(8):: zmelp0(ngc,nt0,ntp0),phase(ngc) , ggitp(ngcgp,ntp0),gggmat(ngcgp,ngp1)
          complex(8):: ggitp_all(ngc, ngp2, ntp0)

          phase(:)=[(exp( -img*tpi*sum((matmul(symope,kvec)+matmul(qlat,ngveccR(:,igc)))*shtv) ),igc=1,ngc)]  !prepared by CPU
          !$acc data copyin(dgeigqk, geigq, phase, ngvecpB1, ngvecpB2, ngveccR, nadd, ggg(1:nggg), nvgcgp2(1:3,1:ngcgp), &
          !$acc             igggi(nxi:nxe,nyi:nye,nzi:nze), igcgp2i(nnxi:nnxe,nnyi:nnye,nnzi:nnze), ppovlinv(1:ngb,1:nbb)) &
          !$acc      create(zmelp0, nn, gggmat, igcgp2i_, ggitp, ggitp_all)
          !$acc kernels loop independent collapse(2)
          do igcgp2 = 1, ngcgp
            do igp1 =1, ngp1
              nn(1:3) = ngvecpB1(1:3,igp1) - nvgcgp2(1:3,igcgp2) - nadd(1:3)
              if(nn(1)<nxi .OR. nxe<nn(1) .OR. nn(2)<nyi .OR. nye<nn(2) .OR. nn(3)<nzi .OR. nze<nn(3)) cycle
              iggg = igggi(nn(1),nn(2),nn(3))
              if(iggg>=0) gggmat(igcgp2,igp1)=ggg(iggg)
            enddo
          enddo
          !$acc end kernels
          !$acc kernels loop independent collapse(2)
          do igc = 1, ngc
           do igp2 =1, ngp2
              nn(1:3) = ngveccR(1:3,igc) + ngvecpB2(1:3,igp2)
              igcgp2i_(igc,igp2) = igcgp2i(nn(1),nn(2),nn(3))
            enddo
          enddo
          !$acc end kernels
          ierr = zmm(gggmat, geigq(1,itq(nqini)), ggitp, ngcgp, ntp0, ngp1, ldB = ngpmx)

          !$acc kernels loop independent collapse(2)
          do igp2 = 1, ngp2
            do igc = 1, ngc
              ggitp_all(igc, igp2, 1:ntp0) = ggitp(igcgp2i_(igc,igp2), 1:ntp0)
            enddo
          enddo
          !$acc end kernels
          ierr = zmm_batch(ggitp_all, dgeigqk(1,nmini), zmelp0, ngc, nt0, ngp2, ntp0, ldB = ngpmx, sameB = .true.)
          !$acc kernels
          do igc = 1, ngc
            zmelp0(igc,1:nt0,1:ntp0) = phase(igc)*zmelp0(igc,1:nt0,1:ntp0)
          enddo
          !$acc end kernels
          allocate(zmelt_d(ngc,nt0,ntp0))
          ierr = zmm(ppovlinv, zmelp0, zmelt_d, ngc, ntp0*nt0, ngc) 
          !$acc kernels
          zmelt(nbloch+1:nbloch+ngc,nctot+1:nctot+nt0,ncc+1:ncc+ntp0) = zmelt_d(1:ngc,1:nt0,1:ntp0)
          !$acc end kernels
          !$acc end data
          deallocate(zmelt_d)
        endblock ZmelIPW
        deallocate(geigq, dgeigqk)
      endif

      allocate(zmel(nbb,ns1:ns2,nqtot))
      !$acc data copyin(ppovlz(1:ngb,1:nbb))
      ! ierr = zmm(mm_op_c, mm_op_n, nbb, nmtot*nqtot, ngb, (1d0,0d0), ppovlz, ngb, zmelt, ngb, (0d0,0d0), zmel, nbb)
      ierr = zmm(ppovlz, zmelt, zmel, nbb, nmtot*nqtot, ngb, opA = m_op_c)
      if(zmelconjg) then
        !$acc kernels
        zmel = dconjg(zmel)
        !$acc end kernels
      endif
      !$acc end data
    endblock ZmelBlock
    call stopwatch_pause(t_sw_zmel)

    call stopwatch_start(t_sw_x0)
    TimeConsumingRcxq: block 
      integer :: jpm, it, itp, ittp, nttp_max
      integer, allocatable :: nttp(:,:),  itw(:,:,:), itpw(:,:,:)
      complex(8), allocatable :: zw(:,:), wzw(:,:)
      real(8), allocatable :: whw(:,:,:)
#ifdef __GPU
      attributes(device) :: zw, wzw
#endif
      allocate(nttp(nwhis,npm), source = 0)
      do icoun = icounkmink, icounkmaxk
        jpm = jpmc(icoun)
        do iw = iwini(icoun), iwend(icoun)
          nttp(iw,jpm) = nttp(iw,jpm) + 1
        enddo
      enddo

      nttp_max = maxval(nttp(1:nwhis,1:npm))
      allocate (itw(nttp_max,nwhis,npm), source = 0)
      allocate (itpw(nttp_max,nwhis,npm), source = 0)
      allocate (whw(nttp_max,nwhis,npm), source = 0d0)

      nttp(1:nwhis,1:npm) = 0
      do icoun = icounkmink, icounkmaxk
        jpm = jpmc(icoun)
        it  = itc (icoun)
        itp = itpc(icoun)
        do iw = iwini(icoun), iwend(icoun)
          nttp(iw,jpm) = nttp(iw,jpm) + 1
          ittp = nttp(iw,jpm)
          itw(ittp,iw,jpm) = it
          itpw(ittp,iw,jpm) = itp
          whw(ittp,iw,jpm) = whwc(iw-iwini(icoun)+icouini(icoun))
        enddo
      enddo
      ! write(stdo, ftox) 'nttp_max, sum(nttp)', nttp_max, sum(nttp(1:nwhis,1:npm))
      ! do iw = 1, nwhis
      !   if(nttp(iw,1:npm) > 0) write(stdo,ftox) 'iw:', iw, nttp(1:nwhis,1:npm)
      ! enddo
      allocate(wzw(nttp_max,npr), zw(nttp_max,npr))
      !$acc host_data use_device(rcxq)
      !$acc data copyin(whw, itw, itpw)
      do jpm = 1, npm
        !!$omp parallel do private(iw, ittp, igb1, it, itp, zw, wzw) shared(rcxq) schedule(dynamic,10) 
        do iw = 1, nwhis
          if (nttp(iw,jpm) < 1) cycle
          !$acc kernels loop independent collapse(2)
          do ittp = 1, nttp(iw,jpm)
            do igb1 = 1, npr
              it  = itw(ittp,iw,jpm)
              itp = itpw(ittp,iw,jpm)
              zw(ittp,igb1) = zmel(igb1,it,itp)
              wzw(ittp,igb1) = zmel(igb1,it,itp)*whw(ittp,iw,jpm)
            enddo
          enddo
          !$acc end kernels
          ierr = zmm(zw, wzw, rcxq(1,1,iw,jpm), npr, npr, nttp(iw,jpm), &
                  &  opA = m_op_c, beta = (1d0,0d0), ldA = nttp_max, ldB = nttp_max)
        enddo
        !!$omp end parallel do
      enddo
      !$acc end data
      !$acc end host_data
      deallocate(itw, itpw, whw, wzw, zw)

    endblock TimeConsumingRcxq
    call stopwatch_pause(t_sw_x0)
    deallocate(zmel)
  endblock GPUblock
end subroutine x0gpu
