!> get zxq and zxqi for given q
module m_x0kf
  use m_lgunit,only: stdo
  use m_keyvalue,only : Getkeyvalue
  use m_pkm4crpa,only : Readpkm4crpa
  use m_zmel,only: get_zmel_init_gemm, zmel !,get_zmel_init1,get_zmel_init2
  use m_freq,only: npm, nwhis
  use m_genallcf_v3,only:  nsp=>nspin ,ndima,nctot, nband
  use m_read_bzdata,only:  nqbz,ginv,nqibz,  rk=>qbz,wk=>wbz
  use m_rdpp,only: nbloch
  use m_readqg,only: ngpmx,ngcmx
  use m_qbze,only: nqbze
!  use m_readhbe,only: nband
  use m_hamindex,only: ngrp
  use m_tetwt,only:  gettetwt,tetdeallocate, whw,ihw,nhw,jhw,n1b,n2b,nbnb,nbnbx,nhwtot
  use m_ftox
  use m_readVcoud,only:   vcousq,zcousq,ngb,ngc
  use m_kind,only: kp => kindrcxq
  use m_mpi,only: ipr
  implicit none
  public:: x0kf_zxq, deallocatezxq, deallocatezxqi
  complex(kind=kp), public, allocatable:: zxqi(:,:,:)   !Not yet protected because of main_hx0fp0
  complex(kind=kp), public, pointer:: zxq(:,:,:) => null()
  complex(kind=kp), allocatable, target:: rcxq(:,:,:)
  integer,public::npr
  private
  
  integer:: ncount,ncoun
  integer,allocatable:: nkmin(:), nkmax(:),nkqmin(:),nkqmax(:),kc(:)
  integer,allocatable:: icounkmin(:),icounkmax(:)
  
  real(8),public,allocatable:: whwc(:)
  integer,allocatable,public:: iwini(:),iwend(:),itc(:),itpc(:),jpmc(:),icouini(:)
  !
  real(8),public::qrk(3),qq(3)
  integer,public::ns1,ns2,ispm,ispq,nqini,nqmax,icounkmink,icounkmaxk
  logical,external:: cmdopt0 
  logical:: debug = .false.
contains
  function X0kf_v4hz_init(job,q,isp_k,isp_kq, iq, crpa, ikbz_in, fkbz_in) result(ierr) !index accumulation. Initialzation for calling x0kf_v4h
    implicit none
    intent(in)::          job,q,isp_k,isp_kq, iq, crpa
    !! Get ncount index to drive x0kf_v4h. Call job=0 and job=1 successively. 
    !! ncount, ngb and nqibz are keys to estimate computational efforts.
    integer:: irot=1,ierr,isp_k,isp_kq, iq, jpm, ibib, iw,igb2,igb1,it,itp,job,icount, ncc,icoun,k
    real(8):: q(3), imagweight, wpw_k, wpw_kq
    logical:: crpa !,showicount=.false.
    integer, intent(in),optional :: ikbz_in, fkbz_in
    integer :: ikbz, fkbz
    ikbz = 1
    fkbz = nqbz
    if(present(ikbz_in) .and. present(fkbz_in)) then
      ikbz = ikbz_in
      fkbz = fkbz_in
    endif
    if(ipr) write(stdo,'(" x0kf_v4hz_init: job q =",i3,3f8.4)') job,q
    ierr=-1
    ncc=merge(0,nctot,npm==1)
    if(job==0) then
      if(allocated(nkmin)) deallocate(nkmin)
      if(allocated(nkqmin)) deallocate(nkqmin)
      if(allocated(nkmax)) deallocate(nkmax)
      if(allocated(nkqmax)) deallocate(nkqmax)
      allocate( nkmin(ikbz:fkbz),nkqmin(ikbz:fkbz),source= 999999)
      allocate( nkmax(ikbz:fkbz),nkqmax(ikbz:fkbz),source=-999999)
    endif
    if(job==1) then
      if(allocated(whwc)) deallocate(whwc)
      if(allocated(kc)) deallocate(kc)
      if(allocated(iwini)) deallocate(iwini)
      if(allocated(iwend)) deallocate(iwend)
      if(allocated(itc)) deallocate(itc)
      if(allocated(itpc)) deallocate(itpc)
      if(allocated(jpmc)) deallocate(jpmc)
      if(allocated(icouini)) deallocate(icouini)
      allocate( whwc(ncount),kc(ncoun),iwini(ncoun),iwend(ncoun),itc(ncoun), itpc(ncoun), jpmc(ncoun),icouini(ncoun))
      if(allocated(icounkmin)) deallocate(icounkmin)
      if(allocated(icounkmax)) deallocate(icounkmax)
      allocate(icounkmin(ikbz:fkbz),icounkmax(ikbz:fkbz))
    endif
    icount=0
    icoun=0
    AccumulateIndex4icount: do 110 k = 1,nqbz
      if(k < ikbz .OR. k > fkbz) cycle
      if(job==1) icounkmin(k)=icoun+1
      if(job==0) then
        do jpm=1,npm 
           do ibib = 1, nbnb(k,jpm)
            nkmin(k)  = min(n1b(ibib,k,jpm),nkmin(k))
            nkqmin(k) = min(n2b(ibib,k,jpm),nkqmin(k))
            if(n1b(ibib,k,jpm)<=nband) nkmax(k)  = max(n1b(ibib,k,jpm),nkmax(k))
            if(n2b(ibib,k,jpm)<=nband) nkqmax(k) = max(n2b(ibib,k,jpm),nkqmax(k))
          enddo
       enddo
     endif
     flush(stdo)
     !     if(ipr) write(stdo,*)'mm111mmmmmmm22222aaa',nqbz,k,nkqmin(k),nkqmax(k),'job=',job,nbnb(k,1)
     
      if(npm==2.AND.nkqmin(k)/=1)call rx( " When npm==2, nkqmin==1 should be.")
      if (job == 1) then
        ! do jpm = 1, npm
        !   print '(A,2I4,2I5,3I7)', 'nhw(min,max), k,jpm:', k, jpm, minval(nhw(:,k,jpm)), maxval(nhw(:,k,jpm)), &
        !   & sum(nhw(1:nbnbx,k,jpm)), nbnbx, sum(nhw(1:nbnbx,k,jpm))/nbnbx
        ! enddo
      endif
      jpmloop: do 1251 jpm  = 1, npm ! nplusminum=1 usually (or =2)
        ibibloop: do 125 ibib = 1, nbnb(k,jpm) !---  ibib loop, ibib is decomposed to band index pair, it and itp
          ! n1b,n2b --> core after valence.  it,itp --> valence after core 
          if(ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1 >nwhis) call rx( "x0kf_v4hz: iw>nwhis") ! asserting. sanity check
          if(n1b(ibib,k,jpm) > nkmax(k) ) cycle
          if(n2b(ibib,k,jpm) > nkqmax(k)) cycle
          it =merge(nctot+n1b(ibib,k,jpm),             n1b(ibib,k,jpm)-nband,             n1b(ibib,k,jpm)<=nband) !bandindex val or core
          itp=merge(ncc  +n2b(ibib,k,jpm)-nkqmin(k)+1, n2b(ibib,k,jpm)-nkqmin(k)+1-nband, n2b(ibib,k,jpm)<=nband) !bandindex val or core
          if(crpa) then ! constraint RPA mode (juelich verison)
            wpw_k =merge(readpkm4crpa(n1b(ibib,k,jpm),   rk(:,k), isp_k), 0d0, n1b(ibib,k,jpm)<=nband)
            wpw_kq=merge(readpkm4crpa(n2b(ibib,k,jpm), q+rk(:,k), isp_kq),0d0, n2b(ibib,k,jpm)<=nband)
          endif
          icoun=icoun+1
          if(job==1) then
            kc   (icoun) = k
            itc  (icoun)= it
            itpc (icoun)= itp
            jpmc (icoun)= jpm
            iwini(icoun) = ihw(ibib,k,jpm)
            iwend(icoun) = ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1
            icouini(icoun)=icount+1 !icountini is the starting index of iw for give icoun
          endif
          iwloop: do iw = ihw(ibib,k,jpm),ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1 !iw is omega index
            imagweight = whw(jhw(ibib,k,jpm)+iw-ihw(ibib,k,jpm))
            if(crpa) imagweight = imagweight*(1d0-wpw_k*wpw_kq) 
            icount = icount+1 ! icount is including iw count
            if(job==1) whwc (icount)= imagweight
            ! if(ipr) write(stdo,ftox)'x0kf_v4hz_init: icount jpm k it itp iw jpm whw=',icount,jpm, k,it,itp,iw, ftof(whwc(icount))
          enddo iwloop
125     enddo ibibloop
1251  enddo jpmloop
      if(job==1) icounkmax(k)=icoun
110 enddo AccumulateIndex4icount
    ncount = icount
    ncoun  = icoun
    if(job==0.and.ipr) write(stdo,"('x0kf_v4hz_init: job=0 ncount ncoun nqibz=',3i8)") ncount, ncoun,nqibz
    ierr=0
  endfunction x0kf_v4hz_init
  
  subroutine x0kf_zxq(realomega,imagomega, q,iq,npr,schi,crpa,chipm,nolfco,q00,zzr)
    use m_readgwinput,only: ecut, ecuts
    use m_dpsion,only: dpsion5, dpsion_init, dpsion_chiq, dpsion_setup_rcxq
    use m_freq,only: nw_i, nw_w=>nw, niwt=>niw
    use m_readeigen,only:readeval
    use m_freq,only: nw_i,nw,niw 
    use m_zmel,only: Setppovlz,Setppovlz_chipm   ! & NOTE: these data set are stored in this module, and used
    use m_stopwatch
    use m_readVcoud, only: ReleaseZcousq
    use m_mpi,only: comm_k, mpi__rank_k, mpi__size_k, &
                    mpi__ipr_col, mpi__npr_col, mpi__rank_b, mpi__root_k, comm_b
#ifdef __MP
    use m_mpi,only: MPI__reduceSum => MPI__reduceSum_kind4
#else
    use m_mpi,only: MPI__reduceSum
#endif
    use m_gpu, only: use_gpu
    implicit none
    intent(in)::      realomega,imagomega, q,iq,npr,schi,crpa,chipm,nolfco,q00,zzr
    logical:: realomega,imagomega,crpa,chipm,nolfco
    integer:: iq, isp_k,isp_kq,ix0,is,isf,kx,ierr,npr,k, ipr_col, npr_col
    real(8),optional:: q00(3)
    complex(8),optional:: zzr(:,:)
    real(8):: q(3),schi,ekxx1(nband,nqbz),ekxx2(nband,nqbz)
    character(10) :: i2char
    logical :: tetwtk = .false.
    real(8) :: zmel_max_size
    type(stopwatch) :: t_sw_zmel, t_sw_x0, t_sw_dpsion
    qq=q

    ipr_col = mpi__ipr_col(mpi__rank_b) ! start index of column on xq for product basis set
    npr_col = mpi__npr_col(mpi__rank_b) ! number of columns on xq

    if(cmdopt0('--tetwtk'))  tetwtk=.true.
    if(cmdopt0('--emptyrun'))  return
    call getkeyvalue("GWinput","zmel_max_size",zmel_max_size,default=1d0) !in GB
    if(zmel_max_size < 0.001d0) zmel_max_size = 1d0
    if(chipm .AND. nolfco) then; call setppovlz_chipm(zzr,npr)
    else;                        call setppovlz(q,matz=.true.,npr=npr)!2024-5-23 obata. A minor bug to consume memory: Set npr=1 for EPSPP0 mode(no lfc)
    endif
    call ReleaseZcousq() !Release zcousq used in Setppovlz
    if(associated(zxq)) nullify(zxq)
    if(allocated(rcxq)) then
      !$acc exit data delete(rcxq)
      deallocate(rcxq)
    endif
    allocate(rcxq(1:npr,1:npr_col,(1-npm)*nwhis:nwhis)) ! rcxq(:,:,0) is empty until Helbert transformation.
    !$acc enter data create(rcxq)
    if(nw_w > nwhis) call rx('nwhis is smaller than nw_w')
    if(mpi__root_k) then
      if(realomega) then
        zxq(1:,1:,nw_i:) => rcxq(1:npr,1:npr_col,nw_i:nw_w) !nw_i = 0 (npm=1) nw_i = -nw_w (npm=2)
        !$acc enter data create(zxq)
      endif
      if(imagomega) then
        allocate(zxqi(npr,npr_col,niw))
        !$acc enter data create(zxqi)
      endif
    endif
    if(ipr) write(stdo,ftox)' size of rcxq:', npr, npr_col, nwhis*npm+1
    call flush(stdo)
    !$acc kernels
    rcxq(:,:,:) = (0d0,0d0)
    !$acc end kernels
    isloop: do 1103 isp_k = 1,nsp
      GETtetrahedronWeight:block
        isp_kq = merge(3-isp_k,isp_k,chipm) 
        do kx = 1, nqbz
          ekxx1(1:nband,kx) = readeval(  rk(:,kx), isp_k ) ! read eigenvalue
          ekxx2(1:nband,kx) = readeval(q+rk(:,kx), isp_kq) !
        enddo
        if(.not.tetwtk) then
          call gettetwt(q,iq,isp_k,isp_kq,ekxx1,ekxx2,nband=nband) ! tetrahedron weight
          ierr=x0kf_v4hz_init(0,q,isp_k,isp_kq,iq, crpa)
          ierr=x0kf_v4hz_init(1,q,isp_k,isp_kq,iq, crpa) 
          call tetdeallocate()
        endif
      endblock GETtetrahedronWeight
      x0kf_v4hz_block: block !call x0kf_v4hz(q,isp_k,isp_kq,iq, npr,q00,chipm,nolfco,zzr,nmbas)
        integer:: k,jpm, ibib, iw,igb2,igb1,it,itp, nkmax1,nkqmax1, ib1, ib2, ngcx,ix,iy,igb
        integer:: izmel,nmtot,nqtot,iwmax,ifi0,icoucold,icoun, icount, kold
        ! integer:: nwj(nwhis,npm),imb, igc ,neibz,icc,ig,ikp,i,j,itimer
        real(8):: imagweight, wpw_k,wpw_kq,qa,q0a 
        complex(8):: img=(0d0,1d0)
!         zmel0mode: if(cmdopt0('--zmel0')) then ! For epsPP0. Use zmel-zmel0 (for subtracting numerical error) for matrix elements.
!           zmel0block : block
!             real(8)::  q1a,q2a,rfac00
!             complex(kind=kp),allocatable:: zmel0(:,:,:)
!             kold = -999 
!             q1a=sum(q00**2)**.5
!             q2a=sum(q**2)**.5
!             rfac00=q2a/(q2a-q1a)

!             if(tetwtk) then
!               do k = 1, nqbz
!                 if(mod(k-1, mpi__size_k) /= mpi__rank_k)  cycle
!                 call gettetwt(q,iq,isp_k,isp_kq,ekxx1,ekxx2,nband=nband, ikbz_in = k, fkbz_in = k) ! tetrahedron weight
!                 ierr=x0kf_v4hz_init(0,q,isp_k,isp_kq,iq, crpa, ikbz_in = k, fkbz_in = k)
!                 ierr=x0kf_v4hz_init(1,q,isp_k,isp_kq,iq, crpa, ikbz_in = k, fkbz_in = k) 
!                 call tetdeallocate()

!                 call x0kf_zmel(q00,k, isp_k,isp_kq)
!                 if(allocated(zmel0)) deallocate(zmel0)
!                 allocate(zmel0,source=zmel)
!                 call x0kf_zmel(q, k, isp_k,isp_kq)
!                 do icoun = icounkmin(k), icounkmax(k)
!                   jpm = jpmc(icoun)
!                   it  = itc (icoun)
!                   itp = itpc(icoun)
!                   do iw=iwini(icoun),iwend(icoun)
!                     icount= icouini(icoun)+iw-iwini(icoun)
!                     if(abs(zmel0(1,it,itp))>1d10) cycle
!                     rcxq(1,1,iw*(3-2*jpm))=rcxq(1,1,iw*(3-2*jpm)) +rfac00**2*(abs(zmel(1,it,itp))-abs(zmel0(1,it,itp)))**2 *whwc(icount)
!                   enddo
!                 enddo
!               enddo
!             else
!             zmel0modeicount: do icoun = 1,ncoun 
!               k   = kc(icoun)
!               it  = itc (icoun) !occ      k
!               itp = itpc(icoun) !unocc    q+k
!               jpm = jpmc(icoun) ! \pm omega. Usual mode is only for jpm=1
!               if(mod(k-1, mpi__size_k) /= mpi__rank_k)  cycle
!               if(kold/=k) then
!                 call x0kf_zmel(q00,k, isp_k,isp_kq)!, GPUTEST=GPUTEST)
!                 if(allocated(zmel0)) deallocate(zmel0)
!                 allocate(zmel0,source=zmel)
!                 call x0kf_zmel(q, k, isp_k,isp_kq)!, GPUTEST=GPUTEST)
!                 kold=k
!                 if(ipr) write(stdo,*) 'k, mpi__rank_k', k, mpi__rank_k
!                 call flush(6)
!               endif
! !              write(*,*)'zzzzzzzzzzzzzmel',shape(zmel)
!               do iw=iwini(icoun),iwend(icoun) !iw  = iwc(icount)  !omega-bin
!                 icount= icouini(icoun)+iw-iwini(icoun)
!                 if(abs(zmel0(1,it,itp))>1d10) cycle !We assume rcxq(1) in this mode
!                 rcxq(1,1,iw*(3-2*jpm))=rcxq(1,1,iw*(3-2*jpm)) +rfac00**2*(abs(zmel(1,it,itp))-abs(zmel0(1,it,itp)))**2 *whwc(icount)
!               enddo
!             enddo zmel0modeicount
!             endif
!             !$acc update device (rcxq)
!           endblock zmel0block
!           goto 2000 
!         endif zmel0mode
        if(cmdopt0('--emptyrun')) goto 1590
        call cputid (0)
!        if(GPUTEST) then
          ! rcxq(ibg1,igb2,iw) = \sum_ibib wwk(iw,ibib)* <M_ibg1(q) psi_it(k)| psi_itp(q+k)> < psi_itp | psi_it M_ibg2 > at q
          call stopwatch_init(t_sw_zmel, 'zmel_'//'gemm') !merge('gpu','ori',mask = GPUTEST))
          call stopwatch_init(t_sw_x0, 'x0_'//'gemm') !merge('gpu','ori',mask = GPUTEST))
          kloop:do 1500 k=1,nqbz !zmel = < M(igb q) phi( rk it occ)|  phi(q+rk itp unocc)>
            if(mod(k-1, mpi__size_k) /= mpi__rank_k)  cycle
            if(tetwtk) then
              call gettetwt(q,iq,isp_k,isp_kq,ekxx1,ekxx2,nband=nband, ikbz_in = k, fkbz_in = k) ! tetrahedron weight
              ierr=x0kf_v4hz_init(0,q,isp_k,isp_kq,iq, crpa, ikbz_in = k, fkbz_in = k)
              ierr=x0kf_v4hz_init(1,q,isp_k,isp_kq,iq, crpa, ikbz_in = k, fkbz_in = k) 
              call tetdeallocate()
            endif
            ! qq   = q;              qrk  = q+rk(:,k)
            ! ispm = isp_k;          ispq = isp_kq
            ! ns1  = nkmin(k)+nctot; ns2  = nkmax(k)+nctot
            ! nqini= nkqmin(k);      nqmax= nkqmax(k)
            icounkmink= icounkmin(k); icounkmaxk= icounkmax(k)
            debug=cmdopt0('--debugzmel')
            if(debug.and.ipr) write(stdo,ftox) 'ggggggggg goto get_zmel_init_gemm',k, nkmin(k),nkmax(k),nctot
            NMBATCH: BLOCK 
              integer :: nsize, nns, ibatch, nbatch, ns12
              integer, allocatable :: ns1lists(:), ns2lists(:)
              nsize = (nkqmax(k)-nkqmin(k))*npr   !
              nns = (nkmax(k) - nkmin(k) + 1) !number of middle states
              nbatch = ceiling(dble(nns)*nsize*16/1000**3/zmel_max_size)
              allocate(ns1lists(nbatch), ns2lists(nbatch))
              ns1 = nkmin(k) + nctot 
              do ibatch = 1, nbatch
                ns12 = (nns + ibatch - 1)/nbatch
                ns1lists(ibatch) = ns1
                ns2lists(ibatch) = ns1 + ns12 - 1
                ns1 = ns2lists(ibatch) + 1
              enddo
              do ibatch = 1, nbatch
                ns1 = ns1lists(ibatch)
                ns2 = ns2lists(ibatch)
                ns12 = ns2 - ns1 + 1
                if(ns12 == 0) cycle
                call stopwatch_start(t_sw_zmel)
                if(ipr) write(stdo,ftox) 'zmel_batch:', ibatch, ns1, ns2, nbatch
                if(use_gpu) then
                  !Currently, mpi version of get_zmel_init_gpu which is available by adding comm argument for MPI communicator,
                  !but, MPI communication is significant bottle-neck in the case where GPUs are used. Therefore, it is only used in without GPU case.
                  call get_zmel_init_gemm(q=q+rk(:,k), kvec=q, irot=1, rkvec=q, ns1=ns1,ns2=ns2, ispm=isp_k, &
                       nqini=nkqmin(k),nqmax=nkqmax(k), ispq=isp_kq,nctot=nctot, ncc=merge(0,nctot,npm==1),iprx=.false., &
                       zmelconjg=.true.)
                else
                  call get_zmel_init_gemm(q=q+rk(:,k), kvec=q, irot=1, rkvec=q, ns1=ns1,ns2=ns2, ispm=isp_k, &
                       nqini=nkqmin(k),nqmax=nkqmax(k), ispq=isp_kq,nctot=nctot, ncc=merge(0,nctot,npm==1),iprx=.false., &
                       zmelconjg=.true., comm = comm_b)
                endif
                call stopwatch_pause(t_sw_zmel)
                call stopwatch_start(t_sw_x0)
                call x0gemm(rcxq, npr, ipr_col, npr_col, nwhis, npm, ns1, ns2)
                call stopwatch_pause(t_sw_x0)
              enddo
              deallocate(ns1lists, ns2lists)
            END BLOCK NMBATCH
            if(ipr) write(stdo,ftox) 'end of k:', k ,' of:',nqbz, 'zmel:', ftof(stopwatch_lap_time(t_sw_zmel),4), '(sec)', &
                                                        ' x0:', ftof(stopwatch_lap_time(t_sw_x0),4), '(sec)'
            call flush(6)
1500      enddo kloop
!         else ! NOTE: kloop10:do 1510 is equivalent to do 1500. 2024-3-25
!           kloop10:do 1510 k=1,nqbz !zmel = < M(igb q) phi( rk it occ)|  phi(q+rk itp unocc)>, where it=nm1:nm2
!             if(mod(k-1, mpi__size_k) /= mpi__rank_k)  cycle
!             call stopwatch_start(t_sw_zmel)
!             if(cmdopt0('--emptyrun')) cycle
!             call get_zmel_init(q=q+rk(:,k), kvec=q, irot=1, rkvec=q, ns1=nkmin(k)+nctot,ns2=nkmax(k)+nctot, ispm=isp_k, &
!                  nqini=nkqmin(k),nqmax=nkqmax(k), ispq=isp_kq,nctot=nctot, ncc=merge(0,nctot,npm==1),iprx=.false.,zmelconjg=.true.)
!             call stopwatch_pause(t_sw_zmel)
!             call stopwatch_start(t_sw_x0)
!             icounloop: do 1000 icoun=icounkmin(k),icounkmax(k)
!               TimeConsumingRcxq: block 
!                 complex(8):: zmelzmel(npr,npr_col)
!                 if(debug) if(ipr) write(stdo,ftox)'icoun: iq k jpm it itp n(iw)=',icoun,iq,k,jpm,it,itp,iwend(icoun)-iwini(icoun)+1
!                 associate( &
!                      jpm => jpmc(icoun),&  !\pm omega 
!                      it  => itc (icoun),&  !occ      at k
!                      itp => itpc(icoun))   !unocc    at q+k
!                   do concurrent(igb1=1:npr,igb2=1:npr_col) 
!                     zmelzmel(igb1,igb2)= dconjg(zmel(igb1,it,itp))*zmel(igb2+ipr_col-1,it,itp) 
!                   enddo
!                   forall(iw=iwini(icoun):iwend(icoun))& !rcxq is hermitian, thus, we can reduce computational time half.
!                        rcxq(:,:,iw,jpm)=rcxq(:,:,iw,jpm)+ whwc(iw-iwini(icoun)+icouini(icoun))* zmelzmel(:,:) !may use zaxpy and symmetrize afterwards
!                   !forall(iw=iwini(icoun):iwend(icoun)) nwj(iw,jpm)=nwj(iw,jpm)+iwend(icoun)-iwini(icoun)+1 !counter check
!                 endassociate
!               endblock TimeConsumingRcxq
! 1000        enddo icounloop
!             call stopwatch_pause(t_sw_x0)
!             if(ipr) write(stdo,ftox) 'end of k:', k ,' of:',nqbz, 'zmel:', ftof(stopwatch_lap_time(t_sw_zmel),4), '(sec)', &
!                                                         ' x0:', ftof(stopwatch_lap_time(t_sw_x0),4), '(sec)'
!             call flush(6)
! 1510      enddo kloop10
!        endif
        call stopwatch_show(t_sw_zmel)
        call stopwatch_show(t_sw_x0)
        call cputid (0)
1590    continue
2000   continue
       if(debug.and.ipr) write(stdo,ftox)"--- x0kf_v4hz: end: sumcheck abs(rcxq)=",sum(abs(rcxq(:,:,:)))
      endblock x0kf_v4hz_block
      deallocate(whwc, kc, iwini,iwend, itc,itpc, jpmc,icouini, nkmin,nkmax,nkqmin,nkqmax,icounkmin,icounkmax)
      HilbertTransformation:if(isp_k==nsp .OR. chipm) then
        !Get real part. When chipm=T, do dpsion5 for every isp_k; When =F, do dpsion5 after rxcq accumulated for spins
        mpi_k_accumulate: block
          !To reduce memory allocation size in MPI__reduceSUM, mpi_reduce operation has been split for each iw and jpm
          integer :: jpm, iw
          if(mpi__size_k > 1) then
            !$acc update host(rcxq)
            do jpm=1, npm
              do iw=1, nwhis
                call MPI__reduceSum(0, rcxq(1,1,iw*(3-2*jpm)), npr*npr_col, communicator = comm_k) 
              enddo
            enddo
            !$acc update device(rcxq)
          endif
        end block mpi_k_accumulate
        if(mpi__root_k) then
          call stopwatch_init(t_sw_dpsion, 'dpsion') !merge('gpu','ori',mask = GPUTEST))
          call stopwatch_start(t_sw_dpsion)
          call dpsion_init(realomega, imagomega, chipm)
          !$acc host_data use_device(rcxq, zxqi)
          call dpsion_chiq(realomega, imagomega, chipm, rcxq, zxqi, npr, npr_col, schi, isp_k, ecut)
          !$acc end host_data
          call stopwatch_pause(t_sw_dpsion)
          call stopwatch_show(t_sw_dpsion)
        endif
        !set zero (if isp_k == 1) or chi+- in rcxq (if isp_k == 2)
        if(chipm) call dpsion_setup_rcxq(rcxq, npr, npr_col, isp_k)
      endif HilbertTransformation 
1103 enddo isloop
  end subroutine x0kf_zxq
  subroutine deallocatezxq()
    !$acc exit data delete(zxq)
    nullify(zxq)
  end subroutine deallocatezxq
  subroutine deallocatezxqi()
    !$acc exit data delete(zxqi)
    deallocate(zxqi)
  end subroutine deallocatezxqi
  subroutine x0kf_zmel( q,k, isp_k,isp_kq)!, GPUTEST) ! Return zmel= <phi phi |M_I> in m_zmel
    use m_mpi, only: comm_b
    intent(in)   ::     q,k, isp_k,isp_kq   
    integer::              k,isp_k,isp_kq 
    real(8)::           q(3)
    debug=cmdopt0('--debugzmel')
    if(debug.and.ipr) write(stdo,ftox) 'ggggggggg goto get_zmel_init_gemm',k, nkmin(k),nkmax(k),nctot
    call get_zmel_init_gemm(q=q+rk(:,k), kvec=q, irot=1, rkvec=q, ns1=nkmin(k)+nctot,ns2=nkmax(k)+nctot, ispm=isp_k, &
         nqini=nkqmin(k),nqmax=nkqmax(k), ispq=isp_kq,nctot=nctot, ncc=merge(0,nctot,npm==1),iprx=.false., zmelconjg=.true.)
       !$acc update host(zmel)
  end subroutine x0kf_zmel

end module m_x0kf

!! === calculate chi0, or chi0_pm === 
!!
!! ppovl= <I|J> = O , V_IJ=<I|v|J>
!! (V_IJ - vcoud_mu O_IJ) Zcousq(J, mu)=0, where Z is normalized with O_IJ.
!! <I|v|J>= \sum_mu ppovl*zcousq(:,mu) v^mu (Zcousq^*(:,mu) ppovl)
!!
!! zmelt contains O^-1=<I|J>^-1 factor. Thus zmelt(phi phi J)= <phi |phi I> O^-1_IJ
!! ppovlz(I, mu) = \sum_J O_IJ Zcousq(J, mu)
!!
!!  rcxq (npr,npr,nwhis,npm): for given q,
!!       rcxq(I,J,iw,ipm) = Im (chi0(omega))= \sum_k <I_q psi_k|psi_(q+k)> <psi_(q+k)|psi_k> \delta(\omega- (e_i-ej))
!!        When npm=2 we calculate negative energy part. (time-reversal asymmetry)
!!
! note: zmel: matrix element <phi phi |M>
!! ndima   = total number of atomic basis functions within MT
!! nqbz    = number of k-points in the 1st BZ

! z1p = <M_ibg1 psi_it | psi_itp> < psi_itp | psi_it M_ibg2 >
!  zxq(iw,ibg1,igb2) = \sum_ibib imgw(iw,ibib)* z1p(ibib, igb1,igb2) !ibib means band pair (occ,unocc)
!
!  zzmel(1:nbloch, ib_k,ib_kq)
!      ib_k =[1:nctot]              core
!      ib_k =[nctot+nkmin:nctot+nkmax]  valence
!      ib_kq =[1:ncc]             core
!      ib_kq =[ncc+nkqmin:ncc+nkqmax]  valence range [nkqmin,nkqmax]
!   If jpm=1, ncc=0.  !   If jpm=2, ncc=ncore. nkqmin(k)=1 should be.
! NOTE:
!  q+rk n2b vec_kq  vec_kq_g geig_kq cphi_kq  ngp_kq ngvecp_kq  isp_kq
!    rk n1b vec_k   vec_k_g  geig_k  cphi_k   ngp_k  ngvecp_k   isp_k
!! -------------------------------------------------------------------------------
! note: for usual correlation mode, I think nctot=0
!!--- For dielectric funciton, we use irot=1 kvec=rkvec=q
!            < MPB      middle   |   end >
!!              q      rkvec     | q + rkvec
!                      nkmin:nt0 | nkqmin:ntp0
!                         occ    | unocc
!                      (nkmin=1)
!                      (cphi_k  | cphi_kq !in x0kf)
!!     rkvec= rk(:,k)   ! <phi(q+rk,nqmax)|phi(rk,nctot+nmmax)  MPB(q,ngb )>
!!     qbz_kr= rk(:,k)  !
!!     qibz_k= rk(:,k)  ! k
!! Get_zmelt in m_zmel gives the matrix element zmel,  ZO^-1 <MPB psi|psi> , where ZO is ppovlz 
!! zmel(ngb, nctot+nt0,  ncc+ntp0) in m_zmel
!            nkmin:nt0, nkqmin:ntp0, where nt0=nkmax-nkmin+1  , ntp0=nkqmax-nkqmin+1
