!> get zxq and zxqi for given q
module m_x0kf_ahc
  use m_lgunit,only: stdo
  use m_keyvalue,only : Getkeyvalue
  use m_pkm4crpa,only : Readpkm4crpa
  use m_zmel,only: get_zmel_init_gemm, zmel !,get_zmel_init1,get_zmel_init2
  use m_freq,only: npm, nwhis
  use m_genallcf_v3,only:  nsp=>nspin,nspx,nspc,ndima,nctot, nband,plat,alat
  use m_read_bzdata,only:  nqbz,ginv,nqibz,qbz,  rk=>qbz,wk=>wbz,wik=>wibz
  use m_rdpp,only: nbloch
  use m_readqg,only: ngpmx,ngcmx
  use m_qbze,only: nqbze
!  use m_readhbe,only: nband
  use m_hamindex,only: ngrp,symops
  use m_tetwt,only:  gettetwt,tetdeallocate, whw,ihw,nhw,jhw,n1b,n2b,nbnb,nbnbx,nhwtot
  use m_ftox
  use m_readVcoud,only:   vcousq,zcousq,ngb,ngc
  use m_kind,only:kindrcxq
  use m_setqibz_lmfham,only: set_qibz,irotg
  implicit none
  public:: x0kf_ahc, deallocatezxq, deallocatezxqi
  complex(8),public,allocatable:: zxq(:,:,:), zxqi(:,:,:)   !Not yet protected because of main_hx0fp0
  complex(kindrcxq),allocatable:: rcxq(:,:,:,:)
  real(8),allocatable::omegaahc_tet(:,:,:), omegaahc_sp(:,:,:), omegasum_tet(:,:,:), omegasum_sp(:,:,:)
!  complex(4),allocatable:: rcxq4(:,:,:,:)
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
  logical:: debug
contains
  function X0kf_v4hz_init(job,q,isp_k,isp_kq, iq, crpa) result(ierr) !index accumulation. Initialzation for calling x0kf_v4h
    implicit none
    intent(in)::          job,q,isp_k,isp_kq, iq, crpa
    !! Get ncount index to drive x0kf_v4h. Call job=0 and job=1 successively. 
    !! ncount, ngb and nqibz are keys to estimate computational efforts.
    integer:: irot=1,ierr,isp_k,isp_kq, iq, jpm, ibib, iw,igb2,igb1,it,itp,job,icount, ncc,icoun,k
    real(8):: q(3), imagweight, wpw_k, wpw_kq
    logical:: crpa !,showicount=.false.
    write(stdo,'(" x0kf_v4hz_init: job q =",i3,3f8.4)') job,q
    ierr=-1
    ncc=merge(0,nctot,npm==1)
    if(job==0) allocate( nkmin(nqbz),nkqmin(nqbz),source= 999999)
    if(job==0) allocate( nkmax(nqbz),nkqmax(nqbz),source=-999999)
    if(job==1) allocate( whwc(ncount),kc(ncoun),iwini(ncoun),iwend(ncoun),itc(ncoun), itpc(ncoun), jpmc(ncoun),icouini(ncoun))
    if(job==1) allocate(icounkmin(nqbz),icounkmax(nqbz))
    icount=0
    icoun=0
    AccumulateIndex4icount: do 110 k = 1,nqbz
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
      if(npm==2.AND.nkqmin(k)/=1)call rx( " When npm==2, nkqmin==1 should be.")
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
            ! write(stdo,ftox)'x0kf_v4hz_init: icount jpm k it itp iw jpm whw=',icount,jpm, k,it,itp,iw, ftof(whwc(icount))
          enddo iwloop
125     enddo ibibloop
1251  enddo jpmloop
      if(job==1) icounkmax(k)=icoun
110 enddo AccumulateIndex4icount
    ncount = icount
    ncoun  = icoun
    if(job==0) write(stdo,"('x0kf_v4hz_init: job=0 ncount ncoun nqibz=',3i8)") ncount, ncoun,nqibz
    ierr=0
  endfunction x0kf_v4hz_init
  
  subroutine x0kf_ahc(realomega,imagomega, q,iq,npr,schi,crpa,chipm,nolfco,q00,zzr)
    use m_readgwinput,only: ecut,ecuts
    use m_dpsion,only: dpsion5
    use m_freq,only: nw_i, nw_w=>nw, niwt=>niw
    use m_readeigen,only:readeval
    use m_freq,only: nw_i,nw,niw 
    use m_zmel,only: Setppovlz,Setppovlz_chipm   ! & NOTE: these data set are stored in this module, and used
    use m_stopwatch
    use m_mpi,only: comm_k, mpi__rank_k, mpi__size_k, MPI__reduceSum, &
         mpi__ipr_col, mpi__npr_col, mpi__rank_b, mpi__root_k, comm_b,&
         MPI__AllreduceSumReal, MPI__AllreduceSumRealSca
    use m_gpu, only: use_gpu
!    use m_data_gpu, only: SetDataGPU_inkx, ExitDataGPU_inkx
    use m_procar, only: read_sdenmat, sdendwgtall
    implicit none
    intent(in)::      realomega,imagomega, q,iq,npr,schi,crpa,chipm,nolfco,q00,zzr
    logical:: realomega,imagomega,crpa,chipm,nolfco
    integer:: iq, isp_k,isp_kq,ix0,is,isf,kx,ierr,npr,k, ipr_col, npr_col,ifx0
    real(8),optional:: q00(3)
    complex(8),optional:: zzr(:,:)
    real(8):: q(3),schi,ekxx1(nband,nqbz),ekxx2(nband,nqbz),frr,wgt,vol,tripl
    character(10) :: i2char
    type(stopwatch) :: t_sw_zmel, t_sw_x0
    logical:: cmdopt0
    real(8),allocatable::ku(:,:),kbu(:,:,:)
    real(8),allocatable::wbb(:),bb(:,:)
    integer::iko_ixs(2),iko_fxs(2),nbb,ifahc,nband_k,nband_kq,ibb1,ibb2,ix1,ix2
    integer,allocatable::ikbidx(:,:)
    complex(8),allocatable::zmelahc_k(:,:,:),zmelahc_kq(:,:,:),uumat_k(:,:,:,:,:),zzahc(:,:,:,:) !,uumat_kq(:,:,:,:,:)
    complex(8),parameter::img=(0d0,1d0)
    vol = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
    !    npr=nprin
    qq=q
!    GPUTEST = .true. !cmdopt0('--gpu')

    ipr_col = mpi__ipr_col(mpi__rank_b) ! start index of column on xq for product basis set
    npr_col = mpi__npr_col(mpi__rank_b) ! number of columns on xq

    if(mpi__root_k) then
      if(realomega) allocate(zxq(npr,npr_col,nw_i:nw),source=(0d0,0d0))
      if(imagomega) allocate(zxqi(npr,npr_col,niw),source=(0d0,0d0))
    endif
    if(cmdopt0('--emptyrun'))  return
    if(chipm .AND. nolfco) then; call setppovlz_chipm(zzr,npr)
    else;                        call setppovlz(q,matz=.true.,npr=npr)!2024-5-23 obata. A minor bug to consume memory: Set npr=1 for EPSPP0 mode(no lfc)
    endif
    if(cmdopt0('--qibzonly')) call set_qibz(plat,qbz,nqbz,symops,ngrp)
    
    ReadinBBVEC: if(cmdopt0('--ahc')) then
       RBBVEC:block
         integer:: nqbz2,nspin2,ifbb,iqr,ib,ibr,iq
         open(newunit=ifbb,file='BBVEC')
         read(ifbb,*)
         read(ifbb,*) nbb, nqbz2
         if(nqbz /= nqbz2) call rx('x0kf_v4h: nqbz /= nqbz2')
         allocate(bb(3,nbb),wbb(nbb),ikbidx(nbb,nqbz))
         do ib = 1,nbb
            read(ifbb,*) bb(1,ib),bb(2,ib),bb(3,ib),wbb(ib)
            write(stdo,ftox)'bbbbb bb wbb=',ftof(bb(:,ib)),wbb(ib)
         enddo
!         allocate(ku(3,nqbz),kbu(3,nbb,nqbz))
         do iq = 1,nqbz
            read(ifbb,*) !iqr, ku(:,iq)
            do ib = 1,nbb
               read(ifbb,*) iqr,ibr,ikbidx(ib,iq)!,kbu(:,ib,iq)
            enddo
         enddo
         read(ifbb,*)
         read(ifbb,*) nspin2
         if(nspx /= nspin2) call rx('x0kv_4h: nsp /= nspin2!')
         do is = 1, nspx
            read(ifbb,*) iko_ixs(is), iko_fxs(is)
         enddo
         close(ifbb)
       endblock RBBVEC
    endif ReadinBBVEC
    
!    call SetDataGPU_inkx(set_ppovlz_in_gpu = .false.)
    allocate(uumat_k(nband,nband,nbb,nqbz,nspx),source=dcmplx(0d0,0d0))!,uumat_kq(nband_kq,nband_kq,nbb,nqbz))
    call readuun(ikbidx,iko_ixs,iko_fxs,nband,nqbz,nbb,nspx, uumat_k)  ! read <u|u> matrix for isp_k
    if(debug)write(stdo,ftox)'uuuuusumabs=',nband,nqbz,nbb,nspx,sum(abs(uumat_k))
    ! if (cmdopt0('--mloahc')) then
    !    ReadCMLO: block
    !      integer::iqibz,iq,ificmlo,nMLO
    !      real(8)::qp(***)
    !      do iqibz=1,nqibz
    !         open(newunit=ificmlo,file='Cmlo'//trim(xt(iqibz))//trim(xt(isp_k)),form='unformatted')
    !         read(ificmlo) iq,qp,nMLO
    !         read(ificmlo) cmlo(1:nMLO)
    !         close(ificmlo)
    !      enddo
    !    endblock ReadCMLO
    !    GetMLOUU: block
    !      integer::i
    !      !!! <u_kn MLO|u_k+bm MLO> = (CMLO +)^-1 <u_kn PMT|u_k+bm PMT> (CMLO)^-1
    !    endblock GetMLOUU
    ! endif
    
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     FecheckSO1nooffSO2: block
!       integer:: im,in,ix
!       ibb1=1
!       ibb2=2
!       k=2
!       if(nspc==2) then
!          im=14
!          in=16
!          isp_k=1
!       else   
!          im=5
!          in=7
!          isp_k=2
!       endif
!       write(stdo,ftox)'sssmmm mat', &!ftof(ekxx0(im)),ftof(ekxx0(in)),'it itp=',im,in,&
!            'abs=',ftod(abs(uumat_k(im,in,ibb1,k,isp_k))),ftod(abs(uumat_k(im,in,ibb2,k,isp_k)))
!       ix=0
!       do in=1,nband
!          if(abs(uumat_k(im,in,ibb2,k,isp_k))>1d-8) then
!             ix=ix+1 
!             write(stdo,ftox)'       ',im,in,' ',ix,ftod(abs(uumat_k(im,in,ibb2,k,isp_k)))
!          endif
!       enddo
!     endblock FecheckSO1nooffSO2
! !!!!!!!!!!!!!!!!!!!!!!!!
        
    nband_k=nband
    isloop: do 1103 isp_k = 1,nspx
       isqloop: do 1104 isp_kq = 1, nspx
          if(isp_k/=isp_kq) cycle
          ReadinUU: if (cmdopt0('--ahc')) then
             nband_k = nband !iko_fxs(isp_k)
!             nband_kq=nband !iko_fxs(isp_kq)
             ! WriteUU: if (cmdopt0('--UUMAT')) then
             !    UUMATblock: block
             !      integer :: iq, ib, ifuu, ibnd1, ibnd2
             !      character(2) :: char
             !      character(1) :: charisp
             !      do iq = 1, nqbz
             !         write(char,'(i2.2)') iq
             !         write(charisp,'(i1.1)') isp_k
             !         open(newunit=ifuu, file='UUMAT-k'//trim(char)//'.isp'//trim(charisp))
             !         write(ifuu,*) "### kp(1:3)  bb(1:3)  i  j  Re(<u|u>)  Im(<u|u>)  |<u|u>| "
             !         do ib = 1, nbb
             !            do ibnd1 = 1, nband_k
             !               do ibnd2 = 1, nband_k
             !                  if (abs(uumat_k(ibnd1,ibnd2,ib,iq,isp_k))>1d10) then
             !                     write(ifuu,'(3f10.5,3x,3f10.5,3x,i3,x,i3,x,3f16.9)') rk(:,iq), bb(:,ib), ibnd1, ibnd2, &
             !                          0d0, 0d0, 0d0
             !                  else
             !                     write(ifuu,'(3f10.5,3x,3f10.5,3x,i3,x,i3,x,3f16.9)') rk(:,iq), bb(:,ib), ibnd1, ibnd2, &
             !                          dreal(uumat_k(ibnd1,ibnd2,ib,iq,isp_k)), &
             !                          dimag(uumat_k(ibnd1,ibnd2,ib,iq,isp_k)), abs(uumat_k(ibnd1,ibnd2,ib,iq,isp_k))
             !                  endif
             !               enddo
             !            enddo
             !         enddo
             !      enddo
             !      close(ifuu)
             !    endblock UUMATBLOCK
             ! endif WriteUU
          endif ReadinUU
          GETtetrahedronWeight:block
            ! isp_kq = merge(3-isp_k,isp_k,chipm.and.nspx==2)
            write(*,*) isp_k, isp_kq
            do kx = 1, nqbz
               ekxx1(1:nband,kx) = readeval(  rk(:,kx), isp_k ) ! read eigenvalue
               ekxx2(1:nband,kx) = readeval(q+rk(:,kx), isp_kq) !
            enddo
            call gettetwt(q,iq,isp_k,isp_kq,ekxx1,ekxx2,nband=nband) ! tetrahedron weight
            ierr=x0kf_v4hz_init(0,q,isp_k,isp_kq,iq, crpa)
            ierr=x0kf_v4hz_init(1,q,isp_k,isp_kq,iq, crpa) 
            call tetdeallocate()
          endblock GETtetrahedronWeight
          x0kf_v4hz_block: block !call x0kf_v4hz(q,isp_k,isp_kq,iq, npr,q00,chipm,nolfco,zzr,nmbas)
            use m_freq,only:frhis
            integer:: k,jpm, ibib, iw,igb2,igb1,it,itp, nkmax1,nkqmax1, ib1, ib2, ngcx,ix,iy,igb
            integer:: izmel,nmtot,nqtot,ierr,iwmax,ifi0,icoucold,icoun
            integer:: nwj(nwhis,npm),imb, igc ,neibz,icc,ig,ikp,i,j,itimer,icount, kold
            integer:: im,in,ib,ic,iq
            real(8):: imagweight, wpw_k,wpw_kq,qa,q0a,efshift=0d0
            complex(8):: img=(0d0,1d0)
            logical :: cmdopt0,cmdopt2,GPUTEST
            character*4:: charnum4
            character*7:: charnum7
            character(8):: charext
            character(256)::strn
            if(.not.allocated(rcxq)) then
               allocate(rcxq(npr,npr_col,nwhis,npm))
               rcxq=0d0
               !           if(GPUTEST) then
               write(stdo,ftox)' size of rcxq:', npr, npr_col, nwhis, npm
               !$acc enter data create(rcxq) 
               !$acc kernels
               rcxq(1:npr,1:npr_col,1:nwhis,1:npm) = (0d0,0d0)
               !$acc end kernels
               !           endif
            endif
            if (.not. allocated(omegaahc_tet)) then
               allocate(omegaahc_tet(3,3,nqbz),source=0d0)
               allocate(omegaahc_sp(3,3,nqbz),source=0d0)
               allocate(omegasum_tet(3,3,nqbz),source=0d0)
               allocate(omegasum_sp(3,3,nqbz),source=0d0)
            endif
            if (cmdopt2('-EfermiShifteV=',strn)) then
               read(strn,*) efshift
            endif
            ! SDENMAT: block
            !   use m_procar,only: m_sden_add
            !   integer::iq
            !   do iq=1,nqbz
            !      call m_sden_add(iq,isp_k,ef0*,evl*,q,nband,evec*,ndima) ! ***
            !   enddo
            ! endblock SDENMAT
            zmel0mode: if(cmdopt0('--zmel0')) then ! For epsPP0. Use zmel-zmel0 (for subtracting numerical error) for matrix elements.
               zmel0block : block
                 real(8)::  q1a,q2a,rfac00
                 complex(8),allocatable:: zmel0(:,:,:)
                 kold = -999 
                 q1a=sum(q00**2)**.5
                 q2a=sum(q**2)**.5
                 rfac00=q2a/(q2a-q1a)
                 zmel0modeicount: do icoun = 1,ncoun 
                    k   = kc(icoun)
                    it  = itc (icoun) !occ      k
                    itp = itpc(icoun) !unocc    q+k
                    jpm = jpmc(icoun) ! \pm omega. Usual mode is only for jpm=1
                    if(mod(k-1, mpi__size_k) /= mpi__rank_k)  cycle
                    if(kold/=k) then
                       call x0kf_zmel(q00,k, isp_k,isp_kq)!, GPUTEST=GPUTEST)
                       if(allocated(zmel0)) deallocate(zmel0)
                       allocate(zmel0,source=zmel)
                       call x0kf_zmel(q, k, isp_k,isp_kq)!, GPUTEST=GPUTEST)
                       kold=k
                       write(6,*) 'k, mpi__rank_k', k, mpi__rank_k
                       call flush(6)
                    endif
                    !              write(*,*)'zzzzzzzzzzzzzmel',shape(zmel)
                    do iw=iwini(icoun),iwend(icoun) !iw  = iwc(icount)  !omega-bin
                       icount= icouini(icoun)+iw-iwini(icoun)
                       if(abs(zmel0(1,it,itp))>1d10) cycle !We assume rcxq(1) in this mode
                       rcxq(1,1,iw,jpm)=rcxq(1,1,iw,jpm) +rfac00**2*(abs(zmel(1,it,itp))-abs(zmel0(1,it,itp)))**2 *whwc(icount)
                    enddo
                 enddo zmel0modeicount
                 !$acc update device (rcxq)
               endblock zmel0block
               goto 2000 
            endif zmel0mode
            if(cmdopt0('--emptyrun')) goto 1590
!             GPUTEST=.false.
!             if(GPUTEST) then
!                ! rcxq(ibg1,igb2,iw) = \sum_ibib wwk(iw,ibib)* <M_ibg1(q) psi_it(k)| psi_itp(q+k)> < psi_itp | psi_it M_ibg2 > at q
!                call stopwatch_init(t_sw_zmel, 'zmel_'//'gpu') !merge('gpu','ori',mask = GPUTEST))
!                call stopwatch_init(t_sw_x0, 'x0_'//'gpu') !merge('gpu','ori',mask = GPUTEST))
!                kloop:do 1500 k=1,nqbz !zmel = < M(igb q) phi( rk it occ)|  phi(q+rk itp unocc)>
!                   if(mod(k-1, mpi__size_k) /= mpi__rank_k)  cycle
!                   ! qq   = q;              qrk  = q+rk(:,k)
!                   ! ispm = isp_k;          ispq = isp_kq
!                   ! ns1  = nkmin(k)+nctot; ns2  = nkmax(k)+nctot
!                   ! nqini= nkqmin(k);      nqmax= nkqmax(k)
!                   icounkmink= icounkmin(k); icounkmaxk= icounkmax(k)
!                   call stopwatch_start(t_sw_zmel)
!                   debug=cmdopt0('--debugzmel')
!                   if(debug) write(stdo,ftox) 'ggggggggg goto get_zmel_init_gemm',k, nkmin(k),nkmax(k),nctot
!                   if(use_gpu) then
!                      !Currently, mpi version of get_zmel_init_gpu which is available by adding comm argument for MPI communicator,
!                      !but, MPI communication is significant bottle-neck in the case where GPUs are used. Therefore, it is only used in without GPU case.
!                      call get_zmel_init_gemm(q=q+rk(:,k), kvec=q, irot=1, rkvec=q, ns1=nkmin(k)+nctot,ns2=nkmax(k)+nctot, ispm=isp_k, &
!                           nqini=nkqmin(k),nqmax=nkqmax(k), ispq=isp_kq,nctot=nctot, ncc=merge(0,nctot,npm==1),iprx=.false., &
!                           zmelconjg=.true.)
!                   else
!                      call get_zmel_init_gemm(q=q+rk(:,k), kvec=q, irot=1, rkvec=q, ns1=nkmin(k)+nctot,ns2=nkmax(k)+nctot, ispm=isp_k, &
!                           nqini=nkqmin(k),nqmax=nkqmax(k), ispq=isp_kq,nctot=nctot, ncc=merge(0,nctot,npm==1),iprx=.false., &
!                           zmelconjg=.true., comm = comm_b)
!                   endif
!                   call stopwatch_pause(t_sw_zmel)
!                   call stopwatch_start(t_sw_x0)
!                   call x0gpu(rcxq, npr, ipr_col, npr_col, nwhis, npm)
!                   call stopwatch_pause(t_sw_x0)
!                   write(6,ftox) 'end of k:', k ,' of:',nqbz, 'zmel:', ftof(stopwatch_lap_time(t_sw_zmel),4), '(sec)', &
!                        ' x0:', ftof(stopwatch_lap_time(t_sw_x0),4), '(sec)'
!                   call flush(6)
! 1500           enddo kloop
!             else 
               !where(abs(uumat_k)>1d6) uumat_k=0d0 !zeroclear padding part (no data region)
              kloop10:do 1510 k=1,nqbz !6,6 !1,nqbz !6,8 !zmel = < M(igb q) phi( rk it occ)|  phi(q+rk itp unocc)>
                if(mod(k-1, mpi__size_k) /= mpi__rank_k)  cycle
                if(cmdopt0('--emptyrun')) cycle
                if(cmdopt0('--x0test')) then
                  continue
                elseif (cmdopt0('--ahc')) then
                  if (cmdopt0('--AHCMAT')) open(newunit=ifahc,file='AHCMAT-q'//trim(charext(k))//'.isp'//trim(charext(isp_k))//trim(charext(isp_kq)))
                  AHCmatrix: block !rk=qbz k is index of iqbz
                    use m_nvfortran,only: findloc
                    real(8),parameter::epsdege=1d-4 !Ry !degeneracy
                    integer,parameter::ndd=10           !degeneracy
                    integer:: in1,in2,inbb(-1:1,nbb,nband_k),ibb,in1i,in1e,in2i,in2e,ii,ee,ierr,inmi,inmx,nww=4
                    complex(8):: facb(3,nbb), bnn(ndd,nband_k,nbb) 
                    complex(8):: uubb(nbb,nbb,nband_k,nband_k)
                    real(8):: ekxx(1:nband_k,nbb), ekxx0(1:nband_k), norm, cutuu=.3d0
                    logical::cmdopt2
                    character(256)::strn
                    forall(ibb=1:nbb) facb(:,ibb) = bb(:,ibb)*wbb(ibb)
                    !                      do ibb=1,nbb
                    !                        write(stdo,ftox)'fffffffffbbbb',k,ibb,facb(:,ibb)
                    !                      enddo  
                    forall(ibb=1:nbb) ekxx(1:nband_k,ibb) = readeval(rk(:,k)+bb(:,ibb), isp_k )
                    ekxx0(1:nband_k)                      = readeval(rk(:,k),           isp_k )
                    ierr=0
                    if (cmdopt2('--nww=',strn)) read(strn,*) nww
                    if (cmdopt2('--cutuu=',strn)) read(strn,*) cutuu
                    !                      do concurrent(ibb=1:nbb, in=1:nband_k) !For in at k, we look for inbb at k+b. inbb(-1),inbb(1) is lowest and highest.
                    do ibb=1,nbb
                      do in=1,nband_k !For in at k, we look for inbb at k+b. inbb(-1),inbb(1) is lowest and highest.
                        if(ekxx0(in)>1d6) cycle    ! 1d6 for eigenvalue is skip data padding (no data for given in).

                        Bandmatching: block !inbb(0,ibb,in):center 
                          real(8),parameter::epsdege=1d-4 !Ry !degeneracy
                          ! integer:: nww=4  !2,5, Maybe smaller for denser mesh
                          ! real(8):: cutuu=.3d0       !.3   Maybe larger(1 is max) for denser mesh.　
                          if(abs(uumat_k(in,in,ibb,k,isp_k))<cutuu) then
                            inmi=max(in-nww,1)    !search range \pm 1
                            inmx=min(in+nww,nband)
                            !center for ibb,in
                            inbb(0,ibb,in)= inmi-1+ maxloc(abs(uumat_k(in,inmi:inmx,ibb,k,isp_k)), dim=1, & !center  mapping in ---> inbb
                                 mask=                     abs(uumat_k(in,inmi:inmx,ibb,k,isp_k))<1d6     &
                                 .and.                     abs(uumat_k(in,inmi:inmx,ibb,k,isp_k))>1d-8)
                          else
                            inbb(0,ibb,in)= in
                          endif
                          if(debug)write(stdo,ftox)'bbbbbbbbmmmmmm',k,isp_k, ibb,in, ftod(ekxx(inbb(0,ibb,in),ibb),2), ftod(ekxx0(in),2)
                          if(ekxx(inbb(0,ibb,in),ibb)>1d6.or.ekxx0(in)>1d6) then
                            inbb(-1,ibb,in)= 99999
                            cycle
                          endif
                          ! Range of degenerated bands inbb-1 to inbb 1. ekxx(inbb0)-epsdege < ekxx(i,ibb) < ekxx(inbb0)
                          inbb(-1,ibb,in)= findloc( ekxx(:,ibb)-ekxx(inbb(0,ibb,in),ibb)>-epsdege,value=.true.,dim=1)            !-1:left end
                          inbb( 1,ibb,in)= findloc( ekxx(:,ibb)-ekxx(inbb(0,ibb,in),ibb)<epsdege,value=.true.,dim=1,back=.true.) ! 1:right end
                        endblock Bandmatching

                        ii = inbb(-1,ibb,in) 
                        ee = inbb( 1,ibb,in)
                        if(ee-ii+1>ndd) ierr=1
                        bnn(1:ee-ii+1,in,ibb) = dconjg(uumat_k(in,ii:ee,ibb,k,isp_k))
                        norm = sum(abs(bnn(1:ee-ii+1,in,ibb))**2)**.5 !normalization
                        if(norm<1d-8) then
                          bnn(1:ee-ii+1,in,ibb)=0d0
                        else   
                          bnn(1:ee-ii+1,in,ibb)= bnn(1:ee-ii+1,in,ibb)/norm !for given in at k, we have 1:ee-ii+1 degenerated bands at k+bb
                        endif
                        if(debug.and.sum(abs(bnn(1:ee-ii+1,in,ibb)))>1d3) write(stdo,ftox)'bbbbbbbbbsssssss', sum(abs(bnn(1:ee-ii+1,in,ibb)))
                      enddo
                    enddo

                    if(ierr==1) call rx('ee-ii+1>ndd: something wrong or enlarge ndd')
                    allocate(zzahc(3,3,nband_k,nband_k),source=(0d0,0d0))
                    !                      iminloop: do concurrent(im=1:nband_k,in=1:nband_k) ! im occupied, in unoccupied
                    do im=1,nband_k
                      do in=1,nband_k ! im occupied, in unoccupied                        
                        if(ekxx0(im)>1d6 .or. ekxx0(in)>1d6) cycle
                        do ibb1=1,nbb
                          do ibb2=1,nbb !do concurrent(ibb1=1:nbb,ibb2=1:nbb) 
                            if(debug)write(stdo,ftox)'iiiiiiiiimmmmmm in=',im,in,ibb1,ibb2,k,isp_k
                            in1i=inbb(-1,ibb1,in)
                            in1e=inbb( 1,ibb1,in)
                            in2i=inbb(-1,ibb2,in)
                            in2e=inbb( 1,ibb2,in)
                            !write(stdo,ftox)'iiiiiiiiiiiinnnnnnnn',in1i,in2i,in1e,in2e
                            if(in1i==99999.or.in2i==99999) then
                              uubb(ibb1,ibb2,im,in)=0d0
                            else
                              uubb(ibb1,ibb2,im,in)=  sum(uumat_k(im,in1i:in1e,ibb1,k,isp_k)*bnn(1:in1e-in1i+1,in,ibb1)) &
                                   *           dconjg(sum(uumat_k(im,in2i:in2e,ibb2,k,isp_k)*bnn(1:in2e-in2i+1,in,ibb2)))

                              !if(abs(uubb(ibb1,ibb2,im,in))>1d3)
                              if(debug.and.(im==8.and.in==9)) then
                                write(stdo,ftox)'uuuubbrange',ibb2,ibb1,im,in,k,isp_k, 'range=',in1i,in2i,in1e,in2e
                                do i=in1i,in1e
                                  write(stdo,ftox)'uuuuuuuubb1=',im,i,k,isp_k,ibb1,ftod(uumat_k(im,i,ibb1,k,isp_k),2),'r:',ibb2,ibb1,in1i,in2i,in1e,in2e
                                enddo
                                do i=in2i,in2e
                                  write(stdo,ftox)'uuuuuuuubb2=',im,i,k,isp_k,ibb2,ftod(uumat_k(im,i,ibb2,k,isp_k),2),'r:',ibb2,ibb1,in1i,in2i,in1e,in2e
                                enddo
                              endif

                              !if(abs(uubb(ibb1,ibb2,im,in))>1d3)
                              !    write(stdo,ftox)'uuuuuuuubb=',ibb1,ibb2,im,in,k,isp_k,abs(uubb(ibb1,ibb2,im,in))
                            endif
                            if(debug) then
                              if(ibb1==1.and.ibb2==2.and.im<16) then
                                write(stdo,ftox)'sssssss mat',ftof(ekxx0(im)),ftof(ekxx0(in)),'isp_k it itp=',isp_k,k,im,in,&
                                     'abs=',ftod(abs(uumat_k(im,in,ibb1,k,isp_k))),ftod(abs(uumat_k(im,in,ibb2,k,isp_k)))
                              endif
                            endif
                          enddo                           !enddo; enddo
                        enddo                           !enddo; enddo
                        zzahc(:,:,im,in) = matmul(facb(:,:),matmul(uubb(:,:,im,in),transpose(facb(:,:)))) !
                        if((abs(zzahc(1,2,im,in)))>.1d0) write(stdo,ftox)'zzzzsssss=',im,in,k,isp_k&
                             ,ftod(zzahc(1,2,im,in)),ftod(zzahc(2,1,im,in))
                      enddo
                      enddo
                        !                      enddo iminloop
                    endblock AHCmatrix
                  else
                     call get_zmel_init_gemm(q=q+rk(:,k), kvec=q, irot=1, rkvec=q, ns1=nkmin(k)+nctot, &
                          ns2=nkmax(k)+nctot, ispm=isp_k, nqini=nkqmin(k),nqmax=nkqmax(k), &
                          ispq=isp_kq,nctot=nctot,ncc=merge(0,nctot,npm==1),iprx=.false.,zmelconjg=.true.)
                  endif
                  icounloop: do 1000 icoun=icounkmin(k),icounkmax(k)
                     ! call get_zmel_init is equivalent to call x0kf_zmel(q, k, isp_k,isp_kq) 
                     TimeConsumingRcxq: block
                       complex(8):: zmelzmel(npr,npr)
                       real(8) :: zmelzmelahc(3,3)
                       integer:: alpha, beta
                       real(8),parameter:: pi=3.141592653589793
                       if(debug) write(stdo,ftox)'icoun: iq k jpm it itp n(iw)=',icoun,iq,k,jpm,it,itp,iwend(icoun)-iwini(icoun)+1
                       associate( &
                            jpm => jpmc(icoun),&  !\pm omega 
                            it  => itc (icoun),&  !occ      at k
                            itp => itpc(icoun))   !unocc    at q+k

!! NOTE since Im(CC)=0 for n=m, we don't need next skip mechanism for it>=itp.
!! if(cmdopt0('--ahc').and. it>=itp) cycle ! we can not skip it>=itp 2025-2-10
                         
                         if(cmdopt0('--x0test')) then
                            forall(iw=iwini(icoun):iwend(icoun)) rcxq(:,:,iw,jpm)=rcxq(:,:,iw,jpm)+ whwc(iw-iwini(icoun)+icouini(icoun))
                         elseif (cmdopt0('--ahc')) then ! AHC: Anomalous Hall conductivity matrix
                            do concurrent(alpha=1:3,beta=1:3)
!                               zmelzmelahc(alpha,beta) = dimag(zmelahc_k(alpha,it,itp)*zmelahc_k(beta,itp,it) &
!                                    -                    zmelahc_k(beta,it,itp)*zmelahc_k(alpha,itp,it))
                               zmelzmelahc(alpha,beta) = -dimag(zzahc(alpha,beta,it,itp)-zzahc(beta,alpha,it,itp))
                            enddo
                            if (sum(abs(uumat_k(it,itp,:,k,isp_k)))>1d6) cycle !.or.sum(abs(uumat_kq(it,itp,:,k)))>1d6) cycle
!                            if (cmdopt0('--AHCMAT')) then
!                               write(ifahc,ftox) it,itp,ftod(zmelzmelahc(1,2),3),&
!                                    ftod(zmelzmelahc(1,3),3),  ftod(zmelzmelahc(2,3),3),'  k=',ftof(rk(:,k),4)
!                            endif
                            wgt = sum(whwc(icouini(icoun):iwend(icoun)-iwini(icoun)+icouini(icoun)))/pi
                            omegaahc_tet(:,:,k) = omegaahc_tet(:,:,k) + wgt *zmelzmelahc(:,:)  ! tetrahedron sum
!                            omegaahc_sp(:,:,k) = omegaahc_sp(:,:,k) + zmelzmelahc(:,:)  ! direct sum
                            omegasum_tet(:,:,k) = omegasum_tet(:,:,k) + wgt
!                            omegasum_sp(:,:,k) = omegasum_sp(:,:,k) + 1d0
                         else
                            do concurrent(igb1=1:npr,igb2=1:npr) 
                               zmelzmel(igb1,igb2)= dconjg(zmel(igb1,it,itp))*zmel(igb2,it,itp) 
                            enddo
                            forall(iw=iwini(icoun):iwend(icoun))& !rcxq is hermitian, thus, we can reduce computational time half.
                                 rcxq(:,:,iw,jpm)=rcxq(:,:,iw,jpm)+ whwc(iw-iwini(icoun)+icouini(icoun))*zmelzmel(:,:) !may use zaxpy and symmetrize afterwards
                         endif
                         !!forall(iw=iwini(icoun):iwend(icoun)) nwj(iw,jpm)=nwj(iw,jpm)+iwend(icoun)-iwini(icoun)+1 !counter check
                       endassociate
                     endblock TimeConsumingRcxq
1000              enddo icounloop

!direct sum                  
                  directsum: block
                    use m_ReadEfermi,only: ef
                    integer:: ibibc
                    real(8):: ekxxx1(1:nband_k)  ,rydberg !                  real(8):: ekxxx2(1:nband_kq) 
                    ekxxx1(1:nband_k)  = readeval(  rk(:,k), isp_k )   !              ekxxx2(1:nband_kq) = readeval(q+rk(:,k), isp_kq)
                    do it=1, nband_k
                       do itp=1,nband_k
                          if( ekxxx1(it) > ef+efshift/rydberg()) cycle
                          if( ekxxx1(itp)< ef+efshift/rydberg()) cycle
                          if( ekxxx1(it) >1d6) cycle
                          if( ekxxx1(itp)>1d6) cycle
                          ! call get_zmel_init is equivalent to call x0kf_zmel(q, k, isp_k,isp_kq) 
                          TTimeConsumingRcxq: block
                            complex(8):: zmelzmel(npr,npr)
                            real(8) :: zmelzmelahc(3,3)
                            integer:: alpha, beta
                            real(8),parameter:: pi=3.141592653589793
                            if (cmdopt0('--ahc')) then ! AHC: Anomalous Hall conductivity matrix
                               do concurrent(alpha=1:3,beta=1:3)
                                  ! zmelzmelahc(alpha,beta) = dimag(zmelahc(alpha,it,itp)*dconjg(zmelahc(beta,it,itp)) &
                                  !      -                    zmelahc(beta,it,itp)*dconjg(zmelahc(alpha,it,itp)))
                                  !                               zmelzmelahc(alpha,beta) = dimag(zmelahc_k(alpha,it,itp)*zmelahc_kq(beta,itp,it) &
                                  !                                    -                    zmelahc_kq(beta,it,itp)*zmelahc_k(alpha,itp,it))

                                  !zmelzmelahc(alpha,beta) = dimag(zmelahc_k(alpha,it,itp)*zmelahc_k(beta,itp,it) &
                                  !     -                    zmelahc_k(beta,it,itp)*zmelahc_k(alpha,itp,it))
                                  zmelzmelahc(alpha,beta) = -dimag(zzahc(alpha,beta,it,itp)-zzahc(beta,alpha,it,itp))
                               enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                               
                               if(debug.and.abs(zzahc(1,2,it,itp))>1d-10) &
                                    write(stdo,ftox)'mmm mat',ftof(ekxxx1(it)),ftof(ekxxx1(itp)),'it itp=',it,itp,&
                                    'zmel12=',ftof(zzahc(1,2,it,itp)),'abs=',abs(zzahc(1,2,it,itp))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                               
                               if (sum(abs(uumat_k(it,itp,:,k,isp_k)))>1d6) cycle !.or.sum(abs(uumat_kq(it,itp,:,k)))>1d6) cycle
                               if (cmdopt0('--AHCMAT').and.abs(zmelzmelahc(1,2))>1d-1) then
                                  write(ifahc,ftox) it,itp,k,ftod(zmelzmelahc(1,2),2), ftof(ekxxx1(it),3),ftof(ekxxx1(itp),3)!,&
!                                       ftod(zmelzmelahc(1,3),3),  ftod(zmelzmelahc(2,3),3),'  k=',ftof(rk(:,k),4)
                               endif
                               omegaahc_sp(:,:,k) = omegaahc_sp(:,:,k) + zmelzmelahc(:,:)  ! direct sum
                               omegasum_sp(:,:,k) = omegasum_sp(:,:,k) + 1d0
                            endif
                          endblock TTimeConsumingRcxq
                       enddo
                    enddo
                  endblock directsum
!----------------------------------
!                  if(cmdopt0('--ahc')) deallocate(zmelahc_k)
                  if(cmdopt0('--ahc')) deallocate(zzahc)
                  if(cmdopt0('--ahc').and.cmdopt0('--AHCMAT')) close(ifahc)
1510           enddo kloop10
!            endif
1590        continue
2000        continue
            if(cmdopt0('--x0test')) then
               open(newunit=ifx0,file='x0test.'//charnum7(iq)//'.dat') !iq is after nqibz
               iwloop: do iw= 1,nwhis
                  write(ifx0,ftox) ftof(q),' ',iw,ftod(2d0*frhis(iw)),ftod(2d0*frhis(iw+1)),' ',rcxq(1,1,iw,1) !frhis in a.u.
               enddo iwloop
            endif
            if (cmdopt0('--ahc')) then
               sigAHC:block ! AHC: calculate anomalous Hall conductivity
                 integer:: i,k,iahc_tet,iahc_sp,isum,isum_tet,isum_sp
                 real(8),parameter:: fac=4.599849969441961d4 ! (e^2/hbar)/a_0 / (Ohm^-1 cm^-1), where a_0=1 bohr
                 real(8)::sigmaahc_tet(3,3),sigmaahc_sp(3,3),omegaahcsum_tet(3,3),omegaahcsum_sp(3,3),kw
                 real(8)::summ_tet, summ_sp
                 logical::lahc,cmdopt0
                 character(1) :: charisp_k, charisp_kq
                 character(256)::dum1,dum2

                 sigmaahc_tet = 0d0
                 sigmaahc_sp = 0d0
                 write(stdo,*) 'AHC: calculate sigma'
                 write(stdo,*) 'Efermi shift (eV): ', efshift
                 omegaahcsum_tet = 0d0
                 omegaahcsum_sp = 0d0
                 summ_tet = 0d0
                 summ_sp = 0d0
                 do k=1,nqbz
                    omegaahcsum_tet(:,:) = omegaahcsum_tet(:,:) + omegaahc_tet(:,:,k)
                    omegaahcsum_sp(:,:) = omegaahcsum_sp(:,:) + wk(k)*omegaahc_sp(:,:,k)
                    summ_tet = summ_tet + omegasum_tet(1,1,k)
                    summ_sp = summ_sp + wk(k)*omegasum_sp(1,1,k)
                 enddo
                 do i=1,3
                    call MPI__AllreduceSumReal(omegaahcsum_tet(i,:),3,communicator=comm_k)
                    call MPI__AllreduceSumReal(omegaahcsum_sp(i,:),3,communicator=comm_k)
                 enddo
                 call MPI__AllreduceSumRealSca(summ_tet,communicator=comm_k)
                 call MPI__AllreduceSumRealSca(summ_sp,communicator=comm_k)
                 sigmaahc_tet(:,:) = -fac * omegaahcsum_tet(:,:) / vol
                 sigmaahc_sp(:,:) = -fac * omegaahcsum_sp(:,:) / vol

                 write(charisp_k,'(i1.1)') isp_k
                 write(charisp_kq,'(i1.1)') isp_kq
                 if(mpi__root_k) then
                    inquire(file='ahc_tet.isp'//trim(charisp_k)//trim(charisp_kq)//'.dat',exist=lahc)
                    if (.not. lahc) then
                       open(newunit=iahc_tet,file='ahc_tet.isp'//trim(charisp_k)//trim(charisp_kq)//'.dat',action='write',status='new')
                       write(iahc_tet,*) '### sigma(alpha,beta) in (Ohm^-1 cm^-1)  alpha,beta=x,y,z '
                       write(iahc_tet,*) '### 1:dEf  2:xx  3:xy  4:xz  5:yx  6:yy  7:yz  8:zx  9:zy  10:zz '
                       open(newunit=iahc_sp,file='ahc_sp.isp'//trim(charisp_k)//trim(charisp_kq)//'.dat',action='write',status='new')
                       write(iahc_sp,*) '### sigma(alpha,beta) in (Ohm^-1 cm^-1)  alpha,beta=x,y,z '
                       write(iahc_sp,*) '### 1:dEf  2:xx  3:xy  4:xz  5:yx  6:yy  7:yz  8:zx  9:zy  10:zz '
                       open(newunit=isum,file='sum.isp'//trim(charisp_k)//trim(charisp_kq)//'.dat',action='write',status='new')
                       write(isum,*) '### sigma(alpha,beta) in (Ohm^-1 cm^-1)  alpha,beta=x,y,z '
                       write(isum,*) '### 1:dEf  2:tet  3:direct '
                    else
                       open(newunit=iahc_tet,file='ahc_tet.isp'//trim(charisp_k)//trim(charisp_kq)//'.dat', action='write',position='append',status='old')
                       open(newunit=iahc_sp,file='ahc_sp.isp'//trim(charisp_k)//trim(charisp_kq)//'.dat',   action='write',position='append',status='old')
                       open(newunit=isum,file='sum.isp'//trim(charisp_k)//trim(charisp_kq)//'.dat',         action='write',position='append',status='old')
                    endif
                    write(iahc_tet,'(f10.7,x,9e17.9)') efshift,sigmaahc_tet(1:3,1:3)
                    close(iahc_tet)
                    write(iahc_sp,'(f10.7,x,9e17.9)') efshift,sigmaahc_sp(1:3,1:3)
                    close(iahc_sp)
                    write(isum,'(f10.7,x,2e17.9)') efshift,summ_tet,summ_sp
                    close(isum)
                 endif
                 if (isp_k==nspx .and. isp_kq==nspx) call rx0('AHC: Finish calculations')
               endblock sigAHC
            endif
            write(stdo,ftox)"--- x0kf_ahc: end",sum(abs(rcxq(:,:,1:nwhis,1:npm))) !,sum((rcxq(:,:,1:nwhis,1:npm)))
          endblock x0kf_v4hz_block
          deallocate(whwc, kc, iwini,iwend, itc,itpc, jpmc,icouini, nkmin,nkmax,nkqmin,nkqmax,icounkmin,icounkmax)
          ! HilbertTransformation:if(isp_k==nspx .OR. chipm) then
          !    !Get real part. When chipm=T, do dpsion5 for every isp_k; When =F, do dpsion5 after rxcq accumulated for spins
          !    !        if(GPUTEST) then
          !    !$acc exit data copyout(rcxq)
          !    !        endif
          !    mpi_kaccumulate: block
          !      integer :: jpm, iw
          !      !To reduce memory allocation size in MPI__reduceSUM, mpi_reduce operation has been split for each iw and jpm
          !      do jpm=1, npm
          !         do iw=1, nwhis
          !            call MPI__reduceSum(0, rcxq(1,1,iw,jpm), npr*npr_col, communicator = comm_k) 
          !         enddo
          !      enddo
          !    end block mpi_kaccumulate
          !    if(mpi__root_k.and..not.cmdopt0('--ahc')) then
          !       call dpsion5(realomega, imagomega, rcxq, npr, npr_col, zxq, zxqi, chipm, schi,isp_k,  ecut,ecuts)
          !    endif
          !    if(cmdopt0('--ahc')) then
          !       deallocate(rcxq,omegaahc_tet,omegaahc_sp,omegasum_tet,omegasum_sp)
          !    else
          !       deallocate(rcxq)
          !    endif
          ! endif HilbertTransformation
1104   enddo isqloop
1103 enddo isloop
!     call ExitDataGPU_inkx()
  end subroutine x0kf_ahc
    
  subroutine deallocatezxq()
    deallocate(zxq)
  end subroutine deallocatezxq
  subroutine deallocatezxqi()
    deallocate(zxqi)
  end subroutine deallocatezxqi
  subroutine x0kf_zmel( q,k, isp_k,isp_kq)!, GPUTEST) ! Return zmel= <phi phi |M_I> in m_zmel
    use m_mpi, only: comm_b
    intent(in)   ::     q,k, isp_k,isp_kq   
    integer::              k,isp_k,isp_kq 
    real(8)::           q(3)
    debug=cmdopt0('--debugzmel')
!    logical, intent(in), optional:: GPUTEST
!    call get_zmel_init(q=q+rk(:,k), kvec=q, irot=1, rkvec=q, nm1=nkmin(k)+nctot, nm2=nkmax(k)+nctot, ispm=isp_k, &
!         nqini=nkqmin(k), nqmax=nkqmax(k), ispq=isp_kq,nctot=nctot, ncc=merge(0,nctot,npm==1), iprx=.false., zmelconjg=.true.)
!    if (present(GPUTEST)) then
    !      if (GPUTEST) then
    if(debug) write(stdo,ftox) 'ggggggggg goto get_zmel_init_gemm',k, nkmin(k),nkmax(k),nctot
    call get_zmel_init_gemm(q=q+rk(:,k), kvec=q, irot=1, rkvec=q, ns1=nkmin(k)+nctot,ns2=nkmax(k)+nctot, ispm=isp_k, &
         nqini=nkqmin(k),nqmax=nkqmax(k), ispq=isp_kq,nctot=nctot, ncc=merge(0,nctot,npm==1),iprx=.false., zmelconjg=.true.)
       !$acc update host(zmel)
!      endif
!    else
!      call get_zmel_init(q=q+rk(:,k), kvec=q, irot=1, rkvec=q, ns1=nkmin(k)+nctot, ns2=nkmax(k)+nctot, ispm=isp_k, &
!           nqini=nkqmin(k), nqmax=nkqmax(k), ispq=isp_kq,nctot=nctot, ncc=merge(0,nctot,npm==1), iprx=.false., zmelconjg=.false.)
!    endif
  end subroutine x0kf_zmel
end module m_x0kf_ahc

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

subroutine readuun(ikbidx,iko_ixs,iko_fxs,nband,nqbz,nbb,nspin,uum)
  use m_lgunit,only:stdo
  use m_ftox
  implicit none
  integer:: iqbz,nband,nqbz,nbb,nspin,ifuu,ibb,isp,iko_ixs(2),iko_fxs(2),ii,ie
  complex(8) :: uum(nband,nband,nbb,nqbz,nspin)
  integer:: ikbidx(nbb,nqbz)
  do isp=1,nspin
    ii=iko_ixs(isp)
    ie=iko_fxs(isp)
    call readuu(isp,ii,ie,ikbidx,nqbz,nbb,uum(ii:ie,ii:ie,:,:,isp))
  enddo
end subroutine readuun
