!> get zxq and zxqi for given q
module m_x0kf
  use m_lgunit,only: stdo
  use m_keyvalue,only : Getkeyvalue
  use m_pkm4crpa,only : Readpkm4crpa
  use m_zmel,only: Get_zmel_init,zmel !,get_zmel_init1,get_zmel_init2
  use m_freq,only: npm, nwhis
  use m_genallcf_v3,only:  nsp=>nspin ,nlmto,nctot
  use m_read_bzdata,only:  nqbz,ginv,nqibz,  rk=>qbz,wk=>wbz
  use m_rdpp,only: nbloch
  use m_readqg,only: ngpmx,ngcmx
  use m_qbze,only: nqbze
  use m_readhbe,only: nband
  use m_hamindex,only: ngrp
  use m_tetwt,only:  gettetwt,tetdeallocate, whw,ihw,nhw,jhw,n1b,n2b,nbnb,nbnbx,nhwtot
  use m_ftox
  use m_readVcoud,only:   vcousq,zcousq,ngb,ngc
  use m_kind,only:kindrcxq
  implicit none
  public:: x0kf_zxq, deallocatezxq, deallocatezxqi
  complex(8),public,allocatable:: zxq(:,:,:), zxqi(:,:,:)   !Not yet protected because of main_hx0fp0
  complex(kindrcxq),allocatable:: rcxq(:,:,:,:)
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
  
  subroutine x0kf_zxq(realomega,imagomega, q,iq,nprin,schi,crpa,chipm,nolfco,q00,zzr)
    use m_readgwinput,only: ecut,ecuts
    use m_dpsion,only: dpsion5
    use m_freq,only: nw_i, nw_w=>nw, niwt=>niw
    use m_readeigen,only:readeval
    use m_freq,only: nw_i,nw,niw 
    use m_zmel,only: Setppovlz,Setppovlz_chipm   ! & NOTE: these data set are stored in this module, and used
    implicit none
    intent(in)::      realomega,imagomega, q,iq,nprin,schi,crpa,chipm,nolfco,q00,zzr
    logical:: realomega,imagomega,crpa,chipm,nolfco
    integer:: iq, isp_k,isp_kq,ix0,is,isf,kx,ierr,nprin,k
    real(8),optional:: q00(3)
    complex(8),optional:: zzr(:,:)
    real(8):: q(3),schi,ekxx1(nband,nqbz),ekxx2(nband,nqbz)
    character(10) :: i2char
    logical:: cmdopt0,GPUTEST
    real(8) :: t1, t2
    npr=nprin
    qq=q
    GPUTEST = cmdopt0('--gpu')
    if(realomega) allocate(zxq(npr,npr,nw_i:nw),source=(0d0,0d0))
    if(imagomega) allocate(zxqi(npr,npr,niw),source=(0d0,0d0))
    if(cmdopt0('--emptyrun'))  return
    if(chipm .AND. nolfco) then;  call setppovlz_chipm(zzr,npr)
    else;                         call Setppovlz(q,matz=.true.)
    endif
    isloop: do 1103 isp_k = 1,nsp
      GETtetrahedronWeight:block
        isp_kq = merge(3-isp_k,isp_k,chipm) 
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
        integer:: k,jpm, ibib, iw,igb2,igb1,it,itp, nkmax1,nkqmax1, ib1, ib2, ngcx,ix,iy,igb
        integer:: izmel,nmtot,nqtot,ierr,iwmax,ifi0,icoucold,icoun
        integer:: nwj(nwhis,npm),imb, igc ,neibz,icc,ig,ikp,i,j,itimer,icount, kold 
        real(8):: imagweight, wpw_k,wpw_kq,qa,q0a 
        complex(8):: img=(0d0,1d0)
        logical,parameter:: debug=.false.
        if(.not.allocated(rcxq)) then
           allocate( rcxq(npr,npr,nwhis,npm))
           rcxq=0d0
           if(GPUTEST) then
             write(stdo,ftox)'size of rcxq:', npr, npr, nwhis, npm
             !$acc enter data create(rcxq) 
             !$acc kernels
             rcxq(1:npr,1:npr,1:nwhis,1:npm) = (0d0,0d0)
             !$acc end kernels
           endif
        endif
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
              if(kold/=k) then
                call x0kf_zmel(q00,k, isp_k,isp_kq)
                if(allocated(zmel0)) deallocate(zmel0)
                allocate(zmel0,source=zmel)
                call x0kf_zmel(q, k, isp_k,isp_kq)
                kold=k
              endif
              do iw=iwini(icoun),iwend(icoun) !iw  = iwc(icount)  !omega-bin
                icount= icouini(icoun)+iw-iwini(icoun)
                if(abs(zmel0(1,it,itp))>1d10) cycle !We assume rcxq(1) in this mode
                rcxq(1,1,iw,jpm)=rcxq(1,1,iw,jpm) +rfac00**2*(abs(zmel(1,it,itp))-abs(zmel0(1,it,itp)))**2 *whwc(icount)
              enddo
            enddo zmel0modeicount
          endblock zmel0block
          goto 2000 
        endif zmel0mode
        if(cmdopt0('--emptyrun')) goto 1590

        call cputid (0)
        if(GPUTEST) then
          ! rcxq(ibg1,igb2,iw) = \sum_ibib wwk(iw,ibib)* <M_ibg1(q) psi_it(k)| psi_itp(q+k)> < psi_itp | psi_it M_ibg2 > at q
          kloop:do 1500 k=1,nqbz !zmel = < M(igb q) phi( rk it occ)|  phi(q+rk itp unocc)>
            qq   = q;              qrk  = q+rk(:,k)
            ispm = isp_k;          ispq = isp_kq
            ns1  = nkmin(k)+nctot; ns2  = nkmax(k)+nctot
            nqini= nkqmin(k);      nqmax= nkqmax(k)
            icounkmink= icounkmin(k); icounkmaxk= icounkmax(k)
            call x0gpu(rcxq,npr,nwhis,npm)
1500      enddo kloop
        else ! NOTE: kloop10:do 1510 is equivalent to do 1500. 2024-3-25
          kloop10:do 1510 k=1,nqbz !zmel = < M(igb q) phi( rk it occ)|  phi(q+rk itp unocc)>
            if(cmdopt0('--emptyrun')) cycle
            call cpu_time(t1)
            call get_zmel_init(q=q+rk(:,k), kvec=q, irot=1, rkvec=q, ns1=nkmin(k)+nctot,ns2=nkmax(k)+nctot, ispm=isp_k, &
                 nqini=nkqmin(k),nqmax=nkqmax(k), ispq=isp_kq,nctot=nctot, ncc=merge(0,nctot,npm==1),iprx=.false.,zmelconjg=.true.)
            call cpu_time(t2)
            print '(1x,A,F10.6)','zmel:cpu' ,(t2-t1)
            call cpu_time(t1)
            icounloop: do 1000 icoun=icounkmin(k),icounkmax(k)
              ! call get_zmel_init is equivalent to call x0kf_zmel(q, k, isp_k,isp_kq) 
              TimeConsumingRcxq: block 
                complex(8):: zmelzmel(npr,npr)
                if(debug) write(stdo,ftox)'icoun: iq k jpm it itp n(iw)=',icoun,iq,k,jpm,it,itp,iwend(icoun)-iwini(icoun)+1
                associate( &
                     jpm => jpmc(icoun),&  !\pm omega 
                     it  => itc (icoun),&  !occ      at k
                     itp => itpc(icoun))   !unocc    at q+k
                  do concurrent(igb1=1:npr,igb2=1:npr) 
                    zmelzmel(igb1,igb2)= dconjg(zmel(igb1,it,itp))*zmel(igb2,it,itp) 
                  enddo
                  forall(iw=iwini(icoun):iwend(icoun))& !rcxq is hermitian, thus, we can reduce computational time half.
                       rcxq(:,:,iw,jpm)=rcxq(:,:,iw,jpm)+ whwc(iw-iwini(icoun)+icouini(icoun))* zmelzmel(:,:) !may use zaxpy and symmetrize afterwards
                  !forall(iw=iwini(icoun):iwend(icoun)) nwj(iw,jpm)=nwj(iw,jpm)+iwend(icoun)-iwini(icoun)+1 !counter check
                endassociate
              endblock TimeConsumingRcxq
1000        enddo icounloop
            call cpu_time(t2)
            print '(1x,A,F10.6)','x0:cpu' ,(t2-t1)
            call flush(stdo)
1510      enddo kloop10
        endif
        call cputid (0)
1590    continue
2000   continue
        write(stdo,ftox)"--- x0kf_v4hz: end",sum(abs(rcxq(:,:,1:nwhis,1:npm))) !,sum((rcxq(:,:,1:nwhis,1:npm)))
      endblock x0kf_v4hz_block
      deallocate(whwc, kc, iwini,iwend, itc,itpc, jpmc,icouini, nkmin,nkmax,nkqmin,nkqmax,icounkmin,icounkmax)
      HilbertTransformation:if(isp_k==nsp .OR. chipm) then
        !Get real part. When chipm=T, do dpsion5 for every isp_k; When =F, do dpsion5 after rxcq accumulated for spins
        if(GPUTEST) then
          !$acc exit data copyout(rcxq)
        endif
        call dpsion5(realomega, imagomega, rcxq, npr,npr, zxq, zxqi, chipm, schi,isp_k,  ecut,ecuts) 
        deallocate(rcxq)
      endif HilbertTransformation 
1103 enddo isloop
  end subroutine x0kf_zxq
  subroutine deallocatezxq()
    deallocate(zxq)
  end subroutine deallocatezxq
  subroutine deallocatezxqi()
    deallocate(zxqi)
  end subroutine deallocatezxqi
  subroutine x0kf_zmel( q,k, isp_k,isp_kq) ! Return zmel= <phi phi |M_I> in m_zmel
    intent(in)   ::     q,k, isp_k,isp_kq   
    integer::              k,isp_k,isp_kq 
    real(8)::           q(3)
    call get_zmel_init(q=q+rk(:,k), kvec=q, irot=1, rkvec=q, ns1=nkmin(k)+nctot, ns2=nkmax(k)+nctot, ispm=isp_k, &
         nqini=nkqmin(k), nqmax=nkqmax(k), ispq=isp_kq,nctot=nctot, ncc=merge(0,nctot,npm==1), iprx=.false., zmelconjg=.true.)
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
!! nlmto   = total number of atomic basis functions within MT
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
