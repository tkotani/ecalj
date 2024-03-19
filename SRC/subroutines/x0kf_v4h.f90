module m_x0kf
  use m_lgunit,only: stdo
  use m_keyvalue,only : Getkeyvalue
  use m_pkm4crpa,only : Readpkm4crpa
  use m_zmel,only: Get_zmel_init,zmel 
  use m_freq,only: npm, nwhis
  use m_genallcf_v3,only:  nsp=>nspin ,nlmto,nctot
  use m_read_bzdata,only:  nqbz,ginv,nqibz,  rk=>qbz,wk=>wbz
  use m_rdpp,only: nbloch
  use m_readqg,only: ngpmx,ngcmx
  use m_qbze,only: nqbze
  use m_readhbe,only: nband
  use m_hamindex,only: ngrp
  use m_tetwt,only:     whw,ihw,nhw,jhw,n1b,n2b,nbnb,nbnbx,nhwtot
  use m_ftox
  use m_readVcoud,only:   Readvcoud,vcousq,zcousq,ngb,ngc
  implicit none
  public:: X0kf_v4hz_init, X0kf_v4hz,x0kf_zmel,X0kf_v4hz_init_write,X0kf_v4hz_init_read,DeallocateRcxq
  integer,public,allocatable,protected:: kc(:)
  integer,public,protected:: ncount,ncoun
  complex(8),allocatable,public:: rcxq(:,:,:,:) !rcxq(:,:,:) !main output owe did not protect, because this is destroyed in dpsion
  private 
  integer,allocatable:: nkmin(:), nkmax(:),nkqmin(:),nkqmax(:) 
  real(8),allocatable:: whwc(:)
  integer,allocatable:: iwini(:),iwend(:),itc(:),itpc(:),jpmc(:),icouini(:)
!  complex(8),allocatable :: zmel0(:,:,:) 
contains
  subroutine DeallocateRcxq()
    deallocate(rcxq)
  end subroutine DeallocateRcxq
  subroutine X0kf_v4hz_init_write(iq,is)
    integer:: iq,is,ix0
    character(10) :: i2char
    open(newunit=ix0,file='x0icount.'//trim(i2char(iq))//'_'//trim(i2char(is)),form='unformatted')
    write(ix0)ncount,ncoun
    write(ix0) whwc, kc, iwini,iwend, itc,itpc, jpmc,icouini
    write(ix0) nkmin,nkmax,nkqmin,nkqmax
    deallocate(nkmin,nkmax,nkqmin,nkqmax)!,skipk)
    deallocate( whwc, kc, iwini,iwend, itc, itpc, jpmc,icouini )
    close(ix0)
  end subroutine X0kf_v4hz_init_write
  subroutine X0kf_v4hz_init_read(iq,is)
    integer:: iq,is,ix0
    character(10) :: i2char
    open(newunit=ix0,file='x0icount.'//trim(i2char(iq))//'_'//trim(i2char(is)),form='unformatted')
    read(ix0)ncount,ncoun
    if(allocated(whwc)) deallocate( whwc, kc, iwini,iwend, itc, itpc, jpmc, nkmin,nkmax,nkqmin,nkqmax)
    allocate( whwc(ncount), kc(ncoun), iwini(ncoun),iwend(ncoun),itc(ncoun),itpc(ncoun), jpmc(ncoun),&
         icouini(ncoun) )
    allocate( nkmin(nqbz),nkmax(nqbz),nkqmin(nqbz),nkqmax(nqbz))
    read(ix0) whwc, kc, iwini,iwend, itc,itpc, jpmc,icouini
    read(ix0) nkmin,nkmax,nkqmin,nkqmax
    close(ix0)
  end subroutine X0kf_v4hz_init_read
  function X0kf_v4hz_init(job,q,isp_k,isp_kq, iq, nmbas, crpa) result(ierr) !index accumulation. Initialzation for calling x0kf_v4h
    intent(in)::          job,q,isp_k,isp_kq, iq, nmbas, crpa
    !! Get ncount index for x0kf_v4h in the private variables
    !! Call job=0 and job=1 successively. 
    !!   ncount, ngb and nqibz are the key to estimate computational efforts.
    !! See x0kf_v4hz. We have OUTPUT to private variables nkmin ... jpmc.
    integer::  irot=1,ierr,k,isp_k,isp_kq, iq, jpm, ibib, iw,igb2,igb1,it,itp,job,icount, nmbas,ncc,icoun
    real(8):: q(3), imagweight, wpw_k, wpw_kq
    logical ::  crpa ,showicount=.false.
    write(stdo,'(" x0kf_v4hz_init: job q =",i3,3f8.4)') job,q
    ncc=merge(0,nctot,npm==1)
    if(job==0) allocate( nkmin(nqbz),nkqmin(nqbz),source= 999999)
    if(job==0) allocate( nkmax(nqbz),nkqmax(nqbz),source=-999999)
    if(job==1) allocate( whwc(ncount),kc(ncoun),iwini(ncoun),iwend(ncoun),itc(ncoun), itpc(ncoun), jpmc(ncoun),&
         icouini(ncoun))
    icount=0
    icoun=0
    AccumulateIndex4icount: do 110 k = 1,nqbz 
       if(job==0) then
          do jpm=1,npm         !npm
             do ibib = 1, nbnb(k,jpm)
                nkmin(k)  = min(n1b(ibib,k,jpm),nkmin(k))
                nkqmin(k) = min(n2b(ibib,k,jpm),nkqmin(k))
                if(n1b(ibib,k,jpm)<=nband) nkmax(k)  = max(n1b(ibib,k,jpm),nkmax(k))
                if(n2b(ibib,k,jpm)<=nband) nkqmax(k) = max(n2b(ibib,k,jpm),nkqmax(k))
             enddo
          enddo
       endif
       if(npm==2.AND.nkqmin(k)/=1)call rx( " When npm==2, nkqmin==1 should be.")! write(stdo,*)' nkqmin nkqmax  nkmin nkmax=',nkqmin,nkqmax,nkmin,nkmax
       jpmloop: do 1251 jpm  = 1, npm ! npm=1 usually (or =2)
          ibibloop: do 125 ibib = 1, nbnb(k,jpm) !---  ibib loop, ibib is decomposed to band index pair, it and itp
             !! n1b,n2b --> core after valence.  it,itp --> valence after core 
             if(ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1 >nwhis) call rx( "x0kf_v4hz: iw>nwhis") ! sanity check
             if(n1b(ibib,k,jpm) > nkmax(k) ) cycle
             if(n2b(ibib,k,jpm) > nkqmax(k)) cycle
             it =merge(nctot +n1b(ibib,k,jpm),             n1b(ibib,k,jpm)-nband,             n1b(ibib,k,jpm)<=nband) !bandindex val or core
             itp=merge(ncc   +n2b(ibib,k,jpm)-nkqmin(k)+1, n2b(ibib,k,jpm)-nkqmin(k)+1-nband, n2b(ibib,k,jpm)<=nband) !bandindex val or core
             if(crpa) then ! constraint RPA mode (juelich verison)
                wpw_k =merge(readpkm4crpa(n1b(ibib,k,jpm),   rk(:,k), isp_k), 0d0, n1b(ibib,k,jpm)<=nband)
                wpw_kq=merge(readpkm4crpa(n2b(ibib,k,jpm), q+rk(:,k), isp_kq),0d0, n2b(ibib,k,jpm)<= nband)
             endif
             icoun=icoun+1
             if(job==1) then
                kc   (icoun) = k
                itc  (icoun)= it
                itpc (icoun)= itp
                jpmc (icoun)= jpm
                iwini(icoun) = ihw(ibib,k,jpm)
                iwend(icoun) = ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1
                icouini(icoun)=icount+1
!                icoend(icoun) =icount+ iwend(icoun)-iwend(icoun)+1
             endif
             iwloop: do iw = ihw(ibib,k,jpm),ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1 !iiww=iw+ihw(ibib,k)-1
                imagweight = whw(jhw(ibib,k,jpm)+iw-ihw(ibib,k,jpm))
                if(crpa) imagweight = imagweight*(1d0-wpw_k*wpw_kq) ! if(eibzmode) imagweight = nwgt(k)*imagweight
                icount = icount+1
                if(job==1) then
                   whwc (icount)= imagweight
!                   iwc  (icount)= iw
!                   icouc(icount)= icoun
                   if(showicount)write(stdo,ftox)'x0kf_v4hz_init: icount jpm k it itp iw jpm whw=',&
                        icount,jpm, k,it,itp,iw, ftof(whwc(icount))
                endif
             enddo iwloop
125       enddo ibibloop
1251   enddo jpmloop
110 enddo AccumulateIndex4icount
    ncount = icount
    ncoun  = icoun
    ierr=0
    if(job==0) write(stdo,"('x0kf_v4hz_init: job=0 ncount ncoun nqibz=',3i8)") ncount, ncoun,nqibz
  endfunction x0kf_v4hz_init
  
  subroutine X0kf_v4hz (q, isp_k,isp_kq, iq, npr,   q00,chipm,nolfco,zzr,nmbas)
    use m_zmel,only: Setppovlz,Setppovlz_chipm   ! & NOTE: these data set are stored in this module, and used
    intent(in)   ::     q, isp_k,isp_kq, iq,       q00,chipm,nolfco,zzr,nmbas !q00 is optional
    logical,optional:: chipm,nolfco
    real(8),optional::q00(3)
    integer,optional::nmbas
    integer:: k,isp_k,isp_kq,iq, jpm, ibib, iw,igb2,igb1,it,itp, npr, nkmax1,nkqmax1, ib1, ib2, ngcx,ix,iy,igb, irot=1
    integer:: izmel,ispold,nmtot,nqtot ,ispold2,ierr,iwmax,ifi0,icoucold,icoun
    integer:: nwj(nwhis,npm),imb, nadd(3), igc ,neibz,icc,ig,ikp,i,j,itimer,icount, kold 
    real(8):: q(3),imagweight, wpw_k,wpw_kq,qa,q0a !    complex(8):: rcxq (npr*(npr+1)/2,nwhis,npm) !only upper-right of rcxq
    complex(8):: img=(0d0,1d0),zmelt2!, zmelzmel(npr*(npr+1)/2)
    complex(8),optional:: zzr(:,:)
    complex(8),allocatable :: zmelt(:,:,:)
    logical :: iww2=.true., cmdopt0,emptyrun,nolfcc, localfieldcorrectionllw
    logical,parameter:: debug=.false.
    emptyrun=cmdopt0('--emptyrun')
    nolfcc=nolfco
    if(chipm .AND. nolfcc) then
       call setppovlz_chipm(zzr,npr) !!!!!!!! zzr??? optional chipm,nolfco zzr(1:nbloch,1:npr)
    else
       call Setppovlz(q,matz=.true.)
    endif
    if(.not. allocated(rcxq)) allocate( rcxq(npr,npr,nwhis,npm),source=(0d0,0d0)) !allocate( rcxq(npr*(npr+1)/2,nwhis,npm),source=(0d0,0d0))
    
    zmel0mode: if(cmdopt0('--zmel0')) then !this is for epsPP0. Use zmel-zmel0 (for subtracting numerical error) for matrix elements.
      zmel0block : block
        real(8)::  q1a,q2a,rfac00
        complex(8),allocatable:: zmel0(:,:,:)
        kold = -999 
        q1a=sum(q00**2)**.5
        q2a=sum(q**2)**.5
        rfac00=q2a/(q2a-q1a)
        zmel0modeicount: do icoun = 1,ncoun 
          k   = kc(icoun)
          it  = itc (icoun)  !occ      k
          itp = itpc(icoun) !unocc    q+k
          jpm = jpmc(icoun) ! \pm omega. Usual mode is only for jpm=1
          if(kold/=k) then
            call x0kf_zmel(q00,k, isp_k,isp_kq)
            if(allocated(zmel0)) deallocate(zmel0)
            allocate(zmel0,source=zmel)
            call x0kf_zmel(q, k, isp_k,isp_kq)
            kold=k
          endif
          do iw=iwini(icoun),iwend(icoun) !!iw  = iwc(icount)  !omega-bin
            icount= icouini(icoun)+iw-iwini(icoun)
            if(abs(zmel0(1,it,itp))>1d10) cycle !We assume rcxq(1) in this mode
            rcxq(1,1,iw,jpm)=rcxq(1,1,iw,jpm) +rfac00**2*(abs(zmel(1,it,itp))-abs(zmel0(1,it,itp)))**2 *whwc(icount)
          enddo
        enddo zmel0modeicount
      endblock zmel0block
      goto 2000 
    endif zmel0mode
    
    kold = -999
    icoucold=-999    ! rcxq(ibg1,igb2,iw) = sum_ibib wwk(iw,ibib)* <M_ibg1(q) psi_it(k)| psi_itp(q+k)> < psi_itp | psi_it M_ibg2 >
    mainloop4rcxqsum: do 1000 icoun=1,ncoun ! = 1,ncount
      if(emptyrun) then
        write(stdo,ftox)'icoun: iq k jpm it itp n(iw)=',icoun,iq,k,jpm,it,itp,iwend(icoun)-iwini(icoun)+1
        cycle
      endif
      k = kc(icoun)
      if(kold/=k) then !get matrix element: zmel for k
        call x0kf_zmel(q, k, isp_k,isp_kq) !Return zmel(igb q,  k it occ,   q+k itp unocc)
        kold=k
      endif
      associate( &
           jpm => jpmc(icoun),&  !\pm omega 
           it  => itc (icoun),&  !occ      at k
           itp => itpc(icoun))   !unocc    at q+k
        block !######### time-consuming part 
          complex(8):: zmelzmel(npr,npr)
          do concurrent(igb1=1:npr,igb2=1:npr) !upper-light block of zmel*zmel
            zmelzmel(igb1,igb2)= dconjg(zmel(igb1,it,itp))*zmel(igb2,it,itp) !right-upper half
          enddo
          forall(iw=iwini(icoun):iwend(icoun)) rcxq(:,:,iw,jpm)=rcxq(:,:,iw,jpm)+whwc(iw-iwini(icoun)+icouini(icoun))*zmelzmel(:,:) !zaxpy
          !  forall(iw=iwini(icoun):iwend(icoun)) nwj(iw,jpm)=nwj(iw,jpm)+iwend(icoun)-iwini(icoun)+1 !onlyfor check counter
        endblock
      endassociate
1000 enddo mainloop4rcxqsum
2000 continue
    deallocate(nkmin,nkmax,nkqmin,nkqmax,whwc,kc,itc,itpc,iwini,iwend,jpmc,icouini)
    write(stdo,ftox) " --- x0kf_v4hz: end" !write(stdo,"(' --- ', 3d13.5)")sum(abs(rcxq(:,:,1:nwhis,1:npm))),sum((rcxq(:,:,1:nwhis,1:npm)))
  end subroutine x0kf_v4hz

  
  subroutine X0kf_zmel( q,k, isp_k,isp_kq) ! Return zmel= <phi phi |M_I> in m_zmel
    intent(in)   ::     q,k, isp_k,isp_kq   
    integer:: k,isp_k,isp_kq,iq 
    real(8):: q(3)
    call get_zmel_init(q=q+rk(:,k), kvec=q, irot=1, rkvec=q, ns1=nkmin(k)+nctot, ns2=nkmax(k)+nctot, ispm=isp_k, &
         nqini=nkqmin(k), nqmax=nkqmax(k), ispq=isp_kq,nctot=nctot, ncc=merge(0,nctot,npm==1), iprx=.false., zmelconjg=.true.)
  end subroutine x0kf_zmel
end module m_x0kf !subroutine x0kf_vhz_symmetrize here removed. See old source code ~2022.

!    intent(out)  ::                                  rcxq
    !! === calculate chi0, or chi0_pm === ! eibzmode, 
    !! We calculate imaginary part of chi0 along real axis.
    !!
    !! NOTE: rcxq is i/o variable for accumulation. We use E_mu basis when chipm=F.
    !!
    !! ppovl= <I|J> = O , V_IJ=<I|v|J>
    !! (V_IJ - vcoud_mu O_IJ) Zcousq(J, mu)=0, where Z is normalized with O_IJ.
    !! <I|v|J>= \sum_mu ppovl*zcousq(:,mu) v^mu (Zcousq^*(:,mu) ppovl)
    !!
    !! zmelt contains O^-1=<I|J>^-1 factor. Thus zmelt(phi phi J)= <phi |phi I> O^-1_IJ
    !! ppovlz(I, mu) = \sum_J O_IJ Zcousq(J, mu)
    !!
    !! OUTPUT:
    !!  rcxq (npr,npr,nwhis,npm): for given q,
    !!       rcxq(I,J,iw,ipm) =
    !!       Im (chi0(omega))= \sum_k <I_q psi_k|psi_(q+k)> <psi_(q+k)|psi_k> \delta(\omega- (e_i-ej))
    !!      When npm=2 we calculate negative energy part. (time-reversal asymmetry)
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


