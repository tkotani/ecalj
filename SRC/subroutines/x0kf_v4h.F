      module m_x0kf
      use m_keyvalue,only : Getkeyvalue
      use m_pkm4crpa,only : Readpkm4crpa
      use m_zmel,only     : Get_zmel_init,Dconjg_zmel,Deallocate_zmel, zmel ,get_zmel_modex0, GramSchmidt_zmel,rwzmel
      use m_freq,only:  npm, nwhis
      use m_genallcf_v3,only:  nsp=>nspin ,nlmto,nctot,ngrp
      use m_read_bzdata,only:  nqbz,qbas,ginv,nqibz,  rk=>qbz,wk=>wbz
      use m_rdpp,only: nbloch
      use m_readqg,only: ngpmx,ngcmx
      use m_qbze,only: nqbze 
      use m_readhbe,only: nband
      use m_readgwinput,only: nbcut,nbcut2
!! q-dependent quantities
      use m_readVcoud,only: zcousq,ngc,ngb
      use m_tetwt,only:     whw,ihw,nhw,jhw,n1b,n2b,nbnb,nbnbx,nhwtot
      implicit none
      
      public:: X0kf_v4hz_init, X0kf_v4hz, X0kf_v4hz_symmetrize,X0kf_v4hz_init_write,X0kf_v4hz_init_read,x0kf_zmel
      private
      integer:: ncc, ncount
      integer,allocatable:: nkmin(:), nkmax(:),nkqmin(:),nkqmax(:)
      real(8),allocatable:: whwc(:)
      integer,allocatable:: kc(:),iwc(:),itc(:),itpc(:),jpmc(:)
      contains
!--------------------------------------------------------------------------      
      subroutine X0kf_v4hz_init_write(iq,is) 
      integer:: iq,is,ix0
      character*8 :: charext
      open(newunit=ix0,file='x0icount.'//trim(charext(iq))//'_'//trim(charext(is)),form='unformatted')
      write(ix0)ncount
      write(ix0) whwc, kc, iwc, itc,itpc, jpmc
      write(ix0) nkmin,nkmax,nkqmin,nkqmax
      deallocate(nkmin,nkmax,nkqmin,nkqmax)!,skipk)
      deallocate( whwc, kc, iwc, itc, itpc, jpmc )
      close(ix0)
      end subroutine

!--------------------------------------------------------------------------      
      subroutine X0kf_v4hz_init_read(iq,is)
      integer:: iq,is,ix0
      character*8 :: charext
      open(newunit=ix0,file='x0icount.'//trim(charext(iq))//'_'//trim(charext(is)),form='unformatted')
      read(ix0)ncount
      if(allocated(whwc)) deallocate( whwc, kc, iwc, itc, itpc, jpmc, nkmin,nkmax,nkqmin,nkqmax)
      allocate( whwc(ncount), kc(ncount), iwc(ncount), itc(ncount), itpc(ncount), jpmc(ncount) )
      allocate( nkmin(nqbz),nkmax(nqbz),nkqmin(nqbz),nkqmax(nqbz))
      read(ix0) whwc, kc, iwc, itc,itpc, jpmc
      read(ix0) nkmin,nkmax,nkqmin,nkqmax
      close(ix0)
      end subroutine

!--------------------------------------------------------------------------      
      function X0kf_v4hz_init(job,q,isp_k,isp_kq, iq, nmbas,eibzmode, nwgt, crpa) result(ierr)
      intent(in)::            job,q,isp_k,isp_kq, iq, nmbas,eibzmode, nwgt, crpa
!! === initialzation for calling x0kf_v4h.
!! Get ncount index for x0kf_v4h in the private variables
!! Call job=0 and job=1 successively.      
!!   ncount, ngb and nqibz are the key to estimate computational efforts.
!! See x0kf_v4hz.
!     ! We have OUTPUT to private variables nkmin ... jpmc.
      integer::  irot=1,ierr         !, ntqxx,nbmax!,nctot
      integer:: k,isp_k,isp_kq, iq, jpm, ibib, iw,igb2,igb1,it,itp,job,icount, nmbas
      integer::  nwgt(nqbz) 
      real(8):: q(3), imagweight, wpw_k, wpw_kq
      logical :: iww2=.true., eibzmode, crpa, exchange=.false.
!!
      write(6,'(" x0kf_v4hz_init: job q=",i3,3f8.4)') job,q
      if(npm==1) then
         ncc=0
      else
         ncc=nctot !stored to private variable
      endif
      if(job==0) then
         allocate( nkmin(nqbz),nkmax(nqbz),nkqmin(nqbz),nkqmax(nqbz)) !,skipk(nqbz) )
      elseif(job==1) then
         allocate( whwc(ncount), kc(ncount), iwc(ncount), itc(ncount), itpc(ncount), jpmc(ncount) )
      endif
      icount=0
      do 110 k = 1,nqbz
        if( (eibzmode.and.nwgt(k)==0) .or. sum(nbnb(k,1:npm))==0   ) then
           cycle
        endif
        if(job==0) then
           nkmin(k) = 999999
           nkmax(k)= -999999
           nkqmin(k)= 999999
           nkqmax(k)=-999999
           do jpm=1,npm         !npm
              do ibib = 1, nbnb(k,jpm)
                 nkmin(k)  = min(n1b(ibib,k,jpm),nkmin(k))
                 nkqmin(k) = min(n2b(ibib,k,jpm),nkqmin(k))
                 if(n1b(ibib,k,jpm)<=nband) nkmax(k)  = max(n1b(ibib,k,jpm),nkmax(k))
                 if(n2b(ibib,k,jpm)<=nband) nkqmax(k) = max(n2b(ibib,k,jpm),nkqmax(k))
              enddo
           enddo
        endif
c     nkqmin(k)  = nkqmin(k)          ! nkqmin = the num of min   n2 =unocc for jpm=1
c     nt0   = nkmax(k)
c     ntp0  = nkqmax(k) - nkqmin(k) +1
!! sanity check
        if( npm==2.and. nkqmin(k)/=1) then
          write(6,*)' npm==2 nkqmin nkqmax  nkmin nkmax=',nkqmin,nkqmax,nkmin,nkmax
          call rx( " When npm==2, nkqmin==1 should be.")
        endif
c        if(nkmin(k)/=1) call rx( " nkmin==1 should be.") !nov 2021 commneted out for --intraband only
!      ib_k =[1:nctot]              core
!      ib_k =[nctot+nkmin:nctot+nkmax]  valence
!      ib_kq =[1:ncc]             core
!      ib_kq =[ncc+nkqmin:ncc+nkqmax]  valence range [nkqmin,nkqmax]
!   If jpm=1, ncc=0.
!   If jpm=2, ncc=ncore. nkqmin(k)=1 should be.
! There is a little confusion. n1b index contains cores are after valence.
! You can see codes to treat the confusion.
! NOTE: 
!  q+rk n2b vec_kq  vec_kq_g geig_kq cphi_kq  ngp_kq ngvecp_kq  isp_kq
!    rk n1b vec_k   vec_k_g  geig_k  cphi_k   ngp_k  ngvecp_k   isp_k
        do 1251 jpm  = 1, npm ! npm=1 usually (or =2)
        do 125 ibib = 1, nbnb(k,jpm) !---  ibib loop, ibib is decomposed to band index pair, it and itp
!! n1b,n2b --> core after valence.  it,itp --> valence after core
c          it_(ibib, jpm,k)= -999
c          itp_(ibib,jpm,k)= -999
          if(n1b(ibib,k,jpm) <= nband) then
            it = nctot + n1b(ibib,k,jpm) !valence
            if(it > nctot + nkmax(k) ) cycle
          else
            it = n1b(ibib,k,jpm) - nband !core
          endif
          if( n2b(ibib,k,jpm) <= nband) then
            itp = ncc + n2b(ibib,k,jpm) - nkqmin(k) + 1 !val
            if(itp > ncc + nkqmax(k)-nkqmin(k)+1 ) cycle
          else
            itp =  n2b(ibib,k,jpm) - nkqmin(k) + 1 - nband !core
         endif
!! nbcut mechanism
          if(nbcut==0.and.nbcut2==999999) then
             continue
          else   
             if(jpm==1) then 
                if( n1b(ibib,k,jpm) <= nbcut .and. nbcut2<n2b(ibib,k,jpm) ) then
                   if(iww2) then
                      write(6,"(' nband_chi0 nbcut nbcut2 n2b n1b=',4i6)") nbcut,n2b(ibib,k,jpm),n1b(ibib,k,jpm)
                      iww2=.false.
                   endif
                   cycle
                endif
             else               !jpm==2 
                if( n2b(ibib,k,jpm) <= nbcut .and. nbcut2<n1b(ibib,k,jpm) ) then
                   if(iww2) then
                      write(6,"(' nband_chi0 nbcut nbcut2 n2b n1b=',4i6)") nbcut,n2b(ibib,k,jpm),n1b(ibib,k,jpm)
                      iww2=.false.
                   endif
                   cycle
                endif
             endif
          endif
!! sanity check
          if (ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1 >nwhis) call rx( "x0kf_v4hz: iw>nwhis")
!! constraint RPA mode (juelich verison)
          if(crpa) then
            if(n1b(ibib,k,jpm) <= nband) then
              wpw_k= readpkm4crpa(n1b(ibib,k,jpm),   rk(:,k), isp_k)
            else
              wpw_k=0d0
            endif  
            if(n2b(ibib,k,jpm) <= nband) then
              wpw_kq= readpkm4crpa(n2b(ibib,k,jpm), q+rk(:,k), isp_kq)
            else
              wpw_kq= 0d0
            endif  
          endif
          do iw = ihw(ibib,k,jpm),ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1 !iiww=iw+ihw(ibib,k)-1
             imagweight = whw(jhw(ibib,k,jpm)+iw-ihw(ibib,k,jpm))
             if(crpa)     imagweight = imagweight*(1d0-wpw_k*wpw_kq)
             if(eibzmode) imagweight = nwgt(k)*imagweight
             icount = icount+1
             if(job==1) then
                whwc (icount)= imagweight
                kc   (icount) = k
                iwc  (icount)= iw
                itc  (icount)= it
                itpc (icount)= itp
                jpmc (icount)= jpm
c        write(6,"(a,6i5,d13.5)")'uuuuu k it itp iw jpm whw=',icount,k,it,itp,iw,jpm,whwc(icount)
             endif
          enddo                 ! iw
 125     continue
 1251    continue
 110  continue
      ncount = icount
      ierr=0
      if(job==0) write(6,"('x0kf_v4hz_init: job=0 ncount ngb nqibz=',3i8)") ncount, ngb, nqibz
      end function x0kf_v4hz_init
!! --------------------------------------------------------------------------------
      subroutine X0kf_zmel ( q, isp_k,isp_kq, iq, nmbas, eibzmode, nwgt, rcxq,epsppmode)
      intent(in)   ::        q, isp_k,isp_kq, iq, nmbas, eibzmode, nwgt,      epsppmode 
!! === calculate chi0, or chi0_pm ===
!! We calculate imaginary part of chi0 along real axis.
!!
!! NOTE: rcxq is i/o variable for accumulation. We use E_mu basis when chipm=F.
!!  
!!
!! ppovl= <I|J> = O , V_IJ=<I|v|J>
!! (V_IJ - vcoud_mu O_IJ) Zcousq(J, mu)=0, where Z is normalized with O_IJ.
!! <I|v|J>= \sum_mu ppovl*zcousq(:,mu) v^mu (Zcousq^*(:,mu) ppovl)     
!!
!! zmelt contains O^-1=<I|J>^-1 factor. Thus zmelt(phi phi J)= <phi |phi I> O^-1_IJ
!! ppovlz(I, mu) = \sum_J O_IJ Zcousq(J, mu)
!!
!! OUTPUT:
!!  rcxq (nmbas,nmbas,nwhis,npm): for given q,
!!       rcxq(I,J,iw,ipm) = 
!!       Im (chi0(omega))= \sum_k <I_q psi_k|psi_(q+k)> <psi_(q+k)|psi_k> \delta(\omega- (e_i-ej))
!!      When npm=2 we calculate negative energy part. (time-reversal asymmetry)
!!
!!note: zmel: matrix element <phi phi |M>
!! nlmto   = total number of atomic basis functions within MT
!! nqbz    = number of k-points in the 1st BZ
!!
      logical,optional ::epsppmode
      integer:: k,isp_k,isp_kq,iq, jpm, ibib, iw,igb2,igb1,it,itp
      real(8):: q(3)
      complex(8):: rcxq (nmbas,nmbas,nwhis,npm)
      complex(8) :: imag=(0d0,1d0),trc,aaa 
      integer::   nadd(3), igc !nband,!ngpmx, ngcmx,nqbze, ngc,
      complex(8),allocatable :: zmelt(:,:,:)
      complex(8),allocatable::  z1p(:,:)
      logical,parameter:: debug=.false.
      real(8) :: imagweight 
      integer:: nmbas !, imb1,imb2, imb verbose,
      logical :: iww2=.true.
      complex(8):: img=(0d0,1d0),zmelt2 !,zzz(ngbb)
      integer ::  nkmax1,nkqmax1, ib1, ib2, ngcx,ix,iy,igb !nkqmin, nkqmax,
      logical :: eibzmode  
      integer::  nwgt(nqbz) 
      logical::  exchange=.false. !eibz4x0,
      real(8)::  wpw_k,wpw_kq
      integer::  irot=1, neibz,icc,ig,eibzmoden,ikp,i,j,itimer,icount,iele
      integer:: ieqbz,kold,nxxxx 
      integer::nkmin_,nkqmin_,nkoff,nkqoff,ispold,izmel,nmini,nqini,nmtot,nqtot!,ispold2
      kold = -999
      do 1000 icount = 1,ncount
        k = kc(icount)
        nkmin_  = nkmin(k)
        nkqmin_ = nkqmin(k)
        nmini= nkmin_ 
        nqini= nkqmin_
        if(kold/=k .or. ispold/=isp_k ) then !.or. ispold2/=isp_kq) then  
           call Get_zmel_modex0(nkmin_,nkqmin_,isp_k,isp_kq)!oct-2021(nkmin(k),nkqmin(k),isp_k,isp_kq)
           call Get_zmel_init(exchange, q+rk(:,k), q, irot, q, nxxxx, nkmax(k),nkqmax(k),nctot,ncc)
           call Dconjg_zmel()   !zmel = dconjg(zmel)
           write(6,*) 'k zmel writing...',iq,k,isp_k
           call rwzmel(iq,k,isp_k,'w')
           kold = k
           ispold=isp_k
        endif
 1000 continue
      call Deallocate_zmel()
      end subroutine x0kf_zmel
      
!! --------------------------------------------------------------------------------
      subroutine X0kf_v4hz ( q, isp_k,isp_kq, iq, nmbas, eibzmode, nwgt, rcxq,epsppmode)
      intent(in)   ::        q, isp_k,isp_kq, iq, nmbas, eibzmode, nwgt,      epsppmode 
!! === calculate chi0, or chi0_pm ===
!! We calculate imaginary part of chi0 along real axis.
!!
!! NOTE: rcxq is i/o variable for accumulation. We use E_mu basis when chipm=F.
!!  
!!
!! ppovl= <I|J> = O , V_IJ=<I|v|J>
!! (V_IJ - vcoud_mu O_IJ) Zcousq(J, mu)=0, where Z is normalized with O_IJ.
!! <I|v|J>= \sum_mu ppovl*zcousq(:,mu) v^mu (Zcousq^*(:,mu) ppovl)     
!!
!! zmelt contains O^-1=<I|J>^-1 factor. Thus zmelt(phi phi J)= <phi |phi I> O^-1_IJ
!! ppovlz(I, mu) = \sum_J O_IJ Zcousq(J, mu)
!!
!! OUTPUT:
!!  rcxq (nmbas,nmbas,nwhis,npm): for given q,
!!       rcxq(I,J,iw,ipm) = 
!!       Im (chi0(omega))= \sum_k <I_q psi_k|psi_(q+k)> <psi_(q+k)|psi_k> \delta(\omega- (e_i-ej))
!!      When npm=2 we calculate negative energy part. (time-reversal asymmetry)
!!
!!note: zmel: matrix element <phi phi |M>
!! nlmto   = total number of atomic basis functions within MT
!! nqbz    = number of k-points in the 1st BZ
!!
      logical,optional ::epsppmode
      integer:: k,isp_k,isp_kq,iq, jpm, ibib, iw,igb2,igb1,it,itp
      real(8):: q(3)
      complex(8):: rcxq (nmbas,nmbas,nwhis,npm)
      complex(8) :: imag=(0d0,1d0),trc,aaa 
      integer::   nadd(3), igc !nband,!ngpmx, ngcmx,nqbze, ngc,
      complex(8),allocatable :: zmelt(:,:,:)
      complex(8),allocatable::  z1p(:,:)
      logical,parameter:: debug=.false.
      real(8) :: imagweight 
      integer:: nmbas !, imb1,imb2, imb verbose,
      logical :: iww2=.true.
      complex(8):: img=(0d0,1d0),zmelt2 !,zzz(ngbb)
      integer ::  nkmax1,nkqmax1, ib1, ib2, ngcx,ix,iy,igb !nkqmin, nkqmax,
      logical :: eibzmode  
      integer::  nwgt(nqbz) 
      logical::  exchange=.false. !eibz4x0,
      real(8)::  wpw_k,wpw_kq!, vec_kcrpa(3),vec_kqcrpa(3)
      integer::  irot=1         !, ntqxx,nbmax!,nctot
      integer::  neibz,icc,ig,eibzmoden,ikp,i,j,itimer,icount,iele!,ngbb !eibzsym(ngrp,-1:1),
      integer:: ieqbz,kold,nxxxx !igx(ngrp*2,nqbz),igxt(ngrp*2,nqbz),
      integer::nkmin_,nkqmin_,izmel,ispold,nmini,nqini,nmtot,nqtot !,ispold2
c      logical:: interbandonly,intrabandonly,cmdopt0
      
c      integer,allocatable:: it_(:,:,:),itp_(:,:,:)
!!-----------------------------------------------------------
!! Main loop over k-points ---------------------------------------------------------
!! z1p = <M_ibg1 psi_it | psi_itp> < psi_itp | psi_it M_ibg2 >
!     ! zxq(iw,ibg1,igb2) = \sum_ibib imgw(iw,ibib)* z1p(ibib, igb1,igb2) !ibib means band pair (occ,unocc)
!      
!  zzmel(1:nbloch, ib_k,ib_kq)
!      ib_k =[1:nctot]              core
!      ib_k =[nctot+nkmin:nctot+nkmax]  valence
!      ib_kq =[1:ncc]             core
!      ib_kq =[ncc+nkqmin:ncc+nkqmax]  valence range [nkqmin,nkqmax]
!   If jpm=1, ncc=0.
!   If jpm=2, ncc=ncore. nkqmin(k)=1 should be.
! There is a little confusion. n1b index contains cores are after valence.
! You can see codes to treat the confusion.
! NOTE: 
!  q+rk n2b vec_kq  vec_kq_g geig_kq cphi_kq  ngp_kq ngvecp_kq  isp_kq
!    rk n1b vec_k   vec_k_g  geig_k  cphi_k   ngp_k  ngvecp_k   isp_k
c
!! -------------------------------------------------------------------------------
!!note: for usual correlation mode, I think nctot=0
!!--- For dielectric funciton, we use irot=1 kvec=rkvec=q
!            < MPB      middle   |   end >
!!              q      rkvec     | q + rkvec  
!                      nkmin:nt0 | nkqmin:ntp0
!                         occ    | unocc     
!                      (nkmin=1)
!                      (cphi_k  | cphi_kq !in x0kf)
!
!!     rkvec= rk(:,k)   ! <phi(q+rk,nqmax)|phi(rk,nctot+nmmax)  MPB(q,ngb )>
!!     qbz_kr= rk(:,k)  ! 
!!     qibz_k= rk(:,k)  ! k
!         
!! Get_zmelt in m_zmel gives the matrix element zmel,  ZO^-1 <MPB psi|psi> , where ZO is ppovlz
!! zmel(ngb, nctot+nt0,ncc+ntp0) in m_zmel
!        nkmin:nt0,           nkqmin:ntp0
!     nt0=nkmax-nkmin+1  , ntp0=nkqmax-nkqmin+1
      kold = -999
      do 1000 icount = 1,ncount
        if(debug) call cputid(0)
        k = kc(icount)
! oct15 2021
        nkmin_  = nkmin(k)
        nkqmin_ = nkqmin(k)
        nmini= nkmin_ 
        nqini= nkqmin_
! This gives matrix element zmel(igb q,  k it occ,   q+k itp unocc)
        if(kold/=k.or.isp_k/=ispold ) then !.or. ispold2/=isp_kq) then  
c         call Get_zmel_modex0(nkmin_,nkqmin_,isp_k,isp_kq)!oct-2021(nkmin(k),nkqmin(k),isp_k,isp_kq)
c         call Get_zmel_init(exchange, q+rk(:,k), q, irot, q, nxxxx, nkmax(k),nkqmax(k),nctot,ncc)
c         call Dconjg_zmel() !zmel = dconjg(zmel)
cc           nmtot = nctot + nkmax(k)  -nmini+1 ! = phi_middle
cc           nqtot = ncc   + nkqmax(k) -nqini+1 ! = phi_end
           call rwzmel(iq,k,isp_k,'r') !,nmbas,nmtot,nqtot)
           kold=k
           ispold=isp_k
        endif
c----------------------------
!!  z1p = <M_ibg1 psi_it | psi_itp> < psi_itp | psi_it M_ibg2 >
!!  zxq(iw,ibg1,igb2) = sum_ibib wwk(iw,ibib)* z1p(ibib, igb1,igb2)
!! n1b,n2b --> core after valence.  it,itp --> valence after core
        it  = itc(icount)  !occ      k
        itp = itpc(icount) !unocc    q+k
        iw  = iwc(icount)  !omega-bin
        jpm = jpmc(icount)      ! \pm omega
        do igb2=1,nmbas  !this part dominates cpu time most time consuming...........
        do igb1=1,igb2
             rcxq(igb1,igb2,iw,jpm) = 
     &       rcxq(igb1,igb2,iw,jpm) + dconjg(zmel(igb1,it,itp))*zmel(igb2,it,itp) * whwc(icount)
!! NOTE                               whwc is the img weight given by tetrahedron method.
        enddo     
        enddo     
 1000 continue
      call Deallocate_zmel()
      deallocate(nkmin,nkmax,nkqmin,nkqmax)!,skipk)
      deallocate(whwc,kc,itc,itpc,iwc,jpmc) !z1p
!! Hermitianize. jun2012takao moved from dpsion5 ====
      do jpm=1,npm
      do iw= 1,nwhis
         do igb2= 1,nmbas       !eibzmode assumes nmbas1=nmbas2
         do igb1= 1,igb2-1
               rcxq(igb2,igb1,iw,jpm) = dconjg(rcxq(igb1,igb2,iw,jpm))
         enddo
         enddo
      enddo
      enddo
 9999 continue
      write(6,"(' --- x0kf_v4hz: end')") !, 3d13.5)")
      if(debug) write(6,"(' --- ', 3d13.5)")
     &     sum(abs(rcxq(1:nmbas,1:nmbas,1:nwhis,1:npm))),sum((rcxq(1:nmbas,1:nmbas,1:nwhis,1:npm)))
      end subroutine x0kf_v4hz
      
!! --------------------------------------------------------------------------------
      subroutine X0kf_v4hz_symmetrize(q, iq, nolfco, zzr,nmbas, chipmzzr, eibzmode, eibzsym, rcxq)
      use m_readqgcou,only: qtt_, nqnum
      use m_rotMPB,only:  rotMPB2 
      intent(in)::                    q, iq, nolfco, zzr,nmbas, chipmzzr, eibzmode, eibzsym
!! === symmetrization for EPIBZ mode ===
      integer(4):: nnc,iq,iatom,nctot,nbmx,iopen 
     &  ,jpm,ibib,itps,nt0,ntp0,ngp_kq,ngp_k,it,itp,iw,igb2,igb1, nn,no,isx,iclose,k 
      real(8):: q(3),ebmx 
      complex(8):: rcxq (nmbas,nmbas,nwhis,npm)
      complex(8):: imag=(0d0,1d0),trc,aaa !phase(natom),
      integer(4):: nadd(3)
      logical,parameter:: debug=.false.
      complex(8) :: zmelt1,zmelt2,zmeltt(ngb)      !...........................sf 21May02
      real(8) :: imagweight !............................sf 21May02
      real(8):: eband(nband)!,ebandr(nband),ebandqr(nband)
      integer(4):: verbose
      logical   :: nolfco !iepsmode
      integer(4):: nmbas, imb1,imb2, imb !nmbas1x !nmbas2,nmbas1,
      real(8):: vec_kq_g(3),vec_k_g(3),vec_kq(3),vec_k(3),quu(3),tolq=1d-8,quu1(3),quu2(3)!tolqu=1d-4,
      integer(4):: nbcut,nbcut2
      logical :: iww1=.true.,iww2=.true.
      complex(8):: img=(0d0,1d0)
      integer(4):: nkmin,  nkmax, nkqmin, nkqmax,nkmax1,nkqmax1
      integer(4):: ib1, ib2,      ngcx,ix,iy
      complex(8),target :: zzr(ngb,nmbas) !ppovlz(ngb,ngb),
      integer:: igb
      logical:: checkbelong,eibzmode, chipmzzr
      complex(8)::  zcousqc(ngb,ngb) !zcousq(ngb,ngb) , 
      integer::  eibzsym(ngrp,-1:1),neibz,icc,ig,eibzmoden,ikp,i,j,itimer,icount,iele
      integer:: irotm,nrotmx,ixx,iyy,itt,ntimer, nccc, nxx,iagain,irotm1,irotm2
      integer,allocatable:: i1(:,:),i2(:,:),nrotm(:)
      complex(8),allocatable:: zrotm(:,:),zrr(:,:),zrrc(:,:),zrr_(:,:,:),zrrc_(:,:,:),zmmm(:), 
     & zrrx(:,:),rcxq_core(:,:), zcousqr(:,:,:),rcxq0(:,:),rcxq00(:,:),rcxq000(:,:),rcxqwww(:,:)
      logical :: exchange=.false.
      integer:: irot=1
      integer:: ntqxx,nbmax
!! == Symmetrizer of EIBZ PRB.81,125102(2010) Eq.(51) july2012takao ==
!! This may be not so effective ---> only for limited cases?
!! --- zrotm(J,J') = <Mbar^k_J| \hat{A}^k_i Mbar^k_J'>. ---
!! We do \sum_i T_alpha_i [ zrotm_i^dagger (I,I') P_I'J' zrom_i(J'J) ] 
!! (exactrly speaking, we insert conversion matrix between Enu basis and M_I basis).
!!
!! input qin = q
!! \hat{A}^k_i  is specified by symops(:,:,igx),and igxt (-1 for time-reversal).
!! Note that k= \hat{A}^k_i(k) (S_A^k)
!! See Eq.(51) around in PRB81 125102(2010)
!! 
!! === zmelt conversion === 
      if(nolfco .and. nmbas==1) then
         write(6,*)' nmbas=1 nolfco=T ---> not need to symmetrize'
         goto 9999
      endif
!!
      if(eibzmode) then
        call iqindx2(q, ginv, qtt_, nqnum, ikp,quu) !to get ikp for timereversal mode
        if(sum(abs(q-quu))>tolq) call rx( 'x0kf_v4h_symmetrize: eibz 111 q/quu')
        neibz = sum(eibzsym(:,1))+sum(eibzsym(:,-1)) 
             !itimer=-1 means time reversal. eibzsym(ig,itimer) where ig: space rotation.
        write(6,"(' --- goto symmetrization --- ikp neibz q=',2i3,3f12.8)")ikp,neibz,q
        call cputid2(' --- x0kf: start symmetrization  ',0)
        ntimer=1
        if(sum(eibzsym(:,-1))>0) ntimer=2 !timereversal case
        allocate(zrotm(ngb,ngb),nrotm(ngrp*2))
!!
!! == Assemble rotantion matrx zrr,zrrc ==
!! Rotation matrix zrrx can be a sparse matrix. 
!! Thus it is stored to  "i1(nrotmx,nccc),i2(nrotmx,nccc),zrr(nrotmx,icc),nrotm(icc)".
!! See folloings: matmul(rcxqwww,zrrx) is given by
!!     do irotm1 = 1,nrotm(icc)
!!       rcxq0(:,i2(irotm1,icc)) = rcxqwww(:,i1(irotm1,icc)) * zrr(irotm1,i2(irotm1,icc))

        allocate(zrrx(nmbas,nmbas))
        nrotmx = 10000 !trial value
        do 1011  !this loop is only in order to to set large enough nrotmx.
        if(allocated(i1)) deallocate(i1,i2,zrr,zrrc)
        nccc=ngrp*2
        allocate(i1(nrotmx,nccc),i2(nrotmx,nccc),zrr(nrotmx,nccc),zrrc(nrotmx,nccc)) 
        i1=-99999
        i2=-99999
        zrr=-99999d0
        zrrc=-99999d0
        call cputid2(' --- x0kf:11111   :',0)
        icc=0
        do itimer=1,-1,-2
          if(ntimer==1.and.itimer==-1) exit
          if(itimer==1 ) itt=1
          if(itimer==-1) itt=2
          do ig=1,ngrp
            if(eibzsym(ig,itimer)==1) then
              icc=icc+1
!! Get rotation matrix zrrx, which can be a sparse matrix. Thus stored to zrr.
              call rotMPB2(nbloch,ngb,q,ig,itimer,ginv,zrotm) 
              if(nolfco.and.chipmzzr) then 
!!   We assume <svec_I | svec_J >= \delta_IJ, In addition, we use fact that we have no IPW parts in svec.
!!   If IPW part exist, we may have to take into account <IPW|IPW> matrix, e.g. as in ppovlz. 
!!   svec --> zzr
               if(itimer==1) then
                  zrrx= matmul(transpose(dconjg(zzr)), matmul(zrotm, zzr)) 
                else  
                  zrrx= matmul(transpose(zzr), matmul(dconjg(zrotm), zzr)) 
                endif 
              elseif(nolfco) then 
                call rx( 'x0kf_v4h_symmetrize: this case is not implemented xxxxxxxxxxxxxx')
              else   
!! zrotm(J,J') is the rotation matrix = <Mbar^k_J| \hat{A}^k_i Mbar^k_J'>
!! See rotMPB2 defined in readeigen.F.
!! zrrx(mu nu)= dconjg(Zcousq(I, mu)) *zrotm(I,J)* Zcousq(J, nu)
!! zrrx is very sparse matrix. Size is \sim ngb or something.

c$$$                if(itimer==1) then
c$$$                  call matmmsparse(zcousqinv,zrotm,zcousq,zrrx,ngb,1d-8,iele)
c$$$                  ! this means zrrx= matmul(zcousqinv,matmul(zrotm, zcousq)) 
c$$$                else  
c$$$                  call matmmsparse(dconjg(zcousqinv),dconjg(zrotm),zcousq,zrrx,ngb,1d-8,iele) 
c$$$                  ! this means zrrx= matmul(dconjg(zcousqinv),matmul(dconjg(zrotm), zcousq)) 
c$$$                endif 
                if(itimer==1) then
                  zrrx=zrotm
c                  call matmmsparse(zcousqinv,zrotm,zcousq,zrrx,ngb,1d-8,iele)
                  ! this means zrrx= matmul(zcousqinv,matmul(zrotm, zcousq)) 
                else  
                  zrrx=dconjg(zrotm)
c                  call matmmsparse(dconjg(zcousqinv),dconjg(zrotm),zcousq,zrrx,ngb,1d-8,iele) 
                  ! this means zrrx= matmul(dconjg(zcousqinv),matmul(dconjg(zrotm), zcousq)) 
                endif 
              endif
              i1(:,icc)=0
              i2(:,icc)=0
              irotm=0
              iagain=0
              do ix=1,ngb
              do iy=1,ngb
                if(abs(zrrx(ix,iy))>1d-8) then
                  irotm=irotm+1
                  if(irotm>nrotmx) then
                    iagain=1
                  endif
                  if(iagain/=1) then
                  i1(irotm,icc)=ix
                  i2(irotm,icc)=iy
                  zrr(irotm,icc) = zrrx(ix,iy)
                  zrrc(irotm,icc)= dconjg(zrr(irotm,icc))
                  endif
                endif   
              enddo
              enddo
              if(iagain==1) then
                nrotmx=irotm !enlarge allocation and do things again.
                write(6,*)' warn:(slow speed) xxxx goto 1011 xxxxxx nrotmx+=nrotmx+10000 again'
                goto 1011
                !enlarge nrotmx ang try it again.
              endif   
              nrotm(icc)=irotm
              if(debug) write(6,*)'ig itimer icc nrotm=',ig,itimer,icc,nrotm(icc) ,iele
            endif
          enddo
        enddo
        exit
 1011   continue !only when nrotmx overflow.

!! === main part to obtain symmetrized rcxq  ===
!! neibz is total number of symmetrization operation. 
!!      rcxq is rotated and accumulated; finally divied by neibz
        zcousqc = dconjg(transpose(zcousq))
        if(debug) call cputid2(' --- x0kf:qqqqq222ini:',0)
!$OMP parallel private(rcxq000,icc,itt,icount,rcxqwww,rcxq00,rcxq0,rcxq_core)
        allocate(rcxq0(ngb,ngb),rcxq00(ngb,ngb),rcxq000(ngb,ngb),rcxqwww(ngb,ngb),rcxq_core(ngb,ngb))
!$OMP master
!$         write(6,'(a,i5,a,i5)') 'OMP parallel nwhis, threads=',omp_get_num_threads(),' nwhis=',nwhis
!$OMP end master
!$OMP do 
        do iw=1,nwhis
        do jpm=1,npm
          rcxq000 = 0d0
          icc=0
          do itimer=1,-1,-2
            if(itimer==1 ) itt=1
            if(itimer==-1) itt=2
            icount=0
            if(itimer==1) then
              rcxqwww = rcxq(:,:,iw,jpm)
            else  
              rcxqwww = transpose(rcxq(:,:,iw,jpm))
            endif  
            rcxq00 = 0d0
            do ig=1,ngrp
              if(eibzsym(ig,itimer)==1) then
               icount=icount+1
               icc=icc+1 
               rcxq0 =0d0

c$$$               if(itimer==1) then
c$$$               do irotm1 = 1,nrotm(icc)
c$$$               do irotm2 = 1,nrotm(icc)
c$$$               rcxq0(i2(irotm2,icc),i2(irotm1,icc)) =rcxq0(i2(irotm2,icc),i2(irotm1,icc)) 
c$$$     &              +    zrrc(irotm2,icc)* rcxq(i1(irotm2,icc),i1(irotm1,icc),iw,jpm)*zrr(irotm1,icc)
c$$$               enddo
c$$$               enddo
c$$$               else
c$$$               do irotm1 = 1,nrotm(icc)
c$$$               do irotm2 = 1,nrotm(icc)
c$$$               rcxq0(i2(irotm1,icc),i2(irotm2,icc)) =rcxq0(i2(irotm1,icc),i2(irotm2,icc)) !transpose
c$$$     &              +    zrrc(irotm2,icc)* rcxq(i2(irotm2,icc),i1(irotm1,icc),iw,jpm)*zrr(irotm1,icc)
c$$$               enddo
c$$$               enddo
c$$$               endif

!!  Followings are equivalent with
!!            rcxq00= rcxq00 + matmul(zrrc_(:,:,icc),matmul(rcxqwww,zrr_(:,:,icc)))
               do irotm1 = 1,nrotm(icc)
c                 if(abs(zrr(irotm1,icc))<1d-8) cycle
                 rcxq0(:,i2(irotm1,icc)) =rcxq0(:,i2(irotm1,icc)) + rcxqwww(:,i1(irotm1,icc)) * zrr(irotm1,icc)
               enddo
               do irotm2 = 1,nrotm(icc)
c                 if(abs(zrrc(irotm2,icc))<1d-8) cycle
                 rcxq00(i2(irotm2,icc),:)= rcxq00(i2(irotm2,icc),:) + zrrc(irotm2,icc) * rcxq0(i1(irotm2,icc),:)
               enddo

c               if(itimer==1) then
c                 rcxq000 = rcxq000 + rcxq00
c               else  
c                 rcxq000 = rcxq000 + transpose(rcxq00)
c               endif
c
c$$$               do irotm = 1,nrotm(icc)
c$$$                iyy = i1(irotm,icc)
c$$$                iy  = i2(irotm,icc)
c$$$                rcxq0(:,iy)= rcxq0(:,iy)+ rcxq(:,iyy,iw,jpm)* zrr(irotm,icc) 
c$$$               enddo  
c$$$               do irotm = 1,nrotm(icc)
c$$$                iyy = i1(irotm,icc)
c$$$                iy  = i2(irotm,icc)
c$$$                rcxq00(iy,:)= rcxq00(iy,:)+ dconjg(zrr(irotm,icc)) * rcxq0(iyy,:)
c$$$               enddo
c$$$
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c               if(iw==1.and.jpm==1) then
c                  write(6,"('bbbbbbb ig icc iw jpm rcxq', 4i3, 13d13.6)")
c     &                 ig,icc,iw,jpm, sum(abs(rcxq00)), rcxq00(1,1),sum(abs(rcxqwww)),sum((rcxqwww))
c               endif   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

              endif
            enddo

c$$$           if(itimer==1) then
c$$$             rcxq000(:,:) = matmul(zcousqc,matmul(rcxq00,zcousq))
c$$$c$$$c               call zgemm("N","N",ngb,ngb,ngb, (1d0,0d0), rcxq00, ngb, zcousq,ngb, (0d0,0d0), rzc,ngb)
c$$$c$$$c               call zgemm("N","N",ngb,ngb,ngb, (1d0,0d0), zcousqc,ngb, rzc,ngb, (0d0,0d0), rcxq000,ngb)
c$$$           elseif(icount>0) then  
c$$$c$$$c           write(6,*)'qqqqq icount=',icount
c$$$c$$$c           rcxq000(:,:) = rcxq000(:,:) + transpose(matmul(transpose(zcousq),matmul(rcxq00,dconjg(zcousq))))
c$$$             rcxq000(:,:) = rcxq000(:,:) +   matmul(matmul(zcousqc,transpose(rcxq00)),zcousq)
c$$$           endif

            if(itimer==1) then
               rcxq000=rcxq00
            else
               rcxq000=rcxq000+rcxq00
            endif   
          enddo
          rcxq_core = rcxq000/neibz
#if 1
!! matmul(rcxq(:,:,iw,jpm),zcousq) fails in ifort 14.0.3.
!! It looks that ifort 14.0.3 has a bug
!! But, zgemm works. So I changed like that. 
          call zgemm('N','N',ngb,ngb,ngb,(1.0d0,0.0d0),rcxq_core,ngb,zcousq, ngb, (0.0d0,0.0d0),rcxq000,ngb)
          call zgemm('N','N',ngb,ngb,ngb,(1.0d0,0.0d0),zcousqc  ,ngb,rcxq000,ngb, (0.0d0,0.0d0),rcxq_core,ngb)
          rcxq(:,:,iw,jpm) = rcxq_core 
#else
          rcxq(:,:,iw,jpm) = matmul(zcousqc,matmul(rcxq_core,zcousq))
#endif
          enddo
          enddo
!$OMP end  do 
        deallocate(rcxq00,rcxq000,rcxq0,rcxqwww)
!$OMP end parallel 
        deallocate(zrotm,i1,i2)

c$$$        allocate(zcousqr(ngb,ngb,neibz),rcxq0(ngb,ngb),rcxq00(ngb,ngb),rcxqtr(ngb,ngb))
c$$$        icc=0
c$$$        do itimer=1,-1,-2
c$$$        do ig=1,ngrp
c$$$          if(eibzsym(ig,itimer)==1) then
c$$$            icc=icc+1
c$$$            if(itimer==1) then
c$$$              call rotMPB(zcousq,nbloch,ngb,q,ig,itimer,ginv,zcousqr(1,1,icc))
c$$$            else
c$$$!! time reversal mapping ---
c$$$              call rotMPB(dconjg(zcousq),nbloch,ngb,q,ig,itimer,ginv,zcousqr(1,1,icc))
c$$$            endif
c$$$          endif
c$$$        enddo
c$$$        enddo
c$$$
c$$$        do iw=1,nwhis
c$$$        do jpm=1,npm
c$$$          rcxq0=0d0
c$$$          icc=0
c$$$c          do itimer=1,1 !1,-1,-2
c$$$          do itimer=1,-1,-2
c$$$          do ig=1,ngrp
c$$$            if(eibzsym(ig,itimer)==1) then
c$$$             icc=icc+1
c$$$             rcxq00(:,:) = matmul(dconjg(transpose(zcousqr(:,:,icc))),
c$$$     &                          matmul(rcxq(:,:,iw,jpm),zcousqr(:,:,icc)))
c$$$!! time reversal mapping ---
c$$$             if(itimer==-1) rcxq00(:,:) = transpose(rcxq00)
c$$$             rcxq0(:,:) = rcxq0(:,:)+ rcxq00(:,:)
c$$$            endif
c$$$          enddo
c$$$          enddo
c$$$          rcxq(:,:,iw,jpm)=rcxq0(:,:)/neibz
c$$$        enddo
c$$$        enddo
c$$$        deallocate(zcousqr,rcxq0,rcxq00,rcxqtr)
        if(debug) call cputid2(' --- qqqqq222end:',0)
      endif  
 9999 continue
      write(6,"(' --- x0kf_v4hz_symmetrize: end')")
      end subroutine x0kf_v4hz_symmetrize
      
      end module m_x0kf

