      subroutine x0kf_v4hz (npm,ncc,
     i            ihw,nhw,jhw,whw,nhwtot, ! tetwt5
     i            n1b,n2b,nbnbx,nbnb,     !  use whw by tetwt5 ,
     i                  q,
     i                  nsp,isp_k,isp_kq,!symmetrize,
     i                  qbas,ginv,rk,wk,  
     d                  nlmto,nqbz,nctot, 
     d                  nbloch,nwt, 
     i      iq, ngbb, ngc, ngpmx,ngcmx,  
     i      nqbze, nband,nqibz, 
     o      rcxq,
     i     nolfco,zzr,nmbas, zcousq,
     i     chipmzzr,eibzmode, 
     i     nwgt,igx,igxt,ngrp,eibzsym, crpa)
      use m_readqg,only   :  readqg
      use m_readeigen,only:  readeval 
      use m_keyvalue,only   :  getkeyvalue
      use m_rotMPB,only   :  rotMPB2 
      use m_readqgcou,only:  
     o qtt_, nqnum 
      use m_pkm4crpa,only : Readpkm4crpa
      use m_zmel,only     : Get_zmelt2,Dconjg_zmel,Deallocate_zmel,
     o zmel 

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
!! when nmbas1=2, this works in a special manner for nolfco=T chipm=F. mar2012takao
!!
!! 
!! OUTPUT:
!!  rcxq (nmbas,nmbas,nwt,npm): for given q,
!!       rcxq(I,J,iw,ipm) = 
!!       Im (chi0(omega))= \sum_k <I_q psi_k|psi_(q+k)> <psi_(q+k)|psi_k> \delta(\omega- (e_i-ej))
!!      When npm=2 we calculate negative energy part. (time-reversal asymmetry)
!!
!!
!!!!!
!! zmel: matrix element <phi phi |M>
!!    q,iq:
!! nlmto   = total number of atomic basis functions within MT
!! nqbz    = number of k-points in the 1st BZ
!!
      implicit none
      integer:: npm,ncc,ngbb,natom,nwt,nsp,isp_k,isp_kq,nlmto 
     &  ,nl,nclass,nnc,nbloch,iq,nqibz,iatom,nctot,nbmx,iopen 
     &  ,jpm,ibib,itps,nt0,ntp0,ngp_kq,ngp_k,it,itp,iw,igb2,igb1,ngb
     &  ,nn,no,isx,iclose,k,nbnbx,nqbz
      real(8):: q(3),qbas(3,3),ginv(3,3),rk(3,nqbz),wk(nqbz),ebmx
      complex   (8):: rcxq (nmbas,nmbas,nwt,npm)
      complex(8) :: imag=(0d0,1d0),trc,aaa 
      complex(8),allocatable:: cphi_k(:,:),cphi_kq(:,:),geig_kq(:,:),geig_k(:,:)
      integer:: ngpmx, ngcmx, nqbze, nband,
     &             ngc,nadd(3), igc
      complex(8),allocatable :: zmelt(:,:,:)
      real(8) :: qbasinv(3,3), det,qdiff(3),add(3),symope(3,3)
     &    ,shtv(3)=(/0d0,0d0,0d0/)
      data symope /1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0/
      complex(8),allocatable :: ttx(:,:)
      complex(8),allocatable::  z1p(:,:)
      integer ::  nbnb(nqbz,npm), 
     &  n1b(nbnbx,nqbz,npm), n2b(nbnbx,nqbz,npm)
      complex(8),allocatable:: zzmel(:,:,:)
      logical,parameter:: debug=.false.
      logical:: hist ,ipr
      integer:: nhwtot,
     &  ihw(nbnbx,nqbz,npm),nhw(nbnbx,nqbz,npm),jhw(nbnbx,nqbz,npm)
      real(8):: whw(nhwtot)
      complex(8) :: zmelt1,zmelt2,zmeltt(ngbb) 
      real(8) :: imagweight 
      integer::verbose,nmbas, imb1,imb2, imb 
      real(8):: eband(nband)
      logical   :: nolfco 
      complex(8):: zq01,zq02,zq0zq0 
      real(8):: vec_kq_g(3),vec_k_g(3),vec_kq(3),vec_k(3),quu(3),tolq=1d-8,quu1(3),quu2(3) 
      integer:: nbcut,nbcut2
      logical :: iww1=.true.,iww2=.true.
      integer:: ifpomat, nbloch_r, ngbo,iqxdummy
      complex(8),allocatable:: pomat(:,:), zmeltn(:,:,:)
      real(8):: q_r(3)
      complex(8):: img=(0d0,1d0),zzz(ngbb)
      integer :: nkmin,  nkmax, nkqmin, nkqmax,nkmax1,nkqmax1, ib1, ib2, ngcx,ix,iy,igb
      logical :: onceww , z1offd 
      complex(8),target :: zzr(ngbb,nmbas) 
      integer:: ngrp,nwgt(nqbz) 
      integer:: igx(ngrp*2,nqbz),igxt(ngrp*2,nqbz),ieqbz
      logical:: checkbelong,eibzmode, chipmzzr
      complex(8):: zcousq(ngbb,ngbb) 
      complex(8),allocatable:: zcousqr(:,:,:),rcxq0(:,:),rcxq00(:,:),rcxq000(:,:),rcxqwww(:,:)
      complex(8):: zmeltx(ngbb),zmelty(ngbb),zcousqrx(ngbb,ngbb) ,zcousqc(ngbb,ngbb)
     & ,rzc(ngbb,ngbb),cmat(ngbb,ngbb) 
      integer::  eibzsym(ngrp,-1:1),neibz,icc,ig,eibzmoden,ikp,i,j,itimer,icount,iele
      integer:: irotm,nrotmx,ixx,iyy,itt,ntimer, nccc, nxx,iagain,irotm1,irotm2
      integer,allocatable:: i1(:,:),i2(:,:),nrotm(:)
      complex(8),allocatable::zrotm(:,:),zrr(:,:),zrrc(:,:),zrr_(:,:,:),zrrc_(:,:,:)
     &  ,zmmm(:),zrrx(:,:),rcxq_core(:,:)
      logical:: eibz4x0, crpa , exchange=.false.
      real(8):: wpw_k,wpw_kq, vec_kcrpa(3),vec_kqcrpa(3)
      integer:: irot=1, ntqxx,nbmax

!TIME0_1001 
      write(6,'(" x0kf_v4hz: q=",3f8.4,$)')q
      call cputid(0)
      call minv33(qbas,qbasinv)
      allocate(cphi_k(nlmto,nband), cphi_kq(nlmto,nband), geig_kq(ngpmx,nband), geig_k(ngpmx,nband))
      call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
      call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )

!TIME1_1001 "beforedo1000"
!! loop over k-points -------------------------------------------------------------------------
      do 1000 k = 1,nqbz
        if(eibzmode.and.nwgt(k)==0 ) cycle
        if(debug) write(6,'("do 1000  k=",i4,3f10.4)')k,rk(:,k)
        ipr=(k<5.or.k==nqbz.or.debug)
        if(sum(nbnb(k,1:npm))==0) cycle
!TIME0_1101
        if(k<=5.or. (mod(k,max(10,nqbz/20))==1.or.k>nqbz-10) ) then
          write(6,"('  x0kf_v4hz: k rk=',i7,3f10.4,$)")k, rk(:,k)
          call cputid(0)
        endif
!! tetra 
        nkmin = 999999
        nkmax= -999999
        nkqmin= 999999
        nkqmax=-999999
        do jpm=1,npm !npm
          do ibib = 1, nbnb(k,jpm)
            nkmin  = min(n1b(ibib,k,jpm),nkmin)
            nkqmin = min(n2b(ibib,k,jpm),nkqmin)
            if(n1b(ibib,k,jpm)<=nband) nkmax  = max(n1b(ibib,k,jpm),nkmax)
            if(n2b(ibib,k,jpm)<=nband)   nkqmax = max(n2b(ibib,k,jpm),nkqmax)
          enddo
        enddo

        itps  = nkqmin          ! nkqmin = the num of min   n2 =unocc for jpm=1
        nt0   = nkmax
        ntp0  = nkqmax - nkqmin +1
        if( npm==2.and. nkqmin/=1) then
          write(6,*)' npm==2 nkqmin nkqmax  nkmin nkmax=',nkqmin,nkqmax,nkmin,nkmax
          call rx( " When npm==2, nkqmin==1 should be.")
        endif
        if(nkmin/=1) call rx( " nkmin==1 should be.")

!  zzmel(1:nbloch, ib_k,ib_kq)
!      ib_k =[1:nctot]              core
!      ib_k =[nctot+nkmin:nctot+nkmax]  valence
!      ib_kq =[1:ncc]             core
!      ib_kq =[ncc+nkqmin:ncc+nkqmax]  valence range [nkqmin,nkqmax]
!   If jpm=1, ncc=0.
!   If jpm=2, ncc=ncore. itps=1 should be.
! There is a little confusion. n1b index contains cores are after valence.
! You can see codes to treat the confusion.
! NOTE: 
!  q+rk n2b vec_kq  vec_kq_g geig_kq cphi_kq  ngp_kq ngvecp_kq  isp_kq
!    rk n1b vec_k   vec_k_g  geig_k  cphi_k   ngp_k  ngvecp_k   isp_k


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
        ngb  = nbloch + ngc 
!!Get the matrix element zmel   ZO^-1 <MPB psi|psi> , where ZO is ppovlz
!!  Output is zmel(ngb, nctot+nt0,ncc+ntp0)  nkmin:nt0, nkqmin:ntp0
! nt0=nkmax-nkmin+1  , ntp0=nkqmax-nkqmin+1
        call get_zmelt2(exchange,  
     &      q,irot,q,ngc,ngb,                   ! MPB
     &        nkmin, nkmax, isp_k,nctot,        
              ! middle state 1:nt0   --> true index of eigen is mkmin:mkmin+nt0-1     + nctot
     &      q+rk(:,k),nkqmin,nkqmax,isp_kq,ncc) 
              ! end state    1:ntp0  --> true index of eigen is nkqmin:nkqmin+ntp0-1  
        call Dconjg_zmel() !zmel = dconjg(zmel) 
        allocate( z1p(ngb,ngb) )
 8829   continue


!TIME1_1101 "before_get_zmelt2"
!TIME0_1201
c----------------------------
!!  z1p = <M_ibg1 psi_it | psi_itp> < psi_itp | psi_it M_ibg2 >
!!  zxq(iw,ibg1,igb2) = sum_ibib wwk(iw,ibib)* z1p(ibib, igb1,igb2)

        do 25 jpm  = 1, npm !
        do 25 ibib = 1, nbnb(k,jpm) !---  ibib loop

!! n1b,n2b --> core after valence.  it,itp --> valence after core
          if(n1b(ibib,k,jpm) <= nband) then
            it = nctot + n1b(ibib,k,jpm) !valence
            if(it > nctot + nkmax ) cycle
          else
            it = n1b(ibib,k,jpm) - nband !core
          endif
          if( n2b(ibib,k,jpm) <= nband) then
            itp = ncc + n2b(ibib,k,jpm) - itps + 1 !val
            if(itp > ncc + nkqmax-itps+1 ) cycle
          else
            itp =  n2b(ibib,k,jpm) - itps + 1 - nband !core
          endif

!! nbcut mechanism
          if(jpm==1) then 
            if( n1b(ibib,k,jpm) <= nbcut .and. nbcut2<n2b(ibib,k,jpm) ) then
             if(iww2) then
              write(6,"(' nband_chi0 nbcut nbcut2 n2b n1b=',4i6)") nbcut,n2b(ibib,k,jpm),n1b(ibib,k,jpm)
              iww2=.false.
             endif
             cycle
            endif
          else !jpm==2 
            if( n2b(ibib,k,jpm) <= nbcut .and. nbcut2<n1b(ibib,k,jpm) ) then
             if(iww2) then
              write(6,"(' nband_chi0 nbcut nbcut2 n2b n1b=',4i6)") nbcut,n2b(ibib,k,jpm),n1b(ibib,k,jpm)
              iww2=.false.
             endif
             cycle
            endif
          endif
!!
          if (ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1 >nwt) call rx( "x0kf_v4hz: iw>nwt")
          do igb2=1, nmbas     !....................................
             zmelt2 = zmel(igb2,it,itp) !zmelt(igb2,it,itp)
             do igb1=1,igb2
                z1p(igb1,igb2) = dconjg(zmel(igb1,it,itp)) * zmelt2
             enddo
          enddo
!! ---------------------------------
          if(crpa) then
            if(n1b(ibib,k,jpm) <= nband) then
              call readpkm4crpa(n1b(ibib,k,jpm),   rk(:,k), isp_k,   wpw_k) !k  n1b
            else
              wpw_k=0d0
            endif  
            if(n2b(ibib,k,jpm) <= nband) then
              call readpkm4crpa(n2b(ibib,k,jpm), q+rk(:,k), isp_kq,  wpw_kq) !kq n2b
            else
              wpw_kq=0d0
            endif  
          endif
c$$$cccccccccccccccccccccccccccccccccc
c$$$ For SrVO3 test
c$$$          if(crpa) then
c$$$             wpw_k=0d0
c$$$             if(15<n1b(ibib,k,jpm).and.n1b(ibib,k,jpm)<19) wpw_k=1d0
c$$$             wpw_kq=0d0
c$$$             if(15<n2b(ibib,k,jpm).and.n2b(ibib,k,jpm)<19) wpw_kq=1d0
c$$$          endif
c$$$ccccccccccccccccccccccccccccccccccc

!$OMP parallel private(imagweight) 
!$OMP master
!$        if (jpm.eq.1 .and. ibib.eq.nbnb(k,1)) then 
!$          write(6,'(a,i5,a,i5)') 'OMP parallel iw, threads=', omp_get_num_threads(),
!$   .     ' nw=',ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-ihw(ibib,k,jpm)
!$        endif
!$OMP end master
!$OMP do 
          do iw = ihw(ibib,k,jpm),ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1 !iiww=iw+ihw(ibib,k)-1
             imagweight = whw(jhw(ibib,k,jpm)+iw-ihw(ibib,k,jpm))
             if(crpa) imagweight = imagweight*(1d0-wpw_k*wpw_kq)
             if(eibzmode) imagweight = nwgt(k)*imagweight  
             do igb2=1,nmbas      !this part dominates cpu time most time consuming...........
                do igb1=1,igb2
                   rcxq(igb1,igb2,iw,jpm) = !here we  sum over ibib (or n, n') and k.
     &                  rcxq(igb1,igb2,iw,jpm) + z1p(igb1,igb2)*imagweight !sum over spin in hx0fp0
                enddo           !igb1
             enddo              !igb2
          enddo                 ! iw
!$OMP end do 
!$OMP end parallel 
 25     continue
!kino 2014-08-13  !$OMP end parallel 
!TIME1_1201 "after_rcxq"
        deallocate(z1p)         !zmelt,z1p)
        call Deallocate_zmel()
        if(debug) call cputid(0)
        if(debug) write(6,*)' end of kloop k jpm=',k,jpm
 1000 continue
!TIME0_1301

!! ==== Hermitianize. jun2012takao moved from dpsion5 ====
        do jpm=1,npm
        do iw= 1,nwt
         do igb2= 1,nmbas !eibzmode assumes nmbas1=nmbas2
         do igb1= 1,igb2-1
           rcxq(igb2,igb1,iw,jpm) = dconjg(rcxq(igb1,igb2,iw,jpm))
         enddo
         enddo
        enddo
        enddo
!TIME1_1301 "before_eibzmode_symmetrization"
 9999 continue
      write(6,"(' --- x0kf_v4hz: end')") !, 3d13.5)")
      if(debug) write(6,"(' --- ', 3d13.5)")
     &     sum(abs(rcxq(1:nmbas,1:nmbas,1:nwt,1:npm))),sum((rcxq(1:nmbas,1:nmbas,1:nwt,1:npm)))
      end subroutine x0kf_v4hz

!! --------------------------------------------------------------------------------
      subroutine x0kf_v4hz_symmetrize (npm, !ncc,
c     i            ihw,nhw,jhw,whw,nhwtot, ! tetwt5
c     i            n1b,n2b,nbnbx,nbnb,     !  use whw by tetwt5 ,
     i                  q,
     i                  nsp,isp_k,isp_kq, !symmetrize,
     i                  qbas,ginv,!rk,wk,  
c     i                  mdim, 
c     d                  nlmto,nqbz,nctot, 
c     d                  natom, 
     d                  nbloch,nwt, 
     i      iq, ngbb, ngc, ngpmx,ngcmx,  
     i      nqbze, nband,nqibz, 
     o      rcxq,
     i     nolfco,zzr,nmbas, zcousq,
     i     chipmzzr,eibzmode, 
     i     ngrp,eibzsym) !, crpa) !nwgt,igx,igxt,
      use m_rotMPB,only   :  rotMPB2 
      use m_readqgcou,only:  
     o qtt_, nqnum 
!! === symmetrization for EPIBZ mode ===
      implicit none
      integer(4):: npm,ncc,ngbb,natom,nwt,nsp,isp_k,isp_kq 
     &  ,nl,nclass,nnc,nbloch,iq,nqibz,iatom,nctot,nbmx,iopen 
     &  ,jpm,ibib,itps,nt0,ntp0,ngp_kq,ngp_k,it,itp,iw,igb2,igb1,ngb
     &  ,nn,no,isx,iclose,k,nbnbx,nqbz
      real(8):: q(3),qbas(3,3),ginv(3,3),ebmx !,rk(3,nqbz),wk(nqbz)
      complex   (8):: rcxq (nmbas,nmbas,nwt,npm)
      complex(8) :: imag=(0d0,1d0),trc,aaa !phase(natom),
      integer(4):: ngpmx, ngcmx, nqbze, nband,
     &              ngc,nadd(3)
c      integer(4) ::  nbnb(nqbz,npm), n1b(nbnbx,nqbz,npm), n2b(nbnbx,nqbz,npm)
      logical,parameter:: debug=.false.
c      integer(4):: nhwtot, ihw(nbnbx,nqbz,npm),nhw(nbnbx,nqbz,npm),jhw(nbnbx,nqbz,npm)
c      real(8):: whw(nhwtot)
      complex(8) :: zmelt1,zmelt2,zmeltt(ngbb)      !...........................sf 21May02
      real(8) :: imagweight !............................sf 21May02
      real(8):: eband(nband)!,ebandr(nband),ebandqr(nband)
      integer(4):: verbose

      logical   :: nolfco !iepsmode
      integer(4):: nmbas, imb1,imb2, imb !nmbas1x !nmbas2,nmbas1,
      real(8):: vec_kq_g(3),vec_k_g(3),vec_kq(3),vec_k(3),quu(3),tolq=1d-8,quu1(3),quu2(3) !tolqu=1d-4,
      integer(4):: nbcut,nbcut2
      logical :: iww1=.true.,iww2=.true.
      complex(8):: img=(0d0,1d0)
      integer(4):: nkmin,  nkmax, nkqmin, nkqmax,nkmax1,nkqmax1
      integer(4):: ib1, ib2,      ngcx,ix,iy
      complex(8),target :: zzr(ngbb,nmbas) !ppovlz(ngbb,ngbb),
      integer:: igb
      integer:: ngrp
      logical:: checkbelong,eibzmode, chipmzzr
      complex(8):: zcousq(ngbb,ngbb) ,  zcousqc(ngbb,ngbb)
      integer::  eibzsym(ngrp,-1:1),neibz,icc,ig,eibzmoden,ikp,i,j,itimer,icount,iele
      integer:: irotm,nrotmx,ixx,iyy,itt,ntimer, nccc, nxx,iagain,irotm1,irotm2
      integer,allocatable:: i1(:,:),i2(:,:),nrotm(:)
      complex(8),allocatable:: zrotm(:,:),zrr(:,:),zrrc(:,:),zrr_(:,:,:),zrrc_(:,:,:),zmmm(:),zrrx(:,:)
      complex(8),allocatable:: rcxq_core(:,:)
      complex(8),allocatable:: zcousqr(:,:,:),rcxq0(:,:),rcxq00(:,:),rcxq000(:,:),rcxqwww(:,:)

      logical :: exchange=.false.
      integer:: irot=1
      integer:: ntqxx,nbmax
    
!! ------------------------------------------------------------------------
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
c 5000 continue
!! === zmelt conversion === 
      if(nolfco .and. nmbas==1) then
         write(6,*)' nmbas=1 nolfco=T ---> not need to symmetrize'
         goto 9999
      endif
!!
      if(eibzmode) then
        ngb  = nbloch + ngc 
        if(ngb/=ngbb) then
           write(6,*)' x0kf_v4h: ngb ngbb=',ngb,ngbb
           call rx( 'x0kf_v4h: ngb/=ngbb')
        endif
!TIME0_1401 
        call iqindx2(q, ginv, qtt_, nqnum, ikp,quu) !to get ikp for timereversal mode
!TIME1_1401 "after_iqindx2"
!TIME0_1501
        if(sum(abs(q-quu))>tolq) call rx( 'x0kf_v2h: eibz 111 q/quu')
        neibz = sum(eibzsym(:,1))+sum(eibzsym(:,-1)) 
             !itimer=-1 means time reversal. eibzsym(ig,itimer) where ig: space rotation.
        write(6,"(' --- goto symmetrization --- ikp neibz q=',2i3,3f12.8)")ikp,neibz,q
        call cputid2(' --- x0kf: start symmetrization  ',0)

c        allocate(rcxq0(ngb,ngb),rcxq00(ngb,ngb),rcxq000(ngb,ngb),rcxqwww(ngb,ngb))
        ntimer=1
        if(sum(eibzsym(:,-1))>0) ntimer=2 !timereversal case
        allocate(zrotm(ngb,ngb),nrotm(ngrp*2))
!!
c        zcousqinv=zcousq
c        call matcinv(ngb,zcousqinv)

!! == Assemble rotantion matrx zrr,zrrc ==
!! Rotation matrix zrrx can be a sparse matrix. 
!! Thus it is stored to  "i1(nrotmx,nccc),i2(nrotmx,nccc),zrr(nrotmx,icc),nrotm(icc)".
!! See folloings: matmul(rcxqwww,zrrx) is given by
!!     do irotm1 = 1,nrotm(icc)
!!       rcxq0(:,i2(irotm1,icc)) = rcxqwww(:,i1(irotm1,icc)) * zrr(irotm1,i2(irotm1,icc))

        allocate(zrrx(nmbas,nmbas))
        nrotmx = 10000 !trial value
!TIME1_1501 "before_1011"
        do 1011  !this loop is only in order to to set large enough nrotmx.
!TIME0_1601 
        if(allocated(i1)) deallocate(i1,i2,zrr,zrrc)!,zrr_,zrrc_)
        nccc=ngrp*2
        allocate(i1(nrotmx,nccc),i2(nrotmx,nccc),zrr(nrotmx,nccc),zrrc(nrotmx,nccc)) !,zrr_(ngb,ngb,nccc),zrrc_(ngb,ngb,nccc))
        i1=-99999
        i2=-99999
        zrr=-99999d0
        zrrc=-99999d0
        call cputid2(' --- x0kf:11111   :',0)
!TIME1_1601 "allocatezrr"
!!
        icc=0
        do itimer=1,-1,-2
          if(ntimer==1.and.itimer==-1) exit
          if(itimer==1 ) itt=1
          if(itimer==-1) itt=2
          do ig=1,ngrp
            if(eibzsym(ig,itimer)==1) then
              icc=icc+1
!TIME0_1701

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
                 call rx( 'x0kf_v4h: this case is not implemented xxxxxxxxxxxxxx')
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
!TIME1_1701 "end_matmmsparse"
!TIME0_1801
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
!TIME1_1801 "before_iagain1"
!TIME0_1901
              if(iagain==1) then
                nrotmx=irotm !enlarge allocation and do things again.
                write(6,*)' warn:(slow speed) xxxx goto 1011 xxxxxx nrotmx+=nrotmx+10000 again'
                goto 1011
                !enlarge nrotmx ang try it again.
              endif   
              nrotm(icc)=irotm
              if(debug) write(6,*)'ig itimer icc nrotm=',ig,itimer,icc,nrotm(icc) ,iele
!TIME1_1901 "end_ig_itimer_icc_nrotm"
            endif
          enddo
        enddo
        exit
 1011   continue !only when nrotmx overflow.
!TIME0_2001

!! === main part to obtain symmetrized rcxq  ===
!! neibz is total number of symmetrization operation. 
!!      rcxq is rotated and accumulated; finally divied by neibz
        zcousqc = dconjg(transpose(zcousq))
        if(debug) call cputid2(' --- x0kf:qqqqq222ini:',0)
!$OMP parallel private(rcxq000,icc,itt,icount,rcxqwww,rcxq00,rcxq0,rcxq_core)
        allocate(rcxq0(ngb,ngb),rcxq00(ngb,ngb),rcxq000(ngb,ngb),rcxqwww(ngb,ngb),rcxq_core(ngb,ngb))
!$OMP master
!$         write(6,'(a,i5,a,i5)') 'OMP parallel nwt, threads=',omp_get_num_threads(),' nwt=',nwt
!$OMP end master
!$OMP do 
        do iw=1,nwt
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
!TIME1_2001 "after_sym_rcxq"
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
c$$$        do iw=1,nwt
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
!kino 2014.08.19 use automatic deallocation,      deallocate(cphi_k,cphi_kq,geig_kq,geig_k)
      write(6,"(' --- x0kf_v4hz_symmetrize: end')")
      end subroutine x0kf_v4hz_symmetrize




C=====================================================================
        subroutine dpsion5  (frhis,nwhis, freqr,nw_w, freqi,niwt, 
     i                                 realomega,   imagomega,      !freqr ->frhis ...sf
     i                    rcxq, npm, nw_i,nmbas1,nmbas2,
     o    zxq,zxqi,
c     i    nolfco,chipm,schi,isp,  rcxqmean,nmbas,  !iepsmode, rcxqmean, ! epsmode
     i    chipm,schi,isp,   !No nolfco mode. Apr2012.
     i    ecut,ecuts)
c     o    x0mean)
C- Calculate W-v zxqi(on the imaginary axis) and zxq(real axis) from sperctum weight rcxq.
Cr v4 works for timereversal=F (npm=2 case).
Cr  See rcxq_zcxq for rcxq, which contains the spectrum weight for given bins along the real-axis.
Cr ! Note that zxq and zxqi are not accumlating
Ci frhis(1:nwhis+1) :: specify histgram bins i-th bin is [frhis(i), frhis(i+1)].
Ci          We suppose "freqr(i)=moddle of i-th bin; freqr(0)=0."
Ci          (I think called routine hilbertmat itself is not limited by this condition).
Ci freqr (0:nw_w) : Calcualte zxq for these real energies.
Ci freqi (1:niwt) : Calcualte zxqi for these imaginary energies.
Ci   realomega  : A switch to calculate zxq or not.
Ci   imagomega: : A switch to calculate zxqi or not.
Ciw rcxq may be altered ---used as work area.
Cio   zxq :  W-v along the real axis on freqr(0:nw_w)
Cio   zxqi:  W-v along the imag axis on freqi(niwt)
C!
C1 Feb2006:  v4 for timereversal=F
C! July2005: v3Add spin chipm mode
C! July2005: This version alter rcxq----it is used as work area.
C! sergey faleev Apr 2002 ; Rebuiled by takao
C------------------------------------------------------------------
        implicit none
        integer(4):: nw_w,niwt,igb1,igb2, iw,iwp,nwhis,ix,npm,ifxx,nmbas1,nmbas2
        real(8) :: freqi(niwt),pi,px,omp,om,om2,om1, !omg2max from hx0fp0
     &  frhis(nwhis+1), freqr(0:nw_w), aaa,d_omg
        logical :: realomega, imagomega
        complex(8):: rcxq(nmbas1,nmbas2, nwhis,npm) !sf 13June
c      logical   :: iepsmode
        logical :: chipm

        integer(4)::isp,ispx !, nmbas
c      complex(8):: rcxqmean(nwhis,npm,nmbas,nmbas)  !takao sep2006 add nmbas
C... ecut mode
        real(8):: ecut,ecuts,wcut,wcutef,dee,schi
        logical ::debug=.false.
        real(8),allocatable :: his_L(:),his_R(:),his_C(:)
        integer(4) it
        real(8):: domega_r,domega_c,domega_l,delta_l,delta_r
        real(8),allocatable ::rmat(:,:,:),rmati(:,:,:),rmatt(:,:,:),imatt(:,:,:)
        complex(8),allocatable :: rmatiC(:,:,:),imattC(:,:,:)
        complex(8) ::beta,wfac
        complex(8):: zz
        complex(8),allocatable :: zxqn(:,:),zxqn1(:,:,:),rx0mean1(:,:,:),rx0mean(:)
        complex(8),allocatable:: rrr(:)

        integer(4)::nw_i,jpm,ipm,verbose,isgi
c      complex(8):: x0mean(nw_i:nw_w,nmbas,nmbas)
        complex(8):: 
     o   zxq (nmbas1,nmbas2, nw_i: nw_w), !iw=0 means omg=0,
        !iw=1:nw_w corresponds to iw's bit of the frequensy histogram 
     o   zxqi(nmbas1,nmbas2,niwt),img !npm), img  !zxqi(...,npm) may2006

        real(8),allocatable:: ebb(:)
        integer(4):: ii,i,ibas1,ibas2
        logical :: evaltest !,testtr

c        if(verbose()>89) debug=.true.
c ------------------------------------------------
        write(6,'(" -- dpsion5: start...   ",$)')
        write(6,"('  nw_w nwhis=',2i5)") nw_w,nwhis
        if(debug) then
          write(6,*)' nmbas1 nmbas2 nwhis npm =',  nmbas1,nmbas2,nwhis,npm
          write(6,*)' sumchk rcxq=', sum(abs(rcxq))
        endif
        pi  = 4d0*datan(1d0)
        img = (0d0,1d0)
        call cputid(0)
        ispx = isp
        if(schi<0) then
          ispx = 3-isp !flip
        endif

!! Check freqr
        if(realomega) then
          if( nwhis <= nw_w ) then
            write(6,*)nwhis,nw_w
            call rx( ' dpsion5: nwhis<=nw_w')
          endif
          if( freqr(0)/=0d0 ) call rx( ' dpsion5: freqr(0)/=0d0')
!! I think current version allows any freqr(iw), independent from frhis.
c$$$          aaa = 0d0
c$$$          if(nw_w>0) then
c$$$            do iw = 1,nw_w
c$$$              aaa = aaa + abs( freqr(iw) - (frhis(iw)+frhis(iw+1))/2d0 )
c$$$              if(debug) write(6,"(' iw freqr frhis_m=',i5,2f13.6)" )
c$$$     &        iw,freqr(iw),  (frhis(iw)+frhis(iw+1))/2d0
c$$$            enddo
c$$$            if(aaa>1d-10)call rx( 'dpsion5:freqr/=frhis_m is not implimented yet')
c$$$          endif
        endif !realomega

C--------------------------------------------------------------
!! Each histogram bins are  [his_Left, his_Right], and  his_Center is middle.
!! his_C(0) is at zero. his_R(0) and his_L(0) are not defined.
        if(debug) write(6,*)' dpsion5: RRR 2222222222 '
        allocate(his_L(-nwhis:nwhis),his_R(-nwhis:nwhis),his_C(-nwhis:nwhis))
        his_L(1:nwhis) = frhis(  1:  nwhis)
        his_R(1:nwhis) = frhis(1+1:1+nwhis)
        his_C(1:nwhis) = (his_L(1:nwhis) + his_R(1:nwhis) )/2d0
        do iw= 1,nwhis
          his_L(-iw) = -his_R(iw)
          his_R(-iw) = -his_L(iw)
          his_C(-iw) = -his_C(iw)
        enddo
        his_C(0) = 0d0; his_R(0)=-999; his_L(0)=-999
C
        if(debug) write(6,*)'sumchk 111 rcxq=', sum(abs(rcxq))

        do iw= 1, nwhis
          if(ecut<1d9) then
            wfac= wcutef(his_C(iw), ecut,ecuts)
          else
            wfac= 1d0
          endif
! rcxq is used as work---> rcxq= Average value of Im chi.
! Note rcxq is "negative" (
          do jpm=1,npm
            call dscal(2*nmbas1*nmbas2, -wfac/(his_r(iw)-his_l(iw)),rcxq(1,1,iw,jpm),1)
          enddo
c         if(debug) write(6,*) 'dpsion5: RRR 7777 iw wfac=',iw,wfac,ecut,ecuts
        enddo
        if(debug) write(6,*)'sumchk 122 rcxq=', sum(abs(rcxq))

C... Temporary. maybe, we will have better procedure...
ctakao moved this to x0kv_v4h.F jun2012takao
!! hermitianize.
c        if(nmbas1==nmbas2) then !Is this required??? apr2012takao
c          do jpm=1,npm
c            do iw= 1, nwhis
c              do igb2= 1, nmbas2
c                do igb1= 1, igb2-1
c                  rcxq(igb2,igb1,iw,jpm) = dconjg(rcxq(igb1,igb2,iw,jpm))
c                enddo
c              enddo
c            enddo
c          enddo
c        endif
cccccccccccccccccccc
        if(debug) write(6,*)'sumchk 222 rcxq=', sum(abs(rcxq))

        if(evaltest().and.nmbas1==nmbas2) then
          write(6,"('hhh --- EigenValues for rcxq --------')")
          allocate(ebb(nmbas1))
          do jpm= 1,npm
            do iw = 1, nwhis
              call diagcvh2(rcxq(:,:,iw,jpm),nmbas1,ebb)
              do ii=1,nmbas1
                write(6,"('hhh1: xxxxxxxxxxxxxxxxx',2i4)") jpm,iw
                if(abs(ebb(ii))>1d-8.and.ebb(ii)>0) 
     &        write(6,"('hhh1: jpm iw eb=',2i4,d13.5)") jpm,iw,ebb(ii)
              enddo
            enddo
          enddo
          deallocate(ebb)
        endif

C--- realomega case
        if(realomega)then
          write(6,*) " --- realomega --- "
          if(npm==1) then
            allocate( rmat(0:nw_w,-nwhis:nwhis,npm), rrr(-nwhis:nwhis))
            rmat  = 0d0
            do it =  0,nw_w
              zz = freqr(it) !his_C(it)
              call hilbertmat(zz, nwhis,his_L,his_C,his_R, rrr)
              rmat(it,:,1) = dreal(rrr)
            enddo ;   if(debug) write(6,*) 'dpsion5: RRR 55555555555'
            allocate( rmatt(0:nw_w,nwhis,npm) )
            if(     chipm.and.ispx==1 ) then
              rmatt(:,:,1) = rmat(:,1:nwhis,1)
            elseif( chipm.and.ispx==2 ) then
              do iw= 1,nwhis
                rmatt(:,iw,1) = -rmat(:,-iw,1)
              enddo
            else
              do iw= 1,nwhis
                rmatt(:,iw,1) = rmat(:,iw,1) - rmat(:,-iw,1)
              enddo
            endif
            deallocate(rmat,rrr)
          else ! npm==2 case -------------------------------------------------
            allocate( rmatt(-nw_w:nw_w,nwhis,npm), rrr(-nwhis:nwhis))
            rmatt = 0d0
            do it  =  -nw_w,nw_w
              if(it<0) then
                zz = -freqr(-it) !his_C(it)
              else
                zz = freqr(it) !his_C(it)
              endif
              call hilbertmat(zz, nwhis,his_L,his_C,his_R, rrr)
              rmatt(it,:,1) =  dreal(rrr  (1:nwhis))
              rmatt(it,:,2) = -dreal(rrr(-1:-nwhis:-1))
            enddo ;   if(debug) write(6,*) 'dpsion5: RRR2 55555555555'
            deallocate(rrr)
          endif
          rmatt = rmatt/pi ; if(debug) write(6,*)'dpsion5: RRR 6666'

!! takao remove if(nolfc) block here.
c          write(6,*) " --- realomega dgemm--- "


!! WARN! I think npm==2.and.chipm does not make sense. apr2012.
!!
          if(npm==2.and.chipm) 
Cstop2rx 2013.08.09 kino     &      stop 'x0kf_v4h:npm==2.and.chipm is not meaningful probably'
     &      call rx( 'x0kf_v4h:npm==2.and.chipm is not meaningful probably')


!! Note rcxq is negative now (converted at the top of this routine !!!
          if(     chipm .and. ispx==2 ) then
            !nothing here
            !Since the range of zxq is nw_i=0, we have no area to store negative energy part of chipm.
          elseif( chipm               ) then
            call zaxpy( nmbas1*nmbas2*nw_w, img, rcxq, 1, zxq(1,1,1), 1)
          else
            zxq = 0d0    ! not accumlating case.
            call zaxpy( nmbas1*nmbas2*nw_w, img, rcxq(1,1,1,1), 1, zxq(1,1,1), 1)
          endif

          if(npm==2) then
            do iw=1,nw_w
              call zaxpy( nmbas1*nmbas2, img, rcxq(1,1,iw,2),1, zxq(:,:,-iw),1)
            enddo
          endif

          if(npm==1) then
            call dgemm('n','t',  2*nmbas1*nmbas2, nw_w+1, nwhis, 1d0, 
     &              rcxq, 2*nmbas1*nmbas2,  rmatt, nw_w+1,
     &               1d0, zxq, 2*nmbas1*nmbas2 )
          elseif(npm==2) then
            call dgemm('n','t',  2*nmbas1*nmbas2,   npm*nw_w+1, nwhis, 1d0, 
     &              rcxq(1,1,1,1), 2*nmbas1*nmbas2, rmatt(:,:,1), npm*nw_w+1,
     &               1d0, zxq, 2*nmbas1*nmbas2 )
            call dgemm('n','t',  2*nmbas1*nmbas2,   npm*nw_w+1, nwhis, 1d0, 
     &              rcxq(1,1,1,2), 2*nmbas1*nmbas2, rmatt(:,:,2), npm*nw_w+1,
     &               1d0, zxq, 2*nmbas1*nmbas2 )
          else
Cstop2rx 2013.08.09 kino            stop 'dpsion5: npm=1 or 2'
            call rx( 'dpsion5: npm=1 or 2')
          endif
          deallocate(rmatt)
        endif

!! === imagomega case      imatt(niwt -->niwt,npm may2005 ===
        if(imagomega) then
          allocate( rrr(-nwhis:nwhis))
          if(npm==1) then
            allocate( rmati (niwt,-nwhis:nwhis,npm))
            rmati= 0d0
          else
            allocate( rmatiC(niwt,-nwhis:nwhis,npm))
            rmatiC = 0d0
          endif ;   if(debug) write(6,*) 'dpsion5: III 111111155555555555'
          do it =  1,niwt
            zz = img*freqi(it) !his_C(it)
            call hilbertmat(zz, nwhis,his_L,his_C,his_R, rrr) !Im(zz)>0
            if(npm==1) then
              rmati (it,:,1) = dreal(rrr)
            else
              rmatiC(it,:,1) = rrr
            endif
          enddo ;   if(debug) write(6,*) 'dpsion5: III 55555555555'
!! ==== npm=1 case ====
          if(npm==1) then
            allocate( imatt(niwt, nwhis,npm) )
            do iw= 1,nwhis
              imatt(:,iw,1) = rmati(:,iw,1) - rmati(:,-iw,1)
            enddo
            deallocate(rmati,rrr)
            imatt = imatt/pi; if(debug) write(6,*) 'dpsion5: III  '
            call dgemm('n','t',  2*nmbas1*nmbas2, niwt, nwhis, 1d0, 
     &              rcxq, 2*nmbas1*nmbas2, imatt, niwt,
     &               0d0, zxqi, 2*nmbas1*nmbas2 )
            deallocate(imatt)
!! ==== npm=2 case ====
          else
            allocate( imattC(niwt, nwhis,npm) )
            do iw= 1,nwhis
              imattC(:,iw,1) =   rmatiC(:, iw,1)
              imattC(:,iw,2) = - rmatiC(:,-iw,1)
            enddo
            deallocate(rmatiC,rrr)
            imattC = imattC/pi; if(debug) write(6,*) 'dpsion5: IIIc '
            call zgemm('n','t',  nmbas1*nmbas2, niwt, nwhis, 1d0, 
     &              rcxq(1,1,1,1), nmbas1*nmbas2, imattC(1,1,1), niwt,
     &              0d0, zxqi,    nmbas1*nmbas2 )
            call zgemm('n','t',  nmbas1*nmbas2, niwt, nwhis, 1d0, 
     &              rcxq(1,1,1,2), nmbas1*nmbas2, imattC(1,1,2), niwt,
     &              1d0, zxqi,    nmbas1*nmbas2 )
            deallocate(imattC)
          endif
        endif
        deallocate(his_L,his_C,his_R)
        write(6,'("         end dpsion5 ",$)')
        call cputid(0)
        end
        logical function checkbelong(qin, qall, nq,ieibz) !ieibz is also returned
        integer:: nq,ieibz
        real(8):: qin(3), qall(3,nq),tolq=1d-8
        checkbelong=.false.
        do i=1,nq
           if(sum(abs(qin-qall(:,i)))<tolq) then
              ieibz=i
              checkbelong=.true.
              return
           endif
        enddo
        end


!!----------------------------------------------------------------------------
      subroutine hilbertmat (zz,nwhis, his_L,his_C,his_R, rmat)
C- Martix for hilbert transformation, rmat.
Cr  zz is real--->  no img*delta function part
Cr   zz is complex (and Im(zz)>0) : includes all contribution when Im(zz)>eps
Co  rmat(-nwhis:nwhis) : rmat(0) is not meaningful.
Ci i-th Histgram bin on real axis are given by [his_L, his_R]. center is his_C.
Cr f(zz) = \int_-x(nwhis)^x(nwhis) f(x)/(zz-x)
Cr       = \sum_{i/=0} rmat(i)*f(i)
Cr     ,where f(i) is the average value at i-th bin.
C!!! 23May2006 I think
C!!! rmat is --------------
C!!! f(zz) = - \int_-x(nwhis)^x(nwhis) f(x)/(zz-x)
C!!!       = - \sum_{i/=0} rmat(i)*f(i)
C I forgot minus sign in the previous note.
C-------------------------------
      implicit none
      integer(4):: iw,nwhis
      complex(8) ::zz,imgepsz
      real(8)    :: his_L(-nwhis:nwhis),his_C(-nwhis:nwhis),his_R(-nwhis:nwhis)
      complex(8) :: rr_fac(-nwhis:nwhis),rl_fac(-nwhis:nwhis),img=(0d0,1d0)
      real(8)::  eps=1d-8, epsz=1d-13,delta_r,delta_l,ddr,ddl
      complex(8):: domega_c,domega_r,domega_l
      complex(8) ::  rmat(-nwhis:nwhis)
      imgepsz =img*epsz
      do iw = -nwhis, nwhis
        if(iw==0) cycle
        domega_r = zz - his_R(iw) + imgepsz
        domega_c = zz - his_C(iw) + imgepsz
        domega_l = zz - his_L(iw) + imgepsz
        if( abs(domega_c)<eps .or. abs(domega_r)<eps ) then
          rr_fac(iw) = 0d0
        else
! rr_fac(his_C(is)) = \int^{his_R}_{his_C} d omega' /(his_C(is) -omega')
c            rr_fac(iw) = log( abs((domega_r/domega_c)) )
          rr_fac(iw) = log( domega_r/domega_c )
        endif
        if( abs(domega_c)<eps .or. abs(domega_l)<eps ) then
          rl_fac(iw) = 0d0
        else
! rl_fac(his_C(is)) = \int^{his_C}^{his_L} d omega' /(his_C(is) -omega')
c            rl_fac(iw) = log( abs((domega_c/domega_l)) )
          rl_fac(iw) = log( domega_c/domega_l)
        endif
      enddo
      rmat=0d0
      do iw = -nwhis, nwhis !symmetric version. iw=0 is meaningless
        if(iw==0) cycle
c          if(debug) print *,' it iw=',it, iw
        domega_c = zz - his_C(iw)
        if(iw==  nwhis) then
          delta_r = his_R(iw)   - his_C(iw)
        elseif(iw== -1) then
          delta_r =   0d0       - his_C(iw)
        else
          delta_r = his_C(iw+1) - his_C(iw)
        endif
!         if(debug) print *,' it iw RRR1'
        if(iw== -nwhis) then
          delta_l = his_C(iw)  - his_L(iw)
        elseif(iw==  1) then
          delta_l = his_C(iw)  - 0d0
        else
          delta_l = his_C(iw)  - his_C(iw-1)
        endif
!         if(debug) print *,' it iw RRR2'
!          ddr = (his_R(iw)-his_C(iw))/delta_r
!          ddl = (his_C(iw)-his_L(iw))/delta_l
        rmat(iw)  = rmat(iw  ) + rr_fac(iw)*( 1d0-domega_C/delta_r) !+ ddr
        if(iw/=nwhis.and.iw/=-1) then
          rmat(iw+1) = rmat(iw+1) + rr_fac(iw)*domega_C/delta_r     !- ddr
        endif
        rmat(iw)  = rmat(iw) + rl_fac(iw)*( 1d0+domega_C/delta_l)   !- ddl
        if(iw/=-nwhis.and. iw/=1) then
          rmat(iw-1) = rmat(iw-1) - rl_fac(iw)*domega_C/delta_l     !+ ddl
        endif
ccccccccccccccccccccccccc
c no-derivarive test
c          rmat(iw)  =  rr_fac(iw) + rl_fac(iw)
ccccccccccccccccccccccccc
      enddo
      end

c$$$      subroutine reducezmel(aold, ngbo,ngb,nx,
c$$$     i       io,   in,   nmat, pmat,
c$$$     i     io_q, in_q, nmat_q, pmat_q,
c$$$     o      anew)
c$$$c   For given q+G basis, we augment the basis within MT.
c$$$c   For given atom and l  prod and prodd at MT boundary (reserved in PPBRD
c$$$      integer(4):: nmat,io(nmat),in(nmat),nmat_q,io_q(nmat),in_q(nmat)
c$$$      complex(8):: aold(ngbo,nx), anew(ngb,nx),pmat(nmat) ,pmat_q(nmat)
c$$$      anew=0d0
c$$$      do ix=1,nmat
c$$$          anew(in(ix), :)
c$$$     &  = anew(in(ix), :)   +  pmat(ix) * aold(io(ix), :)
c$$$      enddo
c$$$      do ix=1,nmat_q
c$$$          anew(in_q(ix), :)
c$$$     &  = anew(in_q(ix), :)  + dconjg(pmat_q(ix)) * aold(io_q(ix), :)
c$$$      enddo
c$$$      end

      real(8) function wcutef(e,ecut,ecuts)
      real(8):: e,ecut,ecuts
c      wcutef = 1d0/( exp((e-ecut)/ecuts)+ 1d0)
      wcutef = exp( -(e/ecut)**2 ) ! ecuts is not used in this case
      end
