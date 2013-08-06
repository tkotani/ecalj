      subroutine x0kf_v4hz(
     .     npm,ncc,
     .     ihw,nhw,jhw,
     .     whw,
     .     nhwtot,
     .     n1b,n2b,
     .     nbnbx,
     .     nbnb,
     .     q,
     .     nsp,isp1,isp2,
     .     symmetrize,
     .     qbas,ginv,rk,wk,
     .     mdim,nlnm,nlnmv,nlnmc,iclass,
     .     ppb,
     .     icore,ncore,
     .     nlmto,nqbz,noccx,noccxv,nl,nclass,natom,nnc,nlnmx,mdimx,
     .     nbloch,nwt,iq,ngbb,ngc,ngpmx,ngcmx,nqbze,nband,nqibz, 
     .     rcxq,
     .     nolfco,
     .     zzr,
     .     nmbas,
     .     ppovlz,
     .     ppovl,zcousq,
     .     chipmzzr,eibzmode,
     .     nwgt,igx,igxt,
     .     ngrp,eibzsym)

C     !---- modules
      use m_mpi
      use m_readqg,only :  readqg
      use m_readeigen,only: readeval,readcphi,readgeig
      use keyvalue,only:  getkeyvalue
      use m_rotMPB,only:  rotMPB2
      use m_readqgcou,only: qtt_,nqnum
      implicit none

C     !---- argments
      integer :: npm,ncc
      integer :: ihw(nbnbx,nqbz,npm),nhw(nbnbx,nqbz,npm)
      integer :: jhw(nbnbx,nqbz,npm)
      real(8) :: whw(nhwtot)
      integer :: nhwtot
      integer :: n1b(nbnbx,nqbz,npm),n2b(nbnbx,nqbz,npm)
      integer :: nbnbx
      integer :: nbnb(nqbz,npm)
      real(8) :: q(3)
      integer :: nsp,isp1,isp2
      logical :: symmetrize
      real(8) :: qbas(3,3),ginv(3,3),rk(3,nqbz),wk(nqbz)
      integer :: mdim(natom),nlnm(*),nlnmv(*),nlnmc(*),iclass(*)
      real(8) :: ppb(*)
      integer :: icore(*),ncore(*)
      integer :: nlmto
      integer :: nqbz,noccx,noccxv,nl,nclass,natom,nnc,nlnmx,mdimx
      integer :: nbloch,nwt,iq,ngbb,ngc,ngpmx,ngcmx,nqbze,nband,nqibz
      complex(8) :: rcxq(nmbas,nmbas,nwt,npm)
      logical :: nolfco
      complex(8),target :: zzr(ngbb,nmbas)
      integer :: nmbas
      complex(8),target :: ppovlz(ngbb,ngbb)
      complex(8) :: ppovl(ngc,ngc),zcousq(ngbb,ngbb) 
      logical :: chipmzzr,eibzmode
      integer :: nwgt(nqbz),igx(ngrp*2,nqbz),igxt(ngrp*2,nqbz)
      integer :: ngrp
      integer :: eibzsym(ngrp,-1:1)


C     !---- functions
      integer :: iopen,iclose
      logical :: checkbelong
      logical :: smbasis
      integer :: nocc

C     !---- local constants
      real(8), parameter :: symope(3,3) =
     .     reshape( (/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/),(/3,3/) )
      real(8), parameter :: shtv(3)=(/0d0,0d0,0d0/)
      complex(8), parameter :: img=(0d0,1d0)
      logical, parameter:: debug=.false.


C     !---- local variables
      integer :: iatom,nctot,nbmx,jpm,ibib,itps,nt0,ntp0
      integer :: ngp1,ngp2,it,itp,iw,igb2,igb1,ngb,nn,no,isx,k
      integer :: nadd(3)
      integer :: ngvecpB1(3,ngpmx),ngvecpB2(3,ngpmx)
      real(8) :: ebmx
      real(8) :: qbasinv(3,3),det,qdiff(3),add(3)
      complex(8) :: phase(natom)
      complex(8), allocatable :: cphik(:,:),cphikq(:,:)
      complex(8), allocatable :: geig1(:,:),geig2(:,:)
      integer :: imdim(natom),iatomp(natom)
      logical :: ipr
      real(8) :: imagweight
      real(8) :: eband(nband)
      real(8) :: qu1(3),qu2(3),qu1x(3),qu2x(3),quu(3),tolqu=1d-4
      integer :: nbcut,nbcut2
      logical :: iww1=.true.,iww2=.true.
      integer :: ifpomat,iqxdummy
      complex(8),allocatable :: pomat(:,:)
      real(8) :: q_r(3)
      integer :: nkmin,nkmax,nkqmin,nkqmax
      real(8) :: qq(3)
      integer :: ix,iy
      logical :: onceww
      complex(8),allocatable:: rcxq0(:,:),rcxq00(:,:),rcxq000(:,:),rcxqwww(:,:)
      complex(8) :: zcousqc(ngbb,ngbb),zcousqinv(ngbb,ngbb)
      integer :: neibz,icc,ig,ikp,i,j,itimer,icount,iele
      integer :: irotm,nrotmx,itt,ntimer,nccc,iagain,irotm1,irotm2
      integer,allocatable:: i1(:,:),i2(:,:),nrotm(:)
      complex(8),allocatable:: zrotm(:,:),zrr(:,:),zrrc(:,:),zrrx(:,:)

C     !---- mpi parallelized arrays
      complex(8),allocatable :: zmat_mpiM(:,:)
      complex(8),allocatable :: zpsi_mpiM(:)
      complex(8),allocatable :: zmelt_mpiM(:)
      complex(8),allocatable :: z1p_mpiM(:,:)
      complex(8),allocatable :: rcxq_mpiMP(:,:,:,:)
      complex(8),allocatable :: rcxq_mpiMPW(:,:,:,:)
      complex(8),allocatable :: rcxq_mpiP(:,:,:,:)

c$$$      real(8) :: time_start, time_end
    

!TIME0 "start xkf_v4hz      
      write(*,'(" x0kf_v4hz: q=",3f8.4,$)')q
      call cputid(0)
      imdim(1) = 1
      do iatom = 1,natom
         iatomp(iatom) = iatom
         if(iatom<natom) imdim(iatom+1)=imdim(iatom)+mdim(iatom)
      end do
      nctot      = noccx - noccxv
      call dinv33(qbas,0,qbasinv,det)
      phase= (1d0,0d0)          !coskt = 1d0; sinkt = 0d0
      allocate( cphik(nlmto,nband),cphikq(nlmto,nband) )
      allocate( geig1(ngpmx,nband),  geig2(ngpmx,nband) )
      call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
      call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
      call getnemx(nbmx,ebmx,7,.true.)


!TIME1 "before if smbasis"
!TIME0 
      if(smbasis()) then        !need to check again, when we will make smbasis on.
         if(ncc/=0) then
            print *,"Timereversal=F(ncc/=0) not yet implemented for smbasis."
            print *," pomat should be generated correctly                  ."
            stop !"Timereversal=F(ncc/=0) not yet implemented for smbasis."
         endif

         ifpomat = iopen('POmat',0,-1,0) !oct2005

         do
            read(ifpomat) q_r,nn,no,iqxdummy !readin reduction matrix pomat
            allocate( pomat(nn,no) )
            read(ifpomat) pomat
            if( sum(abs(q-q_r))<1d-10) then ! .and.kx <= nqibz ) then
               write(6,*) 'ok find the section for give qibz_k'
               exit
            endif
c$$$  deallocate(pomat) ! MIZUHO-IR
         end do
         if( sum(abs(q-q_r))>1d-10 ) then
            write(6,"('q  =',3f12.5)") q
            write(6,"('q_r=',3f12.5)") q_r
            stop !'POmat reading err q/=q_r'
         endif
         isx = iclose('POmat')
      endif
!TIME1 "after if smbasis"

      ngb=nbloch+ngc
      call MPI__x0kf_rankdivider(nbnb,nqbz,npm,ngb,nwt)

      allocate( zmat_mpiM(ngb,mpi__Ms:mpi__Mf) )
      allocate( zpsi_mpiM(mpi__Ms:mpi__Mf) )
      allocate( zmelt_mpiM(mpi__Ms:mpi__Mf) )
      allocate( z1p_mpiM(ngb,mpi__Ms:mpi__Mf) )
      allocate( rcxq_mpiMP(ngb,mpi__Ms:mpi__Mf,mpi__Wnall,mpi__Ps:mpi__Pe) )
      allocate( rcxq_mpiMPW(ngb,mpi__Ms:mpi__Mf,mpi__Ws:mpi__We,mpi__Ps:mpi__Pe) )
      allocate( rcxq_mpiP(ngb,ngb,nwt,mpi__Ps:mpi__Pe) )

!$OMP PARALLEL DO PRIVATE(igb1,igb2,jpm)
      do igb1=1,ngb
         do igb2=mpi__Ms,mpi__Me
            do jpm=mpi__Ps, mpi__Pe
               rcxq_mpiMP(igb1,igb2,:,jpm) = rcxq(igb1,igb2,:,jpm)
            end do
         end do
      end do
!$OMP END PARALLEL

      if(chipmzzr) then
!$OMP PARALLEL DO PRIVATE(igb1,igb2)
         do igb1=1,ngb
            do igb2=mpi__Ms,mpi__Me
               zmat_mpiM(igb1,igb2) =  zzr(igb2,igb1)
            end do
         end do
!$OMP END PARALLEL
      elseif(nolfco .and. nmbas==1) then
!$OMP PARALLEL DO PRIVATE(igb1,igb2)
         do igb1=1,ngb
            do igb2=mpi__Ms,mpi__Me
               zmat_mpiM(igb1,igb2) = ppovlz(igb2,igb1)
            end do
         end do
!$OMP END PARALLEL
      else
!$OMP PARALLEL DO PRIVATE(igb2)
         do igb1=1,ngb
            do igb2=mpi__Ms,mpi__Me
               if( igb1 == igb2 .and. igb2<=nbloch ) then
                  zmat_mpiM(igb1,igb2) = 1.0d0
               else
                  zmat_mpiM(igb1,igb2) = 0.0d0
               end if
            end do
         end do
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(igb1,igb2)
         do igb1=1,ngb
            if( igb1<=nbloch ) cycle
            do igb2=mpi__Ms,mpi__Me
               if( igb2<=nbloch ) cycle
               zmat_mpiM(igb1,igb2) = ppovl(igb2-nbloch,igb1-nbloch)
            end do
         end do
!$OMP END PARALLEL
      endif

c$$$      if( mpi__rootP .and. mpi__rootB .and. mpi__rootM ) then
c$$$         time_start = mpi_wtime()
c$$$      end if

      
      qq=0d0
C     ! parallelized B&K-loop
      do k=1, nqbz
         if( .not. mpi__Btask1(k) ) cycle ! skip, not for this process

         if(eibzmode.and.nwgt(k)==0 ) cycle
         if(debug) write(*,'("do 1000  k=",i4)')k
         ipr=(k<5.or.k==nqbz.or.debug)
         if(sum(nbnb(k,1:npm))==0) cycle
!TIME0 "start do 1000"
         if(k<=5) then
            write(6,"('  x0kf_v4hz: k rk=',i7,3f10.4,$)")k, rk(:,k)
            call cputid(0)
         elseif(mod(k,max(10,nqbz/20))==1.or.k>nqbz-10) then
            write(6,"('  x0kf_v4hz: k rk=',i7,3f10.4,$)")k, rk(:,k)
            call cputid(0)
         endif

         call  readcphi(  rk(:,k)-qq, nlmto,isp1,   quu,cphik) !quu is used q for eigenfunctions.
         call  readcphi(q+rk(:,k)-qq, nlmto,isp2,   quu,cphikq)

         nkmin = 999999
         nkmax= -999999
         nkqmin= 999999
         nkqmax=-999999
         do jpm=1,npm
            do ibib = 1, nbnb(k,jpm)
               nkmin  = min(n1b(ibib,k,jpm),nkmin)
               nkqmin = min(n2b(ibib,k,jpm),nkqmin)
               if(n1b(ibib,k,jpm)<=nband)   nkmax = max(n1b(ibib,k,jpm),nkmax)
               if(n2b(ibib,k,jpm)<=nband)  nkqmax = max(n2b(ibib,k,jpm),nkqmax)
            end do
         end do
         call readeval(q+rk(:,k)-qq,isp2,eband)
         nkqmax = nocc (eband,ebmx,.true.,nband)
         if(npm==2) then
            call readeval(rk(:,k)-qq,isp1,eband)
            nkmax = nocc (eband,ebmx,.true.,nband)
         endif
         itps  = nkqmin
         nt0   = nkmax
         ntp0  = nkqmax - nkqmin +1
         if( npm==2.and. nkqmin/=1) then
            write(6,*)' npm==2 nkqmin nkqmax  nkmin nkmax=',nkqmin,nkqmax,nkmin,nkmax
            stop !" When npm==2, nkqmin==1 should be."
         endif
         if(nkmin/=1) then
            stop !" nkmin==1 should be."
         endif

         if(onceww(5) ) print *,' nctot ncc=',nctot,ncc




         call readqg('QGpsi',q+rk(:,k)-qq,ginv, qu1x, ngp1, ngvecpB1)
         call readqg('QGpsi',  rk(:,k)-qq,ginv, qu2x, ngp2, ngvecpB2)
!     !      ngp1 = ngpn(kp)  ! q+k   ntp0 in FBZ
!     !      ngp2 = ngpn(k)   ! k     np0  in FBZ
!     !      ngc              ! q          in IBZ
         ngb  = nbloch + ngc    ! This is not ngbb for smbasis()=T. oct2005
         if(ngb/=ngbb) then
            print *,' x0kf_v4h: ngb ngbb=',ngb,ngbb
            stop !'x0kf_v4h: ngb/=ngbb'
         endif

!     !                        q       k        q+k

         call readgeig(q+rk(:,k)-qq, ngpmx,isp2, qu1, geig1)
         call readgeig(  rk(:,k)-qq, ngpmx,isp1, qu2, geig2)

         if(sum(abs(qu1-qu1x))>tolqu) then
            write(6,"('qu1 :',3d23.10)") qu1
            write(6,"('qu1x:',3d23.10)") qu1x
            stop !'x0kf_v4hz:qu1/=qu1x'
         endif
         if(sum(abs(qu2-qu2x))>tolqu) then
            write(6,"('qu2 :',3d23.10)") qu2
            write(6,"('qu2x:',3d23.10)") qu2x
            stop !'x0kf_v4hz:qu2/=qu2x'
         endif
         qdiff = q    - qu1    + qu2
                                ! q   - (q+k)   + k  is not zero.
                                ! qc  -  q1     + q2
         add   = matmul(qbasinv, qdiff)
         nadd  = idint( add + dsign(.5d0,add)) !  print *,' qdif=',qdiff,qbkp(:),rk(:,k)
         if(sum(abs(add-nadd))>1d-10) stop !"sexc: abs(add-nadd))>1d-10"

!TIME1 "before psicb_v3"
!TIME0 
C     ! initialize data in advace to calculate zpsi
         call element_psicbA(
     .        nctot,ncc,nt0,ntp0,iclass,phase,
     .        cphik(1,nkmin),cphikq(1,nkqmin), 
     .        ppb,nlnmv,nlnmc,mdim,imdim,iatomp,
     .        mdimx,nlmto,nbloch,ngc,nlnmx,natom,nclass,
     .        icore,ncore,nl,nnc,
     .        zpsi_mpiM(1),0,0,0,.true.)

C     ! initialize data in advace to calculate zpsi
         call element_psicbB(
     .        nctot,ncc,nt0,ntp0,iclass,phase,
     .        cphik(1,nkmin),cphikq(1,nkqmin), 
     .        ppb,nlnmv,nlnmc,mdim,imdim,iatomp,
     .        mdimx,nlmto,nbloch,ngc,nlnmx,natom,nclass,
     .        icore,ncore,nl,nnc,
     .        zpsi_mpiM(1),0,0,0,.true.)

C     ! initialize data in advace to calculate zpsi
         call element_psi2b(
     .        nctot,ncc,nt0,ntp0,iclass,phase,
     .        cphik(1,nkmin),cphikq(1,nkqmin),
     .        ppb,nlnmv, nlnmc,mdim,imdim,iatomp,
     .        mdimx,nlmto,nbloch,ngc,nlnmx,natom,nclass,
     .        zpsi_mpiM(1),0,0,0,.true.)
!TIME1 "after psicb_v3"
!TIME0 

!TIME1 "before melpln2"
!TIME0  
C     ! initialize data in advace to calculate zpsi
         if(ngc/=0) then
            call element_melpln2(
     .           nbloch,ngb,nctot,ncc,ngp1,ngvecpB1,ngp2,ngvecpB2,
     .           ngc,nadd,geig1,geig2,ngpmx,nband,itps,ntp0,nt0,
     .           shtv,q,q,symope,qbas,zpsi_mpiM(1),0,0,0,.true.)
         endif

!TIME1 "after melpln2"
!TIME0 

         if(smbasis()) then
            stop !' not maintained recently need to check again' !takao for safe jun2012
         endif

!TIME1  "after matmul zmelt"
!TIME0 

C     ! parallelized P-loop
         do jpm=mpi__Ps, mpi__Pe
            if( .not. mpi__Btask2(k,jpm) ) cycle ! skip, not for this process

C     ! parallelized B&K-loop
            do ibib = 1, nbnb(k,jpm)
               if( .not. mpi__Btask3(k,jpm,ibib) ) cycle ! skip, not for this process

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

               if(jpm==1) then
                  if(n2b(ibib,k,jpm)>nbmx)  then !nbmx
                     if(iww1) then
                        print *,' nband_chi0 nbmx=',nbmx
                        iww1=.false.
                     endif
                     cycle
                  endif
                  if( n1b(ibib,k,jpm) <= nbcut .and.
     .                 nbcut2<n2b(ibib,k,jpm) ) then
                     if(iww2) then
                        write(6,"(' nband_chi0 nbcut nbcut2 n2b n1b=',4i6)")
     .                       nbcut,n2b(ibib,k,jpm),n1b(ibib,k,jpm)
                        iww2=.false.
                     endif
                     cycle
                  endif

               else
                  if( n1b(ibib,k,jpm) > nbmx) then !nbmx
                     if(iww1) then
                        print *,' nband_chi0 nbmx=',nbmx
                        iww1=.false.
                     endif
                     cycle
                  endif
                  if( n2b(ibib,k,jpm) <= nbcut .and.
     .                 nbcut2<n1b(ibib,k,jpm) ) then
                     if(iww2) then
                        write(6,"(' nband_chi0 nbcut nbcut2 n2b n1b=',4i6)")
     .                       nbcut,n2b(ibib,k,jpm),n1b(ibib,k,jpm)
                        iww2=.false.
                     endif
                     cycle
                  endif
               endif

               if (ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1 >nwt) stop !"x0kf_v4hz: iw>nwt"


C     ! parallelized M-loop
               do igb2=mpi__Ms,mpi__Me

                  if( it<=nctot .and. itp> ncc ) then
                     if( igb2<=nbloch ) then
C     ! zpsi at (igb2,it,itp) is obtained by element_psicbA
                        call element_psicbA(
     .                       nctot,ncc,nt0,ntp0,iclass,phase,
     .                       cphik(1,nkmin),cphikq(1,nkqmin), 
     .                       ppb,nlnmv,nlnmc,mdim,imdim,iatomp,
     .                       mdimx,nlmto,nbloch,ngc,nlnmx,natom,nclass,
     .                       icore,ncore,nl,nnc,
     .                       zpsi_mpiM(igb2),igb2,it,itp,.false.)
                     else
C     ! zpsi at (igb2,it,itp) is zero
                        zpsi_mpiM(igb2) = 0.0d0
                     end if
                  else if( it> nctot .and. itp<=ncc ) then
                     if( igb2<=nbloch ) then
C     ! zpsi at (igb2,it,itp) is obtained by element_psicbB
                        call element_psicbB(
     .                       nctot,ncc,nt0,ntp0,iclass,phase,
     .                       cphik(1,nkmin),cphikq(1,nkqmin), 
     .                       ppb, nlnmv,nlnmc,mdim, imdim,iatomp,
     .                       mdimx,nlmto,nbloch,ngc,nlnmx,natom,nclass,
     .                       icore,ncore,nl,nnc,
     .                       zpsi_mpiM(igb2),igb2,it,itp,.false.)
                     else
C     ! zpsi at (igb2,it,itp) is zero
                        zpsi_mpiM(igb2) = 0.0d0
                     end if
                  else if( it> nctot .and. itp> ncc ) then
                     if( igb2<=nbloch ) then
C     ! zpsi at (igb2,it,itp) is obtained by element_psi2b
                        call element_psi2b(
     .                       nctot,ncc,nt0,ntp0,iclass,phase,
     .                       cphik(1,nkmin), cphikq(1,nkqmin),
     .                       ppb, nlnmv, nlnmc,mdim, imdim,iatomp,
     .                       mdimx,nlmto,nbloch,ngc,nlnmx,natom,nclass,
     .                       zpsi_mpiM(igb2),igb2,it,itp,.false.)
                     else
C     ! zpsi at (igb2,it,itp) is obtained by element_melpln
                        call element_melpln2(
     .                       nbloch,ngb,nctot,ncc,ngp1,ngvecpB1,
     .                       ngp2,ngvecpB2,ngc,nadd,
     .                       geig1,geig2,ngpmx,nband,itps,ntp0,nt0,
     .                       shtv,q,q,symope,qbas,
     .                       zpsi_mpiM(igb2),igb2,it,itp,.false.)
                     end if
                  else
C     ! zpsi at (igb2,it,itp) is zero
C     ! this code is never executed
                     zpsi_mpiM(igb2) = 0.0d0
                  end if
               end do ! igb2

               call MPI__ZGEMV(zmat_mpiM,zpsi_mpiM,zmelt_mpiM)
c$$$               do igb1=1,nmbas
c$$$                  zmelt(igb1) = 0.0d0
c$$$                  do igb2=1,nbloch+ngc
c$$$                     zmelt(igb1) = zmelt(igb1) + zmat(igb1,igb2)*zpsi(igb2)
c$$$                  end do
c$$$               end do

               call MPI__ZGERC(zmelt_mpiM,zmelt_mpiM,z1p_mpiM)
c$$$               do igb2=1, nmbas
c$$$                  do igb1=1, nmbas
c$$$                     z1p(igb1,igb2) = dconjg(zmelt(igb1)) * zmelt(igb2)
c$$$                  end do
c$$$               end do

C     ! not parallelized W-loop
               do iw = ihw(ibib,k,jpm),ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1
                  imagweight = whw(jhw(ibib,k,jpm)+iw-ihw(ibib,k,jpm))
                  if(eibzmode) imagweight = nwgt(k)*imagweight  

                  call MPI__ZAMPM(imagweight,z1p_mpiM,rcxq_mpiMP(:,:,iw,jpm))
c$$$                  do igb2=1,nmbas
c$$$                     do igb1=1, nmbas
c$$$                        rcxq(igb1,igb2,iw,jpm) =
c$$$     &                       rcxq(igb1,igb2,iw,jpm) + z1p(igb1,igb2)*imagweight
c$$$                     end do ! igb1
c$$$                  end do ! igb2
               end do ! iw
            end do ! jbib
         end do ! jpm

!TIME1 "after rcxq"

         if(debug) call cputid(0)
         if(debug) write(6,*)' end of kloop k jpm=',k,jpm
      end do ! k

c$$$      if( mpi__rootP .and. mpi__rootB .and. mpi__rootM ) then
c$$$         time_end = mpi_wtime()
c$$$         write(*,*) "DEBUG: timeA=", mpi__rankQ, time_end - time_start
c$$$      end if

c$$$      if( mpi__rootP .and. mpi__rootB .and. mpi__rootM ) then
c$$$         time_start = mpi_wtime()
c$$$      end if


C     ! parallelized P-loop
      do jpm=mpi__Ps, mpi__Pe
C     ! not parallelized W-loop
        do iw=1, nwt
C     ! reduce B&K parallel
           call MPI__AllreduceMatrixB(rcxq_mpiMP(:,:,iw,jpm))
C     ! gather M parallel
           call MPI__AllgatherMatrixM(rcxq_mpiMP(:,:,iw,jpm),rcxq_mpiP(:,:,iw,jpm))
        end do
      end do

C     ! gather P parallel
      call MPI__AllgatherArrayP(rcxq_mpiP(:,:,:,:),rcxq(:,:,:,:))


!$OMP PARALLEL DO PRIVATE(igb1,igb2,iw,jpm)
      do igb1=1,ngb
         do igb2=mpi__Ms,mpi__Me
            do iw=mpi__Ws, mpi__We
               do jpm=mpi__Ps, mpi__Pe
                  rcxq_mpiMPW(igb1,igb2,iw,jpm) = rcxq_mpiMP(igb1,igb2,iw,jpm)
               end do
            end do
         end do
      end do
!$OMP END PARALLEL

c$$$      if( mpi__rootP .and. mpi__rootB .and. mpi__rootM ) then
c$$$         time_end = mpi_wtime()
c$$$         write(*,*) "DEBUG: timeB=", mpi__rankQ, time_end - time_start
c$$$      end if

      deallocate( zmat_mpiM, zpsi_mpiM, zmelt_mpiM, z1p_mpiM )
      deallocate( rcxq_mpiMP, rcxq_mpiP, rcxq_mpiMPW )

!     ! Not need to be symmetrized
      if(nolfco .and. nmbas==1) then
         print *,' nmbas=1 nolfco=T ---> not need to symmetrize'
         goto 9999
      endif
!TIME0

      write(*,*) "=== DEBUG: rcxq ", sum(rcxq)


!TIME1 "before eibzmode"


!     ! ------------------------------------------------------------------------
!     ! == Symmetrizer of EIBZ PRB.81,125102(2010) Eq.(51) july2012takao ==
!     ! --- zrotm(J,J') = <Mbar^k_J| \hat{A}^k_i Mbar^k_J'>. ---
!     ! We do \sum_i T_alpha_i [ zrotm_i^dagger (I,I') P_I'J' zrom_i(J'J) ] 
!     ! (exactrly speaking, we insert conversion matrix between Enu basis and M_I basis).
!     !
!     ! input qin = q
!     ! \hat{A}^k_i  is specified by symops(:,:,igx),and igxt (-1 for time-reversal).
!     ! Note that k= \hat{A}^k_i(k) (S_A^k)
!     ! See Eq.(51) around in PRB81 125102(2010)
!     ! 

!     ! === zmelt conversion === 
      if(eibzmode.and.symmetrize) then
!TIME0 
         call iqindx2(q, ginv, qtt_, nqnum, ikp,quu) !to get ikp for timereversal mode
!TIME1 "after iqindx2"
!TIME0
         if(sum(abs(q-quu))>1d-8) stop !'x0kf_v2h: eibz 111 q/quu'
         neibz = sum(eibzsym(:,1))+sum(eibzsym(:,-1)) 
                                !itimer=-1 means time reversal. eibzsym(ig,itimer) where ig: space rotation.
         write(6,"(' --- goto symmetrization --- ikp neibz q=',2i3,3f12.8)")ikp,neibz,q
         call cputid2(' --- qqqqq: start symmetrization  ',0)

         ntimer=1
         if(sum(eibzsym(:,-1))>0) ntimer=2 !timereversal case

         allocate(zrotm(ngb,ngb),nrotm(ngrp*2))
         zcousqinv=zcousq
         call matcinv(ngb,zcousqinv)

!     ! == Assemble rotantion matrx zrr,zrrc ==
!     ! Rotation matrix zrrx can be a sparse matrix. 
!     ! Thus it is stored to  "i1(nrotmx,nccc),i2(nrotmx,nccc),zrr(nrotmx,icc),nrotm(icc)".
!     ! See folloings: matmul(rcxqwww,zrrx) is given by
!     !     do irotm1 = 1,nrotm(icc)
!     !       rcxq0(:,i2(irotm1,icc)) = rcxqwww(:,i1(irotm1,icc)) * zrr(irotm1,i2(irotm1,icc))

         allocate(zrrx(nmbas,nmbas))
         nrotmx = 10000         !trial value
!TIME1 "before 1011"
         do 1011                !this loop is only in order to to set large enough nrotmx.
!TIME0 
            if(allocated(i1)) deallocate(i1,i2,zrr,zrrc) !,zrr_,zrrc_)
            nccc=ngrp*2
            allocate(i1(nrotmx,nccc),i2(nrotmx,nccc),zrr(nrotmx,nccc),zrrc(nrotmx,nccc)) !,zrr_(ngb,ngb,nccc),zrrc_(ngb,ngb,nccc))
            i1=-99999
            i2=-99999
            zrr=-99999d0
            zrrc=-99999d0
            call cputid2(' --- qqqqq111   :',0)
!TIME1 "call cputid2 qqqqq11"
            icc=0
            do itimer=1,-1,-2
               if(ntimer==1.and.itimer==-1) exit
               if(itimer==1 ) itt=1
               if(itimer==-1) itt=2
               do ig=1,ngrp
                  if(eibzsym(ig,itimer)==1) then
                     icc=icc+1
!TIME0

!     ! Get rotation matrix zrrx, which can be a sparse matrix. Thus stored to zrr.
                     call rotMPB2(nbloch,ngb,q,ig,itimer,ginv,zrotm) 
                     if(nolfco.and.chipmzzr) then 
!     !   We assume <svec_I | svec_J >= \delta_IJ, In addition, we use fact that we have no IPW parts in svec.
!     !   If IPW part exist, we may have to take into account <IPW|IPW> matrix, e.g. as in ppovlz. 
!     !   svec --> zzr
                        if(itimer==1) then
                           zrrx= matmul(transpose(dconjg(zzr)), matmul(zrotm, zzr)) 
                        else  
                           zrrx= matmul(transpose(zzr), matmul(dconjg(zrotm), zzr)) 
                        endif 
                     elseif(nolfco) then 
                        stop !'x0kf_v4h: this case is not implemented xxxxxxxxxxxxxx'
                     else   
!     ! zrotm(J,J') is the rotation matrix = <Mbar^k_J| \hat{A}^k_i Mbar^k_J'>
!     ! See rotMPB2 defined in readeigen.F.
!     ! zrrx(mu nu)= dconjg(Zcousq(I, mu)) *zrotm(I,J)* Zcousq(J, nu)
!     ! zrrx is very sparse matrix. Size is \sim ngb or something.

                        if(itimer==1) then
                           zrrx=zrotm
                        else  
                           zrrx=dconjg(zrotm)
                        endif 
                     endif
!TIME1 "end matmmsparse"
!TIME0
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
                        end do
                     end do
!TIME1 "before iagain=1"
!TIME0
                     if(iagain==1) then
                        nrotmx=irotm !enlarge allocation and do things again.
                        print *," warn:(slow speed) xxxx goto 1011 ",
     .                       "xxxxxx nrotmx+=nrotmx+10000 again"
                        goto 1011 !enlarge nrotmx ang try it again.
                     endif   
                     nrotm(icc)=irotm
                     print *,'ig itimer icc nrotm=',ig,itimer,icc,nrotm(icc) ,iele
!TIME1 "end ig itimer icc, nrotm"
                  endif
               end do
            end do
!TIME0
            exit
 1011    continue               !only when nrotmx overflow.


!     ! === main part to obtain symmetrized rcxq  ===
!     ! neibz is total number of symmetrization operation. 
!     !      rcxq is rotated and accumulated; finally divied by neibz
         zcousqc=dconjg(transpose(zcousq))
         call cputid2(' --- qqqqq222ini:',0)
!$OMP parallel private(rcxq000,icc,itt,icount,rcxqwww,rcxq00,rcxq0)
         allocate(rcxq0(ngb,ngb),rcxq00(ngb,ngb),rcxq000(ngb,ngb),rcxqwww(ngb,ngb))
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

!     !  Followings are equivalent with
!     !            rcxq00= rcxq00 + matmul(zrrc_(:,:,icc),matmul(rcxqwww,zrr_(:,:,icc)))
                        do irotm1 = 1,nrotm(icc)
                           rcxq0(:,i2(irotm1,icc)) =rcxq0(:,i2(irotm1,icc)) + rcxqwww(:,i1(irotm1,icc)) * zrr(irotm1,icc)
                        end do
                        do irotm2 = 1,nrotm(icc)
                           rcxq00(i2(irotm2,icc),:)= rcxq00(i2(irotm2,icc),:) + zrrc(irotm2,icc) * rcxq0(i1(irotm2,icc),:)
                        end do

                     endif
                  end do

                  if(itimer==1) then
                     rcxq000=rcxq00
                  else
                     rcxq000=rcxq000+rcxq00
                  endif   
               end do
               rcxq(:,:,iw,jpm) = rcxq000/neibz
ccccccccccccccccccc
               rcxq(:,:,iw,jpm) = matmul(rcxq(:,:,iw,jpm),zcousq)
               rcxq(:,:,iw,jpm) = matmul(zcousqc,rcxq(:,:,iw,jpm))
ccccccccccccccccccc
            end do
         end do
!$OMP end  do 
         deallocate(rcxq00,rcxq000,rcxq0,rcxqwww)
!$OMP end parallel 
!TIME1 "after sym rcxq"
         deallocate(zrotm,i1,i2)

         call cputid2(' --- qqqqq222end:',0)
      endif  

 9999 continue
      deallocate(cphik,cphikq,geig1,geig2)
      if(smbasis()) deallocate(pomat)
      write(6,"(' --- x0kf_v4hz: end')")
      end subroutine x0kf_v4hz




      subroutine element_psi2b(
     .     nctot,ncc,nt0,ntp0,iclass,
     .     phase,cphik,cphikq,
     .     ppb,
     .     nlnmv,nlnmc,mdim,imdim,iatomp,
     .     mdimx,nlmto,nbloch,ngc,nlnmx,natom,nclass,
     .     zpsi,target_ib,target_it,target_itp,initialize)
      implicit none
C     !---- arguments
      integer :: nctot,ncc,nt0,ntp0,iclass(natom)
      complex(8) :: phase(natom), cphik(nlmto,*),cphikq(nlmto,ntp0)
      real(8) :: ppb(nlnmx,nlnmx,mdimx,nclass)
      integer :: nlnmv(nclass),nlnmc(nclass),mdim(nclass),imdim(natom),iatomp(natom)
      integer :: mdimx,nlmto,nbloch,ngc,nlnmx,natom,nclass
      complex(8) :: zpsi
      integer :: target_ib, target_it, target_itp
      logical :: initialize

C     !---- local variables
      integer:: ia,ic,nc,iv1,iv2,nv,nc1,ias,iap,icp,i,it,itp,jp,ib
      complex(8) :: zz, zwork
C     !---- saved variables
      integer,allocatable,save :: iasx(:)

      if( initialize ) then
         if(mdimx /= maxval(mdim) ) then
            stop !'psi2b_v3: wrong mdimx'
         end if
         if( sum(mdim(iclass(1:natom)))/= nbloch ) then
            stop !'psi2b_v3: wrong nbloch'
         end if

         if( allocated(iasx) ) deallocate(iasx)
         allocate(iasx(natom))
         ias = 1
         do ia = 1,natom
            iasx(ia) = ias
            ias = ias + nlnmv(iclass(ia))
         end do
         if(ias-1/=nlmto) stop !' psi2b_v3:sum(nlnmv)/= nlmto'

         return
      end if

      itp = target_itp-ncc
c$$$  do itp=1, ntp0
c$$$  if( ncc+itp /= target_itp ) cycle

      it = target_it-nctot
c$$$  do it=1, nt0
c$$$  if( nctot+it /= target_it ) cycle

      do ia=1,natom
         ic = iclass(ia)
         nc = nlnmc(ic)
         nv = nlnmv(ic)
         nc1 = nc + 1
         if (nc+ nlnmv(ic) > nlnmx) stop !'psi2b_v3: nlnmx exceeded'
         
         ias = iasx(ia)
         iap = iatomp(ia)
         icp = iclass(iap)
         
         do i=1,mdim(icp)
            ib=imdim(iap)-1+i
            if( ib /= target_ib ) cycle ! skip, not for this process

            zwork = 0.0d0
            do iv1=1, nv
               zz = 0.0d0
!$OMP PARALLEL DO PRIVATE(iv2) REDUCTION(+:zz)
               do iv2=1, nv
                  zz = zz + ppb(nc+iv2,nc+iv1,i,icp)
     .                 *cphikq(ias-1+iv2,itp)
               end do
!$OMP END PARALLEL
               zwork = zwork + cphik(ias-1+iv1,it) * dconjg(zz)
            end do
            zpsi = zwork * phase(ia) 
         end do
      end do
c$$$  end do
c$$$  end do

      end subroutine element_psi2b


      subroutine element_psicbA(
     .     nctot,ncc,nt0,ntp0,iclass,
     i     phase,cphik,cphikq, 
     i     ppb,
     i     nlnmv,nlnmc,mdim,imdim,iatomp,
     i     mdimx,nlmto,nbloch,ngc,nlnmx,natom,nclass,
     i     icore,ncore,nl,nnc,  
     o     zpsi,target_ib,target_it,target_itp,initialize)
      implicit none
C     !---- arguments

      integer :: nctot,ncc,nt0,ntp0,iclass(natom)
      complex(8) :: phase(natom),cphik(nlmto,*),cphikq(nlmto,*)
      real(8) :: ppb(nlnmx,nlnmx,mdimx,nclass)
      integer :: nlnmv(nclass),nlnmc(nclass),mdim(nclass)
      integer :: imdim(natom),iatomp(natom)
      integer :: mdimx,nlmto,nbloch,ngc,nlnmx,natom,nclass
      integer :: icore(nl*nl*nnc,nclass),ncore(nclass),nl,nnc
      complex(8) :: zpsi
      integer :: target_ib, target_it, target_itp
      logical :: initialize
C     !---- local variables
      integer :: ib,ias,ics,ia,ic,iv,nc,nv,iap,icp,i,itp,it,icr
      complex(8) :: cphik_ppb

      if( initialize ) then
         if(ncc/=0 .and. ncc/=nctot) then
            stop !"psicb_v3: ncc/=0 and ncc/=ncctot"
         endif

         if(sum(ncore(iclass(1:natom)))  /= nctot) then
            stop !"psicb_v3:sum(ncore) wrong"
         end if

         return
      end if

      itp = target_itp-ncc
c$$$  do itp=1, ntp0
c$$$  if( ncc+itp /= target_itp ) cycle

      ias = 1
      ics = 0
      do ia=1, natom
         ic = iclass(ia)
         nc = nlnmc(ic)
         nv = nlnmv(ic)
         iap = iatomp(ia)
         icp = iclass(iap)
         ib = imdim(iap)-1

         do it=1, ncore(ic)
            if( ics+it /= target_it ) cycle ! skip, not for this process

            do i=1, mdim(icp)
               ib=ib+1
               if( ib /= target_ib ) cycle ! skip, not for this process
               
               icr=icore(it,ic)

               cphik_ppb = 0.0d0
!$OMP PARALLEL DO PRIVATE(iv) REDUCTION(+:cphik_ppb)
               do iv=0, nv-1
                  cphik_ppb = cphik_ppb 
     .                 + cphikq(ias+iv,itp)*ppb(nc+1+iv,icr,i,icp)
               end do
!$OMP END PARALLEL
               zpsi = phase(ia)*dconjg(cphik_ppb)
            end do

            ias = ias + nlnmv(ic)
            ics = ics + ncore(ic)
         end do
      end do
c$$$  end do

      end subroutine element_psicbA


      subroutine element_psicbB(
     .     nctot,ncc,nt0,ntp0,iclass,
     i     phase,cphik,cphikq, 
     i     ppb,
     i     nlnmv,nlnmc,mdim,imdim,iatomp,
     i     mdimx,nlmto,nbloch,ngc,nlnmx,natom,nclass,
     i     icore,ncore,nl,nnc,  
     o     zpsi,target_ib,target_it,target_itp,initialize)
      implicit none
C     !---- arguments

      integer :: nctot,ncc,nt0,ntp0,iclass(natom)
      complex(8) :: phase(natom),cphik(nlmto,*),cphikq(nlmto,*)
      real(8) :: ppb(nlnmx,nlnmx,mdimx,nclass)
      integer :: nlnmv(nclass),nlnmc(nclass),mdim(nclass)
      integer :: imdim(natom),iatomp(natom)
      integer :: mdimx,nlmto,nbloch,ngc,nlnmx,natom,nclass
      integer :: icore(nl*nl*nnc,nclass),ncore(nclass),nl,nnc
      complex(8) :: zpsi
      integer :: target_ib, target_it, target_itp
      logical :: initialize
C     !---- local variables
      integer :: ib,ias,ics,ia,ic,nc,iv,nv,iap,icp,i,itp,it,icr
      complex(8) :: cphik_ppb

      if( initialize ) then
         if(ncc/=0 .and. ncc/=nctot) then
            stop !"psicb_v3: ncc/=0 and ncc/=ncctot"
         endif

         if(sum(ncore(iclass(1:natom)))  /= nctot) then
            stop !"psicb_v3:sum(ncore) wrong"
         end if

         return
      end if

      if(ncc == 0 ) return

      it = target_it-nctot
c$$$  do it=1, nt0
c$$$  if( nctot+it /= target_it ) cycle

      ias = 1
      ics = 0
      do ia=1, natom
         ic = iclass(ia)
         nc = nlnmc(ic)
         nv = nlnmv(ic)
         iap = iatomp(ia)
         icp = iclass(iap)
         ib = imdim(iap)-1

         do itp=1, ncore(ic)
            if( ics+itp /= target_itp ) cycle ! skip, not for this process

            do i=1, mdim(icp)
               ib=ib+1
               if( ib /= target_ib ) cycle ! skip, not for this process
               
               icr=icore(itp,ic)

               cphik_ppb = 0.0d0
!$OMP PARALLEL DO PRIVATE(iv) REDUCTION(+:cphik_ppb)
               do iv=0, nv-1
                  cphik_ppb = cphik_ppb 
     .                 + cphik(ias+iv,itp)*ppb(nc+1+iv,icr,i,icp)
               end do
!$OMP END PARALLEL
               zpsi = dconjg(phase(ia))*(cphik_ppb)
            end do

            ias = ias + nlnmv(ic)
            ics = ics + ncore(ic)
         end do
      end do
c$$$  end do

      end subroutine element_psicbB


      subroutine element_melpln2(
     .     nbloch,ngb,nctot,ncc,
     .     ngp1,ngvecp1,ngp2,ngvecp2,ngc,nadd,
     .     geig1,geig2,ngpmx,nband,itps,ntp0,nt0,
     .     shtv,q,qi,symope,qbas,
     .     zpsi,target_ib,target_it,target_itp,initialize)
C     !---- modules
      use m_read_ppovl,only:getppx, ppx,ngc2,ngvecc2
      implicit none
C     !---- arguments

      integer :: nbloch,ngb,nctot,ncc
      integer :: ngp1,ngvecp1(3,ngp1),ngp2,ngvecp2(3,ngp2),ngc,nadd(3)
      complex(8):: geig1(ngpmx,nband),geig2(ngpmx,nband),ge2
      integer:: ngpmx,nband,itps,itp,ntp0,nt0
      real(8) :: shtv(3),q(3),qi(3),symope(3,3),qbas(3,3)
      complex(8) :: zpsi
      integer :: target_ib,target_it,target_itp
      logical :: initialize
C     !---- local variables
      integer :: it,nnum,inum,ig1,ig2,igc
      integer :: iqi,igc2,ngcx,nn(1:3)
      integer :: nxmin,nymin,nzmin,nxmax,nymax,nzmax
      real(8) :: qg(3)
      complex(8),allocatable :: gg(:)
C     !---- saved variables
      integer,allocatable,save :: ipp(:,:,:)
      integer,allocatable,save :: igcs(:,:)
      complex(8),allocatable,save :: phas(:,:)
C     !---- local constants
      complex(8),parameter :: img=(0d0,1d0)
      real(8), parameter :: pi=3.1415926535897932D0

      if(ngc==0) return

      if( initialize ) then
         call getppx(qbas,ngcx,qi,symope)
         ppx=dconjg(ppx)

         if(ngc/=ngcx) then
            print *, qi, ngcx, ngc
            stop !'melpln2: ngc/= ngcx by getppx:PPOVL'
         endif

         nxmin=minval(ngvecc2(1,1:ngc2))
         nxmax=maxval(ngvecc2(1,1:ngc2))
         nymin=minval(ngvecc2(2,1:ngc2))
         nymax=maxval(ngvecc2(2,1:ngc2))
         nzmin=minval(ngvecc2(3,1:ngc2))
         nzmax=maxval(ngvecc2(3,1:ngc2))

         if( allocated(ipp) ) deallocate(ipp)
         allocate(ipp(nxmin:nxmax,nymin:nymax,nzmin:nzmax))

         ipp = -10000000
         do igc2 =1,ngc2
            ipp(ngvecc2(1,igc2),ngvecc2(2,igc2),ngvecc2(3,igc2))=igc2
         end do

         do ig1 = 1,ngp1
            do ig2 = 1,ngp2
               nn(1:3) = ngvecp1(1:3,ig1)-ngvecp2(1:3,ig2) - nadd
               if(nn(1)<nxmin.or.nn(1)>nxmax) then
                  stop !"melpln2:nn(1) out of range"
               endif   
               if(nn(2)<nymin.or.nn(2)>nymax) stop !"melpln2:nn(2) out of range"
               if(nn(3)<nzmin.or.nn(3)>nzmax) stop !"melpln2:nn(3) out of range"
            end do
         end do

         if( allocated(phas) ) deallocate(phas)
         if( allocated(igcs) ) deallocate(igcs)

         allocate(phas(ngp2,ngp1))
         allocate(igcs(ngp2,ngp1))

         do ig1 = 1,ngp1
!$OMP PARALLEL DO PRIVATE(ig2,nn,qg)
            do ig2 = 1,ngp2
               nn(1:3) = ngvecp1(1:3,ig1)-ngvecp2(1:3,ig2) - nadd
               igcs(ig2,ig1)    = ipp(nn(1),nn(2),nn(3))
               qg(1:3) = 2d0*pi * (q(1:3)+ matmul(qbas, nn(1:3)) )
               phas(ig2,ig1) = exp(+img*sum(qg*shtv))
            end do
!$OMP END PARALLEL
         end do

         return
      end if

      allocate(gg(ngc2))

      itp = target_itp-ncc
c$$$  do itp=1,ntp0
c$$$  if( ncc+itp /= target_itp ) cycle

      it = target_it-nctot
c$$$  do it =1,nt0
c$$$  if( nctot+it /= target_it ) cycle

      gg=0d0
      do ig1=1, ngp1
!$OMP PARALLEL DO PRIVATE(ig2,ge2,igc2)
         do ig2=1, ngp2
            ge2=geig2(ig2,it)*phas(ig2,ig1)
            igc2=igcs(ig2,ig1)
            gg(igc2)=gg(igc2)+dconjg(geig1(ig1,itps-1+itp))*ge2
         end do
!$OMP END PARALLEL
      end do

      do igc=1, ngc
         if( nbloch+igc /= target_ib ) cycle ! skip, not for this process

         zpsi=0.0d0
!$OMP PARALLEL DO PRIVATE(igc,igc2) REDUCTION(+:zpsi)
         do igc2=1, ngc2
            zpsi=zpsi+ppx(igc,igc2)*gg(igc2)
         end do
!$OMP END PARALLEL
      end do
c$$$  end do ! it
c$$$  end do ! itp

      deallocate(gg)

      end subroutine element_melpln2




C=====================================================================
      subroutine dpsion5(frhis,nwhis, freqr,nw_w, freqi,niwt, 
     i     realomega,   imagomega, !freqr ->frhis ...sf
     i     rcxq, npm, nw_i,nmbas1,nmbas2,
     o     zxq,zxqi,
c     i    nolfco,chipm,schi,isp,  rcxqmean,nmbas,  !iepsmode, rcxqmean, ! epsmode
     i     chipm,schi,isp,      !No nolfco mode. Apr2012.
     i     ecut,ecuts)
c     o    x0mean)
C     - Calculate W-v zxqi(on the imaginary axis) and zxq(real axis) from sperctum weight rcxq.
C     r v4 works for timereversal=F (npm=2 case).
C     r  See rcxq_zcxq for rcxq, which contains the spectrum weight for given bins along the real-axis.
C     r ! Note that zxq and zxqi are not accumlating
C     i frhis(1:nwhis+1) :: specify histgram bins i-th bin is [frhis(i), frhis(i+1)].
C     i          We suppose "freqr(i)=moddle of i-th bin; freqr(0)=0."
C     i          (I think called routine hilbertmat itself is not limited by this condition).
C     i freqr (0:nw_w) : Calcualte zxq for these real energies.
C     i freqi (1:niwt) : Calcualte zxqi for these imaginary energies.
C     i   realomega  : A switch to calculate zxq or not.
C     i   imagomega: : A switch to calculate zxqi or not.
C     iw rcxq may be altered ---used as work area.
C     io   zxq :  W-v along the real axis on freqr(0:nw_w)
C     io   zxqi:  W-v along the imag axis on freqi(niwt)
C     !
C     1 Feb2006:  v4 for timereversal=F
C     ! July2005: v3Add spin chipm mode
C     ! July2005: This version alter rcxq----it is used as work area.
C     ! sergey faleev Apr 2002 ; Rebuiled by takao
C------------------------------------------------------------------
      implicit none
      integer:: nw_w,niwt,igb1,igb2, iw,iwp,nwhis,ix,npm,ifxx,nmbas1,nmbas2
      real(8) :: freqi(niwt),pi,px,omp,om,om2,om1, !omg2max from hx0fp0
     &     frhis(nwhis+1), freqr(0:nw_w), aaa,d_omg
      logical :: realomega, imagomega
      complex(8):: rcxq(nmbas1,nmbas2, nwhis,npm) !sf 13June
c     logical   :: iepsmode
      logical :: chipm

      integer::isp,ispx      !, nmbas
c     complex(8):: rcxqmean(nwhis,npm,nmbas,nmbas)  !takao sep2006 add nmbas
C...  ecut mode
      real(8):: ecut,ecuts,wcut,wcutef,dee,schi
      logical ::debug=.false.
      real(8),allocatable :: his_L(:),his_R(:),his_C(:)
      integer it
      real(8):: domega_r,domega_c,domega_l,delta_l,delta_r
      real(8),allocatable ::rmat(:,:,:),rmati(:,:,:),rmatt(:,:,:),imatt(:,:,:)
      complex(8),allocatable :: rmatiC(:,:,:),imattC(:,:,:)
      complex(8) ::wfac
      complex(8):: zz
      complex(8),allocatable :: zxqn(:,:),zxqn1(:,:,:),rx0mean1(:,:,:),rx0mean(:)
      complex(8),allocatable:: rrr(:)

      integer::nw_i,jpm,ipm,verbose,isgi
c     complex(8):: x0mean(nw_i:nw_w,nmbas,nmbas)
      complex(8):: 
     o     zxq (nmbas1,nmbas2, nw_i: nw_w), !iw=0 means omg=0,
                                !iw=1:nw_w corresponds to iw's bit of the frequensy histogram 
     o     zxqi(nmbas1,nmbas2,niwt),img !npm), img  !zxqi(...,npm) may2006

      real(8),allocatable:: ebb(:)
      integer:: ii,i,ibas1,ibas2
      logical :: evaltest       !,testtr

c     if(verbose()>89) debug=.true.
c     ------------------------------------------------
      write(*,'(" -- dpsion5: start...   ",$)')
      write(*,"('  nw_w nwhis=',2i5)") nw_w,nwhis
      if(debug) then
         write(6,*)' nmbas1 nmbas2 nwhis npm =',  nmbas1,nmbas2,nwhis,npm
         write(6,*)' sumchk rcxq=', sum(abs(rcxq))
      endif
      pi  = 4d0*datan(1d0)
      img = (0d0,1d0)
      call cputid(0)
      ispx = isp
      if(schi<0) then
         ispx = 3-isp           !flip
      endif

C...  Check freqr=frhis_m.
C...  But I think now this is not necessary. You can supply any freqr(iw). But be careful.
      if(realomega) then
         if( nwhis <= nw_w ) then
            print *,nwhis,nw_w
            stop !' dpsion5: nwhis<=nw_w'
         endif
         if( freqr(0)/=0d0 ) stop !' dpsion5: freqr(0)/=0d0'
         aaa = 0d0
         if(nw_w>0) then
            do iw = 1,nw_w
               aaa = aaa + abs( freqr(iw) - (frhis(iw)+frhis(iw+1))/2d0 )
               if(debug) write(6,"(' iw freqr frhis_m=',i5,2f13.6)" )
     &              iw,freqr(iw),  (frhis(iw)+frhis(iw+1))/2d0
            end do
            if(aaa>1d-10)stop !'dpsion5:freqr/=frhis_m is not implimented yet'
         endif
      endif                     !realomega

C--------------------------------------------------------------
!     Each bin  [his_Left, his_Right]  his_Center is middle.
!     his_C(0) is at zero. his_R(0) and his_L(0) are not defined.
      if(debug) print *,' dpsion5: RRR 2222222222 '
      allocate(his_L(-nwhis:nwhis),his_R(-nwhis:nwhis),his_C(-nwhis:nwhis))
      his_L(1:nwhis) = frhis(  1:  nwhis)
      his_R(1:nwhis) = frhis(1+1:1+nwhis)
      his_C(1:nwhis) = (his_L(1:nwhis) + his_R(1:nwhis) )/2d0
      do iw= 1,nwhis
         his_L(-iw) = -his_R(iw)
         his_R(-iw) = -his_L(iw)
         his_C(-iw) = -his_C(iw)
      end do
      his_C(0) = 0d0; his_R(0)=-999; his_L(0)=-999
C     
      if(debug) write(6,*)'sumchk 111 rcxq=', sum(abs(rcxq))

      do iw= 1, nwhis
         if(ecut<1d9) then
            wfac= wcutef(his_C(iw), ecut,ecuts)
         else
            wfac= 1d0
         endif
!     rcxq is used as work---> rcxq= Average value of Im chi.
!     Note rcxq is "negative" (
         do jpm=1,npm
            call dscal(2*nmbas1*nmbas2, -wfac/(his_r(iw)-his_l(iw)),rcxq(1,1,iw,jpm),1)
         end do
c     if(debug) print *, 'dpsion5: RRR 7777 iw wfac=',iw,wfac,ecut,ecuts
      end do
      if(debug) write(6,*)'sumchk 122 rcxq=', sum(abs(rcxq))

C...  Temporary. maybe, we will have better procedure...
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
     &                 write(6,"('hhh1: jpm iw eb=',2i4,d13.5)") jpm,iw,ebb(ii)
               end do
            end do
         end do
         deallocate(ebb)
      endif

C---  realomega case
      if(realomega)then
         write(6,*) " --- realomega --- "
         if(npm==1) then
            allocate( rmat(0:nw_w,-nwhis:nwhis,npm), rrr(-nwhis:nwhis))
            rmat  = 0d0
            do it =  0,nw_w
               zz = freqr(it)   !his_C(it)
               call hilbertmat(zz, nwhis,his_L,his_C,his_R, rrr)
               rmat(it,:,1) = dreal(rrr)
            end do ;   if(debug) print *, 'dpsion5: RRR 55555555555'
            allocate( rmatt(0:nw_w,nwhis,npm) )
            if(     chipm.and.ispx==1 ) then
               rmatt(:,:,1) = rmat(:,1:nwhis,1)
            elseif( chipm.and.ispx==2 ) then
               do iw= 1,nwhis
                  rmatt(:,iw,1) = -rmat(:,-iw,1)
               end do
            else
               do iw= 1,nwhis
                  rmatt(:,iw,1) = rmat(:,iw,1) - rmat(:,-iw,1)
               end do
            endif
            deallocate(rmat,rrr)
         else                   ! npm==2 case -------------------------------------------------
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
            end do ;   if(debug) print *, 'dpsion5: RRR2 55555555555'
            deallocate(rrr)
         endif
         rmatt = rmatt/pi ; if(debug) print *,'dpsion5: RRR 6666'

!     ! takao remove if(nolfc) block here.
c     write(6,*) " --- realomega dgemm--- "


!     ! WARN! I think npm==2.and.chipm does not make sense. apr2012.
!     !
         if(npm==2.and.chipm) 
     &        stop !'x0kf_v4h:npm==2.and.chipm is not meaningful probably'


!     ! Note rcxq is negative now (converted at the top of this routine !!!
         if(     chipm .and. ispx==2 ) then
                                !nothing here
                                !Since the range of zxq is nw_i=0, we have no area to store negative energy part of chipm.
         elseif( chipm               ) then
            call zaxpy( nmbas1*nmbas2*nw_w, img, rcxq, 1, zxq(1,1,1), 1)
         else
            zxq = 0d0           ! not accumlating case.
            call zaxpy( nmbas1*nmbas2*nw_w, img, rcxq(1,1,1,1), 1, zxq(1,1,1), 1)
         endif

         if(npm==2) then
            do iw=1,nw_w
               call zaxpy( nmbas1*nmbas2, img, rcxq(1,1,iw,2),1, zxq(:,:,-iw),1)
            end do
         endif

         if(npm==1) then
            call dgemm('n','t',  2*nmbas1*nmbas2, nw_w+1, nwhis, 1d0, 
     &           rcxq, 2*nmbas1*nmbas2,  rmatt, nw_w+1,
     &           1d0, zxq, 2*nmbas1*nmbas2 )
         elseif(npm==2) then
            call dgemm('n','t',  2*nmbas1*nmbas2,   npm*nw_w+1, nwhis, 1d0, 
     &           rcxq(1,1,1,1), 2*nmbas1*nmbas2, rmatt(:,:,1), npm*nw_w+1,
     &           1d0, zxq, 2*nmbas1*nmbas2 )
            call dgemm('n','t',  2*nmbas1*nmbas2,   npm*nw_w+1, nwhis, 1d0, 
     &           rcxq(1,1,1,2), 2*nmbas1*nmbas2, rmatt(:,:,2), npm*nw_w+1,
     &           1d0, zxq, 2*nmbas1*nmbas2 )
         else
            stop !'dpsion5: npm=1 or 2'
         endif
         deallocate(rmatt)
      endif

!     ! === imagomega case      imatt(niwt -->niwt,npm may2005 ===
      if(imagomega) then
         allocate( rrr(-nwhis:nwhis))
         if(npm==1) then
            allocate( rmati (niwt,-nwhis:nwhis,npm))
            rmati= 0d0
         else
            allocate( rmatiC(niwt,-nwhis:nwhis,npm))
            rmatiC = 0d0
         endif ;   if(debug) print *, 'dpsion5: III 111111155555555555'
         do it =  1,niwt
            zz = img*freqi(it)  !his_C(it)
            call hilbertmat(zz, nwhis,his_L,his_C,his_R, rrr) !Im(zz)>0
            if(npm==1) then
               rmati (it,:,1) = dreal(rrr)
            else
               rmatiC(it,:,1) = rrr
            endif
         end do ;   if(debug) print *, 'dpsion5: III 55555555555'
!     ! ==== npm=1 case ====
         if(npm==1) then
            allocate( imatt(niwt, nwhis,npm) )
            do iw= 1,nwhis
               imatt(:,iw,1) = rmati(:,iw,1) - rmati(:,-iw,1)
            end do
            deallocate(rmati,rrr)
            imatt = imatt/pi; if(debug) print *, 'dpsion5: III  '
            call dgemm('n','t',  2*nmbas1*nmbas2, niwt, nwhis, 1d0, 
     &           rcxq, 2*nmbas1*nmbas2, imatt, niwt,
     &           0d0, zxqi, 2*nmbas1*nmbas2 )
            deallocate(imatt)
!     ! ==== npm=2 case ====
         else
            allocate( imattC(niwt, nwhis,npm) )
            do iw= 1,nwhis
               imattC(:,iw,1) =   rmatiC(:, iw,1)
               imattC(:,iw,2) = - rmatiC(:,-iw,1)
            end do
            deallocate(rmatiC,rrr)
            imattC = imattC/pi; if(debug) print *, 'dpsion5: IIIc '
            call zgemm('n','t',  nmbas1*nmbas2, niwt, nwhis, 1d0, 
     &           rcxq(1,1,1,1), nmbas1*nmbas2, imattC(1,1,1), niwt,
     &           0d0, zxqi,    nmbas1*nmbas2 )
            call zgemm('n','t',  nmbas1*nmbas2, niwt, nwhis, 1d0, 
     &           rcxq(1,1,1,2), nmbas1*nmbas2, imattC(1,1,2), niwt,
     &           1d0, zxqi,    nmbas1*nmbas2 )
            deallocate(imattC)
         endif
      endif
      deallocate(his_L,his_C,his_R)
      write(*,'("         end dpsion5 ",$)')
      call cputid(0)
      end subroutine dpsion5


      logical function checkbelong(qin, qall, nq,ieibz) !ieibz is also returned
      integer:: nq,ieibz
      real(8):: qin(3), qall(3,nq)
      checkbelong=.false.
      do i=1,nq
         if(sum(abs(qin-qall(:,i)))<1d-8) then
            ieibz=i
            checkbelong=.true.
            return
         endif
      end do
      end function checkbelong



