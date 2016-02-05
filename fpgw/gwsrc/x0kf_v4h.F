      subroutine x0kf_v4hz (npm,ncc,
     i            ihw,nhw,jhw,whw,nhwtot, ! tetwt5
     i            n1b,n2b,nbnbx,nbnb,     !  use whw by tetwt5 ,
     i                  q,
     i                  nsp,isp_k,isp_kq,symmetrize,
     i                  qbas,ginv,rk,wk,  
c     i                  mdim, 
     d                  nlmto,nqbz,nctot, 
c     d                  natom, 
     d                  nbloch,nwt, 
     i      iq, ngbb, ngc, ngpmx,ngcmx,  
     i      nqbze, nband,nqibz, 
     o      rcxq,
     i     nolfco,zzr,nmbas, zcousq,
     i     chipmzzr,eibzmode, 
     i     nwgt,igx,igxt,ngrp,eibzsym, crpa)
      use m_readqg,only   :  readqg
      use m_readeigen,only:  readeval 
      use keyvalue,only   :  getkeyvalue
      use m_rotMPB,only   :  rotMPB2 
      use m_readqgcou,only:  
     o qtt_, nqnum 
      use m_pkm4crpa,only : readpkm4crpa
      use m_zmel,only     : get_zmelt2,
     o zmel !,ppbir ,ppovlz

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
c              See also tetwt5. and check weight mode=4 of hx0fp0 and (mode=5,6).
c
c- takao kotani Apr 2002   This originated from Ferdi's x0k.
cr daxpy dominates the cpu time
c
c
c x0(i,j)(q,w) = S[k=FBZ] S[t=occ] S[t'=unocc]
c  <M(q,i) psi(k,t) |psi(k+q,t')> <psi(k+q,t')| psi(k,t) M(q,j)>
c  { 1/[w-e(k+q,t')+e(k,t)+i*delta] - 1/[w+e(k+q,t')-e(k,t)-i*delta] }
c  ; w is real. x0 is stored into rcxq.
c
c zzmel =  <psi(k+q,t') | psi(k,t) B(R,i)>
c zmel  =  <psi(k+q,t') | psi(k,t) M(R,i)>
c rcxq  =  zeroth order response function along the positive real axis.
c          Note this is accmulating variable. Equivalnet with zxq. See rcxq2zxq below.
c
c q       = q-vector in x(q,iw)
c ifchi   = direct access unit file for cphi, the coefficient of eigenfunction for argumentation wave.
c qbas    = base reciprocal lattice vectors
c ginv    = inverse of qbas s. indxrk.f
c
c ppb     = <phi(RLn) phi(RL'n') B(R,i)>
c
c iclass  = given an atom, tells the class
c iindxk  = index for k-points in the FBZ
c rk      = k-points in the 1st BZ
c wk      = weight for each k-point in the 1st BZ
c freq    = frequency points along positive imaginary axis
c
c
c mdim    = dimension of B(R,i) for each atom R
c nlnmx   = maximum number of l,n,m
c nlmto   = total number of LMTO basis functions
c nqbz    = number of k-points in the 1st BZ
c n1,n2,n3= divisions along base reciprocal lattice vectors
c natom   = number of atoms
c noccx   = maximum number of occupied states
c noccxv  = maximum number of occupied valence states
c nbloch  = total number of Bloch basis functions
c
c cphi_k cphi_kq:  b(k) and b(k+q)
c  : coefficients of eigenfunctions for argumentation waves in each MT
c
      implicit none
      integer(4):: npm,ncc,ngbb,natom,nwt,nsp,isp_k,isp_kq,nlmto !,noccx,noccxv
     &  ,nl,nclass,nnc,nlnmx,nbloch,iq,nqibz,iatom,nctot,nbmx,iopen !mdimx,
     &  ,jpm,ibib,itps,nt0,ntp0,ngp_kq,ngp_k,it,itp,iw,igb2,igb1,ngb
     &  ,nn,no,isx,iclose,k,nbnbx,nqbz
      real(8):: q(3),qbas(3,3),ginv(3,3),rk(3,nqbz),wk(nqbz),ebmx
c      complex   (8):: rcxq (ngbb,ngbb, nwt,npm),aaa
      complex   (8):: rcxq (nmbas,nmbas,nwt,npm)
      complex(8) :: imag=(0d0,1d0),trc,aaa !phase(natom),
      complex(8),allocatable:: cphi_k(:,:),cphi_kq(:,:),geig_kq(:,:),geig_k(:,:)
      integer(4):: ngpmx, ngcmx, nqbze, nband,
     &              ngc,nadd(3), !ngvecpB(3,ngpmx,nqbze), ngpn(nqbze),
     &             igc, !ngveccB(3,ngcmx),
     &             ngvecp_kq(3,ngpmx),ngvecp_k(3,ngpmx)
      complex(8),allocatable :: zmelt(:,:,:)!,zmelt1(:,:,:)
      real(8) :: qbasinv(3,3), det,qdiff(3),add(3),symope(3,3)
     &    ,shtv(3)=(/0d0,0d0,0d0/)
      data symope /1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0/
c      real(8) :: ppb_unused(*)

c      integer(4) :: mdim(natom)
      
      complex(8),allocatable :: ttx(:,:)
      complex(8),allocatable::  z1p(:,:)
      integer(4) ::  nbnb(nqbz,npm), 
     &  n1b(nbnbx,nqbz,npm), n2b(nbnbx,nqbz,npm)
      complex(8),allocatable:: zzmel(:,:,:)
c      integer(4)::imdim(natom),iatomp(natom)
      logical,parameter:: debug=.false.
c---tetwt5
      logical:: hist      ,ipr
      integer(4):: nhwtot,
     &  ihw(nbnbx,nqbz,npm),nhw(nbnbx,nqbz,npm),jhw(nbnbx,nqbz,npm)
      real(8):: whw(nhwtot)
      complex(8) :: zmelt1,zmelt2,zmeltt(ngbb)      !...........................sf 21May02
c         complex(8), allocatable :: zxq_(:,:,:) !...........sf 21May02
      real(8) :: imagweight !............................sf 21May02
c         logical :: takao=.false. !.........................sf 21May02
c         allocate( zxq_( nbloch + ngc,nbloch + ngc,nwt)) !..sf 21May02
      integer(4)::nocc
      real(8):: eband(nband)
c      integer(4):: n_index_qbz
c      integer(4):: index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)
c-----
c      integer(4):: nlnm(*),nlnmv(*),nlnmc(*)!,iclass(*)!,icore(*),ncore(*)
      integer(4):: verbose
c---for iepsmode
      logical   :: nolfco !iepsmode
      integer(4):: nmbas, imb1,imb2, imb !nmbas1x !nmbas2,nmbas1,
      complex(8):: zq01,zq02
c      real(8)   :: zq0zq0
      complex(8)   :: zq0zq0 !This is a bug for the case of two atoms per cell!!! oct2006
c      complex(8):: rcxqmean(nwt,npm,nmbas1,nmbas2)

      real(8):: vec_kq_g(3),vec_k_g(3),vec_kq(3),vec_k(3),quu(3),tolqu=1d-4,quu1(3),quu2(3)

      integer(4):: nbcut,nbcut2
      logical :: iww1=.true.,iww2=.true.

      logical:: smbasis
      integer(4):: ifpomat, nbloch_r, ngbo,iqxdummy
c      integer(4),allocatable:: io(:),in(:),io_q(:),in_q(:)
      complex(8),allocatable:: pomat(:,:), zmeltn(:,:,:)
      real(8):: q_r(3)
      complex(8):: img=(0d0,1d0),zzz(ngbb)

      integer(4):: nkmin,  nkmax, nkqmin, nkqmax
c      real(8):: qq(3)
      integer(4):: ib1, ib2,      ngcx,ix,iy

      integer(4),allocatable:: ngvecc(:,:)
      logical :: onceww !testtr,negative_testtr

!! takao apr2012
      logical :: z1offd !, z1stcol
      complex(8),target :: zzr(ngbb,nmbas) !ppovlz(ngbb,ngbb),
      integer:: igb
      logical:: symmetrize

!! jun2012takao
c      real(8):: qeibz(3,nqbz),  !    aik(3,3,ngrpt) 
      integer:: ngrp,nwgt(nqbz) !,ngrpt, aiktimereversal(ngrpt),nwgtieibz,ieibz
      integer:: igx(ngrp*2,nqbz),igxt(ngrp*2,nqbz),ieqbz
!! nwgt(neibz
      logical:: checkbelong,eibzmode, chipmzzr
      complex(8):: zcousq(ngbb,ngbb)  !ppovl(ngc,ngc) ,
      complex(8),allocatable:: zcousqr(:,:,:),rcxq0(:,:),rcxq00(:,:),rcxq000(:,:),rcxqwww(:,:)
c      complex(8):: zcousqrsum(ngbb,ngbb,2),    zmeltx(ngbb),zmelty(ngbb),zcousqrx(ngbb,ngbb)  
      complex(8):: zmeltx(ngbb),zmelty(ngbb),zcousqrx(ngbb,ngbb) ,zcousqc(ngbb,ngbb)
     & ,rzc(ngbb,ngbb),zcousqinv(ngbb,ngbb),cmat(ngbb,ngbb)
      integer::  eibzsym(ngrp,-1:1),neibz,icc,ig,eibzmoden,ikp,i,j,itimer,icount,iele
      integer:: irotm,nrotmx,ixx,iyy,itt,ntimer, nccc, nxx,iagain,irotm1,irotm2
      integer,allocatable:: i1(:,:),i2(:,:),nrotm(:)
      complex(8),allocatable:: zrotm(:,:),zrr(:,:),zrrc(:,:),zrr_(:,:,:),zrrc_(:,:,:),zmmm(:),zrrx(:,:)
      complex(8),pointer:: zmat(:,:)
c      complex(8),allocatable,target:: ppovl_(:,:)
c#ifdef USE_GEMM_FOR_SUM
c      complex(8),allocatable :: zmelt_tmp(:,:,:)
c#endif
      complex(8),allocatable:: rcxq_core(:,:)
!$    integer:: omp_get_num_threads
      logical:: eibz4x0
      logical :: crpa 
      real(8):: wpw_k,wpw_kq
      real(8):: vec_kcrpa(3),vec_kqcrpa(3)

      logical :: exchange=.false.
      integer:: irot=1
      integer:: ntqxx,nbmax
      

      if (symmetrize) goto 5000 

c ------------------------------------------------
!TIME0_1001 
      write(6,'(" x0kf_v4hz: q=",3f8.4,$)')q
      call cputid(0)
c$$$!! check eibzmode assumes nmbas1=nmbs2
c$$$      if(eibzmode) then
c$$$         if(nmbas1/=nmbas2) then
c$$$            write(6,*)'x0kf_v4h: eibzmode=T only allow nmbas1=nmbas2.',nmbas1,nmbas2,nmbas
c$$$            stop 'x0kf_v4h: eibzmode=T only allow nmbas1=nmbas2.'
c$$$         endif
c$$$      endif
!!
c$$$      imdim(1) = 1
c$$$      do iatom = 1,natom
c$$$        iatomp(iatom) = iatom
c$$$        if(iatom<natom) imdim(iatom+1)=imdim(iatom)+mdim(iatom)
c$$$      enddo
c      nctot      = noccx - noccxv
      call dinv33(qbas,0,qbasinv,det)
c      phase= (1d0,0d0) !coskt = 1d0; sinkt = 0d0
      allocate(cphi_k(nlmto,nband),cphi_kq(nlmto,nband), geig_kq(ngpmx,nband),  geig_k(ngpmx,nband) )
      call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
      call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
      call getnemx(nbmx,ebmx,7,.true.)
c$$$!TIME0 
c$$$      if(smbasis()) then !need to check again, when we will make smbasis on.
c$$$cccccccccccccccccccccccccccccccc
c$$$        if(ncc/=0) then
c$$$          write(6,*)"Timereversal=F(ncc/=0) not yet implemented for smbasis."
c$$$          write(6,*)" pomat should be generated correctly                  ."
c$$$Cstop2rx 2013.08.09 kino          stop "Timereversal=F(ncc/=0) not yet implemented for smbasis."
c$$$          call rx( "Timereversal=F(ncc/=0) not yet implemented for smbasis.")
c$$$        endif
c$$$ccccccccccccccccccccccccccccccccc
c$$$
c$$$        ifpomat = iopen('POmat',0,-1,0) !oct2005
c$$$C... smoothed mixed basis !oct2005
c$$$C This replace original zmelt with new zmelt based on smoothed mixed basis.
c$$$        do
c$$$          read(ifpomat) q_r,nn,no,iqxdummy !readin reduction matrix pomat
c$$$c          write(6,"('ttt: q  =',3f12.5)") q
c$$$c          write(6,"('ttt: q_r=',3f12.5)") q_r
c$$$          allocate( pomat(nn,no) )
c$$$          read(ifpomat) pomat
c$$$          if( sum(abs(q-q_r))<1d-10) then ! .and.kx <= nqibz ) then
c$$$            write(6,*) 'ok find the section for give qibz_k'
c$$$            exit
c$$$!         elseif (kx >nqibz ) then
c$$$!           exit
c$$$          endif
c$$$          deallocate(pomat)
c$$$        enddo
c$$$        if( sum(abs(q-q_r))>1d-10 ) then
c$$$          write(6,"('q  =',3f12.5)") q
c$$$          write(6,"('q_r=',3f12.5)") q_r
c$$$Cstop2rx 2013.08.09 kino          stop 'POmat reading err q/=q_r'
c$$$          call rx( 'POmat reading err q/=q_r')
c$$$        endif
c$$$        isx = iclose('POmat')
c$$$      endif
c$$$!TIME1 "after if smbasis"

ckino
!KINO      it=0
!KINO      do k=1,nqbz
!KINO        if(eibzmode.and.nwgt(k)==0 ) cycle
!KINO        if(sum(nbnb(k,1:npm))==0) cycle
!KINO        it=it+1
!KINO      enddo
!KINO      write(6,'(a,i3,1x,a,i4)')'iq=',iq,'active-k-points=',it
ckinoend

!TIME1_1001 "beforedo1000"
!! loop over k-points -------------------------------------------------------------------------
c      qq=0d0
      do 1000 k = 1,nqbz
        if(eibzmode.and.nwgt(k)==0 ) cycle
        if(debug) write(6,'("do 1000  k=",i4)')k
        ipr=(k<5.or.k==nqbz.or.debug)
        if(sum(nbnb(k,1:npm))==0) cycle
!TIME0_1101
        if(k<=5.or. (mod(k,max(10,nqbz/20))==1.or.k>nqbz-10) ) then
          write(6,"('  x0kf_v4hz: k rk=',i7,3f10.4,$)")k, rk(:,k)
          call cputid(0)
        endif

!! --- tetra ------ override nt0 itps ntp0 ---------------
        nkmin = 999999
        nkmax= -999999
        nkqmin= 999999
        nkqmax=-999999
        do jpm=1,npm !npm
          do ibib = 1, nbnb(k,jpm)
            nkmin  = min(n1b(ibib,k,jpm),nkmin)
            nkqmin = min(n2b(ibib,k,jpm),nkqmin)
            if(n1b(ibib,k,jpm)<=nband)   nkmax = max(n1b(ibib,k,jpm),nkmax)
            if(n2b(ibib,k,jpm)<=nband)  nkqmax = max(n2b(ibib,k,jpm),nkqmax)
          enddo
        enddo
        call readeval(q+rk(:,k),isp_kq,eband)
        nkqmax = nocc (eband,ebmx,.true.,nband)
        if(npm==2) then
          call readeval(rk(:,k),isp_k,eband)
          nkmax = nocc (eband,ebmx,.true.,nband)
        endif
        itps  = nkqmin          ! nkqmin = the num of min   n2 =unocc for jpm=1
        nt0   = nkmax
        ntp0  = nkqmax - nkqmin +1
        if( npm==2.and. nkqmin/=1) then
          write(6,*)' npm==2 nkqmin nkqmax  nkmin nkmax=',nkqmin,nkqmax,nkmin,nkmax
          call rx( " When npm==2, nkqmin==1 should be.")
        endif
        if(nkmin/=1) then
          call rx( " nkmin==1 should be.")
        endif

!... feb2006
!  zzmel(1:nbloch, ib_k,ib_kq)
!      ib_k =[1:nctot]              core
!      ib_k =[nctot+1:nctot+nkmax]  valence
!      ib_kq =[1:nctot]             core
!      ib_kq =[ncc+1:ncc+nkqmax - nkqmin] valence range [nkqmin,nkqmax]
!   If jpm=1, ncc=0.
!   If jpm=2, ncc=ncore. itps=1 should be.
! There is a little confusion. n1b index contains cores are after valence.
! You can see codes to treat the confusion.
! NOTE: 
!  q+rk n2b vec_kq  vec_kq_g geig_kq cphi_kq  ngp_kq ngvecp_kq  isp_kq
!    rk n1b vec_k   vec_k_g  geig_k  cphi_k   ngp_k  ngvecp_k   isp_k

c        if(ipr) then
c          write(6,"(' nkRange  nkqRange=',2i6,2x,2i6)") nkmin,nkmax,nkqmin,nkqmax
c        endif




c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$        goto 8828
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$!TIME1 "before readphi"
c$$$!TIME0 
c$$$!!--- calculate the matrix elements <psi(k+q,t') | psi(k,t) B(R,i)>
c$$$!! Read cphi part of eigenfunctions for k and k+q
c$$$        call  readcphi(q+rk(:,k)-qq, nlmto,isp_kq, quu2,cphi_kq)
c$$$        call  readcphi(  rk(:,k)-qq, nlmto,isp_k, quu1,cphi_k) !quu is used q for eigenfunctions.
c$$$!! ... core
c$$$        if(debug)  call cputid(0)
c$$$!!        allocate( zzmel(nbloch,noccx, ntp0) )
c$$$!!                         q        k    q+k
c$$$        if(debug)  write(6,*)nbloch,nctot,nt0,ncc,ntp0
c$$$        allocate( zzmel(nbloch, nctot+nt0, ncc+ntp0) )
c$$$        if(debug)  write(6,"('zzzw nkmin nkqmin=',2i5)") nkmin,nkqmin
c$$$        if(onceww(5) ) write(6,*)' nctot ncc=',nctot,ncc
c$$$c      allocate( ppb(nlnmx*nlnmx*mdimx*nclass))
c$$$c      ppb=ppbir(:,irot,ispq)
c$$$        call psicb_v3  ( nctot,ncc,nt0,ntp0,iclass,phase,
c$$$     i              cphi_k (1,nkmin),  
c$$$     i              cphi_kq(1,nkqmin), 
c$$$     i              ppbir(:,irot,isp_k),!ppb,            
c$$$     i              nlnmv,nlnmc,mdim, 
c$$$     i              imdim,iatomp,
c$$$     i              mdimx,nlmto,nbloch,nlnmx,natom,nclass,
c$$$     i              icore,ncore,nl,nnc,
c$$$     o              zzmel)
c$$$C ... valence
c$$$      if(debug) write(6,'("4 zzzqqbbb222 ",3d13.5)') sum(abs(zzmel)),sum(zzmel)
c$$$        !call cputid(0)
c$$$        call psi2b_v3  ( nctot,ncc,nt0,ntp0,iclass,phase,
c$$$     i              cphi_k(1,nkmin), 
c$$$     i              cphi_kq(1,nkqmin),
c$$$     i              ppbir(:,irot,isp_k),! ppb,         
c$$$     i              nlnmv, nlnmc,mdim,
c$$$     i              imdim,iatomp,
c$$$     d              mdimx,nlmto,nbloch,nlnmx,natom,nclass,
c$$$     o              zzmel)
c$$$        if(debug) call cputid(0)
c$$$        if(debug) write(6,'("4 zzzqqbbb222 ",3d13.5)') sum(abs(zzmel)),sum(zzmel)
c$$$!TIME1 "after psicb_v3"
c$$$!TIME0 
c$$$!! --- IPW set
c$$$        call readqg('QGpsi',q+rk(:,k)-qq,ginv, vec_kq, ngp_kq, ngvecp_kq)
c$$$        call readqg('QGpsi',  rk(:,k)-qq,ginv, vec_k,  ngp_k,  ngvecp_k)
c$$$!!      ngp_kq = ngpn(kp)  ! q+k   ntp0 in FBZ
c$$$!!      ngp_k = ngpn(k)   ! k     np0  in FBZ
c$$$!!      ngc              ! q          in IBZ
c$$$        ngb  = nbloch + ngc ! This is not ngbb for smbasis()=T. oct2005
c$$$        if(ngb/=ngbb) then
c$$$           write(6,*)' x0kf_v4h: ngb ngbb=',ngb,ngbb
c$$$           call rx( 'x0kf_v4h: ngb/=ngbb')
c$$$        endif
c$$$!!                        q       k        q+k
c$$$        allocate( zmel(ngb,nctot+nt0,ncc+ntp0) )
c$$$        allocate( z1p(ngb,ngb) )
c$$$!! ... read plane wave part of eigenfunction. (note isp is opposite).
c$$$        call readgeig(q+rk(:,k)-qq, ngpmx,isp_kq, vec_kq_g, geig_kq)
c$$$        call readgeig(  rk(:,k)-qq, ngpmx,isp_k,  vec_k_g,  geig_k)
c$$$        if(sum(abs(vec_kq_g-vec_kq))>tolqu) then
c$$$          write(6,"('vec_kq_g :',3d23.10)") vec_kq_g
c$$$          write(6,"('vec_kq  :',3d23.10)") vec_kq
c$$$          call rx( 'x0kf_v4hz:vec_kq_g/=vec_kq')
c$$$        endif
c$$$        if(sum(abs(vec_k_g-vec_k))>tolqu) then
c$$$          write(6,"('vec_k_g :',3d23.10)") vec_k_g
c$$$          write(6,"('vec_k:  ',3d23.10)") vec_k
c$$$          call rx( 'x0kf_v4hz:vec_k_g/=vec_k')
c$$$        endif
c$$$!!      qdiff = q   - qbkp(:) + rk(:,k)
c$$$        qdiff = q    - vec_kq    + vec_k
c$$$        ! q   - (q+k)   + k  is not zero.
c$$$        ! qc  -  q1     + q2
c$$$        add   = matmul(qbasinv, qdiff)
c$$$        nadd  = idint( add + dsign(.5d0,add))  !  write(6,*)' qdif=',qdiff,qbkp(:),rk(:,k)
c$$$        if(sum(abs(add-nadd))>1d-10) call rx( "sexc: abs(add-nadd))>1d-10")
c$$$        zmel = 0d0
c$$$!TIME1 "before melpln2t"
c$$$!TIME0  
c$$$!! <Bq Pq2|Pq1> = < Bq Pqu2|  Pqu1> *exp(i2pi nadd)
c$$$        if(ngc/=0) then !Aug2005
c$$$          call melpln2t(ngp_kq, ngvecp_kq  ! q+k  ; kp ngp_kq 1:ntp0 q-point
c$$$     &             , ngp_k, ngvecp_k       ! k    ; k  ngp_k 1:nt0  occupied
c$$$     &           , ngc,  nadd,
c$$$     &       geig_kq(1:ngp_kq, itps:itps+ntp0-1), ntp0, ! q+k  ; kq
c$$$     &       geig_k(1:ngp_k, 1:nt0     ),  nt0, !    ; k
c$$$     i       shtv,  q, q,  symope,qbas,
c$$$     i       vec_kq, !qt oct2013
c$$$     o       zmel (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,ncc+1:ncc+ntp0))
c$$$        endif
c$$$!! == zmelt contain O^-1=<I|J>^-1 factor. zzmel(J,it,itp)= \sum_I <phi phi|I> O^-1_IJ ==
c$$$        zmel(1:nbloch, 1:nctot+nt0, 1:ncc+ntp0) = 
c$$$     &  zzmel(1:nbloch, 1:nctot+nt0, 1:ncc+ntp0)
c$$$!!                        k            q+k
c$$$        deallocate(zzmel)
c$$$         if(debug) write(6,'("4 zzzppp222bbb ",3d13.5)') sum(abs(zmel)),sum(zmel)
c$$$        if(debug) call cputid(0)
c$$$!TIME1 "after melpln2t"
c$$$!TIME0 
c$$$!! === zmelt conversion on different basis. 
c$$$        allocate(zmmm(nmbas))   ! this is also obsolete if USE_GEMM_FOR_SUM
c$$$        if(chipmzzr) then       !spin moment basis.
c$$$          zmat => zzr
c$$$        elseif(nolfco .and. nmbas==1) then !for <e^iqr|x0|e^iqr>
c$$$          zmat => ppovlz
c$$$        else                    !may2013  this removes O^-1 factor from zmelt
c$$$          allocate(ppovl_(ngb,ngb))
c$$$          ppovl_=0d0
c$$$          do i=1,nbloch
c$$$            ppovl_(i,i)=1d0
c$$$          enddo
c$$$          ppovl_(nbloch+1:nbloch+ngc,nbloch+1:nbloch+ngc)=ppovl
c$$$          if(.not.eibz4x0()) then !sep2014 added for eibz4x0=F
c$$$            ppovl_=matmul(ppovl_,zcousq)
c$$$          endif
c$$$          zmat => ppovl_
c$$$        endif
c$$$!! :: zmelt conversion muplitpled by zzr.
c$$$        if(verbose()>39) write(6,*)'info: USE GEMM FOR SUM (zmelt = zmelt*zmat) in x0kf_v4h.F'
c$$$        allocate( zmelt_tmp(ngb,nctot+nt0,ncc+ntp0) )
c$$$        call zcopy(ngb*(nctot+nt0)*(ncc+ntp0),zmel,1,zmelt_tmp,1)
c$$$        call zgemm('T','N',ngb,(nctot+nt0)*(ncc+ntp0),ngb,(1d0,0d0),
c$$$     .   zmat,ngb,zmelt_tmp,ngb,(0d0,0d0),zmel,ngb)
c$$$        deallocate(zmelt_tmp)
c$$$        deallocate(zmmm)
c$$$        if(allocated(ppovl_)) deallocate(ppovl_)
c$$$!TIME1  "after matmul zmel"
c$$$c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$        print *,'xxxxxxxxxxxx 8829 xxxxxxxxxxxxxx'
c$$$         goto 8829
c$$$c        if(debug) write(6,'("4 zzzppp 111 ",3d15.6)') sum(abs(zmel)),sum(zmel)
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$c$$$!!march2013--> this if branch of nolfco is now unified to do 25 
c$$$c$$$!! No LocalFieldCorrection mode 
c$$$c$$$        if(nolfco) then !iepsmode==202) then ! just for <exp(iq r)|x0(q,\omega)|exp(iq r)>
c$$$c$$$          do jpm  = 1,npm !
c$$$c$$$            do ibib = 1, nbnb(k,jpm) !---  ibib loop
c$$$c$$$              if(jpm==1) then
c$$$c$$$                if( n1b(ibib,k,jpm) <= nbcut.and. n2b(ibib,k,jpm)>nbcut2) then !oct2005
c$$$c$$$                  if(iww2) then
c$$$c$$$                    write(6,"(' nband_chi0 nbcut nbcut2 n2b n1b=',4i6)") nbcut,n2b(ibib,k,jpm),n1b(ibib,k,jpm)
c$$$c$$$                    iww2=.false.
c$$$c$$$                  endif
c$$$c$$$                  cycle
c$$$c$$$                endif
c$$$c$$$              else
c$$$c$$$                if( n2b(ibib,k,jpm) <= nbcut.and. n1b(ibib,k,jpm)>nbcut2) then !oct2005
c$$$c$$$                  if(iww2) then
c$$$c$$$                    write(6,"(' nband_chi0 nbcut nbcut2 n2b n1b=',4i6)") nbcut,n2b(ibib,k,jpm),n1b(ibib,k,jpm)
c$$$c$$$                    iww2=.false.
c$$$c$$$                  endif
c$$$c$$$                  cycle
c$$$c$$$                endif
c$$$c$$$              endif
c$$$c$$$
c$$$c$$$              if( jpm==1.and.n2b(ibib,k,jpm) > nbmx) cycle
c$$$c$$$              if( jpm==2.and.n1b(ibib,k,jpm) > nbmx) cycle
c$$$c$$$
c$$$c$$$              if( n1b(ibib,k,jpm) <= nband) then
c$$$c$$$                it = nctot + n1b(ibib,k,jpm) !valence
c$$$c$$$              else
c$$$c$$$                it = n1b(ibib,k,jpm) - nband !core
c$$$c$$$              endif
c$$$c$$$
c$$$c$$$              if( n2b(ibib,k,jpm) <= nband) then
c$$$c$$$                itp = ncc + n2b(ibib,k,jpm) - itps + 1 !val
c$$$c$$$                if(itp > ncc + nkqmax-itps+1 ) cycle
c$$$c$$$              else
c$$$c$$$                itp =       n2b(ibib,k,jpm) - itps + 1 - nband !core
c$$$c$$$              endif
c$$$c$$$
c$$$c$$$              do imb2=1,nmbas
c$$$c$$$                zq02 = zmelt(imb2,it,itp)
c$$$c$$$                do imb1=1,imb2
c$$$c$$$                  zq01 = zmelt(imb1,it,itp)
c$$$c$$$                  zq0zq0 = dconjg(zq01)*zq02
c$$$c$$$                  do iw = ihw(ibib,k,jpm),ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1
c$$$c$$$                    if (iw .gt. nwt) stop "x0kf_v4hz: iw > nwt"
c$$$c$$$                    !iiww=iw+ihw(ibib,k)-1
c$$$c$$$                    imagweight   = whw(jhw(ibib,k,jpm)+iw-ihw(ibib,k,jpm))
c$$$c$$$                    if(eibzmode) imagweight = nwgt(k)*imagweight  
c$$$c$$$c                    rcxqmean(iw,jpm,imb1,imb2) =  ! here we  sum over ibib (or n, n') and k.
c$$$c$$$c     &              rcxqmean(iw,jpm,imb1,imb2) + zq0zq0*imagweight
c$$$c$$$                    rcxq(imb1,imb2,iw,jpm) =  ! here we  sum over ibib (or n, n') and k.
c$$$c$$$     &              rcxq(imb1,imb2,iw,jpm) + zq0zq0*imagweight !sum over spin in hx0fp0
c$$$c$$$                  enddo ! iw
c$$$c$$$                enddo ! imb1
c$$$c$$$              enddo ! imb2
c$$$c$$$            enddo ! ----- ibib loop
c$$$c$$$          enddo ! ----- jpm loop
c$$$c$$$          deallocate(zmelt,z1p)
c$$$c$$$          cycle !cycye do 1000 here
c$$$c$$$        endif  
c$$$c$$$!TIME1 "before jpm ibib loop"
c$$$c$$$!TIME0
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$ 8828   continue


c$$$ this (ppovlz generation) is moved to hx0fp0 and/or hx0fp0_sc.
c$$$!! === zmelt conversion on different basis. 
c$$$        if(chipmzzr) then       !spin moment basis.
c$$$c          if(allocated(ppovlz)) deallocate(ppovlz)
c$$$c          allocate(ppovlz(ngb,nmbas))
c$$$c          ppovlz= zzr
c$$$        elseif(nolfco .and. nmbas==1) then !for <e^iqr|x0|e^iqr>
c$$$          continue
c$$$        else                    !may2013  this removes O^-1 factor from zmelt
c$$$          allocate(ppovl_(ngb,ngb))
c$$$          ppovl_=0d0
c$$$          do i=1,nbloch
c$$$            ppovl_(i,i)=1d0
c$$$          enddo
c$$$          ppovl_(nbloch+1:nbloch+ngc,nbloch+1:nbloch+ngc)=ppovl
c$$$          if(.not.eibz4x0()) then !sep2014 added for eibz4x0=F
c$$$            ppovl_=matmul(ppovl_,zcousq)
c$$$          endif
c$$$          ppovlz = ppovl_
c$$$          deallocate(ppovl_)
c$$$        endif
c$$$        if(allocated(zmel)) deallocate(zmel)
c        nbmax= nctot+nt0
c        ntqxx= ncc+ntp0
cc        call get_zmelt(exchange,q,kx,qibz_k,irot,qbz_kr,kr,isp,
cc     &       ngc,ngb,nbmax,ntqxx,isp_k,isp_kq)
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
        ngb  = nbloch + ngc 
!!Get the matrix element zmel   ZO^-1 <MPB psi|psi> , where ZO is ppovlz
!!  Output is zmel(ngb, nctot+nt0,ncc+ntp0)  nkmin:nt0, nkqmin:ntp0
! nt0=nkmax-nkmin+1  , ntp0=nkqmax-nkqmin+1
        call get_zmelt2(exchange,  
     &      q,irot,q,ngc,ngb,                       !MPB
     &                   nkmin, nkmax, isp_k,nctot, !middle state 1:nt0   --> true index of eigen is mkmin:mkmin+nt0-1     + nctot
     &      q+rk(:,k),nkqmin,nkqmax,isp_kq,ncc)  !end state    1:ntp0  -->                     is nkqmin:nkqmin+ntp0-1  
        zmel = dconjg(zmel)
        allocate( z1p(ngb,ngb) )

 8829   continue
        if(debug) write(6,'("4 zzzppp 222 ",3d15.6)') sum(abs(zmel)),sum(zmel)

!TIME1_1101 "before_get_zmelt2"
!TIME0_1201
c----------------------------
!!  z1p = <M_ibg1 psi_it | psi_itp> < psi_itp | psi_it M_ibg2 >
!!  zxq(iw,ibg1,igb2) = sum_ibib wwk(iw,ibib)* z1p(ibib, igb1,igb2)
!KINO        write(6,'(a,i4)')'kino: npm=',npm
!kino 2014-08-13  !$OMP parallel private(it,itp,iww1,iww2,zmelt2,imagweight)
        do 25 jpm  = 1, npm !
        do 25 ibib = 1, nbnb(k,jpm) !---  ibib loop
!KINO          write(6,'(a,5i8)')'kino: ngb,hw(ibib,k,jpm),ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1=',
!KINO&      ngb,ihw(ibib,k,jpm),ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1
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
            if(n2b(ibib,k,jpm)>nbmx)  then  !nbmx
              if(iww1) then
                write(6,*)' nband_chi0 nbmx=',nbmx
                iww1=.false.
              endif
              cycle
            endif
            if( n1b(ibib,k,jpm) <= nbcut .and. nbcut2<n2b(ibib,k,jpm) ) then
              if(iww2) then
                write(6,"(' nband_chi0 nbcut nbcut2 n2b n1b=',4i6)") nbcut,n2b(ibib,k,jpm),n1b(ibib,k,jpm)
                iww2=.false.
              endif
              cycle
            endif

          else !jpm==2 ------------------------------
            if( n1b(ibib,k,jpm) > nbmx) then  !nbmx
              if(iww1) then
                write(6,*)' nband_chi0 nbmx=',nbmx
                iww1=.false.
              endif
              cycle
            endif
            if( n2b(ibib,k,jpm) <= nbcut .and. nbcut2<n1b(ibib,k,jpm) ) then
              if(iww2) then
                write(6,"(' nband_chi0 nbcut nbcut2 n2b n1b=',4i6)") nbcut,n2b(ibib,k,jpm),n1b(ibib,k,jpm)
                iww2=.false.
              endif
              cycle
            endif
          endif

cccccccccccccccccccccccccccccccccc takao variant begin
cc         if(takao) then
cc
cc         do ic = 1,ngb
cc           z1p(1:ngb,ic) =
cc     &     zmelt(ic,it,itp)*dconjg(zmelt(1:ngb,it,itp))
cc         end do
cc         ihww = ihw(ibib,k)
cc
c1ini------------
cc         do iw = 1, nhw(ibib,k)
cc           rviw = whw(jhw(ibib,k)+iw-1)
cC ... this part dominates the cpu time --------------------!
c!          call zaddr_(zxq(1,1,ihww+iw-1),rviw,z1p,ngb**2)
cc           call daxpy(ngb**2*2,rviw,z1p,1,
cc     &                 zxq(1,1,ihww+iw-1),1)
cc         enddo
c1end---------
c2ini --------
cc         call  rcxq_zxq(rc1p,z1p,ngb,-1)
cc         do iw = 1, nhw(ibib,k)
cc           rviw = whw(jhw(ibib,k)+iw-1)
C ... this part dominates the cpu time --------------------!
!          call zaddr_(rcxq(1,1,ihww+iw-1),rviw,z1p,ngb**2)
cc           call daxpy(ngb**2,rviw,rc1p,1,
cc     &                rcxq(1,1,ihww+iw-1),1)
cc         enddo
c2end --------------
cc         else
ccccccccccccccccccccccccc takao variant end

c$$$          if(newanisox.and.eibzmoden==1) then ! This is slow.
c$$$            zmeltx =  zmelt(:,it,itp)
c$$$            z1p=0d0
c$$$            do ieqbz =1, nwgt(k) !equivalent points for ieibz
c$$$                                 !igx,igxt specifies space-group operation (including ID)
c$$$              call rotMPB(zcousq,nbloch,ngbb,q,igx(ieqbz,k),igxt(ieqbz,k),ginv,zcousqrx) !zcousqr=Rotate_igx(zcousq) 
c$$$              zmelty = matmul(zmeltx,zcousqrx)
c$$$              do igb2=1, ngb !....................................
c$$$                zmelt2 = zmelty(igb2)
c$$$              do igb1=1,igb2
c$$$                z1p(igb1,igb2) = z1p(igb1,igb2) + dconjg(zmelty(igb1)) * zmelt2
c$$$              enddo
c$$$              enddo
c$$$            enddo  
c$$$          else
Cstop2rx 2013.08.09 kino          if (ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1 >nwt) stop "x0kf_v4hz: iw>nwt"
          if (ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1 >nwt) call rx( "x0kf_v4hz: iw>nwt")


!kino 2014-08-13  !$OMP do private(zmelt2)
          do igb2=1, nmbas     !....................................
             zmelt2 = zmel(igb2,it,itp) !zmelt(igb2,it,itp)
             do igb1=1,igb2
c                z1p(igb1,igb2) = dconjg(zmelt(igb1,it,itp)) * zmelt2
                z1p(igb1,igb2) = dconjg(zmel(igb1,it,itp)) * zmelt2
             enddo
          enddo

!! ---------------------------------
          if(crpa) then
c            print *,'readout readqkm init'
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
c            print *,'readout readqkm=',wpw_k,wpw_kq
          endif
          
ccccccccccccccccccccccccccccccc
!kino 2014-08-13  !$OMP end do
c$$$  endif

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

c$$$cccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c        if(ipr) then
c$$$          do jpm=1,npm
c$$$            write(6,"(' k jpm sum(rcxq) ngb ngbb=',2i5,2d23.15,2i8)")
c$$$     &      k,jpm,sum(rcxq(:,:,:,jpm)),ngb,ngbb
c$$$          enddo
c$$$          do ib1 =1,ngbb
c$$$             if(ib1<4) then
c$$$             elseif(ib1>ngbb-3) then
c$$$             else
c$$$                cycle
c$$$             endif
c$$$          do ib2 =1,ngbb
c$$$             if(ib2<4) then
c$$$             elseif(ib2>ngbb-3) then
c$$$             else
c$$$                cycle
c$$$             endif
c$$$          do iw =1,nwt
c$$$            write(6,"('uuu: k iw ib1 ib2 sum(rcxq)=',4i5,4d23.15)")
c$$$     &       k,iw,ib1,ib2,(rcxq(ib1,ib2,iw,1)), (rcxq(ib1,ib2,iw,2))
c$$$          enddo
c$$$          enddo
c$$$          enddo
c$$$c        endif
c$$$cccccccccccccccccccccccccccccccccccccccccccccc
        deallocate(z1p,zmel) !zmelt,z1p)
        if(debug) call cputid(0)
        if(debug) write(6,*)' end of kloop k jpm=',k,jpm
 1000 continue

!! Not need to be symmetrized
      if(nolfco .and. nmbas==1) then
         write(6,*)' nmbas=1 nolfco=T ---> not need to symmetrize'
         goto 9999
      endif
!TIME0_1301


!! ==== Hermitianize. jun2012takao moved from dpsion5 ====
c      if(eibzmode) then !comment out sep2014. 
        do jpm=1,npm
        do iw= 1,nwt
         do igb2= 1,nmbas !eibzmode assumes nmbas1=nmbas2
         do igb1= 1,igb2-1
           rcxq(igb2,igb1,iw,jpm) = dconjg(rcxq(igb1,igb2,iw,jpm))
         enddo
         enddo
        enddo
        enddo
c       endif
!TIME1_1301 "before_eibzmode_symmetrization"

!! == End of eibzmode=F (no symmetrization required). ==
      goto 9999 ! finalize



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
 5000 continue
!! === zmelt conversion === 
      if(nolfco .and. nmbas==1) then
         write(6,*)' nmbas=1 nolfco=T ---> not need to symmetrize'
         goto 9999
      endif
!!
      if(eibzmode.and.symmetrize) then
        ngb  = nbloch + ngc ! This is not ngbb for smbasis()=T. oct2005
        if(ngb/=ngbb) then
           write(6,*)' x0kf_v4h: ngb ngbb=',ngb,ngbb
           call rx( 'x0kf_v4h: ngb/=ngbb')
        endif

!TIME0_1401 
        call iqindx2(q, ginv, qtt_, nqnum, ikp,quu) !to get ikp for timereversal mode
!TIME1_1401 "after_iqindx2"
!TIME0_1501
        if(sum(abs(q-quu))>1d-8) call rx( 'x0kf_v2h: eibz 111 q/quu')
        neibz = sum(eibzsym(:,1))+sum(eibzsym(:,-1)) 
             !itimer=-1 means time reversal. eibzsym(ig,itimer) where ig: space rotation.
        write(6,"(' --- goto symmetrization --- ikp neibz q=',2i3,3f12.8)")ikp,neibz,q
        call cputid2(' --- x0kf: start symmetrization  ',0)

c        allocate(rcxq0(ngb,ngb),rcxq00(ngb,ngb),rcxq000(ngb,ngb),rcxqwww(ngb,ngb))
        ntimer=1
        if(sum(eibzsym(:,-1))>0) ntimer=2 !timereversal case
        allocate(zrotm(ngb,ngb),nrotm(ngrp*2))
!!
        zcousqinv=zcousq
        call matcinv(ngb,zcousqinv)

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
        call cputid2(' --- x0kf:qqqqq222ini:',0)
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
        call cputid2(' --- qqqqq222end:',0)
      endif  

 9999 continue
!kino 2014.08.19 use automatic deallocation,      deallocate(cphi_k,cphi_kq,geig_kq,geig_k)
c      if(smbasis()) deallocate(pomat)
      write(6,"(' --- x0kf_v4hz: end')")
      end subroutine x0kf_v4hz




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

C... Check freqr=frhis_m.
C... But I think now this is not necessary. You can supply any freqr(iw). But be careful.
        if(realomega) then
          if( nwhis <= nw_w ) then
            write(6,*)nwhis,nw_w
Cstop2rx 2013.08.09 kino            stop ' dpsion5: nwhis<=nw_w'
            call rx( ' dpsion5: nwhis<=nw_w')
          endif
Cstop2rx 2013.08.09 kino          if( freqr(0)/=0d0 ) stop ' dpsion5: freqr(0)/=0d0'
          if( freqr(0)/=0d0 ) call rx( ' dpsion5: freqr(0)/=0d0')
          aaa = 0d0
          if(nw_w>0) then
            do iw = 1,nw_w
              aaa = aaa + abs( freqr(iw) - (frhis(iw)+frhis(iw+1))/2d0 )
              if(debug) write(6,"(' iw freqr frhis_m=',i5,2f13.6)" )
     &     iw,freqr(iw),  (frhis(iw)+frhis(iw+1))/2d0
            enddo
Cstop2rx 2013.08.09 kino            if(aaa>1d-10)stop 'dpsion5:freqr/=frhis_m is not implimented yet'
            if(aaa>1d-10)call rx( 'dpsion5:freqr/=frhis_m is not implimented yet')
          endif
        endif !realomega

C--------------------------------------------------------------
! Each bin  [his_Left, his_Right]  his_Center is middle.
! his_C(0) is at zero. his_R(0) and his_L(0) are not defined.
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
        real(8):: qin(3), qall(3,nq)
        checkbelong=.false.
        do i=1,nq
           if(sum(abs(qin-qall(:,i)))<1d-8) then
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
