c------------------------------------------
      subroutine x0kf_v4hz (npm,ncc,
     i            ihw,nhw,jhw,whw,nhwtot, ! tetwt5
     i            n1b,n2b,nbnbx,nbnb,     !  use whw by tetwt5 ,
     i                  q,
     i                  nsp,isp1,isp2,      !ifcphi,ifrb,ifcb,ifrhb,ifchb,
     i                  qbas,ginv,rk,wk,  !,iindxk,
     i                  mdim, nlnm, nlnmv,nlnmc, iclass, ! not the w(pointer)
     i                  ppb, 
     i                  icore,ncore, ! not the w(pointer)
     d                  nlmto,nqbz,noccx,noccxv,
     d                  nl,nclass,natom,nnc,
     d                  nlnmx,mdimx,nbloch,nwt,
c
c     i      iq,  ngpn,ngc, ngpmx,ngcmx,  geigB, jan2004
c     i       ngvecpB,ngveccB,
     i      iq, ngbb, ngc, ngpmx,ngcmx,  !oct2005 ngbb is not ngc + nbloch for smbasis=T.
     i      nqbze, nband,nqibz, 
     o      rcxq,
c------------------------------------------
     i     nolfco,zzr,nmbas,nmbas1,nmbas2, ppovlz, ppovl,zcousq, !iepsmode,zzr, nmbas added Sep2006
     i     newanisox,z1offd,
     o     rcxqmean,
     i     nwgt,igx,igxt,ngrp,eibzsym)!, zcousqrsum) !aik,aiktimereversal,ngrpt)
      use m_readqg,only :  readqg
      use m_readeigen,only: readeval,readcphi,readgeig
      use keyvalue,only:  getkeyvalue
      use m_rotMPB,only:  rotMPB2 !rotMPB,
      use m_readqgcou,only: qtt_,nqnum !itimermap,

!! when nmbas1=2, this works in a special manner for nolfco=T and newanisox=T. mar2012takao

C- Only the spectrum function (imaginary part) along the real axis.
Cr This version x0kf_v2hx is for timereversl=F (npm=2). dec2006
c   zxq(i,j,iw) = (-1/pi) \int_{freqr(iw)}^{freqr(iw+1)} d \omega Img<M(i)|x0(q,\omega)|M(j)>
c              for a given q and spin. Stored into rcxq equivalent with zxq.
c              See rcxq_zxq below.
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
c cphik cphikq:  b(k) and b(k+q)
c  : coefficients of eigenfunctions for argumentation waves in each MT
c
      implicit none
      integer(4):: npm,ncc,ngbb,natom,nwt,nsp,isp1,isp2,nlmto,noccx,noccxv
     &  ,nl,nclass,nnc,nlnmx,mdimx,nbloch,iq,nqibz,iatom,nctot,nbmx,iopen
     &  ,jpm,ibib,itps,nt0,ntp0,ngp1,ngp2,it,itp,iw,igb2,igb1,ngb
     &  ,nn,no,isx,iclose,k,nbnbx,nqbz
      real(8):: q(3),qbas(3,3),ginv(3,3),rk(3,nqbz),wk(nqbz),ebmx
      complex   (8):: rcxq (ngbb,ngbb, nwt,npm),aaa
      complex(8) :: imag=(0d0,1d0),trc,phase(natom)
      complex(8),allocatable:: cphik(:,:),cphikq(:,:),geig1(:,:),geig2(:,:)
      integer(4):: ngpmx, ngcmx, nqbze, nband,
     &              ngc,nadd(3), !ngvecpB(3,ngpmx,nqbze), ngpn(nqbze),
     &             igc, !ngveccB(3,ngcmx),
     &             ngvecpB1(3,ngpmx),ngvecpB2(3,ngpmx)
      complex(8),allocatable :: zmelt(:,:,:)!,zmelt1(:,:,:)
      real(8) :: qbasinv(3,3), det,qdiff(3),add(3),symope(3,3)
     &    ,shtv(3)=(/0d0,0d0,0d0/)
      data symope /1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0/
      real(8) :: ppb(*)

      integer(4) :: mdim(natom)
      
      complex(8),allocatable :: ttx(:,:)
      complex(8),allocatable::  z1p(:,:)
      integer(4) ::  nbnb(nqbz,npm), 
     &  n1b(nbnbx,nqbz,npm), n2b(nbnbx,nqbz,npm)
      complex(8),allocatable:: zzmel(:,:,:)
      integer(4)::imdim(natom),iatomp(natom)
      logical,parameter:: debug=.true.
c---tetwt5
      logical:: hist      ,ipr
      integer(4):: nhwtot,
     &  ihw(nbnbx,nqbz,npm),nhw(nbnbx,nqbz,npm),jhw(nbnbx,nqbz,npm)
      real(8):: whw(nhwtot)
      complex(8) :: zmelt2      !...........................sf 21May02
c         complex(8), allocatable :: zxq_(:,:,:) !...........sf 21May02
      real(8) :: imagweight !............................sf 21May02
c         logical :: takao=.false. !.........................sf 21May02
c         allocate( zxq_( nbloch + ngc,nbloch + ngc,nwt)) !..sf 21May02
      integer(4)::nocc
      real(8):: eband(nband)
c      integer(4):: n_index_qbz
c      integer(4):: index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)
c-----
      integer(4):: nlnm(*),nlnmv(*),nlnmc(*),iclass(*),icore(*),ncore(*)
      integer(4):: verbose
c---for iepsmode
      logical   :: nolfco !iepsmode
      integer(4):: nmbas, imb1,imb2, nmbas1,nmbas2,imb
      complex(8):: zq01,zq02,zzr(ngbb,nmbas)
c      real(8)   :: zq0zq0
      complex(8)   :: zq0zq0 !This is a bug for the case of two atoms per cell!!! oct2006
      complex(8):: rcxqmean(nwt,npm,nmbas1,nmbas2)

      real(8):: qu1(3),qu2(3),qu1x(3),qu2x(3),quu(3),tolqu=1d-4

      integer(4):: nbcut,nbcut2
      logical :: iww1=.true.,iww2=.true.

      logical:: smbasis
      integer(4):: ifpomat, nbloch_r, ngbo,iqxdummy
c      integer(4),allocatable:: io(:),in(:),io_q(:),in_q(:)
      complex(8),allocatable:: pomat(:,:), zmeltn(:,:,:)
      real(8):: q_r(3)
      complex(8):: img=(0d0,1d0)

      integer(4):: nkmin,  nkmax, nkqmin, nkqmax
      real(8):: qq(3)
      integer(4):: ib1, ib2,      ngcx,ix,iy

      integer(4),allocatable:: ngvecc(:,:)
      logical :: onceww !testtr,negative_testtr

!! takao apr2012
      logical :: newanisox, z1stcol,z1offd
      complex(8) :: ppovlz(ngbb,ngbb),zzz(ngbb)
      integer:: igb

!! jun2012takao
c      real(8):: qeibz(3,nqbz),  !    aik(3,3,ngrpt) 
      integer:: ngrp,nwgt(nqbz) !,ngrpt, aiktimereversal(ngrpt),nwgtieibz,ieibz
      integer:: igx(ngrp*2,nqbz),igxt(ngrp*2,nqbz),ieqbz
!! nwgt(neibz
      logical:: eibz4x0,checkbelong,eibzmode
      complex(8):: ppovl(ngc,ngc) ,zcousq(ngbb,ngbb) 
      complex(8),allocatable:: zcousqr(:,:,:),rcxq0(:,:),rcxq00(:,:),rcxq000(:,:),rcxqwww(:,:)
c      complex(8):: zcousqrsum(ngbb,ngbb,2),    zmeltx(ngbb),zmelty(ngbb),zcousqrx(ngbb,ngbb)  
      complex(8):: zmeltx(ngbb),zmelty(ngbb),zcousqrx(ngbb,ngbb),zcousqc(ngbb,ngbb)
     & ,rzc(ngbb,ngbb),zrrx(ngbb,ngbb),zcousqinv(ngbb,ngbb),cmat(ngbb,ngbb)
      integer::  eibzsym(ngrp,-1:1),neibz,icc,ig,eibzmoden,ikp,i,j,itimer,icount,iele
      integer:: irotm,nrotmx,ixx,iyy,itt,ntimer, nccc, nxx,iagain,irotm1,irotm2
      integer,allocatable:: i1(:,:),i2(:,:),nrotm(:)
      complex(8),allocatable:: zrotm(:,:),zrr(:,:),zrrc(:,:)
c ------------------------------------------------


cooo ccccccccccccccccccccccccccccc
c      complex(8),allocatable :: ppovlzinv(:,:)
c      allocate(ppovlzinv(ngbb,ngbb))
c      ppovlzinv = ppovlz
c      call matcinv(ngbb,ppovlzinv)
cooo end cccccccccccccccccccccccccccccccccccccccccccc




      write(*,'(" x0kf_v4hz: q=",3f8.4,$)')q
      call cputid(0)
c      if(verbose()>=80) debug=.true.

      imdim(1) = 1
      do iatom = 1,natom
        iatomp(iatom) = iatom
        if(iatom<natom) imdim(iatom+1)=imdim(iatom)+mdim(iatom)
      enddo
c dimensions and constants !     nlmtoband  = nlmto*nband
      nctot      = noccx - noccxv
      call dinv33(qbas,0,qbasinv,det)
      phase= (1d0,0d0) !coskt = 1d0; sinkt = 0d0
      allocate(cphik (nlmto,nband),cphikq(nlmto,nband),
     &         geig1(ngpmx,nband),  geig2(ngpmx,nband) )

C ---- overlaid ntp0 ---------- kq 21Jun2003 takao
      call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
      call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
      call getnemx(nbmx,ebmx,7,.true.)

      if(smbasis()) then
ccccccccccccccccccccccccccccccccc
        if(ncc/=0) then
          print *,"Timereversal=F(ncc/=0) not yet implemented for smbasis."
          print *," pomat should be generated correctly                  ."
          stop "Timereversal=F(ncc/=0) not yet implemented for smbasis."
        endif
ccccccccccccccccccccccccccccccccc
        ifpomat = iopen('POmat',0,-1,0) !oct2005
C... smoothed mixed basis !oct2005
C This replace original zmelt with new zmelt based on smoothed mixed basis.
        do
          read(ifpomat) q_r,nn,no,iqxdummy !readin reduction matrix pomat
c          write(6,"('ttt: q  =',3f12.5)") q
c          write(6,"('ttt: q_r=',3f12.5)") q_r
          allocate( pomat(nn,no) )
          read(ifpomat) pomat
          if( sum(abs(q-q_r))<1d-10) then ! .and.kx <= nqibz ) then
            write(6,*) 'ok find the section for give qibz_k'
            exit
!         elseif (kx >nqibz ) then
!           exit
          endif
          deallocate(pomat)
        enddo
        if( sum(abs(q-q_r))>1d-10 ) then
          write(6,"('q  =',3f12.5)") q
          write(6,"('q_r=',3f12.5)") q_r
          stop 'POmat reading err q/=q_r'
        endif
        isx = iclose('POmat')
      endif

c$$$ccccccccccccccccccccccccccccccccccccccc
c$$$      if(testtr()) then
c$$$        allocate(ngvecc(3,ngc))
c$$$        call readqg('QGcou', q,ginv, qu2x, ngcx, ngvecc)
c$$$        if(ngc/=ngcx) then
c$$$           write(6,*) 'ngc ngcx =',ngc,ngcx
c$$$           stop ' ngc/=ngcx'
c$$$        endif
c$$$        write(9005,*) ngc
c$$$        do igc=1,ngc
c$$$          write(9005,"(3i5)") ngvecc(1:3,igc)
c$$$        enddo
c$$$      endif
c$$$ccccccccccccccccccccccccccccccccccccccc

      eibzmode = eibz4x0().and.(.not.nolfco)
      eibzmoden = 2
      if(.not.eibzmode) eibzmoden=0

C... loop over k-points
      qq=0d0
      do 1000 k = 1,nqbz
        if(eibzmode.and.nwgt(k)==0 ) cycle

        if(debug) write(*,'("do 1000  k=",i4)')k
        ipr=(k<5.or.k==nqbz.or.debug)
        if(sum(nbnb(k,1:npm))==0) cycle
        if(k<=5) then
          write(6,"(' x0kf_v4hz: k rk=',i4,3f10.4,$)")k, rk(:,k)
          call cputid(0)
        elseif(mod(k,10)==1.or.k==nqbz) then
          write(6,"(' x0kf_v4hz: k rk=',i4,3f10.4)")k, rk(:,k)
        endif
c     


C... Read cphi part of eigenfunctions for k and k+q
c        kisp= isp + nsp*(k-1)
c        call rwdd1   (ifcphi, kisp, 2*nlmto*noccxv, cphik)
c        kpisp= isp + nsp*(kp-1)
c        call rwdd1   (ifcphi, kpisp, 2*nlmtoband,  cphikq)
        call  readcphi(  rk(:,k)-qq, nlmto,isp1,   quu,cphik) !quu is used q for eigenfunctions.
        call  readcphi(q+rk(:,k)-qq, nlmto,isp2,   quu,cphikq)

C --- tetra ------ override nt0 itps ntp0 ---------------
c        nkmax  = 0;  nkqmin = nband + nctot   ! for sanity check
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
          enddo
        enddo
        call readeval(q+rk(:,k)-qq,isp2,eband)
        nkqmax = nocc (eband,ebmx,.true.,nband)
        if(npm==2) then
          call readeval(rk(:,k)-qq,isp1,eband)
          nkmax = nocc (eband,ebmx,.true.,nband)
        endif

        itps  = nkqmin          ! nkqmin = the num of min   n2 =unocc for jpm=1
        nt0   = nkmax
        ntp0  = nkqmax - nkqmin +1
        if( npm==2.and. nkqmin/=1) then
          write(6,*)' npm==2 nkqmin nkqmax  nkmin nkmax=',nkqmin,nkqmax,nkmin,nkmax
          stop " When npm==2, nkqmin==1 should be."
        endif
        if(nkmin/=1) then
          stop " nkmin==1 should be."
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

        if(ipr) then
          write(6,"(' nkRange  nkqRange=',2i6,2x,2i6)") nkmin,nkmax,nkqmin,nkqmax
        endif

C --- calculate the matrix elements <psi(k+q,t') | psi(k,t) B(R,i)>
C ... core
c      write (*,*) 'matrix elements'
        if(debug)  call cputid(0)
c        allocate( zzmel(nbloch,noccx,ntp0) )
c                        q        k          q+k
        if(debug)  write(6,*)nbloch,nctot,nt0,ncc,ntp0
        allocate( zzmel(nbloch, nctot+nt0, ncc+ntp0) )
        if(debug)  then
          write(6,"('zzzw nkmin nkqmin=',2i5)") nkmin,nkqmin
        endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(onceww(5) ) print *,' nctot ncc=',nctot,ncc
c        stop 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        call psicb_v3  ( nctot,ncc,nt0,ntp0,iclass,phase,
     i              cphik (1,nkmin),  
     i              cphikq(1,nkqmin), 
     i              ppb,            
     i              nlnmv,nlnmc,mdim, 
     i              imdim,iatomp,
     i              mdimx,nlmto,nbloch,nlnmx,natom,nclass,
     i              icore,ncore,nl,nnc,
     o              zzmel)
C ... valence
        !call cputid(0)
        call psi2b_v3  ( nctot,ncc,nt0,ntp0,iclass,phase,
     i              cphik(1,nkmin), 
     i              cphikq(1,nkqmin),
     i              ppb,         
     i              nlnmv, nlnmc,mdim,
     i              imdim,iatomp,
     d              mdimx,nlmto,nbloch,nlnmx,natom,nclass,
     o              zzmel)
        if(debug) call cputid(0)

C --- IPW set
        call readqg('QGpsi',q+rk(:,k)-qq,ginv, qu1x, ngp1, ngvecpB1)
        call readqg('QGpsi',  rk(:,k)-qq,ginv, qu2x, ngp2, ngvecpB2)

c      ngp1 = ngpn(kp)  ! q+k   ntp0 in FBZ
c      ngp2 = ngpn(k)   ! k     np0  in FBZ
!     ngc              ! q          in IBZ
        ngb  = nbloch + ngc ! This is not ngbb for smbasis()=T. oct2005

c                        q       k        q+k
        allocate( zmelt(ngb,nctot+nt0,ncc+ntp0) )
cooo
c        allocate( zmelt1(ngb,nctot+nt0,ncc+ntp0) )
        allocate( z1p(ngb,ngb) )
C ... read plane wave part of eigenfunction.
        call readgeig(q+rk(:,k)-qq, ngpmx,isp2, qu1, geig1)
        call readgeig(  rk(:,k)-qq, ngpmx,isp1, qu2, geig2)

cccccccccccccc
c      write(6,"('uuu sumtest g1a=',6f8.4,2i4,2d13.5)")
c     &  q+rk(:,k), qu1,itps,itps+ntp0-1, sum(geig1(1:ngp1, itps:itps+ntp0-1))
cccccccccccc
c      print *,'uuu sumtest g1a=',sum(geig1(1:ngp1, itps:itps+ntp0-1))
c      geig1(1:ngp1, itps:itps+ntp0-1)= geigB(1:ngp1, itps:itps+ntp0-1,kp)
c      print *,' sumtest g1b=',sum(geig1(1:ngp1, itps:itps+ntp0-1))
c      print *,'uuu sumtest g2a=',sum(geig2(1:ngp2, 1:nt0) )
c      geig2(1:ngp2, 1:nt0     )= geigB(1:ngp2, 1:nt0,    k )
c      print *,' sumtest g2b=',sum(geig2(1:ngp2, 1:nt0) )
        if(sum(abs(qu1-qu1x))>tolqu) then
          write(6,"('qu1 :',3d23.10)") qu1
          write(6,"('qu1x:',3d23.10)") qu1x
          stop 'x0kf_v4hz:qu1/=qu1x'
        endif
        if(sum(abs(qu2-qu2x))>tolqu) then
          write(6,"('qu2 :',3d23.10)") qu2
          write(6,"('qu2x:',3d23.10)") qu2x
          stop 'x0kf_v4hz:qu2/=qu2x'
        endif

c     qdiff = q   - qbkp(:) + rk(:,k)
        qdiff = q    - qu1    + qu2
        ! q   - (q+k)   + k  is not zero.
        ! qc  -  q1     + q2
        add   = matmul(qbasinv, qdiff)
        nadd  = idint( add + dsign(.5d0,add))  !  print *,' qdif=',qdiff,qbkp(:),rk(:,k)
        if(sum(abs(add-nadd))>1d-10) stop "sexc: abs(add-nadd))>1d-10"
        zmelt = 0d0
        if(ngc/=0) then !Aug2005
          call melpln2(  ngp1, ngvecpB1  ! q1=q+k  ; kp ngp1 1:ntp0 q-point
     &             , ngp2, ngvecpB2  ! q2=k    ; k  ngp2 1:nt0  occupied
c     &           , ngc,  ngveccB,  nadd,
     &           , ngc,  nadd,
     &       geig1(1:ngp1, itps:itps+ntp0-1), ntp0, ! q1=q+k  ; kp
     &       geig2(1:ngp2, 1:nt0     ),  nt0, ! q2=k    ; k
     i       shtv,  q, q,  symope,qbas,
C... extensiton to nbloch+ngc       
     o       zmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,ncc+1:ncc+ntp0))
        endif
        zmelt(1:nbloch, 1:nctot+nt0, 1:ncc+ntp0) = 
     &  zzmel(1:nbloch, 1:nctot+nt0, 1:ncc+ntp0)
!                         k            q+k
        deallocate(zzmel)
        if(debug) write(*,'("4 z1pp definitions begin",3d13.6,$)') sum(abs(zmelt))
        if(debug) call cputid(0)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$cc      nt0=4
c$$$cc      ntp0=22
c$$$cc      itps=5
c      write(6,"('uuu scheck= ',6f9.4, 8i4, 4d13.6)")
c     & q+rk(:,k),rk(:,k),
c     &  ngb,ngc,nctot,nt0,ntp0,ngp1,ngp2,itps
c     & ,sum(abs(geig1(1:ngp1, itps:itps+ntp0-1)))
c     & ,sum(abs(geig2(1:ngp2, 1:nt0     )))
c     & ,sum(abs(zmelt(1:ngb,1:nctot+nt0,1:ntp0)))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C...Feb2006 We have gotten zmelt until here. -------------------------------
!   zmelt(1:nbloch, 1:nctot+nt0, 1:ncc+ntp0)
!   ,where ncc=0(npm=1) or ncc=nctot(npm=2)
!    Note
!          1:nctot are  for core
!         nctot+1:nctot+nt0 are for valence
!    n1b,n2b has different ordering---core is after valence.
!... --- end of matrix elements --------------------------------------------
        if(smbasis()) then !
           stop ' not maintained recently need to check again' !takao for safe jun2012
C... smoothed mixed basis !oct2005
C This replace original zmelt with new zmelt based on smoothed mixed basis.
          if(nn /= ngbb) stop 'x0kf_v4h: ngb is not ngbb'
          if(no /= ngb)  stop 'x0kf_v4h: POmat err no/=ngb'
          ngb = nn       ! Renew ngb
          allocate( zmeltn(ngb, 1:nctot+nt0, 1:ntp0) )
          call matm( pomat, zmelt, zmeltn, nn,no,(nctot+nt0)*ntp0)
          deallocate(zmelt)
          allocate( zmelt(ngb, 1:nctot+nt0, 1:ntp0) )
          zmelt = zmeltn
          deallocate(zmeltn)
c        write(6,"('smbasis: ngb ngbb nn no=',4i5)") ngb,ngbb,nn,no
        endif

!! === zmelt conversion === \mu-basis representation.
c        if(newanisox.and.eibzmode) then
c          if(ngb/=ngbb) then
c            print *,' cx0kf_v4h: ngb ngbb=',ngb,ngbb
c            stop 'x0kf_v4h: ngb/=ngbb'
c          endif
c          do it=1,nctot+nt0
c          do itp=1,ncc+ntp0  !! See sxcf_fal2.F L933 around. zmelt1
c            zmelt(nbloch+1:nbloch+ngc,it,itp)=matmul(zmelt(nbloch+1:nbloch+ngc,it,itp),ppovl)
c          enddo
c          enddo
c        elseif(newanisox) then !.and.iq==1) then
        if(newanisox) then !.and.iq==1) then
          if(ngb/=ngbb) then
            print *,' cx0kf_v4h: ngb ngbb=',ngb,ngbb
            stop 'x0kf_v4h: ngb/=ngbb'
          endif
          do it=1,nctot+nt0
          do itp=1,ncc+ntp0  !! See sxcf_fal2.F L933 around. zmelt1
            zmelt(:,it,itp) =  matmul( zmelt(:,it,itp),ppovlz(:,:) )
          enddo
          enddo
        endif


c$$$cooo test mode
c$$$          do it=1,nctot+nt0
c$$$          do itp=1,ncc+ntp0
c$$$            do imb=1,ngb !max(nmbas1,nmbas2)
c$$$              zzz(imb) =  sum( zmelt(1:ngb,it,itp)*ppovlzinv(1:ngb,imb) )
c$$$            enddo
c$$$            zmelt(:,it,itp) =  zzz(:)
c$$$c            zmelt(:,it,itp) =  matmul( zmelt(:,it,itp),ppovlzinv(:,:) )
c$$$
c$$$c            do imb=1,ngb !max(nmbas1,nmbas2)
c$$$c              zzz(imb)=0d0
c$$$c              do igb=1,ngb
c$$$c              zzz(imb)=zzz(imb)+zmelt(igb,it,itp)*ppovlzinv(igb,imb)
c$$$c              enddo
c$$$c            enddo
c$$$c            zmelt(:,it,itp) =  zzz(:)
c$$$          enddo
c$$$          enddo

!! New section Dec 2003.
        if(nolfco) then !iepsmode==202) then ! just for <exp(iq r)|x0(q,\omega)|exp(iq r)>
          if(.not.newanisox) then
            do it=1,nctot+nt0
              do itp=1,ncc+ntp0
                do imb=1,max(nmbas1,nmbas2)
                  zmelt(imb,it,itp) =  sum(zmelt(:,it,itp)*zzr(:,imb) )
                enddo
              enddo
            enddo
          endif

          do jpm  = 1,npm !
            do ibib = 1, nbnb(k,jpm) !---  ibib loop
              if(jpm==1) then
                if( n1b(ibib,k,jpm) <= nbcut.and. n2b(ibib,k,jpm)>nbcut2) then !oct2005
                  if(iww2) then
                    write(6,"(' nband_chi0 nbcut nbcut2 n2b n1b=',4i6)") nbcut,n2b(ibib,k,jpm),n1b(ibib,k,jpm)
                    iww2=.false.
                  endif
                  cycle
                endif
              else
                if( n2b(ibib,k,jpm) <= nbcut.and. n1b(ibib,k,jpm)>nbcut2) then !oct2005
                  if(iww2) then
                    write(6,"(' nband_chi0 nbcut nbcut2 n2b n1b=',4i6)") nbcut,n2b(ibib,k,jpm),n1b(ibib,k,jpm)
                    iww2=.false.
                  endif
                  cycle
                endif
              endif

              if( jpm==1.and.n2b(ibib,k,jpm) > nbmx) cycle
              if( jpm==2.and.n1b(ibib,k,jpm) > nbmx) cycle

              if( n1b(ibib,k,jpm) <= nband) then
                it = nctot + n1b(ibib,k,jpm) !valence
              else
                it = n1b(ibib,k,jpm) - nband !core
              endif

              if( n2b(ibib,k,jpm) <= nband) then
                itp = ncc + n2b(ibib,k,jpm) - itps + 1 !val
                if(itp > ncc + nkqmax-itps+1 ) cycle
              else
                itp =       n2b(ibib,k,jpm) - itps + 1 - nband !core
              endif

              do imb1=1,nmbas1
                zq01 = zmelt(imb1,it,itp)
                z1stcol=(nmbas1==2.and.imb1==2.and.z1offd) !a specical treatment 
                if(z1stcol) zq01 = zmelt(1,it,itp)
                do imb2=1,nmbas2
                  zq02 = zmelt(imb2,it,itp)
                  zq0zq0 = dconjg(zq01)*zq02
                  if(z1stcol) zq0zq0 = dconjg(zq0zq0)
                  do iw=ihw(ibib,k,jpm),ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1!1,nhw(ibib,k) !iw is in range 1:nwhis
                    if (iw .gt. nwt) stop "x0kf_v4hz: iw > nwt"
                    !iiww=iw+ihw(ibib,k)-1
                    imagweight   = whw(jhw(ibib,k,jpm)+iw-ihw(ibib,k,jpm))
                    rcxqmean(iw,jpm,imb1,imb2) =  ! here we  sum over ibib (or n, n') and k.
     &              rcxqmean(iw,jpm,imb1,imb2) 
     &              + zq0zq0*imagweight !sum over spin in hx0fp0
c                    print *,'rrrrr: ',k,iw,imagweight,zq0zq0
                  enddo ! iw
                enddo ! imb2
              enddo ! imb1
            enddo ! ----- ibib loop
          enddo ! ----- jpm loop
          deallocate(zmelt,z1p)
          cycle
        endif

c----------------------------
c  z1p = <M_ibg1 psi_it | psi_itp> < psi_itp | psi_it M_ibg2 >
c  zxq(iw,ibg1,igb2) = sum_ibib wwk(iw,ibib)* z1p(ibib, igb1,igb2)
c        write(6,*)" nkmax ",nkmax,nkqmax,itps
        do 25 jpm  = 1, npm !
        do 25 ibib = 1, nbnb(k,jpm) !---  ibib loop
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
                print *,' nband_chi0 nbmx=',nbmx
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
                print *,' nband_chi0 nbmx=',nbmx
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
            do igb2=1, ngb !....................................
              zmelt2 = zmelt(igb2,it,itp)
            do igb1=1,igb2
              z1p(igb1,igb2) = dconjg(zmelt(igb1,it,itp)) * zmelt2
            enddo
            enddo
c$$$        endif

          if (ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1 >nwt) stop "x0kf_v4hz: iw>nwt"
!$OMP parallel do  private(imagweight) shared(rcxq)
          do iw = ihw(ibib,k,jpm),ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1   !iiww=iw+ihw(ibib,k)-1
            imagweight = whw(jhw(ibib,k,jpm)+iw-ihw(ibib,k,jpm))
            if(eibzmoden==2) imagweight = nwgt(k)*imagweight  
            do igb2=1,ngb  !this part dominates cpu time most time consuming...........
            do igb1=1,igb2
              rcxq(igb1,igb2,iw,jpm) =  !here we  sum over ibib (or n, n') and k.
     &        rcxq(igb1,igb2,iw,jpm) + z1p(igb1,igb2)*imagweight !sum over spin in hx0fp0
             enddo !igb1
             enddo !igb2
          enddo ! iw
!$OMP end parallel do

 25     continue
        if(debug) call cputid(0)
        deallocate(zmelt,z1p)

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

c       if(ipr) then
c         do iw= 1,nwt
c           aaa = 0d0
c         do ix= 1,ngb
c           aaa= aaa + rcxq (ix,ix, iw)
c         enddo
c         write(6,"(' diag sum iw rcxq=',i5,2d13.5)") iw,aaa
c         enddo
c       endif
        if(debug) write(6,*)' end of kloop k jpm=',k,jpm
        if(ipr) write(6,*)
        !write(*,'("2 after ibib",$)') !..............sf 21May02
        !    call cputid(0)
 1000 continue

!! ==== hermitianize. jun2012takao moved from dpsion5 ====
      if(.not.nolfco) then
        do jpm=1,npm
        do iw= 1,nwt
         do igb2= 1,ngb
         do igb1= 1,igb2-1
           rcxq(igb2,igb1,iw,jpm) = dconjg(rcxq(igb1,igb2,iw,jpm))
         enddo
         enddo
        enddo
        enddo
      endif

!! ==symmetrizer of EIBZ PRB.81,125102(2010) Eq.(51) july2012takao==
      if(eibzmoden==2) then
        print *,'eibzmoden=2 goto iqindx2'
        call iqindx2(q, ginv, qtt_, nqnum, ikp,quu) !to get ikp for timereversal mode
        if(sum(abs(q-quu))>1d-8) stop 'x0kf_v2h: eibz 111 q/quu'
        neibz=sum(eibzsym(:,1))+sum(eibzsym(:,-1)) !itimer=-1 means time reversal. eibzsym(ig,itimer) where ig: space rotation.
        write(6,"(' --- goto symmetrization --- ikp neibz q=',2i3,3f12.8)")ikp,neibz,q
        call cputid2(' --- qqqqq: start symmetrization  ',0)

        allocate(rcxq0(ngb,ngb),rcxq00(ngb,ngb),rcxq000(ngb,ngb),rcxqwww(ngb,ngb))
        ntimer=1
        if(sum(eibzsym(:,-1))>0) ntimer=2 !timereversal case
        allocate(zrotm(ngb,ngb),nrotm(ngrp*2))
!!
        zcousqinv=zcousq
        call matcinv(ngb,zcousqinv)
        nrotmx = 10000 !trial value
 1011   continue !only when nrotmx overflow.
        if(allocated(i1)) deallocate(i1,i2,zrr,zrrc)
        nccc=ngrp*2
        allocate(i1(nrotmx,nccc),i2(nrotmx,nccc),zrr(nrotmx,nccc),zrrc(nrotmx,nccc))
        i1=-99999
        i2=-99999
        zrr=-99999d0
        zrrc=-99999d0
        call cputid2(' --- qqqqq111   :',0)
!!
        icc=0
        do itimer=1,-1,-2
          if(ntimer==1.and.itimer==-1) exit
          if(itimer==1 ) itt=1
          if(itimer==-1) itt=2
          do ig=1,ngrp
            if(eibzsym(ig,itimer)==1) then
              icc=icc+1
              call rotMPB2(nbloch,ngb,q,ig,itimer,ginv,zrotm) 
!!this is time-consuming. we may reduce compuatational time by sparse matrix multiplication.
!! takao tried in a Si case. but matmmsparse was not so efficient...
!! zrrx is the rotation matrix for \mu based representaion of mixed basis. PRB.81,125102(2010) Eq.(51)
              if(itimer==1) then
                call matmmsparse(zcousqinv,zrotm,zcousq,zrrx,ngb,1d-8,iele)
c                zrrx= matmul(zcousqinv,matmul(zrotm, zcousq)) 
              else  
                call matmmsparse(dconjg(zcousqinv),dconjg(zrotm),zcousq,zrrx,ngb,1d-8,iele)
c                zrrx= matmul(dconjg(zcousqinv),matmul(dconjg(zrotm), zcousq)) 
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
                nrotmx=irotm+10000 !enlarge allocation and do things again.
                print *,' warn:(slow speed) xxxx goto 1011 xxxxxx nrotmx+=nrotmx+10000 again'
                goto 1011
              endif   
              nrotm(icc)=irotm
              print *,'ig itimer icc nrotm=',ig,itimer,icc,nrotm(icc) ,iele
            endif
          enddo
        enddo

!! === main part to obtain symmetrized rcxq  ===
c        call cputid(0)
        zcousqc=dconjg(transpose(zcousq))
        call cputid2(' --- qqqqq222ini:',0)
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

               do irotm1 = 1,nrotm(icc)
                 if(abs(zrr(irotm1,icc))<1d-8) cycle
                 rcxq0(:,i2(irotm1,icc)) =rcxq0(:,i2(irotm1,icc)) + rcxqwww(:,i1(irotm1,icc)) * zrr(irotm1,icc)
               enddo
               do irotm2 = 1,nrotm(icc)
                 if(abs(zrrc(irotm2,icc))<1d-8) cycle
                 rcxq00(i2(irotm2,icc),:)= rcxq00(i2(irotm2,icc),:) + zrrc(irotm2,icc) * rcxq0(i1(irotm2,icc),:)
               enddo

c               if(itimer==1) then
c                 rcxq000 = rcxq000 + rcxq00
c               else  
c                 rcxq000 = rcxq000 + transpose(rcxq00)
c               endif

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
              endif
            enddo

c$$$           if(itimer==1) then
c$$$             rcxq000(:,:) = matmul(zcousqc,matmul(rcxq00,zcousq))
c$$$c$$$c               call zgemm("N","N",ngb,ngb,ngb, (1d0,0d0), rcxq00, ngb, zcousq,ngb, (0d0,0d0), rzc,ngb)
c$$$c$$$c               call zgemm("N","N",ngb,ngb,ngb, (1d0,0d0), zcousqc,ngb, rzc,ngb, (0d0,0d0), rcxq000,ngb)
c$$$           elseif(icount>0) then  
c$$$c$$$c           print *,'qqqqq icount=',icount
c$$$c$$$c           rcxq000(:,:) = rcxq000(:,:) + transpose(matmul(transpose(zcousq),matmul(rcxq00,dconjg(zcousq))))
c$$$             rcxq000(:,:) = rcxq000(:,:) +   matmul(matmul(zcousqc,transpose(rcxq00)),zcousq)
c$$$           endif

            if(itimer==1) then
               rcxq000=rcxq00
            else
               rcxq000=rcxq000+rcxq00
            endif   
          enddo
          rcxq(:,:,iw,jpm) = rcxq000/neibz
        enddo
        enddo
        deallocate(rcxq00,rcxq000,rcxq0,zrotm,i1,i2,rcxqwww)
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
      deallocate(cphik,cphikq,geig1,geig2)
      if(smbasis()) deallocate(pomat)
      print *,' x0kf_v4hz: end'
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
        write(*,'(" dpsion5: start...   ",$)')
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
          ispx = 3-isp !flip
        endif

C... Check freqr=frhis_m.
C... But I think now this is not necessary. You can supply any freqr(iw). But be careful.
        if(realomega) then
          if( nwhis <= nw_w ) then
            print *,nwhis,nw_w
            stop ' dpsion5: nwhis<=nw_w'
          endif
          if( freqr(0)/=0d0 ) stop ' dpsion5: freqr(0)/=0d0'
          aaa = 0d0
          if(nw_w>0) then
            do iw = 1,nw_w
              aaa = aaa + abs( freqr(iw) - (frhis(iw)+frhis(iw+1))/2d0 )
              if(debug) write(6,"(' iw freqr frhis_m=',i5,2f13.6)" )
     &     iw,freqr(iw),  (frhis(iw)+frhis(iw+1))/2d0
            enddo
            if(aaa>1d-10)stop 'dpsion5:freqr/=frhis_m is not implimented yet'
          endif
        endif !realomega

C--------------------------------------------------------------
! Each bin  [his_Left, his_Right]  his_Center is middle.
! his_C(0) is at zero. his_R(0) and his_L(0) are not defined.
        if(debug) print *,' dpsion5: RRR 2222222222 '
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
c         if(debug) print *, 'dpsion5: RRR 7777 iw wfac=',iw,wfac,ecut,ecuts
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
            enddo ;   if(debug) print *, 'dpsion5: RRR 55555555555'
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
            enddo ;   if(debug) print *, 'dpsion5: RRR2 55555555555'
            deallocate(rrr)
          endif
          rmatt = rmatt/pi ; if(debug) print *,'dpsion5: RRR 6666'

!! takao remove if(nolfc) block here.
          write(6,*) " --- realomega dgemm--- "


!! WARN! I think npm==2.and.chipm does not make sense. apr2012.
!!
          if(npm==2.and.chipm) 
     &      stop 'x0kf_v4h:npm==2.and.chipm is not meaningful probably'


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
            stop 'dpsion5: npm=1 or 2'
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
          endif ;   if(debug) print *, 'dpsion5: III 111111155555555555'
          do it =  1,niwt
            zz = img*freqi(it) !his_C(it)
            call hilbertmat(zz, nwhis,his_L,his_C,his_R, rrr) !Im(zz)>0
            if(npm==1) then
              rmati (it,:,1) = dreal(rrr)
            else
              rmatiC(it,:,1) = rrr
            endif
          enddo ;   if(debug) print *, 'dpsion5: III 55555555555'
!! ==== npm=1 case ====
          if(npm==1) then
            allocate( imatt(niwt, nwhis,npm) )
            do iw= 1,nwhis
              imatt(:,iw,1) = rmati(:,iw,1) - rmati(:,-iw,1)
            enddo
            deallocate(rmati,rrr)
            imatt = imatt/pi; if(debug) print *, 'dpsion5: III  '
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
            imattC = imattC/pi; if(debug) print *, 'dpsion5: IIIc '
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
        write(*,'("         end dpsion5 ",$)')
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

