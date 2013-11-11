      subroutine x0kf_v2hz (
     i            ihw,nhw,jhw,whw,nhwtot, ! tetwt5
     i            n1b,n2b,nbnbx,nbnb,     !  use whw by tetwt5 ,
     i                  q,
c     i                  ifev, 
     i                  nsp,isp1,isp2,      !ifcphi,ifrb,ifcb,ifrhb,ifchb,
     i                  qbas,ginv,rk,wk,  !,iindxk,
c     i            index_qbz,n_index_qbz,
c     i                  mdim,inlnm,inlnmv,inlnmc,iiclass,
     i                  mdim, nlnm, nlnmv,nlnmc, iclass, ! not the w(pointer)
     i                  ppb, !pdb,dpb,ddb,
c     i                  freqr,
c     i                  iicore,incore,
     i                  icore,ncore, ! not the w(pointer)
c     &                  ecore,
c
     d                  nlmto,nqbz,noccx,noccxv,
     d                  nl,nclass,natom,nnc,
     d                  nlnmx,mdimx,nbloch,nwt,
c
c     i      iq,  ngpn,ngc, ngpmx_in,ngcmx,  geigB, jan2004
c     i       ngvecpB,ngveccB,
     i      iq, ngbb, ngc, ngpmx_in,ngcmx,  !oct2005 ngbb is not ngc + nbloch for smbasis=T.
c     i      ngveccB,
     i      nqbze, nband,nqibz,
c
     o      rcxq,
c------------------------------------------
     i     nolfco,zzr, !iepsmode,zzr,
     o     rcxqmean)
      use m_readqg, only: readqg
      use m_readeigen, only: readeval,readcphi,readgeig
      use keyvalue
c     o     zxq)
c- Only the spectrum function (imaginary part) along the real axis.
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
c
      implicit real*8(a-h,o-z)
      dimension q(3),qbas(3,3),ginv(3,3),rk(3,nqbz),wk(nqbz)
c     i          freqr(nwt),
c     w          ek(nband),ekq(nband),ekc(noccx), ecore(noccx - noccxv)
c     complex(8):: zxq (nbloch+ngc,nbloch+ngc, nwt)
c      complex   (8):: rcxq (nbloch+ngc,nbloch+ngc, nwt)

      integer(4):: ngbb
      complex   (8):: rcxq (ngbb,ngbb, nwt),aaa

c      dimension qk(3) !,qbkp(3)
      complex(8) :: imag=(0d0,1d0),trc
      complex(8),allocatable:: cphik(:,:),cphikq(:,:),geig1(:,:),geig2(:,:)

c      integer w
c      common /w/ w(1)

      integer(4):: ngpmx_in, ngcmx, nqbze, nband,
     &              ngc,nadd(3), !ngvecpB(3,ngpmx_in,nqbze), ngpn(nqbze),
     &             igc, !ngveccB(3,ngcmx),
     &             ngvecpB1(3,ngpmx_in),ngvecpB2(3,ngpmx_in)

c      complex(8):: geigB  (ngpmx_in,nband,nqbze)
c      complex(8):: geigB  (1,1,1)

      complex(8),allocatable :: zmelt(:,:,:)
c
      real(8) :: qbasinv(3,3), det,qdiff(3),add(3),symope(3,3),shtv(3)

      data symope /1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0/
      data shtv   /0d0,0d0,0d0/

      real(8) :: ppb(*)
      real(8) :: coskt(natom), sinkt(natom),qx(3)
      integer(4) :: mdim(natom)
      complex(8),allocatable :: ppovl(:,:),ttx(:,:)

      complex(8),allocatable::  z1p(:,:)
c      real(8)   ,allocatable:: rc1p(:,:)

c      complex(8) :: zviw

c      complex(8) ::  wwk(nbnbx,nqbz,nwt)
      integer(4) ::  nbnb(nqbz), 
     &  n1b(nbnbx,nqbz), n2b(nbnbx,nqbz)

      complex(8),allocatable:: zzmel(:,:,:)
      integer(4)::imdim(natom),iatomp(natom)
      logical:: debug=.false.,ipr

c---tetwt5
      logical:: hist
      integer(4):: nbnbx,nqbz,nhwtot,
     &  ihw(nbnbx,nqbz),nhw(nbnbx,nqbz),jhw(nbnbx,nqbz)
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
      complex(8):: zq0,zzr(ngbb)
      real(8)   :: zq0zq0
      complex(8):: rcxqmean(nwt)

c      complex(8),allocatable:: geigBk(:,:,:),geigbkp(:,:,:)
c      logical:: ifgeigb,iffk,iffkp
c      character*3 charnum3

      real(8):: qu1(3),qu2(3),qu1x(3),qu2x(3),quu(3)
c      complex(8)::   geig1(ngpmx_in,nband),  geig2(ngpmx_in,nband)

      integer(4):: nbcut,nbcut2
      logical :: iww1=.true.,iww2=.true.

      logical:: smbasis
      integer(4):: ifpomat, nbloch_r, ngbo,iqx
c      integer(4),allocatable:: io(:),in(:),io_q(:),in_q(:)
      complex(8),allocatable:: pomat(:,:), zmeltn(:,:,:)
      real(8):: q_r(3)

c ------------------------------------------------
      write(*,'("---  begin x0kf_v2hz xxx--- ",$)')
      if(debug) write(6,*) "x0kf_v2hz aaaa0"
      call cputid(0)
      if(verbose()>=80) debug=.true.

      imdim(1) = 1
      do iatom = 1,natom
        iatomp(iatom) = iatom
        if(iatom<natom) imdim(iatom+1)=imdim(iatom)+mdim(iatom)
      enddo
      if(debug) write(6,*) "x0kf_v2hz aaaa1"
c dimensions and constants
      nlmtoband  = nlmto*nband
      nctot      = noccx - noccxv
      call dinv33(qbas,0,qbasinv,det)
      if(debug) write(6,*) "x0kf_v2hz aaaa2"
c work arrays
      coskt = 1d0
      sinkt = 0d0
      allocate(cphik (nlmto,nband),cphikq(nlmto,nband),
     &         geig1(ngpmx,nband),  geig2(ngpmx,nband) )
      if(debug) write(6,*) "x0kf_v2hz aaaa3"
C ---- overlaid ntp0 ---------- kq 21Jun2003 takao
      call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
      call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
      call getnemx(nbmx,ebmx,7,.true.)

      if(smbasis()) then
        ifpomat = iopen('POmat',0,-1,0) !oct2005
C... smoothed mixed basis !oct2005
C This replace original zmelt with new zmelt based on smoothed mixed basis.
        do
          read(ifpomat) q_r,nn,no,iqx !readin reduction matrix pomat
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
Cstop2rx 2013.08.09 kino          stop 'POmat reading err q/=q_r'
          call rx( 'POmat reading err q/=q_r')
        endif
        isx = iclose('POmat')
      endif

C... loop over k-points
      do 1000 k = 1,nqbz
        if(debug) write(*,'("1  begin k-cycle",$)')
        if(debug) call cputid(0)
        ipr=(k<5.or.k==nqbz.or.debug)

        if(nbnb(k)==0) cycle

        if(k<=5) then
          write(6,"(' x0kf_v2hz: k rk=',i4,3f10.4,$)")k, rk(:,k)
          call cputid(0)
        elseif(mod(k,10)==1.or.k==nqbz) then
          write(6,"(' x0kf_v2hz: k rk=',i4,3f10.4)")k, rk(:,k)
        endif

C... Read cphi part of eigenfunctions for k and k+q
c        kisp= isp + nsp*(k-1)
c        call rwdd1   (ifcphi, kisp, 2*nlmto*noccxv, cphik)
c        kpisp= isp + nsp*(kp-1)
c        call rwdd1   (ifcphi, kpisp, 2*nlmtoband,  cphikq)
        call  readcphi(  rk(:,k),nlmto,isp1,   quu,cphik)
        call  readcphi(q+rk(:,k),nlmto,isp2,   quu,cphikq)

C --- tetra ------ override nt0 itps ntp0 ---------------
        n1bmax  = 0; n2bmin = nband + nctot +1000000  ! for sanity check
        do ibib = 1, nbnb(k)
          if( n1b(ibib,k)<=nband) n1bmax = max(n1b(ibib,k),n1bmax)
          n2bmin = min(n2b(ibib,k),n2bmin)
        enddo
        nt0  = n1bmax           ! n1bmax = the num of maxmum occ
        itps = n2bmin           ! n2bmin = the num of min unocc
        ntp0 = nband - (itps -1)

c       call rwdd1  (ifev,kp,    nband,eband)
        call readeval(q+rk(:,k),isp2,eband)
        n2bmax = nocc (eband,ebmx,.true.,nband)
        ntp0   = n2bmax - (itps -1)
        if(ipr) write(6,"(' n1bmax n2bmin n2bmax=',3i10)")n1bmax,n2bmin,n2bmax
c        do ib = 1, nband
c          write(6,"(i3,f13.5)") ib,eband(ib)
c        enddo

C --- calculate the matrix elements <psi(k+q,t') | psi(k,t) B(R,i)>
C ... core
c      write (*,*) 'matrix elements'
        if(debug)  write(*,'("2  matrix elements",$)')
        if(debug)  call cputid(0)
        allocate( zzmel(nbloch,noccx,ntp0) )
        call psicb_v2  (icore,ncore,ntp0,iclass,
     i              coskt,sinkt,
     i              cphikq(1,itps), !rbkq(1,itps),cbkq(1,itps),rhbkq(1,itps),chbkq(1,itps),
     i              ppb,            !,pdb,dpb,ddb,
     i              nlnmv,nlnmc,mdim, 
     i              imdim,iatomp,
     d              mdimx,nlmto,nbloch,nlnmx,noccx,ntp0,natom,nclass,
     d              nl,nnc,
     o              zzmel)
C ... valence
        !write(*,'("2.1 2nd matrix elements",$)')
        !call cputid(0)
        call psi2b_v2  (nt0,ntp0,iclass,coskt,sinkt,
     i              cphik,         ! rbk,cbk,rhbk,chbk,
     i              cphikq(1,itps),! rbkq(1,itps),cbkq(1,itps),rhbkq(1,itps),chbkq(1,itps),
     i              ppb,           !,pdb,dpb,ddb,
     i              nlnmv, nlnmc,mdim,nctot,
     i              imdim,iatomp,
     d              mdimx,nlmto,nbloch,nlnmx,noccxv,noccx,ntp0,
     d              natom,nclass,
     o              zzmel)
        if(debug) write(*,'("3 after matrix elements",$)')
        if(debug) call cputid(0)

C --- IPW set
        call readqg('QGpsi',q+rk(:,k),ginv, qu1x, ngp1, ngvecpB1)
        call readqg('QGpsi',  rk(:,k),ginv, qu2x, ngp2, ngvecpB2)
c      ngp1 = ngpn(kp)  ! q+k   ntp0 in FBZ
c      ngp2 = ngpn(k)   ! k     np0  in FBZ
!     ngc              ! q          in IBZ
        ngb  = nbloch + ngc ! This is not ngbb for smbasis()=T. oct2005

        allocate( zmelt(ngb,nctot+nt0,ntp0) )
        allocate( z1p(ngb,ngb) )

C ... read plane wave part of eigenfunction.
        call readgeig(q+rk(:,k), ngpmx,isp2, qu1, geig1)
cccccccccccccc
c      write(6,"('uuu sumtest g1a=',6f8.4,2i4,2d13.5)")
c     &  q+rk(:,k), qu1,itps,itps+ntp0-1, sum(geig1(1:ngp1, itps:itps+ntp0-1))
cccccccccccc
c      print *,'uuu sumtest g1a=',sum(geig1(1:ngp1, itps:itps+ntp0-1))
c      geig1(1:ngp1, itps:itps+ntp0-1)= geigB(1:ngp1, itps:itps+ntp0-1,kp)
c      print *,' sumtest g1b=',sum(geig1(1:ngp1, itps:itps+ntp0-1))
        call readgeig(  rk(:,k), ngpmx,isp1, qu2, geig2)
c      print *,'uuu sumtest g2a=',sum(geig2(1:ngp2, 1:nt0) )
c      geig2(1:ngp2, 1:nt0     )= geigB(1:ngp2, 1:nt0,    k )
c      print *,' sumtest g2b=',sum(geig2(1:ngp2, 1:nt0) )
Cstop2rx 2013.08.09 kino        if(sum(abs(qu1-qu1x))>1d-10) stop 'x0kf_v2hz:qu1/=qu1x'
        if(sum(abs(qu1-qu1x))>1d-10) call rx( 'x0kf_v2hz:qu1/=qu1x')
Cstop2rx 2013.08.09 kino        if(sum(abs(qu2-qu2x))>1d-10) stop 'x0kf_v2hz:qu2/=qu2x'
        if(sum(abs(qu2-qu2x))>1d-10) call rx( 'x0kf_v2hz:qu2/=qu2x')

c     qdiff = q   - qbkp(:) + rk(:,k)
        qdiff = q    - qu1    + qu2
        ! q   - (q+k)   + k  is not zero.
        ! qc  -  q1     + q2
        add   = matmul(qbasinv, qdiff)
        nadd  = idint( add + dsign(.5d0,add))  !  print *,' qdif=',qdiff,qbkp(:),rk(:,k)
Cstop2rx 2013.08.09 kino        if(sum(abs(add-nadd))>1d-10) stop "sexc: abs(add-nadd))>1d-10"
        if(sum(abs(add-nadd))>1d-10) call rx( "sexc: abs(add-nadd))>1d-10")
        zmelt = 0d0
        if(debug) write(*,'("333 goto melpln2t",$)')
        if(ngc/=0) then !Aug2005
          call melpln2t(  ngp1, ngvecpB1  ! q1=q+k  ; kp ngp1 1:ntp0 q-point
     &             , ngp2, ngvecpB2  ! q2=k    ; k  ngp2 1:nt0  occupied
c     &           , ngc,  ngveccB,  nadd,
     &           , ngc,  nadd,
     &       geig1(1:ngp1, itps:itps+ntp0-1), ntp0, ! q1=q+k  ; kp
     &       geig2(1:ngp2, 1:nt0     ),  nt0, ! q2=k    ; k
     i       shtv,  q, q,  symope,qbas,
     i       q, !qt oct2013
C... extensiton to nbloch+ngc       
     o       zmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0))
        endif
        zmelt  (1:nbloch, 1:nctot+nt0, 1:ntp0) = 
     &  zzmel(1:nbloch, 1:nctot+nt0, 1:ntp0)
!                         k            q+k
        deallocate(zzmel)
        if(debug) write(*,'("4 z1pp definitions begin",$)')
        if(debug) call cputid(0)

ccccccccccccccccccccccccccccccccccccccccc
c      write(6,"('uuu scheck= ',6f9.4, 7i4, 3d13.6)")
c     & q+rk(:,k),rk(:,k),
c     &  ngb,ngc,nctot,nt0,ntp0,ngp1,ngp2
c     & ,sum(abs(geig1(1:ngp1, itps:itps+ntp0-1)))
c     &, sum(abs(geig2(1:ngp2, 1:nt0     )))
c     & ,sum(abs(zmelt(1:ngb,1:nctot+nt0,1:ntp0)))
ccccccccccccccccccccccccccccccccccccccccc

!... --- end of matrix elements --------------------------------------------------

        if(smbasis()) then !
C... smoothed mixed basis !oct2005
C This replace original zmelt with new zmelt based on smoothed mixed basis.
Cstop2rx 2013.08.09 kino          if(nn /= ngbb) stop 'x0kf_v2h: ngb is not ngbb'
          if(nn /= ngbb) call rx( 'x0kf_v2h: ngb is not ngbb')
Cstop2rx 2013.08.09 kino          if(no /= ngb)  stop 'x0kf_v2h: POmat err no/=ngb'
          if(no /= ngb)  call rx( 'x0kf_v2h: POmat err no/=ngb')
          ngb = nn       ! Renew ngb
          allocate( zmeltn(ngb, 1:nctot+nt0, 1:ntp0) )
          call matm( pomat, zmelt, zmeltn, nn,no,(nctot+nt0)*ntp0)
          deallocate(zmelt)
          allocate( zmelt(ngb, 1:nctot+nt0, 1:ntp0) )
          zmelt = zmeltn
          deallocate(zmeltn)
c       write(6,"('smbasis: ngb ngbb nn no=',4i5)") ngb,ngbb,nn,no
        endif

cccccccccccccccccccccccccccccccccccccccccccccccccc
ctakao test
c        zmelt(1:nbloch, nctot+1:nctot+6, 1:ntp0) = 0d0
cccccccccccccccccccccccccccccccccccccccccccccccccc

!New section Dec 2003.
        if(nolfco) then !iepsmode==202) then ! just for <exp(iq r)|x0(q,\omega)|exp(iq r)>
          do ibib = 1, nbnb(k) !---  ibib loop

            if( n1b(ibib,k) <= nbcut.and. n2b(ibib,k)>nbcut2) then !oct2005
              if(iww2) then
                write(6,"(' nband_chi0 nbcut nbcut2 n2b n1b=',4i6)") nbcut,n2b(ibib,k),n1b(ibib,k)
                iww2=.false.
              endif
              cycle
            endif

            if( n2b(ibib,k) > nbmx) cycle
            if( n1b(ibib,k) <= nband) then
              it = nctot + n1b(ibib,k) !valence
            else
              it = n1b(ibib,k) - nband !core
            endif
            itp = n2b(ibib,k) - itps + 1
            if(itp> n2bmax-itps+1) cycle
            zq0 = sum(zmelt(1:ngb,it,itp)*zzr(1:ngb))
            zq0zq0 = dconjg(zq0)*zq0
c$$$cccccccccccccccccccccccccccc
c$$$            if(k/=10) cycle
c$$$            write(6,"('uuu',14i5)")k,ibib,n1b(ibib,k),n2b(ibib,k), it,itp
c$$$            write(6,"('uuu zq0zq0 sumzzr=',2d24.16)")zq0zq0,sum(abs(zzr(1:ngb)))
c$$$            write(6,"('uuu sumzmelt =',2d24.16)") sum(zmelt(1:ngb,it,itp))
c$$$            write(6,"('uuu sumzmelt1=',2d24.16)") sum(zmelt(1:nbloch,it,itp))
c$$$            write(6,"('uuu sumzmelt2=',2d24.16)") sum(zmelt(nbloch+1:ngb,it,itp))
c$$$           write(6,"('uuu ng chk=',18i5)")
c$$$     &       ngp1,  sum(abs(ngvecpB1(:,1:ngp1))), ngp2,sum(abs(ngvecpB2(:,1:ngp2)))
c$$$     &       ,  ngc,   nadd,  ntp0, nt0
c$$$           write(6,"('uuu q chk=',3f9.4,9f9.4)") q(:), qbas(:,:)
c$$$           write(6,"('uuu itps ntp0=',18i5)") itps,itps+ntp0-1
c$$$           write(6,"('uuu geig=',d23.16)")
c$$$     &        sum(abs(geig1(1:ngp1, itps:itps+ntp0-1))),
c$$$     &        sum(abs(geig2(1:ngp2, 1:nt0)))
c$$$cccccccccccccccccccccccccccc
            do iw=ihw(ibib,k),ihw(ibib,k)+nhw(ibib,k)-1!1,nhw(ibib,k) !iw is in range 1:nwhis
Cstop2rx 2013.08.09 kino              if (iw .gt. nwt) stop "x0kf_v2hz: iw > nwt"
              if (iw .gt. nwt) call rx( "x0kf_v2hz: iw > nwt")
              !iiww=iw+ihw(ibib,k)-1
              imagweight   = whw(jhw(ibib,k)+iw-ihw(ibib,k))
              rcxqmean(iw) =  ! here we  sum over ibib (or n, n') and k.
     &        rcxqmean(iw) + zq0zq0*imagweight !sum over spin in hx0fp0
            enddo ! iw
          enddo ! ----- ibib loop
          deallocate(zmelt,z1p)
          cycle
        endif

c----------------------------
c  z1p = <M_ibg1 psi_it | psi_itp> < psi_itp | psi_it M_ibg2 >
c  zxq(iw,ibg1,igb2) = sum_ibib wwk(iw,ibib)* z1p(ibib, igb1,igb2)
c
        do 25 ibib = 1, nbnb(k) !---  ibib loop
          if( n2b(ibib,k) > nbmx) then  !nbmx
            if(iww1) then
              print *,' nband_chi0 nbmx=',nbmx
              iww1=.false.
            endif
            cycle
          endif
          if( n1b(ibib,k) <= nbcut.and. n2b(ibib,k)>nbcut2) then
            if(iww2) then
              write(6,"(' nband_chi0 nbcut nbcut2 n2b n1b=',4i6)") nbcut,n2b(ibib,k),n1b(ibib,k)
              iww2=.false.
            endif
            cycle
          endif

          if(n1b(ibib,k) <= nband) then
            it = nctot + n1b(ibib,k) !valence
          else
            it = n1b(ibib,k) - nband !core
          endif
          itp  = n2b(ibib,k) - itps + 1
c-------
          if(itp> n2bmax-itps+1) cycle !takao 21June2003

cccccccccccccccccccccccccccccccccc takao variant begin
cc         if(takao) then
cc
cc         do ic = 1,ngb
cc           z1p(1:ngb,ic) =
cc     &     zmelt(ic,it,itp)*dconjg(zmelt(1:ngb,it,itp))
cc         end do
cc         ihww = ihw(ibib,k)

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

          !.....................................................................sf 21May02
ccccccccccccccccccccccccccccccc faleev variant begin,it about 15% faster for MnO case

          do igb2=1, ngb !....................................
            zmelt2 = zmelt(igb2,it,itp)
            ! if (abs(zmelt2) .lt. 1d-4) cycle
            do igb1=1,igb2
              z1p(igb1,igb2) = 
     &       dconjg(zmelt(igb1,it,itp)) * zmelt2
            enddo
          enddo

Cstop2rx 2013.08.09 kino          if (ihw(ibib,k)+nhw(ibib,k)-1 >nwt) stop "x0kf_v2hz: iw>nwt"
          if (ihw(ibib,k)+nhw(ibib,k)-1 >nwt) call rx( "x0kf_v2hz: iw>nwt")
CCC!$OMP parallel do  private(imagweight) shared(rcxq)
          do iw=ihw(ibib,k),ihw(ibib,k)+nhw(ibib,k)-1
            !iiww=iw+ihw(ibib,k)-1
            imagweight=whw(jhw(ibib,k)+iw-ihw(ibib,k))
            do igb2=1,ngb  !this part dominates cpu time...........
              do igb1=1,igb2
c     call zdaxpy(igb2,imagweight,z1p(1,igb2),zxq_(1,igb2,iiww))
                rcxq(igb1,igb2,iw)=  !here we  sum over ibib (or n, n') and k.
     &         rcxq(igb1,igb2,iw)+z1p(igb1,igb2)*imagweight !sum over spin in hx0fp0
              enddo !igb1
            enddo !igb2
          enddo ! iw
CCC!$OMP end parallel do
c           do iw=ihw(ibib,k),ihw(ibib,k)+nhw(ibib,k)-1!1,nhw(ibib,k) !iw is in range 1:nwhis
c             if (iw .gt. nwt) stop "x0kf: iw > nwt"
c              !iiww=iw+ihw(ibib,k)-1
c             imagweight=whw(jhw(ibib,k)+iw-ihw(ibib,k))
c             do igb2=1,ngb  !this part dominates cpu time...........
c                   ! if (abs(zmelt(igb2,it,itp)) .lt. 1d-4) cycle
c             do igb1=1,igb2
cc     call zdaxpy(igb2,imagweight,z1p(1,igb2),zxq_(1,igb2,iiww))
c               rcxq(igb1,igb2,iw)=  !here we  sum over ibib (or n, n') and k.
c     &         rcxq(igb1,igb2,iw)+z1p(igb1,igb2)*imagweight !sum over spin in hx0fp0
c             enddo !igb1
c             enddo  !igb2
c           enddo ! iw
cc            endif
cccccccccccccccccccccccccccccccccccccccccc faleev variant end
 25     continue
        if(debug) write(*,'("6 mine cycles end",$)')
        if(debug) call cputid(0)

        !write(*,'("2 after ibib",$)') !............................sf 21May02
        ! call cputid(0)
cc        if(.not.takao) then !......................sf 21May02
c
c        do iw=1,nwt
c           write(*,*)iw,nwt,ngb
c          !!!call  rcxq_zxq(rc1p,zxq_(1:ngb,1:ngb,iw),ngb,-1)
c        do igb2=1,ngb
c          do igb1=1,igb2
c             rcxq(igb1,igb2,iw)=rcxq(igb1,igb2,iw)+ zxq_(igb1,igb2,iw)
c                                             !rc1p(:,igb2)
c           enddo
c        enddo
c        enddo !....................................sf 21May02
cc        endif
        deallocate(zmelt,z1p)
        if(ipr) write(6,"(' sum rcxq=', 2d23.15,2i8)") sum(rcxq),ngb,ngbb
c       if(ipr) then
c         do iw= 1,nwt
c           aaa = 0d0
c         do ix= 1,ngb
c           aaa= aaa + rcxq (ix,ix, iw)
c         enddo
c         write(6,"(' diag sum iw rcxq=',i5,2d13.5)") iw,aaa
c         enddo
c       endif
        if(debug) write(6,*)' end of kloop k=',k
        if(ipr) write(6,*)
        !write(*,'("2 after ibib",$)') !..............sf 21May02
        !    call cputid(0)

 1000 continue
c      print *,' x0kf_v2hz: 111'
c      deallocate(zxq_) !...........................sf 21May02
      deallocate(cphik,cphikq,geig1,geig2)
      if(smbasis()) deallocate(pomat)
      print *,' x0kf_v2hz: end'
cccccccccccccccccccccccccccccccccccccccccccc
c        print *,' ---- trace chck ---- '
c        do iw= 1,nwt
c        do ix=1,5
c          write(6,"(a,2i4,3d13.5)")' ttt=',iw,ix,
c     &    rcxq(ix,ix+1,iw),rcxq(ix,ix+1,iw),rcxq(ix+1,ix,iw)
c        enddo
c        enddo
c        stop 'test zzzzzzzzzzz1'
cccccccccccccccccccccccccccccccccccccccccccc
      end
